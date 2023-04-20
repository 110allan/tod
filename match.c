#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include "ketopt.h" // command-line argument parser
#include "kthread.h" // multi-threading models: pipeline and multi-threaded for loop
#include "kstring.h"
#include "kvec.h"
#include "khashl.h" // hash table
KHASHL_SET_INIT(static, uint64_set_t, uint64_set, uint64_t, kh_hash_uint64, kh_eq_generic)

KHASHL_MAP_INIT(static, uint64_map_t, uint64_map, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)


#include "kseq.h" // FASTA/Q parser

KSEQ_INIT(gzFile, gzread)

#include "misc.h"
;

typedef struct {
    int block_len;
    int n_thread;
    int only_count;
    FILE *out;
    int matrix_out;
	int text_out;
} tod_opt_t;


static inline void tod_opt_init(tod_opt_t *opt)
{
    memset(opt, 0, sizeof(tod_opt_t));
    opt->out = stdout;
    opt->block_len = 100000000;
    opt->n_thread = 12;
    opt->only_count = 0;
    opt->matrix_out = 0;
}

/**************************************************************  MATCH  ************************************************************************/
typedef struct {
    const tod_opt_t *p;
    int k;
    uint64_t mask, shift;
    kseq_t *ks;
    uint64_set_t *h;
} match_pl_t;

typedef struct { // data structure for each step in kt_pipeline()
    match_pl_t *p;
    int n, m, sum_len;
    int *len;
    char **seq;
    char **seqname;
    char **ret;
    int *retlen;
} match_step_t;

static uint64_set_t   *_load_kmers_hash(const char *fn, int *kmer_size)
{
    uint64_set_t *h = 0;
	gzFile fp;
	if((fp = (fn && strcmp(fn, "-") ? gzopen(fn, "r") : gzdopen(fileno(stdin), "r"))) == 0) {
        return NULL;
    }
    
    kstream_t *ks = ks_init(fp);
    kstring_t str = {0, 0, 0};
    ks_resize(&str, 256);
    if(kmer_size) {
        *kmer_size = 0;
    }
    if(ks_getuntil(ks, KS_SEP_SPACE, &str, 0) >= 0) {
        int k = str.l;
        int i, l;
        uint64_t kmer, krev;
        if(kmer_size) {
            *kmer_size = k;
        }
        if(*kmer_size > 32) {
            return 0;
        }

        h = uint64_set_init();
        uint64_t mask = (1ULL << k * 2) - 1, shift = (k - 1) * 2;
        do {
            for(i = l = 0, kmer = krev = 0; i < k; ++i) {
                int c = seq_nt4_table[(uint8_t)str.s[i]];
                int absent;
                if(c < 4) {
                    kmer = (kmer << 2 | c) & mask;                  // forward strand
                    krev = krev >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
                    if(++l >= k) {  // we find a k-mer
                        uint64_set_put(h,
                                       (kmer < krev ? kmer : krev),
                                       &absent); // only add one strand!
                    }
                } else {
                    l = 0, kmer = krev = 0;
                }

            }
        } while(ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0);
    }

    ks_destroy(ks);
    gzclose(fp);
    free(str.s);
    return h;
}


static void _match_worker_for(void *data, long idx, int tid) // callback for kt_for()
{
    match_step_t *s = (match_step_t *)data;
    char *seq = s->seq[idx];
    int len = s->len[idx];
    char *seqname = s->seqname[idx];
    int match_nk = 0;
    kvec_t(int) match_pos;
    kv_init(match_pos);
    kv_resize(int, match_pos, 1024);
    uint64_set_t *uniq_kmers = uint64_set_init();

    uint64_set_t *h = s->p->h;
    int k = s->p->k;
    uint64_t mask = s->p->mask, shift = s->p->shift;
    int i, l;
    uint64_t kmer, krev;
    for(i = l = 0, kmer = krev = 0; i < len; ++i) {
        int c = seq_nt4_table[(uint8_t)seq[i]];
        if(c < 4) {
            kmer = (kmer << 2 | c) & mask;                  // forward strand
            krev = krev >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
            if(++l >= k) {  // we find a k-mer
                khint_t itr;
                itr = uint64_set_get(h, kmer < krev ? kmer : krev);
                if(itr != kh_end(h)) {
                    kv_push(int, match_pos, i + 1 - k);
                    match_nk++;
                    int absent;
                    uint64_set_put(uniq_kmers, itr, &absent);
                }
            }
        } else {
            l = 0, kmer = krev = 0;
        }

    }
    // output
    kstring_t ret_str = {0, 0, 0};
    ks_resize(&ret_str, 1024);
    kputs(seqname, &ret_str);
    kputc('\t', &ret_str);
    kputw(len, &ret_str);
    kputc('\t', &ret_str);
    kputw(kh_size(uniq_kmers), &ret_str);
    kputc('\t', &ret_str);
    kputw(match_nk, &ret_str);
    if((!s->p->p->only_count) && match_nk) {
        kputc('\t', &ret_str);
        kputw(kv_A(match_pos, 0), &ret_str);
        for(int i = 1; i < match_nk; ++i) {
            kputc(',', &ret_str);
            kputw(kv_A(match_pos, i), &ret_str);
        }
    }
    kv_destroy(match_pos);
    uint64_set_destroy(uniq_kmers);
    MALLOC(s->ret[idx], ret_str.l + 1);
    memcpy(s->ret[idx], ret_str.s, ret_str.l + 1);
    s->retlen[idx] = ret_str.l;
    free(ret_str.s);
}

static void *_match_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
    match_pl_t *p = (match_pl_t *)data;
    if(step == 0) {  // step 1: read a block of sequences
        int ret;
        match_step_t *s;
        CALLOC(s, 1);
        s->p = p;

        while((ret = kseq_read(p->ks)) >= 0) {
            int l = p->ks->seq.l;
            if(l < p->k) {
                continue;
            }
            if(s->n == s->m) {
                s->m = s->m < 256 ? 256 : s->m + (s->n >> 1);
                REALLOC(s->len, s->m);
                REALLOC(s->seq, s->m);

                REALLOC(s->seqname, s->m);
                REALLOC(s->ret, s->m);
                REALLOC(s->retlen, s->m);
            }

            //seq
            MALLOC(s->seq[s->n], l);
            memcpy(s->seq[s->n], p->ks->seq.s, l);
            //seqname
            MALLOC(s->seqname[s->n], p->ks->name.l + 1);
            memcpy(s->seqname[s->n], p->ks->name.s, p->ks->name.l + 1);
            s->len[s->n++] = l;
            s->sum_len += l;

            if(s->sum_len >= p->p->block_len) {
                break;
            }
        }
        if(s->sum_len == 0) {
            free(s);
        } else {
            return s;
        }
    } else if(step == 1) {  // step 2: match k-mers to hash table
        match_step_t *s = (match_step_t *)in;
        kt_for(p->p->n_thread, _match_worker_for, s, s->n);
        return s;
    } else if(step == 2) {
        match_step_t *s = (match_step_t *)in;
        //free
        for(int i = 0; i < s->n; ++i) {
            //fprintf(s->p->out, "%s\n", s->ret[i]);
            const int len = s->retlen[i];
            s->ret[i][len] = '\n';
            fwrite(s->ret[i], len + 1, 1, s->p->p->out);
            free(s->seq[i]);
            free(s->seqname[i]);
            free(s->ret[i]);
        }
        free(s->seq);
        free(s->len);
        free(s->seqname);
        free(s->ret);
        free(s->retlen);
        free(s);
    }
    return 0;
}


static void match(const char *kmer_file_fn, const char *fn, const tod_opt_t *opt)
{
    int k = 0;
    uint64_set_t *h = _load_kmers_hash(kmer_file_fn, &k);
    if(!h) {
        fprintf(stderr, "load kmers error \n");
        return;
    } else {
        fprintf(stderr, "load %d kmers (k=%d) \n", kh_size(h), k);
    }

    match_pl_t pl;
    gzFile fp;

	if((fp=(fn && strcmp(fn, "-") ? gzopen(fn, "r") : gzdopen(fileno(stdin), "r"))) == 0) {
        return;
    }

    pl.ks = kseq_init(fp);
    pl.k = k;
    pl.mask = (1ULL << k * 2) - 1;
    pl.shift = (k - 1) * 2;
    pl.h = h;
    pl.p = opt;
    kt_pipeline(3, _match_pipeline, &pl, 3);
    uint64_set_destroy(h);
    kseq_destroy(pl.ks);
    gzclose(fp);
    return;
}

/**************************************************************  END MATCH  ************************************************************************/


/**************************************************************  MATCH2  ************************************************************************/
typedef kvec_t(uint32_t) vec_uint32_t;

typedef struct {
    const tod_opt_t *p;
    int k;
    uint64_t mask, shift;
    kseq_t *ks;
    uint64_t *kmer_table;
	uint32_t kmer_table_size;
	uint32_t dump_id;
	vec_uint32_t dump_cols;
} match2_pl_t;

typedef struct { // data structure for each step in kt_pipeline()
    match2_pl_t *p;
    int n, m, sum_len;
    int *len;
    char **seq;
    char **seqname;
    uint16_t* cnts;//    number_of_kmers rows and nseqs cols
} match2_step_t;

typedef kvec_t(uint64_t) vec_uint64_t;

static uint64_t  *_load_kmers_vec(const char *fn, int *kmer_size, uint32_t *size)
{

    vec_uint64_t vec;
    kv_init(vec);
    kv_resize(uint64_t, vec, 16384);
    if(size) {
        *size = 0;
    }
    gzFile fp ;
	if((fp=(fn && strcmp(fn, "-") ? gzopen(fn, "r") : gzdopen(fileno(stdin), "r"))) == 0) {
        return NULL;
    }

    kstream_t *ks = ks_init(fp);
    kstring_t str = {0, 0, 0};
    ks_resize(&str, 256);
    if(kmer_size) {
        *kmer_size = 0;
    }
    if(ks_getuntil(ks, KS_SEP_SPACE, &str, 0) >= 0) {
        int k = str.l;
        int i, l;
        uint64_t kmer, krev;
        if(kmer_size) {
            *kmer_size = k;
        }
        if(*kmer_size > 32) {
            return 0;
        }


        uint64_t mask = (1ULL << k * 2) - 1, shift = (k - 1) * 2;
        do {
            for(i = l = 0, kmer = krev = 0; i < k; ++i) {
                int c = seq_nt4_table[(uint8_t)str.s[i]];
                int absent;
                if(c < 4) {
                    kmer = (kmer << 2 | c) & mask;                  // forward strand
                    krev = krev >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
                    if(++l >= k) {  // we find a k-mer
                        kv_push(uint64_t, vec, (kmer < krev ? kmer : krev));
                    }
                } else {
                    l = 0, kmer = krev = 0;
                }

            }
        } while(ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0);
    }

    ks_destroy(ks);
    gzclose(fp);
    free(str.s);
    if(size) {
        *size = kv_size(vec);
    }
    return vec.a;
}


static void _match2_worker_for(void *data, long idx, int tid) // callback for kt_for()
{
    match2_step_t *s = (match2_step_t *)data;
    char *seq = s->seq[idx];
    int len = s->len[idx];
    //char *seqname = s->seqname[idx];
    
    //get kmer cnt
    uint64_map_t *kmer_cnts = uint64_map_init();
    int k = s->p->k;
    uint64_t mask = s->p->mask, shift = s->p->shift;
    int i, l;
    uint64_t kmer, krev;
    for(i = l = 0, kmer = krev = 0; i < len; ++i) {
        int c = seq_nt4_table[(uint8_t)seq[i]];
        if(c < 4) {
            kmer = (kmer << 2 | c) & mask;                  // forward strand
            krev = krev >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
            if(++l >= k) {  // we find a k-mer
                khint_t itr;
                int absent;
                itr = uint64_map_put(kmer_cnts, (kmer < krev ? kmer : krev), &absent);  // only add one strand!
                if(absent) {
                    kh_val(kmer_cnts, itr) = 0;
                }
                ++kh_val(kmer_cnts, itr);
            }
        } else {
            l = 0, kmer = krev = 0;
        }
    }


    // output
    uint32_t v;
    for(int i = 0; i < s->p->kmer_table_size; ++i) {
        khint_t itr;
        itr = uint64_map_get(kmer_cnts, s->p->kmer_table[i]);
        if(itr != kh_end(kmer_cnts)) {
			v=kh_val(kmer_cnts, itr);
			s->cnts[i*s->n+idx]=v>UINT16_MAX?UINT16_MAX:v;
        } 
    }

	uint64_map_destroy(kmer_cnts);
}

static void *_match2_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
    match2_pl_t *p = (match2_pl_t *)data;
    if(step == 0) {  // step 1: read a block of sequences
        int ret;
        match2_step_t *s;
        CALLOC(s, 1);
        s->p = p;

        while((ret = kseq_read(p->ks)) >= 0) {
            int l = p->ks->seq.l;
            if(l < p->k) {
                continue;
            }
            if(s->n == s->m) {
                s->m = s->m < 256 ? 256 : s->m + (s->n >> 1);
                REALLOC(s->len, s->m);
                REALLOC(s->seq, s->m);

                REALLOC(s->seqname, s->m);
            }

            //seq
            MALLOC(s->seq[s->n], l);
            memcpy(s->seq[s->n], p->ks->seq.s, l);
            //seqname
            MALLOC(s->seqname[s->n], p->ks->name.l + 1);
            memcpy(s->seqname[s->n], p->ks->name.s, p->ks->name.l + 1);
			
            s->len[s->n++] = l;
            s->sum_len += l;

            if(s->sum_len >= p->p->block_len) {
                break;
            }
        }
        if(s->sum_len == 0) {
            free(s);
        } else {
        	//ret_str
        	s->cnts = calloc(p->kmer_table_size*s->n,sizeof(uint16_t));
			kv_push(uint32_t,p->dump_cols,s->n);
            return s;
        }
    } else if(step == 1) {  // step 2: match k-mers to hash table
        match2_step_t *s = (match2_step_t *)in;
        kt_for(p->p->n_thread, _match2_worker_for, s, s->n);
        return s;
    } else if(step == 2) {
        match2_step_t *s = (match2_step_t *)in;
		char dump_prefix[1024]={0};
		sprintf(dump_prefix, "tod.match.matrix.tmp_file.pid%d.%d", getpid(),p->dump_id++);
		FILE* file;
        if((file = fopen(dump_prefix, "wb")) == NULL){
            fprintf(stderr, " can not write tmp file %s \n", dump_prefix);
            abort();
        }
		fwrite(s->cnts,sizeof(uint16_t), p->kmer_table_size*s->n, file);
		fclose(file);
        //free
        for(int i = 0; i < s->n; ++i) {
            //fprintf(s->p->out, "%s\n", s->ret[i]);
   			
            free(s->seq[i]);
            free(s->seqname[i]);         
        }
		
        free(s->seq);
        free(s->len);
        free(s->seqname);
        free(s->cnts);
        free(s);
    }
    return 0;
}


static void match2(const char *kmer_file_fn, const char *fn, const tod_opt_t *opt)
{
    int k = 0;
    uint32_t size;
    uint64_t *kmer_table = _load_kmers_vec(kmer_file_fn, &k, &size);
    if(!kmer_table) {
        fprintf(stderr, "load kmers error \n");
        return;
    } else {
        fprintf(stderr, "load %d kmers (k=%d) \n", size, k);
    }

    match2_pl_t pl;
    gzFile fp;
	if((fp=(fn && strcmp(fn, "-") ? gzopen(fn, "r") : gzdopen(fileno(stdin), "r"))) == 0) {
        return;
    }

    pl.ks = kseq_init(fp);
    pl.k = k;
    pl.mask = (1ULL << k * 2) - 1;
    pl.shift = (k - 1) * 2;
    pl.kmer_table = kmer_table;
	pl.kmer_table_size  = size;
    pl.p = opt;
	pl.dump_id=0;
	kv_init(pl.dump_cols);
    kt_pipeline(3, _match2_pipeline, &pl, 3);
	//paste temp files
	uint32_t cnt_number=0;
	for(int i=0;i<kv_size(pl.dump_cols);++i ) cnt_number+=kv_A(pl.dump_cols,i);
	if(size && cnt_number) {
		// write header
		if(!opt->text_out) {
			uint32_t t[3] = {k, cnt_number, size};
			fwrite("TOD\2", 1, 4,opt->out);
			fwrite(t, 4, 3, opt->out);
		}
		char dump_prefix[1024]={0};
		FILE* fps[pl.dump_id];
		for(int i=0;i<pl.dump_id;++i) {
			sprintf(dump_prefix, "tod.match.matrix.tmp_file.pid%d.%d", getpid(),i);
			fps[i] = fopen(dump_prefix,"rb");
			if(fps[i]==NULL) {
				fprintf(stderr, " can not open tmp file %s \n", dump_prefix);
            	abort();
			}
		}

		uint16_t buffer[cnt_number];
		kstring_t str={0,0,0};
		ks_resize(&str, 32+1+6*cnt_number);
		for(int i=0;i<size;++i) {
			for(int j=0,offset=0;j<pl.dump_id;++j) {
				fread(buffer+offset, sizeof(uint16_t), pl.dump_cols.a[j], fps[j]);
				offset+=pl.dump_cols.a[j];
			}
			if(!opt->text_out) {
				fwrite(&kmer_table[i], 8, 1, opt->out);
				fwrite(buffer, sizeof(uint16_t), cnt_number  , opt->out);
			} else {
				kmer2seq(str.s, k, kmer_table[i]);
				str.l= k, str.s[str.l] = 0;
		
				for(int j = 0; j < cnt_number; ++j) {
					char c = (j == 0) ? '\t' : ',';
					kputc(c, &str);
					kputw(buffer[j], &str);
				}
				
				kputc('\n', &str);
				
				fwrite(str.s, str.l, 1, opt->out);
			}
		}
		free(str.s);
		for(int i=0;i<pl.dump_id;++i) {
			fclose(fps[i]);
			sprintf(dump_prefix, "tod.match.matrix.tmp_file.pid%d.%d", getpid(),i);
			remove(dump_prefix);
		}
	}
    free(kmer_table);
    kseq_destroy(pl.ks);
	kv_destroy(pl.dump_cols);
    gzclose(fp);
    return;
}

/**************************************************************  END MATCH2  ************************************************************************/


int kmer_match(int argc, char *argv[])
{
    int parse_ok = 1;
    tod_opt_t todopt;
    tod_opt_init(&todopt);
    ketopt_t o = KETOPT_INIT;
    int c;
    static ko_longopt_t longopts[] = {
        { "matrix-out", ko_no_argument, 300 }, //output matrix
        { NULL, 0, 0 }
    };
    while((c = ketopt(&o, argc, argv, 1, "b:t:co:Th", longopts)) >= 0) {
        switch(c) {
            case 'b':
                todopt.block_len = atoi(o.arg);
                break;
            case 'c':
                todopt.only_count = 1;
                break;
			case 'T':
                todopt.text_out = 1;
                break;
            case 't':
                todopt.n_thread = atoi(o.arg);
                break;
            case 'o': {
                if((todopt.out = (strcmp(o.arg, "-")) ? fopen(o.arg, "wb") : stdout) == 0) {
                    fprintf(stderr, "Invalid out file %s \n", o.arg);
                    parse_ok = 0;
                }
                break;
            }
            case 'h': {
                parse_ok = 0;
                break;
            }
            case 300:
                todopt.matrix_out = 1;
                break;

            default: {
                fprintf(stderr, "Invalid option %c \n", c);
                parse_ok = 0;
                break;
            }

        }

    }
    if((parse_ok == 0) || (argc - o.ind != 2)) {
        fprintf(stderr, "Usage: tod map [options] <in.fa> <kmer file> \n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -b INT           block size [%d]\n", todopt.block_len);
        fprintf(stderr, "  -t INT           number of worker threads [%d]\n", todopt.n_thread);
        fprintf(stderr, "  -c               print only a count of sequences\n");
        fprintf(stderr, "  -o STR           output file[stdout]\n");
		fprintf(stderr, "  -T               output in text format instead of binary format \n");
		fprintf(stderr, "  --matrix-out     output file in matrix format\n");		
        return 1;
    }

    if(!todopt.matrix_out) {
        match(argv[o.ind + 1], argv[o.ind], &todopt);
    } else {
        match2(argv[o.ind + 1], argv[o.ind], &todopt);
    }
    fclose(todopt.out);
    return 0;
}

