#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

#include <zlib.h>
#include <stdint.h>
#include "ketopt.h" // command-line argument parser
#include "kvec.h"
#include "kstring.h"
#include "kseq.h"
#include "ksort.h"
#include "misc.h"




// output format, collapse,add(sum),sub,min,max and file_id
// output format, collapse,add(sum),sub,min,max and file_id
#define OUT_NULL 0
#define OUT_COLLAPSE 1
#define OUT_SUM 2
#define OUT_ADD OUT_SUM
#define OUT_SUB 3
#define OUT_MIN  4
#define OUT_MAX 5
#define OUT_VALID_LIMIT OUT_MAX


#define SET_GROUP 0
#define SET_INTERSECT 1
#define SET_UNION 2
#define SET_DIFF 3
#define SET_OP_COUNT SET_DIFF+1
static const char str_set_op[][10] = {"group", "intersect", "union", "diff"};

KSTREAM_INIT(gzFile, gzread, 32768)
;


typedef kvec_t(gzFile) vec_file_t;
typedef kvec_t(int) vec_int_t;
typedef kvec_t(uint16_t) vec_uint16_t;

typedef struct {
	int first;
	int second;
} int_pair_t;
typedef kvec_t(int_pair_t) vec_int_pair_t;

typedef struct {
	int_pair_t min_max_number;// [0,INT_MAX],use negtive value for invalid value
	vec_int_pair_t  min_max_counts;
    uint16_t fill_cnt;
    vec_file_t vec_fps;

    FILE *out;
    int outtype;
    int set_type;
    int text_out;
    uint8_t *present_files, *absent_files,files_cap;
} tod_opt_t;

#define MAX_KMER_SIZE 32

typedef struct {
    uint64_t kmer;
    uint16_t cnts[];
} kc_t;

typedef struct {
    uint32_t idx,nbyte;
	uint8_t *buf;
} heap_t;

#define heap_lt(a, b) (*(uint64_t*)((a).buf) > *(uint64_t*)((b).buf))

KSORT_INIT(kmc, heap_t, heap_lt)

;

#define N_BLOCK_KMER 16384
#define N_BLOCK_MAX 3
#define N_BLOCK_END (N_BLOCK_KMER*N_BLOCK_MAX)

typedef struct {
    kstream_t *ks;
    int fileid, adjust_file_id;
    int ksize, cnt_number;
    int binary_mode;
    //buffer of the file
    kstring_t str;
    bool cache_valid;
    uint8_t *cached_buf;
    bool iseof;
} infile_t;

static inline int _tfile_kc_init(infile_t *infile)
{
    if(infile->str.l && infile->str.s[infile->str.l - 1] == '\n') {
        --infile->str.l;
        if(infile->str.s[infile->str.l - 1] == '\r') {
            --infile->str.l;
        }
        infile->str.s[infile->str.l] = '\0';
    } else {
        if(ks_getuntil2(infile->ks, KS_SEP_LINE, &infile->str, 0, 1) <= 0) {
            fprintf(stderr, "[set] empty file \n");
            return -1;
        }
    }
    if(strrchr(infile->str.s, '\n')) {
        fprintf(stderr, "[set] unknown file \n");
    }

    kstring_t *str = &infile->str;

    int *fields, n;
    fields = ksplit(str, 0, &n);
    int ksize = strlen(str->s + fields[0]);

    if(ksize > MAX_KMER_SIZE) {
        if(fields) {
            free(fields);
        }
        fprintf(stderr, "[set] files with kmer size >%d not allowed \n", MAX_KMER_SIZE);
        return -1;
    }


    vec_uint16_t vec_cnt;
    kv_init(vec_cnt);
    if(n > 1) {
        char *p = str->s + fields[1], *p_last = p;
        uint32_t v;
        for(; *p; ++p) {
            if(*p == ',') {
                *p = 0;
                v = atoi(p_last);
                if(v > UINT16_MAX) {
                    v = UINT16_MAX;
                }
                kv_push(uint16_t, vec_cnt, v);
                p_last = p + 1;
            }
        }

        v = atoi(p_last);
        if(v > UINT16_MAX) {
            v = UINT16_MAX;
        }
        kv_push(uint16_t, vec_cnt, v);
    }

    infile->cnt_number = kv_size(vec_cnt);
    infile->ksize = ksize;
    infile->iseof = false;
    infile->cache_valid = true;
    infile->cached_buf = malloc(8 + sizeof(uint16_t) * infile->cnt_number);
    *(uint64_t *)infile->cached_buf = seq2kmer(str->s + fields[0], ksize);
    for(int i = 0; i < infile->cnt_number; ++i) {
        *((uint16_t *)(infile->cached_buf + 8) + i) = kv_A(vec_cnt, i);
    }
    free(fields);
    kv_destroy(vec_cnt);
    return ksize;
}


static inline int _bfile_kc_init(infile_t *infile)
{
    kstring_t *str = &infile->str;
    int c;
    while(str->l < 16 && ((c = ks_getc(infile->ks)) >= 0)) {
        str->s[str->l++] = c;
    }
    str->s[str->l] = '\0';

    if((str->l != 16) || (strncmp(str->s + 0, "TOD\2", 4))) {
        fprintf(stderr, "[set] file error,not tod binary file \n");
        return -1;
    }

    // kmer_size(32bit) nsample(32bit) number_of_kmers(32bit)

    int ksize = *(uint32_t *)(str->s + 4);
    int cnt_number = *(uint32_t *)(str->s + 8);
    //int nk = *(uint32_t *)(str->s + 12);

    if(ksize > MAX_KMER_SIZE) {
        fprintf(stderr, "[set] files with kmer size >%d not allowed \n", MAX_KMER_SIZE);
        return -1;
    }

    const uint64_t mask = (1ULL << (ksize * 2)) - 1;

    str->l = 0;
    while(str->l < (8 + cnt_number * 2) && ((c = ks_getc(infile->ks)) >= 0)) {
        str->s[str->l++] = c;
    }
    str->s[str->l] = '\0';

    if(str->l != 8 + cnt_number * 2) {
        fprintf(stderr, "[set] empty/trucate file \n");
        return -1;
    }

    vec_uint16_t vec_cnt;
    kv_init(vec_cnt);
    for(int i = 0; i < cnt_number; ++i) {
        kv_push(uint16_t, vec_cnt, *(uint16_t *)(str->s + 8 + i * 2));
    }

    infile->iseof = false;
    infile->ksize = ksize;
    infile->cnt_number = kv_size(vec_cnt);
    infile->cache_valid = true;
    infile->cached_buf = malloc(8 + sizeof(uint16_t) * infile->cnt_number);

    *(uint64_t *)infile->cached_buf = (*(uint64_t *)(str->s + 0)) & mask;
    for(int i = 0; i < infile->cnt_number; ++i) {
        *((uint16_t *)(infile->cached_buf + 8) + i) = kv_A(vec_cnt, i);
    }
    kv_destroy(vec_cnt);
    return ksize;
}

static inline int _tfile_kc_next(infile_t *infile, uint8_t *buf, size_t nkmer)
{
    const int ksize = infile->ksize;
    
    kstring_t *str = &infile->str;
    const int nbyte = sizeof(uint64_t) + infile->cnt_number * sizeof(uint16_t);
    int read_cnt = 0;
    if(nkmer && infile->cache_valid) {
        infile->cache_valid = false;
        memcpy(buf, infile->cached_buf, nbyte);
        nkmer--;
        buf += nbyte;
        read_cnt++;
    }
    while(nkmer--) {
        if(ks_getuntil(infile->ks, KS_SEP_LINE, str, 0) < 0) {
            return read_cnt;
        }

        if(ks_len(str) == 0) {
            return read_cnt;
        }

        *(uint64_t *)buf = seq2kmer(str->s, ksize), buf += sizeof(uint64_t);
        int cnt_number = infile->cnt_number;
        for(char *p = str->s + ksize + 1, *p_last = p; cnt_number; ++p) {
            if((*p == ',') || isspace(*p) || (*p == '\0')) {
                bool skip = *p != ',';
                *p = 0;
                uint32_t v = atoi(p_last);
                if(v > UINT16_MAX) {
                    v = UINT16_MAX;
                }

                *(uint16_t *)buf = v, buf += sizeof(uint16_t);
                cnt_number--;
                p_last = p + 1;
                if(skip) {
                    break;
                }
            }
        }
        read_cnt++;
    }
    return read_cnt;
}


static inline int _bfile_kc_next(infile_t *infile, uint8_t *buf, size_t nkmer)
{
    //const int ksize = infile->ksize;
    const int cnt_number = infile->cnt_number;
    kstring_t *str = &infile->str;
    const int nbyte = sizeof(uint64_t) + cnt_number * sizeof(uint16_t);
    int read_cnt = 0;
    if(nkmer && infile->cache_valid) {
        infile->cache_valid = false;
        memcpy(buf, infile->cached_buf, nbyte);
        nkmer--;
        buf += nbyte;
        read_cnt++;
    }

    int c;
    while(nkmer--) {
        str->l = 0;
        while(str->l < ( sizeof(uint64_t) + cnt_number * sizeof(uint16_t)) && ((c = ks_getc(infile->ks)) >= 0)) {
            str->s[str->l++] = c;
        }
        str->s[str->l] = '\0';
        if(ks_len(str) == 0) {
            return read_cnt;
        }

        if(ks_len(str) != sizeof(uint64_t) + cnt_number * sizeof(uint16_t)) {
            fprintf(stderr, "[set] file trucate \n");
            return -1;
        }

        memcpy(buf, str->s, nbyte), buf += nbyte;
        read_cnt++;
    }
    return read_cnt;
}

static inline int file_kc_init(infile_t *infile)
{
    if(!(infile && infile->ks)) {
        return -1;
    }
    kstring_t *str = &infile->str;
    str->s = 0, str->l = 0, str->m = 0;
    ks_resize(str, 1024);

    int c;

    str->l = 0;
    while(str->l < 4 && ((c = ks_getc(infile->ks)) >= 0)) {
        str->s[str->l++] = c;
    }
    str->s[str->l] = '\0';

    if((str->l != 4) || (strncmp(str->s + 0, "TOD\2", 4))) {
        infile->binary_mode = false;
    } else {
        infile->binary_mode = true;
    }

    return infile->binary_mode ? _bfile_kc_init(infile) : _tfile_kc_init(infile);
}

static inline int file_kc_next(infile_t *infile, uint8_t *buf, size_t nkmer)
{
    return (infile->iseof || nkmer == 0) ? 0 : (infile->binary_mode ? _bfile_kc_next(infile, buf, nkmer) : _tfile_kc_next(infile, buf, nkmer));
}

static inline void file_kc_destroy(infile_t *infile)
{
    free(infile->cached_buf);
    free(infile->str.s);
}

typedef struct st_bufs_info {
    uint32_t valid_i, cap_i, nbyte; // [0,valid)  [0,cap)
    uint8_t *buf;// per file N_BLOCK_KMER*N_BLOCK_MAX size
} buf_t;

typedef struct {
    const tod_opt_t *opt;
    int total_extend_cnt_number;
    int out_cnt_number;
    infile_t *infiles;
    int infile_number;
    uint32_t header_init;
    buf_t *bufs;
} set_op_pl_t;

typedef struct {
    const set_op_pl_t *p;

    uint32_t *begs, *ends;//[beg,end)
    
	uint8_t *out_buf;
	uint32_t out_buf_n;
	uint32_t out_buf_nbyte;
} set_op_step_t;

set_op_step_t *set_op_step_init(const set_op_pl_t *p)
{
    set_op_step_t *s;
    CALLOC(s, 1);//all the in_cnts/out_cnts set 0
    s->p = p;
   
    s->begs = calloc(p->infile_number, sizeof(uint32_t));
    s->ends = calloc(p->infile_number, sizeof(uint32_t));
    return s;
}

void set_op_step_destory(set_op_step_t **ss)
{
    free((*ss)->begs);
    free((*ss)->ends);
    free(*ss);
}

void output_kmc(uint8_t *des, uint64_t kmer, uint16_t *src, uint32_t len, const set_op_pl_t *p)
{
    //output
    int cnt = 0;
    uint16_t *in_cnts = src;
    uint16_t *out_cnts = (uint16_t*)(des + sizeof(uint64_t));
    *(uint64_t *)des = kmer;

    switch(p->opt->outtype) {
        case OUT_ADD://case OUT_SUM:

            for(uint16_t *cnt_p = in_cnts; cnt_p != src + len; cnt_p++) {
                cnt += *cnt_p;
            }
            *out_cnts = cnt > UINT16_MAX ? UINT16_MAX : cnt;
            break;
        case OUT_MIN:

            cnt = *in_cnts;
            for(uint16_t *cnt_p = in_cnts + 1; cnt_p != src + len; cnt_p++) {
                if(cnt > *cnt_p) {
                    cnt = *cnt_p;
                }
            }
            *out_cnts = cnt > UINT16_MAX ? UINT16_MAX : cnt;
            break;
        case OUT_MAX:

            cnt = *in_cnts;
            for(uint16_t *cnt_p = in_cnts + 1; cnt_p != src + len; cnt_p++) {
                if(cnt < *cnt_p) {
                    cnt = *cnt_p;
                }
            }
            *out_cnts = cnt > UINT16_MAX ? UINT16_MAX : cnt;
            break;
        case OUT_SUB:

            cnt = *in_cnts;
            for(uint16_t *cnt_p = in_cnts + 1; cnt_p != src + len; cnt_p++) {
                cnt -= *cnt_p;
            }
            *out_cnts = (cnt > UINT16_MAX) ? UINT16_MAX : (cnt < 0 ? 0 : cnt);
            break;
        case OUT_COLLAPSE:
            memcpy(out_cnts, in_cnts, len * sizeof(uint16_t));
            break;
        case OUT_NULL:
            break;
        default: {
            infile_t *infile = p->infiles + p->opt->outtype - OUT_VALID_LIMIT - 1;
            for(int i = 0; i < infile->cnt_number; ++i) {
                *(out_cnts + i) = in_cnts[infile->adjust_file_id + i];
            }
            break;
        }
    }

}

static void *set_op_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
    set_op_pl_t *p = (set_op_pl_t *)data;
    if(step == 0) {
        set_op_step_t *s = set_op_step_init(p);
        const int file_cnt = p->infile_number;

        // read N_BLOCK_KMER kmers and update cap_i, kmers may not be used in the last block
        for(int i = 0; i < file_cnt; ++i) {
            //fprintf(stderr, "[---] file id valid_i cap_i=(%u %u %u)->", i, p->bufs[i].valid_i, p->bufs[i].cap_i);

            if(p->infiles[i].iseof || (p->bufs[i].cap_i != p->bufs[i].valid_i)) {
                continue;
            }

            uint8_t *buf = p->bufs[i].buf + p->bufs[i].cap_i * p->bufs[i].nbyte;
            int read_cnt = file_kc_next(p->infiles + i, buf, N_BLOCK_KMER);
            if(read_cnt < N_BLOCK_KMER) {
                p->infiles[i].iseof = true;
            }
            if(read_cnt > 0) {
                p->bufs[i].cap_i += read_cnt;
            }

            if(p->bufs[i].cap_i == N_BLOCK_END) {
                p->bufs[i].cap_i = 0;
            }

            //fprintf(stderr, "(%u %u %u)\n", i, p->bufs[i].valid_i, p->bufs[i].cap_i);
        }

        // get minimal  kmer
        uint64_t minimal_kmer = UINT64_MAX, kmer;
        for(int i = 0; i < file_cnt; ++i) {
            if(p->bufs[i].cap_i != p->bufs[i].valid_i) {
                const int rh_circular = p->bufs[i].cap_i ? (p->bufs[i].cap_i - 1) : (p->bufs[i].cap_i - 1 + N_BLOCK_END);
                kmer = *(uint64_t *)(p->bufs[i].buf + rh_circular * p->bufs[i].nbyte);
                if(kmer < minimal_kmer) {
                    minimal_kmer = kmer;
                }
            }
        }

        // update valid_i
        for(int i = 0; i < file_cnt; ++i) {
            s->begs[i] = p->bufs[i].valid_i;
            // search minimal kmer in last block
            int lh = p->bufs[i].valid_i;
            int rh = p->bufs[i].cap_i - 1;
            if(rh < lh) {
                rh += N_BLOCK_END;
            }

            kmer = *(uint64_t *)(p->bufs[i].buf + (rh % (N_BLOCK_END)) * p->bufs[i].nbyte);
            if(kmer == minimal_kmer) {  // check right boundary
                s->ends[i] = p->bufs[i].cap_i;
                p->bufs[i].valid_i = s->ends[i] ;
                //fprintf(stderr, "step file %u =[%u %u)\n", i, s->begs[i], s->ends[i]);
                continue;
            }

            while(lh <= rh) {  //binary search
                int mid = (lh + rh) / 2;
                int mid_circular = mid % (N_BLOCK_END);

				__builtin_prefetch ((uint64_t*)(p->bufs[i].buf+(mid + 1 + rh)/2), 0, 1);
	        	__builtin_prefetch ((uint64_t*)(p->bufs[i].buf+(lh + mid - 1)/2), 0, 1);

                kmer =  *(uint64_t *)(p->bufs[i].buf + mid_circular * p->bufs[i].nbyte);
                if(kmer >  minimal_kmer) {
                    rh = mid - 1;
                } else {
                    lh = mid + 1;
                }
            }
            if(rh < 0) {
                s->ends[i] = p->bufs[i].valid_i;
                //fprintf(stderr, "step file %u =[%u %u)\n", i, s->begs[i], s->ends[i]);
                continue;
            }

            s->ends[i] = (rh + 1) % (N_BLOCK_END);
            p->bufs[i].valid_i = s->ends[i] ;
            //fprintf(stderr, "step file %u =[%u %u)\n", i, s->begs[i], s->ends[i]);
        }

        bool gotany = false;
        for(int i = 0; i < file_cnt; ++i) {
            if(s->begs[i] != s->ends[i]) {
                gotany = true;
                break;
            }
        }
        if(!gotany) {
            set_op_step_destory(&s);
        } else {
            //fprintf(stderr, "read block num=%d kmers  first kmer=%lu \n", s->cnt, s->kmer_buffer[0]);
            return s;
        }

    } else if(step == 1) {
        set_op_step_t *s = (set_op_step_t *)in;
        uint8_t *pav_files = (uint8_t *)calloc(bit_array_size(p->infile_number), 1);
        int valid_file_cnt = 0;
        heap_t heap[p->infile_number];
        for(int i = 0; i < p->infile_number; ++i) {
            if(s->begs[i] != s->ends[i]) {
                heap[valid_file_cnt].buf = p->bufs[i].buf + s->begs[i] * p->bufs[i].nbyte;
                heap[valid_file_cnt].nbyte = p->bufs[i].nbyte;
                heap[valid_file_cnt].idx = i;
                valid_file_cnt++;
            }
        }

        ks_heapmake(kmc, valid_file_cnt, heap);
        
		s->out_buf_nbyte = (p->out_cnt_number * sizeof(uint16_t) + sizeof(uint64_t));
        s->out_buf = malloc(s->out_buf_nbyte * N_BLOCK_KMER * p->infile_number);
        s->out_buf_n = 0;
		
		kc_t* kc = calloc(1,sizeof(kc_t) + sizeof(uint16_t) * p->total_extend_cnt_number);
		kc->kmer = *(uint64_t *)heap[0].buf;
		
		__builtin_prefetch((int_pair_t*)(s->p->opt->min_max_counts.a), 0, 3);
		__builtin_prefetch((int_pair_t*)(s->p->opt->min_max_counts.a)+1, 0, 3);
		__builtin_prefetch((uint64_t*)kc, 1, 3);
		__builtin_prefetch((uint64_t*)kc+1, 1, 3);

		int pass_number = 0;
		const int min_number=s->p->opt->min_max_number.first,max_number=s->p->opt->min_max_number.second;
        while(true) {
            const uint8_t *buf = heap[0].buf;
            uint64_t kmer = *(uint64_t *)buf;
            if(kmer != kc->kmer || valid_file_cnt == 0) {
                // output
                // check present
                bool check_ok = true;
                for(int i = 0; i < bit_array_size(p->infile_number); ++i) {
                    if(((p->opt->present_files[i] & pav_files[i]) != p->opt->present_files[i])
                       || (p->opt->absent_files[i] & pav_files[i])) {
                        check_ok = false;
                        break;
                    }
                }
                if(check_ok) {
                    if(pass_number < min_number || pass_number > max_number) {
                        check_ok = false;
                    }
                }
                if(check_ok) {
                    output_kmc(s->out_buf + s->out_buf_nbyte * s->out_buf_n,kc->kmer, 
						kc->cnts, p->total_extend_cnt_number, p);
                    s->out_buf_n++;
                }
                if(valid_file_cnt == 0) {
                    break;
                }
                pass_number = 0;
                memset(kc->cnts, 0, sizeof(uint16_t)*p->total_extend_cnt_number);
                memset(pav_files, 0, bit_array_size(p->infile_number));
                kc->kmer = kmer;
            }
			
            const int fileid = heap[0].idx;
            const int nbyte = heap[0].nbyte;
            const infile_t *infile = s->p->infiles + fileid;
            const uint32_t min_count = kv_A(s->p->opt->min_max_counts, fileid).first,
                      max_count = kv_A(s->p->opt->min_max_counts, fileid).second;

			__builtin_prefetch((uint64_t*)(p->bufs[fileid].buf + nbyte), 0, 1);
			
            if(infile->cnt_number > 0) {
                memcpy(kc->cnts + infile->adjust_file_id, buf + 8, nbyte - 8);
            } else {
                kc->cnts[infile->adjust_file_id] = p->opt->fill_cnt;
            }
            
            int all_pass = 1;
            for(int j = 0;
                j < (infile->cnt_number > 0 ? infile->cnt_number : 1) ;
                ++j) {
                const uint32_t cnt = kc->cnts[infile->adjust_file_id + j];
                if(cnt < min_count || (cnt > max_count)) {
                    all_pass = 0;
                    break;
                }
				
            }
            if(all_pass) {
                pass_number++;
				bit_array_set(pav_files, infile->fileid);
            }

            // get next record
            heap[0].buf += heap[0].nbyte;
            if(heap[0].buf == p->bufs[fileid].buf + N_BLOCK_END * heap[0].nbyte) {
                heap[0].buf = p->bufs[fileid].buf;
            }

            if(heap[0].buf == p->bufs[fileid].buf + s->ends[fileid]*heap[0].nbyte) {
                heap[0] = heap[valid_file_cnt-- - 1];
            }

            // adjust heap
            if(valid_file_cnt > 1) {
                ks_heapadjust(kmc, 0, valid_file_cnt, heap);
            }

        }
		free(kc);
        free(pav_files);
		return s;
    }else if(step == 2) {
		set_op_step_t *s = (set_op_step_t *)in;
        int ksize = p->infile_number ? p->infiles[0].ksize : 0;
		if(!s->p->opt->text_out) {
			// write binary header
			if(p->header_init == 0) {
				p->header_init = 1;
				uint32_t t[3] = {ksize, p->out_cnt_number, 0};
				fwrite("TOD\2", 1, 4, p->opt->out);
				fwrite(t, 4, 3, p->opt->out);
			}
		
			if(s->out_buf_n) {
				fwrite(s->out_buf, s->out_buf_nbyte, s->out_buf_n, p->opt->out);
			}
		
		} else {
			kstring_t str = {0, 0, 0};
			ks_resize(&str, s->out_buf_n * 256);
			for(uint8_t *b = s->out_buf; b != s->out_buf + s->out_buf_n * s->out_buf_nbyte; b += s->out_buf_nbyte) {
				ks_resize(&str, str.l + ksize);
				kmer2seq(str.s + str.l, ksize, *(uint64_t *)b);
				str.l += ksize, str.s[str.l] = 0;
		
				for(int i = 0; i < p->out_cnt_number; ++i) {
					char c = (i == 0) ? '\t' : ',';
					kputc(c, &str);
					kputw(*(uint16_t *)(b + sizeof(uint64_t) + i * sizeof(uint16_t)), &str);
				}
		
				kputc('\n', &str);
			}
		
			if(str.l) {
				fwrite(str.s, sizeof(char), str.l, p->opt->out);
			}
			free(str.s);
		}
		if(s->out_buf) free(s->out_buf);
		set_op_step_destory(&s);

	}

    return 0;
}


void tod_opt_init(tod_opt_t *opt)
{
    memset(opt, 0, sizeof(tod_opt_t));
    opt->min_max_number.first = 1;
    opt->min_max_number.second = INT_MAX;
    opt->out = stdout;
    opt->outtype = OUT_NULL;
    kv_init(opt->vec_fps);
	kv_init(opt->min_max_counts);
	kv_resize(int_pair_t,opt->min_max_counts,32);
	
	opt->files_cap = 8;
	opt->absent_files = (uint8_t *)calloc(bit_array_size(opt->files_cap), 1);
    opt->present_files = (uint8_t *)calloc(bit_array_size(opt->files_cap), 1);
    opt->fill_cnt = 1;
}

static inline void tod_opt_update(tod_opt_t *opt)
{
    // auto adjust the min/max kmer count
    const int invalid_value=-1;
    if(kv_size(opt->min_max_counts) < kv_size(opt->vec_fps)) {
        for(int i = kv_size(opt->min_max_counts); i < kv_size(opt->vec_fps); ++i) {
			int_pair_t tmp={invalid_value,invalid_value};
			kv_push(int_pair_t, opt->min_max_counts, tmp);
        }
    } else {
        kv_size(opt->min_max_counts) = kv_size(opt->vec_fps);
    }
	
	int min_cnt=invalid_value,max_cnt=invalid_value;
	for(int i=0;i<kv_size(opt->min_max_counts);++i) {
		int_pair_t* min_max=&kv_A(opt->min_max_counts,i);
		if(min_max->first!=invalid_value) min_cnt=min_max->first;
		else if(min_cnt==invalid_value) {min_cnt=1;min_max->first=min_cnt;}
		else min_max->first=min_cnt;
		
		if(min_max->second!=invalid_value) max_cnt=min_max->second;
		else if(max_cnt==invalid_value)  {max_cnt=INT_MAX;min_max->second=max_cnt; } 
		else min_max->second=max_cnt;
	}

	int max_fileid = kv_size(opt->vec_fps);
	if(opt->files_cap <max_fileid) {
		opt->absent_files = realloc(opt->absent_files,bit_array_size(max_fileid));
		opt->present_files = realloc(opt->present_files,bit_array_size(max_fileid));
		memset(opt->absent_files+bit_array_size(opt->files_cap),
				0,
				bit_array_size(max_fileid)-bit_array_size(opt->files_cap));
		memset(opt->present_files+bit_array_size(opt->files_cap),
				0,
				bit_array_size(max_fileid)-bit_array_size(opt->files_cap));
		opt->files_cap = max_fileid;	
	} 
	
    switch(opt->set_type) {
        case SET_INTERSECT:
            opt->min_max_number.first = kv_size(opt->vec_fps);
            opt->min_max_number.second = kv_size(opt->vec_fps);
            for(int i = 0; i < kv_size(opt->vec_fps); ++i) {
                bit_array_set(opt->present_files, i);
            }
            break;
        case SET_UNION:
            opt->min_max_number.first = 1;
            opt->min_max_number.second = INT_MAX;
            break;
        case SET_DIFF:
            opt->min_max_number.first = 1;
            opt->min_max_number.second = 1;
            bit_array_set(opt->present_files, 0);
            for(int i = 1; i < kv_size(opt->vec_fps); ++i) {
                bit_array_set(opt->absent_files, i);
            }
            break;
        default:
            break;
    }

}


static inline void tod_opt_destroy(tod_opt_t *opt)
{
    if(opt->vec_fps.a) {
        kv_destroy(opt->vec_fps);
    }
    if(opt->min_max_counts.a) {
        kv_destroy(opt->min_max_counts);
    }

    if(opt->present_files) {
        free(opt->present_files);
    }

    if(opt->absent_files) {
        free(opt->absent_files);
    }

    fclose(opt->out);
}




static int kmer_set_op(int type, int argc, char *argv[])
{
    int parse_ok = 1;
    ketopt_t o = KETOPT_INIT;
    int c;
    tod_opt_t todopt;
    tod_opt_init(&todopt);
    todopt.set_type = type;

    static ko_longopt_t longopts[] = {
        { "fill", ko_required_argument, 300 }, //default value for file with only kmer column
        { "force", ko_required_argument, 301 }, //Kmer is forced to exist on these files
        { NULL, 0, 0 }
    };

    while((c = ketopt(&o, argc, argv, 1, "c:C:n:N:O:o:Th", longopts)) >= 0) {
        switch(c) {
            case 'c':
            case 'C': {
                int *fields, n;
                kstring_t s = {0};
                kputs(o.arg, &s);
                fields = ksplit(&s, ',', &n);
				
				while(kv_size(todopt.min_max_counts)<n) {
					int_pair_t tmp={1,INT_MAX};
					kv_push(int_pair_t,todopt.min_max_counts , tmp);
				}
                for(int i = 0; i < n; ++i) {
					int_pair_t* min_max = &kv_A(todopt.min_max_counts,i);
                    if(c == 'c') {
						min_max->first = (atoi(s.s + fields[i])<0)?0:(atoi(s.s + fields[i]));						
                    } else {
						min_max->second = (atoi(s.s + fields[i])<0)?INT_MAX:(atoi(s.s + fields[i]));	          
                    }
                }
                free(s.s);
                free(fields);
                break;
            }

            case 'n':
            case 'N': {
                if(c == 'n') {
					todopt.min_max_number.first = (atoi(o.arg)<0)?0:atoi(o.arg);
                } else {
                    todopt.min_max_number.second = (atoi(o.arg)<0)?INT_MAX:atoi(o.arg);
                }
                break;
            }

            case 'o': {
                if((todopt.out = (strcmp(o.arg, "-")) ? fopen(o.arg, "wb") : stdout) == 0) {
                    fprintf(stderr, "Invalid out file %s \n", o.arg);
                    parse_ok = 0;
                }
                break;
            }
            case 'O': {
                for(char *p = o.arg; *p != '\0'; ++p) {
                    *p = toupper(*p);
                }
                if(strcmp(o.arg, "COLLAPSE") == 0) {
                    todopt.outtype = OUT_COLLAPSE;
                }  else if(strcmp(o.arg, "SUM") == 0) {
                    todopt.outtype = OUT_SUM;
                }  else if(strcmp(o.arg, "ADD") == 0) {
                    todopt.outtype = OUT_ADD;
                } else if(strcmp(o.arg, "SUB") == 0) {
                    todopt.outtype = OUT_SUB;
                } else if(strcmp(o.arg, "MIN") == 0) {
                    todopt.outtype = OUT_MIN;
                } else if(strcmp(o.arg, "MAX") == 0) {
                    todopt.outtype = OUT_MAX;
                } else if(atoi(o.arg) > 0) {
                    todopt.outtype = OUT_VALID_LIMIT + atoi(o.arg);
                } else {
                    parse_ok = 0;
                }
                break;
            }
            case 'h': {
                parse_ok = 0;
                break;
            }
            case 'T': {
                todopt.text_out = 1;
                break;
            }
            case 300: {
                todopt.fill_cnt = atoi(o.arg);
                break;
            }
            case 301: {
                int *fields, n;
                kstring_t s = {0};
                kputs(o.arg, &s);
                fields = ksplit(&s, ',', &n);
                int max_fileid = 0;
                for(int i = 0; i < n; ++i) {
                    int fileid = s.s[fields[i]] == '^' ? (atoi(s.s + fields[i] + 1)) : (atoi(s.s + fields[i]));
                    if(fileid > max_fileid) {
                        max_fileid = fileid;
                    }
                }
                if(todopt.files_cap <max_fileid) {
					todopt.absent_files = realloc(todopt.absent_files,bit_array_size(max_fileid));
					todopt.present_files = realloc(todopt.present_files,bit_array_size(max_fileid));
					memset(todopt.absent_files+bit_array_size(todopt.files_cap),
							0,
							bit_array_size(max_fileid)-bit_array_size(todopt.files_cap));
					memset(todopt.present_files+bit_array_size(todopt.files_cap),
							0,
							bit_array_size(max_fileid)-bit_array_size(todopt.files_cap));
					todopt.files_cap = max_fileid;	
				} 
				
                for(int i = 0; i < n; ++i) {
                    int fileid;
                    if(s.s[fields[i]] == '^') {
                        fileid = atoi(s.s + fields[i] + 1)-1;
                        bit_array_set(todopt.absent_files, fileid);
                    } else {
						fileid = atoi(s.s + fields[i])-1;
                        bit_array_set(todopt.present_files, fileid);
                    }
                }
                free(s.s);
                free(fields);
                break;
            }
            default: {
                fprintf(stderr, "Invalid option %c \n", c);
                parse_ok = 0;
                break;
            }
        }
    }

    if((argc - o.ind < 1) || (parse_ok == 0)) {


        fprintf(stderr, "Usage: tod %s [options]  <kmer files ... > \n", str_set_op[type]);
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -c       INT     minimum kmer count present in file,separated by commas in multiple files [1]\n");
        fprintf(stderr, "  -C       INT     maximum kmer count present in file,separated by commas in multiple files [MAX_INT]\n");
        fprintf(stderr, "  -O       FILE    report counts of each kmer in a user-friendly manner,can be collapse,sum,sub,min,max,FILE-ID []\n");
        fprintf(stderr, "  -o       FILE    output file[-] \n");
        fprintf(stderr, "  -T               output in text format instead of binary format \n");
        if(type == SET_GROUP) {
            fprintf(stderr, "  -n       INT     minimum number of files containing the specified kmer count [%d]\n", todopt.min_max_number.first);
            fprintf(stderr, "  -N       INT     maximum number of files containing the specified kmer count [%d]\n", todopt.min_max_number.second);
            fprintf(stderr, "  --force  STR     a comma-separated list to indicate the presence/absence(^) of kmer in specified FILE-ID files[]\n");
        }
        fprintf(stderr, "  --fill   INT     default value for file with only kmer column[%d] \n", todopt.fill_cnt);

        tod_opt_destroy(&todopt);
        return 1;
    }


    do {
        // open files
        int i = o.ind;
        for(; i < argc; ++i) {
            gzFile fp = strcmp(argv[i], "-") ? gzopen(argv[i], "r") : gzdopen(fileno(stdin), "r");
            if(fp == 0) {
                fprintf(stderr, "[set] open file error %s \n", argv[i]);
                break;
            }

            kv_push(gzFile, todopt.vec_fps, fp);

        }
        if(i < argc) {
            break;
        }

        if((todopt.outtype > OUT_VALID_LIMIT) && todopt.outtype - OUT_VALID_LIMIT > kv_size(todopt.vec_fps))  {
            fprintf(stderr, "The outype set exceeds the number of files.\n");
            break;
        }

        // auto adjust the min/max kmer count
        tod_opt_update(&todopt);

        // print parameters
        fprintf(stderr, "[set] %s ,%lu files in total.\n", str_set_op[type], kv_size(todopt.vec_fps));
        if(type == SET_GROUP) {
            fprintf(stderr, "[set] select kmers present in [%d,%d] out of the %lu files\n", todopt.min_max_number.first, todopt.min_max_number.second, kv_size(todopt.vec_fps));
        }
        fprintf(stderr, "[set] kmer is considered to be present only when the count of kmers meets the specified range: \n");
        for(int i = o.ind; i < argc; ++i) {
            fprintf(stderr, "[set] ---------- [%d,%d] in %s \n", kv_A(todopt.min_max_counts, i - o.ind).first, kv_A(todopt.min_max_counts, i - o.ind).second,   argv[i]);
        }

        set_op_pl_t pl;
        pl.opt = &todopt;
        pl.header_init = 0;
        pl.infiles = malloc(sizeof(infile_t) * kv_size(todopt.vec_fps));
        pl.infile_number = kv_size(todopt.vec_fps);
        pl.total_extend_cnt_number = 0;

        // infile init
        for(int i = 0; i < kv_size(todopt.vec_fps); ++i) {
            infile_t *infile = pl.infiles + i;
            infile->ks = ks_init(kv_A(todopt.vec_fps, i));
            infile->fileid = i;

            // ksize and cnt_number
            file_kc_init(infile);

            infile->adjust_file_id = pl.total_extend_cnt_number;
            pl.total_extend_cnt_number += infile->cnt_number < 1 ? 1 : infile->cnt_number;
        }


        switch(pl.opt->outtype) {
            case OUT_COLLAPSE:
                pl.out_cnt_number = pl.total_extend_cnt_number;
                break;
            case OUT_NULL:
                pl.out_cnt_number = 0;
                break;
            case OUT_ADD:
            case OUT_SUB:
            case OUT_MAX:
            case OUT_MIN:
                pl.out_cnt_number = 1;
                break;
            default:
                pl.out_cnt_number = pl.infiles[pl.opt->outtype - OUT_VALID_LIMIT - 1].cnt_number;
                break;
        }

        pl.bufs = malloc(sizeof(buf_t) * pl.infile_number);
        for(int i = 0; i < pl.infile_number; ++i) {
            pl.bufs[i].valid_i = 0;
            pl.bufs[i].cap_i = 0;
            const int nbyte = sizeof(uint64_t) + sizeof(uint16_t) * pl.infiles[i].cnt_number;
            pl.bufs[i].buf =  calloc(N_BLOCK_END, nbyte);
            pl.bufs[i].nbyte = nbyte;
        }
        kt_pipeline(1, set_op_pipeline, &pl, 3);

        for(int i = 0; i < kv_size(todopt.vec_fps); ++i) {
            ks_destroy(pl.infiles[i].ks);
            file_kc_destroy(&pl.infiles[i]);
        }
        free(pl.infiles);
        for(int i = 0; i < pl.infile_number; ++i) {
            free(pl.bufs[i].buf);
        }
        free(pl.bufs);
    } while(0);

    // clean
    for(int i = 0; i < kv_size(todopt.vec_fps); ++i) {
        gzclose(kv_A(todopt.vec_fps, i));
    }
    tod_opt_destroy(&todopt);
    return 0;
}

int kmer_set_group(int argc, char *argv[])
{
    return kmer_set_op(SET_GROUP, argc, argv);
}


int kmer_set_intersect(int argc, char *argv[], int type)
{
    return kmer_set_op(SET_INTERSECT, argc, argv);
}

int kmer_set_union(int argc, char *argv[])
{
    return kmer_set_op(SET_UNION, argc, argv);
}

int kmer_set_diff(int argc, char *argv[])
{
    return kmer_set_op(SET_DIFF, argc, argv);
}


