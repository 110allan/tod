#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include "ketopt.h"
#include "kstring.h"
#include "kthread.h"
#include "kseq.h"
#include "kvec.h"
#include "ksort.h"

#include "misc.h"

typedef kvec_t(int) vec_int_t;
// output format, collapse,add(sum),sub,min,max and file_id
#define OUT_NULL 0
#define OUT_COLLAPSE 1
#define OUT_SUM 2
#define OUT_ADD OUT_SUM
#define OUT_SUB 3
#define OUT_MIN  4
#define OUT_MAX 5
#define OUT_VALID_LIMIT OUT_MAX


KSEQ_INIT(gzFile, gzread)

#include "khashl.h" // hash table
#define KC_BITS 16
#define KC_MAX ((1<<KC_BITS) - 1)

// ckmer: counting-contained kmer k{64-KC_BITS}c{KC_BITS}
#define ckmer_eq(a, b) ((a)>>KC_BITS == (b)>>KC_BITS)
#define ckmer_hash(a) ((a)>>KC_BITS)
KHASHL_SET_INIT(, ckmer_set_t, ckmer_set, uint64_t, ckmer_hash, ckmer_eq)
;
typedef kvec_t(uint64_t) vec_uint64_t;

typedef struct {
    ckmer_set_t *h;
    uint64_t kmer;
} ckmer_hashset_t;

#define ckmer_hashset_lt(a, b) ((a).kmer < (b).kmer)
KSORT_INIT(ckmer_hashset, ckmer_hashset_t, ckmer_hashset_lt)
;
typedef kvec_t(ckmer_hashset_t) vec_ckmer_hashset_t;
;
static vec_ckmer_hashset_t *vec_ckmer_set_init(int p)
{
    vec_ckmer_hashset_t *h = calloc(1, sizeof(vec_ckmer_hashset_t));
    kv_init(*h);
    kv_resize(ckmer_hashset_t, *h, 1 << p);
    kv_size(*h) = 1 << p;

    for(int i = 0; i < 1 << p; ++i) {
        kv_A(*h, i).kmer = i;
        kv_A(*h, i).h = ckmer_set_init();
    }
    return h;
}

static void vec_ckmer_set_release(vec_ckmer_hashset_t *h)
{
    for(int i = 0; i < kv_size(*h); ++i) {
        ckmer_set_destroy(kv_A(*h, i).h);
    }
    kv_destroy(*h);
}

static void add_kmer_buf(const char *seq, int len, int k, int p, int compressed, vec_uint64_t *buf) // extract k-mers in $seq to buffer vector
{
    int i, l;
    uint64_t x[2], mask = (1ULL << k * 2) - 1, shift = (k - 1) * 2;
    for(i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
        if(compressed && i && seq[i] == seq[i - 1]) {
            continue;
        }
        int c = seq_nt4_table[(uint8_t)seq[i]];
        if(c < 4) {  // not an "N" base
            x[0] = (x[0] << 2 | c) & mask;                  // forward strand
            x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
            if(++l >= k) {  // we find a k-mer
                // 64        2k  2k-p   0
                // ********* *** ********
                // [0:2k-p]: low bits for kmer in ckmer;   [2k-p,2k]: high bits for hash id
                const uint64_t kmer = (x[0] < x[1] ? x[0] : x[1])&mask;
                kv_push(uint64_t, buf[ kmer << (64 - k * 2) >> (64 - p)],  kmer);
            }
        } else {
            l = 0, x[0] = x[1] = 0;    // if there is an "N", restart
        }
    }
}

typedef struct { // global data structure for kt_pipeline()
    int k, block_len, n_thread, p, homo_compressed;
    kseq_t *ks;
    vec_ckmer_hashset_t *h;
	uint64_t sum_len;
	uint64_t kmer_cnt;
} pldat_t;

typedef struct { // data structure for each step in kt_pipeline()
    pldat_t *p;
    int n, m, sum_len, nk;
    int *len;
    char **seq;
    vec_uint64_t *buf;
} stepdat_t;

static void worker_for(void *data, long i, int tid) // callback for kt_for()
{
    stepdat_t *s = (stepdat_t *)data;
    vec_uint64_t *b = &s->buf[i];
    ckmer_set_t *h = kv_A(*s->p->h, i).h;
    int p = s->p->p;
    for(int j = 0; j < kv_size(*b); ++j) {
        khint_t k;
        int absent;
        k = ckmer_set_put(h, (kv_A(*b, j) & ((1ULL << (s->p->k * 2 - p)) - 1)) << KC_BITS, &absent);
        if((kh_key(h, k)&KC_MAX) < KC_MAX) {
            ++kh_key(h, k);
        }
    }
}

static void *worker_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
    pldat_t *p = (pldat_t *)data;
    if(step == 0) {  // step 1: read a block of sequences
        int ret;
        stepdat_t *s;
        CALLOC(s, 1);
        s->p = p;
        while((ret = kseq_read(p->ks)) >= 0) {
            int l = p->ks->seq.l;
            if(l < p->k) {
                continue;
            }
            if(s->n == s->m) {
                s->m = s->m < 16 ? 16 : s->m + (s->n >> 1);
                REALLOC(s->len, s->m);
                REALLOC(s->seq, s->m);
            }
            MALLOC(s->seq[s->n], l);
            memcpy(s->seq[s->n], p->ks->seq.s, l);
            s->len[s->n++] = l;
            s->sum_len += l;
            s->nk += l - p->k + 1;
            if(s->sum_len >= p->block_len) {
                break;
            }
        }
        if(s->sum_len == 0) {
            free(s);
        } else {
            return s;
        }
    } else if(step == 1) {  // step 2: extract k-mers
        stepdat_t *s = (stepdat_t *)in;
        int i, n = 1 << p->p, m;
        CALLOC(s->buf, n);
        m = (int)(s->nk * 1.2 / n) + 1;
        for(i = 0; i < n; ++i) {
            kv_resize(uint64_t, s->buf[i], m > 8 ? m : 8);
        }
        for(i = 0; i < s->n; ++i) {
            add_kmer_buf(s->seq[i], s->len[i], p->k, p->p, s->p->homo_compressed, s->buf);
            free(s->seq[i]);
        }

		p->sum_len += s->sum_len;
		for(i = 0; i < n; ++i) {
			p->kmer_cnt += kv_size(s->buf[i]);
        }
				
        free(s->seq);
        free(s->len);

        return s;
    } else if(step == 2) {  // step 3: insert k-mers to hash table
        stepdat_t *s = (stepdat_t *)in;
        int i, n = 1 << p->p;
        kt_for(p->n_thread, worker_for, s, n);
        for(i = 0; i < n; ++i) {
            kv_destroy(s->buf[i]);
        }
        free(s->buf);
        free(s);
    }
    return 0;
}

static vec_ckmer_hashset_t *count_file(const char *fn, int k, int p, int block_size, int n_thread, int homo_compressed, vec_ckmer_hashset_t *h,uint64_t* sumlen,uint64_t* kmercnt)
{
    pldat_t pl;
    gzFile fp;
    fp = fn && strcmp(fn, "-") ? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
    if(fp == 0) {
        return 0;
    }
    pl.ks = kseq_init(fp);
    pl.k = k;
    pl.n_thread = n_thread;
    pl.h = h ? h : vec_ckmer_set_init(p);
    pl.block_len = block_size;
    pl.p = p;
    pl.homo_compressed = homo_compressed;
	pl.sum_len = 0;
	pl.kmer_cnt = 0;
    kt_pipeline(3, worker_pipeline, &pl, 3);
	
	*sumlen = pl.sum_len;
	*kmercnt = pl.kmer_cnt;
    kseq_destroy(pl.ks);
    gzclose(fp);
    return pl.h;
}

/*hist*/
typedef struct {
    uint32_t c[KC_MAX + 1];
} buf_cnt_t;

typedef struct {
    const vec_ckmer_hashset_t *h;
    buf_cnt_t *cnt;/*every thread*/
} hist_aux_t;

static void worker_dump_hist(void *data, long i, int tid) // callback for kt_for()
{
    hist_aux_t *a = (hist_aux_t *)data;
    uint32_t *cnt = a->cnt[tid].c;
    ckmer_set_t *h = kv_A(*a->h, i).h;
    khint_t k;
    for(k = 0; k < kh_end(h); ++k)
        if(kh_exist(h, k)) {
            int c = kh_key(h, k) & KC_MAX;
            ++cnt[c < KC_MAX ? c : KC_MAX];
        }
}

static int dump_hist(const vec_ckmer_hashset_t *h, int n_thread, const char *fn)
{
    FILE *fp;
    hist_aux_t a;
    uint32_t cnt[KC_MAX];
    int i, j;
    a.h = h;
    if((fp = strcmp(fn, "-") ? fopen(fn, "wb") : stdout) == 0) {
        return -1;
    }
    CALLOC(a.cnt, n_thread);
    kt_for(n_thread, worker_dump_hist, &a, kv_size(*h));
    for(i = 0; i < KC_MAX; ++i) {
        cnt[i] = 0;
    }
    for(j = 0; j < n_thread; ++j)
        for(i = 0; i < KC_MAX; ++i) {
            cnt[i] += a.cnt[j].c[i];
        }
    free(a.cnt);
    for(i = 1; i < KC_MAX; ++i) {
        if(cnt[i]) {
            fprintf(fp, "%d\t%u\n", i, cnt[i]);
        }
    }
    fclose(fp);
    return 0;
}
/*tsih*/

#pragma pack (2)
typedef struct {
    uint64_t kmer;
    uint16_t cnt;
} kmc_pair_t;
#pragma pack ()

typedef kvec_t(kmc_pair_t) vec_kmc_pair_t;

#define kmc_pair_lt(a, b) ((a).kmer < (b).kmer)

KSORT_INIT(kmc_pair, kmc_pair_t, kmc_pair_lt)

;
typedef struct {
    const vec_ckmer_hashset_t *h;
    vec_kmc_pair_t *vec_ckmer;// every hashid
    int ksize, p;
} kmerdump_aux_t;

static void worker_dump_kmer(void *data, long i, int tid) // callback for kt_for()
{
    kmerdump_aux_t *a = (kmerdump_aux_t *)data;
    int ksize = a->ksize, p = a->p;
    ckmer_set_t *h = kv_A(*a->h, i).h;
    if(!kh_size(h)) {
        return;
    }

    kv_init(a->vec_ckmer[i]);
    kv_resize(kmc_pair_t, a->vec_ckmer[i], kh_size(h));

    for(khint_t k = 0; k < kh_end(h); ++k) {
        if(kh_exist(h, k)) {
            const uint64_t key = kh_key(h, k);
            const uint16_t cnt = key & KC_MAX;
            const uint64_t kmer = ((key >> KC_BITS) | (kv_A(*a->h, i).kmer << (ksize * 2 - p))) & ((1ULL << ksize * 2) - 1);
            const kmc_pair_t v = {kmer, cnt};
            kv_push(kmc_pair_t, a->vec_ckmer[i], v);
        }
    }

    ks_introsort(kmc_pair, kv_size(a->vec_ckmer[i]), a->vec_ckmer[i].a);
    ckmer_set_destroy(kv_A(*a->h, i).h);
    kv_A(*a->h, i).h = NULL;//to avoid double free in ckmer_set_destroy
}


static int dump_kmer(const vec_ckmer_hashset_t *hashsets, int n_thread, int ksize, int p, const char *fn, int text_out)
{
    FILE *fp;
    if((fp = (fn && strcmp(fn, "-")) ? fopen(fn, "wb") : stdout) == 0) {
        return -1;
    }

    ks_introsort(ckmer_hashset, kv_size(*hashsets), hashsets->a);
    kmerdump_aux_t a;
    a.h = hashsets;
    a.ksize = ksize;
    a.p = p;
    a.vec_ckmer = calloc(kv_size(*hashsets), sizeof(vec_kmc_pair_t));

    kt_for(n_thread, worker_dump_kmer, &a, kv_size(*hashsets));

    if(!text_out) {
        uint32_t nk = 0;
        uint32_t t[3];
        for(int i = 0; i < kv_size(*hashsets); ++i) {
            nk += kv_size(a.vec_ckmer[i]);
        }
        fwrite("TOD\2", 1, 4, fp);
        t[0] = ksize, t[1] = 1, t[2] = nk;
        fwrite(t, 4, 3, fp);
    }

    for(int i = 0; i < kv_size(*hashsets); ++i) {
        if(!text_out) {
            if(kv_size(a.vec_ckmer[i])) {
                fwrite(a.vec_ckmer[i].a, sizeof(kmc_pair_t), kv_size(a.vec_ckmer[i]), fp);
            }
        } else {
            char seqs[ksize + 1];
            seqs[ksize] = '\0';
            for(int idx = 0; idx < kv_size(a.vec_ckmer[i]); ++idx) {
                kmer2seq(seqs, ksize, a.vec_ckmer[i].a[idx].kmer);
                fprintf(fp, "%s\t%u \n", seqs, a.vec_ckmer[i].a[idx].cnt);
            }
        }
    }

    for(int i = 0; i < kv_size(*hashsets); ++i) {
        kv_destroy(a.vec_ckmer[i]);
    }
    free(a.vec_ckmer);

    fclose(fp);
    return 0;
}

int kmer_count(int argc, char *argv[])
{
    vec_ckmer_hashset_t *h = 0;
    int c, k = 21, p = KC_BITS, block_size = 100000000, n_thread = 4;
    int hompolymer_compressed = 0;
    int text_out = 0;
    char *outfile = "-";
    char *hist_file = NULL;

    ketopt_t o = KETOPT_INIT;
    int parse_ok = 1;
    while((c = ketopt(&o, argc, argv, 1, "k:p:Tt:H:o:c", 0)) >= 0) {
        switch(c) {
            case 'k':
                k = atoi(o.arg);
                break;
            case 'p':
                p = atoi(o.arg);
                break;
            case 'T':
                text_out = 1;
                break;
            case 't':
                n_thread = atoi(o.arg);
                break;
            case 'H':
                hist_file = o.arg;
                break;
            case 'o':
                outfile = o.arg;
                break;
            case 'c':
                hompolymer_compressed = 1;
                break;
            default:
                fprintf(stderr, "Invalid option %c \n", c);
                parse_ok = 0;
        }
    }

    if((argc - o.ind < 1) || (!parse_ok)) {
        fprintf(stderr, "Usage: tod count [options] <in.fa>\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -k INT	  k-mer size [%d]\n", k);
        fprintf(stderr, "  -c 		  count kmer in hompolymer compressed mode \n");
        fprintf(stderr, "  -t INT	  number of worker threads [%d]\n", n_thread);
        fprintf(stderr, "  -H FILE	  report a histogram of count for each kmer A\n");
        fprintf(stderr, "  -o FILE	  output file, report counts of each kmer [%s] \n", outfile);
        fprintf(stderr, "  -T 		  output in text format instead of binary format \n");
        return -1;
    }

    h = 0;
	uint64_t sum_len,kmer_cnt,total_sum_len=0,total_kmer_cnt=0;
    for(int i = o.ind; i < argc; ++i) {
        h = count_file(argv[i], k, p, block_size, n_thread, hompolymer_compressed, h,&sum_len,&kmer_cnt);
		total_sum_len+=sum_len,total_kmer_cnt+=kmer_cnt;
    }
	fprintf(stderr, "[stats] total sequence len = %lu ,kmer count = %lu \n",total_sum_len,total_kmer_cnt);
    if(hist_file) {
        dump_hist(h, n_thread, hist_file);
    }

    dump_kmer(h, n_thread, k, p, outfile, text_out);
    vec_ckmer_set_release(h);
    free(h);
    return 0;
}

int kmer_base_dump(FILE *in, FILE *out, int min_number, int max_number,  vec_int_t *min_counts, vec_int_t *max_counts, int outtype)
{


    uint32_t ksize, cnt_number, nk;
    //open files for read and write
    int file_error = 0;
    do {
        char header[4];
        if((fread(header, 1, 4, in) != 4) || (strncmp(header, "TOD\2", 4) != 0)) {
            fprintf(stderr, "ERROR: wrong file header.\n");
            file_error = 1;
            break;
        }

        uint32_t t[3];
        if(fread(t, 4, 3, in) != 3) {
            fprintf(stderr, "ERROR: wrong TOD input file.\n");
            file_error = 1;
            break;
        }
        ksize = t[0], cnt_number = t[1], nk = t[2];
        if(ksize > 32) {
            fprintf(stderr, ">32 kmer not allowed.\n");
            file_error = 1;
            break;
        }

        if(cnt_number == 0) {
            outtype = OUT_NULL;
        }
        if((outtype > OUT_VALID_LIMIT) && outtype - OUT_VALID_LIMIT > cnt_number)  {
            fprintf(stderr, "The outtype set exceeds the number of columns in the input file.\n");
            file_error = 1;
            break;
        }
        // auto adjust the min/max kmer count
        size_t extend_number = cnt_number;
        if(kv_size(*min_counts) < extend_number) {
            for(int i = kv_size(*min_counts); i < extend_number; ++i) {
                const int last_v = kv_A(*min_counts, kv_size(*min_counts) - 1);
                kv_push(int, *min_counts, last_v);

            }
        } else {
            kv_size(*min_counts) = extend_number;
        }

        if(kv_size(*max_counts) < extend_number) {
            for(int i = kv_size(*max_counts); i < extend_number; ++i) {
                const int last_v = kv_A(*max_counts, kv_size(*max_counts) - 1);
                kv_push(int, *max_counts, last_v);
            }
        } else {
            kv_size(*max_counts) = extend_number;
        }
    } while(0);

    if(file_error) {
        return -1;
    }

    //read from file
    uint64_t kmer;
    uint16_t result_cnt[cnt_number];

    kstring_t str = {0, 0, 0};
    ks_resize(&str, 1024);
    while(1) {

        if(fread(&kmer, 8, 1, in) != 1) {
            file_error = 1;
            break;
        }
        kmer2seq(str.s, ksize, kmer);
        str.l = ksize, str.s[str.l] = 0;

        if(cnt_number) {
            int pass_number = 0;
            for(int i = 0; i < cnt_number; ++i) {
                if(fread(&result_cnt[i], 2, 1, in) != 1) {
                    return;
                }

                if(!(result_cnt[i] < min_counts->a[i] || (max_counts->a[i] > 0 && result_cnt[i] > max_counts->a[i]))) {
                    pass_number++;
                }
            }

            if(pass_number < min_number || (max_number > 0 && pass_number > max_number)) {
                continue;
            }
        }
        //output
        kmer2seq(str.s, ksize, kmer);
        str.l = ksize, str.s[str.l] = 0;

        int cnt = result_cnt[0];
        switch(outtype) {
            case OUT_ADD://case OUT_SUM:
                kputc('\t', &str);
                for(int i = 1; i < cnt_number; ++i) {
                    cnt += result_cnt[i];
                }

                kputw(cnt, &str);
                break;
            case OUT_MIN:
                kputc('\t', &str);
                for(int i = 1; i < cnt_number; ++i)
                    if(cnt > result_cnt[i]) {
                        cnt = result_cnt[i];
                    }

                kputw(cnt, &str);
                break;
            case OUT_MAX:
                kputc('\t', &str);
                for(int i = 1; i < cnt_number; ++i)
                    if(cnt < result_cnt[i]) {
                        cnt = result_cnt[i];
                    }

                kputw(cnt, &str);
                break;
            case OUT_SUB:
                kputc('\t', &str);
                for(int i = 1; i < cnt_number; ++i) {
                    cnt -= result_cnt[i];
                }
                if(cnt < 0) {
                    cnt = 0;
                }
                kputw(cnt, &str);
                break;
            case OUT_COLLAPSE:
                kputc('\t', &str);
                kputw(cnt, &str);
                for(int i = 1; i < cnt_number; ++i) {
                    kputc(',', &str), kputw(result_cnt[i], &str);
                }
                break;
            case OUT_NULL:
                break;
            default:
                kputc('\t', &str);
                kputw(result_cnt[outtype - OUT_VALID_LIMIT - 1], &str);
                break;

        }
        kputc('\n', &str);
        fwrite(str.s, sizeof(char), str.l, out);
    }
    free(str.s);
    return file_error;
}


int kmer_dump(int argc, char *argv[])
{
    int parse_ok = 1;
    int c;

    int min_number = 1, max_number = -1;
    int outtype = OUT_COLLAPSE;
    vec_int_t min_counts;
    vec_int_t max_counts;
    kv_init(min_counts);
    kv_init(max_counts);

    FILE *in = NULL, *out = stdout;
    ketopt_t o = KETOPT_INIT;
    while((c = ketopt(&o, argc, argv, 1, "c:C:n:N:O:o:h", 0)) >= 0) {
        switch(c) {

            case 'c':
            case 'C': {
                int *fields, n;
                kstring_t s = {0};
                kputs(o.arg, &s);
                fields = ksplit(&s, ',', &n);
                for(int i = 0; i < n; ++i)
                    if(c == 'c') {
                        kv_push(int, min_counts, atoi(s.s + fields[i]));
                    } else {
                        kv_push(int, max_counts, atoi(s.s + fields[i]));
                    }
                free(s.s);
                free(fields);
                break;
            }

            case 'n':
            case 'N': {
                if(c == 'n') {
                    min_number = atoi(o.arg);
                } else {
                    max_number = atoi(o.arg);
                }
                break;
            }

            case 'o': {
                if((out = (strcmp(o.arg, "-")) ? fopen(o.arg, "wb") : stdout) == 0) {
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
                    outtype = OUT_COLLAPSE;
                }  else if(strcmp(o.arg, "SUM") == 0) {
                    outtype = OUT_SUM;
                }  else if(strcmp(o.arg, "ADD") == 0) {
                    outtype = OUT_ADD;
                } else if(strcmp(o.arg, "SUB") == 0) {
                    outtype = OUT_SUB;
                } else if(strcmp(o.arg, "MIN") == 0) {
                    outtype = OUT_MIN;
                } else if(strcmp(o.arg, "MAX") == 0) {
                    outtype = OUT_MAX;
                } else if(atoi(o.arg) > 0) {
                    outtype = OUT_VALID_LIMIT + atoi(o.arg);
                } else {
                    parse_ok = 0;
                }
                break;
            }
            case 'h': {
                parse_ok = 0;
                break;
            }

            default: {
                fprintf(stderr, "Invalid option %c \n", c);
                parse_ok = 0;
                break;
            }
        }
    }


    if(!kv_size(min_counts)) {
        kv_push(int, min_counts, 1);
    }

    if(!kv_size(max_counts)) {
        kv_push(int, max_counts, -1);
    }

    if(!(parse_ok && (argc - o.ind <= 1))) {

        fprintf(stderr, "Usage: tod dump [options] <in.kc>\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -c INT     minimum kmer count present in file,separated by commas in multiple files [%d]\n", kv_A(min_counts, 0));
        fprintf(stderr, "  -C INT     maximum kmer count present in file,separated by commas in multiple files [%d]\n", kv_A(max_counts, 0));
        fprintf(stderr, "  -n INT     minimum number of files containing the specified kmer count [%d]\n", min_number);
        fprintf(stderr, "  -N INT     maximum number of files containing the specified kmer count [%d]\n", max_number);
        fprintf(stderr, "  -O STR     report counts of each kmer in a user-friendly manner,can be collapse,sum,sub,min,max,FILE-ID [collapse]\n");
        fprintf(stderr, "  -o FILE    output file [-] \n");

        kv_destroy(min_counts);
        kv_destroy(max_counts);
        return -1;
    }


    if((argc - o.ind == 0) || (strcmp(argv[o.ind], "-") == 0)) {
        in = stdin;
    } else {
        in = fopen(argv[o.ind], "rb") ;
    }


    if(kmer_base_dump(in, out, min_number, max_number, &min_counts, &max_counts, outtype) < 0) {
        fprintf(stderr, "ERROR: read in file.\n");
    }

    kv_destroy(min_counts);
    kv_destroy(max_counts);

    if(in) {
        fclose(in);
    }
    if(out) {
        fclose(out);
    }
    return 0;
}

