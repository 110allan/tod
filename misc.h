#ifndef MISC_H
#define MISC_H


#define bit_array_size(n) ((n)/8+1)
#define bit_array_set(a,i)   ((a)[(i)/8] |=   1 << ((i)%8))
#define bit_array_clear(a,i) ((a)[(i)/8] &= ~(1 << ((i)%8)))
#define bit_array_test(a,i)  ((a)[(i)/8] &   (1 << ((i)%8)))

#define CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

// extern unsigned char seq_nt4_table[256];
static unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};



static inline uint64_t seq2kmer(char *seq, int k)
{
    uint64_t kmer = 0, krev = 0;
    const uint64_t shift = (k - 1) * 2;
    for(int i = 0; i < k; i++) {
        int c = seq_nt4_table[(uint8_t)seq[i]] & 0x3 ; //^ACGT->A
        kmer = (kmer << 2) | c;
        krev = krev >> 2 | (uint64_t)(~c) << shift;
    }
    return kmer < krev ? kmer : krev;
}

static inline void kmer2seq(char * seq, int ksize, uint64_t kmer)
{
    for(int i = 0; i < ksize; i++) {
        seq[i] = "ACGTN"[(kmer >> ((ksize - 1 - i) * 2)) & 0x03];
    }
}

static inline uint64_t hash(uint64_t key, uint64_t mask) // invertible integer hash function
{
	//return key&mask;
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

// The inversion of hash_64(). Modified from <https://naml.us/blog/tag/invertible>
static inline uint64_t unhash(uint64_t key, uint64_t mask)
{
	//return key&mask;
	uint64_t tmp;

	// Invert key = key + (key << 31)
	tmp = (key - (key << 31));
	key = (key - (tmp << 31)) & mask;

	// Invert key = key ^ (key >> 28)
	tmp = key ^ key >> 28;
	key = key ^ tmp >> 28;

	// Invert key *= 21
	key = (key * 14933078535860113213ull) & mask;

	// Invert key = key ^ (key >> 14)
	tmp = key ^ key >> 14;
	tmp = key ^ tmp >> 14;
	tmp = key ^ tmp >> 14;
	key = key ^ tmp >> 14;

	// Invert key *= 265
	key = (key * 15244667743933553977ull) & mask;

	// Invert key = key ^ (key >> 24)
	tmp = key ^ key >> 24;
	key = key ^ tmp >> 24;

	// Invert key = (~key) + (key << 21)
	tmp = ~key;
	tmp = ~(key - (tmp << 21));
	tmp = ~(key - (tmp << 21));
	key = ~(key - (tmp << 21)) & mask;

	return key;
}


static inline  uint64_t hash2( uint64_t x )  {
    x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
    x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
    x = ( ( x >> 32 ) ^ x );
    return x;
}

static inline  uint64_t unhash2 (uint64_t x ) {
    x = ( ( x >> 32 ) ^ x ) * 0xCFEE444D8B59A89B;
    x = ( ( x >> 32 ) ^ x ) * 0xCFEE444D8B59A89B;
    x = ( ( x >> 32 ) ^ x );
    return x;
}



#endif
