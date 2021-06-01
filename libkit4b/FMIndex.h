#pragma once

const int cMinL2BlockSize  = 1;						// minimum accepted level 2 block size in 1K increments
const int cDfltL2BlockSize = 1;						// default level 2 block size in 1K increments
const int cMaxL2BlockSize  = 4;						// maximum accepted level 2 block size in 1K increments
const int cMinL1BlockSize  = cMinL2BlockSize * 4;	// minimum accepted level 1 block size (in 1K increments)
const int cDfltL1BlockSize = cMinL1BlockSize * 4;	// default level 1 block size (in 1K increments)
const int cMaxL1BlockSize  = 0x7fff;				// maximum accepted level 1 block size (in 1K increments)
const double cDfltMarkerFreq = 0.02;				// marker frequency (0.0 <= Freq <= 1.0)

const int BZ_N_GROUPS = 6;
const int BZ_MAX_CODE_LEN = 23;
const int BZ_MAX_ALPHA_SIZE = 258;
const int ALPHASIZE = 256;

		
#define EXT ".fmi"



class CFMIndex
{
	typedef struct TAG_sBucket_lev1 {
	  uint32_t *occ;             /* occ chars of compact alph in prev. superbuc */
	  uint16_t alpha_size;
	  uint8_t *bool_char_map;   /* boolean map of chars occurring in this superbucket */
	} bucket_lev1;
	 
	typedef struct TAG_sFm_index {
	  uint8_t *text;					/* input text */
	  uint8_t *oldtext;				/* A copy of the input text */
	  uint8_t *compress;				/* compress text */
	  uint8_t *bwt;					/* BWT of input text */
	  uint32_t *lf;					/* lf-mapping or suffix-array*/
	  uint32_t compress_size;			/* size of compressed file */
	  uint32_t text_size;				/* size of text */
	  bucket_lev1 *buclist_lev1; 	/* array of num_bucs buckets */
		
	  /* Info readed/writed on index prologue */
	  uint32_t bucket_size_lev1;		/* size of level 1 buckets */	
	  uint32_t bucket_size_lev2;		/* size of level 2 buckets */
	  uint32_t num_bucs_lev1;			/* number buckets lev 1 */
	  uint32_t num_bucs_lev2;			/* number buckets lev 2 */
	  uint32_t bwt_eof_pos;			/* position of EOF within BWT */
	  uint16_t alpha_size;				/* actual size of alphabet in input text */
	  uint16_t type_compression;		/* buckets lev 1 type of compression */
	  double freq;					/* frequency of marked chars */
	  uint16_t owner;					/* == 0 frees the text and alloc a new with overshoot */	
	  uint16_t compress_owner;         /* == 1 frees the compress */
	  uint16_t smalltext;				/* If == 1 stores plain text without compression */

	  /* Starting position (in byte) of each group of info */
	  uint32_t start_prologue_info_sb;	/* byte inizio info sui superbuckets */
	  uint32_t start_prologue_info_b;	/* byte inizio posizioni buckets */
	  uint8_t *start_prologue_occ;	/* byte inizio posizioni marcate */
	  uint32_t start_positions;        /* byte inizio posizioni marcate per build */
	  uint32_t *start_lev2;			/* starting position of each buckets in compr file */

	  /* Chars remap info */
	  uint16_t bool_char_map[ALPHASIZE];/* is 1 if char i appears on text */
	  uint8_t char_map[ALPHASIZE];	 /* cm[i]=j say that chars i is remapped on j */
	  uint8_t inv_char_map[ALPHASIZE]; /* icm[i] = j say that j was remapped on i */

	  /* Multiple locate. Chars substitution */
	  uint8_t specialchar;			/* carattere speciale che indica marcamento */
	  uint8_t subchar;				/* carattere sostituito dal carattere speciale */

	  // Da spostare in una struct working_space
	  /* Running temp info of actual superbucket and bucket */
	  uint32_t pfx_char_occ[ALPHASIZE];	/* i stores # of occ of chars 0.. i-1 in the text */
	  uint16_t	bool_map_sb[ALPHASIZE]; 	/* info alphabet for superbucket to be read */
	  uint8_t inv_map_sb[ALPHASIZE]; 	/* inverse map for the current superbucket */
	  uint16_t alpha_size_sb;	  					/* current superbucket alphasize */
	  
	  uint8_t	bool_map_b[ALPHASIZE];		/* info alphabet for the bucket to be read */
	  uint8_t inv_map_b[ALPHASIZE];		/* inverse map for the current superbucket */
	  uint16_t alpha_size_b;     					/* actual size of alphabet in bucket */
	    
	  uint8_t mtf[ALPHASIZE];  		/* stores MTF-picture of bucket to be decompressed */
	  uint8_t *mtf_seq;				/* store bucket decompressed */
	  uint32_t occ_bucket[ALPHASIZE];  /* number chars occurences in the actual bucket needed by Mtf2 */
	  uint16_t int_dec_bits; 			/* log2(log2(text_size)) */
	  uint16_t log2textsize; 			/* int_log2(s.text_size-1) */
	  uint16_t var_byte_rappr;   		/* variable byte-length repr. (log2textsize+7)/8;*/
	  uint32_t num_marked_rows;		/* number of marked rows */
	  uint32_t bwt_occ[ALPHASIZE];     /* entry i stores # of occ of chars 0..i-1 */
	  uint32_t char_occ[ALPHASIZE];	/* (useful for mtf2) entry i stores # of occ of char i == bwt_occ[i]-bwt_occ[i-1] */
	  uint16_t skip;			/* 0 no marked pos, 1 all pos are marked, 2 only a char is marked*/
	  uint16_t sb_bitmap_size;			/* size in bytes of the bitmap of superbuckets */	
	  uint32_t occcharinf; 			/* Number occs of the chars  < index->specialchar */
	  uint32_t occcharsup; 			/* Number occs of the chars  <= index->specialchar */
	  
	  /* Needed by fm_build */
	  uint32_t *loc_occ;				/* Positions of marked rows */
	  
	} fm_index;





	/* Report rows from mulri_count */
	typedef struct TYPE_sMulti_count {
		uint32_t first_row; 			/* riga inizio occorrenza */
		uint32_t elements;   			/* numero occorrenze */	
		} multi_count;


	int __Fm_Verbose;

	fm_index *m_pIndex;							// contains complete context for compressed index

	uint32_t * __Num_Bytes;				/* numero byte letti/scritti */
	uint8_t * __MemAddress;				/* indirizzo della memoria dove scrivere */
	int __Bit_buffer_size;						/* number of unread/unwritten bits in Bit_buffer */
	uint32_t __pos_read;
	uint32_t __Bit_buffer;

	uint8_t *pmtf_start;
	int gmtflen;

	int allocated;
	int used;	/* var usate dalla count multipla */
	multi_count *lista;


	int count_row_mu (uint8_t * pattern, uint32_t len, uint32_t sp,uint32_t ep);
	inline void get_pos (uint32_t first_row, uint32_t element, uint16_t step, uint32_t * pos);
	int multi_locate (uint32_t sp, uint32_t element, uint32_t * positions);
	void preBmBc(uint8_t *x, int m, int bmBc[]);
	void suffixes(uint8_t *x, int m, int *suff);
	void preBmGs(uint8_t *x, int m, int bmGs[]);
	int fm_boyermoore(uint8_t * pattern, uint32_t length, uint32_t ** occ, uint32_t * numocc);
	int open_file(char * filename, uint8_t ** file, uint32_t * size);
	int fm_read_basic_prologue (void);
	int fm_multi_count (uint8_t * pattern, uint32_t len,multi_count ** list);

	int occ_all (uint32_t sp, uint32_t ep, uint32_t * occsp, uint32_t * occep,uint8_t * char_in);
	int get_info_sb (uint32_t pos, uint32_t * occ);
	int get_info_b (uint8_t ch, uint32_t pos, uint32_t *occ, int flag);
	uint8_t get_b_multihuf(uint32_t k, uint32_t * occ,int is_odd);
	inline void unmtf_unmap (uint8_t * mtf_seq, int len_mtf);
	int compress_bucket(uint8_t *in, uint32_t len, uint16_t alphasize);
	void mtf_string(uint8_t *in, uint8_t *out, uint32_t len, uint16_t mtflen);
	int fm_bwt_compress(void);
	void fm_unmtf(uint8_t *in, uint8_t *out, int length);
	int fm_bwt_uncompress(void);

	int fm_multihuf_compr (uint8_t * in, int len, int alpha_size);
	inline int decode_unary(void);
	int fm_multihuf_decompr (uint8_t * dest, int alpha_size, int limit);
	int fm_uncompress_bucket_multihuf (uint8_t * dest, int len, int alpha_size);

	uint8_t huf_len[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];	// coding and decoding
	int huf_code[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];	// coding
	int rfreq[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];	// coding
	int mtf_freq[BZ_MAX_ALPHA_SIZE];	// coding

	uint8_t huf_minLens[BZ_N_GROUPS];	// decoding
	int huf_limit[BZ_N_GROUPS][BZ_MAX_CODE_LEN];	// decoding
	int huf_base[BZ_N_GROUPS][BZ_MAX_CODE_LEN];	// decoding
	int huf_perm[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];	// decoding

	void fm_init_bit_writer (uint8_t * mem, uint32_t * pos_mem);
	void fm_bit_write (int n, uint32_t vv);
	int fm_bit_read (int n);
	void fm_bit_write24(int bits, uint32_t num);
	void fm_uint_write (uint32_t uu);
	uint32_t fm_uint_read (void);
	void fm_init_bit_reader (uint8_t * mem);
	void fm_bit_flush (void);
	uint32_t fm_integer_decode (unsigned short int headbits);

	int parse_options(char *optionz);
	int build_index(uint8_t *text, uint32_t length, char *build_options);
	int build_sa(void);
	void count_occ(void);
	int build_bwt(void);
	int compute_locations(void);
	int compute_info_superbuckets(void);
	int compute_info_buckets(void);
	void write_prologue(void);
	int compress_superbucket(uint32_t num);
	void write_susp_infos(void);
	void write_locations(void);
	int select_subchar(void);
	void dealloc(void);
	void dealloc_bucketinfo(void);
	int errore(int error);
	int save_index(char *filename);
	int fm_read_file(char *filename,		// file to read from 
				 uint8_t **textt,			// returned buffer allocated using malloc() containing contents of file + space for suffix sorting 
				 uint32_t *length);

	const char *error_index(int e);
	int int_log2(int u);
	int int_pow2(int u);

	uint32_t go_back(uint32_t row, uint32_t len, uint8_t *dest);
	uint32_t go_forw(uint32_t row, uint32_t len, uint8_t *dest);
	uint32_t fl_map(uint32_t row, uint8_t ch);
	uint8_t get_firstcolumn_char(uint32_t row);

	int fm_snippet(uint32_t row, uint32_t plen, uint32_t clen, uint8_t *dest, 
			   uint32_t *snippet_length);
	int read_prologue(void);
	int uncompress_data(void);
	int uncompress_superbucket(uint32_t numsb, uint8_t *out);
	int fm_compute_lf(void);
	int fm_invert_bwt(void);
	void free_unbuild_mem(void);
	int fm_unbuild(uint8_t ** text, uint32_t *length);

	void hbMakeCodeLengths(uint8_t *len,uint32_t *freq,uint32_t alphaSize,uint32_t maxLen);
	void hbAssignCodes(uint32_t *code,uint8_t *length,uint32_t minLen,uint32_t maxLen,uint32_t alphaSize);
	void hbCreateDecodeTables(uint32_t *limit,uint32_t *base,uint32_t *perm,uint8_t *length,uint32_t minLen,uint32_t maxLen,uint32_t alphaSize);
	

public:
	CFMIndex(void);
	~CFMIndex(void);
	int CreateIndex(char *pszInFile,		// create from contents of this file
					char *pszOutFile,		// write index into this file
					unsigned int *pText_len = NULL, // returned file content length before index created	
					unsigned int *pIndex_len = NULL, // created index length
					int bsl1 = cDfltL1BlockSize,	// level 1 block size (in 1K increments) will be forced to be multiple of bsl2
					int bsl2 = cDfltL2BlockSize,	// level 2 block size (in byte increments) will be forced to be multiple of 256
					double Freq = cDfltMarkerFreq);	// marker frequency (0.0 <= Freq <= 1.0)

	int load_index (char * filename);
	int extract(uint32_t from, uint32_t to, uint8_t **dest,uint32_t *snippet_length);
	int free_index (void);
	int display(uint8_t *pattern, uint32_t length, uint32_t nums, uint32_t *numocc, 
			uint8_t **snippet_text, uint32_t **snippet_len);
	int locate (uint8_t * pattern, uint32_t length, uint32_t ** occ,uint32_t * numocc);
	int count (uint8_t * pattern, uint32_t length, uint32_t * numocc);
	int get_length (uint32_t * length);

	int fm_build_config(double freq, uint32_t bsl1, uint32_t bsl2, uint16_t owner);
	int fm_build(uint8_t *text, uint32_t length);
	int index_size(uint32_t *size);


	int load_index_mem(uint8_t *compress, uint32_t size);		// loads compressed content from user supplied memory
	int save_index_mem(uint8_t *compress);					// saves compressed content into user supplied memory	

	char *GetErrText(int error);

};

