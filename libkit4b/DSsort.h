#pragma once

const int cMaxComparedLCP = 300000;					// maximum length LCP compared, if longer then treat LCP+1 as mismatch

const int cDSSFREESIZE=5000;
const int cDSSMax_thresh=30;



class CDSsort
{
/* ------- node of blind trie -------- */ 
typedef struct nodex {
  int32_t skip;
  uint8_t key;  
  struct nodex  *down;      // first child
  struct nodex *right;      // next brother
	} node;


		void *m_pFreeArray[cDSSFREESIZE];
		node *m_pBufn;
		int32_t m_bufn_num;
		int32_t m_free_num;
		int32_t m_Aux_written;
		int32_t *m_pAux;
		node **m_ppStack;
		int m_Stack_size;
		int32_t m_Cmp_done;
		int32_t  m_Text_size;           // size of input string 
		uint8_t  *m_pText;        // input string+ overshoot
		int32_t  *m_pSuffixArray;       // suffix array
		uint8_t  *m_pEndOfText;   // m_pText+m_Text_size
		int32_t  *m_pAnchorRank;        // rank (in the sorted suffixes of the  
									// anchor points (-1 if rank is unknown))
		uint16_t  *m_pAnchorOffset;     // offset (wrt to the anchor) of the suffix
									// whose rank is in m_pAnchorRank. 
		int32_t m_NumAnchorPts;           // number of anchor points
		int32_t m_FtabArray[65537];   
		int32_t m_RunningOrderArray[256];

		uint8_t m_BucketRankedArray[65536];

		int32_t m_AnchorDist;                // distance between anchors
		int32_t m_DsVerbose;                // how verbose it the algorithm?
		int32_t m_DsWordSize;              // # of bytes in word in mkqs
		int32_t m_MkQsThresh;               // recursion limit for mk quicksort:
		int32_t m_MaxPseudoAnchorOffset; // maximum offset considered when 
										// searching a pseudo anchor
		int32_t m_B2gRratio;                // maximum ratio bucket_size/group_size
										// accepted for pseudo anchor_sorting
		int32_t m_UpdateAnchorRanks;      // if!=0 update anchor ranks when determining
										// rank for pseudo-sorting
		int32_t m_BlindSortRatio;         // blind sort is used for groups of size 
										// <= m_Text_size/m_BlindSortRatio

		int32_t m_Calls2HelpedSort;     
		int32_t m_Calls2AnchorSortForw;     
		int32_t m_Calls2AnchorSortBbackw;    
		int32_t m_Calls2PseudoAnchorSortForw;      
		int32_t m_Calls2DeepSort;     
		int32_t m_ShallowLimit;						// Max depth for shallow sorting
		uint8_t *m_pShallowTextLimit;         // m_pText+m_ShallowLimit
		int32_t m_CmpLeft;
		int32_t m_LCPAuxArray[1+cDSSMax_thresh];
		int32_t *m_pLCP; 



	void blind_ssort(int32_t *a, int32_t n, int32_t depth);
	node *find_companion(node *head, uint8_t *s);
	node *get_leaf(node *head);
	inline node *new_node__blind_ssort(void);
	void insert_suffix(node *h, int32_t suf, int n, uint8_t mmchar);
	void traverse_trie(node *h);
	inline int32_t get_lcp_unrolled(uint8_t *b1, uint8_t *b2, int32_t cmp_limit);
	int32_t compare_suffixes(int32_t suf1, int32_t suf2, int32_t depth);
	static int neg_integer_cmp(const void *a, const void *b);
	static int integer_cmp(const void *a, const void *b);
	inline int32_t cmp_unrolled_lcp(uint8_t *b1, uint8_t *b2);
	void qs_unrolled_lcp(int32_t *a, int n, int depth, int blind_limit);
	void deep_sort(int32_t *a, int32_t n, int32_t depth);
	void calc_running_order(void);
	void set_global_variables(void);
	int compute_overshoot(void);
	void pretty_putchar(int c);
	int scmp3(uint8_t *p, uint8_t *q, int *l, int maxl);
	void helped_sort(int32_t *a, int n, int depth);
	void pseudo_or_deep_sort(int32_t *a, int32_t n, int32_t depth);
	void pseudo_anchor_sort(int32_t *a,int32_t n,int32_t pseudo_anchor_pos, int32_t offset);
	void general_anchor_sort(int32_t *a, int32_t n,int32_t anchor_pos, int32_t anchor_rank, int32_t offset);
	int32_t get_rank(int32_t pos);
	int32_t get_rank_update_anchors(int32_t pos);
	void update_anchors(int32_t *a, int32_t n);
	int32_t split_group(int32_t *a, int n, int depth,int offset,int32_t pivot,int *first);
	void shallow_sort(int32_t *a, int n);
	inline void vecswap2(int32_t *a, int32_t *b, int n);
	inline int32_t *med3func(int32_t *a, int32_t *b, int32_t *c, uint8_t *text_depth);
	void shallow_mkq(int32_t *a, int n, uint8_t *text_depth);
	void shallow_mkq16(int32_t *a, int n, uint8_t *text_depth);
	void shallow_mkq32(int32_t *a, int n, uint8_t *text_depth);
	int check_global_variables(void);
	inline int32_t cmp_unrolled_shallow_lcp(uint8_t *b1, uint8_t *b2);
	void shallow_inssort_lcp(int32_t *a, int32_t n, uint8_t *text_depth);
	void free_node_mem(void);
public:
	CDSsort(void);
	~CDSsort(void);
	void ds_ssort(uint8_t *pToSort,	// string to sort suffixes on (NOTE: must have additional allocated memory past end for holding overshoot)
				  int32_t *pSuffixArray,	// sort into this prealloc'd suffix array
				  int32_t ToSortLen);		// string length (NOTE: excludes additional memory at end of pToSort allocated for overshoot)
	int init_ds_ssort(int adist = 500, int bs_ratio = 2000); // returns required overshoot memory size
};
