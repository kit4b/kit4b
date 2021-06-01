/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#pragma once

const uint32_t cMaxEdges = 250;						// limit any vertex to have at most this many inbound or outbound edges

													// Note: normalised overlap edge scores are independently calculated using 1 match, -1 mismatch, -3 gap open and -1 gap extension
const int cMinMin1kScore = 900;						 // when processing the multialignment overlaps for graph edges then only accept edge if normalised overlap score at least this per 1Kbp of overlap
const int cDfltMin1kScore = 980;					 // when processing the multialignment overlaps for graph edges then only accept edge if normalised overlap score at least this per 1Kbp of overlap
const int cMaxMin1kScore = 1000;					 // when processing the multialignment overlaps for graph edges then perfect match required

const uint32_t cInitialAllocVertices =500000;		// initially allocate for this many vertices, will be realloc'd if required
const double cReallocVertices   =   0.3;	    // then, as may be required, realloc in increments of this proportion of existing vertices
const uint32_t cInitialAllocEdges =  (cInitialAllocVertices * 10);	// initially allocate for this many outgoing edges
const double cReallocEdges      =   0.3;		// then, as may be required, realloc in increments of this proportion of existing outgoing edges
const uint32_t cInitalComponentsAlloc   = 50000;	// not expecting too many assemblies with more than 50K components but will realloc if required
const double cReallocComponents  = 0.3;			// realloc as may be required in this proportion of existing components

const uint64_t cMaxGraphEdges = 0x07fffffff;		// allowing at most this many edges in graph

const uint32_t cTransitStackAlloc   =    1000000;		// initially allocate for transition stack to hold this many entries
const double cReallocTransitStack =    0.3;		// realloc by this proportion of existing stack entries

const uint32_t cInitialAllocTraceBacks=100000;	// initially allocate for this many vertices, will be realloc'd if required
const double cReallocTraceBacks   =   0.3;	    // then, as may be required, realloc in increments of this proportion of existing tracebacks

const uint32_t cMaxDiscRemaps = 1000;					// remap disconnected graph identifiers list limit

typedef enum TAG_eVerticesSortOrder {
	eVSOUnsorted = 0,	//  unsorted or sort order indeterminate
	eVSOVertexID,		// sorted by vertex identifier ascending
	eVSOSeqID,			// sorted by sequence identifier ascending
	eVSOVertexComponentID // sorted by ComponentID ascending
	} eVerticesSortOrder;

typedef enum TAG_eOverlapClass {
	eOLCOverlapping = 0,	// probe classified as overlapping target, either 5' or 3'
	eOLCcontains,			// probe completely contains the target
	eOLCcontained,			// probe is completely contained within the target
	eOLCartefact			// probe contains a subsequence of target, classifying as an artefact overlap and not further processed
} eOverlapClass;

#pragma pack(1)

// graph consists of vertices (representing sequences) and connecting edges (overlaying sequences) between adjacent vertices
// edges are of two types - forward edges represent overlaying sequences (vertexA is overlaying vertexB), and overlaid (vertexA is overlaid by vertexB)
// it is also important to note that the graph is likely to contain many - perhaps millions - of disconnected components 

// an instance of a Vertex which is used to represent a read sequence
typedef struct TAG_sGraphVertex {
	    tVertID VertexID;		// monotonically ascending (1..V) identifier uniquely identifies this vertex
		tEdgeID OutEdgeID;		// m_pGraphOutEdges[OutEdgeID-1] at which outgoing edges from this vertex start, 0 if none
		tEdgeID InEdgeID;		// m_pGraphInEdges[InEdgeID-1] identifying incoming edges to this vertex start, 0 if none
		tSeqID SeqID;			// identifies the sequence being represented by this vertex
		uint32_t SeqLen;			//  sequence is this length
		uint32_t RecurseDepth;		// recurse depth at which this vertex is processed - used to determine circular reference
		tComponentID ComponentID; // vertex is a member of this disconnected component
		uint8_t DegreeOut;		// number of outgoing edges from this vertex to any adjacent vertices, clamped to max cMaxEdges (currently 250)
		uint8_t DegreeIn;			// number of incoming edges from any adjacent vertices, clamped to max cMaxEdges (currently 250)
		uint8_t flgEmitted:1;		// sequence has been emitted as part of a fragment;	
		uint8_t flgRmvEdges:1;	// sequence may contain SMRTbell hairpin as this read aligns at least twice to the same other read
		uint8_t flgPathAccepted:1;	// this vertex has been accepted as part of a highest scoring path
		uint8_t flgPathScored:1;		// set if highest scoring path already processed
		uint8_t flgPathTerm:1;		// vertex classified as path terminating
		uint64_t PathScore;			// highest score for any path originating from this vertex
		tEdgeID PathScoreEdgeID;	// highest scoring path starts with this outgoing edge
	} tsGraphVertex;


// an instance of an outgoing edge, from 'FromVertexID' to 'ToVertexID' 
// normally will be sorted in FromVertexID.ToVertexID ascending order
typedef struct TAG_sGraphOutEdge {
	tVertID FromVertexID;	// edge is outgoing from this vertex 
	tVertID ToVertexID;		// edge is incoming to this Vertex
	uint32_t FromSeqLen;		// 'From' vertex sequence length is this length - more efficent to duplicate copy of length with edges as it saves a lookup in graph vertices when length is needed
	uint32_t ToSeqLen;		// 'To' vertex sequence length is this length - more efficent to duplicate copy of length with edges as it saves a lookup in graph vertices when length is needed
	uint32_t FromSeq5Ofs;		// overlap of FromVertexID onto ToVertexID starts at this base relative to FromVertexID 5' start
	uint32_t FromSeq3Ofs;		// overlap of FromVertexID onto ToVertexID ends at this base relative to FromVertexID 5' start
	uint32_t ToSeq5Ofs;		// overlap onto ToVertexID from FromVertexID starts at this base relative to ToVertexID 5' start
	uint32_t ToSeq3Ofs;		// overlap onto ToVertexID from FromVertexID ends at this base relative to ToVertexID 5' start
	uint32_t Score;			// score associated with this overlap, higher scores represent higher confidence in the overlap
	uint32_t ScoreAlignLen;	// scored for this alignment length: (1 + ProbeAlignLen + TargAlignLen) / 2
	uint16_t flgFromAntisense:1;// set if From was revcpl as antisense probe 
	uint16_t flgToAntisense:1;// set if To was revcpl as antisense target 
	uint16_t flgRemove:1;		// set if this edge marked for removal
	uint16_t flgTravFwd:1;	// set after this edge has been traversed from FromVertexID to ToVertexID
	uint16_t flgTravRev:1;	// set after this edge has been traversed from ToVertexID to FromVertexID
	uint16_t flgInfBackEdge:1;// set if this is an inferred back edge, if 0 then edge was result of actual alignment
	} tsGraphOutEdge;

typedef struct TAG_sOverlappedSeq {
				tSeqID FromSeqID;			// identifies the 'From' overlapping (FromVertexID) sequence
				tSeqID ToSeqID;			   // identifies the 'To' overlapped (ToVertexID) sequence or a completely contained sequence
				uint32_t FromSeqLen;			// 'From' sequence length is this length
				uint32_t ToSeqLen;			// 'To' sequence length is this length
				uint32_t FromSeq5Ofs;			// overlap of FromVertexID onto ToVertexID starts at this base relative to FromVertexID 5' start
				uint32_t FromSeq3Ofs;			// overlap of FromVertexID onto ToVertexID ends at this base relative to FromVertexID 5' start
				uint32_t ToSeq5Ofs;			// overlap onto ToVertexID from FromVertexID starts at this base relative to ToVertexID 5' start
				uint32_t ToSeq3Ofs;			// overlap onto ToVertexID from FromVertexID ends at this base relative to ToVertexID 5' start
				uint32_t Score;				// score associated with this overlap, higher scores represent higher confidence in the overlap
				uint32_t ScoreAlignLen;		// scored for this alignment length: (1 + ProbeAlignLen + TargAlignLen) / 2
				uint8_t flgFromAntisense:1;	// set if From was revcpl as antisense probe 
				uint8_t flgToAntisense:1;		// set if To was revcpl as antisense target 
				eOverlapClass OverlapClass;	// classification of overlap from OverlappingSeqID onto OverlappedSeqID
	} tsOverlappedSeq;

const int cAllocPathEdges = 10000;			// initally allocate for this many path edges per component; realloc as needed

typedef struct TAG_sPathTraceBack {
				 tComponentID ComponentID;      // path is for this component
				 tVertID VertexID;				// path includes this vertex sequence
				 uint32_t SeqLen;					// vertex sequence length
				 uint32_t PathOfs;				// vertex sequence is starting at this path offset (0..CurrentPathLength)
				 uint32_t Off5;					// vertex sequence to include in path starts from this sequence 5' offset (0..SeqLen-1)
				 bool bRevCpl;					// if true then reverse complement the vertex sequence before applying Off5 and Off3
	} tsPathTraceBack;

// disconnected graph components
typedef struct TAG_sComponent {
	tComponentID ComponentID;			// identifies this component
	tVertID VertexID;					// one of the vertices which is in this component
	uint32_t NumVertices;					// there are this many vertices in this component
	tVertID PathStartVertexID;			// highest scoring path starts from this vertex
	uint64_t PathScore;					// highest scoring path has this score
	uint32_t PathLength;					// highest scoring path is this length
	uint32_t NumPathEdges;				// number of edges in path
	uint32_t NumTraceBacks;				// number of tracebacks in path
	uint32_t StartTraceBackID;			// path starts with this traceback
} tsComponent;


typedef struct TAG_sRemapComponentID {
	tComponentID From;		// map from	
	tComponentID To;		// map to
} tsRemapComponentID;

#pragma pack()

class CAssembGraph
{
	CMTqsort m_MTqsort;				// multithreaded sorting

	bool m_bTerminate;				// if set true then all threads should terminate processing

	uint32_t m_MinScaffScoreThres;      // edges (overlaps) must be at least this score to be accepted for processing

	eVerticesSortOrder m_VerticesSortOrder; // current graph vertex sort order
	bool m_bVertexEdgeSet;				// true if InEdgeID/OutEdgeID have been initialised
	uint32_t m_UsedGraphVertices;			// number of graph vertices currently used
	uint32_t m_AllocGraphVertices;		// number of graph vertices allocated
	tsGraphVertex *m_pGraphVertices;   // allocated to hold array of graph vertices

	bool m_bSenseOnlyOvlps;				// if false then both sense/sense and sense/antisense overlaps will be accepted and processed, otherwise sense/sense only overlaps accepted and processed
	bool m_bOutEdgeSorted;				// true if m_pGraphOutEdges has been sorted in ascending FromVertexID.ToVertexOrder
	uint32_t m_UsedGraphOutEdges;			// number of forward graph edges currently used
	uint32_t m_AllocGraphOutEdges;		// number of forward graph edges allocated
	tsGraphOutEdge *m_pGraphOutEdges;	// allocated to hold array of forward edges

	bool m_bInEdgeSorted;				// true if m_pGraphInEdges has been sorted in ascending ToVertexOrder.FromVertexID
	uint32_t m_UsedGraphInEdges;			// number of inbound graph edges currently used
	uint32_t m_AllocGraphInEdges;			// number of inbound graph edges allocated
	tEdgeID *m_pGraphInEdges;			// index onto m_pGraphOutEdges which is sorted in ToVertexID.FwdVertexID ascending order

	uint32_t m_NumComponents;				// number of components
	uint32_t m_AllocComponents;			// number of components allocated
	tsComponent *m_pComponents;			// allocated to hold array of identified components

	uint32_t m_CurTransitDepth;			// current transition depth
	uint32_t m_MaxTransitDepth;			// deepest transition required
	uint32_t m_AllocTransitStack;			// allocation is for this many transition entries 
	tVertID *m_pTransitStack;			// allocated to hold transition entries

	uint32_t m_UsedTraceBacks;			// currently using this many tracebacks
	uint32_t m_AllocdTraceBacks;			// allocd to hold this many tracebacks
	tsPathTraceBack *m_pPathTraceBacks; // to hold all path tracebacks


	uint32_t m_NumDiscRemaps;				// number of disconnected graph identifiers requiring remaps
	tsRemapComponentID m_DiscRemaps[cMaxDiscRemaps];	// to hold disconnected graph identifiers requiring remaps
	uint32_t m_NumReducts;			// number of graph node reductions
	bool m_bReduceEdges;			// if true then attempt to reduce extraneous graph edges when current graph is full and more node memory needs to be allocated  

	bool m_bAcceptOrphanSeqs;		// if true then report also report sequences which have no overlap with any other sequence

	int m_NumThreads;				// use at most this number threads for graph processing
	bool m_bMutexesCreated;			// set true if mutexes and rwlocks created/initialised

#ifdef WIN32
	alignas(4)	volatile unsigned int m_CASSerialise; // used with synchronous compare and swap (CAS) for serialising access - replaces AcquireSerialise() as much more efficient
	alignas(4)	volatile unsigned int m_CASLock; // used with synchronous compare and swap (CAS) for serialising access -  - replaces AcquireLock(True) as much more efficient
#else
	__attribute__((aligned(4))) volatile unsigned int m_CASSerialise; // used with synchronous compare and swap (CAS) for serialising access - replaces AcquireSerialise() as much more efficient
	__attribute__((aligned(4))) volatile unsigned int m_CASLock; // used with synchronous compare and swap (CAS) for serialising access -  - replaces AcquireLock(True) as much more efficient
#endif
	void AcquireCASSerialise(void);
	void ReleaseCASSerialise(void);
	void AcquireCASLock(void);
	void ReleaseCASLock(void);


	static int SortOutEdgeFromVertexID(const void *arg1, const void *arg2);
	static int SortOutEdgeFromVertexIDSeqOfs(const void *arg1, const void *arg2);
	static int SortInEdgesToVertexID(const void *arg1, const void *arg2);
	static int SortOutEdgeToVertexIDSeqOfs(const void *arg1, const void *arg2);

	static int SortVerticesSeqID(const void *arg1, const void *arg2);
	static int SortVerticesVertexID(const void *arg1, const void *arg2);
	static int SortVerticesComponentID(const void *arg1, const void *arg2);
	static int SortComponentNumVertices(const void *arg1, const void *arg2);
	static int SortComponentID(const void *arg1, const void *arg2);

	tsGraphOutEdge *					// ptr to edge with matching FromVertexID and ToVertexID, or NULL if unable to locate				
		LocateFromToEdge(tVertID FromVertexID,	// match edge with this FromVertexID which is
				  tVertID ToVertexID);			// to this ToVertexID

	tsGraphOutEdge *							// ptr to edge with matching FromVertexID and ToVertexID, or NULL if unable to locate				
		LocateToFromEdge(tVertID ToVertexID,	// match edge with this ToVertexID which is 
				  tVertID FromVertexID);		// from this FromVertexID

	uint32_t			// index+1 in m_pFwdGraphEdges of first matching FromVertexID, or 0 if non matching				
		LocateFirstFwdEdgeID(tVertID FromVertexID);	// find first matching 
	
	uint32_t			// index+1 in m_pFwdGraphEdges of first matching DnSeqID, or 0 if non matching				
		LocateFirstDnSeqID(tSeqID DnSeqID);		// find first matching 


	tsGraphVertex *			// ptr to vertex corresponding to SeqID , or NULL if unable to locate				
			LocateVertexSeqID(tSeqID SeqID);		// match this SeqID

	tsGraphVertex *			// ptr to vertex corresponding to VertixID , or NULL if unable to locate				
			LocateVertex(tVertID VertixID);		// match this VertixID

	uint32_t					// number of replacements
		ReplaceComponentID(tComponentID ToReplaceID,		// replace existing disconnected graph identifiers 
							tComponentID ReplaceID);		// with this identifier

	uint32_t					// number of replacements
		ReplaceComponentIDs(void);

	int CreateMutexes(void);
	void DeleteMutexes(void);

	int			// stack depth or < 1 if errors
		PushTransitStack(tVertID VertexID);		// push VertexID onto stack 
	tVertID PopTransitStack(void);				// popped VertexID or 0 if stack was empty
	tVertID PeekTransitStack(void);				// peeked VertexID or 0 if stack empty
	void ClearTransitStack(void);					// remove all entries from stack

	uint32_t							// number of vertices with both inbound and outbound edges
		VertexConnections(void);		// identify and mark vertices which have multiple inbound edges
	uint32_t GenSeqFragment(tsGraphVertex *pVertex);		// initial seed vertex
	uint32_t  TransitIdentDiscGraph(tVertID VertexID,tComponentID ComponentID);
	uint32_t  ClearEdgeTravFwdRevs(void);
	uint32_t	ClearDiscCompIDs(void);

public:
	CAssembGraph(void);
	~CAssembGraph(void);

	void Reset(void);								// reset and free any allocated resources
	teBSFrsltCodes		// initialise with ScaffScoreThres and  maxThreads
		Init(bool bSenseOnlyOvlps = false,				// if false then both sense/sense and sense/antisense overlaps will be accepted and processed, otherwise sense/sense only overlaps accepted and processed
			 int ScaffScoreThres = cDfltMin1kScore,		// accepted edges must be of at least this overlap score
			 bool bAcceptOrphanSeqs = false,			// also report sequences which are not overlapped or overlapping any other sequence
						int MaxThreads = 8);			// set number of threads

	uint32_t								// returned vertex identifier
		AddVertex(uint32_t SeqLen,		// sequence length
				  tSeqID SeqID);		// allocates and initialises a new graph vertex which references this sequence

	uint32_t								// 0 if errors else number of vertices finalised
		FinaliseVertices(void);			// added vertices are sorted ready for graph edges to be added

	uint32_t								// number of edges removed	 
		ReduceEdges(void);				// reduce graph by detecting and removing extraneous edges


	uint32_t									// returns total number of edges, including this edge if accepted, thus far accepted 
		AddEdge(tSeqID FromSeqID,			// identifies the 'From' or overlapping sequence
				tSeqID ToSeqID,				// identifies the 'To' overlapped sequence
				uint32_t FromSeqLen,			// 'From' sequence length is this length
				uint32_t ToSeqLen,			// 'To' sequence length is this length
				uint32_t Score,				// score associated with this overlap, higher scores represent higher confidence in the overlap
				uint32_t ScoreAlignLen,		// scored for this alignment length: (1 + ProbeAlignLen + TargAlignLen) / 2
				uint32_t FromSeq5Ofs,			// overlap of FromSeqID onto ToSeqID starts at this base relative to FromSeqID 5' start
				uint32_t FromSeq3Ofs,			// overlap of FromSeqID onto ToSeqID ends at this base relative to FromSeqID 5' start
				uint32_t ToSeq5Ofs,			// overlap onto ToSeqID from FromSeqID starts at this base relative to ToSeqID 5' start
				uint32_t ToSeq3Ofs,			// overlap onto ToSeqID from FromSeqID ends at this base relative to ToSeqID 5' start
				eOverlapClass OverlapClass,	// classification of overlap from FromSeqID onto ToSeqID, note that classification must be eOLCOverlapping
				bool bSenseOvlpAnti);	// false: 'From' sense overlaps 'To' sense, true: 'From' sense overlaps 'To' antisense 

	uint32_t									// returns total number of edges , including any of these edges, if accepted, thus far accepted 
		AddEdges(uint32_t NumSeqs,				// number of overlapped sequences
				tsOverlappedSeq *pOverlappedSeqs); // overlapped sequences

	uint32_t									// returns total number of edges accepted after finalisation
		FinaliseEdges(void);
	
	uint32_t GetNumGraphVertices(void);		// returns current number of graph vertices

	uint32_t GetNumGraphOutEdges(void);		// returns current number of graph forward edges

	uint32_t GetNumReducts(void);				// returns current number of edge reductions

	uint32_t IdentifyDiscComponent(tVertID VertexID, // start component traversal from this vertex
				tComponentID ComponentID);			 // mark all traversed vertices as members of this component

	int32_t											// returned From sequence extension; -1 if no sequence extension
	OverlapAcceptable(tsGraphOutEdge *pEdge,		// overlap edge
				uint8_t FromOvlpClass = 0,	// From vertex overlap classification; bit 0 set if From vertex evaluated as antisense in current path
				uint8_t *pToOvlpClass = NULL);		// returned To vertex overlap classification; bit 0 set if To vertex is evaluated as antisense in current path

	uint64_t
		ScorePaths(tVertID VertexID);			// score all paths starting with outgoing edges from this vertex

	uint64_t										// highest scoring of any path from pEdge
		ScorePath(uint32_t Depth,					// current recursive depth - used to detect circular paths
				  tsGraphOutEdge *pEdge,		// score paths starting with this edge
				uint8_t FromOvlpClass = 0,	// From vertex overlap classification; bit 0 set if From vertex evaluated as antisense in current path
				uint8_t *pToOvlpClass = NULL);	// returned To vertex overlap classification; bit 0 set if To vertex is evaluated as antisense in current path



	int										 // eBSFSuccess or otherwise
		FindHighestScoringPaths(void);		 // score all possible paths and record highest scoring path for each component

	int											 // eBSFSuccess or otherwise
			ReportVerticesEdgesGEXF(char *pszOutFile,		 // report as GEXF format all vertices and their edges for each component to this file
								CSeqStore *pSeqStore);	// holds sequences used to assemble contig

	int												 // eBSFSuccess or otherwise
			ReportVerticesEdgesGraphML(char *pszOutFile,		 // report as GraphML on all vertices and their edges for each component to this file
										   CSeqStore *pSeqStore);	// holds sequences used to assemble contig

	int											// returned CurrentPathLength
		AddTraceBackPath(bool bFirst,			// true if 1st vertex in path (starting a new path)
				 tComponentID ComponentID,      // path is for this component
				 tVertID VertexID,				// path includes this vertex sequence
				 uint32_t SeqLen,					// vertex sequence length
				 uint32_t PathOfs,				// vertex sequence is starting at this path offset (0..CurrentPathLength)
				 uint32_t Off5,					// vertex sequence to include in path starts from this sequence 5' offset (0..SeqLen-1)
				 bool bRevCpl);					// if true then reverse complement the vertex sequence before applying Off5

	int										// eBSFSuccess or otherwise
		GenTraceBackPath(tsComponent *pComponent); // generate traceback path for this component

	uint32_t	IdentifyDiscComponents(void);

	int WriteContigSeqs(char *pszOutFile,CSeqStore *pSeqStore);  // write assmbled PacBio contig sequences to this output file
};



