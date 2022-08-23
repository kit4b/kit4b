#pragma once

#include "BKScommon.h"

const int cListenBacklog = 5;					// maximum length of the queue of pending connections

const int cMaxConcurrentRequests = min(4095,(int)(cMaxServiceInsts * cMaxNumSessions));	// can process at most this many concurrent service requests over all session instances independent of service type
const int cMaxReqID = cMaxConcurrentRequests;	// request identifiers will range from 1..cMaxConcurrentRequests

// when negotiating with potential service providers then minimal buffer tx/rx buffer sizes are allocated
const int cMinTxRxBuffSize = (cMaxServiceTypes * sizeof(tsServiceDetail)) + sizeof(tsBKSReqServices) * 3;	// always allocate at least this sized TxdBuff/RxdBuffs - ensures negotiation frames fit!

// Session service providers will be in one of these exclusive states
typedef enum TAG_eBKSEPProvState
{
	eBKSPSUndefined = 0,				// service provider state is undefined waiting on a potential service provider to connect
	eBKSPSRegisteringNeg,				// initiated negotiation of capabilities with connected service provider
	eBKSPSRegisteredActv,				// completed capabilities negotiation and provider can now actively process service requests
	eBKSPSRegisteredTerm,				// currently terminating service provider connection
	eBKSPSPlaceHolder
} teBKSPEPProvState;

typedef enum TAG_eSessEstabState
{
	eSESnone = 0,			// not yet started
	eSESTxReqServices,		   // send required service type details to potential provider
	eSESRxOfferedService,	   // expecting to receive potential providers offered service
	eSESTxAcceptService,   // sending acceptance of offered service to provider
	eSESTxRejectService,   // sending rejection of offered service to provider
	eSESPlaceHolder
} teSessEstabState;


#pragma pack(1)

typedef int64_t tJobIDEx;			// >0 extended job identifier; 0 if none assigned, <0 if errors

typedef struct TAG_sTxdRxd
{
	uint32_t SessionID;			// uniquely identifying this session between service requester and provider, requester specified
	uint16_t flgSelMonExcept : 1; // socket to be monitored for exceptions
	uint16_t flgSelMonRead : 1;	// socket to be monitored for received data available to be read
	uint16_t flgSelMonWrite : 1; // socket to be monitored for write completed
	uint16_t flgKeepAliveReq: 1; // set if a keepalive is to be sent
	uint16_t flgRxCplt : 1;	// set when complete frame received 	
	uint16_t flgTxCplt : 1;	// set when complete frame sent
	uint16_t flgErr : 1;       // set on any unrecoverable error
	uint16_t flgErrReason : 4;	// holds reason for socket level error flag set
	socket_t  Socket;		// assumed connected socket
	time_t PacRxdAtSecs;	// the time at which a frame was last received, used for determining if session still active
	time_t PacTxdAtSecs;	// the time at which a frame was last sent, used for keep alive generation

	uint8_t RxdRxFrameID;		// RxFrameID of last received frame - will be 0 if no frame yet received by session peer
	uint8_t RxdTxFrameID;		// TxFrameID of last received frame - will be 0 if no frame yet received from session peer
	uint8_t TxFrameID;		// FrameID of next frame to send, will wraparound back to 1 when incremented past 0x07f 
	uint32_t TotRxd;			// pRcvBuff contains this total number of received bytes
	uint32_t CurPacRxd;		// number of bytes in currently received packet
	uint32_t TotTxd;			// total number of bytes in pTxdBuff to be sent
	uint32_t CurTxd;			// number of bytes currently sent from pTxdBuff
	uint32_t AllocdRxdBuff;	// pRxdBuff allocated to hold at most this many received bytes
	uint32_t AllocdTxdBuff;	// pTxdBuff allocated to hold at most this many bytes
	uint8_t *pRxdBuff;		// receiving data into this buffer
	uint8_t *pTxdBuff;		// sending data from this buffer
    SOCKADDR_STORAGE  IPaddress;	// remote IP address + port of endpoint service provider (IPv4 or IPv6)
} tsTxdRxd;


typedef struct TAG_sBKSSessEstab
{
	uint8_t SEState;			// current session establishment state - eSESnone, eSESTxServices, eSESRxOfferService, eSESTxAcceptService 
	time_t StartSecs;		// started session establishment at this time
	teBKSPType BKSPType;	// confirmation of service type being accepted or eBKSPTUndefined if offered type not accepted 
	uint32_t MaxInstances;	// max number of instances supported
	uint32_t MaxClassInstances;// Session can instantiate at most this many class instances
	tsTxdRxd TxdRxd;		// holding low level send/receive buffers + connected socket
} tsBKSSessEstab;

// Session service providers are registered into the following
typedef struct TAG_sBKSRegSession
{
	uint32_t SessionID;					// uniquely identifying this session between service requester and provider
	uint32_t TypeSessionID;				// used as an index within the type to reference this session
	uint8_t BKSPType;						// Session is providing this teBKSType service
	uint8_t BKSPState;					// Session is currently in this teBKSPEPProvState registration state 
	uint32_t MaxInstances;				// Session can process at most this many service instances
	uint32_t MaxClassInstances;			// Session can instantiate at most this many class instances
	uint32_t NumBusy;						// total number of instances currently committed to requests, processing, or with responses
	uint32_t NumReqs;						// number of instances with service request ready to send to a service provider
	uint32_t NumProcs;					// number of instances with service requests currently being processed by service provider
	uint32_t NumCpltd;					// number of instances with service requests finished processing with response available
	uint64_t TotNumCpltd;					// total number of completed requests in this session
} tsBKSRegSession;

typedef struct TAG_sReqRespInst
{
	int64_t JobIDEx;		// unique job identifier as was generated when job was accepted for submission
	uint32_t SubmitAt;	// when job was accepted for submission
	uint32_t CpltdAt;		// when job was returned as completed
	uint16_t ReqID;		// instance specific identifier
	uint16_t FlgReq:1;	// this job instance has been initialised and is ready to be sent for processing
	uint16_t FlgProc:1;   // this job instance is currently being processed by service provider
	uint16_t FlgCpltd:1;  // service provider has completed processing and returned results are available
	uint32_t JobRslt;			// service provider completion result
	uint64_t ClassInstanceID;	// class instance referenced
	uint32_t ClassMethodID;	// identifies class method
	uint32_t ParamSize;		// instance specific parameter size
	uint32_t InDataSize;		// instance specific input data size
	uint32_t OutDataSize;		// instance specific result data size
	uint8_t Data[1];		// when service requested then parameters followed by input data, if service response then response result data
} tsReqRespInst;

typedef struct TAG_sBKSRegSessionEx
{
	struct TAG_sBKSRegSessionEx *pNext;	// Sessions are linked as will be dynamically allocated
	tsBKSRegSession Session;			// Session

	uint32_t LastChkdReqIdx;				// start checking for service requests to send from this instance index
	uint32_t NumClassInstances;			// number of class instance identifiers currently in ClassInstanceIDs
	uint64_t ClassInstanceIDs[cMaxClassInsts];	// holds all instantiated class instance identifiers for this session
	uint32_t AllocdReqResp;				// allocation size for pReqResp
	uint8_t *pReqResp;					// allocation for NumInstances of requests and associated responses
	tsTxdRxd TxdRxd;					// holding low level send/receive buffers + connected socket
} tsBKSRegSessionEx;


typedef struct TAG_sBKSType
{
	tsServiceDetail Detail;				// service type detail
	uint32_t MaxSessions;					// allowing at most this many sessions of specified service type
	uint32_t  NumSessions;				// currently this number of registered sessions of this type
	uint32_t NumInstances;				// totaling this many service instances
	uint32_t ReqRespInstSize;				// each service instance requires this much memory
	uint32_t FlgReq:1;					// at least 1 session has a service request ready to send to a service provider
	tsBKSRegSessionEx *pFirstSession;	// ptr to first session of this type
	tsBKSRegSessionEx *pSessions[cMaxNumSessions]; // ptr to each session as indexed by TypeSessvIdx  
} tsBKSType;


typedef struct TAG_sChkPtReq {	// individual check pointed request
	uint32_t ClassMethodID;	// identifies class method
	uint32_t ParamSize;		// instance specific parameter size
	uint32_t InDataSize;		// instance specific input data size
	uint32_t MaxOutDataSize;	// instance specific max expected output data for this method
	uint8_t Data[1];			// parameters followed by input data
	} tsChkPtReq;

typedef struct TAG_sChkPtReqs { // list of all check pointed requests for a specific class instance
	uint32_t SessionID;			// uniquely identifying this session between service requester and provider (set to 0 if this class instance can be reused)
	uint32_t TypeSessionID;		// used as an index within the type to reference this session
	uint64_t ClassInstanceID;		// check pointing this class instance
	uint32_t JobRslt;				// last job processing result as returned by service provider
	uint32_t UsedRespData;		// pRespData currently holds this sized response for last job submitted
	uint32_t AllocRespData;		// pRespData allocated to this size
    uint8_t *pRespData;			// response data for last check pointed method successfully processed by service provider
	uint32_t NumChkPtReqs;		// number of check pointed requests in ChkPtReqs[]
	size_t AllocMem;			// total memory allocated to hold this check pointed list (includes this header)
	size_t UsedMem;				// currently using this sized memory to hold check pointed list (includes this header)
	tsChkPtReq ChkPtReqs[1];	// 1st check pointed request
	} tsChkPtReqs;

typedef struct TAG_sRequesterThread {
	int ThreadIdx;					// uniquely identifies this thread
	void *pThis;					// will be initialised to pt to class instance
#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	uint32_t threadID;				// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int Rslt;						// processing result
} tsRequesterThread;

#pragma pack()

class CBKSRequester
{
	char m_szHostName[cMaxHostNameLen+1];				// host on which to listen for connections
	char m_szServiceName[cMaxServiceNameLen+1];			// listening on this port

	tsRequesterThread m_RequesterThread;				// requester thread parameters

#ifdef WIN32
#ifdef _ChkLockDepth_			// if checking that serialisation locks are actually working!
	alignas(4) volatile	uint32_t m_LockDepth;
	alignas(4) volatile	uint32_t m_CASSerialiseDepth;
#endif
	alignas(4) volatile	uint32_t m_CreatedMutexes;					// set to 1 after serialisation locks/mutexes initialised
	alignas(4) volatile	uint32_t m_ThreadActive;						// set to 1 after requester thread has started and initialised
	alignas(4) volatile uint32_t m_Terminate;							// set to 1 if requester thread is required to self terminate
	alignas(4) volatile uint32_t m_LockEnabled;						// set to 1 if AcquireLock has been initialised and can be used for serialisation
	alignas(4) volatile uint32_t m_CASSerialise;						// used with synchronous compare and swap (CAS) for serialising access 
	alignas(4) volatile uint32_t m_NumPendingReqs;					// number of outstanding pending requests from RMI worker threads to be processed
	alignas(4) volatile uint32_t m_NumPendingResps;					// number of outstanding pending response checks from RMI worker threads to be processed
	alignas(4) volatile LONG m_TotRespsAvail;						// current total number of responses available to RMI worker threads over all sessions

#else
#ifdef _ChkLockDepth_		// if checking that serialisation locks are actually working!
	__attribute__((aligned(4))) volatile	uint32_t m_LockDepth;
	__attribute__((aligned(4))) volatile	uint32_t m_CASSerialiseDepth;
#endif
	__attribute__((aligned(4))) volatile uint32_t m_CreatedMutexes;	// set to 1 after serialisation locks/mutexes initialised
	__attribute__((aligned(4))) volatile uint32_t m_ThreadActive;		// set to 1 after requester thread has started and initialised
	__attribute__((aligned(4))) volatile uint32_t m_Terminate;		// set to 1 if requester thread is required to self terminate
	__attribute__((aligned(4))) volatile uint32_t m_LockEnabled;		// set to 1 if AcquireLock has been initialised and can be used for serialisation
	__attribute__((aligned(4))) volatile uint32_t m_CASSerialise;		// used with synchronous compare and swap (CAS) for serialising access 
	__attribute__((aligned(4))) volatile uint32_t m_NumPendingReqs;	// number of outstanding pending requests from RMI worker threads to be processed
	__attribute__((aligned(4))) volatile uint32_t m_NumPendingResps;	// number of outstanding pending response checks from RMI worker threads to be processed
	__attribute__((aligned(4))) volatile uint32_t m_TotRespsAvail;	// current total number of responses available to RMI worker threads over all sessions
#endif

	bool m_bNotifiedReqs;								// set TRUE if control message has been sent to control socket 2

	bool m_bSessionTermReq;								// true whilst a Session  is being deleted

	uint32_t m_NumChkPtReqs;								// m_pChkPtReqs contains this many currently active check pointed class request ptrs
	uint32_t m_MaxChkPtReqs;								// m_pChkPtReqs allocated to hold at most this many check pointed class request ptrs
	tsChkPtReqs **m_ppChkPtReqs;							// allocated to hold check pointed class request ptrs

	uint32_t m_NumInitTypes;								// number of service types initialised in m_pBKSTypes
	tsBKSType *m_pBKSTypes;								// pts to allocated array of entries for each potentially supported service type indexed by Type-1

	uint32_t m_NumSessEstabs;								// number of sessions currently in pBKSSessEstab being negotiated with
	tsBKSSessEstab *m_pBKSSessEstabs;					// pts to allocated array of session peer endpoints currently being negotiated with
 
	uint32_t m_NumSessions;								// currently there are this many registered sessions

	uint32_t m_SessionIDVect[(cMaxConcurrentSessions+31)/32];	// bit vector containing all possible session identifiers, if bit set then that session identifier has been allocated and is in use
	uint32_t m_ReqIDVect[(cMaxReqID+31) / 32];				// bit vector containing all possible request identifiers, if bit set then that request identifier has been allocated and is in use

	socket_t m_ListenerSock;							// listening on this socket for Session connections

	tsTxdRxd m_Ctrl[2];									// Ctrl[0] written to by threads needing to signal select() processing thread, select() processing thread monitors m_Ctrl[1]

	uint32_t												// returned request identifier or 0 if all identifiers have already been allocated
			AllocReqID(void);							// returns next available unused request identifier and sets that identifier in m_ReqIDVect[] as now allocated

	bool												// true if ReqID was previously allocated (used), false if previously unallocated (unused)
			UnallocReqID(uint32_t ReqID);					// sets the ReqID in m_ReqIDVect as unallocated and available to be reused

	bool												// true if ReqID is currently allocated (used), false if unallocated (unused)
			IsAllocReqID(uint32_t ReqID);					// checking this ReqID in m_ReqIDVect

	uint32_t												// returned session identifier or 0 if all identifiers have already been allocated
		AllocSessionID(void);							// returns next available unused session identifier and sets that identifier in m_SessionIDVect as allocated

	bool				// true if SessionID was previously allocated (used), false if previously unallocated (unused)
		UnallocSessionID(uint32_t SessionID);   // sets the SessionID in m_SessionIDVect as unallocated and available to be reused

	bool												// true if SessionID is currently allocated (used), false if unallocated (unused)
			IsAllocSessionID(uint32_t SessionID);			// checking this SessionID in m_SessionIDVect

	tsBKSRegSessionEx *LocateSession(uint32_t SessionID);	// locates an established session which is identified by SessionID

	void
		ResetSessEstab(tsBKSSessEstab *pSessEstab,			// reset this tsBKSSessEstab instance
									  bool bKeepAllocs = true,// if true then retain any existing buffer allocations
									  bool bKeepSessionID = false);	// if true then retain SessionID as session has been accepted
	bool
		StartSessEstab(uint32_t SessionID,								// Session identifier for this potential Session
									socket_t Socket,					// communicating with Session over this connected socket
									SOCKADDR_STORAGE  *pIPaddress);		// peer is at this socket network address

	bool
		ProgressSessEstab(tsBKSSessEstab *pSessEstab,
										bool bCpltdWrite);	// false if frame received, true if frame sent

	bool
		AcceptFullSession(tsBKSSessEstab *pSessEstab);		// accepting session being established as full session ready for normal payload request/responses

	int	TerminateAllSessions(void);

	int				// returns number of sessions deleted
		DeleteAllSessionsInState(teBKSPEPProvState TermState = eBKSPSRegisteredTerm);	// terminate and delete any session in specific state (normally eBKSPSRegisteredTerm) or with an invalid socket

	int	CreateMutexes(void);				// create mutexes used in access serialisations
	int	DeleteMutexes(void);				// delete mutexes used in access serialisations

	bool AcquireLock(bool bExclusive);	// lock for serialised access by multiple concurrent reader threads (bExclusive == false), or serialised access by single thread (bExclusive == true)
	bool ReleaseLock(bool bExclusive);	// release serialised access lock

	void AcquireCASSerialise(bool bPriority = false);	// if bPriority true then backoff time is reduced relative to if bPriority is false, increasing the probability of acquiring the serialisation lock if there is contention
	void ReleaseCASSerialise(void);

#ifdef _WIN32
	SRWLOCK m_hRWLock;					// serialising multiple reader access but single writer
#else
	pthread_rwlock_t m_hRWLock;			// serialising multiple reader access but single writer
#endif

	bool InitialiseListener(const char* pszHost, char *pszService);  // initialise listening socket
	bool InitialiseCtrlSocks(void); // initialise the control sockets in m_Ctrl[]
	int	AcceptConnections(void);	 // start accepting connections
	int	SendRequestFrames(void);			// iterate all sessions and if any frames ready to send and room to accept the frame in TxdBuff then initiate the sending
	int  // 0: accepted frame, -1: JobIDEx errors, -2 Session or type errors, -3 mismatch between instance JobIDEx's, ClassInstanceID mismatch
			ProcessResponseFrame(tsBKSRegSessionEx *pSession);			// process a received response frame on this session
	bool RxData(tsTxdRxd *pRxd);
	bool TxData(tsTxdRxd *pTxd);

	bool NotifyCtrl(uint8_t Msg = 0);			// notify via control sockets that, default, there is at least 1 request available to be sent to a session instance provider
	bool SendCtrlMsg(int Len,				// number of control message bytes to send ptd at by pCtrlMsg, can be 0
								uint8_t *pCtrlMsg);	// pCtrlMsg pts to control message bytes
	int								// received message is this length, 0 if no outstanding messages, < 0 if errors
			RcvCtrl(int BuffLen,			// available buffer into which the control message can be copied 
		uint8_t *pBuff);			// copy control message into this buffer
	bool ProcessCtrlMsg(int MsgLen,uint8_t *pMsg);			// process a message received by control socket m_Ctrl[1] - note: currently these messages are simply discarded

	uint32_t									// identifies created and initialised check point requests list
		CreateChkPtReqs(uint32_t SessionID,	// uniquely identifying this session between service requester and provider
				uint32_t TypeSessionID,		// used as an index within the type to reference this session
				uint64_t ClassInstanceID,		// check pointing this class instance
			    uint32_t ClassMethodID,		// identifies class method
				uint32_t ParamSize,			// instance specific parameter size
				uint32_t InDataSize,			// instance specific input data size
				uint32_t MaxOutDataSize,		// instance specific max expected output data for this method
				uint8_t *pData);				// request parameters followed by input data

	bool
			DeleteChkPtReqs(uint32_t ChkPtID);  // identifies check point list to be deleted

	bool			// true if request was successfully check pointed
		AddChkPtReq(uint32_t ChkPtID,  // identifies check point list to extend with this request
				uint32_t SessionID,			// uniquely identifying this session between service requester and provider
			    uint32_t ClassMethodID,		// identifies class method
				uint32_t ParamSize,			// instance specific parameter size
				uint32_t InDataSize,			// instance specific input data size
				uint32_t MaxOutDataSize,		// instance specific max expected output data for this method
				uint8_t *pData);				// request parameters followed by input data

	bool			// true if response was successfully updated		
		AddChkPtResp(uint32_t ChkPtID,  // identifies check point list to update with this response
					uint32_t SessionID,			// uniquely identifying this session between service requester and provider
					uint32_t JobRslt,				// job processing result as returned by service provider
					uint32_t RespSize,			// response size
					uint8_t *pData);				// response data

	bool ShutdownConnection(socket_t *pSocket);

	int			// returns 0 if no sockets to be monitored with select(), on windows the total number of monitored sockets, on linux the highest socket file descriptor plus 1
		SetupFDSets(fd_set& ReadFDs,			// select() read available socket descriptor set  
					fd_set& WriteFDs,			// select() write accepted socket descriptor set
					fd_set& ExceptFDs);			// select() exceptions descriptor set


	tJobIDEx										// packed job identifier or 0 if range errors
				PackIntoJobIDEx(uint32_t ReqID,							// must be in the range 1..16777215 (24bits)
									   uint32_t SessionID,				// service provider session identifier, must be in the range 1..131071 (17bits)
									   uint32_t InstanceID,				// service instance in the service providers session SessionID, must be in the range 1..511 (9bits)
									   uint32_t TypeID,					// identifies service being provided by the service provider, must be in the range 1..15 (4bits)
									   uint32_t TypeSessionID);			// index of this session in the service type, must be in the range 1..511 (9bits)

	bool				// false if any range errors whilst unpacking
			UnpackFromJobIDEx(tJobIDEx JobIDEx,								// unpack from this extended job identifier
										 uint32_t *pReqID,				    // returned ReqID, will be in the range 1..16777215
										 uint32_t *pSessionID,				// returned SessionID, will be in the range 1..131071
										 uint32_t *pInstanceID,				// returned InstanceID, will be in the range 1..511
										 uint32_t *pTypeID,					// returned service TypeID, will be in the range 1..15
										 uint32_t *pTypeSessionID);			// returned index of this session in the service type, will be in the range 1..511

#ifdef WIN32	
	const char *WSAGetLastErrorMessage(const char* pcMessagePrefix, int nErrorID = 0);
#endif
	
	int	Reset(bool bReleaseLock);			// locked mutex requires releasing immediately prior to deleting the mutexes

	

public:
	CBKSRequester();
	~CBKSRequester();

	int StartThread(tsRequesterThread *pPars);


						// NOTE: defaults chosen are representative for Smith-Waterman service processing
	int					// returns total number of registered service types or teBSFrsltCodes error code if any parameterisation errors or already registered type
		RegServiceType(teBKSPType BKSPType = eBKSPTSmithWaterman,							// registering this service type
									  uint32_t MinProviderVersion = cMinProviderVersion,		// service provider version must be at least this software version
									  uint32_t MaxProviderVersion = cMaxProviderVersion,		// service provider version must be no more than this software version
									  uint32_t KeepaliveSecs = cDfltKeepaliveSecs,			// expecting packet activity from session peer with periodicity of no more than this number of seconds
									  uint32_t MaxSessions = cDfltNumSessions,				// allowing at most this many concurrent sessions of specified service type
									  uint32_t MinServiceInsts = cMinServiceInsts,			// any session to support at least this minimum number of service instances
									  uint32_t MaxServiceInsts = cDfltServiceInsts,			// limit any session to support a maximum of this many service instances
									  uint32_t MaxProcSecs = cMaxSWProcSecs,                   // expecting any service provider to take no more than this number of seconds to complete processing a given request
									  uint32_t MaxParamLen = cMaxSWParamLen,					// service job parameters can be up to this length
  									  uint32_t MaxQuerySeqLen = cMaxSWQuerySeqLen,			// query sequences can be up to this length
									  uint32_t MaxTargSeqLen = cMaxSWTargSeqLen,				// target sequences can be up to this length
								      uint32_t MaxReqPayloadSize = cMaxSWReqPayloadSize,		// request payloads to the service provider, including framing, can be up to this size (UINT8s),
									  uint32_t MaxRespPayloadSize = cMaxSWRespPayloadSize);	// response payloads from the service provider, including framing, can be up to this  size (UINT8s)

	int						// returns number of registered sessions which are in requested state
		GetNumSessions(teBKSPType BKSPType = eBKSPTSmithWaterman,					// may request number of sessions providing this specific service type, or if eBKSPTUndefined then the total
					 teBKSPEPProvState BKSPState = eBKSPSRegisteredActv);			// sessions must be in this specific state

	int						// returns total number of service instances in all registered sessions which are in requested state
		GetNumInstances(teBKSPType BKSPType = eBKSPTSmithWaterman,					// may request number of instances providing this specific service type, or if eBKSPTUndefined then the total
					   teBKSPEPProvState BKSPState = eBKSPSRegisteredActv);			// sessions providing service instances must be in this specific state

	int													// actual number of sessions returned
		GetSessions(teBKSPType BKSPType,				// return sessions for this service type,  or if eBKSPTUndefined then all service types registered
					 teBKSPEPProvState BKSPState,		// sessions must be in this state
				int MaxSessions,						// return at most this many sessions in pSessions[]
				tsBKSRegSession *pSessions);			// caller has preallocated to hold returned array snapshot of sessions

	int 
			Initialise(char* pszHost = NULL,			// listening on this host/IP address; NULL to use first INET IP local to this machine
						char *pszService = NULL,			// listening on this service/port; NULL to use default port 
						teBKSPType BKSPType = eBKSPTSmithWaterman,					// registering this service type
									  uint32_t MinProviderVersion = cMinProviderVersion,		// service provider version must be at least this software version
									  uint32_t MaxProviderVersion = cMaxProviderVersion,		// service provider version must be no more than this software version
									  uint32_t KeepaliveSecs = cDfltKeepaliveSecs,			// expecting packet activity from session peer with periodicity of no more than this number of seconds
									  uint32_t MaxSessions = cDfltNumSessions,				// allowing at most this many concurrent sessions of specified service type
									  uint32_t MinServiceInsts = cMinServiceInsts,			// any session to support at least this minimum number of service instances
									  uint32_t MaxServiceInsts = cDfltServiceInsts,			// limit any session to support a maximum of this many service instances
									  uint32_t MaxProcSecs = cMaxSWProcSecs,                   // expecting any service provider to take no more than this number of seconds to complete processing a given request
									  uint32_t MaxParamLen = cMaxSWParamLen,					// service job parameters can be up to this length
  									  uint32_t MaxQuerySeqLen = cMaxSWQuerySeqLen,			// query sequences can be up to this length
									  uint32_t MaxTargSeqLen = cMaxSWTargSeqLen,				// target sequences can be up to this length
								      uint32_t MaxReqPayloadSize = cMaxSWReqPayloadSize,		// request payloads to the service provider, including framing, can be up to this size (UINT8s),
									  uint32_t MaxRespPayloadSize = cMaxSWRespPayloadSize);	// response payloads from the service provider, including framing, can be up to this  size (UINT8s)

	bool Run(int Secs=180);					// false if unable to begin thread execution or if executing thread took too long (default 3min) to initialise ready to accept connections

	void Terminate(int Secs = 120);			// allow this many seconds for thread to self-terminate before forcing termination

	int												// total number of classes 
		GetNumClassInstances(teBKSPType TypeID = eBKSPTSmithWaterman,		// service required
						uint32_t *pCommited = NULL,			// returned number of class instances currently committed or instantiated 
						uint32_t *pUncommited = NULL);    // returned number of class instances currently not committed and available to be instantiated

	int				// -2: parameter errors, -1: class instance no longer exists, 0: currently no available service instance 1: if job accepted
			AddJobRequest(  tJobIDEx *pJobID,					// returned unique job identifier by which job can later be referenced
							teBKSPType TypeID,					// service type required
							  uint64_t ClassInstanceID,			// class instance on which job method is to be applied
							  uint32_t ClassMethodID,				// class method to apply on the class instance
									 uint32_t ParamsSize = 0,		// processing parameters are this total size in bytes
									 void *pParams = NULL,		// service processing parameters
									 uint32_t InDataSize = 0,		// service processing input data is this total size in bytes
									 void *pInData = NULL,		// service processing input data
									 uint32_t SessionID = 0);     // if 0 then use session as specified by ClassInstanceID, otherwise use session corresponding to specific session identifier

	int				// < 0 if job no longer exists, 0 if job still being processed, > 0 if job completed
		GetJobResponse(tJobIDEx	JobID,			// unique job identifier returned when job was submitted
				  uint64_t *pClassInstanceID,		// returned class instance on which job method was applied
				  uint32_t *pClassMethodID,		// returned class method applied
					   uint32_t *pJobRslt,		// job processing result as returned by service provider
						uint32_t *pOutDataSize,	// (IN) service processing output results expected to be at most this total length, [OUT] bytes of response data copied into pOutData 
						void *pOutData,			// service processing output results data
						bool bRetain=false);		// true if job response is to be retained and not deleted; subsequent call with bRetain==false will delete this response


};

