#pragma once

#include "BKScommon.h"

const uint32_t cServiceProviderVersion = 1;			// service provider is at this version
const uint32_t cMaxServiceProviderInsts = cMaxServiceInsts;	    // limited to support a maximum of this many service instances
// when negotiating with potential service requesters then minimal buffer tx/rx buffer sizes are allocated
const int cMinTxRxBuffSize = (cMaxServiceTypes * sizeof(tsServiceDetail)) + sizeof(tsBKSReqServices) * 3;	// always allocate at least this sized TxdBuff/RxdBuffs - ensures negotiation frames fit!

const int cMaxReqDataSize =  cMaxSWReqPayloadSize;		// each worker thread allocates to process up to this much request data
const int cMaxReqParamSize = cMaxSWParamLen;			// each worker thread allocates to process up to this much parameterisation data
const int cMaxRespDataSize = cMaxSWRespPayloadSize;		// each worker thread allocates to return up to this much response data
const int cMaxMFABuffSize =  cMaxSWMAFBuffSize;			// each worker thread allocates to hold at most this sized MAlignCols2fasta/MAlignCols2MFA alignments plus row descriptor prefixes

// service providers will be in one of these exclusive states
typedef enum TAG_eBKSPProvState
{
	eBKSPSUndefined = 0,				// service provider state is undefined, yet to connect with a server
	eBKSPSWaitReqServices,				// connected with server and waiting for server to send list of required services
	eBKSPSSendOfferedService,			// sending server the offered service
	eBKSPSWaitAcceptService,			// waiting for server to accept offered service
	eBKSPSAcceptedServiceActv,			// server has accepted offered service so can now actively process service requests
	eBKSPSAcceptedServiceTerm,			// session is in the process of terminating
	eBKSPSPlaceHolder
} teBKSPProvState;

#pragma pack(1)
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

typedef struct TAG_sReqResp {
	int64_t JobIDEx;						// request identifier as received from server for this request instance
	uint64_t ClassInstanceID;				// class instance referenced
	uint32_t ClassMethodID;				// identifies class method
	uint32_t InstanceID;					// service instance specific identifier
	uint32_t InstanceIDEx;				// combination of both the InstanceID (in bits 0..9) and ......
	uint32_t flgReqAvail : 1;				// this service instance is available for processing 
	uint32_t flgProc: 1;					// this service instance is currently being processed
	uint32_t flgCpltd: 1;					// service processing has completed and resultset can be sent back to service requester
	uint32_t JobRslt;						// completion result
	uint32_t ParamSize;					// instance specific parameter size
	uint32_t InDataSize;					// instance specific input data size
	uint32_t OutDataSize;					// instance specific result data size
	uint8_t Data[1];						// when service requested then parameters followed by input data, if service response then response result data
} tsReqResp;

typedef struct TAG_sBKSRegSessionEx
{
	uint8_t BKSPType;						// session is providing this teBKSType service
	uint8_t BKSPState;					// session is currently in this teBKSPProvState registration state 
	uint32_t NumInstances;				// session can process at most this many service instances
	uint32_t MaxReqPayloadSize;			// request payloads from the service requester, including framing, can be up to this size (UINT8s),
	uint32_t MaxRespPayloadSize;			// response payloads from the service provider to requester, including framing, can be up to this  size (UINT8s)
	uint32_t ReqRespInstSize;				// each service instance requires this much memory
	uint32_t KeepaliveSecs;				// expecting packet activity from endpoint with periodicity of no more than this number of seconds
	uint32_t InstancesBusy;				// currently this many requests being serviced in this session (total of instances ready to be processed, being processed, and processing has completed)
	uint32_t InstancesReqAvail;			// current number of requests ready to be processed
	uint32_t InstancesProc;				// current number of requests currently being processed
	uint32_t InstancesCpltd;				// current number of instances completed processing and ready for resultsets to be sent back to requester
	uint32_t TotNumRequests;				// total number of requests processed by this service provider in the current session
	uint32_t AllocdReqResp;				// allocation size for pReqResp
	uint8_t *pReqResp;					// allocation for NumInstances of requests and associated responses (tsReqResp's)
	tsTxdRxd TxdRxd;					// holding low level send/receive buffers + connected socket
} tsBKSRegSessionEx;


typedef struct TAG_sBKSType
{
	teBKSPType BKSPType;					// registering this service type
	uint32_t ProviderVersion;					// service provider version 
	uint32_t MaxServiceInsts;					// limited to support a maximum of this many service instances
	uint32_t MaxQuerySeqLen;					// accepting query sequences up to this length
	uint32_t MaxTargSeqLen;					// accepting target sequences up to this length
	uint32_t MaxReqPayloadSize;				// request payloads from the service requester, including framing, can be up to this size (UINT8s),
	uint32_t MaxRespPayloadSize;				// response payloads from the service provider to requester, including framing, can be up to this  size (UINT8s)

	double ScalePayload;					// ScaledPayload = (MaxReqPayloadSize + MaxTargSeqLen) * ScalePayload
	double ScaleTargLen;					// ScaledTargLen = MaxTargSeqLen * ScaleTargLen
	double ScaleQueryLen;					// ScaledQueryLen = MaxTargSeqLen * ScaleQueryLen
	double ScaleScaledTargQuery;			// ScaledScaledTargQuery = (ScaledTargLen * ScaledQueryLen) * ScaleScaledTargQuery;
	double AvailResources;					// OfferedInstances = AvailResources / (ScaledScaledTargQuery + ScalePayload);
	uint32_t(*pfnCostNumInstances)(teBKSPType BKSPType,			// costing function for determining number of service instances to offer for this service type
										uint32_t MaxServiceInsts,			// limited to support a maximum of this many service instances
										uint32_t MaxQuerySeqLen,			// accepted query sequences can be up to this length
										uint32_t MaxTargSeqLen,			// accepted target sequences can be up to this length
										uint32_t MaxReqPayloadSize,		// request payloads from the service requester, including framing, can be up to this size (UINT8s),
										uint32_t MaxRespPayloadSize,		// response payloads from the service provider to requester, including framing, can be up to this  size (UINT8s)
										double ScalePayload,			// ScaledPayload = (MaxReqPayloadSize + MaxTargSeqLen) * ScalePayload
										double ScaleTargLen,			// ScaledTargLen = MaxTargSeqLen * ScaleTargLen
										double ScaleQueryLen,			// ScaledQueryLen = MaxTargSeqLen * ScaleQueryLen
										double ScaleScaledTargQuery,	// ScaledScaledTargQuery = (ScaledTargLen * ScaledQueryLen) * ScaleScaledTargQuery
										double AvailResources);			// OfferedInstances = AvailResources / (ScaledScaledTargQuery + ScalePayload)
} tsBKSType;

typedef struct TAG_sWorkerInstance {
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
	uint8_t *pReqData;				// malloc'd (cMaxReqDataSize) to hold request data
	uint8_t *pParamData;              // malloc'd (cMaxReqParamSize) to hold any paramertisations
	uint8_t *pRespData;               // malloc'd (cMaxRespDataSize) to hold response data
	char *pszBuffer;                // malloc'd )(cMaxMFABuffSize) to hold any textual alignment sequences
} tsWorkerInstance;

typedef struct TAG_sClassInstance {
	uint64_t ClassInstanceID;				// monotonically incremented 
	time_t LastAccessed;				// when this class instance was last accessed, used to expire unused class instances
	CSSW *pClass;						// if non-null then pts to instantiated class instance	
	} tsClassInstance;

#pragma pack()


class CBKSProvider
{

	int m_MaxConnWait;									// wait for at most this many minutes for connection
	int m_MaxServInsts;									// max number of service instances supported
	int m_MaxServMemGB;									// max allocatable memory available for all service instances
	bool m_bCreatedMutexes;								// set true after serialisation locks/mutexes initialised
	bool m_bTermConnectionReq;							// connection termination has been requested and is either being progressed or termination has completed
	bool m_bNotifiedReqs;								// set TRUE if control message has been sent to control socket 2
	uint32_t m_NumInitTypes;								// number of service types initialised in m_BKSTypes[]
	tsBKSType m_BKSTypes[eBKSPTPlaceHolder - 1];		// entry for each potentially supported service type indexed by Type-1

	tsBKSRegSessionEx m_BKSConnection;							// server connection

	uint32_t m_ReqNumWorkerInsts;							// number of worker instance threads to start
#ifdef WIN32
	alignas(4) volatile uint32_t m_NumPendCpltd;			// count of threads waiting to gain lock in JobResponse() to notify job completed response 
	alignas(4) volatile uint32_t  m_NumWorkerInsts;				// number of worker instance threads actually started
	alignas(4) volatile uint32_t  m_NumSWAlignReqs;				// number of SW alignments requested
	alignas(4) volatile uint32_t m_TermAllThreads;                // will be set to 1 if all worker threads are to terminate
#else
	__attribute__((aligned(4))) volatile uint32_t m_NumPendCpltd;			// count of threads waiting to gain lock in JobResponse() to notify job completed response 
	__attribute__((aligned(4))) volatile uint32_t  m_NumWorkerInsts;				// number of worker instance threads actually started
	__attribute__((aligned(4)))  volatile uint32_t  m_NumSWAlignReqs;				// number of SW alignments requested
	__attribute__((aligned(4))) volatile uint32_t m_TermAllThreads;                  // will be set to 1 if all worker threads are to terminate
#endif
	tsWorkerInstance m_WorkerInstances[cMaxServiceInsts];	// to hold all worker instance thread parameters

	uint32_t m_MaxClassInsts;									// at most this many class instances can be instantiated
	uint32_t m_NumClassInsts;									// currently this many class instances have been instantiated in m_ClassInstances[];
	uint64_t m_HiClassInstanceID;								// highest class instance identifier thus far allocated							
	tsClassInstance m_ClassInstances[cMaxClassInsts];		// all possible class instances

	tsTxdRxd m_Ctrl[2];									// Ctrl[0] written to by threads needing to signal select() processing thread, select() processing thread monitors m_Ctrl[1]

	teBSFrsltCodes										// cBSFSuccess if no errors and registration process is continued, cBSFSocketErr if any errors and connection has been terminated, eBSFerrMem if unable to allocate memory
		ProcessSessEstab(bool bCpltdWrite);				// false if frame received, true if frame sent

	int													// number of received requests allocated to service instances
		HandleServiceRequests(void);					// handler for received requests for allocation of a service instance to process the service request

	int                                         // number of responses assembled into m_BKSConnection.TxdRxd.pTxdBuff ready to be sent back to requester
		HandleServiceResponses(void);			// locate those service instances with responses ready to be sent to the requester and assemble responses into m_BKSConnection.TxdRxd.pTxdBuf


	int	TerminateConnection(bool bFreeMem,		  // true if any allocated memory to be free'd
							bool bCloseCtrl,	  // close control socket pair
							bool bResetInitTypes); // reset all service types

	int	CreateMutexes(void);				// create mutexes used in access serialisations
	int	DeleteMutexes(void);				// delete mutexes used in access serialisations

	void AcquireLock(bool bExclusive);	// lock for serialised access by multiple concurrent reader threads (bExclusive == false), or serialised access by single thread (bExclusive == true)
	void ReleaseLock(bool bExclusive);	// release serialised access lock

#ifdef _WIN32
	SRWLOCK m_hRWLock;					// serialising multiple reader access but single writer
#else
	pthread_rwlock_t m_hRWLock;			// serialising multiple reader access but single writer
#endif
#ifdef WIN32
	alignas(4) volatile unsigned int m_CASSerialise; // used with synchronous compare and swap (CAS) for serialising access 
#else
	__attribute__((aligned(4))) volatile unsigned int m_CASSerialise; // used with synchronous compare and swap (CAS) for serialising access
#endif
	void AcquireCASSerialise(bool bPriority = false);	// if bPriority true then backoff time is reduced relative to if bPriority is false, increasing the probability of acquiring the serialisation lock if there is contention
	void ReleaseCASSerialise(void);


	int				// 1 if timed out attempting to connect to server, 0 if terminated because requested to terminate, -1 if socket level errors
		ConnectServer(int MaxConnWait,							 // wait for at most this many minutes for connection
						const char* pszHost, char *pszService);  // connect to server at pszHost:pszService

	bool InitialiseCtrlSocks(void); // initialise the control sockets in m_Ctrl[]
																// InitialiseConnect 
	bool					// true if connecting socket initialised
		InitialiseConnect(int MaxConnWait,			// wait for at most this many minutes for connection
							const char* pszHost,	// connecting to this server host/IP address; NULL to use first INET IP local to this machine
							char *pszService);		// connecting to this server service/port; NULL to use default port 

	bool RxData(tsTxdRxd *pRxd);
	bool TxData(tsTxdRxd *pTxd);

	bool NotifyCtrl(uint8_t Msg = 0);			// notify via control sockets that, default, there is at least 1 response available to be sent to service requester
	bool SendCtrlMsg(int Len,				// number of control message bytes to send ptd at by pCtrlMsg, can be 0
								uint8_t *pCtrlMsg);	// pCtrlMsg pts to control message bytes
	int								// received message is this length, 0 if no outstanding messages, < 0 if errors
			RcvCtrl(int BuffLen,			// available buffer into which the control message can be copied 
		uint8_t *pBuff);			// copy control message into this buffer
	bool ProcessCtrlMsg(int MsgLen,uint8_t *pMsg);			// process a message received by control socket m_Ctrl[1] - note: currently these messages are simply discarded

	bool ShutdownConnection(socket_t *pSocket);


	int			// returns 0 if no sockets to be monitored with select(), on windows the total number of monitored sockets, on linux the highest socket file descriptor plus 1
		SetupFDSets(fd_set& ReadFDs,			// select() read available socket descriptor set  
					fd_set& WriteFDs,			// select() write accepted socket descriptor set
					fd_set& ExceptFDs);			// select() exceptions descriptor set
#ifdef WIN32
	const char *WSAGetLastErrorMessage(const char* pcMessagePrefix,int nErrorID = 0);
#endif
	int	Reset(void);

	int Initialise(int MaxConnWait,				// wait for at most this many minutes for connection
					int MaxServInsts,			// max number of service instances supported
						int MaxClassInsts,		// max number of class instances supported (may differ from MaxServInsts)
						int MaxServMemGB);		// max allocatable memory available for all service instances

	int StartWorkerThreads(uint32_t NumInstances,	// this number of worker threads required
						   uint8_t BKSPType);		// workers are providing this service type

	int TerminateWorkerThreads(void);			// terminate all worker threads

	int											// marshaled parameter required this many bytes
		MarshalResp(uint8_t *pInto,				// marshal into this list
				teRMIParamType Type,				// parameter type
				void *pValue,					// parameter value
				uint32_t ValLen);					// length of parameter ptd to by pValue, only used if parameter type is pUint8

	int
		UnmarshalReq(uint32_t DataLen,
					uint8_t *pFrom,		// unmarshal from this marshalled parameter list
					void *pValue);

public:
	CBKSProvider();
	~CBKSProvider();

	int ProcWorkerThread(tsWorkerInstance *pThreadPar);  // worker thread startup entry for processing requests

	int					// returns total number of registered service types or teBSFrsltCodes error code if any parameterisation errors or already registered type
		RegServiceType(teBKSPType BKSPType = eBKSPTSmithWaterman,			// registering this service type
					   uint32_t ProviderVersion = cServiceProviderVersion,	// service provider version 
					   uint32_t MaxServiceInsts = cMaxServiceProviderInsts,	// limited to support a maximum of this many service instances
					   uint32_t MaxQuerySeqLen = cMaxSWQuerySeqLen,			// accepted query sequences can be up to this length
					   uint32_t MaxTargSeqLen = cMaxSWTargSeqLen,				// accepted target sequences can be up to this length
					   uint32_t MaxReqPayloadSize = cMaxSWReqPayloadSize,		// request payloads from the service requester, including framing, can be up to this size (UINT8s),
					   uint32_t MaxRespPayloadSize = cMaxSWRespPayloadSize,	// response payloads from the service provider to requester, including framing, can be up to this  size (UINT8s)
						double ScalePayload = 1.0,					// ScaledPayload = (MaxReqPayloadSize + MaxTargSeqLen) * ScalePayload
						double ScaleTargLen = 1.0,					// ScaledTargLen = MaxTargSeqLen * ScaleTargLen
						double ScaleQueryLen = 1.0,					// ScaledQueryLen = MaxTargSeqLen * ScaleQueryLen
						double ScaleScaledTargQuery = 0.05,			// ScaledScaledTargQuery = (ScaledTargLen * ScaledQueryLen) * ScaleScaledTargQuery;
						double AvailResources = 256000000000.0,				// OfferedInstances = AvailResources / (ScaledScaledTargQuery + ScalePayload);
					   uint32_t  (*pfnCostNumInstances)(teBKSPType BKSPType,				// determine number of service instances to offer for this service type
													  uint32_t MaxServiceInsts,			// limited to support a maximum of this many service instances
													  uint32_t MaxQuerySeqLen,			// accepted query sequences can be up to this length
													  uint32_t MaxTargSeqLen,			// accepted target sequences can be up to this length
													  uint32_t MaxReqPayloadSize,		// request payloads from the service requester, including framing, can be up to this size (UINT8s),
													  uint32_t MaxRespPayloadSize,		// response payloads from the service provider to requester, including framing, can be up to this  size (UINT8s)
													  double ScalePayload,					// ScaledPayload = (MaxReqPayloadSize + MaxTargSeqLen) * ScalePayload
													  double ScaleTargLen,					// ScaledTargLen = MaxTargSeqLen * ScaleTargLen
													  double ScaleQueryLen,				// ScaledQueryLen = MaxTargSeqLen * ScaleQueryLen
													  double ScaleScaledTargQuery,			// ScaledScaledTargQuery = (ScaledTargLen * ScaledQueryLen) * ScaleScaledTargQuery
													  double AvailResources) = NULL);				// OfferedInstances = AvailResources / (ScaledScaledTargQuery + ScalePayload)

	
	// after all service types offered have been registered then Process() can be invoked													
	int 
		Process(int MaxConnWait,						// wait for at most this many minutes for connection 
				int MaxServInsts,						// max number of service instances supported
				int MaxClassInsts,						// max number of class instances supported (may differ from MaxServInsts)
				int MaxServMemGB,						// max allocatable memory available for all service instances 
				const char* pszHost = NULL,				// connect to this host/IP address; NULL to use first INET IP local to this machine
				char *pszService = NULL);				// connect to this service/port; NULL to use default port 


	int			// 0 no service request job available, 1 job details are being returned, -1 job is available but MaxParamSize < required, -2 job available but MaxRequestData < required, -3 if session being terminated
			GetJobToProcess(uint32_t *pInstanceID,		// service instance identifier, the JobResponse() must use this identifier
						uint64_t *pClassInstanceID,	// returned class instance to apply job to
						uint32_t *pClassMethodID,		// returned class method to apply on class instance
					    uint32_t *pMaxParamsSize,		// on input, max sized parameter block accepted, on return then the actual size of the parameter block
						uint8_t *pParams,				// returned parameter block
						uint32_t *pMaxRequestData,	// on input the max sized request data block accepted, on return the actual size of the request data block
						uint8_t *pRequestData);		// returned request data block


	int			// 0 if response accepted, -1 if job does not exist or parameterisation errors, -3 if session terminating 
			JobResponse(int32_t InstanceID,			// service instance identifier returned by GetJobToProcess
					uint64_t ClassInstanceID,		// response is for this class instance
				    uint32_t ProcRslt,			// service request processing result
					uint32_t ResponseSize,		// response data block size
					uint8_t *pResponseData);		// response data


	tsClassInstance *LocateClassInstance(uint64_t ClassInstanceID);
	tsClassInstance *AllocClassInstance(void); // allocate a new tsClassInstance and initialise with ClassInstanceID 
	bool FreeClassInstance(uint64_t ClassInstanceID); // free a previously allocated class instance

};


