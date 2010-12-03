//====================================================================
//$$
//====================================================================

#ifndef MPI_CHANNEL
#define MPI_CHANNEL

using namespace std;

#define MPI_CHANNEL_DEFAULT_BUFFER_SIZE 100000

class CMPI_channel {
private:
  long int max_MPIbuffer_size;

  char *sendBuffer,**recvBuffer;
  long int sendptr,*recvptr;
  
  long int *RecvDataLength;
  int sendThread,TotalThreadsNumber,ThisThread;

  void exit(long int nline,const char* msg=NULL) {
    char str[2000];
    FILE* errorlog=fopen("error.log","a+");

    time_t TimeValue=time(0);
    tm *ct=localtime(&TimeValue);

    if (msg==NULL) sprintf(str," exit: line=%ld, file=%s\n",nline,__FILE__);
    else sprintf(str," exit: line=%ld, file=%s, message=%s\n",nline,__FILE__,msg);

    fprintf(errorlog,"Thread=%i: (%i/%i %i:%i:%i)\n",ThisThread,ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);
    fprintf(errorlog,"%s\n\n",msg);

    fclose(errorlog);
    std::exit(0);
  }

public:

  void init(long int buffersize) {
    max_MPIbuffer_size=buffersize;
    sendBuffer=NULL,sendptr=0,recvptr=NULL;
    sendThread=-1;

    #ifdef MPI_ON
    MPI_Comm_rank(MPI_COMM_WORLD,&ThisThread);
    MPI_Comm_size(MPI_COMM_WORLD,&TotalThreadsNumber);
    #else
    ThisThread=0,TotalThreadsNumber=1;
    #endif

    try {
      recvBuffer=new char* [TotalThreadsNumber];
      recvptr=new long int[TotalThreadsNumber];
      RecvDataLength=new long int[TotalThreadsNumber];
    }
    catch (bad_alloc) {exit(__LINE__,"Error: CMPI_channel::CMPI_channel(long) cannot allocate buffer");} 
 
    for (int thread=0;thread<TotalThreadsNumber;thread++) recvBuffer[thread]=NULL,recvptr[thread]=0,RecvDataLength[thread]=0; 
  };

  CMPI_channel(long int buffersize) {
    init(buffersize);
  }

  CMPI_channel() {
    max_MPIbuffer_size=0;
    sendBuffer=NULL,sendptr=0,recvptr=NULL;
    sendThread=-1;

    ThisThread=0,TotalThreadsNumber=0;

    recvBuffer=NULL;
    recvptr=NULL;
    RecvDataLength=NULL;
  };

  void remove() {
    int thread;

    if (sendBuffer!=NULL) {
      delete [] sendBuffer;
      sendBuffer=NULL;
    }

    if (recvBuffer!=NULL) {
      for (thread=0;thread<TotalThreadsNumber;thread++) if (recvBuffer[thread]!=NULL) delete [] recvBuffer[thread];
      delete [] recvBuffer;
      delete [] recvptr;
      delete [] RecvDataLength;

      recvBuffer=NULL,recvptr=NULL,RecvDataLength=NULL;
    }
  }; 

  ~CMPI_channel() {
     remove();
  } 


  void openSend(int send) {
    if (max_MPIbuffer_size==0) init(MPI_CHANNEL_DEFAULT_BUFFER_SIZE);
    if (sendBuffer!=NULL) exit(__LINE__,"CMPI_channel::openSend reinitialization of the send buffer"); 

    try {
      sendBuffer=new char[max_MPIbuffer_size];
    }
    catch (bad_alloc) {
      char msg[200]; 

      sprintf(msg,"Error: CMPI_channel::openSend cannot allocate send buffer, need %ld bytes\n", max_MPIbuffer_size);
      exit(__LINE__,msg);
    }

    sendptr=0,sendThread=send; 
  };

  void flush() {
    if (sendptr!=0) {

#ifdef MPI_ON
      MPI_Send(&sendptr,1,MPI_LONG,sendThread,0,MPI_COMM_WORLD);
      MPI_Send(sendBuffer,sendptr,MPI_CHAR,sendThread,0,MPI_COMM_WORLD);
#endif

      sendptr=0;
    }
  };


  void closeSend() {
    if (sendptr!=0) {
      flush();

      delete [] sendBuffer;

      sendBuffer=NULL;
      sendThread=-1;
    }
  };

  template<class T> inline void send(T data) {
    register long int i,length=sizeof(T); 
    register char *ptr;

    if (sendptr+length>=max_MPIbuffer_size) {

#ifdef MPI_ON
      MPI_Send(&sendptr,1,MPI_LONG,sendThread,0,MPI_COMM_WORLD); 
      MPI_Send(sendBuffer,sendptr,MPI_CHAR,sendThread,0,MPI_COMM_WORLD);
#endif

      sendptr=0;
    } 
     
    for (i=0,ptr=(char*)&data;i<length;i++,ptr++) sendBuffer[sendptr++]=*ptr; 
  }


  template<class T> void inline send(T* data,long int nsend) {
    register long int i;
    register T* dataptr; 

    for (i=0,dataptr=data;i<nsend;i++,dataptr++) send(*dataptr);
  }

  void openRecv(int thread) { 
    if (max_MPIbuffer_size==0) init(MPI_CHANNEL_DEFAULT_BUFFER_SIZE);
    if (recvBuffer[thread]!=NULL) exit(__LINE__,"CMPI_channel::openRecv reinitialization of the recv buffer");   

    try {
      recvBuffer[thread]=new char[max_MPIbuffer_size];
    }
    catch (bad_alloc) {
      char msg[200];

      sprintf(msg,"Error: CMPI_channel::openRecv cannot allocate recv buffer, need %ld bytes\n", max_MPIbuffer_size);
      exit(__LINE__,msg);
    }

    recvptr[thread]=0;
  };

  void openRecvAll() {
    int thread;

    for (thread=0;thread<TotalThreadsNumber;thread++) if (thread!=ThisThread) openRecv(thread);
  };
 

  void closeRecv(int thread) {
    if (recvptr[thread]!=RecvDataLength[thread]) {
      char msg[200];

      sprintf(msg,"Error: closeRecv: (recvptr[thread]!=RecvDataLength[thread]) for thread=%i\n",thread);
      exit(__LINE__,msg);
    }

    delete [] recvBuffer[thread];
    recvBuffer[thread]=NULL,recvptr[thread]=0,RecvDataLength[thread]=0; 
  };

  void closeRecvAll() {
    int thread;

    for (thread=0;thread<TotalThreadsNumber;thread++) if (thread!=ThisThread) closeRecv(thread);
  };

  template<class T> inline void recv(T& data,int thread) {
    register long int i,length=sizeof(T);
    register T t;
    register char *ptr;
   
    if (recvptr[thread]==RecvDataLength[thread]) {

#ifdef MPI_ON
      MPI_Status status;

      MPI_Recv(RecvDataLength+thread,1,MPI_LONG,thread,0,MPI_COMM_WORLD,&status);
      MPI_Recv(recvBuffer[thread],RecvDataLength[thread],MPI_CHAR,thread,0,MPI_COMM_WORLD,&status);
#endif

      recvptr[thread]=0;
    }

    for (i=0,ptr=(char*)&t;i<length;i++,ptr++) *ptr=recvBuffer[thread][recvptr[thread]++];    
    data=t;  
  }

  template<class T> inline void recv(T* data,long int nrecv,int thread) {
    register long int i;
    register T *dataptr;

    for (i=0,dataptr=data;i<nrecv;i++,dataptr++) recv(*dataptr,thread);
  } 

  template<class T> inline T recv(int thread) {
    T data;
   
    recv(data,thread); 
    return data;
  }
};

#endif 
