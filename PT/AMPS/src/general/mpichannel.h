//====================================================================
//$$
//====================================================================

#ifndef MPI_CHANNEL
#define MPI_CHANNEL

#include "global.h"

using namespace std;

#define MPI_CHANNEL_DEFAULT_BUFFER_SIZE 100000

#define _MPI_CHANNEL_MODE_BCAST_    0
#define _MPI_CHANNEL_MODE_SENDRECV_ 1

class CMPI_channel {
private:
  long int max_MPIbuffer_size;

  char *sendBuffer,**recvBuffer;
  long int sendptr,*recvptr;
  
  long int *RecvDataLength;
  int TotalThreadsNumber,BcastThread;

  int ChannelMode;

public:
  int ThisThread,sendThread;


  void exit(long int nline,const char* msg=NULL) {
    char str[2000];
    FILE* errorlog=fopen("$ERRORLOG","a+");

    time_t TimeValue=time(0);
    tm *ct=localtime(&TimeValue);

    if (msg==NULL) sprintf(str," exit: line=%ld, file=%s\n",nline,__FILE__);
    sprintf(str," exit: line=%ld, file=%s, message=%s\n",nline,__FILE__,msg);

    fprintf(errorlog,"Thread=%i: (%i/%i %i:%i:%i)\n",ThisThread,ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);
    fprintf(errorlog,"%s\n\n",msg);

    printf("$PREFIX:Thread=%i: (%i/%i %i:%i:%i)\n",ThisThread,ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);
    printf("$PREFIX:%s\n\n",msg);

    fclose(errorlog);
    ::exit(0);
  }

public:



  void RedirectSendBuffer(int newSendThread) {
    if (sendThread==-1) exit(__LINE__,"Error: need to init first");
    flush();
    sendThread=newSendThread;
  }

  void RedirectRecvBuffer(int newRecvThread,int oldRecvThread=-1) {

    if (recvBuffer[newRecvThread]!=NULL) return;

//    if (recvBuffer[newRecvThread]!=NULL) exit(__LINE__,"Error: the recv direction is already initialized");

    if (oldRecvThread==-1) { //search for an allocated buffer
      for (int i=0;i<TotalThreadsNumber;i++) if (recvBuffer[i]!=NULL) {
        oldRecvThread=i;
        break;
      }

      if (oldRecvThread==-1) exit(__LINE__,"Error: need to allocate the recv buffer first");
    }

    if ((recvBuffer[oldRecvThread]==NULL)||(recvptr[oldRecvThread]!=RecvDataLength[oldRecvThread])) exit(__LINE__,"Error: cannot redirect the recv buffer");

    recvBuffer[newRecvThread]=recvBuffer[oldRecvThread],recvptr[newRecvThread]=0,RecvDataLength[newRecvThread]=0;
    recvBuffer[oldRecvThread]=NULL,recvptr[oldRecvThread]=0,RecvDataLength[oldRecvThread]=0;
  }

  void init(long int buffersize) {
    max_MPIbuffer_size=buffersize;
    sendBuffer=NULL,sendptr=0,recvptr=NULL;
    sendThread=-1;

    #ifdef MPI_ON
    MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&ThisThread);
    MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR,&TotalThreadsNumber);
    #else
    ThisThread=0,TotalThreadsNumber=1;
    #endif

    try {
      recvBuffer=new char* [TotalThreadsNumber];
      recvptr=new long int[TotalThreadsNumber];
      RecvDataLength=new long int[TotalThreadsNumber];
    }
    catch (bad_alloc &ba) {exit(__LINE__,"Error: CMPI_channel::CMPI_channel(long) cannot allocate buffer");}
 
    for (int thread=0;thread<TotalThreadsNumber;thread++) recvBuffer[thread]=NULL,recvptr[thread]=0,RecvDataLength[thread]=0; 
  };

  CMPI_channel(long int buffersize) {
    ChannelMode=_MPI_CHANNEL_MODE_SENDRECV_,BcastThread=-1;

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

    ChannelMode=_MPI_CHANNEL_MODE_SENDRECV_,BcastThread=-1;
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
    catch (bad_alloc &ba) {
      char msg[200]; 

      sprintf(msg,"Error: CMPI_channel::openSend cannot allocate send buffer, need %ld bytes\n", max_MPIbuffer_size);
      exit(__LINE__,msg);
    }

    sendptr=0,sendThread=send; 
    ChannelMode=_MPI_CHANNEL_MODE_SENDRECV_;
  };

  void flush() {
    if (sendptr!=0) {

#ifdef MPI_ON
      if (ChannelMode==_MPI_CHANNEL_MODE_SENDRECV_) {
        MPI_Send(&sendptr,1,MPI_LONG,sendThread,0,MPI_GLOBAL_COMMUNICATOR);
        MPI_Send(sendBuffer,sendptr,MPI_CHAR,sendThread,0,MPI_GLOBAL_COMMUNICATOR);
      }
      else {
        MPI_Bcast(&sendptr,1,MPI_LONG,sendThread,MPI_GLOBAL_COMMUNICATOR);
        MPI_Bcast(sendBuffer,sendptr,MPI_CHAR,sendThread,MPI_GLOBAL_COMMUNICATOR);
      }
#endif

      sendptr=0;
    }
  };


  void closeSend() {
    if (sendBuffer!=NULL) {
      flush();

      delete [] sendBuffer;

      sendBuffer=NULL;
      sendThread=-1,sendptr=0;
    }
  };

  template<class T> inline void send(T data) {
    register long int i,length=sizeof(T); 
    register char *ptr;

    if (sendptr+length>=max_MPIbuffer_size) {

#ifdef MPI_ON
      if (ChannelMode==_MPI_CHANNEL_MODE_SENDRECV_) {
        MPI_Send(&sendptr,1,MPI_LONG,sendThread,0,MPI_GLOBAL_COMMUNICATOR);
        MPI_Send(sendBuffer,sendptr,MPI_CHAR,sendThread,0,MPI_GLOBAL_COMMUNICATOR);
      }
      else {
        MPI_Bcast(&sendptr,1,MPI_LONG,sendThread,MPI_GLOBAL_COMMUNICATOR);
        MPI_Bcast(sendBuffer,sendptr,MPI_CHAR,sendThread,MPI_GLOBAL_COMMUNICATOR);
      }
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
    catch (bad_alloc &ba) {
      char msg[200];

      sprintf(msg,"Error: CMPI_channel::openRecv cannot allocate recv buffer, need %ld bytes\n", max_MPIbuffer_size);
      exit(__LINE__,msg);
    }

    recvptr[thread]=0;
    ChannelMode=_MPI_CHANNEL_MODE_SENDRECV_;
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
   
    if (recvptr[thread]>=RecvDataLength[thread]) {

#ifdef MPI_ON
      if (ChannelMode==_MPI_CHANNEL_MODE_SENDRECV_) {
        MPI_Status status;

        MPI_Recv(RecvDataLength+thread,1,MPI_LONG,thread,0,MPI_GLOBAL_COMMUNICATOR,&status);
        MPI_Recv(recvBuffer[thread],RecvDataLength[thread],MPI_CHAR,thread,0,MPI_GLOBAL_COMMUNICATOR,&status);
      }
      else {
        MPI_Bcast(RecvDataLength+thread,1,MPI_LONG,thread,MPI_GLOBAL_COMMUNICATOR);
        MPI_Bcast(recvBuffer[thread],RecvDataLength[thread],MPI_CHAR,thread,MPI_GLOBAL_COMMUNICATOR);
      }
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

  //open/close procedures in the Bcast mode
  void openBcast(int thread) {
    if (ThisThread==thread) openSend(thread);
    else openRecv(thread);

    ChannelMode=_MPI_CHANNEL_MODE_BCAST_,BcastThread=thread;
  }

  void closeBcast() {
    if (ChannelMode!=_MPI_CHANNEL_MODE_BCAST_) exit(__LINE__,"wrong option");

    if (ThisThread==BcastThread) closeSend();
    else closeRecv(BcastThread);

    BcastThread=-1;
  }


};

#endif 
