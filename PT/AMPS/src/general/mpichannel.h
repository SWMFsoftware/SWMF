//====================================================================
//$Id$
//====================================================================

#ifndef MPI_CHANNEL
#define MPI_CHANNEL

class CMPI_channel {
private:
  long int max_MPIbuffer_size;

  char *sendBuffer,**recvBuffer;
  long int sendptr,*recvptr;
  
  long int *RecvDataLength;
  int sendThread;

public:

  CMPI_channel(long int buffersize) {
    max_MPIbuffer_size=buffersize;
    sendBuffer=NULL,sendptr=0,recvptr=NULL;
    sendThread=-1;

    try {
      recvBuffer=new char* [TotalThreadsNumber];
      recvptr=new long int[TotalThreadsNumber];
      RecvDataLength=new long int[TotalThreadsNumber];
    }
    catch (bad_alloc) {
      printf("Error: CMPI_channel::CMPI_channel(long) cannot allocate buffer\n");
      exit(__LINE__,__FILE__);
    }

    for (int thread=0;thread<TotalThreadsNumber;thread++) recvBuffer[thread]=NULL,recvptr[thread]=0,RecvDataLength[thread]=0; 
  };

 ~CMPI_channel() {
    int thread;

    if (sendBuffer!=NULL) delete [] sendBuffer;

    if (recvBuffer!=NULL) {
      for (thread=0;thread<TotalThreadsNumber;thread++) if (recvBuffer[thread]!=NULL) delete [] recvBuffer[thread];
      delete [] recvBuffer;
      delete [] recvptr;
      delete [] RecvDataLength;
    }
  }; 

  void openSend(int send) {
    if (sendBuffer!=NULL) {
      printf("CMPI_channel::openSend reinitialization of the send buffer\n");
      exit(__LINE__,__FILE__);
    }  

    try {
      sendBuffer=new char[max_MPIbuffer_size];
    }
    catch (bad_alloc) {
      printf("Error: CMPI_channel::openSend cannot allocate send buffer, need %ld bytes\n", max_MPIbuffer_size);
      exit(__LINE__,__FILE__);
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
    if (recvBuffer[thread]!=NULL) {  
      printf("CMPI_channel::openRecv reinitialization of the recv buffer\n");
      exit(__LINE__,__FILE__);
    }

    try {
      recvBuffer[thread]=new char[max_MPIbuffer_size];
    }
    catch (bad_alloc) {
      printf("Error: CMPI_channel::openRecv cannot allocate recv buffer, need %ld bytes\n", max_MPIbuffer_size);
      exit(__LINE__,__FILE__);
    }

    recvptr[thread]=0;
  };

  void openRecvAll() {
    int thread;

    for (thread=0;thread<TotalThreadsNumber;thread++) if (thread!=ThisThread) openRecv(thread);
  };
 

  void closeRecv(int thread) {
    if (recvptr[thread]!=RecvDataLength[thread]) {
      printf("Error: closeRecv: (recvptr[thread]!=RecvDataLength[thread]) for thread=%i\n",thread);
      exit(__LINE__,__FILE__);
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
