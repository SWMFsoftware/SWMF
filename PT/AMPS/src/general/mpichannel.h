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
      exit(0);
    }

    for (int thread=0;thread<TotalThreadsNumber;thread++) recvBuffer[thread]=NULL,recvptr[thread]=0,RecvDataLength[thread]=0; 
  };

 ~CMPI_channel() {
    int thread;

    if (sendBuffer!=NULL) delete [] sendBuffer;

    if (recvBuffer!=NULL) {
      for (thread=0;thread<TotalThreadsNumber;thread++) if (recvBuffer[thread]!=NULL) delete [] recvBuffer[thread];
      delete [] recvBuffer,recvptr,RecvDataLength;
    }
  }; 

  void openSend(int send) {
    if (sendBuffer!=NULL) {
      printf("CMPI_channel::openSend reinitialization of the send buffer\n");
      exit(0);
    }  

    try {
      sendBuffer=new char[max_MPIbuffer_size];
    }
    catch (bad_alloc) {
      printf("Error: CMPI_channel::openSend cannot allocate send buffer, need %i bytes\n", max_MPIbuffer_size);
      exit(0);
    }

    sendptr=0,sendThread=send; 
  };

  void flush() {
    if (sendptr!=0) {
      MPI_Status status;
      MPI_Request request;

      MPI_Send(&sendptr,1,MPI_LONG,sendThread,0,MPI_COMM_WORLD);

      MPI_Isend(sendBuffer,sendptr,MPI_CHAR,sendThread,0,MPI_COMM_WORLD,&request);
      MPI_Wait(&request,&status);

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

  template<class T> void send(T* data,long int nsend) {
    long int i,length; 
    char *ptr;

    length=sizeof(T)*nsend;

    if (sendptr+length>=max_MPIbuffer_size) {
      MPI_Status status;
      MPI_Request request;

      MPI_Send(&sendptr,1,MPI_LONG,sendThread,0,MPI_COMM_WORLD); 

      MPI_Isend(sendBuffer,sendptr,MPI_CHAR,sendThread,0,MPI_COMM_WORLD,&request);
      MPI_Wait(&request,&status);

      sendptr=0;
    } 
     
    for (i=0,ptr=(char*)data;i<length;i++,ptr++) sendBuffer[sendptr++]=*ptr; 
  };


  template<class T> void send(T data) {
    send(&data,1);
  };

  void openRecv(int thread) { 
    if (recvBuffer[thread]!=NULL) {  
      printf("CMPI_channel::openRecv reinitialization of the recv buffer\n");
      exit(0);
    }

    try {
      recvBuffer[thread]=new char[max_MPIbuffer_size];
    }
    catch (bad_alloc) {
      printf("Error: CMPI_channel::openRecv cannot allocate recv buffer, need %i bytes\n", max_MPIbuffer_size);
      exit(0);
    }

    recvptr[thread]=0;
  };

  void openRecvAll() {
    int thread;

    for (thread=0;thread<TotalThreadsNumber;thread++) if (thread!=ThisThread) openRecv(thread);
  };
 

  void closeRecv(int thread) {
    delete [] recvBuffer[thread];
   
    recvBuffer[thread]=NULL,recvptr[thread]=0; 
  };

  void closeRecvAll() {
    int thread;

    for (thread=0;thread<TotalThreadsNumber;thread++) if (thread!=ThisThread) closeRecv(thread);
  };

  template<class T> void recv(T* data,long int nrecv,int thread) {
    long int i,length;
    char *ptr;

    length=sizeof(T)*nrecv; 

    if (recvptr[thread]==RecvDataLength[thread]) {
      MPI_Status status;
      MPI_Request request;

      MPI_Recv(RecvDataLength+thread,1,MPI_LONG,thread,0,MPI_COMM_WORLD,&status);

      MPI_Irecv(recvBuffer[thread],RecvDataLength[thread],MPI_CHAR,thread,0,MPI_COMM_WORLD,&request);
      MPI_Wait(&request,&status);

      recvptr[thread]=0;
    }

    for (i=0,ptr=(char*)data;i<length;i++,ptr++) *ptr=recvBuffer[thread][recvptr[thread]++];    
  };

  template<class T> void recv(T& data ,int thread) {
    recv(&data,1,thread);
  }; 

  template<class T> T recv(int thread) {
    T data;
   
    recv(&data,1,thread); 
    return data;
  };
};

#endif 
