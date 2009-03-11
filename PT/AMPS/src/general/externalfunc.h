//===================================================
//$Id$
//===================================================

#ifndef EXTERNALFUNC 
#define EXTERNALFUNC 

#include "specfunc.h"

template<class T>
class CExternalFunction {
private:
  int nfuncs;

public:
  T *func;

  void init(int n) {
    if (nfuncs==0) {
       nfuncs=n;
       func=new T[n];

       for (int i=0;i<n;i++) func[i]=NULL;
    }
    else {
       printf("CExternalFunction: nfuncs!=0 <=> repeated initialization\n");
       exit(__LINE__,__FILE__);
    }
  };

  void reinit(int n) {
    if (n<nfuncs) return;
    else if (nfuncs==0) init(n);
    else {
      T* oldfunc=func;
      int i,oldnfuncs=nfuncs;

      nfuncs=n,func=new T[n]; 
      for (i=0;i<n;i++) func[i]=(i<oldnfuncs) ? oldfunc[i] : NULL;

      delete [] oldfunc;
    }
  };

  int size() {
    return nfuncs;
  };

  void SetFunc(int fnum,T newfunc) {
    if ((func!=NULL)&&(nfuncs>fnum)) func[fnum]=newfunc;
    else { 
      reinit(fnum+1);
      func[fnum]=newfunc;
    }
  };


  void SetFunc(T newfunc) {
    static int nattachedfuncs=0;

    if ((func!=NULL)&&(nfuncs>nattachedfuncs)) func[nattachedfuncs++]=newfunc;
    else {
      reinit(nattachedfuncs+1);
      func[nattachedfuncs++]=newfunc;
    }
  };

  T GetFunc(int n) {
    if ((func==NULL)||(n>=nfuncs)) {
       char str[200];
       sprintf(str,"template<class T> class CExternalFunction::GetFunc: Error: cannot find function for n=%i\n",n);
       PrintErrorLog(str);

       exit(__LINE__,__FILE__);
    }

    return func[n];
  };

  T operator[] (int n) {
    return (n<nfuncs) ? func[n] : NULL;
  };

  bool initialized(int n) {
    if ((func!=NULL)&&(n<nfuncs)) if (func[n]!=NULL) return true;

    return false;
  };

  CExternalFunction() {nfuncs=0,func=NULL;};
};

#endif
