//===================================================
//$Id$
//===================================================

#ifndef EXTERNALFUNC 
#define EXTERNALFUNC 

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
       exit(0);
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

  void SetFunc(int fnum,T newfunc) {
    if ((func!=NULL)&&(nfuncs>fnum))
      func[fnum]=newfunc;
    else if (func==NULL) {
      printf("CExternalFunction::SetFunc: func=NULL\n");
      exit(0);
    }
    else {
      printf("CExternalFunction::SetFunc: func buffer overflow\n");
      exit(0);
    }
  };


  void SetFunc(T newfunc) {
    static int nattachedfuncs=0;

    if ((func!=NULL)&&(nfuncs>nattachedfuncs))
      func[nattachedfuncs++]=newfunc;
    else if (func==NULL) {
      printf("CExternalFunction::SetFunc: func=NULL\n");
      exit(0);
    }
    else {
      printf("CExternalFunction::SetFunc: func buffer overflow\n");
      exit(0);
    }
  };

  T GetFunc(int n) {
    T res=NULL;

    if ((func!=NULL)&&(n<nfuncs)) if (func[n]!=NULL) res=func[n];

    if (res==NULL) {
       printf("Cbc::CExternalFaceBCFunction::GetFunc: Error: cannot find function for n=%i\n",n);
       exit(0);
    }

    return res;
  };

  CExternalFunction() {nfuncs=0,func=NULL;};
};

#endif
