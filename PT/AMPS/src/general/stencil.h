//functionality used for constructing of stencils

/*
 * stencil.h
 *
 *  Created on: Dec 9, 2017
 *      Author: vtenishe
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <list>


#ifndef _GENERAL_STENCIL_H_
#define _GENERAL_STENCIL_H_

class cStencil {
public:

  struct cStencilElement {
    int i,j,k;
    double a; //coefficient used in the stencil
  };

  list<cStencilElement> StencilData;
  cStencilElement *Stencil;
  int StencilLength,AllocatedStencilLength;

  cStencil() {
    Stencil=NULL;
    StencilLength=0,AllocatedStencilLength=0;
  }

  ~cStencil() {
    if (Stencil!=NULL) delete [] Stencil;
    StencilData.clear();
  }

  //remove elements of the list that have the same combination of i,j, and k
  void Simplify() {
    int i,k,j;
    list<cStencilElement>::iterator l0,l1;

    for (l0=StencilData.begin();l0!=StencilData.end();l0++) {
      list<list<cStencilElement>::iterator> CleaningList;

      i=l0->i,j=l0->j,k=l0->k;
      l1=l0;

      for (l1++;l1!=StencilData.end();l1++) {
        if ((i==l1->i)&&(j==l1->j)&&(k==l1->k)) {
          //found another element in the list that point to the point with the same combination of i,j, and k
          //combine l0 and l1
          l0->a+=l1->a;

          CleaningList.push_back(l1);
        }
      }

      for (list<list<cStencilElement>::iterator>::iterator t=CleaningList.begin();t!=CleaningList.end();t++) StencilData.erase(*t);

      CleaningList.clear();
    }
  }

  //in order to speed up usage of the stencil the data are saved in an aray rather than probided as a list
  cStencilElement* GetStencil() {
    Simplify();
    StencilLength=StencilData.size();

    //allocate the data buffer
    if (Stencil!=NULL) {
      if (AllocatedStencilLength<=StencilLength) {
        delete [] Stencil;
        Stencil=NULL;
      }
    }

    if (Stencil==NULL) {
      AllocatedStencilLength=(int)(1.2*StencilLength);
      Stencil=new cStencilElement[AllocatedStencilLength];
    }

    //copy the contant of the list into the array
    int i;
    std::list<cStencilElement>::iterator l0;

    for (i=0,l0=StencilData.begin();l0!=StencilData.end();i++,l0++) Stencil[i]=*l0;

    return this->Stencil;
  }

  //multiply the entire stencil by a constatnt
  friend cStencil& operator *= (cStencil &v1,const double t) {
    for (std::list<cStencilElement>::iterator l0=v1.StencilData.begin();l0!=v1.StencilData.end();l0++) l0->a*=t;
    return v1;
  };

  //add a new element to the stencil
  void add(double a,int i,int j,int k) {
    cStencilElement NewElement;

    NewElement.a=a;
    NewElement.i=i;
    NewElement.j=j;
    NewElement.k=k;

    StencilData.push_back(NewElement);
    Simplify();
  }

  //shift the entire stencil
  void shift(int di,int dj,int dk) {
    for (std::list<cStencilElement>::iterator l0=StencilData.begin();l0!=StencilData.end();l0++) {
      l0->i+=di;
      l0->j+=dj;
      l0->k+=dk;
    }
  }

  //copy the stencil
  cStencil& operator = (const cStencil& v) {
    for (list<cStencilElement>::const_iterator l0=v.StencilData.begin();l0!=v.StencilData.end();l0++) StencilData.push_back(*l0);
    return *this;
  };

  //add and substract another stencil
  friend cStencil& operator += (cStencil &v1,const cStencil &v2) {
    cStencilElement NewElement;

    for (list<cStencilElement>::const_iterator l0=v2.StencilData.begin();l0!=v2.StencilData.end();l0++) {
      NewElement=*l0;
      v1.StencilData.push_back(NewElement);
    }

    return v1;
  };

  friend cStencil& operator -= (cStencil &v1,const cStencil &v2) {
    cStencilElement NewElement;

    for (list<cStencilElement>::const_iterator l0=v2.StencilData.begin();l0!=v2.StencilData.end();l0++) {
      NewElement=*l0;
      NewElement.a*=-1.0;
      v1.StencilData.push_back(NewElement);
    }

    return v1;
  };

  void AddShifled(cStencil& v,int di,int dj,int dk,double c=1.0) {
    cStencilElement NewElement;

    for (std::list<cStencilElement>::iterator l0=v.StencilData.begin();l0!=v.StencilData.end();l0++) {
      NewElement=*l0;

      NewElement.a*=c;
      NewElement.i+=di,NewElement.j+=dj,NewElement.k+=dk;

      StencilData.push_back(NewElement);
    }
  }

  void SubstractShifted(cStencil& v,int di,int dj,int dk,double c=1.0) {
    AddShifled(v,di,dj,dk,-c);
  }

};




#endif /* _GENERAL_STENCIL_H_ */
