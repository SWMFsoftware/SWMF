/*
 * VIRTIS-M.h
 *
 *  Created on: Jan 20, 2016
 *      Author: vtenishe
 */

//geomentry related to VIRTIS-M calculations

#ifndef SRCCG_POSTPROCESS_VIRTIS_M_H_
#define SRCCG_POSTPROCESS_VIRTIS_M_H_

#include "SpiceUsr.h"
#include "specfunc.h"
#include "global.h"
#include "constants.h"
#include "meshAMRcutcell.h"
#include "PostProcess3D.h"

class cVirtisM {
public:

  static const double AngularFieldOfView;
  static const int nFieldOfViewPixels;
  static const double dAnglePixel;



  //frame of reference e0 -> toward the nucleus, e1 -> along the direction of flight, e2 -> the right handed frame of reference
  //xRosetta -> location of the spacecraft in the frame of reference related to the nucleus
  double e0[3],e1[3],e2[3],xRosetta[3];

  void SetFrameAxis(SpiceDouble et);
  void GetNucleusPixelCoordinate(int &e0PixelNucleusCoordinate,int &e1PixelNucleusCoordinate,SpiceDouble et);


  class cBlockNucleus {
  public:
    cVirtisM* Virtis;

    //the range of the VIRTIS map for column integraion
    int nxVirtisBlockPixelRange,nyVirtisBlockPixelRange;
    int AdditionalBlockedPixels;

    //coordinates of the nucleus
    int iCenter,jCenter;
    int  **VirtisMask,**NucleusMask;
    double ***InstrumentPointing;

    void SetPixelLimits(int dxPix,int dyPix); //set the block limit
    void SetBlock(SpiceDouble et,int nNucleusSurfaceFaces,CutCell::cTriangleFace *NucleusSurfaceFaces); //set blocking of the nucleus

    //get column integrals
    void GetColumnIntegralMap(const char *fname,cPostProcess3D::cColumnIntegral::cColumnIntegrationSet* IntegrationSet,cPostProcess3D* PostProcessor);

    cBlockNucleus() {
      nxVirtisBlockPixelRange=0,nyVirtisBlockPixelRange=0;
      AdditionalBlockedPixels=7;

      //init the mask
      int i,j,k,cnt;

      VirtisMask=new int* [cVirtisM::nFieldOfViewPixels];
      VirtisMask[0]=new int[cVirtisM::nFieldOfViewPixels*cVirtisM::nFieldOfViewPixels];

      NucleusMask=new int* [cVirtisM::nFieldOfViewPixels];
      NucleusMask[0]=new int[cVirtisM::nFieldOfViewPixels*cVirtisM::nFieldOfViewPixels];

      for (i=1;i<cVirtisM::nFieldOfViewPixels;i++) {
        VirtisMask[i]=VirtisMask[i-1]+cVirtisM::nFieldOfViewPixels;
        NucleusMask[i]=NucleusMask[i-1]+cVirtisM::nFieldOfViewPixels;
      }

      for (i=0;i<cVirtisM::nFieldOfViewPixels;i++) for (j=0;j<cVirtisM::nFieldOfViewPixels;j++) VirtisMask[i][j]=false,NucleusMask[i][j]=false;

      //init the instrument pointing matrix
      InstrumentPointing=new double **[cVirtisM::nFieldOfViewPixels];
      InstrumentPointing[0]=new double *[cVirtisM::nFieldOfViewPixels*cVirtisM::nFieldOfViewPixels];
      InstrumentPointing[0][0]=new double [3*cVirtisM::nFieldOfViewPixels*cVirtisM::nFieldOfViewPixels];

      for (i=1;i<cVirtisM::nFieldOfViewPixels;i++) InstrumentPointing[i]=InstrumentPointing[i-1]+cVirtisM::nFieldOfViewPixels;
      for (cnt=0,i=0;i<cVirtisM::nFieldOfViewPixels;i++) for (j=0;j<cVirtisM::nFieldOfViewPixels;j++) {
        InstrumentPointing[i][j]=InstrumentPointing[0][0]+cnt;
        cnt+=3;
      }

      for (cnt=0,i=0;i<cVirtisM::nFieldOfViewPixels;i++) for (j=0;j<cVirtisM::nFieldOfViewPixels;j++) for (k=0;k<3;k++) InstrumentPointing[i][j][k]=0.0;
    }
  };

  cBlockNucleus BlockNucleus;

  cVirtisM() {
    BlockNucleus.Virtis=this;

  }
};



#endif /* SRCCG_POSTPROCESS_VIRTIS_M_H_ */
