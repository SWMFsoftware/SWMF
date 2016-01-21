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

  class cBlockNucleus {
  public:
    cVirtisM* Virtis;

    class cMaskPixel {
    public:
      int **mask;

      void SetMask(int nx,int ny);
    };

    cMaskPixel VistisMask,NucleusMask;

    static const double AngularFieldOfView;
    static const int nFieldOfViewPixels;
    static const double dAnglePixel;

    //the range of the VIRTIS map for column integraion
    int nxVirtisBlockPixelRange,nyVirtisBlockPixelRange;
    int AdditionalBlockedPixels;

    //frame of reference e0 -> toward the nucleus, e1 -> along the direction of flight, e2 -> the right handed frame of reference
    //xRosetta -> location of the spacecraft in the frame of reference related to the nucleus
    double e0[3],e1[3],e2[3],xRosetta[3];

    void SetPixelLimits(int dxPix,int dyPix); //set the block limit
    void SetBlock(SpiceDouble et,int nNucleusSurfaceFaces,CutCell::cTriangleFace *NucleusSurfaceFaces); //set blocking of the nucleus

    //get column integrals
    void GetColumnIntegralMap(const char *fname,cPostProcess3D::cColumnIntegral::cColumnIntegrationSet* IntegrationSet,cPostProcess3D* PostProcessor);

    cBlockNucleus() {
      VistisMask.SetMask(nFieldOfViewPixels,nFieldOfViewPixels);
      NucleusMask.SetMask(nFieldOfViewPixels,nFieldOfViewPixels);
      nxVirtisBlockPixelRange=0,nyVirtisBlockPixelRange=0;
      AdditionalBlockedPixels=7;
    }
  };

  cBlockNucleus BlockNucleus;
};



#endif /* SRCCG_POSTPROCESS_VIRTIS_M_H_ */
