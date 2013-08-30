#include "constants.h"
#include "pic.h"


#include "tvis.Kaguya.2009_02_06T23_51_52.h"
#include "tvis.Kaguya.2009_03_27T20_32_29.h"
#include "tvis.Kaguya.2009_04_01T18_50_41.h"
#include "tvis.Kaguya.2009_04_05T13_09_19.h"
#include "tvis.Kaguya.2009_04_08T14_26_05.h"


Moon::Sampling::Kaguya::TVIS::cTvisOrientation Moon::Sampling::Kaguya::TVIS::individualPointsTvisOrientation[]={
{1238626280,  "2009-04-01T18:51:20",  -7.2000005, 283.58658/180.0*Pi, -23.701941/180.0*Pi, {0.0,0.0,0.0},{0.0,0.0,0.0}},
{1238951387,  "2009-04-05T13:09:47",  -0.099999458, 106.75014/180.0*Pi, 22.397675/180.0*Pi, {0.0,0.0,0.0},{0.0,0.0,0.0}},
{1239215182,  "2009-04-08T14:26:22",  -0.95000021,  283.85913/180.0*Pi, -21.983335/180.0*Pi, {0.0,0.0,0.0},{0.0,0.0,0.0}},
};
int Moon::Sampling::Kaguya::TVIS::nTotalIndividualPointsTvisOrientationElements=3;

vector<Moon::Sampling::Kaguya::TVIS::cTvisOrientationListElement> Moon::Sampling::Kaguya::TVIS::TvisOrientationVector;


//init the vector of TVIS orientations
void Moon::Sampling::Kaguya::TVIS::Init() {
  TVIS::TvisOrientationVector.resize(5);

  TVIS::TvisOrientationVector[0].TvisOrientation=TvisOrientation__2009_02_06T23_51_52;
  TVIS::TvisOrientationVector[0].nTotalTvisOrientationElements=nTotalTvisOrientationElements__2009_02_06T23_51_52;

  TVIS::TvisOrientationVector[1].TvisOrientation=TvisOrientation__2009_03_27T20_32_29;
  TVIS::TvisOrientationVector[1].nTotalTvisOrientationElements=nTotalTvisOrientationElements__2009_03_27T20_32_29;

  TVIS::TvisOrientationVector[2].TvisOrientation=TvisOrientation__2009_04_01T18_50_41;
  TVIS::TvisOrientationVector[2].nTotalTvisOrientationElements=nTotalTvisOrientationElements__2009_04_01T18_50_41;

  TVIS::TvisOrientationVector[3].TvisOrientation=TvisOrientation__2009_04_05T13_09_19;
  TVIS::TvisOrientationVector[3].nTotalTvisOrientationElements=nTotalTvisOrientationElements__2009_04_05T13_09_19;

  TVIS::TvisOrientationVector[4].TvisOrientation=TvisOrientation__2009_04_08T14_26_05;
  TVIS::TvisOrientationVector[4].nTotalTvisOrientationElements=nTotalTvisOrientationElements__2009_04_08T14_26_05;
}

