//=========================================================================
//$Id$
//=========================================================================
//The storage for global variables used in the AMR mesh
#include "meshAMRdef.h"
#include "meshAMRinternalSurface.h"

double _MESH_AMR_XMAX_[3]={0.0,0.0,0.0},_MESH_AMR_XMIN_[3]={0.0,0.0,0.0};

 //the static data of spherical internal boundaries
long int cInternalSphericalData::nAzimuthalSurfaceElements=30,cInternalSphericalData::nZenithSurfaceElements=20;
double cInternalSphericalData::dAzimuthalAngle=2.0*Pi/cInternalSphericalData::nAzimuthalSurfaceElements;

#if  _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ == _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_COSINE_DISTRIBUTION_
double cInternalSphericalData::dCosZenithAngle=2.0/cInternalSphericalData::nZenithSurfaceElements;
#elif _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ == _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_UNIFORM_DISTRIBUTION_
double cInternalSphericalData::dZenithAngle=Pi/cInternalSphericalData::nZenithSurfaceElements;
#endif

//the static data for the circle internal boundaries
long int cInternalCircleData::nPolarSurfaceElements=30;
#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_AXIAL_SYMMETRY_
double cInternalCircleData::dPolarAngle=Pi/nPolarSurfaceElements;
#else
double cInternalCircleData::dPolarAngle=2.0*Pi/nPolarSurfaceElements;
#endif


