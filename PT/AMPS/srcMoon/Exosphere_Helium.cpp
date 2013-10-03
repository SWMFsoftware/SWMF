/*
 * Mercury.h
 *
 *  Created on: Feb 10, 2012
 *      Author: vtenishe
 */

//$Id$


#include "pic.h"
#include "Exosphere_Helium.h"




        double HeliumDesorption::SourceRate[PIC::nTotalSpecies],HeliumDesorption::maxLocalSourceRate[PIC::nTotalSpecies];
        cSingleVariableDiscreteDistribution<int> HeliumDesorption::SurfaceInjectionDistribution[PIC::nTotalSpecies];
//        double HeliumDesorption::GetSurfaceElementProductionRate(int nElement,int *spec);


        double HeliumDesorption::GetSurfaceElementProductionRate(int spec,int SurfaceElement,void *SphereDataPointer) {
                double res=0.0,norm[3],sun[3],HeliumSurfaceElementPopulation;

                HeliumSurfaceElementPopulation=((cInternalSphericalData*)SphereDataPointer)->SurfaceElementPopulation[spec][SurfaceElement];
                if (HeliumSurfaceElementPopulation<0.0) return 0.0;

                memcpy(sun,Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT,3*sizeof(double));
                memcpy(norm,(((cInternalSphericalData*)SphereDataPointer)->SurfaceElementExternalNormal+SurfaceElement)->norm,3*sizeof(double));

                //dayside production of helium
                for (int idim=0;idim<DIM;idim++) res+=sun[idim]*norm[idim];
                      res*=(res>0.0) ? sw_Helium_Ion_Flux_1AU*pow(_AU_/Exosphere::xObjectRadial,2.0) : 0.0;
                 return res;

        }


               bool HeliumDesorption::GenerateParticleProperties(int spec,double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double* v_IAU_OBJECT, double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere) {
                double ExternalNormal[3],CosSubSolarAngle;

                //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_OBJECT'
                double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];
                int nZenithElement,nAzimuthalElement;
                unsigned int el;

                el=SurfaceInjectionDistribution[spec].DistributeVariable();
                Exosphere::Planet->GetSurfaceElementIndex(nZenithElement,nAzimuthalElement,el);
                Exosphere::Planet->GetSurfaceElementRandomDirection(ExternalNormal,nZenithElement,nAzimuthalElement);

                CosSubSolarAngle=(Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT[0]*ExternalNormal[0])+
                    (Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT[1]*ExternalNormal[1])+
                    (Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT[2]*ExternalNormal[2]);

                x_LOCAL_IAU_OBJECT[0]=sphereRadius*ExternalNormal[0];
                x_LOCAL_IAU_OBJECT[1]=sphereRadius*ExternalNormal[1];
                x_LOCAL_IAU_OBJECT[2]=sphereRadius*ExternalNormal[2];

                ExternalNormal[0]*=-1.0;
                ExternalNormal[1]*=-1.0;
                ExternalNormal[2]*=-1.0;

                //transfer the position into the coordinate frame related to the rotating coordinate frame 'MSGR_SO'
                x_LOCAL_SO_OBJECT[0]=
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[0][0]*x_LOCAL_IAU_OBJECT[0])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[0][1]*x_LOCAL_IAU_OBJECT[1])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[0][2]*x_LOCAL_IAU_OBJECT[2]);

                x_LOCAL_SO_OBJECT[1]=
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[1][0]*x_LOCAL_IAU_OBJECT[0])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[1][1]*x_LOCAL_IAU_OBJECT[1])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[1][2]*x_LOCAL_IAU_OBJECT[2]);

                x_LOCAL_SO_OBJECT[2]=
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[2][0]*x_LOCAL_IAU_OBJECT[0])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[2][1]*x_LOCAL_IAU_OBJECT[1])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[2][2]*x_LOCAL_IAU_OBJECT[2]);


                //determine if the particle belongs to this processor
                startNode=PIC::Mesh::mesh.findTreeNode(x_LOCAL_SO_OBJECT,startNode);
                if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;

                //generate particle's velocity vector in the coordinate frame related to the planet 'IAU_OBJECT'
                double SurfaceTemperature,vbulk[3]={0.0,0.0,0.0};

                SurfaceTemperature=Exosphere::GetSurfaceTemeprature(CosSubSolarAngle,x_LOCAL_SO_OBJECT);
                PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);


                //transform the velocity vector to the coordinate frame 'MSGR_SO'
                v_LOCAL_SO_OBJECT[0]=
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[3][0]*x_LOCAL_IAU_OBJECT[0])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[3][1]*x_LOCAL_IAU_OBJECT[1])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[3][2]*x_LOCAL_IAU_OBJECT[2])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[3][3]*v_LOCAL_IAU_OBJECT[0])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[3][4]*v_LOCAL_IAU_OBJECT[1])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[3][5]*v_LOCAL_IAU_OBJECT[2]);

                v_LOCAL_SO_OBJECT[1]=
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[4][0]*x_LOCAL_IAU_OBJECT[0])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[4][1]*x_LOCAL_IAU_OBJECT[1])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[4][2]*x_LOCAL_IAU_OBJECT[2])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[4][3]*v_LOCAL_IAU_OBJECT[0])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[4][4]*v_LOCAL_IAU_OBJECT[1])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[4][5]*v_LOCAL_IAU_OBJECT[2]);

                v_LOCAL_SO_OBJECT[2]=
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[5][0]*x_LOCAL_IAU_OBJECT[0])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[5][1]*x_LOCAL_IAU_OBJECT[1])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[5][2]*x_LOCAL_IAU_OBJECT[2])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[5][3]*v_LOCAL_IAU_OBJECT[0])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[5][4]*v_LOCAL_IAU_OBJECT[1])+
                    (Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[5][5]*v_LOCAL_IAU_OBJECT[2]);

                memcpy(x_SO_OBJECT,x_LOCAL_SO_OBJECT,3*sizeof(double));
                memcpy(x_IAU_OBJECT,x_LOCAL_IAU_OBJECT,3*sizeof(double));
                memcpy(v_SO_OBJECT,v_LOCAL_SO_OBJECT,3*sizeof(double));
                memcpy(v_IAU_OBJECT,v_LOCAL_IAU_OBJECT,3*sizeof(double));

                return true;
              }



