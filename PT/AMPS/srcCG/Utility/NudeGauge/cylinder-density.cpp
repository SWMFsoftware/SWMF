

#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>
#include <vector>

#define _DO_NOT_LOAD_GLOBAL_H_
#define _COMPILATION_MODE_ _COMPILATION_MODE__MPI_

#include "../../../src/general/rnd.h"
#include "../../../src/general/specfunc.h"

//const double Pi=3.141592654;

namespace Box {
  double xmin[3]={-1.0,-1.0,0.0};
  double xmax[3]={1.0,1.0,1.0};
}

namespace Gauge {
  double Radius=0.1;
  double Height=0.1; 
}

int nTotalTests=100000000;
double AveragingPhysicalTime=100.0;

double BackgroundDensity=100.0;
double BackgroundVelocity[3]={0.0,0.0,-100.0};

void EvaluateCase(double& Density) {
  int idim,iTest,nface,nd0,nd1,nd2;
  double ParticleWeight,x[3],v[3],norm[3],x0[3],e0[3],e1[3],ParticleTimeCounter=0.0;

  static const int nFaceNodes=4;
  static const int FaceNodeMap[2*3][nFaceNodes]={ {0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7}, {0,1,2,3}, {4,5,6,7}};  //the array dimension: FaceNodeMap[2*_MESH_DIMENSION_][nFaceNodes]
  static const int CornerNodesMap[8][3]={ {0,0,0}, {1,0,0}, {0,1,0}, {1,1,0},   {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}};

  
  for (nface=0;nface<6;nface++) {
    //get the geometry of each face 
    nd0=FaceNodeMap[nface][0];
    nd1=FaceNodeMap[nface][1];
    nd2=FaceNodeMap[nface][2];

    for (idim=0;idim<3;idim++) {
      x0[idim]=(CornerNodesMap[nd0][idim]==1) ? Box::xmax[idim] : Box::xmin[idim];

      e0[idim]=((CornerNodesMap[nd1][idim]==1) ? Box::xmax[idim] : Box::xmin[idim]) - x0[idim];
      e1[idim]=((CornerNodesMap[nd2][idim]==1) ? Box::xmax[idim] : Box::xmin[idim]) - x0[idim];
    }

    switch(nface) {
    case 0: case 1:
      norm[0]=(nface==0) ? -1.0 : 1.0;
      norm[1]=0.0,norm[2]=0.0;
      break;

    case 2: case 3:
      norm[1]=(nface==2) ? -1.0 : 1.0;
      norm[0]=0.0,norm[2]=0.0;
      break;

    case 4: case 5:
      norm[2]=(nface==4) ? -1.0 : 1.0;
      norm[0]=0.0,norm[1]=0.0;
      break;
    }

    //get the numerical flux through the face 
    double c=-Vector3D::DotProduct(norm,BackgroundVelocity); 
    double NumericalFlux=(c>0.0) ? c*Vector3D::Length(e0)*Vector3D::Length(e1) : 0.0;

    //if the numerical flux is positive -> determine the contribution of the flux to the density inside the gauge 
    double t1,t2,tmin,tmax,tStart;

    ParticleWeight=NumericalFlux*AveragingPhysicalTime/nTotalTests;
    
    if (ParticleWeight>0.0) for (iTest=0;iTest<nTotalTests;iTest++) {
      //generate a new particle 
      double c0=rnd(),c1=rnd(),a=0.0,b=0.0,c=0.0,d;
  
      for (idim=0;idim<3;idim++) {
        x[idim]=x0[idim]+c0*e0[idim]+c1*e1[idim];
        v[idim]=BackgroundVelocity[idim];
      }

      //correct the particle velocity such that it always has a radial velocity component
      if (v[0]*v[0]+v[1]*v[1]<1.0E-15*v[2]*v[2]) {
        double r,phi;

        phi=2.0*Pi*rnd();
        r=1.0E-15*fabs(v[2]);

        v[0]=r*sin(phi);
        v[1]=r*cos(phi);
      }

      //determine wether the particle can enter the gauge limits
      // 2 as the limit of the index is CORRECT!!!!!!!!!! 
      for (idim=0;idim<2;idim++) {
        a+=pow(v[idim],2);
        b+=2.0*v[idim]*x[idim];
        c+=pow(x[idim],2);
      }

      c-=pow(Gauge::Radius,2);
      d=pow(b,2)-4.0*a*c;
   
      if (d>0.0) {
        d=sqrt(d);
        t1=(-b+d)/(2.0*a);
        t2=(-b-d)/(2.0*a);

        tmin=std::min(t1,t2);
        tmax=std::max(t1,t2);

        if (tmax>0.0) {
          //the particle can enter the gauge 
          //determine the time that the particle will spend in the gauge 
          double tBottomPlaneIntersection,tTopPlaneIntersection;
          double tEnter,tExit,tCylinderExit;
          bool BottomPlaneIntersectionFlag=false;
          bool ZeroVzFlag=false; 

          //determine wether particle has vz non-zero
          ZeroVzFlag=(tmax*fabs(v[2])>1.0E-6*(Box::xmax[2]-Box::xmin[2])) ? false : true; //false -> vz is non-zero; true -> vz is zero

          if (ZeroVzFlag==false) {
            //vz componen of the partilce velocity is non-zero
            //determine the time of intersection of the particle trajectories with the planes defined by the top and the botton of the gauges
            tTopPlaneIntersection=(Gauge::Height-x[2])/v[2];
            tBottomPlaneIntersection=-x[2]/v[2];



            //the partice can anter into the gauge only if tmax>tTopPlaneIntersection
            if (tmax>tTopPlaneIntersection) {
              //the particle can enter into the gauge
              //determine the enterence time 


              //check wether the initial location is already within the projection of the cylinder
              if ((x[0]*x[0]+x[1]*x[1]<Gauge::Radius*Gauge::Radius)) {
                //the particle can enter only when tmax>tTopPlaneIntersection

                if (tmax>tTopPlaneIntersection) {
                  tEnter=tTopPlaneIntersection;

                  if (tmax>tBottomPlaneIntersection) {
                    tExit=tBottomPlaneIntersection;
                    BottomPlaneIntersectionFlag=true;
                  }
                  else {
                    tExit=tmax;
                    BottomPlaneIntersectionFlag=false;
                  }
                }
                else {
                  //the particle exit the culinder befor entering intothe gauge
                  BottomPlaneIntersectionFlag=false;
                  continue;
                }

              }
              else {
                //projection of the initial locaiton is the particle is outdile of the cylinder
                if (tmin>tBottomPlaneIntersection) {
                  //the particle cannot enter into the gauge
                  BottomPlaneIntersectionFlag=false;
                  continue;
                }
                else {
                  tEnter=std::max(tTopPlaneIntersection,tmin);

                  if (tmax>tBottomPlaneIntersection) {
                    //the particle extits through the bottom of the gauge
                    tExit=tBottomPlaneIntersection;
                    BottomPlaneIntersectionFlag=true;
                  }
                  else {
                    //the particle exits through the side of the gauge
                    tExit=tmax;
                    BottomPlaneIntersectionFlag=false;
                  }
                }
              }

              //register the time that the particle has spent in the gauge
              ParticleTimeCounter+=(tExit-tEnter)*ParticleWeight;

              //generate a scattered particles if needed
              if (BottomPlaneIntersectionFlag==true) {
                for (idim=0;idim<3;idim++) x[idim]+=v[idim]*tBottomPlaneIntersection;
  
                //generate the new particle velocity, and determine the time that the particle will spent in the gauge 
                double SpeedInit=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

                do {
                  Vector3D::Distribution::Uniform(v);
                }
                while (v[2]<0.0);

                for (idim=0;idim<3;idim++) v[idim]*=SpeedInit;

                //correct the particle velocity such that it always has a radial velocity component
                if (v[0]*v[0]+v[1]*v[1]<1.0E-15*v[2]*v[2]) {
                  double r,phi;

                  phi=2.0*Pi*rnd();
                  r=1.0E-15*fabs(v[2]);

                  v[0]=r*sin(phi);
                  v[1]=r*cos(phi);
                }

                for (a=0.0,b=0.0,c=0.0,idim=0;idim<3;idim++) {
                  a+=pow(v[idim],2);
                  b+=2.0*v[idim]*x[idim];
                  c+=pow(x[idim],2);
                }

                c-=pow(Gauge::Radius,2);
                d=pow(b,2)-4.0*a*c;

                d=sqrt(d);
                t1=(-b+d)/(2.0*a);
                t2=(-b-d)/(2.0*a);

                tmin=std::min(t1,t2);
                tmax=std::max(t1,t2);

                tCylinderExit=(tmin>0.0) ? tmin : tmax; 
                tTopPlaneIntersection=Gauge::Height/v[2];

                //register the trajectory interval to the bottom plane
                ParticleTimeCounter+=std::min(tCylinderExit,tTopPlaneIntersection)*ParticleWeight;
              }
           }
         }
         else {
          //vz componet is zero
          if (x[2]<Gauge::Height) ParticleTimeCounter+=(tmax-tmin)*ParticleWeight;
         }

       }
     }
   }
 }


//  std::cout << "Density=" << ParticleTimeCounter/AveragingPhysicalTime/(Gauge::Height*Pi*pow(Gauge::Radius,2)) << ", Volume=" << Gauge::Height*Pi*pow(Gauge::Radius,2) << std::endl;

  Density=ParticleTimeCounter/AveragingPhysicalTime/(Gauge::Height*Pi*pow(Gauge::Radius,2));
}


int main() {
  int i;
  double phi,Density;

  const double Speed=100.0;

  std::cout << "Angle\tDensity\t1+sin(phi)\n";

  for (i=0;i<90;i++) {
    phi=i*Pi/180.0;

    BackgroundVelocity[0]=-Speed*cos(phi);
    BackgroundVelocity[1]=0.0;
    BackgroundVelocity[2]=-Speed*sin(phi);

    EvaluateCase(Density);
    std::cout << phi/Pi*180.0 << "\t  " <<  Density << "\t  " << 1.0+sin(phi) << std::endl;
  }
  
  return 1;
}





 

