//$Id$

//implementation of methods in namespace FieldLine

#include "pic.h"

cStack<PIC::FieldLine::cFieldLineVertex>  PIC::FieldLine::VerticesAll;
cStack<PIC::FieldLine::cFieldLineSegment> PIC::FieldLine::SegmentsAll;
PIC::FieldLine::cFieldLine *PIC::FieldLine::FieldLinesAll = NULL;

long int PIC::FieldLine::nFieldLine=0;

namespace PIC{
  namespace FieldLine{
    //=========================================================================
    //set whole state vector on vertex
    void cFieldLineVertex::SetStateVector(double* ElectricFieldIn,
					  double* MagneticFieldIn,
					  double* PlasmaVelocityIn,
					  double  PlasmaDensityIn,
					  double  PlasmaTemperatureIn,
					  double  PlasmaPressureIn){
      SetElectricField(ElectricFieldIn);
      SetMagneticField(MagneticFieldIn);
      SetPlasmaVelocity(PlasmaVelocityIn);
      SetPlasmaDensity(PlasmaDensityIn);
      SetPlasmaTemperature(PlasmaTemperatureIn);
      SetPlasmaPressure(PlasmaPressureIn);
    }
    //=========================================================================
    //get whole state vector on vertex
    void cFieldLineVertex::GetStateVector(double* ElectricFieldOut,
					  double* MagneticFieldOut,
					  double* PlasmaVelocityOut,
					  double& PlasmaDensityOut,
					  double& PlasmaTemperatureOut,
					  double& PlasmaPressureOut){
      GetElectricField(ElectricFieldOut);
      GetMagneticField(MagneticFieldOut);
      GetPlasmaVelocity(PlasmaVelocityOut);
      GetPlasmaDensity(PlasmaDensityOut);
      GetPlasmaTemperature(PlasmaTemperatureOut);
      GetPlasmaPressure(PlasmaPressureOut);
    }
    //=========================================================================
    //get whole state vector on segment
    void cFieldLineSegment::GetStateVector(double  S, //position 0<=S<=1
					   double* ElectricFieldOut,
					   double* MagneticFieldOut,
					   double* PlasmaVelocityOut,
					   double& PlasmaDensityOut,
					   double& PlasmaTemperatureOut,
					   double& PlasmaPressureOut){
      GetElectricField(    S, ElectricFieldOut);
      GetMagneticField(    S, MagneticFieldOut);
      GetPlasmaVelocity(   S, PlasmaVelocityOut);
      GetPlasmaDensity(    S, PlasmaDensityOut);
      GetPlasmaTemperature(S, PlasmaTemperatureOut);
      GetPlasmaPressure(   S, PlasmaPressureOut);
    }
    //=========================================================================
    bool cFieldLine::is_broken(){
      int count;
      cFieldLineSegment* Segment;
      //  cFieldLineVertex*  Vertex; 
      
      //check connectivity forward
      count = 1, Segment = FirstSegment;//, Vertex  = FirstVertex;
      for(int iSegment=0; iSegment<nSegment; iSegment++){
	if(Segment->GetNext()==NULL //|| Vertex->GetNext()==NULL) break;
	   )break;
	Segment = Segment->GetNext();
	//    Vertex  = Vertex->GetNext();
	count++;
      }
      if(count<nSegment || Segment != LastSegment)// || Vertex != LastVertex)
	return true;
      
      //check connectivity backward
      count = 1, Segment = LastSegment;//, Vertex  = LastVertex;
      for(int iSegment=0; iSegment<nSegment; iSegment++){
	if(Segment->GetPrev()==NULL// || Vertex->GetPrev()==NULL) break;
	   )break;
	Segment = Segment->GetPrev();
	//    Vertex  = Vertex->GetPrev();
	count++;
      }
      if(count<nSegment || Segment != FirstSegment)// || Vertex != FirstVertex)
	return true;
      
      //the line is fully connected
      return false;
    }
    
    //=========================================================================
    void cFieldLine::Add(double *xIn){
      // check if field lineis unset
      if(IsSet == 0){
	if(FirstVertex==NULL){
	  //allocate the first vertex
	  LastVertex = (FirstVertex = VerticesAll.newElement());
	  FirstVertex->SetX(xIn);
	  //the first segment can't be allocated yet
	  nSegment = 0; TotalLength = 0.0; IsSet = 0;
	}
	else{
	  //allocate the second vertex 
	  LastVertex = VerticesAll.newElement();
	  LastVertex->SetX(xIn);
	  LastVertex->SetPrev(FirstVertex);
	  FirstVertex->SetNext(LastVertex);
	  //allocate the first segment
	  LastSegment = (FirstSegment = SegmentsAll.newElement());
	  LastSegment->SetVertices(LastVertex->GetPrev(), LastVertex);
	  nSegment++; TotalLength+= LastSegment->GetLength(); IsSet = 1;
	}
	return;
      }
      
      //allocate vertex
      cFieldLineVertex* newVertex;
      newVertex = VerticesAll.newElement();
      newVertex->SetX(xIn);
      
      //connect it to the last vertex in the field line
      LastVertex->SetNext(newVertex);
      newVertex->SetPrev(LastVertex);
      
      //now the new vertex is the last one
      LastVertex = newVertex;
      
      //allocate segment
      cFieldLineSegment* newSegment;
      newSegment = SegmentsAll.newElement();
      newSegment->SetVertices(LastVertex->GetPrev(), LastVertex);
      
      //connect it to the last segment in the line
      LastSegment->SetNext(newSegment);
      newSegment->SetPrev(LastSegment);
      
      //now the new segment is the last one
      LastSegment = newSegment;
      
      //update house-keeping data
      nSegment++;
      TotalLength += LastSegment->GetLength();
    }

    //=========================================================================
    void cFieldLine::SetMagneticField(double *BIn, int iVertex){
      cFieldLineVertex *Vertex;
      if(iVertex > 0.5*nSegment){
	Vertex = LastVertex;
	for(int i=nSegment; i>iVertex; i--)
	  Vertex = Vertex->GetPrev();
      }
      else{
	Vertex = FirstVertex;
	for(int i=0; i<iVertex; i++)
	  Vertex = Vertex->GetNext();
      }
      Vertex->SetMagneticField(BIn);
    }
    
    //=========================================================================
    void cFieldLine::Output(FILE* fout, bool GeometryOnly=false){
      
      if(!GeometryOnly)
	exit(__LINE__, __FILE__,"Not implemented for multiple processors");
      
      cFieldLineVertex *Vertex=FirstVertex;
      for(int iVertex=0; iVertex<=nSegment; iVertex++){

	//print coordinates
	double x[DIM];
	Vertex->GetX(x);
	for(int idim=0; idim<DIM; idim++) {fprintf(fout, "%e ", x[idim]);}

	//print magnetic field
	double B[DIM];
	Vertex->GetMagneticField(B);
	for(int idim=0; idim<DIM; idim++) {fprintf(fout, "%e ", B[idim]);}

	fprintf(fout,"\n");
	Vertex = Vertex->GetNext();
      }
    }
    
    //=========================================================================
    void Init(){
      
      if(PIC::nTotalThreads > 1)
	exit(__LINE__, __FILE__,"Not implemented for multiple processors");
      
      //allocate container for field lines
      FieldLinesAll = new cFieldLine [nFieldLineMax];
      
    }
    
    //=========================================================================
    void Output(char* fname, bool GeometryOnly){
      
      FILE* fout;
      fout = fopen(fname,"w");
      
      fprintf(fout, "TITLE=\"Field line geometry\"");
      
#if DIM == 3 
      fprintf(fout,"VARIABLES=\"x\",\"y\",\"z\",\"Bx\",\"By\",\"Bz\"\n");
#else 
      exit(__LINE__,__FILE__,"not implemented");
#endif
      
      for(int iFieldLine=0; iFieldLine<nFieldLine; iFieldLine++){
	fprintf(fout,"ZONE T=\"Field-line %i\" F=POINT\n",iFieldLine);
	FieldLinesAll[iFieldLine].Output(fout, true);
      }
      
      fclose(fout);
      
    }
    
    //=========================================================================
    void InitSimpleParkerSpiral(double *xStart){
      
#if DIM != 3
      exit(__LINE__,__FILE__,"Implemetned only for 3D case");
#endif

      if(nFieldLine == nFieldLineMax)
	exit(__LINE__,__FILE__,"ERROR: reached limit for field line number");

      //increase counter of field lines
      nFieldLine++;
      
      //spatial step and # of segments
      const double Ds=1E9;
      int nSegment = 1000;
           
      //Parker spiral starts @ 10 Sun radii
      const double R0 = 10 *_SUN__RADIUS_;
      //magnetic field at this location
      const double B0 = 1.83E-6;
      //velocity of solar wind
      const double SpeedSolarWind = 4.0E5;
      //rate of solar rotation at equator (rad/sec)
      const double Omega = 2.97211E-6;      

      //position of the Sun
      double xSun[DIM]={0};
      //relative position to Sun
      double x2Sun[DIM]={0};
      //distance to Sun
      double R2Sun=0;
      //distance to sun in equatorial plane
      double r2Sun=0;

      //find original latitude
      for(int idim=0; idim<DIM; idim++) x2Sun[idim] = xStart[idim] - xSun[idim];
      R2Sun = sqrt(x2Sun[0]*x2Sun[0]+x2Sun[1]*x2Sun[1]+x2Sun[2]*x2Sun[2]);

      //original latitude 
      double costheta = x2Sun[2] / R2Sun;
      double sintheta = sqrt(1 - costheta*costheta); // >=0 always     

      //start position
      double x[DIM]={xStart[0],xStart[1],xStart[2]};


      //if too close to the Sun => move to distance R0
      if(R2Sun < R0)
	for(int idim=0; idim<DIM; idim++) 
	  x[idim] = xSun[idim] + x2Sun[idim] * R0/R2Sun;
      
      
      for(int iSegment=0;iSegment<nSegment; iSegment++){

	//distances to sun
	for(int idim=0; idim<DIM; idim++) x2Sun[idim]=x[idim]-xSun[idim];
	R2Sun = sqrt(x2Sun[0]*x2Sun[0]+x2Sun[1]*x2Sun[1]+x2Sun[2]*x2Sun[2]);
	r2Sun = sqrt(x2Sun[0]*x2Sun[0]+x2Sun[1]*x2Sun[1]);

	//angle in the equatorial plane
	double cosphi = x2Sun[0]/r2Sun, sinphi=x2Sun[1]/r2Sun;

	double misc1 = B0 * pow(R0/R2Sun, 2);
	double misc2 = (R2Sun - R0) * Omega/SpeedSolarWind * sintheta;

	//magnetic field at location x
	double B[DIM];
	B[0] = misc1 * (sintheta*cosphi + misc2 * sinphi);
	B[1] = misc1 * (sintheta*sinphi - misc2 * cosphi);
	B[2] = misc1 *  costheta;

	//add vertex to the spiral
	FieldLinesAll[nFieldLine-1].Add(x);
	FieldLinesAll[nFieldLine-1].SetMagneticField(B,nSegment);

	//magnitude of the magnetic field
	double AbsB = sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);
	
	//direction of the magnetic field
	double b[DIM] = {B[0]/AbsB, B[1]/AbsB, B[2]/AbsB};
	
	// next vertex of the spiral
	x[0]+=b[0]*Ds; x[1]+=b[1]*Ds; x[2]+=b[2]*Ds;
      }
      
    }
    
  }
}

