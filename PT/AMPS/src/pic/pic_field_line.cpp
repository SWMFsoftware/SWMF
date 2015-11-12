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
    void cFieldLine::ResetSegmentWeights(){
      cFieldLineSegment *Segment = FirstSegment;
      double w[PIC::nTotalSpecies];
      //compute weights and normalize them
	for(int spec=0; spec < PIC::nTotalSpecies; spec++)
	  TotalWeight[spec] = 0;
      for(int iSegment=0; iSegment<nSegment; iSegment++){
	//compute weights
	_FIELDLINE_SEGMENT_WEIGHT_(w, Segment);
	//set segment's weight
	Segment->SetWeight(w);
	for(int spec=0; spec < PIC::nTotalSpecies; spec++)
	  TotalWeight[spec] += w[spec];
	//get the next one
	Segment = Segment->GetNext();
      }
    }

    //=========================================================================
    void cFieldLine::GetSegmentRandom(cFieldLineSegment *SegmentOut,
				      double& WeightCorrectionFactor,
				      int spec){
      // choose a random segment on field line and return pointer SegmentOut;
      // since statistical weights may vary, compute a correction factor
      // for a particle to be injected on this segment
      //-----------------------------------------------------------------------
      // choose uniformly a segment, i.e. weight_uniform = 1.0 / nSegment
      int iSegmentChoice = (int)(nSegment * rnd());
      // cycle through segment until get the chosen one
      SegmentOut = FirstSegment;
      for(int iSegment=0; iSegment < iSegmentChoice; iSegment++)
	SegmentOut = SegmentOut->GetNext();

      // SegmentOut now is a pointer to a uniformly chosen segment;
      // find the correction factor: weight_segment / weight_uniform
      //-----------------------------------------------------------------------
      WeightCorrectionFactor = 
	SegmentOut->GetWeight(spec) / TotalWeight[spec] * nSegment;
    }
    
    //=========================================================================
    void cFieldLine::SetMagneticField(double *BIn, int iVertex){
      cFieldLineVertex *Vertex;
      if(iVertex == -1)
	Vertex = LastVertex;
      else if(iVertex > 0.5*nSegment && iVertex <= nSegment){
	Vertex = LastVertex;
	for(int i=nSegment; i>iVertex; i--)
	  Vertex = Vertex->GetPrev();
      }
      else if(iVertex >= 0){
	Vertex = FirstVertex;
	for(int i=0; i<iVertex; i++)
	  Vertex = Vertex->GetNext();
      }
      else
	exit(__LINE__, __FILE__, "ERROR: invalid index of vertex");
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
    void FieldLineWeight_Uniform(double* Weight, cFieldLineSegment* Segment){
      for(int spec=0; spec<PIC::nTotalSpecies; spec++)
	Weight[spec] = 1.0;
    }

    //=========================================================================
    void InitLoop2D(double *xStart,   //start loop here
		    double DArc //increment in the angle of arc 
		    ){
      // in 2D case (e.g. cylindrical symmetry) generate a magnetic field-line
      // based on the background data
#if _PIC_SYMMETRY_MODE_ != _PIC_SYMMETRY_MODE__AXIAL_
      exit(__LINE__, __FILE__,"ERROR: implemented only for axial symmetry!");
#endif

      if(nFieldLine == nFieldLineMax)
	exit(__LINE__,__FILE__,"ERROR: reached limit for field line number");

      //increase counter of field lines
      nFieldLine++;

      double x[3] = {xStart[0], xStart[1], xStart[2]};
#if _PIC_SYMMETRY_MODE_ == _PIC_SYMMETRY_MODE__AXIAL_
      //rotate to the y=0 plane
      x[1] = 0;
      x[0] = pow(x[0]*x[0] + x[1]*x[1], 0.5) * ( (x[0]>0) ? 1 : -1);
#endif
      double xFirst[3] = {x[0],x[1],x[2]};

      // min value resolved for magnetic field (squared)
      const double epsB2 = 1E-30;
      const double epsB  = 1E-15;

      //magnetic field
      double B[3] = {0.0, 0.0, 0.0};
      double b[3] = {0.0, 0.0, 0.0};
      double absB =  0.0;
      CPLR::InitInterpolationStencil(x);
      CPLR::GetBackgroundMagneticField(B);

#if _PIC_SYMMETRY_MODE_ == _PIC_SYMMETRY_MODE__AXIAL_
      // need to have non-zero components in x-z plane
      if(B[0]*B[0] + B[2]*B[2] < epsB2)
	exit(__LINE__,__FILE__, 
	     "ERROR: magnetic field magnitude is below min resolved value");
      absB = pow(B[0]*B[0]+B[2]*B[2],0.5);
      b[0] = B[0]/absB; b[1] = 0; b[2] = B[2]/absB;
#endif
      
      //add the initial vertex
      FieldLinesAll[nFieldLine-1].Add(x);
      FieldLinesAll[nFieldLine-1].SetMagneticField(B,0);

      //housekeeping for controlling loop's generation
      // Arc = curvature * Length
      double Arc = 0.0;
      // direction of new segment
      double Dir[3] = {0.0, 0.0, 0.0};
      // angle swiped by the loop so far
      double Angle = 0.0;
      // estimate for the next segment's length
      double Length = 1E3;
      // position and magnetic field and next candidate vertex
      double xNew[3] = {0.0,0.0,0.0};
      double BNew[3] = {0.0,0.0,0.0};
      double bNew[3] = {0.0,0.0,0.0};
      double absBNew;
      //generate the loop; new vertices are added until:
      // - 2*Pi angle is swiped AND 
      //   candidate vertex is spatially close to the first one ANR
      //   directions of the first and last segments are close => SUCCESS
      // - 4*Pi angle is swiped => FAILURE
      // - exited the domain => FAILURE
      
      while(Angle < 4*Pi){

	// to control inner while loop
	static int countIn;
	const  int countInMax = 100;
	countIn = 0;
	do{
	  countIn++;

	  // new location - predictor
	  for(int idim=0; idim<DIM; idim++)
	    xNew[idim] = x[idim] + Length * b[idim];
	  
	  // get magnetic field at this location
	  CPLR::InitInterpolationStencil(xNew);
	  CPLR::GetBackgroundMagneticField(BNew);
#if _PIC_SYMMETRY_MODE_ == _PIC_SYMMETRY_MODE__AXIAL_
	  // need to have non-zero components in x-z plane
	  if(BNew[0]*BNew[0] + BNew[2]*BNew[2] < epsB2)
	    exit(__LINE__,__FILE__, 
		 "ERROR: magnetic field magnitude is below min resolved");
	  absBNew = pow(BNew[0]*BNew[0]+BNew[2]*BNew[2],0.5);
	  bNew[0] = BNew[0]/absBNew; bNew[1] = 0; bNew[2] = BNew[2]/absBNew;
#endif
	  
	  //find Arc = (curvature * Length)
	  Arc = 0.0;
	  for(int idim=0; idim<DIM; idim++)
	    Arc+= pow(BNew[idim]-B[idim], 2);
	  Arc = pow(Arc, 0.5);
	  Arc/= 0.5 * (absB + absBNew);
	  
	  //check condition to exit loop
	  if(Arc < 0.5 * DArc || Arc > DArc)
	    // too small or too large step; factor -> 1 as countIn grows
	    // to avoid reaching stationary points and infinite loop
	    Length *= 1. + (0.75 * DArc / Arc - 1.) / (1 + countIn);
	  else
	    break;
	  
	  if(countIn > countInMax)
	    exit(__LINE__, __FILE__,"ERROR: can't generate field-line loop ");
	  
	}while(true);
	
	//compute the final direction
#if _PIC_SYMMETRY_MODE_ == _PIC_SYMMETRY_MODE__AXIAL_
	Dir[0] = 0.5*(BNew[0]+B[0]);
	Dir[1] = 0;
	Dir[2] = 0.5*(BNew[2]+B[2]);
	double misc = pow(Dir[0]*Dir[0] + Dir[2]*Dir[2], 0.5);
	Dir[0]/= misc; Dir[2]/= misc;
#endif
	//new vertex location
	for(int idim=0; idim<DIM; idim++)
	  xNew[idim] = x[idim] +  Length * Dir[idim];
	
	//get magnetic field at this location
	CPLR::InitInterpolationStencil(xNew);
	CPLR::GetBackgroundMagneticField(BNew);
#if _PIC_SYMMETRY_MODE_ == _PIC_SYMMETRY_MODE__AXIAL_
	// need to have non-zero components in x-z plane
	if(BNew[0]*BNew[0] + BNew[2]*BNew[2] < epsB2)
	  exit(__LINE__,__FILE__, 
	       "ERROR: magnetic field magnitude is below min resolved");
	absBNew = pow(BNew[0]*BNew[0]+BNew[2]*BNew[2],0.5);
	bNew[0] = BNew[0]/absBNew; bNew[1] = 0; bNew[2] = BNew[2]/absBNew;
#endif

	//housekeeping
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// INCORRECT!!!!! FIX THIS!!!!!
	Angle += acos(b[0]*bNew[0] + b[1]*bNew[1] + b[2]*bNew[2]);
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	//add this vertex
	FieldLinesAll[nFieldLine-1].Add(xNew);
	FieldLinesAll[nFieldLine-1].SetMagneticField(BNew);

	//check conditions for completing the loop
	cFieldLineSegment* First=FieldLinesAll[nFieldLine-1].GetFirstSegment();
	cFieldLineSegment* Last =FieldLinesAll[nFieldLine-1].GetLastSegment();
	double Dist = 0.0; //distance between last and first verticies
	double DAngle = 0.0; //angle between last and first segments
	double DirFirst[3], DirLast[3];
	First->GetDir(DirFirst); Last->GetDir(DirLast);
	for(int idim=0; idim<DIM; idim++){
	  Dist  += pow(xNew[idim]-xFirst[idim], 2);
	  DAngle+= DirLast[idim] * DirFirst[idim];
	}
	Dist   = pow(Dist, 0.5);
	DAngle = acos(DAngle);
	if( fabs( 0.5*Angle/Pi - 1) < 0.1 &&
	    Dist < Length && DAngle < Pi/12){
	  cFieldLineVertex* LastEnd = Last->GetEnd();
	  LastEnd->SetX(xFirst);
	  Last->SetNext(First);
	  return;
	}

	//save current location and magnetic field
	x[0] = xNew[0]; x[1] = xNew[1]; x[2] = xNew[2];
	b[0] = bNew[0]; b[1] = bNew[1]; b[2] = bNew[2];
	B[0] = BNew[0]; B[1] = BNew[1]; B[2] = BNew[2];
	absB = absBNew;

	  
      }
      
      //      exit(__LINE__, __FILE__,"ERROR: can't generate field-line loop ");
      
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
	FieldLinesAll[nFieldLine-1].SetMagneticField(B,iSegment);

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

