//$Id$

//implementation of methods in namespace FieldLine

#include "pic.h"

//static variables of class cFieldLineVertex
int PIC::FieldLine::cFieldLineVertex::totalAssociatedDataLength=-1;
int PIC::FieldLine::cFieldLineVertex::sampleDataLength=-1;
int PIC::FieldLine::cFieldLineVertex::CollectingSamplingOffset=-1;
int PIC::FieldLine::cFieldLineVertex::CompletedSamplingOffset=-1;

namespace PIC{
  namespace FieldLine{

    cDatumStored DatumAtVertexElectricField(3,"\"Ex [V/m]\",\"Ey [V/m]\",\"Ez [V/m]\",",false);
    cDatumStored DatumAtVertexMagneticField(3,"\"Bx [nT]\",\"By [nT]\",\"Bz [nT]\",",true);
    cDatumStored DatumAtVertexPlasmaVelocity(3,"\"Plasma Vx [m/s]\",\"Plasma Vy [m/s]\",\"Plasma Vz [m/s]\"",true);
    cDatumStored DatumAtVertexPlasmaDensity(1,"\"Plasma number density [1/m^3]\"", false);
    cDatumStored DatumAtVertexPlasmaTemperature(1,"\"Plasma Temperature [K]\"", false);
    cDatumStored DatumAtVertexPlasmaPressure(1,"\"Plasma pressure [Pa]\"", false);
    cDatumStored DatumAtVertexMagneticFluxFunction(1,"\"MagneticFluxFunction [nT*m]\"", true);

    cDatumTimed DatumAtVertexParticleWeight(1,"\"Particle Weight\"",false);
    cDatumTimed DatumAtVertexParticleNumber(1,"\"Particle Number\"",true);
    cDatumTimed DatumAtVertexNumberDensity(1,"\"Number Density[1/m^3]\"",true);
    cDatumWeighted DatumAtVertexParticleEnergy(1,"\"Kinetic energy [J]\"",true);

    cDatumWeighted DatumAtGridParticleEnergy(1,"\"Kinetic energy [J]\"",true);

    vector<cDatumStored*> DataStoredAtVertex;
    vector<cDatumSampled*> DataSampledAtVertex;

    cAssociatedDataAMRstack<cFieldLineVertex>  VerticesAll;
    cAssociatedDataAMRstack<cFieldLineSegment> SegmentsAll;
    cFieldLine *FieldLinesAll = NULL;
    
    long int nFieldLine=0;

    double TimeLastUpdate = -1;

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
    cFieldLineVertex* cFieldLine::Add(double *xIn){
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
	return LastVertex;
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

      return LastVertex;
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
    void cFieldLine::GetSegmentRandom(int& iSegment,//cFieldLineSegment *SegmentOut,
				      double& WeightCorrectionFactor,
				      int spec){
      // choose a random segment on field line and return pointer SegmentOut;
      // since statistical weights may vary, compute a correction factor
      // for a particle to be injected on this segment
      //-----------------------------------------------------------------------
      // choose uniformly a segment, i.e. weight_uniform = 1.0 / nSegment
      int iSegmentChoice = (int)(nSegment * rnd());
      iSegment = iSegmentChoice;
      // cycle through segment until get the chosen one
      //      SegmentOut = FirstSegment;
      //      for(int iSegment=0; iSegment < iSegmentChoice; iSegment++)
      //	SegmentOut = SegmentOut->GetNext();

      // SegmentOut now is a pointer to a uniformly chosen segment;
      // find the correction factor: weight_segment / weight_uniform
      //-----------------------------------------------------------------------
      //      WeightCorrectionFactor = 
      //      	SegmentOut->GetWeight(spec) / TotalWeight[spec] * nSegment;
    }
    
    //=========================================================================
    long int InjectParticle(int spec){
      // this is a wrapper that can call either the default injection procedure
      // or a user-defined procedure
      return InjectParticle_default(spec);
    }
    
    //=========================================================================
    long int InjectParticle_default(int spec){
      //namespace alias
      namespace PB = PIC::ParticleBuffer;
      // pointer particle to the particle to be injected
      long int ptr;
      PB::byte* ptrData;
      ptrData = new PB::byte[PB::GetParticleDataLength()];
      
      //default settings
      PB::SetIndividualStatWeightCorrection(1.0,ptrData);
      
      
      PB::SetI(spec, ptrData);
      
      // pick a random field line
      int iFieldLine = (int)(nFieldLine * rnd());
      PB::SetFieldLineId(iFieldLine, ptrData);
      // inject particle onto this field line
      
      int iSegment;
      double WeightCorrection;
      //    FieldLinesAll[iFieldLine].GetSegmentRandom(iSegment,
      //					       WeightCorrection, spec);
      
      //Inject at the beginning of the field line FOR PARKER SPIRAL
      iSegment = 0;
      FieldLinesAll[iFieldLine].GetSegmentRandom(iSegment,
						 WeightCorrection, spec);
      
      
      cFieldLineSegment* Segment=FieldLinesAll[iFieldLine].GetSegment(iSegment);

      double S = iSegment + rnd();
      PB::SetFieldLineCoord(S, ptrData);
      double x[3], v[3];
      FieldLinesAll[iFieldLine].GetCartesian(x, S);
      PB::SetX(x, ptrData);
      
      // fix kintic energy and pitch angle far now
      //    double cosPhi=pow(2,-0.5);
      //    double vpar = 1e6 * cosPhi;
      //    double KinEnergyPerp=5E-17*(1-cosPhi*cosPhi);
      //direction of velocity: INJECT ALONG HE FIELD LINE FOR PARKER SPIRAL
      //    vpar*=(rnd()<0.5)?1:-1;
      
      // inject power law, 
      // see Numerical Studies of the Solar Energetic Particle 
      // Transport and Acceleration, Tenishev et al.
      int q=3;
      double vmin=1e6, vmax=1e6;
      double pvmin = pow(vmin, 1-q), pvmax = pow(vmax, 1-q);
      double r= rnd();
      double absv = pow( (1-r)*pvmin + r*pvmax, 1.0/(1-q));
      //direction of velocity: INJECT ALONG HE FIELD LINE FOR PARKER SPIRAL
      double cosPhi = 1 - rnd();
      double vpar = absv*cosPhi;
      double KinEnergyPerp = 9.1E-31 * 0.5 * (absv*absv - vpar*vpar);
      
      
      //velocity is paral to the field line
      Segment->GetDir(v);
      for(int i=0; i<3; i++)v[i]*=vpar;
      PB::SetV(v, ptrData);
      
      //magnetic field
      double B[3], AbsB;
      Segment->GetMagneticField(S, B);
      AbsB = pow(B[0]*B[0]+B[1]*B[1]+B[2]*B[2], 0.5)+1E-15;
      
      //magnetic moment
      PB::SetMagneticMoment(KinEnergyPerp/AbsB, ptrData);

      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
      node=PIC::Mesh::mesh.findTreeNode(x);
      
      long int res=PB::InitiateParticle(x,NULL,NULL,NULL,
					ptrData,
					_PIC_INIT_PARTICLE_MODE__ADD2LIST_,
					(void*)node);
      delete [] ptrData;
      double misc = PB::GetFieldLineCoord(res);
      return res;
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
    void cFieldLine::GetMagneticField(double* BOut, double S){
      // check correctness
      if(S < 0.0 || S > nSegment)
	exit(__LINE__,__FILE__,
	     "ERROR: trying to get magnetic field at an invalid location");

      // interpolate the magnetic field at the location S:
      //  floor(S) is the number of the segment,
      //  S - floor(S) is the location along segment (between 0 & 1)
      //-----------------------------------------------------------------------
      // number of the begin vertex
      int iSegment = (int) S;
      cFieldLineSegment* Segment = GetSegment(iSegment);
      Segment->GetMagneticField(S - iSegment, BOut);
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

	double Value[3];

	vector<cDatumStored*>::iterator itrDatumStored;
	for(itrDatumStored = DataStoredAtVertex.begin();
	    itrDatumStored!= DataStoredAtVertex.end();  itrDatumStored++)
	  if((*itrDatumStored)->doPrint){
	    Vertex->GetDatum(*(*itrDatumStored), Value);
	    for(int i=0; i<(*itrDatumStored)->length; i++){
	      fprintf(fout, "%e ", Value[i]);}
	  }
      
	vector<cDatumSampled*>::iterator itrDatum;
	for(itrDatum = DataSampledAtVertex.begin();
	    itrDatum!= DataSampledAtVertex.end();  itrDatum++)
	  if((*itrDatum)->doPrint){
	    cDatumTimed*    ptrDatumTimed;
	    cDatumWeighted* ptrDatumWeighted;
	    if((*itrDatum)->type == PIC::Datum::cDatumSampled::Timed_){
	      ptrDatumTimed = static_cast<cDatumTimed*> ((*itrDatum));
	      Vertex->GetDatumAverage(*ptrDatumTimed, Value, 0);}
	    else {
	      ptrDatumWeighted = static_cast<cDatumWeighted*> ((*itrDatum));
	      Vertex->GetDatumAverage(*ptrDatumWeighted, Value, 0);}
	    for(int i=0; i<(*itrDatum)->length; i++){
	      fprintf(fout, "%e ", Value[i]);}
	  }
	

	fprintf(fout,"\n");
	//flush completed sampling buffer
	Vertex->flushCollectingSamplingBuffer();
	Vertex = Vertex->GetNext();
      }
    }
    
    //=========================================================================
    void Init(){
      
      if(PIC::nTotalThreads > 1)
	exit(__LINE__, __FILE__,"Not implemented for multiple processors");
      
      // allocate container for field lines
      FieldLinesAll = new cFieldLine [nFieldLineMax];

      // activate data storage
      long int Offset = 0;
      // activate data that are stored but NOT sampled
      DatumAtVertexMagneticField.       activate(Offset, &DataStoredAtVertex);
      DatumAtVertexElectricField.       activate(Offset, &DataStoredAtVertex);
      DatumAtVertexPlasmaVelocity.      activate(Offset, &DataStoredAtVertex);
      DatumAtVertexPlasmaDensity.       activate(Offset, &DataStoredAtVertex);
      DatumAtVertexPlasmaTemperature.   activate(Offset, &DataStoredAtVertex);
      DatumAtVertexPlasmaPressure.      activate(Offset, &DataStoredAtVertex);
      DatumAtVertexMagneticFluxFunction.activate(Offset, &DataStoredAtVertex);
      // activate data that is sampled
      long int SamplingOffset = Offset; Offset = 0;
      DatumAtVertexParticleWeight.activate(Offset, &DataSampledAtVertex);
      DatumAtVertexParticleNumber.activate(Offset, &DataSampledAtVertex);
      DatumAtVertexNumberDensity. activate(Offset, &DataSampledAtVertex);
      DatumAtVertexParticleEnergy.activate(Offset, &DataSampledAtVertex);

      // assign offsets and data length
      cFieldLineVertex::SetDataOffsets(SamplingOffset, Offset);
      
      PIC::IndividualModelSampling::DataSampledList.push_back(&DatumAtGridParticleEnergy);
    }
    
    //=========================================================================
    void Output(char* fname, bool GeometryOnly){

      //swap sampling offsets
      cFieldLineVertex::swapSamplingBuffers();
      
      FILE* fout;
      fout = fopen(fname,"w");
      
      fprintf(fout, "TITLE=\"Field line geometry\"");
      
#if DIM == 3 
      fprintf(fout,"VARIABLES=\"x\",\"y\",\"z\"");//,\"Bx\",\"By\",\"Bz\"");
      vector<cDatumStored*>::iterator itrDatumStored;
      vector<cDatumSampled*>::iterator itrDatum;
      //      double* Value;
      for(itrDatumStored = DataStoredAtVertex.begin();
	  itrDatumStored!= DataStoredAtVertex.end();  itrDatumStored++)
	if((*itrDatumStored)->doPrint)
	  (*itrDatumStored)->PrintName(fout);
      for(itrDatum = DataSampledAtVertex.begin();
	  itrDatum!= DataSampledAtVertex.end();  itrDatum++)
	if((*itrDatum)->doPrint)
	  (*itrDatum)->PrintName(fout);
      fprintf(fout,"\n");
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
    void Sampling(long int ptr, double Weight, char* CellSamplingBuffer){
      // namespace alias
      namespace PB = PIC::ParticleBuffer;

      // find the field line and location along it
      int iFieldLine = PB::GetFieldLineId(ptr);
      double S = PB::GetFieldLineCoord(ptr);

      // weights of vertices
      double w = S - (int)S;

      //magnetic field
      double B[3];
      FieldLinesAll[iFieldLine].GetMagneticField(B,S);
      double AbsB = pow(B[0]*B[0]+B[1]*B[1]+B[2]*B[2], 0.5);

      //magnetic moment, mass, velocity and energy
      double mu = PB::GetMagneticMoment(ptr);
      int spec = PB::GetI(ptr);
      double v[3];
      double x[3];
      PB::GetV(v, ptr);
      PB::GetX(x, ptr);
      double m0= PIC::MolecularData::GetMass(spec);
      double absv=pow(v[0]*v[0]+v[1]*v[1]+v[2]*v[2], 0.5);
      double E = mu*AbsB + 0.5 * m0 * (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

      // volume
      double volume = FieldLinesAll[iFieldLine].GetSegment(S)->GetLength();
      const double R0 = 10 *_SUN__RADIUS_;
      const double B0 = 1.83E-6;





      // FIX LATER=====================================
      // take gyroradius of electron with 10^7 m/s speed at B0 (~30 m)
      volume *= 3.14 * 900 * (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]) / (R0*R0);
      // FIX LATER=====================================





      //sample to vertices
      cFieldLineVertex* V=FieldLinesAll[iFieldLine].GetSegment(S)->GetBegin();
      V->SampleDatum(DatumAtVertexNumberDensity,Weight/volume, (1-w));
      V->SampleDatum(DatumAtVertexParticleWeight,Weight,(1-w));
      V->SampleDatum(DatumAtVertexParticleEnergy,Weight*E,(1-w));
      //      V->SampleDatum(EnergyFlux(Weight*E*absv/volume*(1-w));
      V = V->GetNext();
      V->SampleDatum(DatumAtVertexNumberDensity,Weight/volume, (w));
      V->SampleDatum(DatumAtVertexParticleWeight,Weight,(w));
      V->SampleDatum(DatumAtVertexParticleEnergy,Weight*E,w);
//
//      V->SampleNumberDensity(Weight/volume * w);
//      V->SampleParticleWeight(Weight*w);
//      V->SampleParticleEnergy(E*Weight*w);
      //      V->SampleEnergyFlux(E*Weight*absv/volume*w);      

      //.........................
      *((double*)(CellSamplingBuffer+DatumAtGridParticleEnergy.offset))+=Weight*E;
	

    }

    //=========================================================================
    void FieldLineWeight_Uniform(double* Weight, cFieldLineSegment* Segment){
      for(int spec=0; spec<PIC::nTotalSpecies; spec++)
	Weight[spec] = 1.0;
    }

    //=========================================================================
    void InitLoop2D(double *xStart,  //start loop here
		    double DArc      //increment in the angle of arc 
		    ){
      // in 2D case (e.g. cylindrical symmetry) generate a magnetic field line;
      // magnetic flux function, psi, is utilized:
      // loop is a line of constant value
      //-----------------------------------------------------------------------
      // check correctness
#if _PIC_SYMMETRY_MODE_ != _PIC_SYMMETRY_MODE__AXIAL_
      exit(__LINE__, __FILE__,"ERROR: implemented only for axial symmetry!");
#endif
      // check for limit of number of field lines allowed
      if(nFieldLine == nFieldLineMax)
	exit(__LINE__,__FILE__,"ERROR: reached limit for field line number");

      // mark time of extraction
      TimeLastUpdate = PIC::SimulationTime::Get();

      // list of variables utilized
      //.......................................................................
      // current location
      double x[3] = {xStart[0], xStart[1], xStart[2]};
      // candidate location
      //      double xNew[3] = {0.0, 0.0, 0.0};
      // the first location; used to close the loop
      double xFirst[3] = {xStart[0], xStart[1], xStart[2]};;
      // min value resolved for magnetic field (squared)
      const double epsB2 = 1E-30;
      const double epsB  = 1E-15;
      // plasma velocity, stored at vertex
      double V[3] = {0.0, 0.0, 0.0};
      // magnetic field vector, unit vector and magnitude
      // NOTE: y-component is ignored!!!
      double B[3] = {0.0, 0.0, 0.0};
      double b[3] = {0.0, 0.0, 0.0};
      double absB =  0.0;
      // flux function
      double Psi0 = 0.0, Psi = 0.0;
      //housekeeping for controlling loop's generation
      // Arc = curvature * Length
      double Arc = 0.0;
      // direction of new segment
      double Dir[3] = {0.0, 0.0, 0.0};
      // angle swiped by the loop so far
      double Angle = 0.0;
      // estimate for the next segment's length
      static double Length = 1E3;
      // position and magnetic field and next candidate vertex
      double xNew[3] = {0.0,0.0,0.0};
      double BNew[3] = {0.0,0.0,0.0};
      double bNew[3] = {0.0,0.0,0.0};
      double absBNew;
      // new vertex
      cFieldLineVertex *Vertex;
      //.......................................................................

      //increase counter of field lines
      nFieldLine++;

      //rotate to the y=0 plane
      x[0] = pow(x[0]*x[0] + x[1]*x[1], 0.5) * ( (x[0]>0) ? 1 : -1);
      x[1] = 0.0;

      // get magnetic field at the current location
      CPLR::InitInterpolationStencil(x);
      CPLR::GetBackgroundMagneticField(B);
      absB = pow(B[0]*B[0]+B[2]*B[2],0.5);
      b[0] = B[0]/absB; b[1] = 0; b[2] = B[2]/absB;

      // target value of the flux function
      Psi0 = CPLR::GetBackgroundMagneticFluxFunction();

      // check correctness: need to have non-zero components in x-z plane
      if(B[0]*B[0] + B[2]*B[2] < epsB2)
	exit(__LINE__,__FILE__, 
	     "ERROR: magnetic field magnitude is below min resolved value");
      
      //add the initial vertex
      Vertex = FieldLinesAll[nFieldLine-1].Add(x);
      //      FieldLinesAll[nFieldLine-1].SetMagneticField(B,0);
      Vertex -> SetMagneticField(B);
      Vertex -> SetDatum(DatumAtVertexMagneticFluxFunction, Psi0);
      CPLR::GetBackgroundPlasmaVelocity(V);
      Vertex -> SetDatum(DatumAtVertexPlasmaVelocity, V);

      //generate the loop
      //......................................................................
      // new vertices are added until:
      // - 2*Pi angle is swiped AND 
      //   candidate vertex is spatially close to the first one => SUCCESS
      // - 4*Pi angle is swiped => FAILURE
      // - exited the domain => FAILURE
      //......................................................................
      // new point is generated in 2 steps
      // 1) make the first guess, x_0, in the direction of B
      // 2) iterate towards point with same value of FLUX via gradient descent
      //......................................................................
      
      while(fabs(Angle) < 4*Pi){

	// get the first guess
	//------------------------------------------------------------------
	// control inner while loop
	int count;
	const  int countMax = 100;
	count = 0;
	//predictor
	while(true){
	  count++;
	  if(count > countMax)
	    exit(__LINE__, __FILE__,"ERROR: can't generate field-line loop ");
	  // get new location
	  xNew[0] = x[0] + Length * b[0];
	  xNew[2] = x[2] + Length * b[2];

	  // get magnetic field at this location
	  CPLR::InitInterpolationStencil(xNew);
	  CPLR::GetBackgroundMagneticField(BNew);
	  absBNew = pow(BNew[0]*BNew[0]+BNew[2]*BNew[2],0.5);
	  // need to have non-zero field 
	  if(absBNew < epsB)
	    exit(__LINE__,__FILE__, 
		 "ERROR: magnetic field magnitude is below min resolved");
	  bNew[0] = BNew[0]/absBNew; bNew[1] = 0; bNew[2] = BNew[2]/absBNew;
	  // find angle of arc to the new location
	  Arc = fabs(b[0]*bNew[2] - b[2]*bNew[0]);
	  //check condition to exit loop
	  if(Arc < 0.5 * DArc || Arc > DArc)
	    // too small or too large step; factor -> 1 as count grows
	    // to avoid reaching stationary points and infinite loop
	    Length *= 1. + (0.75 * DArc / max(Arc, 0.075*DArc) - 1.) / count;
	  else
	    break;
	}
	
	// corrector
	{
	  Dir[0] = 0.5*(BNew[0]+B[0]);
	  Dir[1] = 0;
	  Dir[2] = 0.5*(BNew[2]+B[2]);
	  double misc = pow(Dir[0]*Dir[0] + Dir[2]*Dir[2], 0.5);
	  Dir[0]/= misc; Dir[2]/= misc;
	  //new vertex location
	  xNew[0] = x[0] +  Length * Dir[0];
	  xNew[2] = x[2] +  Length * Dir[2];
	}	
	//------------------------------------------------------------------

	// iterate towards location with needed value of flux function
	//------------------------------------------------------------------
	//get magnetic field at the candidate location
	CPLR::InitInterpolationStencil(xNew);
	CPLR::GetBackgroundMagneticField(BNew);
	// need to have non-zero field
	if(BNew[0]*BNew[0]+BNew[2]*BNew[2] < epsB2)
	  exit(__LINE__,__FILE__, 
	       "ERROR: magnetic field magnitude is below min resolved");
	// flux function
	Psi = CPLR::GetBackgroundMagneticFluxFunction();

	// correct the initial guess until flux function is close enough
	// to the original value
	//------------------------------------------------------------------
	while(fabs(2*(Psi-Psi0)/(Psi+Psi0)) > 1E-5){
	  double misc = BNew[0]*BNew[0] + BNew[2]*BNew[2];
	  xNew[0]+=-0.001 * (Psi-Psi0) * BNew[2] / misc;
	  xNew[2]+= 0.001 * (Psi-Psi0) * BNew[0] / misc;
	  CPLR::InitInterpolationStencil(xNew);
	  CPLR::GetBackgroundMagneticField(BNew);
	  Psi = CPLR::GetBackgroundMagneticFluxFunction();
	}
	//------------------------------------------------------------------

	// next vertex is found;
	// add it to the field line and update information
	absBNew = pow(BNew[0]*BNew[0]+BNew[2]*BNew[2],0.5);
	bNew[0] = BNew[0]/absBNew; bNew[1] = 0; bNew[2] = BNew[2]/absBNew;

	//housekeeping
	Angle += asin(b[0]*bNew[2] - b[2]*bNew[0]);
	
	//add this vertex
	Vertex = FieldLinesAll[nFieldLine-1].Add(xNew);
	//	FieldLinesAll[nFieldLine-1].SetMagneticField(BNew);
	Vertex -> SetMagneticField(BNew);
	Vertex -> SetDatum(DatumAtVertexMagneticFluxFunction, Psi0);
	CPLR::GetBackgroundPlasmaVelocity(V);
	Vertex -> SetDatum(DatumAtVertexPlasmaVelocity, V);

	//check conditions for completing the loop
	cFieldLineSegment* First=FieldLinesAll[nFieldLine-1].GetFirstSegment();
	cFieldLineSegment* Last =FieldLinesAll[nFieldLine-1].GetLastSegment();
	double Dist = 0.0; //distance between last and first verticies
	double DirFirst[3], DirLast[3];
	First->GetDir(DirFirst); Last->GetDir(DirLast);
	for(int idim=0; idim<DIM; idim++){
	  Dist  += pow(xNew[idim]-xFirst[idim], 2);
	}
	Dist   = pow(Dist, 0.5);
	if( fabs( 0.5*fabs(Angle)/Pi - 1) < 0.1 &&
	    Dist < Length){
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
      exit(__LINE__, __FILE__,"ERROR: can't generate field-line loop ");
    }

    //=========================================================================
    // update field lines 
    void Update(){
      cFieldLineVertex*  Vertex;
      cFieldLineSegment* Segment;
      double dt = PIC::SimulationTime::Get() - TimeLastUpdate;
      TimeLastUpdate += dt;
      for(int iFieldLine=0; iFieldLine<nFieldLine; iFieldLine++){
	int cnt = 0;
	bool done = false;
	for(Vertex = FieldLinesAll[iFieldLine].GetFirstVertex() ;
	    !done;
	    Vertex = Vertex->GetNext(), cnt++){
	    done = Vertex==FieldLinesAll[iFieldLine].GetLastVertex() && cnt!=0;
	  double X[3],B[3],V[3];
	  Vertex->GetX(X);
	  Vertex->GetPlasmaVelocity(V);
	  for(int i=0; i<DIM; i++) X[i] += V[i] * dt;
         #if _PIC_SYMMETRY_MODE_ == _PIC_SYMMETRY_MODE__AXIAL_
	  // rotate to the y=0 plane
	  X[0] = pow(X[0]*X[0]+X[1]*X[1], 0.5);
	  X[1] = 0.0;
         #endif //_PIC_SYMMETRY_MODE_ == _PIC_SYMMETRY_MODE__AXIAL_ 
         #if _PIC_DATAFILE__TIME_INTERPOLATION_MODE_ == _PIC_MODE_ON_
	  // flux function
	  CPLR::InitInterpolationStencil(X);
	  CPLR::GetBackgroundMagneticField(B);
	  double Psi = CPLR::GetBackgroundMagneticFluxFunction();
	  double Psi0;
	  Vertex->GetDatum(DatumAtVertexMagneticFluxFunction, Psi0);

	  // correct the initial guess until flux function is close enough
	  // to the original value
	  //------------------------------------------------------------------
	  while(fabs(2*(Psi-Psi0)/(Psi+Psi0)) > 1E-5){
	    double misc = B[0]*B[0] + B[2]*B[2];
	    X[0]+=-0.001 * (Psi-Psi0) * B[2] / misc;
	    X[2]+= 0.001 * (Psi-Psi0) * B[0] / misc;
	    CPLR::InitInterpolationStencil(X);
	    CPLR::GetBackgroundMagneticField(B);
	    Psi = CPLR::GetBackgroundMagneticFluxFunction();
	  }
	  //------------------------------------------------------------------
         #endif//_PIC_DATAFILE__TIME_INTERPOLATION_MODE_ == _PIC_MODE_ON_
	  Vertex->SetX(X);
	  // update background data
	  PIC::CPLR::InitInterpolationStencil(X);
	  PIC::CPLR::GetBackgroundMagneticField(B);
	  PIC::CPLR::GetBackgroundPlasmaVelocity(V);
	  Vertex->SetDatum(DatumAtVertexMagneticField, B);
	  Vertex->SetDatum(DatumAtVertexPlasmaVelocity,V);
	}
	cnt = 0; done = false;
	for(Segment = FieldLinesAll[iFieldLine].GetFirstSegment();
	    !done; 
	    Segment = Segment->GetNext(), cnt++){
	  done = Segment==FieldLinesAll[iFieldLine].GetLastSegment()&&cnt!=0;
	  Segment->SetVertices(Segment->GetBegin(),
			       Segment->GetEnd());
	}
      }
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

  //---------------------------------------------------------------------------

namespace Mover{
namespace FieldLine{
  // procedure that returns parameters of the guiding center motion
  void GuidingCenterMotion(double& ForceParal, 
			   double& AbsB,
			   int spec,long int ptr,
			   int   iFieldLine,
			   double FieldLineCoord){
    /* function returns guiding center velocity in direction perpendicular
     * to the magnetic field and the force parallel to it
     * for Lorentz force (considered here)
     *
     * v_{guide_perp} = 
     *  E\times B / B^2 + 
     *  \mu/(q\gamma) B\times\nabla B / B^2 +
     *  (p_{\parallel}^2)/(q\gamma m_0) B\times(B\cdot\nabla)B/B^4
     *
     * dp_{\parallel}/dt=
     *  q E_{\parallel} - \mu/\gamma \nabla_{\parallel}B
     *
     * \mu = p_{\perp}^2/(2m_0B)
     *
     * \gamma = 1/\sqrt{1-v^2/c^2}
     *********************************************************************/

    // namespace alias
    namespace FL = PIC::FieldLine;

    double ForceParal_LOC = 0.0;
    double AbsB_LOC       = 0.0;
    
    //#if _FORCE_LORENTZ_MODE_ == _PIC_MODE_ON_
    
    double E[3],B[3],B0[3],B1[3], AbsBDeriv;
    double mu     = PIC::ParticleBuffer::GetMagneticMoment(ptr);
    double q      = PIC::MolecularData::GetElectricCharge(spec);
    FL::FieldLinesAll[iFieldLine].GetMagneticField(B0, (int)FieldLineCoord);
    FL::FieldLinesAll[iFieldLine].GetMagneticField(B,       FieldLineCoord);
    FL::FieldLinesAll[iFieldLine].GetMagneticField(B1, (int)FieldLineCoord+1-1E-7);
    AbsB_LOC   = pow(B[0]*B[0] + B[1]*B[1] + B[2]*B[2], 0.5);
    AbsBDeriv = (pow(B1[0]*B1[0] + B1[1]*B1[1] + B1[2]*B1[2], 0.5) -
		 pow(B0[0]*B0[0] + B0[1]*B0[1] + B0[2]*B0[2], 0.5)) / 
      FL::FieldLinesAll[iFieldLine].GetSegmentLength(FieldLineCoord);
    
    //parallel force
#if _PIC__IDEAL_MHD_MODE_ == _PIC_MODE_ON_
    // in this case E = - V \cross B => E_{\paral} = E*b = 0
    ForceParal_LOC = - mu * AbsBDeriv;
#else
    exit(__LINE__, __FILE__, "not implemented");
#endif//_PIC__IDEAL_MHD_MODE_ == _PIC_MODE_ON_
    
    
    //#endif//_FORCE_LORENTZ_MODE_
    
    ForceParal = ForceParal_LOC;
    AbsB       = AbsB_LOC;
    
  }
  
  // mover itself
  int Mover_SecondOrder(long int ptr, double dtTotal, 
			cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode){
    //aliases
    namespace FL = PIC::FieldLine;
    namespace PB = PIC::ParticleBuffer;


    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *newNode=NULL;
    double dtTemp;
    PIC::ParticleBuffer::byte *ParticleData;
    double AbsBInit=0.0, bInit[3]={0.0};
    double vInit[  3]={0.0}, pInit  =0.0, xInit[  3]={0.0}, gammaInit;
    double AbsBMiddle=0.0, bMiddle[3]={0.0};
    double vMiddle[3]={0.0}, pMiddle=0.0, xMiddle[3]={0.0}, gammaMiddle;
    double vFinal[ 3]={0.0}, pFinal =0.0, xFinal[ 3]={0.0}, gammaFinal;
    double c2 = SpeedOfLight*SpeedOfLight;
    double xminBlock[3],xmaxBlock[3];
    int   iFieldLine = -1;
    double FieldLineCoordInit   = -1.0;
    double FieldLineCoordMiddle = -1.0;
    double FieldLineCoordFinal  = -1.0;
    int idim;
    
    double ForceParalInit=0.0, ForceParalMiddle=0.0;
    double dirInit[3]={0.0,0.0,0.0};
    double vparInit=0.0, vparMiddle=0.0, vparFinal=0.0;

    long int LocalCellNumber;
    int i,j,k,spec;
    
    PIC::Mesh::cDataCenterNode *cell;
    bool MovingTimeFinished=false;
    
    double xmin,xmax;
    
    double misc;
	
    ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
    PIC::ParticleBuffer::GetV(vInit,ParticleData);
    PIC::ParticleBuffer::GetX(xInit,ParticleData);
    spec=PIC::ParticleBuffer::GetI(ParticleData);
    double m0 = PIC::MolecularData::GetMass(spec);
    double mu = PIC::ParticleBuffer::GetMagneticMoment(ptr);
    iFieldLine = PIC::ParticleBuffer::GetFieldLineId(ptr);
    FieldLineCoordInit = PIC::ParticleBuffer::GetFieldLineCoord(ptr);
    static long int nCall=0;
    nCall++;
    
    //  memcpy(xminBlock,startNode->xmin,DIM*sizeof(double));
    //  memcpy(xmaxBlock,startNode->xmax,DIM*sizeof(double));
    
    MovingTimeFinished=true;
    
    // predictor step
#if _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE_ == _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE__ON_
    
    GuidingCenterMotion(ForceParalInit,
			AbsBInit,
			spec,ptr,iFieldLine,FieldLineCoordInit);
    
#endif
    
    FL::FieldLinesAll[iFieldLine].GetSegmentDirection(dirInit, FieldLineCoordInit);
    vparInit = vInit[0]*dirInit[0]+vInit[1]*dirInit[1]+vInit[2]*dirInit[2];
    
    
    dtTemp=dtTotal/2.0;
    // advance coordinates half-step
    FieldLineCoordMiddle = FieldLineCoordInit + 
      dtTemp * vparInit  / FL::FieldLinesAll[iFieldLine].GetSegmentLength(FieldLineCoordInit);
    // advance momentum half-step
    vparMiddle = vparInit + 
      dtTemp * ForceParalInit/m0; 
    
    
    
    // check if a particle has left the domain
    if (FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoordMiddle)==NULL) { 
      //the particle left the computational domain
      int code=_PARTICLE_DELETED_ON_THE_FACE_;
      
      //call the function to processes particles that left the domain
      switch(code){
      case _PARTICLE_DELETED_ON_THE_FACE_:
	PIC::ParticleBuffer::DeleteParticle(ptr);
	return _PARTICLE_LEFT_THE_DOMAIN_;
	break;
      default:
	exit(__LINE__,__FILE__,"Error: not implemented");
      }
      
    }
    
    
    // corrector step
    
#if _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE_ == _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE__ON_
    GuidingCenterMotion(ForceParalMiddle,
			AbsBMiddle,
			spec,ptr,iFieldLine,FieldLineCoordMiddle);
#endif
    

    //    FL::FieldLinesAll[iFieldLine].GetSegmentDirection(dirInit, FieldLineCoordMiddle);
    
    



  // advance coordinates full-step
    FieldLineCoordFinal = FieldLineCoordInit + 
      dtTotal * vparMiddle / FL::FieldLinesAll[iFieldLine].GetSegmentLength(FieldLineCoordMiddle);
    // advance momentum full-step
    vparFinal =vparInit + dtTotal * ForceParalMiddle/m0; 

   
    //advance the particle's position and velocity
    //interaction with the faces of the block and internal surfaces
  
    if (FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoordFinal)==NULL) {
      
      //the particle left the computational domain
      int code=_PARTICLE_DELETED_ON_THE_FACE_;
      
      //call the function that process particles that leaved the coputational domain
      switch(code){
      case _PARTICLE_DELETED_ON_THE_FACE_:
	PIC::ParticleBuffer::DeleteParticle(ptr);
	return _PARTICLE_LEFT_THE_DOMAIN_;
      break;
      default:
      exit(__LINE__,__FILE__,"Error: not implemented");
      }
      
    }

    FL::FieldLinesAll[iFieldLine].GetSegmentDirection(vFinal, FieldLineCoordFinal);
    vFinal[0]*=vparFinal;vFinal[1]*=vparFinal;vFinal[2]*=vparFinal;
    FL::FieldLinesAll[iFieldLine].GetCartesian(xFinal, FieldLineCoordFinal);
    
    newNode=PIC::Mesh::Search::FindBlock(xFinal);
    

    PIC::ParticleBuffer::SetV(vFinal,ParticleData);
    PIC::ParticleBuffer::SetX(xFinal,ParticleData);
    PB::SetFieldLineCoord(FieldLineCoordFinal, ParticleData);
  
    //adjust the value of 'startNode'
    startNode=newNode;
    
    if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(xFinal,i,j,k,startNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cell where the particle is located");
    
    
#if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
    PIC::Mesh::cDataBlockAMR *block=startNode->block;
    long int tempFirstCellParticle=block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
    
    PIC::ParticleBuffer::SetNext(tempFirstCellParticle,ParticleData);
    PIC::ParticleBuffer::SetPrev(-1,ParticleData);
    
    if (tempFirstCellParticle!=-1) PIC::ParticleBuffer::SetPrev(ptr,tempFirstCellParticle);
    block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=ptr;
#else
    exit(__LINE__,__FILE__,"Error: the block need to be generalized to the case when _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_");
#endif
    
    
    return _PARTICLE_MOTION_FINISHED_;
  }

}
}
  //namespace Mover -----------------------------------------------------------

}

