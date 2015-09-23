//$Id$

//implementation of methods in namespace FieldLine

#include "pic.h"

cStack<PIC::FieldLine::cFieldLineVertex>  PIC::FieldLine::VerticesAll;
cStack<PIC::FieldLine::cFieldLineSegment> PIC::FieldLine::SegmentsAll;
PIC::FieldLine::cFieldLine *PIC::FieldLine::FieldLinesAll = NULL;

long int PIC::FieldLine::nFieldLine=0;

bool PIC::FieldLine::cFieldLine::is_broken(){
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

//=============================================================================

void PIC::FieldLine::cFieldLine::Add(double *xIn){
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

//=============================================================================

void PIC::FieldLine::cFieldLine::Output(FILE* fout, bool GeometryOnly=false){
  
  if(!GeometryOnly)
    exit(__LINE__, __FILE__,"Not implemented for multiple processors");

  cFieldLineVertex *Vertex=FirstVertex;
  for(int iVertex=0; iVertex<=nSegment; iVertex++){
    double *x;
    Vertex->GetX(x);
    for(int idim=0; idim<DIM; idim++) {
      fprintf(fout, "%e ", x[idim]);
    }
    fprintf(fout,"\n");
    Vertex = Vertex->GetNext();
  }
}

//=============================================================================

void PIC::FieldLine::Init(){

  if(PIC::nTotalThreads > 1)
    exit(__LINE__, __FILE__,"Not implemented for multiple processors");

  //allocate container for field lines
  FieldLinesAll = new cFieldLine [nFieldLineMax];

}

//=============================================================================

void PIC::FieldLine::Output(char* fname, bool GeometryOnly){

  FILE* fout;
  fout = fopen(fname,"w");

  fprintf(fout, "TITLE=\"Field line geometry\"");
  
#if DIM == 3 
  fprintf(fout,"VARIABLES=\"x\",\"y\",\"z\"\n");
#else 
  exit(__LINE__,__FILE__,"not implemented");
#endif

  for(int iFieldLine=0; iFieldLine<nFieldLine; iFieldLine++){
    fprintf(fout,"ZONE T=\"Field-line %i\" F=POINT\n",iFieldLine);
    FieldLinesAll[iFieldLine].Output(fout, true);
  }

  fclose(fout);

}

//=============================================================================

void PIC::FieldLine::InitSimpleParkerSpiral(double *xStart){

#if DIM != 3
  exit(__LINE__,__FILE__,"Implemetned only for 3D case");
#endif

  //position of the Sun
  double xSun[DIM]={0};
  //relative position to Sun
  double xRel[DIM]={xStart[0]-xSun[0],xStart[1]-xSun[1],xStart[2]-xSun[2]};
  //distance to Sun
  double R2Sun=0;
  //distance to sun in equatorial plane
  double r2Sun=0;

  //velocity of solar wind
  const double SpeedSolarWind = 4.0E5;
  //rate of solar rotation at equator (rad/sec)
  const double Omega = 2.97211E-6;
  //spatial step
  const double Ds=1E9;
  int nSegment = 1000;

  if(nFieldLine == nFieldLineMax)
    exit(__LINE__,__FILE__,"ERROR: reached limit for number of field lines");
  
  //increase counter of field lines
  nFieldLine++;

  //start position
  double x[DIM]={xStart[0],xStart[1],xStart[2]};

  //original latitude 
  R2Sun = pow(xRel[0]*xRel[0]+xRel[1]*xRel[1]+xRel[2]*xRel[2], 0.5);
  double costheta = xRel[2] / R2Sun;
  double sintheta = pow(1 - costheta*costheta,0.5); // >=0 always

  //direction
  double a[DIM]={0,0,costheta};

  FieldLinesAll[nFieldLine-1].Add(xStart);

  for(int iSegment=1;iSegment<nSegment; iSegment++){
    //distances to sun
    xRel[0]=x[0]-xSun[0];xRel[1]=x[1]-xSun[1];xRel[2]=x[2]-xSun[2];
    R2Sun = pow(xRel[0]*xRel[0]+xRel[1]*xRel[1]+xRel[2]*xRel[2], 0.5);
    r2Sun = pow(xRel[0]*xRel[0]+xRel[1]*xRel[1], 0.5);
    //angle in the equatorial plane
    double cosphi = xRel[0]/r2Sun, sinphi=xRel[1]/r2Sun;
    // equatorial components of vector a
    a[0] = cosphi + (Omega/SpeedSolarWind)*R2Sun*sinphi;
    a[1] = sinphi - (Omega/SpeedSolarWind)*R2Sun*cosphi;

    //normalize vector a
    double alen = pow(a[0]*a[0]+a[1]*a[1],0.5);
    a[0] *= sintheta/alen; a[1] *= sintheta/alen;
    
    // next vertex of the spiral
    x[0]+=a[0]*Ds; x[1]+=a[1]*Ds; x[2]+=a[2]*Ds;
    FieldLinesAll[nFieldLine-1].Add(x);
  }
  
}
