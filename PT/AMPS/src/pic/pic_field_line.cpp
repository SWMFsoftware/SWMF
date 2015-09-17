//$Id$

//implementation of methods in namespace FieldLine

#include "pic.h"

cStack<PIC::FieldLine::cFieldLineVertex>  PIC::FieldLine::VerticesAll;
cStack<PIC::FieldLine::cFieldLineSegment> PIC::FieldLine::SegmentsAll;
PIC::FieldLine::cFieldLine *PIC::FieldLine::FieldLinesAll = NULL;


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

  for(int iFieldLine=0; iFieldLine<nFieldLineMax; iFieldLine++){
    fprintf(fout,"ZONE T=\"Field-line %i\" F=POINT\n",iFieldLine);
    FieldLinesAll[iFieldLine].Output(fout, true);
  }

  fclose(fout);

}
