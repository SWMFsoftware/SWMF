//$Id$
//header of the functions used in the test

double TheoreticalLifeTime(double *x,int spec,long int ptr,bool &PhotolyticReactionAllowedFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
int ReactionProcessor(double *xInit,double *xFinal,double *vFinal,long int ptr,int &spec,PIC::ParticleBuffer::byte *ParticleData, double ReactionTimeInterval,double ReactionTimeIntervalStarts,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);

