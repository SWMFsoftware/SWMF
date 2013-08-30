
#include "pic.h"

//int PIC::nTotalSpecies=0;
int PIC::ThisThread=0,PIC::nTotalThreads=1;

//the list containing the functions used to exchange the run time execution statistics
vector<PIC::fExchangeExecutionStatistics> PIC::ExchangeExecutionStatisticsFunctions;

//the list contains the functions used for user defined sampling procedures
vector<PIC::IndividualModelSampling::fRequestSamplingData> PIC::IndividualModelSampling::RequestSamplingData;
vector<PIC::IndividualModelSampling::fSamplingProcedure> PIC::IndividualModelSampling::SamplingProcedure;
vector<PIC::IndividualModelSampling::fPrintVariableList> PIC::IndividualModelSampling::PrintVariableList;
vector<PIC::IndividualModelSampling::fInterpolateCenterNodeData> PIC::IndividualModelSampling::InterpolateCenterNodeData;
vector<PIC::IndividualModelSampling::fPrintSampledData> PIC::IndividualModelSampling::PrintSampledData;
vector<PIC::IndividualModelSampling::fRequestStaticCellData> PIC::IndividualModelSampling::RequestStaticCellData;

//generic particle transformation
//PIC::ChemicalReactions::GenericParticleTranformation::fTransformationIndicator *PIC::ChemicalReactions::GenericParticleTranformation::TransformationIndicator=NULL;
//PIC::ChemicalReactions::GenericParticleTranformation::fTransformationProcessor *PIC::ChemicalReactions::GenericParticleTranformation::TransformationProcessor=NULL;

//execution alarm
bool PIC::Alarm::AlarmInitialized=false,PIC::Alarm::WallTimeExeedsLimit=false;
double PIC::Alarm::StartTime=0.0,PIC::Alarm::RequestedExecutionWallTime=0.0;

//the file descriptor and the prefix for output of the diagnostic
FILE* PIC::DiagnospticMessageStream=stdout;
char PIC::DiagnospticMessageStreamName[_MAX_STRING_LENGTH_PIC_]="stdout";

//the directory for output files
char PIC::OutputDataFileDirectory[_MAX_STRING_LENGTH_PIC_]=".";
