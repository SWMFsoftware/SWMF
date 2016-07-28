!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine PW_upper_heat_conduction
  use ModCommonVariables
      use ModCommonPlanet,ONLY: HLPion1,HLPion2,HLPion3,HLPE

  HeatCon_GI(nDim+1,Ion1_)=&
       HLPion1*(State_GV(nDim+1,RhoH3_)/State_GV(nDim+1,RhoE_))*State_GV(nDim+1,Th3_)**2.5
  HeatCon_GI(nDim+1,nIon)=&
       HLPE*State_GV(nDim+1,Te_)**2.5
  HeatCon_GI(nDim+1,Ion2_)=&
       HLPion2*(State_GV(nDim+1,RhoH_)/State_GV(nDim+1,RhoE_))*State_GV(nDim+1,Th_)**2.5
  HeatCon_GI(nDim+1,Ion3_)=&
       HLPion3*(State_GV(nDim+1,RhoH2_)/State_GV(nDim+1,RhoE_))*State_GV(nDim+1,Th2_)**2.5


end subroutine PW_upper_heat_conduction
