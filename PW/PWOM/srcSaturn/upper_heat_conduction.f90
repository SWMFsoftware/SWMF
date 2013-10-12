!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine PW_upper_heat_conduction
  use ModCommonVariables

  HeatCon_GI(nDim+1,Ion1_)=&
       HLPO*(State_GV(nDim+1,RhoO_)/State_GV(nDim+1,RhoE_))*State_GV(nDim+1,To_)**2.5
  HeatCon_GI(nDim+1,nIon)=&
       HLPE*State_GV(nDim+1,Te_)**2.5
  HeatCon_GI(nDim+1,Ion2_)=&
       HLPH*(State_GV(nDim+1,RhoH_)/State_GV(nDim+1,RhoE_))*State_GV(nDim+1,Th_)**2.5


end subroutine PW_upper_heat_conduction
