!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine PW_upper_heat_conduction

  use ModCommonVariables
  
  HeatCon_GI(nDim+1,Ion1_)=&
       HLPO*(State_GV(nDim+1,RhoO_)/State_GV(nDim+1,RhoE_))*State_GV(nDim+1,To_)**2.5
  HeatCon_GI(nDim+1,nIon)=HLPE*State_GV(nDim+1,Te_)**2.5
  HeatCon_GI(nDim+1,Ion2_)=&
       HLPH*(State_GV(nDim+1,RhoH_)/State_GV(nDim+1,RhoE_))*State_GV(nDim+1,Th_)**2.5
  THLP=SQRT(State_GV(nDim+1,Th_))*(1.-0.047*ALOG10(State_GV(nDim+1,Th_)))**2
  TKLP=State_GV(nDim+1,RhoO_)/(State_GV(nDim+1,To_)+16.*State_GV(nDim+1,Th_))**1.5
  HeatCon_GI(nDim+1,Ion2_)=&
       HeatCon_GI(nDim+1,Ion2_)/(1.+(0.7692*(CZHN2+CZHO2)+1.0962*CZHO*THLP)/CZHOX &
       /TKLP)
  THLP=0.5*(XTNMAX+State_GV(nDim+1,The_))
  THLP=SQRT(THLP)*(1.-0.093*ALOG10(THLP))**2
  TKLP=State_GV(nDim+1,RhoO_)/(State_GV(nDim+1,To_)+4.*State_GV(nDim+1,The_))**1.5
  TLLP=State_GV(nDim+1,RhoH_)/(State_GV(nDim+1,The_)+4.*State_GV(nDim+1,Th_))**1.5
  HeatCon_GI(nDim+1,Ion3_)=&
       HLPHE*(State_GV(nDim+1,RhoHe_)/Mass_I(Ion3_))*State_GV(nDim+1,The_)
  HeatCon_GI(nDim+1,Ion3_)=&
       HeatCon_GI(nDim+1,Ion3_)/(0.99*CZHEN2+0.99*CZHEO2+1.02*CZHEO+1.48*CZHEHE*&
       THLP+2.22*CZHEH+1.21*CZHEOX*TKLP+2.23*CZHEHD*TLLP)      
end subroutine PW_upper_heat_conduction
