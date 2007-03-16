subroutine PW_write_restart(&
     nAlt,RAD,GmLat,GmLong,Time,DT,nStep,NameRestart,     &
     dOxyg_C, uOxyg_C, pOxyg_C, TOxyg_C,         &
     dHel_C , uHel_C , pHel_C , THel_C ,         &
     dHyd_C , uHyd_C , pHyd_C , THyd_C ,         &
     dElect_C, uElect_C, pElect_C, TElect_C)

  !Description: Write out restart file for PW

  use ModIoUnit, ONLY: UnitTmp_
  implicit none
  integer, intent(in) :: nAlt,nStep
  real   , intent(in) :: GmLat,GmLong,Time,DT
  real   , intent(in) :: RAD(nAlt)
  real   , intent(in) :: &
       dOxyg_C(nAlt), uOxyg_C(nAlt), pOxyg_C(nAlt), TOxyg_C(nAlt), &
       dHel_C(nAlt) , uHel_C(nAlt) , pHel_C(nAlt) , THel_C(nAlt) , &
       dHyd_C(nAlt) , uHyd_C(nAlt) , pHyd_C(nAlt) , THyd_C(nAlt) , &
       dElect_C(nAlt), uElect_C(nAlt), pElect_C(nAlt), TElect_C(nAlt)
  
  character*100,intent(in)   :: NameRestart
  
  integer :: K
  !_____________________________________________________________________________
  
  open(UnitTmp_, FILE=NameRestart)
  
  WRITE (UnitTmp_,*) TIME,DT,NSTEP
  Write (UnitTmp_,*) GMLAT, GMLONG
  
  WRITE (UnitTmp_,2002)&
       (RAD(K),uOxyg_C(K) ,pOxyg_C(K) ,dOxyg_C(K) ,TOxyg_C(K) ,K=1,nAlt)        
  WRITE (UnitTmp_,2002)&
       (RAD(K),uHel_C(K)  ,pHel_C(K)  ,dHel_C(K)  ,THel_C(K)  ,K=1,nAlt)         
  WRITE (UnitTmp_,2002)&
       (RAD(K),uHyd_C(K)  ,pHyd_C(K)  ,dHyd_C(K)  ,THyd_C(K)  ,K=1,nAlt)        
  WRITE (UnitTmp_,2002)&
       (RAD(K),uELECT_C(K),pELECT_C(K),dELECT_C(K),TELECT_C(K),K=1,nAlt)
  close(UnitTmp_)
  
!2001 format(2(1PE16.6),I10)
2002 format(5(1PE25.16))


end subroutine PW_write_restart
