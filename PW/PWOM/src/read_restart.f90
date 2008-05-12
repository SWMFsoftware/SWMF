
subroutine PW_read_restart
  use ModNumConst,    ONLY: cDegToRad
  use ModIoUnit,      ONLY: io_unit_new,UnitTmp_
  use ModPwom,        ONLY: nLine, Time, GeoMagLat_I, GeoMagLon_I, nStep, &
                            ThetaLine_I, PhiLine_I , State_CVI, nAlt, &
                            NameRestartIn, Dt_I
  use ModCommonPlanet,ONLY: nIon,iRho_I,iU_I,iP_I,iT_I,nVar
  use ModCommonVariables,ONLY: DrBnd

  implicit none
  real   :: ddt1, xxx,Dx, DrFile, Loc,Altmin
  real,allocatable :: StateFile_CVI(:,:,:),Alt_C(:)
  integer:: iError, iIon, iLine,i, nAltFile,iAlt,iMin,iMax

!------------------------------------------------------------------------------
  if (.not.allocated(StateFile_CVI)) &
       allocate(StateFile_CVI(nAlt,nVar,nLine))
  if (.not.allocated(Alt_C)) &
       allocate(Alt_C(nAlt))
  !read in restart data for each line
  do iLine=1,nLine
     OPEN(UNIT=UnitTmp_, FILE=NameRestartIn(iLine), STATUS='OLD')
     READ (UnitTmp_,*) Time,Dt_I(iLine),nAltFile, nStep
     READ (UnitTmp_,*) GeoMagLat_I(iLine),GeoMagLon_I(iLine)
     
     ThetaLine_I (iLine) = (90.0-GeoMagLat_I(iLine)) * cDegToRad
     PhiLine_I   (iLine) = GeoMagLon_I(iLine)        * cDegToRad
     
     do iIon=1,nIon
        READ (UnitTmp_,*) &
             (Alt_C(i),StateFile_CVI(i,iU_I(iIon),iLine),&
             StateFile_CVI(i,iP_I(iIon),iLine),&
             StateFile_CVI(i,iRho_I(iIon),iLine),StateFile_CVI(i,iT_I(iIon),iLine),&
             i=1,nAltFile)
     enddo
     
     CLOSE(UNIT=UnitTmp_)
  enddo

  !check that restart data has same resolution as called for in the simulation
  DrFile = Alt_C(2)-Alt_C(1)
  Altmin = Alt_C(1)
  if (DrFile >= DrBnd - 1.0e-5 .and. DrFile <= DrBnd + 1.0e-5 .and.& 
       nAltFile == nAlt) then
     !Restart File and simulation use same grid.  
     State_CVI(:,:,:) = StateFile_CVI(:,:,:)
  else
     !Restart File and simulation use different grid. Interpolate restart 
     !solution onto simulation grid
     do iLine=1,nLine
        do iAlt=1,nAlt
           Loc  = (DrBnd*iAlt)/DrFile+1.0
           iMin = floor(Loc)
           iMax = ceiling(Loc)
           Dx   = Loc-iMin 
           
           if (iMax <= nAltFile) then
              State_CVI(iAlt,1:nVar,iLine) = &
                   Dx * StateFile_CVI(iMax,1:nVar,iLine) &
                   + (1.0 - Dx) * StateFile_CVI(iMin,1:nVar,iLine)
           else
              State_CVI(iAlt,1:nVar,iLine) = &
                   StateFile_CVI(nAltFile,1:nVar,iLine)
           endif
        enddo
     enddo
  endif
     
  deallocate(StateFile_CVI)
  deallocate(Alt_C)
end subroutine PW_read_restart
