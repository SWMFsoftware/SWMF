module ModCrcmPlot

  implicit none

  private ! except
  public :: Crcm_plot, Crcm_plot_fls, crcm_plot_log
  character(len=5),  public    :: TypePlot   = 'ascii'
  logical,           public    :: DoSavePlot = .false.
  logical,           public    :: DoSaveFlux = .false.
  logical,           public    :: DoSaveLog = .false.
  logical,           public    :: UseSeparatePlotFiles = .false.
  real,              public    :: DtOutput   = 10.0
  real,              public    :: DtLogout   = 10.0
  

  character(len=*), parameter :: NameHeader = 'CRCM output'

contains
  subroutine Crcm_plot(nLat,nLon, X_C,Y_C,Pressure_IC,PressurePar_IC,&
       PressureHot_IC,PparHot_IC, Den_IC, Beq_C,Volume_C,Potential_C,FAC_C, &
       Time,Dt)

    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModCrcmRestart, ONLY: IsRestart
    use ModCrcmPlanet,   ONLY: nspec,NamePlotVar,iPplot_I,iPparplot_I, &
                               iPhotplot_I,iPparhotplot_I, iNplot_I,   &
                               Beq_,Vol_,Pot_,FAC_,nVar
    use ModCrcmGrid,   ONLY: PhiIono_C => phi, LatIono_C => xlatr
    use ModFieldTrace, ONLY: iba
    integer, intent(in) :: nLat, nLon
    real,    intent(in) :: X_C(nLat,nLon), Y_C(nLat,nLon), Time, Dt
    real,    intent(in) :: Pressure_IC(nspec,nLat,nLon), &
                           PressurePar_IC(nspec,nLat,nLon), &
                           Den_IC(nspec,nLat,nLon), & 
                           Beq_C(nLat,nLon),Volume_C(nLat,nLon),   &
                           Potential_C(nLat,nLon), &
                           PressureHot_IC(nspec,nLat,nLon), &
                           PparHot_IC(nspec,nLat,nLon), &
                           FAC_C(nLat,nLon)
    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)
    real, allocatable   :: CoordIono_DII(:,:,:)
    integer             :: iLat,iLon,iSpecies
    integer, parameter  :: x_=1, y_=2, nDim=2
    real                :: Theta, Phi
    character(len=20)   :: NamePlotEq  = 'IM/plots/CRCMeq.outs'
    character(len=22)   :: NamePlotIono= 'IM/plots/CRCMiono.outs'
    character(len=6)    :: TypePosition  ! 'rewind' or 'append'
    real, parameter     :: Gamma = 5./3., rBody = 1.0
    logical             :: IsFirstCall = .true.
    !--------------------------------------------------------------------------
    
    allocate(Coord_DII(nDim,nLat,nLon+1), CoordIono_DII(nDim,nLat,nLon+1), &
         PlotState_IIV(nLat,nLon+1,nVar))

    PlotState_IIV = 0.0
    Coord_DII     = 0.0
    CoordIono_DII = 0.0
    !Set Coords
    Coord_DII(x_,:, 1:nLon) = X_C(:,1:nLon)
    Coord_DII(y_,:, 1:nLon) = Y_C(:,1:nLon)
    
    do iLon = 1,nLon
       do iLat = 1,nLat
          CoordIono_DII(x_,iLat, iLon) = &
               cos(LatIono_C(iLat))*cos(PhiIono_C(iLon))
          CoordIono_DII(y_,iLat, iLon) = &
               cos(LatIono_C(iLat))*sin(PhiIono_C(iLon))
       enddo
    enddo

    !fill ghost cells of Coords
    Coord_DII(x_,:, nLon+1) = X_C(:,1)
    Coord_DII(y_,:, nLon+1) = Y_C(:,1)
    
    CoordIono_DII(x_,:, nLon+1) = CoordIono_DII(x_,:, 1)
    CoordIono_DII(y_,:, nLon+1) = CoordIono_DII(y_,:, 1)

    !Set plot data
    do iSpecies = 1,nspec
       do iLon=1,nLon
       PlotState_IIV(1:iba(iLon),iLon,iPplot_I(iSpecies+1))= &
            Pressure_IC(iSpecies,1:iba(iLon),iLon) 
       PlotState_IIV(1:iba(iLon),iLon,iPparplot_I(iSpecies+1))= &
            PressurePar_IC(iSpecies,1:iba(iLon),iLon) 
       PlotState_IIV(1:iba(iLon),iLon,iPhotplot_I(iSpecies+1))= &
            PressureHot_IC(iSpecies,1:iba(iLon),iLon)
       PlotState_IIV(1:iba(iLon),iLon,iPparhotplot_I(iSpecies+1))= &
            PparHot_IC(iSpecies,1:iba(iLon),iLon)
       PlotState_IIV(1:iba(iLon),iLon,iNplot_I(iSpecies+1))= &
            Den_IC(iSpecies,1:iba(iLon),iLon)  
       PlotState_IIV(1:iba(iLon),iLon,iPplot_I(1))= &
            PlotState_IIV(1:iba(iLon),iLon,iPplot_I(1))&
            +Pressure_IC(iSpecies,1:iba(iLon),iLon) 
       PlotState_IIV(1:iba(iLon),iLon,iPparplot_I(1))= &
            PlotState_IIV(1:iba(iLon),iLon,iPparplot_I(1))&
            +PressurePar_IC(iSpecies,1:iba(iLon),iLon) 
       PlotState_IIV(1:iba(iLon),iLon,iPhotplot_I(1))= &
            PlotState_IIV(1:iba(iLon),iLon,iPhotplot_I(1))&
            + PressureHot_IC(iSpecies,1:iba(iLon),iLon)
       PlotState_IIV(1:iba(iLon),iLon,iPparhotplot_I(1))= &
            PlotState_IIV(1:iba(iLon),iLon,iPparhotplot_I(1))&
            + PparHot_IC(iSpecies,1:iba(iLon),iLon)
       PlotState_IIV(1:iba(iLon),iLon,iNplot_I(1))= &
            PlotState_IIV(1:iba(iLon),iLon,iNplot_I(1))+Den_IC(iSpecies,1:iba(iLon),iLon)  
    end do
 end do
    do iLon=1,nLon
       PlotState_IIV(1:iba(iLon),iLon,Beq_) = Beq_C   (1:iba(iLon),iLon)    
       PlotState_IIV(1:iba(iLon),iLon,Vol_) = Volume_C(1:iba(iLon),iLon)    
       PlotState_IIV(1:iba(iLon),iLon,Pot_) = Potential_C(1:iba(iLon),iLon)    
       PlotState_IIV(1:iba(iLon),iLon,FAC_) = FAC_C   (1:iba(iLon),iLon)    
    enddo
    !fill ghost cells of plot data
    PlotState_IIV(:,nLon+1,iPplot_I(1))= &
         PlotState_IIV(:,1,iPplot_I(1))
    PlotState_IIV(:,nLon+1,iPparplot_I(1))= &
         PlotState_IIV(:,1,iPparplot_I(1))
    PlotState_IIV(:,nLon+1,iPhotplot_I(1))= &
         PlotState_IIV(:,1,iPhotplot_I(1))
    PlotState_IIV(:,nLon+1,iPparhotplot_I(1))= &
         PlotState_IIV(:,1,iPparhotplot_I(1))
    PlotState_IIV(:,nLon+1,iNplot_I(1))= &
         PlotState_IIV(:,1,iNplot_I(1))
    do iSpecies = 1,nspec
       PlotState_IIV(:,nLon+1,iPplot_I(iSpecies+1))= &
            PlotState_IIV(:,1,iPplot_I(iSpecies+1))
       PlotState_IIV(:,nLon+1,iPparplot_I(iSpecies+1))= &
            PlotState_IIV(:,1,iPparplot_I(iSpecies+1))
       PlotState_IIV(:,nLon+1,iPhotplot_I(iSpecies+1))= &
            PlotState_IIV(:,1,iPhotplot_I(iSpecies+1))
       PlotState_IIV(:,nLon+1,iPparhotplot_I(iSpecies+1))= &
            PlotState_IIV(:,1,iPparhotplot_I(iSpecies+1))
       PlotState_IIV(:,nLon+1,iNplot_I(iSpecies+1))= &
            PlotState_IIV(:,1,iNplot_I(iSpecies+1))  
    end do
    PlotState_IIV(:,nLon+1,Beq_) = Beq_C   (:,1)    
    PlotState_IIV(:,nLon+1,Vol_) = Volume_C(:,1)    
    PlotState_IIV(:,nLon+1,Pot_) = Potential_C(:,1)    
    PlotState_IIV(:,nLon+1,FAC_) = FAC_C(:,1)    

    TypePosition = 'append'
    if(IsFirstCall .and. .not. IsRestart) TypePosition = 'rewind'
    IsFirstCall = .false.

    !equatorial plot
    call save_plot_file(NamePlotEq, TypePositionIn=TypePosition,          &
         TypeFileIn=TypePlot,StringHeaderIn = NameHeader,                 &
         NameVarIn = NamePlotVar, nStepIn= nint(Time/Dt),TimeIn=Time,     &
         nDimIn=2,CoordIn_DII=Coord_DII,                                  &
         VarIn_IIV = PlotState_IIV, ParamIn_I = (/Gamma, rBody/))

    ! ionospheric plot
    call save_plot_file(NamePlotIono, TypePositionIn=TypePosition,        &
         TypeFileIn=TypePlot,StringHeaderIn = NameHeader,                 &
         NameVarIn = NamePlotVar, nStepIn= nint(Time/Dt),TimeIn=Time,     &
         nDimIn=2,CoordIn_DII=CoordIono_DII,                              &
         VarIn_IIV = PlotState_IIV, ParamIn_I = (/Gamma, rBody/))
    
    deallocate(Coord_DII, CoordIono_DII, PlotState_IIV)

  end subroutine Crcm_plot
  !============================================================================

  subroutine Crcm_plot_fls(rc,flux,time)
    use ModIoUnit,    ONLY: UnitTmp_
    use ModCrcmGrid,  ONLY:nLat=>np, nLon=>nt, nEnergy=>neng, nPitchAng=>npit,&
                           energy,sinAo,xlat,xmlt,Ebound
    use ModCrcmPlanet,ONLY:nSpecies=>nspec
    use ModFieldTrace,ONLY:ro,bo,xmlto,irm
    use ModCrcmRestart, ONLY: IsRestart
    use ModImTime,    ONLY:iCurrentTime_I
    real, intent(in) :: rc,flux(nSpecies,nLat,nLon,nEnergy,nPitchAng),time
    
    real          :: parmod(1:10)=0.0,lat,ro1,xmlt1,bo1
    integer       :: iLat,iLon,k,m,n,i,nprint
    logical, save :: IsFirstCall = .true.
    character(len=13):: outnameSep
    !--------------------------------------------------------------------------
    nprint=ifix(time/DtOutput)
    write(outnameSep,"(i4.4,i2.2,i2.2,a,i2.2,i2.2)") & 
         iCurrentTime_I(1),iCurrentTime_I(2),iCurrentTime_I(3),'_',&
         iCurrentTime_I(4),iCurrentTime_I(5)
    
    
    do n=1,nSpecies
       If (UseSeparatePlotFiles) then
          
          if (n==1) &
               open(unit=UnitTmp_,file='IM/plots/'//outnameSep//'_h.fls',&
               status='unknown')
          if (n==2 .and. n /= nSpecies) &
               open(unit=UnitTmp_,file='IM/plots/'//outnameSep//'_o.fls',&
               status='unknown')
          if (n==3 .and. n /= nSpecies) &
               open(unit=UnitTmp_,file='IM/plots/'//outnameSep//'_he.fls',&
               status='unknown')
          if (n==nSpecies) &
               open(unit=UnitTmp_,file='IM/plots/'//outnameSep//'_e.fls',&
               status='unknown')
          write(UnitTmp_,"(f10.6,5i6,6x,'! rc in Re,nr,ip,je,ig,ntime')") &
               rc,nLat-1,nLon,nEnergy,nPitchAng,nprint
          write(UnitTmp_,'(6f9.3)') (energy(k),k=1,nEnergy)
          !write(UnitTmp_,'(6f9.3)') (Ebound(k),k=1,nEnergy+1)
          write(UnitTmp_,'(6f9.5)') (sinAo(m),m=1,nPitchAng)
          write(UnitTmp_,'(10f8.3)') (xlat(i),i=2,nLat)
       else
          if (IsFirstCall .and. .not. IsRestart) then
             if (n==1) &
                  open(unit=UnitTmp_,file='IM/plots/CrcmFlux_h.fls',&
                  status='unknown')
             if (n==2 .and. n /= nSpecies) &
                  open(unit=UnitTmp_,file='IM/plots/CrcmFlux_o.fls',&
                  status='unknown')
             if (n==3 .and. n /= nSpecies) &
                  open(unit=UnitTmp_,file='IM/plots/CrcmFlux_he.fls',&
                  status='unknown')
             if (n==nSpecies) &
                  open(unit=UnitTmp_,file='IM/plots/CrcmFlux_e.fls',&
                  status='unknown')
             write(UnitTmp_,"(f10.6,5i6,6x,'! rc in Re,nr,ip,je,ig,ntime')") &
                  rc,nLat-1,nLon,nEnergy,nPitchAng,nprint
             write(UnitTmp_,'(6f9.3)') (energy(k),k=1,nEnergy)
             !write(UnitTmp_,'(6f9.3)') (Ebound(k),k=1,nEnergy+1)
             write(UnitTmp_,'(6f9.5)') (sinAo(m),m=1,nPitchAng)
             write(UnitTmp_,'(10f8.3)') (xlat(i),i=2,nLat)
          else
             if (n==1) &
                  open(unit=UnitTmp_,file='IM/plots/CrcmFlux_h.fls',&
                  status='old', position='append')
             if (n==2 .and. n /= nSpecies)&
                  open(unit=UnitTmp_,file='IM/plots/CrcmFlux_o.fls',&
                  status='old', position='append')
             if (n==3 .and. n /= nSpecies)&
                  open(unit=UnitTmp_,file='IM/plots/CrcmFlux_he.fls',&
                  status='old', position='append')
             if (n==nSpecies)&
                  open(unit=UnitTmp_,file='IM/plots/CrcmFlux_e.fls', &
                  status='old', position='append')
          endif
       endif
       write(UnitTmp_,'(f8.3,10f9.2,"    ! hour,  parmod")') &
            time/3600.,parmod(1:10)
       do iLat=2,nLat             ! Write fluxes @ fixed E & y grids
          do iLon=1,nLon
             lat=xlat(iLat)
             if (iLat.gt.irm(iLon)) lat=xlat(irm(iLon))
             ro1=ro(iLat,iLon)
             if (iLat.gt.irm(iLon)) ro1=ro(irm(iLon),iLon)
             bo1=bo(iLat,iLon)
             if (iLat.gt.irm(iLon)) bo1=bo(irm(iLon),iLon)
             xmlt1=xmlto(iLat,iLon)
             if (iLat.gt.irm(iLon)) xmlt1=xmlto(irm(iLon),iLon)
             write(UnitTmp_,'(f7.2,f6.1,2f8.3,1pe11.3)') &
                  lat,xmlt(iLon),ro1,xmlt1,bo1
             do k=1,nEnergy
                write(UnitTmp_,'(1p,12e11.3)') &
                     (flux(n,iLat,iLon,k,m),m=1,nPitchAng)
             enddo
          enddo
       enddo
       close(UnitTmp_)
    enddo
    IsFirstCall=.false.
    

  end subroutine Crcm_plot_fls

  !============================================================================
  subroutine crcm_plot_log(Time)
    use ModCrcmPlanet,  ONLY:nSpecies=>nspec,NamePlotVarLog
    use ModCrcmRestart, ONLY: IsRestart
    use ModIoUnit,      ONLY: UnitTmp_
    use ModCRCM,        ONLY: nOperator, eChangeOperator_IV,driftin,driftout,&
                              rbsumglobal,dt
    implicit none
    
    real, intent(in) :: Time
    integer, parameter :: nLogVars = 8
    logical, save :: IsFirstCall=.true.
    integer       :: iSpecies,iOperator
    !--------------------------------------------------------------------------

    ! Open file and write header if no restart on first call
    If (IsFirstCall .and. .not.IsRestart) then
       open(unit=UnitTmp_,file='IM/plots/CRCM.log',&
                  status='unknown')
       write(UnitTmp_,*) 'CRCM Logfile'
       write(UnitTmp_,*) NamePlotVarLog
       IsFirstCall = .false.
    else
       open(unit=UnitTmp_,file='IM/plots/CRCM.log',&
                  status='old', position='append')      
       IsFirstCall = .false.
    endif

    ! Write out iteration number (just use time in seconds)
    write(UnitTmp_,'(i7)',ADVANCE='NO') nint(Time/dt)

    ! write time 
    write(UnitTmp_,'(1es13.5)',ADVANCE='NO') time
    
  ! write out the operator changes
    do iSpecies=1,nSpecies
       if (iSpecies < nSpecies) then
          write(UnitTmp_,'(8es13.5)',ADVANCE='NO') & 
               rbsumglobal(iSpecies),eChangeOperator_IV(iSpecies,1:nOperator), &
               driftin(iSpecies),driftout(iSpecies)
       else
          write(UnitTmp_,'(8es13.5)') & 
               rbsumglobal(iSpecies),eChangeOperator_IV(iSpecies,1:nOperator),&
               driftin(iSpecies),driftout(iSpecies)
       endif
    enddo

    close(UnitTmp_)

  end subroutine crcm_plot_log
end module ModCrcmPlot

