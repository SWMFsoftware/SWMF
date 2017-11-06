Module ModImSat
  implicit none
  
  logical :: DoWriteSats  = .false., ReadRestartSat = .false.
  logical,            allocatable :: IsFirstWrite(:,:)
  integer :: nImSats = 0, iStartTime = 0
  real    :: DtSatOut = 1.0 !Default write sat output every minute
  real,               allocatable :: SatLoc_3I(:,:,:), SatFlux_3I(:,:,:)
  real,               allocatable :: BfieldEq2_G(:,:)
  real,               allocatable :: LonGrid_G(:), Flux_G(:,:,:,:,:)
  character(len=100), allocatable :: NameSat_I(:)

  ! Variables for the Prerunsat output.
  real :: DtReadSat=60.
  logical :: DoWritePrerunSat = .false., UsePrerunSat = .false.

contains
    !===========================================================================
  subroutine write_im_sat(iSatIn, nLat,nLon,nEnergy,nAngle,Flux_C)
    ! Write solution interpolated to satellite position to satellite files.
    use ModIoUnit,      ONLY: io_unit_new, UnitTmp_
    use ModInterpolate, ONLY: bilinear, trilinear
    use ModNumConst,    ONLY: cDegToRad,cRadToDeg
    use ModCimi,        ONLY: t=>time
    use ModCimiPlanet,  ONLY: nSpecies=>nspec
    use ModCimiGrid,    ONLY: LonGrid_I=>phi, LatGrid_I=>xlat, &
                              AngleGrid_I=>sinAo
    use ModFieldTrace,  ONLY: BfieldEq_C => bo, iba
    use ModImTime,      ONLY: iCurrentTime_I

    implicit none
    
    ! Arguments
    integer, intent(in) :: iSatIn
    integer,intent(in) :: nLat,nLon,nEnergy,nAngle
    real,   intent(in) :: Flux_C(nSpecies,nLat,nLon,nEnergy,nAngle)
    
    ! Name of subroutine
    character(len=*), parameter :: NameSubSub = 'write_im_sat'

    !internal variables
    integer             :: iError
    logical             :: IsExist
    character(len=100)  :: NameSatFile, StringTime
    character(len=3000) :: HeadVar
    real                :: SatLat, SatLon, SatAng, &
                          LatSatGen,LonSatGen, AngSatGen
    real                :: SatB2, RatioBeqBsat
    integer             :: iSatLat, iSatLonMin,iSatLonMax, &
         iSatAng, iAngle, iEnergy
    integer             :: iSpecies
    character(len=2)    :: NameSpecies
    character(len=8)    :: NameChannel
    character(len=3)    :: numChannels
    !-------------------------------------------------------------------------
    ! Allocate array for satellite flux

    if (.not. allocated(IsFirstWrite)) then
       allocate( IsFirstWrite(3,nImSats) )
       IsFirstWrite = .true.
    end if
    
    if(.not. allocated(SatFlux_3I))  &
         allocate(SatFlux_3I(nSpecies,nEnergy,nAngle))
    if(.not. allocated(BfieldEq2_G)) allocate(BfieldEq2_G(nLat,0:nLon+1))

    ! Allocate Lat, Lon and Flux grids with Ghost Cells
    if(.not. allocated(LonGrid_G)) then
       allocate(LonGrid_G(0:nLon+1))
       LonGrid_G(1:nLon)=LonGrid_I(1:nLon)*cRadToDeg
       LonGrid_G(nLon+1)=360.0
       LonGrid_G(0)=LonGrid_I(nLon)-360.0
    endif
    
    if(.not. allocated(Flux_G)) &
         allocate(Flux_G(nSpecies,nLat,0:nLon+1,nEnergy,nAngle))
    if(iSatIn == 1) then
       Flux_G(nSpecies,:,:,:,:)=0.0
       Flux_G(1:nSpecies,1:nLat,1:nLon,1:nEnergy,1:nAngle) = &
            Flux_C(1:nSpecies,1:nLat,1:nLon,1:nEnergy,1:nAngle)
       Flux_G(1:nSpecies,1:nLat,nLon+1,1:nEnergy,1:nAngle) = &
            Flux_G(1:nSpecies,1:nLat,1,1:nEnergy,1:nAngle)
       Flux_G(1:nSpecies,1:nLat,0,1:nEnergy,1:nAngle) = &
            Flux_G(1:nSpecies,1:nLat,nLon,1:nEnergy,1:nAngle)
    endif

    ! Set BfieldEq^2 on first sat call
    if (iSatIn == 1) then 
       BfieldEq2_G(:,1:nLon) = (BfieldEq_C(:,1:nLon))**2.0
       BfieldEq2_G(:,nLon+1) =  BfieldEq2_G(:,1)
       BfieldEq2_G(:,0)      =  BfieldEq2_G(:,nLon)
    endif

    ! Collect variables.
    SPECIES_LOOP: do iSpecies=1,nSpecies
       
       ! Get satellite location in generalized coordinates.
       ! and interpolate flux to satellite location
       SatLat=SatLoc_3I(1,2,iSatIn)
       SatLon=mod(SatLoc_3I(2,2,iSatIn),360.0)
       if (SatLoc_3I(3,2,iSatIn) == 3 .and. SatLat <= maxval(LatGrid_I)) then
          !get generallized sat lat
          call locate1IM(LatGrid_I, nLat, SatLat, iSatLat)
          LatSatGen = iSatLat &
               + (LatGrid_I(iSatLat+1) - SatLat) &
               / (LatGrid_I(iSatLat+1)-LatGrid_I(iSatLat))
          
          !get generallized sat lon
          call locate1IM(LonGrid_G, nLon+2, SatLon, iSatLonMax)
          LonSatGen=SatLon/(360.0/nLon)+1
          iSatLonMin=floor(LonSatGen)
          select case (iSatLonMax)
          case (48)
             iSatLonMax=iSatLonMax
          case DEFAULT
             iSatLonMax=mod(ceiling(LonSatGen),nLon)
          end select
          
          !get B^2 at sat
          SatB2 = SatLoc_3I(4,2,iSatIn)
          RatioBeqBsat  = sqrt(bilinear(BfieldEq2_G,1,nLat,0,nLon+1,& 
               (/ LatSatGen, LonSatGen /) ) / SatB2)
          
          ! Beq must be minimum so ratio can not be larger than 1
          if (RatioBeqBsat > 1.0) RatioBeqBsat=1.0
          
          !check that sat is in iba domain
          if (iSatLat >= iba(iSatLonMin) .or. &
               iSatLat >= iba(iSatLonMax)) then
             !treat as open
             LatSatGen=-1.0
             LonSatGen=-1.0
             SatFlux_3I(iSpecies,:,:) = 0.0
             SatLoc_3I(3,2,iSatIn) = -1
          else
             !Interpolate solution to sat location
             do iAngle=1,nAngle
                !Get generallized sat pitch-angle
                SatAng = AngleGrid_I(iAngle)*RatioBeqBsat
                call locate1IM(AngleGrid_I, nAngle, SatAng, iSatAng)
                ! Angle grid does not go all the way to zero, if angle
                ! falls below min angle in grid then set it to min angle
                if (iSatAng < 1) then 
                   iSatAng   = 1
                   AngSatGen = 1.0
                else
                   AngSatGen =  iSatAng &
                        + (AngleGrid_I(iSatAng+1) - SatAng) &
                        / (AngleGrid_I(iSatAng+1)-AngleGrid_I(iSatAng))
                endif
                do iEnergy=1,nEnergy
                   SatFlux_3I(iSpecies,iEnergy,iAngle) = &
                        trilinear( Flux_G(iSpecies,:,:,iEnergy,:), &
                        1,nLat,0,nLon+1,1,nAngle, &
                        (/ LatSatGen, LonSatGen, AngSatGen/) )
                enddo
             enddo
          endif
       else
          ! When satellite is on open field lines or outside domain 
          ! set flux to zero
          LatSatGen=-1.0
          LonSatGen=-1.0
          SatFlux_3I(iSpecies,:,:) = 0.0
       end if
       
       ! Write satellite output
       
       ! Build file name; 
       if (IsFirstWrite(iSpecies,iSatIn)) iStartTime = int(t)
       
       !set the file name for each species
       SELECT CASE (iSpecies)
       CASE (1)
          if (iSpecies==nSpecies) then
             NameSpecies='_e'
          else
             NameSpecies='_h'
          endif
       CASE (2)
          if (iSpecies==nSpecies) then
             NameSpecies='_e'
          else
             NameSpecies='_o'
          endif
       CASE DEFAULT
          NameSpecies='_e'
       END SELECT
       
       write(NameSatFile, '(a, i6.6, a)')                    &
            'IM/plots/'//'sat_'//trim(NameSat_I(iSatIn))// &
             NameSpecies//'flux_t',iStartTime,'.sat'
       
       ! Open file in appropriate mode.  Write header if necessary.
       inquire(file=NameSatFile, exist=IsExist)
       
       if ( (.not.IsExist) .or. IsFirstWrite(iSpecies,iSatIn) ) then
          IsFirstWrite(iSpecies,iSatIn) = .false.
          HeadVar = 'it year mo dy hr mn sc msc X Y Z'
          open(unit=UnitTmp_, file=trim(NameSatFile), &
               status='replace', iostat=iError)
          if(iError /= 0) call CON_stop &
               (NameSubSub//' Error opening file '//NameSatFile)
          do iAngle=1,nAngle
             do iEnergy=1,nEnergy
                write(NameChannel,"(a,i2.2,a,i2.2)") &
                     ' E',iEnergy,'@A',iAngle
                HeadVar = trim(HeadVar) // NameChannel
             enddo
          enddo
          ! Write header
          write(UnitTmp_, '(2a)')'IM results for SWMF trajectory file ', &
               trim(NameSat_I(iSatIn))//NameSpecies
          write(UnitTmp_,"(A)") trim(HeadVar)
       else
          open(unit=UnitTmp_, file=trim(NameSatFile), status='OLD',&
               position='append', iostat=iError)
          if(iError /= 0) &
               call CON_stop(NameSubSub//' Error opening file '//NameSatFile)
       end if

       write(UnitTmp_,'(i7)',ADVANCE='NO') int(t)
       write(UnitTmp_,'(i5,5(1X,i2.2),1X,i3.3)',ADVANCE='NO') &
            iCurrentTime_I
       write(UnitTmp_,'(3es13.5)',ADVANCE='NO') SatLoc_3I(1:3,1,iSatIn)
       write(numChannels,'(I3)') nEnergy*nAngle
       write(UnitTmp_,'('//numChannels//'es13.5)') &
            SatFlux_3I(iSpecies,:,:)!SatVar_I
       
       close(UnitTmp_)
       
    end do SPECIES_LOOP

    ! Deallocate Flux_G after last satellite to save memory
    if(iSatIn == nImSats) deallocate(Flux_G)

    if (iSatIn == nImSats .and. DoWritePrerunSat) then
       call save_prerun_sat(t)
    end if
        
  end subroutine write_im_sat

  subroutine save_prerun_sat(tSimulation)

    use ModGmCimi
    use ModIoUnit,  ONLY: UnitTmp_
    real, intent(in) :: tSimulation
    integer          :: iTimeOut,iSat,iRow
    Character(len=100) :: NameFile             ! input file name

    iTimeOut=int(tSimulation)
    
    write(NameFile,"(a,i8.8,a)") &
         'IM/PrerunSat_',iTimeOut,'.dat'
    open(UnitTmp_,file=NameFile,status="replace", form="unformatted")

    write(UnitTmp_) nImSats

    do iSat=1,nImSats
       
       write(UnitTmp_) NameSat_I(iSat)

    end do

    do iSat=1,nImSats

       do iRow=1,2

          write(UnitTmp_) SatLoc_3I(1:4,iRow,iSat)

       end do
       
    end do
    
    close(UnitTmp_)    
    
  end subroutine save_prerun_sat

  subroutine read_prerun_sat(tSimulation)

    use ModGmCimi
    use ModIoUnit,  ONLY: UnitTmp_
    real, intent(in) :: tSimulation
    integer          :: iTimeOut,iSat,iRow
    integer, save    :: iTimeOutPrev = -1
    Character(len=100) :: NameFile             ! input file name

    iTimeOut=int(floor(tSimulation/DtReadSat) * DtReadSat)

    if(iTimeOut == iTimeOutPrev) then
       return
    else
       iTimeOutPrev =iTimeOut
    end if

    write(NameFile,"(a,i8.8,a)") &
         'IM/PrerunSat_',iTimeOut,'.dat'   

    open(UnitTmp_,file=NameFile,status="old", form="unformatted")

    read(UnitTmp_) nImSats

    if (.not. allocated(NameSat_I)) allocate(NameSat_I(nImSats))
    if (.not. allocated(SatLoc_3I)) allocate(SatLoc_3I(4,2,nImSats))
    
    do iSat=1,nImSats
       
       read(UnitTmp_) NameSat_I(iSat)

    end do

    do iSat=1,nImSats

       do iRow=1,2

          read(UnitTmp_) SatLoc_3I(1:4,iRow,iSat)

       end do
       
    end do

    close(UnitTmp_)
    
  end subroutine read_prerun_sat
  
end Module ModImSat

