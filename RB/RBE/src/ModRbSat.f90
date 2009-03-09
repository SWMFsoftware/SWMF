Module ModRbSat
  implicit none
  
  logical :: DoWriteSats  = .false.  !!!DTW 2007                                                                                      
  logical :: IsFirstWrite = .true.
  integer :: nRbSats = 0, iStartIter = 0
  real,               allocatable :: SatLoc_3I(:,:,:), SatFlux_II(:,:)
  real,               allocatable :: BfieldEq2_C(:,:),SatVar_I(:)
  character(len=100), allocatable :: NameSat_I(:)
  
contains
    !===========================================================================
  subroutine write_rb_sat(iSatIn, nLat,nLon,nEnergy,nAngle,Flux_C)
    ! Write solution interpolated to satellite position to satellite files.
    use ModIoUnit,      ONLY: io_unit_new, UnitTmp_
    use ModInterpolate, ONLY: bilinear, trilinear
    use ModNumConst,    ONLY: cDegToRad
    use rbe_cread2,     ONLY: tint
    use rbe_time,       ONLY: t
    use rbe_cgrid,      ONLY: LonGrid_I=>xmltd, LatGrid_I=>xlati, &
                              AngleGrid_I=>gridy
    use rbe_cfield,     ONLY: BfieldEq_C => bo
    use ModRbTime,      ONLY: iCurrentTime_I
    use ModGmRb,        ONLY: StateIntegral_IIV
    implicit none
    
    ! Arguments
    integer, intent(in) :: iSatIn
    
    ! Name of subroutine
    character(len=*), parameter :: NameSubSub = 'write_rb_sat'

    integer,intent(in) :: nLat,nLon,nEnergy,nAngle
    real,   intent(in) :: Flux_C(nLat,nLon,nEnergy,nAngle)
    !internal variables
    integer            :: iError
    logical            :: IsExist
    character(len=100) :: NameSatFile, StringTime
    character(len=200) :: HeadVar
    real               :: SatLat, SatLon, SatAng, &
                          LatSatGen,LonSatGen, AngSatGen
    real               :: SatB2, RatioBeqBsat
    integer            :: iSatLat, iSatLon, iSatAng, iAngle, iEnergy
    character(len=8)  :: NameChannel
    !-------------------------------------------------------------------------
    
    !Allocate array for satellite flux
    
    if(.not. allocated(SatFlux_II))  allocate(SatFlux_II(nEnergy,nAngle))
    if(.not. allocated(BfieldEq2_C)) allocate(BfieldEq2_C(nLat,nLon))
    if(.not. allocated(SatVar_I))    allocate(SatVar_I(2*nEnergy))
    
    !Set BfieldEq^2 on first sat call
    if (iSatIn == 1) BfieldEq2_C = BfieldEq_C**2.0

    ! Build file name; 
    if (IsFirstWrite) iStartIter = int(t/tint)

    write(NameSatFile, '(a, i6.6, a)')                    &
         'RB/plots/'//'sat_'//trim(NameSat_I(iSatIn))// &
         '_n',iStartIter, '.sat'

    ! Open file in appropriate mode.  Write header if necessary.
    inquire(file=NameSatFile, exist=IsExist)

    if ( (.not.IsExist) .or. IsFirstWrite ) then
       IsFirstWrite = .false.
       HeadVar = 'it year mo dy hr mn sc msc X Y Z'
       open(unit=UnitTmp_, file=trim(NameSatFile), status='replace', iostat=iError)
       if(iError /= 0) call CON_stop &
            (NameSubSub//' Error opening file '//NameSatFile)
       do iEnergy=1,nEnergy
          write(NameChannel,"(a,i2.2,a,i2.2)") &
               ' E',iEnergy,' A',iEnergy
          HeadVar = trim(HeadVar) // NameChannel
       enddo
       
       ! Write header
       write(UnitTmp_, '(2a)')'RB results for SWMF trajectory file ', &
            trim(NameSat_I(iSatIn))
       write(UnitTmp_,*) trim(HeadVar)
    else
       open(unit=UnitTmp_, file=trim(NameSatFile), status='OLD',&
            position='append', iostat=iError)
       if(iError /= 0) &
            call CON_stop(NameSubSub//' Error opening file '//NameSatFile)
    end if

    ! Collect variables.


    ! Get satellite location in generalized coordinates.
    ! and interpolate flux to satellite location
    write(*,*) 'RB:Lat Lon:', SatLoc_3I(1:2,2,iSatIn)
    if (SatLoc_3I(3,2,iSatIn) == 3) then
       SatLat=SatLoc_3I(1,2,iSatIn)*cDegToRad
       SatLon=mod(SatLoc_3I(2,2,iSatIn)+180.0,360.0)

       !get generallized sat lat
       call locate1(LatGrid_I, nLat, SatLat, iSatLat)
       LatSatGen = iSatLat &
            + (LatGrid_I(iSatLat+1) - SatLat) &
            / (LatGrid_I(iSatLat+1)-LatGrid_I(iSatLat))

       !get generallized sat lon
       call locate1(LonGrid_I, nLon, SatLon, iSatLon)
       LonSatGen=SatLon/(360.0/nLon)+1
       iSatLon=floor(LonSatGen)
       LonSatGen = iSatLon &
            + (LonGrid_I(iSatLon+1) - SatLon) &
            / (LonGrid_I(iSatLon+1)-LonGrid_I(iSatLon))

       !get B^2 at sat
       SatB2 = SatLoc_3I(4,2,iSatIn)
       RatioBeqBsat  = sqrt(bilinear(BfieldEq2_C,1,nLat,1,nLon,& 
            (/ LatSatGen, LonSatGen /) ) / SatB2)
       
       ! Beq must be minimum so ratio can not be larger than 1
       if (RatioBeqBsat > 1.0) RatioBeqBsat=1.0

       !Interpolate solution to sat location
       do iAngle=1,nAngle
          !Get generallized sat pitch-angle
          SatAng = AngleGrid_I(iAngle)*RatioBeqBsat
          call locate1(AngleGrid_I, nAngle, SatAng, iSatAng)
          ! Angle grid does not go all the way to zero, if angle falls below min
          ! angle in grid then set it to min angle
          if (iSatAng < 1) then 
             iSatAng   = 1
             AngSatGen = 1.0
          else
             AngSatGen =  iSatAng &
            + (AngleGrid_I(iSatAng+1) - SatAng) &
            / (AngleGrid_I(iSatAng+1)-AngleGrid_I(iSatAng))
          endif
          do iEnergy=1,nEnergy
             SatFlux_II(iEnergy,iAngle) = &
                  trilinear( Flux_C(:,:,iEnergy,:),1,nLat,1,nLon, &
                  1,nAngle,(/ LatSatGen, LonSatGen, AngSatGen/) )
          enddo
       enddo
    else
       ! When satellite is on open field lines set flux to zero
       LatSatGen=-1.0
       LonSatGen=-1.0
       SatFlux_II(:,:) = 0.0
    endif
    
    ! Format output
    if (SatLoc_3I(3,2,iSatIn) == 3) then
       SatVar_I(:)=0.0
       do iEnergy=1,nEnergy
          SatVar_I(2*iEnergy-1)=sum(SatFlux_II(iEnergy,:))
          do iAngle=1,nAngle
             SatVar_I(2*iEnergy)=&
                  SatVar_I(2*iEnergy) &
                  + SatFlux_II(iEnergy,iAngle)*AngleGrid_I(iAngle)
          enddo
          SatVar_I(2*iEnergy)=&
               SatVar_I(2*iEnergy)/SatVar_I(2*iEnergy-1)
       enddo
    else
       SatVar_I(:)=0.0
    endif
    ! Write satellite output
    
    write(UnitTmp_,'(i7)',ADVANCE='NO') int(t/tint)
    write(UnitTmp_,'(i5,5(1X,i2.2),1X,i3.3)',ADVANCE='NO') &
         iCurrentTime_I
    write(UnitTmp_,'(3es13.5)',ADVANCE='NO') SatLoc_3I(1:3,1,iSatIn)
    write(UnitTmp_,'(100es13.5)') SatVar_I
    
    close(UnitTmp_)

  end subroutine write_rb_sat

end Module ModRbSat

