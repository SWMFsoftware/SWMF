Module ModPrerunField
  implicit none
  real    :: DtRead
  logical :: DoWritePrerun = .false., UsePrerun=.false.
contains
  !=============================================================================
  subroutine save_prerun(tSimulation)
    use ModGmRb
    use rbe_cread2, ONLY: xnswa, vswa,bxw,byw,bzw 
    use ModIoUnit,  ONLY: UnitTmp_
    real, intent(in) :: tSimulation
    Character(len=100) :: NameFile             !output file name
    integer :: iTimeOut, iLat, iLon, iVar
    !---------------------------------------------------------------------------
    
    ! Create Filename and open file
    iTimeOut=int(tSimulation)
    write(NameFile,"(a,i8.8,a)") &
         'RB/PrerunField_',iTimeOut,'.dat'   
    open(UnitTmp_,file=NameFile,status="replace", form="unformatted")

    ! Write out SW values
    write(UnitTmp_) xnswa(1), vswa(1), bxw(1), byw(1), bzw(1)
    
    ! Write out nPoint and nIntegral
    write(UnitTmp_) nPoint, nIntegral

    ! Write out StateIntegral_IIV
    do iLon = 1,nLon
       do iLat = 1,nLat
          write(UnitTmp_) StateIntegral_IIV(iLat,iLon,1:nIntegral)
       end do
    end do
    
    ! Write out StateLine_VI
    do iVar = 1,nVar
       write(UnitTmp_) StateLine_VI(iVar, 1:nPoint)
    end do
    
    close(UnitTmp_)

  end subroutine save_prerun
    
  !=============================================================================
  subroutine read_prerun(tSimulation)
    use ModGmRb
    use rbe_cread2, ONLY: xnswa, vswa, bxw, byw, bzw, UseGm
    use ModIoUnit,  ONLY: UnitTmp_
    real, intent(in) :: tSimulation
    integer          :: iTimeOut
    integer,save     :: iTimeOutPrev = -1
    integer          :: n, iLat, iLon, iVar
    Logical, save    :: IsFirstCall =.true.
    Character(len=100) :: NameFile             ! input file name
    !---------------------------------------------------------------------------

    ! Set filename for reading
    iTimeOut=int(floor(tSimulation/DtRead) * DtRead)
    
    if(iTimeOut == iTimeOutPrev) then
       return
    else
       iTimeOutPrev =iTimeOut
    end if 
    
    write(NameFile,"(a,i8.8,a)") &
         'RB/PrerunField_',iTimeOut,'.dat'   
    open(UnitTmp_,file=NameFile,status="old", form="unformatted")

    ! read SW values
    read(UnitTmp_) xnswa(1), vswa(1), bxw(1), byw(1), bzw(1)
    
    !  read nPoint and nIntegral
    read(UnitTmp_) nPoint, nIntegral

    ! Allocate StateLine and StateIntegral
    if (allocated(StateLine_VI)) then
       deallocate(StateLine_VI,StateIntegral_IIV)
    endif
    if (.not.allocated(StateLine_VI)) then
       allocate(StateLine_VI(nVar,nPoint),&
            StateIntegral_IIV(nLat,nLon,nIntegral))
    endif

    ! read StateIntegral_IIV
    do iLon = 1,nLon
       do iLat = 1,nLat
          read(UnitTmp_) StateIntegral_IIV(iLat,iLon,1:nIntegral)
       end do
    end do
    
    ! read StateLine_VI
    do iVar = 1,nVar
       read(UnitTmp_) StateLine_VI(iVar, 1:nPoint)
    end do
    
    close(UnitTmp_)
    
    ! create an index array on the first call
    if (IsFirstCall) then
       n = 0
       do iLon = 1, nLon
          do iLat = 1, nLat
           n = n+1
           iLineIndex_II(iLon,iLat) = n
        end do
     end do
     IsFirstCall = .false.
     UseGm = .true.
  endif

    
  end subroutine read_prerun

end Module ModPrerunField
