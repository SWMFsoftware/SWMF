Module Mod3DAtmos
  implicit none

  private

  public :: get_jupiter_atmos

  integer, parameter :: nNeutralVar = 6 ! number of variables to read in
  
  integer :: nLon,nLat,nAtmos              ! number of grid points in J-GITM array
  real, allocatable :: lon_I(:),lat_I(:),alt_I(:) ! Coordinates for J-GITM array
  real, allocatable :: AtmosArray_VIII(:,:,:,:)   ! 3D J-GITM array
  real :: ReadArray_II(2,33) ! Atreya array
  
contains

  subroutine get_jupiter_atmos(mLat,mLon,nSpecies,NeutralDens_IC,NeutralTemp_C)

    use ModInterpolate,  ONLY : trilinear,linear
    use ModSeGrid,       ONLY : Alt_C, nAlt, IsVerbose
!    use ModSeBackground, ONLY : mLat,mLon
    
    integer,            intent(in) :: nSpecies     ! 4
    real,               intent(in) :: mLat,mLon
    real,               intent(out) :: NeutralDens_IC(nSpecies,nAlt)
    real,               intent(out) :: NeutralTemp_C(nAlt)
    
    integer :: iIono, iNeutralVar
    
    real    :: param,H2A
    real, allocatable :: param3d_III(:,:,:)
    real    :: Coords(3)
    logical :: DoExtrapolate

    write(*,*) 'starting get_jupiter_atmos'

    Coords(2) = mLat
    !    Coords(3) = mLon
    Coords(3) = 20
    
    ! Calculate trilinear interpolation of param3d_III at position Coords

    CALL read_Atreya
    CALL read_JGITM_3D

    if(IsVerbose) write(*,*) 'Neutral Array:',ReadArray_II
    
    if(.not.allocated(param3d_III)) allocate(param3d_III(nAtmos,nLat,nLon))

    do iNeutralVar = 1,nNeutralVar
       param3d_III = AtmosArray_VIII(iNeutralVar,:,:,:)
       do iIono=1,nAlt
          Coords(1) = Alt_C(iIono)
          if (Coords(1) > alt_I(nAtmos)) then
             param = 0.0
          else
             param = trilinear(param3d_III, &
                  1,nAtmos,1,nLat,1,nLon, Coords, alt_I(:),lat_I(:),lon_I(:),.false.)
          endif
          if (iNeutralVar.eq.2) then ! for H2
             H2A = 0.0
             H2A = linear(ReadArray_II(1,:),1,33, &
                  Coords(1),ReadArray_II(2,:)*1.e5,.false.)
             write(*,*) H2A*1.e6, param, Coords(1)
! *** Use Atreya neutrals ***
!             if (H2A*1.e6.GT.param.AND. &
!                  Coords(1).LT.maxval(ReadArray_II(2,:)*1.e5)) then
!                write(*,*) 'overwriting param'
!                param = H2A*1.e6
!             endif
          endif
          ! are iCell_D, Dist_D important?
          if (iNeutralVar > 1 .and. iNeutralVar < nNeutralVar) then
             ! store density variables, convert to cm^-3
             NeutralDens_IC(iNeutralVar-1,iIono) = param/1.0e6
          else if (iNeutralVar == nNeutralVar) then
             ! store temperature (final) variable along field line
             NeutralTemp_C(iIono) = param             
          end if
       end do
    end do

!    call con_stop('')
    
    write(*,*) 'done get_jupiter_atmos'
  end subroutine get_jupiter_atmos

  subroutine read_JGITM_3D
    use ModIoUnit,     ONLY: UnitTmp_
    
    real, allocatable :: ReadArray_VIII(:,:,:,:)   ! 3D J-GITM array
    character (len=189) :: tmpline
    character (len=50)  :: filenames(nNeutralVar)
    integer :: iNeutralVar,line
    
    filenames(1) = 'PW/3D_GITM_AuroralMask.txt'
    filenames(2) = 'PW/3D_GITM_H2Density.txt'
    filenames(3) = 'PW/3D_GITM_HeDensity.txt'
    filenames(4) = 'PW/3D_GITM_HDensity.txt'
    filenames(5) = 'PW/3D_GITM_CH4_Density.txt'
    filenames(6) = 'PW/3D_JGITM_Temperature.txt'

    do iNeutralVar = 1,nNeutralVar
       write(*,*) iNeutralVar
       open(UnitTmp_,FILE=filenames(iNeutralVar),STATUS='OLD')
       do line = 1,6
          read(UnitTmp_,'(a)') tmpline
       end do

       read(UnitTmp_,'(a)') tmpline
       read (tmpline(len_trim(tmpline)-2:),*) nLat
       read(UnitTmp_,'(a)') tmpline
       read (tmpline(len_trim(tmpline)-2:),*) nLon
       read(UnitTmp_,'(a)') tmpline
       read (tmpline(len_trim(tmpline)-2:),*) nAtmos
       
       do line = 1,8
          read(UnitTmp_,'(a)') tmpline
       end do

       if(.not.allocated(ReadArray_VIII)) &
            allocate(ReadArray_VIII(4,nAtmos,nLat,nLon))
       if(.not.allocated(AtmosArray_VIII)) &
            allocate(AtmosArray_VIII(nNeutralVar,nAtmos,nLat,nLon))
       if(.not.allocated(lon_I)) allocate(lon_I(nLon))
       if(.not.allocated(lat_I)) allocate(lat_I(nLat))
       if(.not.allocated(alt_I)) allocate(alt_I(nAtmos))

       read(UnitTmp_,*) ReadArray_VIII(:,:,:,:)

       close(UnitTmp_)

       AtmosArray_VIII(iNeutralVar,:,:,:) = ReadArray_VIII(4,:,:,:)
       lon_I(:) = ReadArray_VIII(3,1,1,:)
       lat_I(:) = ReadArray_VIII(2,1,:,1)
       alt_I(:) = ReadArray_VIII(1,:,1,1)*1.0e5 ! convert alt to cm

    end do

  end subroutine read_JGITM_3D

  subroutine read_Atreya
    use ModIoUnit,     ONLY: UnitTmp_
    
    character (len=189) :: tmpline ! 22-ish?
    character (len=50)  :: filename='PW/nH2_Atreya1981.dat'
    integer :: iNeutralVar,line

    open(UnitTmp_,FILE=filename,STATUS='OLD')

    read(UnitTmp_,'(a)') tmpline
    
    read(UnitTmp_,*) ReadArray_II(:,:)    

  end subroutine read_Atreya
    
end Module Mod3DAtmos
