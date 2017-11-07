!                                                                      C
!**********************************************************************C
!                                                                      C
!  Neutral atmosphere using JGITM 3-D output
!                                                                      C
!**********************************************************************C
!
Module ModJupiterAtmos
  implicit none

  private

  public :: JupiterAtmos

  integer :: nLon,nLat,nAtmos              ! number of grid points in J-GITM array
  real, allocatable :: lon_I(:),lat_I(:),alt_I(:) ! Coordinates for J-GITM array
  real, allocatable :: AtmosArray_VIII(:,:,:,:)   ! 3D J-GITM array
  
contains

  SUBROUTINE JupiterAtmos (mLat,mLon,nH2,nH,nH2O,nCH4,Temp)
    use ModParameters
    use ModInterpolate,  ONLY : trilinear
    use ModCommonPlanet, ONLY : nNeutral
    use ModCommonVariables, ONLY: GRAVTY,ALTD,ALTMIN,ALTMAX,NDIM
    implicit none
    
    real,               intent(in) :: mLat,mLon
    real, intent(out) :: nH2(0:MaxGrid+1),nH(0:MaxGrid+1),nH2O(0:MaxGrid+1)
    real, intent(out) :: nCH4(0:MaxGrid+1),Temp(0:MaxGrid+1)

    real :: Alt(0:NDIM+1) ! cm
    real :: nHe(0:NDIM+1)
    integer :: iAlt,maxJGITMAlt
    real :: Scaleheight_H2,Scaleheight_H,Scaleheight_H2O,Scaleheight_CH4

    real    :: param
    real, allocatable :: param3d_III(:,:,:)
    real    :: Coords(3)

    write(*,*) 'starting get_jupiter_atmos'

    Coords(2) = mLat
    Coords(3) = mLon
    
    Alt(0) = ALTMIN
    Alt(NDIM+1) = ALTMAX
    Alt(1:NDIM) = ALTD(1:nDim)
    write(*,*) AltMin,AltMax, AltD(1)
  ! Calculate trilinear interpolation of param3d_III at position Coords

    CALL read_JGITM_3D

    if(.not.allocated(param3d_III)) allocate(param3d_III(nAtmos,nLat,nLon))

    do iAlt=0,NDIM+1
       Coords(1) = Alt(iAlt)
       if (Coords(1).LE.alt_I(nAtmos)) then

          param3d_III = AtmosArray_VIII(2,:,:,:)
          param = trilinear(param3d_III, &
               1,nAtmos,1,nLat,1,nLon, Coords, alt_I(:),lat_I(:),lon_I(:),.false.)
          nH2(iAlt) = param/1.0e6
          
          param3d_III = AtmosArray_VIII(3,:,:,:)
          param = trilinear(param3d_III, &
               1,nAtmos,1,nLat,1,nLon, Coords, alt_I(:),lat_I(:),lon_I(:),.false.)
          nHe(iAlt) = param/1.0e6
          
          param3d_III = AtmosArray_VIII(4,:,:,:)
          param = trilinear(param3d_III, &
               1,nAtmos,1,nLat,1,nLon, Coords, alt_I(:),lat_I(:),lon_I(:),.false.)
          nH(iAlt) = param/1.0e6
          
          param3d_III = AtmosArray_VIII(5,:,:,:)
          param = trilinear(param3d_III, &
               1,nAtmos,1,nLat,1,nLon, Coords, alt_I(:),lat_I(:),lon_I(:),.false.)
          nCH4(iAlt) = param/1.0e6
          
          param3d_III = AtmosArray_VIII(6,:,:,:)
          param = trilinear(param3d_III, &
               1,nAtmos,1,nLat,1,nLon, Coords, alt_I(:),lat_I(:),lon_I(:),.false.)
          Temp(iAlt) = param

          nH2O(iAlt) = 0.0
          
          maxJGITMAlt = iAlt

     else
        Temp(iAlt) = Temp(maxJGITMAlt)   ! constant temp
        
        ! Scale height in cm
        if (iAlt==nDim+1)then
           Scaleheight_H2 =abs(1.380658e-19*Temp(iAlt)/(3.3452462e-27*GRAVTY(iAlt-1)))
           Scaleheight_H  =abs(1.380658e-19*Temp(iAlt)/(1.6726231e-27*GRAVTY(iAlt-1)))
           Scaleheight_H2O=Scaleheight_H2
           Scaleheight_CH4=abs(1.380658e-19*Temp(iAlt)/(2.65686432e-26*GRAVTY(iAlt-1)))
           !        Scaleheight_H2O=1.380658e-19*Temp(iAlt)/(2.99e-26*GRAVTY(iAlt))
        else
           Scaleheight_H2 =abs(1.380658e-19*Temp(iAlt)/(3.3452462e-27*GRAVTY(iAlt)))
           Scaleheight_H  =abs(1.380658e-19*Temp(iAlt)/(1.6726231e-27*GRAVTY(iAlt)))
           Scaleheight_H2O=Scaleheight_H2
           Scaleheight_CH4=abs(1.380658e-19*Temp(iAlt)/(2.65686432e-26*GRAVTY(iAlt)))
           !        Scaleheight_H2O=1.380658e-19*Temp(iAlt)/(2.99e-26*GRAVTY(iAlt))
        endif
        nH2(iAlt)  = nH2(iAlt-1) *exp(-(Alt(iAlt)-Alt(iAlt-1))/Scaleheight_H2)
        nH(iAlt)   = nH(iAlt-1)  *exp(-(Alt(iAlt)-Alt(iAlt-1))/Scaleheight_H)
        nH2O(iAlt) = nH2O(iAlt-1)*exp(-(Alt(iAlt)-Alt(iAlt-1))/Scaleheight_H2O)
        nCH4(iAlt) = nCH4(iAlt-1)*exp(-(Alt(iAlt)-Alt(iAlt-1))/Scaleheight_CH4)

     endif
  end do

  !call con_stop('')
    
    write(*,*) 'done get_jupiter_atmos'

  END SUBROUTINE JupiterAtmos

  subroutine read_JGITM_3D
    use ModCommonPlanet, ONLY : nNeutral
    use ModIoUnit,     ONLY: UnitTmp_
    
    real, allocatable :: ReadArray_VIII(:,:,:,:)   ! 3D J-GITM array
    character (len=189) :: tmpline
    character (len=50)  :: filenames(nNeutral+2)
    integer :: iNeutralVar,line
    
    filenames(1) = 'PW/3D_GITM_AuroralMask.txt'
    filenames(2) = 'PW/3D_GITM_H2Density.txt'
    filenames(3) = 'PW/3D_GITM_HeDensity.txt'
    filenames(4) = 'PW/3D_GITM_HDensity.txt'
    filenames(5) = 'PW/3D_GITM_CH4_Density.txt'
    filenames(6) = 'PW/3D_JGITM_Temperature.txt'

    do iNeutralVar = 1,nNeutral+2
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
            allocate(AtmosArray_VIII(nNeutral+2,nAtmos,nLat,nLon))
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

  end Module ModJupiterAtmos
