program test_w05sc_mod
!
! To compile with pgf90 (e.g., hao linux) (note 8-byte reals):
! pgf90 -o w05sc -r8 read_data.F90 w05sc.F90 test_w05sc.F90
!
!
  use w05sc,only: setmodel,epotval,mpfac
  implicit none
  real :: by,bz,tilt,swvel,swden
  ! points inside boundary is recognized by fill
!  real,parameter :: fill=1.e36
  real,parameter :: fill  =1.e-6
! integer,parameter :: nlat=4, nlon=3
! integer,parameter :: nlat=97, nlon=25
!  integer,parameter :: nlat=97, nlon=1+24*4
!  integer,parameter :: nlat=40, nlon=1+24*4
!  integer,parameter :: nlat=101, nlon=1+24*4

! here nlat = nf + 1 and nlon = nlt + 1

!  integer,parameter :: nlat=125, nlon=1+24*4
  integer,parameter :: nlat=121, nlon=1+24*4

  real :: mlat(nlat),mlon(nlon)
  real :: lon(nlon),mlt(nlon)

  real :: lat_sami(nlon,nlat)
  real :: lon_sami(nlon,nlat)
  real :: pot_sami(nlon,nlat)

  real :: phi_weimer(nlat,nlon)
  
  real :: epot(nlon,nlat),fac(nlon,nlat)
  integer :: i,j,lu
  integer :: appid,iwk_ps,ier
  character(len=24) :: datfile
  character(len=*),parameter ::  file_path="./"
  character(len=80)  :: line

  integer :: ntime,year,doy,itime
  real    :: hrutw

  ! file_path="/Users/danielweimer/Work_Files/FORTRAN/W05/"

  open(101,file='phi_weimer.inp',form='unformatted')

  open(121,file='weimer_grid.dat',form='formatted')
  read(121,*) mlat,mlon
  close(121)

  open(122,file='weimer_input.inp',form='formatted')
  read(122,*) line
  read(122,*) line
  read(122,*) line
  read(122,*) line
  read(122,*) ntime,year,doy
  print *,'ntime,year,doy',ntime,year,doy

!  open(123,file='weimer_time.dat',form='formatted')

  do itime = 1,ntime
!  do itime = 1,22

    print *,'itime',itime

  read(122,*) hrutw,by,bz,swvel,swden,tilt

!     print *,hrutw,by,bz,swvel,swden,tilt
!  write(123,*) hrutw

!
!  by = 0.
!  bz = -10.
!  tilt = 0.
!  swvel = 450.
!  swden = 4.
  call setmodel(by,bz,tilt,swvel,swden,file_path,'epot')
!

! mlat = (/-90., -88.1238292398491, -86.2386359278657, -84.3344382773342,&
!    -82.4013318763435, -80.4295344892688, -78.4094552099168,&
!    -76.331796630125, -74.1876988925388, -71.9689341802758,&
!    -69.6681589022773, -67.2792279882741, -64.7975706790533,&
!    -62.2206194320588, -59.5482728298363, -56.7833601290164,&
!    -53.9320608459732, -51.0042204168578, -48.0134966005524,&
!    -44.9772754602266, -41.916313892128, -38.8540980954293,&
!    -35.8159497801506, -32.8279553674349, -29.9158266703621,&
!    -27.1038148776609, -24.4137889090065, -21.8645574169981,&
!    -19.4714697638694, -17.2462861630082, -15.1972697734841,&
!    -13.3294282264571, -11.6448185129562, -10.142824406667,&
!    -8.82031765103987, -7.67162666281269, -6.68827297583048,&
!    -5.85851734698832, -5.16689314460211, -4.5940469432968,&
!    -4.11722526306697, -3.71151170575937, -3.35148255039153,&
!    -3.01257883277328, -2.67136426606314, -2.3036287214954,&
!    -1.87754943767857, -1.32687203939232, -7.72840966450717e-08,&
!    1.32687203939232, 1.87754943767857, 2.3036287214954, 2.67136426606314,&
!    3.01257883277328, 3.35148255039153, 3.71151170575936, 4.11722526306697,&
!    4.59404694329679, 5.16689314460211, 5.85851734698832, 6.68827297583048,&
!    7.67162666281268, 8.82031765103987, 10.142824406667, 11.6448185129562,&
!    13.3294282264571, 15.1972697734841, 17.2462861630082, 19.4714697638694,&
!    21.8645574169981, 24.4137889090064, 27.1038148776609, 29.9158266703621,&
!    32.8279553674348, 35.8159497801506, 38.8540980954293, 41.916313892128,&
!    44.9772754602266, 48.0134966005524, 51.0042204168578, 53.9320608459731,&
!    56.7833601290163, 59.5482728298363, 62.2206194320588, 64.7975706790533,&
!    67.2792279882741, 69.6681589022773, 71.9689341802758, 74.1876988925387,&
!    76.331796630125, 78.4094552099168, 80.4295344892687, 82.4013318763434,&
!    84.3344382773342, 86.2386359278657, 88.123829239849, 90./)


!  do j=1,nlat
!    mlat(j) = 50. + float(j) 
!  enddo

! mlat = (/50.,60.,70.,80./)
! mlt =  (/0.,6.,12./)

  do i=1,nlon
!    !mlt(i) = float(i-1)
!    !lon(i) = mlt(i)*15.
     mlt(i) = mlon(i)/15.
!!     print *,i,mlon(i),mlt(i)
!    mlt(i) = float(i-1)/float(nlon-1)*24.
!    lon(i) = mlt(i)*360./float(nlon-1)
  enddo
  do j=1,nlat
    do i=1,nlon
      call epotval(mlat(j),mlt(i),fill,epot(i,j))
!      if (itime.eq.1 .and.j.eq.nlat) print *,i,epot(i,j)
      !if (epot(i,j) /= fill) &
      !write(6,"('test_w05sc: mlat=',f8.3,' mlt=',f8.3,' epot=',1pe12.4)") &
      !      mlat(j),mlt(i),epot(i,j)
    enddo
  enddo

 
  do i = 1,nlon
    j  = nlat - 2
    epot(i,j) = 0.5 * (epot(i,j-1)+epot(i,j+1))
  enddo

!    i = nlon
!    j  = nlat - 2
!    epot(i,j) = 0.5 * (epot(i,j-1)+epot(i,j+1))


  if ( itime .eq. 1 ) then
    open(141,file='epot.dat',form='unformatted')
    write(141) epot
    close(141)
  endif



!
! Reset model for magnetic potential:
!
  call setmodel(by,bz,tilt,swvel,swden,file_path,'bpot')
!
! Call mpfac to calculate field-aligned current from magnetic potential model:
!
  do j=1,nlat
    do i=1,nlon
      call mpfac(mlat(j),mlt(i),fill,fac(i,j))
!     if (fac(i,j) /= fill) &
!     write(6,"('test_w05sc: mlat=',f8.3,' mlt=',f8.3,' fac=',1pe12.4)") &
!       mlat(j),mlt(i),fac(i,j)
    enddo
  enddo

!
! Save epot to ascii file, to be read and plotted by IDL:
!
!!  lu = 20
!!  datfile = 'wei05sc_epot_f90.dat '
!!  open(lu,file=datfile,status='replace')
!!  write(lu,"(2i4)") nlon,nlat
!!  write(lu,"(6f8.2)") mlt
!!  write(lu,"(6f18.13)") mlat
!!  write(lu,"(6e12.4)") epot
!!  close(lu)
!!  write(6,"('Wrote ascii file ',a)") datfile

  do j = 1,nlon
    do i = 1,nlat
      lat_sami(j,i) = mlat(i)
      lon_sami(j,i) = mlt(j) * 15.
      pot_sami(j,i) = epot(j,i)
      phi_weimer(i,j) = epot(j,i)
    enddo
  enddo

!  open(101,file='lat_grid.dat',form='formatted')
!  write(101,*) lat_sami
!  close(101)

!  open(101,file='lon_grid.dat',form='formatted')
!  write(101,*) lon_sami
!  close(101)

!  open(101,file='weimer_pot.dat',form='formatted')
!  write(101,*) pot_sami
!  close(101)

  write(101) hrutw
  write(101) phi_weimer * 1.e3 / 300. ! convert to statvolt/cm for sami3

  enddo  

  close(122)
!  close(123)

!
! Save fac to ascii file, to be read and plotted by IDL:
!

! skip this for now

!  datfile = 'wei05sc_fac_f90.dat '
!  open(lu,file=datfile,status='replace')
!  write(lu,"(2i4)") nlon,nlat
!  write(lu,"(6f8.2)") mlt
!  write(lu,"(6f18.13)") mlat
!  write(lu,"(6e12.4)") fac
!  close(lu)
!  write(6,"('Wrote ascii file ',a)") datfile


  close(101)

!  open(101,file='phi_weimer.inp',form='unformatted')
!  read(101) hrutw
!  read(101) phi_weimer
!  print *,'hrutw....',hrutw
!  print *,'phi_weimer',phi_weimer
!  read(101) hrutw
!  print *,'hrutw....1',hrutw
!  close(101)


end program test_w05sc_mod

