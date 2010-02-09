  !----------------------------------------------------------------------------
  !
  !                              crcm.f90
  !
  ! A test program to call subroutine crcm
  !
  ! The CRCM subroutine called from the OpenGGCM MHD model. 
  ! The first call (ijob=0), the CRCM reports its grid dimensions (np,nt).
  ! The second call (ijob=1), the CRCM does its initialization.
  ! On Each subsequent call (ijob.ge.2), the CRCM performs calculation and output
  !   ring current phase space density and pressure at the end of the time
  !   interval: time+delta_t.
  !
  ! This code must be compiled with default double precision (i.e., -r8)
  !
  ! Created on 29 December 2006 by Mei-Ching Fok, Code 673, NASA GSFC.
  !-----------------------------------------------------------------------------
  
  implicit none
  
  integer,parameter :: npmax=500    ! max number of points in tracing a field line
  ! this number can be changed
  integer ijob,np,nt,neng,npit,nspec,i,j
  integer,allocatable,dimension(:) :: kspec
  real delta_t,time,Hiono
  real pi,halfPCP,pcb,shield,cosp,cosl,pot1
  real,allocatable,dimension(:) :: xlat,xmlt,energy,sinAo
  real,allocatable,dimension(:,:) :: pot,rrio,ttio,fac,phot,brad,bphi,beq,ftv
  integer,allocatable,dimension(:,:) :: npt
  real,allocatable,dimension(:,:,:) :: bb
  real,allocatable,dimension(:,:,:,:) :: gb
  real,allocatable,dimension(:,:,:,:,:) :: flux

  ! These parameters can be changed
  delta_t=300.            ! in second. OpenGGCM and CRCM talk every delta_t
  Hiono=120.              ! ionosphere altitude in km in OpenGGCM 

  ! Initial calls of CRCM (ijob=0,1)
  do ijob=0,1  ! CRCM reports grid dimensions (ijob=0), initialization (ijob=1)
     call crcm(ijob,np,nt,neng,npit,nspec,delta_t,Hiono, &
          npmax,npt,gb,bb,brad,bphi,beq,ftv,pot,rrio,ttio,  &
          kspec,xlat,xmlt,energy,sinAo,fac,phot,flux)
     if (ijob.eq.0) then  ! query np,nt,neng,npit,nspec and allocate arrays
        allocate (kspec(nspec),xlat(np),xmlt(nt),energy(neng),sinAo(npit))
        allocate (pot(np,nt),rrio(np,nt),ttio(np,nt),fac(np,nt),phot(np,nt))
        allocate (brad(np,nt),bphi(np,nt),beq(np,nt),ftv(np,nt))
        allocate (npt(np,nt),gb(np,nt,npmax,3),bb(np,nt,npmax))
        allocate (flux(nspec,np,nt,neng,npit))
     endif
  enddo

  ! Start the time loop. For this test, just do one-hour run 
  time=0.
  do ijob=2,13
     ! update magnetic field input to the CRCM: npt,gb,bb,brad,bphi,beq,ftv
     open(unit=11,file='Bfield.dat',status='old')
     read(11,*) npt     ! npt(np,nt)
     read(11,*) gb      ! gb(np,nt,npmax,3), SM coordinates in RE
     read(11,*) bb      ! bb(np,nt,npmax), B at gb in nT
     read(11,*) brad    ! brad(np,nt), equatorial distance in RE
     read(11,*) bphi    ! bphi(np,nt), magnetic local time in hour
     read(11,*) beq     ! beq(np,nt), B at magnetic equator in nT
     read(11,*) ftv     ! ftv(np,nt), flux tube volume per unit flux in RE**3/Wb)
     close(11)

     ! Update ionospheric potential in Volt, pot. In this test run,
     ! Volland-Stern type potential is assumed
     pi=acos(-1.)
     halfPCP=30000.        
     pcb=70.           ! polar cap boundary latitude in degree
     shield=2.         ! shielding factor
     cosp=cos(pcb*pi/180.)
     do i=1,np
        cosl=cos(xlat(i)*pi/180.)
        pot1=halfPCP*(cosp*cosp/cosl/cosl)**shield
        do j=1,nt
           pot(i,j)=pot1*sin(xmlt(j)*pi/180.)
        enddo
     enddo

     ! Update density and temperature mapped to the ionosphere at the CRCM grid.
     ! In this test run, they are assumed to be constants
     rrio=5.e5       ! number density in m^3 mapped to the ionosphere
     ttio=5.e3       ! temperature in eV mapped to the ionosphere

     call crcm(ijob,np,nt,neng,npit,nspec,delta_t,Hiono, &
          npmax,npt,gb,bb,brad,bphi,beq,ftv,pot,rrio,ttio,  &
          kspec,xlat,xmlt,energy,sinAo,fac,phot,flux)
     time=time+delta_t
  enddo

  ! Print results
  open(unit=2,file='crcm.output')
  write(2,*) time
  write(2,*) xlat
  write(2,*) xmlt
  write(2,*) energy
  write(2,*) sinAo
  write(2,*) brad
  write(2,*) bphi
  write(2,*) beq 
  write(2,*) rrio
  write(2,*) ttio
  write(2,*) fac 
  write(2,*) phot
  write(2,*) flux
  close(2)

end program     ! end of the test program
!=============================================================================
subroutine CON_stop(StringError)
  implicit none
  character (len=*), intent(in) :: StringError

  ! Local variables:
  integer :: iError,nError
  !----------------------------------------------------------------------------

  write(*,*) 'Stopping execution with msg:'
  write(*,*) StringError
  stop

end subroutine CON_stop
!=============================================================================
subroutine CON_set_do_test(String, DoTest, DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe

  DoTest = .false.; DoTestMe = .false.

end subroutine CON_set_do_test

