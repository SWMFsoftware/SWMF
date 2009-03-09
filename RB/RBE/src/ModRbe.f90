module rbe_io_unit

  integer :: iUnit1, iUnit2

end module rbe_io_unit
!=============================================================================
module rbe_constant
  implicit none

  real, parameter ::  &
       xme=7.9e15,    &	            ! magnetic dipole moment of the Earth
       re=6.371e6,    &	            ! earth's radius (m)
       xmp=1.673e-27, &             ! mass of H+ in kg
       q=1.6e-19,     &	            ! electron charge
       c=2.998e8                    ! speed of light (m/s)

end module rbe_constant
!=============================================================================
module rbe_grid
  implicit none

  ! max no. of SW-IMF, dst, Kp data pts
  integer, parameter :: nswmax=40000, ndstmax=45000, nKpmax=250  

  ! ir=no. of grids in latitude in ionosphere, ip=no. of grids in
  ! local time, iw = no. of grids in magnetic moment, ik = no. of
  ! grids in invariant K, ns = no. of species, je = no. of 
  ! pts. in the fixed E grid, ig = no. of pts. in the fixed y grid.
  !integer, parameter :: ns=2, ir=51, ip=48, iw=29, ik=29, je=12, ig=12


  ! full version:
  integer, parameter :: ns=2, ir=51, ip=48, iw=29, ik=29, je=12, ig=12
  
  ! quick version:
  ! integer, parameter :: ns=2, ir=51, ip=48, iw=15, ik=15, je=9,  ig=6     
  
  ! dimension of Meredith's data:
  integer,parameter  :: irw=7, ipw=24
  ! dimension of Horne's data:
  integer,parameter  :: irc=5,  ipe=5, iwc=6, ipa=91

end module rbe_grid
!=============================================================================
module rbe_cread1

  implicit none

  character(len=8):: outname
  logical         :: UseSeparatePlotFiles = .false.
  character(len=2):: st2

end module rbe_cread1
!=============================================================================
module rbe_cread2

  use rbe_grid

  implicit none

  real,parameter :: dtmax=3.0
  real :: &
       tstart,tmax,trans,tint,tf,hlosscone,dsth(ndstmax),&
       tdst(ndstmax),timf(nswmax),bxw(nswmax),byw(nswmax),bzw(nswmax),&
       tsw(nswmax),xnswa(nswmax),vswa(nswmax),rc,tKp(nKpmax),&
       xKph(nKpmax)
  integer :: &
       nimf,nsw,iyear,iday,js,itype,nstep,nstept,ndst,nKp,&
       ires,ismo,imod,iprint,ntime,iconvect,init,il,ie,idfa,idfe,iplsp

  character (len=8)::  storm

  logical :: IsStandAlone=.false.,UseGm=.false.,UseIE=.false., &
       UseSplitting = .false.,DoSaveIe = .false.

  logical :: UseMcLimiter = .false.,UseCentralDiff = .false.
  real    :: BetaLimiter  = 2.0

end module rbe_cread2
!=============================================================================
module rbe_cgrid
  use rbe_grid
  implicit none
  real :: xlati(ir),dlati(ir),phi(ip),dphi,xmlt(ip),&
       xmass(ns),si(0:ik+1),ds(ik),rsi,w(0:iw+1),dw(iw),rw,&
       d4(ir,iw,ik),xjac(ir,iw),gride(je),gridp(je),gridy(ig),colat(ir),&
       xmltd(ip)
end module rbe_cgrid
!=============================================================================
module rbe_cfield
  use rbe_grid
  implicit none
  real :: bo(ir,ip),ro(ir,ip),xmlto(ir,ip),y(ir,ip,0:ik+1),&
       Hdens(ir,ip,ik),p(ir,ip,iw,ik),v(ir,ip,iw,ik),ekev(ir,ip,iw,ik),&
       rmir(ir,ip,ik),tcone(ir,ip,iw,ik),tanA2(ir,ip,0:ik+1),&
       volume(ir,ip),bm(ir,ip,ik),gamma(ir,ip,iw,ik),parmod(10),rb,&
       xo(ir,ip),yo(ir,ip),tya(ir,ip,0:ik+1),gridoc(ir,ip)

  integer :: irm(ip),irm0(ip),iba(ip)

end module rbe_cfield
!=============================================================================
module rbe_ccepara
  use rbe_grid
  implicit none
  ! Charge exchange loss rate
  real :: achar(ir,ip,iw,ik)
end module rbe_ccepara
!=============================================================================
module rbe_convect
  use rbe_grid
  implicit none
  real :: potent(ir,ip),xnsw,vsw,Bx,By,Bz
end module rbe_convect
!=============================================================================
module rbe_cVdrift
  use rbe_grid
  implicit none
  ! Drift velocity (lon, poloidal)
  real:: vl(ir,ip,iw,ik),vp(ir,ip,iw,ik)
end module rbe_cVdrift
!=============================================================================
module rbe_cinitial
  use rbe_grid
  implicit none
  real :: f2(ir,ip,iw,ik),xnsw0,vsw0,Bx0,By0,&
       Bz0,vswb0,xnswb0,elb,eub,e_l(ir),ecbf(ir),ecdt(ir),&
       eclc(ir),ecce(ir)
  integer :: iw1(ik),iw2(ik)
end module rbe_cinitial
!=============================================================================
module rbe_cboundary
  use rbe_grid
  implicit none
  real :: fb(ip,iw,ik),vswb,xnswb
end module rbe_cboundary
!=============================================================================
module rbe_plasmasphere
  use rbe_grid
  implicit none
  ! plasmasphere density(m^-3) from pbo_2.f
  real :: par(2),density(ir,ip) 
end module rbe_plasmasphere
!=============================================================================
module rbe_time
  real :: t, dt
  integer :: istep
end module rbe_time
!=============================================================================
Module ModGmRb
  use rbe_grid,ONLY: nLat => ir, nLon => ip
  real, allocatable :: StateLine_VI(:,:),StateIntegral_IIV(:,:,:)
  integer :: iLineIndex_II(nLon,1:nLat),nPoint
  
end Module ModGmRb
!=============================================================================
Module ModChorusIntensity
  use rbe_grid, ONLY: irw, ipw
  real :: wLshell(irw),wmlt(ipw),chorusI(irw,ipw,3)
end Module ModChorusIntensity
!=============================================================================
Module ModChorusDiffCoef
  use rbe_grid, ONLY: irc, ipe, iwc, ipa
  real :: cLshell(irc),ompea(ipe),ckeV(iwc),cPA(ipa),&
          cDEE(irc,ipe,iwc,ipa),cDaa(irc,ipe,iwc,ipa)
end Module ModChorusDiffCoef
!=============================================================================
Module ModWpower
use rbe_grid, ONLY: ir,ip

real :: CHpower(ir,ip),ompe(ir,ip)

end Module ModWpower
!=============================================================================
!=============================================================================
Module ModRbTime
  use ModKind,only:real8_

  real (kind=real8_)   ::  StartTime,CurrentTime !needed for SWMF stuff
  integer              ::  iStartTime_I(7)  =(/1976,6,28,0,0,0,0/)
  integer              ::  iCurrentTime_I(7)=(/1976,6,28,0,0,0,0/)
!  integer              ::  iDOY,IYD
end Module ModRbTime
!=============================================================================
