!The data which are accepted by SEP module from the corona or
!inner heliosphere
!with the identifier old: storage to save the data got in the
!course of previous coupling
!
!rx,ry.rz - the coordinates of the Lagrangian mesh, in au
!
!vx,vy,vz - three components of the velocity in
!the Lagrangian point, in au/hour
!
!bx,by,bz - three components of the magnetic field in 
!the Lagrangian point, in nT
!
!dd - density in nucleons/cm^3
!
!pp - pressure in erg/cm^3

      real rx,ry,rz,vx,vy,vz,bx,by,bz,dens,pres
      include 'param.h'
      integer  iMax,nMax
      PARAMETER(nMax=nRMax)
      common /SP_ihcoord/ rx(nMax),ry(nMax),rz(nMax)
      common /SP_ihvel/   vx(nMax), vy(nMax), vz(nMax)
      common /SP_ihmagf/ bx(nMax), by(nMax), bz(nMax)
      common /SP_ihpdi/     dens(nMax), pres(nMax), iMax
      real vxOld,vyOld,vzOld
      common /SP_oldvel/     vxOld(nMax), vyOld(nMax), vzOld(nMax)
      real  algbb(0:nMax)
      real  algll(0:nMax)
      real  algnn(0:nMax)
      real  vr(0:nMax)
      common /SP_plasma/ algbb,algll,algnn,vr
      integer Old_,New_
      parameter(Old_=1,New_=2)
      integer iShock,iShockOld
      common/SP_ishock/iShock,iShockOld
      real Smooth_VII
      common/SP_smooth/Smooth_VII(11,nMax,Old_:New_)
      integer nResolution
      PARAMETER(nResolution=10)
      real rTransient
      integer iTransient
      common/SP_trans/rTransient,iTransient
      logical DoRestart
      common/SP_restart/DoRestart
