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

      real rx,ry,rz,vx,vy,vz,bx,by,bz,dd,pp
      integer  iMax,nMax
      PARAMETER(nMax=2500)
      common /ihcoord/ rx(nMax),ry(nMax),rz(nMax)
      common /ihvel/   vx(nMax), vy(nMax), vz(nMax)
      common /ihmagf/ bx(nMax), by(nMax), bz(nMax)
      common /ihpdi/     dd(nMax), pp(nMax), iMax
      real vxOld,vyOld,vzOld
      common /oldvel/     vxOld(nMax), vyOld(nMax), vzOld(nMax)
      integer Old_,New_
      parameter(Old_=1,New_=2)
      integer iShock,iShockOld
      common/ishock/iShock,iShockOld
      real Smooth_VII
      common/smooth/Smooth_VII(11,nMax,Old_:New_)
      integer nResolution
      PARAMETER(nResolution=10)
