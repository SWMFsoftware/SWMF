c-----=================================================================
c     =================================================================
      subroutine claw2ez(FirstTime,T1,T2,xlower,xupper,ylower,yupper,
     &     iterations,mbcIN,mxIN,myIN,
     &     qIN,aIN,bIN,cIN)
c     =================================================================
c
c     An easy-to-use clawpack driver routine for simple applications
c
c     Author: Randall J. LeVeque
c     Version of August, 1999 --  CLAWPACK Version 4.0
c
c     Modified: Darren De Zeeuw
c       changed to work as advection solver for use with RCM
c       August 2000
c
      implicit double precision (a-h,o-z)
      external bc2,rpn2,rpt2,src2,b4step2

      common /compsi/ pi
      common /comvt/ tperiod,pi2

      real, dimension(1-mbcIN:mxIN+mbcIN, 1-mbcIN:myIN+mbcIN) :: 
     &     qIN,aIN,bIN,cIN


      parameter (maxmx  =   360)
      parameter (maxmy  =   120)
      parameter (mwork  = 60000)
      parameter (mbc    =     2)
      parameter (meqn   =     1)
      parameter (mwaves =     1)
      parameter (maux   =     3)

      dimension  aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)
      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension work(mwork)
      dimension mthlim(mwaves)
      dimension method(7),dtv(5),cflv(4),nv(2),mthbc(4)

      logical FirstTime

!     save q,aux,work,mx,my,tstart,tend,dtv,
!    &     cflv,nv,method,mthlim,mthbc,info
!
      SAVE
c
c     ----------------------------------------
c-----| Check for correctness of input values:
c     ----------------------------------------
c
      if (mbcIN.ne.mbc) then
        write(6,*) '*** ERROR ***  mbcIN'
        write(6,*) 'require mbcIN = mbc,  mbcIN=',mbcIN,'  mbc=',mbc
        stop 
      endif
      if (mxIN.gt.maxmx)then
        write(6,*) '*** ERROR ***  mxIN'
        write(6,*) 'require mxIN <= maxmx,  mxIN=',mxIN,'  maxmx=',maxmx
        stop 
      endif
      if (myIN.gt.maxmy)then
        write(6,*) '*** ERROR ***  myIN'
        write(6,*) 'require myIN <= maxmy,  mxIN=',mxIN,'  maxmx=',maxmx
        stop 
      endif

      mx = mxIN                 !cells in x direction
      my = myIN                 !cells in y direction

c
c     --------------------------------
c-----| First time in initializations:
c     --------------------------------
c
      if (FirstTime) then
c
c       --------------------------
c-------| Set the input variables:
c       --------------------------
c
c-------compute values for common blocks:
c       pi used in psi.f
        pi = 4.d0 * datan(1.d0)
c       2*pi used in b4step2:
        pi2 = 2.d0*pi

c-------timestepping variables
!       dtv(1)  =   1.0d0       !initial dt (used in all steps if method(1)=0)
        dtv(1)  =   T2-T1       !initial dt (used in all steps if method(1)=0)
        dtv(2)  =   T2-T1       !max allowable dt
        cflv(1) =   1.0d0       !max allowable Courant number
        cflv(2) =   0.9d0       !desired Courant number
        nv(1)   =   500        !max number of time steps per call to claw2

c-------input parameters for clawpack routines
        method(1) = 1           !1 for variable dt,   = 0 for fixed dt
        method(2) = 2           !order
        method(3) = 2           !transverse order
        method(4) = 0           !verbosity of output
        method(5) = 1           !source term splitting
        method(6) = 0           !mcapa
        method(7) = maux        !maux

        mthlim(1:mwaves) =2     !limiter for each wave  (mw=1,mwaves)

        mthbc(1) = 1            !type of boundary conditions at left
        mthbc(2) = 1            !type of boundary conditions at right
        mthbc(3) = 2            !type of boundary conditions at bottom
        mthbc(4) = 2            !type of boundary conditions at top

!        tperiod = 2.0d0
        tperiod = 0.0d0

c-------check for correctness of input values
        if ((mthbc(1).eq.2 .and. mthbc(2).ne.2) .or.
     &       (mthbc(2).eq.2 .and. mthbc(1).ne.2)) then
          write(6,*) '*** ERROR ***  periodic boundary conditions'
          write(6,*) 'require mthbc(1) and mthbc(2) BOTH be set to 2'
          stop 
        endif
        
        if ((mthbc(3).eq.2 .and. mthbc(4).ne.2) .or.
     &       (mthbc(4).eq.2 .and. mthbc(3).ne.2)) then
          write(6,*) '*** ERROR ***  periodic boundary conditions'
          write(6,*) 'require mthbc(3) and mthbc(4) BOTH be set to 2'
          stop 
        endif

c-------check that enough storage has been allocated:
        if (method(5).lt.2) then
          narray = 1            !# only need one qwork array
        else
          narray = 2            !# need two qwork arrays for Strang splitting
        endif
        
        maxm = max0(maxmx, maxmy)
        mwork1 = (maxm+2*mbc)*(10*meqn + mwaves + meqn*mwaves 
     &       + 3*maux + 2) 
     &       + narray * (maxmx + 2*mbc) * (maxmy + 2*mbc) * meqn   

        if (mx.gt.maxmx .or. my.gt.maxmy .or. mwork.lt.mwork1) then
c---------# insufficient storage
          maxmx1 = max0(mx,maxmx)
          maxmy1 = max0(my,maxmy)
          maxm1 = max0(maxmx1,maxmy1)

          mwork1 = (maxm1+2*mbc)*(10*meqn + mwaves + meqn*mwaves
     &         + 3*maux + 2)
     &         + narray * (maxmx1 + 2*mbc) * (maxmy1 + 2*mbc) * meqn

          write(6,*) ' '
          write(6,*) '*** ERROR *** Insufficient storage allocated'
          write(6,*) 'Recompile after increasing values in claw2ex.f:'
          write(6,611) maxmx1
          write(6,612) maxmy1
          write(6,613) mwork1
 611      format(/,'parameter (maxmx = ',i5,')')
 612      format('parameter (maxmy = ',i5,')')
 613      format('parameter (mwork = ',i7,')',/)
          stop
        endif

!        write(6,*) ' '
!        write(6,*) ' finished with CLAW initializations ...'
!        write(6,*) ' '

      end if
c
c     --------------------------
c-----| Copy arrays in to solve:
c     --------------------------
c
c-----copy q
      q(1-mbc:mx+mbc,1-mbc:my+mbc,1) = qIN(1-mbc:mx+mbc,1-mbc:my+mbc)
c
c-----set auxiliary arrays 
c     #   aux(i,j,1) is edge velocity at "left" boundary of grid point (i,j)
c     #   aux(i,j,2) is edge velocity at "bottom" boundary of grid point (i,j)
      DO i=1,mx
      DO j=1,my
          aux(i,j,1) = aIN(i,j)
          aux(i,j,2) = bIN(i,j)
          aux(i,j,3) = cIN(i,j)
      END DO
      END DO

c-----grid spacing
      dx = (xupper - xlower) / float(mx)
      dy = (yupper - ylower) / float(my)
      
c
c     ------------
c-----| Main loop:
c     ------------
c
      tstart = T1
      tend   = T2
c
      call claw2(maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &     q,aux,xlower,ylower,dx,dy,tstart,tend,dtv,
     &     cflv,nv,method,mthlim,mthbc,
     &     work,mwork,info,bc2,rpn2,rpt2,src2,b4step2)

c-----use final dt as starting value on next call
      dtv(1) = dtv(5)

c-----check to see if an error occured:
      if (info .ne. 0) then
        write(6,*) '*** ERROR in claw2 ***  info =',info
        STOP
      endif

c-----write out information about this call to claw:
!      write(6,1010) info,nv(2),tend,
!     &     dtv(3),dtv(4),dtv(5),
!     &     cflv(3),cflv(4)
! 1010 format('CLAW Summary: info=',i2,', steps=',i4,', tend=',d9.3,/,
!     &     '    min  dt=',d9.3,',  max  dt=',d9.3,', last dt=',d9.3,/,
!     &     '    max cfl=',d9.3,', last cfl=',d9.3)

!      write(6,1020) nv(2),cflv(4),dtv(3),tend
! 1020 format('CLAW Summary:  steps=',i4,
!     &     ', last cfl=',d9.3,', min dt=',d9.3,', tend=',d9.3)

c
c     ----------------------------
c-----| Copy array back to return:
c     ----------------------------
c
      qIN(1-mbc:mx+mbc,1-mbc:my+mbc) = q(1-mbc:mx+mbc,1-mbc:my+mbc,1)
      iterations = nv(2)
!
      RETURN 
      END SUBROUTINE Claw2ez
