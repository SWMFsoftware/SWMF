c**************************************************************
c     as of  February 17, 2004   3:05 p.                      *
c     SIMILARITY SOLUTION with B_fi ADDED                     *
c     CONCEPT: insert one M-config and assume similarity      *
c              lambda + time prescribed as wished  --         *
c              variable step-size to speed up - add later??   *
c     this one employs VARIABLE omega,B_r                     *
c     Change so that dbdt, dldt etc do mean that they do mean *
c     Take care of B definition - so that have real line      *
c     This version loses particles at outer spatial boundary  *
c     Concept : no inward flux at i=0 (outw irelevant)        *
c               no inward flux at i=nr+1, keep i=nr           *
c               no inward flux at i=nr+1, keep i=nr           *
c     PEEL-OFF  Routine added - but not used                  *
c**************************************************************
      subroutine sp_set_defaults
      implicit none
      include 'coupler.h'
      integer nr,nmu,nw
      real dim1
      common /SP_size  / nr,nmu,nw, dim1
      real wghtl,wghtmu,wghtw
      common /SP_suly /  wghtl,wghtmu,wghtw
      real  time,tmax,dlnt0,dlnt1,dta
      integer kfriss,kacc
      common /SP_times/  time,tmax,dlnt0,dlnt1,dta,kfriss,kacc
      real emin,emax,eein,ee
      common /SP_energy/ emin,emax,eein,ee(0:nPMax)
      integer kinj
      real einj,pinj,qqp(0:nPMax),qqx(0:nRMax)
      common /SP_quelle/ kinj,einj,pinj,qqp,qqx
      integer iz,massa
      real ekpp,xlmbda0,omega,xscatt1
      common /SP_partid/ iz,massa,ekpp,xlmbda0
!      common /scphys/ wind,omega,xscatt1 !'wind' is undefined
      common /SP_scphys/ omega,xscatt1
      real qex,cmu(nMuMax),scmu(nMuMax)
      real wsc(0:nPMax),xsc(0:nRMax)
      common /SP_scatti/ qex,cmu,scmu,wsc,xsc
      integer krMax,krObs,kEMax,kEObs
      real rObs,eObs
      common /SP_obsrad/ krmax,krobs(5),robs(5)
      common /SP_obserg/ kemax,keobs(5),eobs(5)
      real period
      !-----------------------------------------------------
      nr=nRMax      !1000
      nMu=nMuMax    !10
      nW=nPMax      !150
      wghtl=1.0
      wghtmu=1.0
      wghtw=1.0
      tMax=0.0      
      emin=0.001
      emax=1000.
      einj=0.005
      kFriss=5
      kAcc=5
      period=26.0 !Days
      iz=1
      massa=1
      ekpp=0.33
      xlmbda0=0.4
      qex=0.0
      krMax=4
      kEMax=5
      rObs(1)=0.05
      rObs(2)=0.1
      rObs(3)=0.2
      rObs(4)=1.0
      kEObs(1)=25
      kEObs(2)=50
      kEObs(3)=75
      kEObs(4)=100
      kEObs(5)=125
      if (period.gt.1.e3) then
         omega = 0.
      else
         omega = 4.*asin(1.)/period/24.
      end if
   
      if (xlmbda0.gt.1.e3) then
         xscatt1 = 0.
      else
         xscatt1 = ((3-qex)*(1-qex)/3.)/xlmbda0
      end if
      rTransient = 1.0 ![RSun]
      end subroutine sp_set_defaults
      
      subroutine sp_init
      implicit none
      include 'param.h'
      integer::nr,nmu,nw,kfriss,kacc,nn
      real:: wghtl,wghtmu,wghtw,time,tmax,dlnt0,dt1,dta
      real:: slamb,tblast,tblst1,rshck1,dlnt, dim1
      real:: rmin,rshock,rmax,r(0:nRMax),wind,omega,xscatt1
      real::  ggamma,bbrad,vvmin,vvmax,ddmin,ddmax,
     1                ccmin,ccmax,aamin,aamax,bbmin,bbmax
      common /SP_size  / nr,nmu,nw, dim1        
      common /SP_suly  / wghtl,wghtmu,wghtw
      common /SP_times/  time,tmax,dlnt0,dt1,dta,kfriss,kacc
      common /SP_blast/  slamb,tblast,tblst1,rshck1,dlnt
      common /SP_radio / nn,rmin,rshock,rmax,r
!      common /scphys/ wind,omega,xscatt1 !'wind' is undefined
      common /SP_scphys/ omega,xscatt1
      common /SP_gazdi/  ggamma,bbrad,vvmin,vvmax,ddmin,ddmax,
     1                ccmin,ccmax,aamin,aamax,bbmin,bbmax
      include 'stdout.h'
      integer::jnext,jstep,jsep,kstep
      integer::iStep
      common /SP_spmain/jnext,jsep,jstep,istep
c ----------------------------------------------------------
      jnext = 2
       if(DoWriteAll)write(iStdout,*)prefix,
     1 'first j-stop at ',jnext
      call SP_cool   
      call SP_pangle  
c ----------------------------------------------------------
c     SP_initial conditions
      jstep = 0
      iStep=jStep*kfriss
      if(DoWriteAll)write(iStdout,*)prefix,
     1 'want to start from negative j-step : ',jstep
      jsep = 1
      if(DoWriteAll)write(iStdout,*)prefix,
     1 'Particles be calculated from j-sep : ',jsep
      call SP_initial
      call SP_peeloff    
      call SP_opentime
c ----------------------------------------------------------
      end subroutine sp_init
!===========================================================!
!Data from heliosphere and/or solar corona model are
!obtained through the common blocks which are put into
!the header file coupler.h
!==================coupler.h begin:==========================!
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
!dens - density in nucleons/cm^3
!
!pres - pressure in erg/cm^3
!
!      real rx,ry,rz,vx,vy,vz,bx,by,bz,dens,pres
!      integer iMax,nMax
!      PARAMETER(nMax=2500)
!      common /ihcoord/ rx(nMax),ry(nMax),rz(nMax)
!      common /ihvel/   vx(nMax), vy(nMax), vz(nMax)
!      common /ihmagf/ bx(nMax), by(nMax), bz(nMax)
!      common /ihpdi/     dens(nMax), pres(nMax), iMax
!      real vxOld,vyOld,vzOld
!      common /oldvel/     vxOld(nMax), vyOld(nMax), vzOld(nMax)
!      integer Old_,New_
!      parameter(Old_=1,New_=2)
!      integer iShock,iShockOld
!      common/ishock/ishock/
!      real Smooth_VII
!      common/smooth/Smooth_VII(11,nMax,New_)
!      integer nResolution
!      PARAMETER(nResolution=10)
!=======================end coupler.h======================!
      subroutine SP_clean_coupler
      implicit none
      include 'coupler.h'
      include 'stdout.h'
      integer iR,iV
      real cZero
      parameter(cZero=0.00000000000000000000000000000000)
      do iR=1,nMax
         rx(iR)=cZero
         ry(iR)=cZero
         rz(iR)=cZero
         dens(iR)=cZero
         vx(iR)=cZero
         vy(iR)=cZero
         vz(iR)=cZero
         bx(iR)=cZero
         by(iR)=cZero
         bz(iR)=cZero
         pres(iR)=cZero
         vxOld(iR)=cZero
         vyOld(iR)=cZero
         vzOld(iR)=cZero
         do iV=1,11
            Smooth_VII(iV,iR,Old_)=cZero
            Smooth_VII(iV,iR,New_)=cZero
         end do
      end do
      if(DoWriteAll)write(iStdOut,*)prefix,'Nullify the coupler'
      end subroutine SP_clean_coupler

!-----------------------------------------------------------!
      subroutine SP_MASTER(tSimulation,tFinal)
      implicit none
      real tSimulation 

! intent(inout)
!"In" - the start value of the physical time
!"Out" - the final value of the physical time

      real tFinal
! intent (in)
! The value of the physical time to stop the run
! After the physical time exceeds this value the run stops.
! tFinal DIFFERS from tmax, which is the value of the INTERNAL
! selfsimilar time to stop the run

      real tSimulationStart
      real tSimilarityStart
      include 'param.h'
      integer nr,nmu,nw,kfriss,kacc,nn
      real  wghtl,wghtmu,wghtw,time,tmax,dlnt0,dt1,dta
      real  slamb,tblast,tblst1,rshck1,dlnt, dim1
      real  rmin,rshock,rmax,r(0:nRMax),wind,omega,xscatt1
      real   ggamma,bbrad,vvmin,vvmax,ddmin,ddmax,
     1                ccmin,ccmax,aamin,aamax,bbmin,bbmax
      common /SP_size  / nr,nmu,nw, dim1        
      common /SP_suly  / wghtl,wghtmu,wghtw
      common /SP_times/  time,tmax,dlnt0,dt1,dta,kfriss,kacc
      common /SP_blast/  slamb,tblast,tblst1,rshck1,dlnt
      common /SP_radio / nn,rmin,rshock,rmax,r
!      common /scphys/ wind,omega,xscatt1 !'wind' is undefined
      common /SP_scphys/ omega,xscatt1
      common /SP_gazdi/  ggamma,bbrad,vvmin,vvmax,ddmin,ddmax,
     1                ccmin,ccmax,aamin,aamax,bbmin,bbmax
      logical UseSelfSimilarity,UseRefresh
      common/SP_log/UseSelfSimilarity,UseRefresh
      include 'stdout.h'
      integer::jnext,jstep,jsep,Misc,kstep,iStep
      common /SP_spmain/jnext,jsep,jstep,istep
c ----------------------------------------------------------
      tSimulationStart=tSimulation
      if(UseSelfSimilarity)then
         Misc = kfriss
         call SP_simile(Misc)
         Misc = jstep*kfriss
         tblst1= tblast*exp(float(Misc)*dlnt)
         time = tblst1 
         tSimilarityStart=time
      else
         time=tSimulationStart/3.60e+3
         dt1=(tFinal-tSimulation)/3.60e+3
         dta=1.000001*dt1/float(kacc)
      end if      
c ---------------------------------------------------------- 
      ! Time step loop
 
      do while((time.lt.tmax.or.iStep.ne.(iStep/kfriss)*kfriss.or.
     1    .not.UseSelfSimilarity).and.tSimulation.lt.tFinal)  

!     Do this while internal time less than tMax, 
!     and physical time less than tFinal
!     The first way to stop the run should be disabled in
!     the framework.
!     The second way at the time (05.10.04) is disabled in the 
!     stand alone version
         iStep=iStep+1
         kstep=iStep-kfriss*((iStep-1)/kfriss)
         if(UseSelfSimilarity)then
            tblst1= tblast*exp(float(istep)*dlnt)
            rshck1= rshock*(tblst1/tblast)**slamb
            time = tblst1
            
            tSimulation=tSimulationStart +
     1           (Time-tSimilarityStart)*3.60d3
            
            dt1 = tblst1*dlnt
            dta = dt1/float(kacc)
         end if
         write(iStdOut,900)prefix,
     1        jstep,kstep,istep,time
 900     format(a,5x,'step:',3i10,5x,'t[hour]: ',f12.6)
c     ------------------------------------------ do dynamical step here
         call SP_helios(time,dt1,jstep,kstep)
ccc   if (mod(jstep,200).eq.0.and.kstep.eq.1) read(*,*) lull
         
         if (jstep.ge.jsep) then
         
         ! Loop for acceleration of particles on a static grid
            do  Misc=1,kacc
               
               call SP_source(dta)
               call SP_deltal(dta)
               call SP_lcycle(dta,wghtl)
               call SP_update
               
               call SP_subsub(dta)
               if(.not.UseSelfSimilarity)
     1             tSimulation=tSimulation+dta*3.60e+3
            end do
         end if
c ------------------------------------------    dynamical step done
         jStep=iStep/kFriss
         if(iStep.eq.jStep*kFriss)then
            if (jstep.gt.0) call SP_timevar(jstep,time)
            if (jstep.eq.jnext) then
               write(iStdout,*)prefix,
     1           jstep,'  j-step done -- time: ',time
               if (time.le.1.) Misc = 11
               if (time.gt.1..and.time.le.2.)   Misc=12
               if (time.gt.2..and.time.le.6.)   Misc=13
               if (time.gt.6..and.time.le.12.)  Misc=14
               if (time.gt.12..and.time.le.48.) Misc=15
               if (time.gt.48.) Misc=16
               call SP_alla(Misc,jstep,time)
               call SP_csilla(Misc,istep,time)
               write(iStdout,*)prefix,
     1              'next j-stop ??'
c??   read(*,*)  jnext
               jnext = jnext+24
               write(iStdout,*)prefix,
     1              jnext,'   auto - selected'
            end if
            
            if(UseRefresh)call SP_refresh
         end if
      end do
      ! End of main time step loop
      end subroutine SP_MASTER
!---------------------------------------------------------------!
      subroutine SP_sharpen_and_run(tSimulation,tFinal)
      implicit none
      include 'coupler.h'
      include 'stdout.h'
      real  tSimulation
      real  tFinal
      integer jnext,jstep,jsep,Misc,kstep,iStep
      common /SP_spmain/jnext,jsep,jstep,istep
      integer nr,nmu,nw
      real dim1    
      common /SP_size  / nr,nmu,nw, dim1  
      integer iCoupling,nCoupling,iR
      real DtCoupling
      do iR=1,iMax
         rx(iR)=Smooth_VII( 1,iR,New_) 
         ry(iR)=Smooth_VII( 2,iR,New_)
         rz(iR)=Smooth_VII( 3,iR,New_)
         dens(iR)=Smooth_VII( 4,iR,New_)
         vx(iR)=Smooth_VII( 5,iR,New_)
         vy(iR)=Smooth_VII( 6,iR,New_)
         vz(iR)=Smooth_VII( 7,iR,New_)
         bx(iR)=Smooth_VII( 8,iR,New_)
         by(iR)=Smooth_VII( 9,iR,New_)
         bz(iR)=Smooth_VII(10,iR,New_)
         pres(iR)=Smooth_VII(11,iR,New_)
      end do
      iMax=min(iMax,nR+1)
      nR=iMax-1
      if(DoWriteAll)write(iStdOut,*)prefix,'sharpen and run: ',
     1   'nR,tSimulation,tFinal: ',nR,tSimulation,tFinal
      call SP_get_ishock
      if(tFinal-tSimulation.le.0.00001)then
         !Save the values of log n, log b and so on
         !Do do this, call SP_helios with artificial positive dt
         call SP_helios(0.0,1.0,1,1)
         tSimulation=tFinal
         iMax=0
         return
      end if
      if(iShock<nResolution+iTransient)then
         !Do not sharpen the shock wave
            call SP_get_rshock
            call SP_MASTER(tSimulation,tFinal)
      else
         nCoupling=max(1,iShock-iShockOld)
         if(iShockOld<nResolution+iTransient)then
             nCoupling=1    !No sharpening happened before
         end if
         if(nCoupling.eq.1)then
            call SP_sharpen_profile
            call SP_get_rshock
            call SP_MASTER(tSimulation,tFinal)
            tSimulation=tFinal
         else
            DtCoupling=(tFinal-tSimulation)/nCoupling
            do iCoupling=1,nCoupling
               do iR=1,iMax
                  rx(iR)=(iCoupling*Smooth_VII(1,iR,New_)+
     1                 (nCoupling-iCoupling)*Smooth_VII(1,iR,Old_))/
     2                 nCoupling
                  ry(iR)=(iCoupling*Smooth_VII(2,iR,New_)+
     1                 (nCoupling-iCoupling)*Smooth_VII(2,iR,Old_))/
     2                 nCoupling
                  rz(iR)=(iCoupling*Smooth_VII(3,iR,New_)+
     1                 (nCoupling-iCoupling)*Smooth_VII(3,iR,Old_))/
     2                 nCoupling
                  dens(iR)=(iCoupling*Smooth_VII(4,iR,New_)+
     1                 (nCoupling-iCoupling)*Smooth_VII(4,iR,Old_))/
     2                 nCoupling
                  vx(iR)=(iCoupling*Smooth_VII(5,iR,New_)+
     1                 (nCoupling-iCoupling)*Smooth_VII(5,iR,Old_))/
     2                 nCoupling
                  vy(iR)=(iCoupling*Smooth_VII(6,iR,New_)+
     1                 (nCoupling-iCoupling)*Smooth_VII(6,iR,Old_))/
     2                 nCoupling
                  vz(iR)=(iCoupling*Smooth_VII(7,iR,New_)+
     1                 (nCoupling-iCoupling)*Smooth_VII(7,iR,Old_))/
     2                 nCoupling
                  bx(iR)=(iCoupling*Smooth_VII(8,iR,New_)+
     1                 (nCoupling-iCoupling)*Smooth_VII(8,iR,Old_))/
     2                 nCoupling
                  by(iR)=((iCoupling)*Smooth_VII(9,iR,New_)+
     1                 (nCoupling-iCoupling)*Smooth_VII(9,iR,Old_))/
     2                 nCoupling
                  bz(iR)=(iCoupling*Smooth_VII(10,iR,New_)+
     1                 (nCoupling-iCoupling)*Smooth_VII(10,iR,Old_))/
     2                 nCoupling
                  pres(iR)=(iCoupling*Smooth_VII(11,iR,New_)+
     1                 (nCoupling-iCoupling)*Smooth_VII(11,iR,Old_))/
     2                 nCoupling
               end do
               iShock=iShockOld+iCoupling
               call SP_sharpen_profile
               call SP_get_rshock
               call SP_MASTER(tSimulation,tSimulation+DtCoupling)
            end do
         end if
      end if
      tSimulation=tFinal !To eliminate the error accumulation.
      iMax=0             !To be prepared for the next coupling
      end subroutine sp_sharpen_and_run
!==================================================================
      subroutine SP_get_ishock
      implicit none
      include 'coupler.h'
      include 'stdout.h'
      real DMax,DDens,DDens0,Diff,Rad2,Dist2
      integer i
      real cRsunPerAU,R2
      parameter(cRsunPerAU= 149.59787000000000000/0.696)

      !set iTransient at the last point at which r<rTransient
      i=0;R2=(rTransient/cRSunPerAU)**2
      do while(rx(i+1)**2+ry(i+1)**2+rz(i+1)**2<R2)
         i=i+1
      end do
      iTransient=i
      if(DoWriteAll)write(iStdOut,*)prefix,'iTransient=',iTransient
      iShockOld=iShock
      iShock =iTransient+1
      DMax=0
      DDens = dens(iTransient+1)
      
      do i=iTransient+2,iMax
 

         Dist2 = sqrt((rx(i)-rx(i-1))**2+
     1        (ry(i)-ry(i-1))**2+
     2        (rz(i)-rz(i-1))**2)
         
         if(Dist2==0.0)then
            write(iStdout,*) prefix,'i...=',
     1           i,rx(i),rx(i-1),ry(i),ry(i-1),rz(i),rz(i-1)
         else
            
            DDens0 = DDens
            DDens = dens(i)

C           Find maximum of -(r/rho)*(drho/ds) which is proportional
C           to -1<dr/ds<1 for rho proportional to 1/r**alpha, alpha>=2. 
C           The formula should have a large peak where the field line
C           crosses the shock and drho/ds is a large negative number
C           Note that the field line starts on the Sun and ends at the Earth
            Diff = (DDens0-DDens)/(DDens0+DDens)/Dist2 
     &            * sqrt(rx(i)**2+ry(i)**2+rz(i)**2)
            
            if (Diff.gt.DMax) then
               iShock = i-1
               DMax = Diff
            endif
         endif
      end do
      if(DoWriteAll)then
         write(iStdOut,*)prefix,'iShockOld,iShock=',iShockOld,iShock
      end if
      end subroutine SP_get_ishock
!==================================================================
      subroutine SP_get_rshock
      implicit none
      include 'coupler.h'
      include 'stdout.h'
      real cRsunPerAU
      parameter(cRsunPerAU= 149.59787000000000000/0.696)
      real slamb,tblast,tblst1,rshck1,dlnt
      common /SP_blast/  slamb,tblast,tblst1,rshck1,dlnt
      rshck1=sqrt(rx(iShock)**2+    
     1     ry(iShock)**2+    
     2     rz(iShock)**2)
      write(iStdOut,*)'Shock position is at ',rshck1*cRsunPerAU,' Rs'
      end subroutine SP_get_rshock
!==================================================================    
      subroutine SP_sharpen_profile
      implicit none
      include 'coupler.h'
      integer::nRes2,iUpStream,iDownstream,i
  !-------------------------------------------------
      nRes2=nResolution/2
      iDownstream=iShock-nRes2
      if(iDownStream<1)return
      iUpStream=iShock+nRes2+1
      if(iUpStream>iMax)return
      do i=iDownStream+1,iShock
         dens(i)=dens(iDownStream)
         vx(i)=vx(iDownStream)
         vy(i)=vy(iDownStream)
         vz(i)=vz(iDownStream)
         bx(i)=bx(iDownStream)
         by(i)=by(iDownStream)
         bz(i)=bz(iDownStream)
         pres(i)=pres(iDownStream)
      end do
      do i=iShock+1,iUpStream-1
         dens(i)=dens(iUpStream)
         vx(i)=vx(iUpStream)
         vy(i)=vy(iUpStream)
         vz(i)=vz(iUpStream)
         bx(i)=bx(iUpStream)
         by(i)=by(iUpStream)
         bz(i)=bz(iUpStream)
         pres(i)=pres(iUpStream)
      end do
      end subroutine SP_sharpen_profile 
!=============================================================== 

c ********************* end routine SP_MASTER  *********************

      subroutine SP_refresh
      include 'coupler.h'
      common /SP_size  / nr,nmu,nw, dim1
      common /SP_radio / nn,rmin,rshock,rmax,r(0:nRMax)
      common /SP_times/  time,tmax,dt0,dt1,dta,kfriss,kacc
      common /SP_blast/  slamb,tblast,tblst1,rshck1,dlnt
!      common /scphys/ wind,omega,xscatt1 !'wind' is undefined
      common /SP_scphys/ omega,xscatt1
      common /SP_elem /  zr(0:nRMax),zv(0:nRMax),zp(0:nRMax),zn(0:nRMax)
      common /SP_magia/  qb(0:nRMax),zb(0:nRMax),zt(0:nRMax) 
      common /SP_solutn/ f(0:nRMax,0:nMuMax,0:nPMax),
     1     df(0:nRMax,0:nMuMax,0:nPMax)
      common /SP_smile/  mp,ni
     1       ,eta(0:6000),exx(0:6000),ezz(0:6000),efi(0:6000)
     2       ,evr(0:6000),evx(0:6000),evz(0:6000)
     3       ,edsm(0:6000),edfi(0:6000)
     4       ,ebm(0:6000),ebr(0:6000),ebt(0:6000),ebfi(0:6000)
     5       ,ebx(0:6000),ebz(0:6000),edd(0:6000)
     6       ,dbmds(0:6000),dbfids(0:6000)

c **********************************  concept : *****************
c     shifts back to regions to be considered                   *
c ********************************** ---------- *****************
      
      do i=0,nr-1
         i2=i+1

         algbb(i) = algbb(i2)
         algnn(i) = algnn(i2)
         algll(i) = algll(i2)
         vr(i) = vr(i2)
         do j=0,nmu
            do k=0,nw
               f(i,j,k) = f(i2,j,k)
            enddo
         enddo
      enddo
     
      i = nr
      r(nr) = rshck1*eta(ni)
      vr(nr) = r(nr)/tblst1*evr(ni)
        
      bbr   = ebr(ni)/r(i)**2
      bbm   = ebm(ni)/r(i)**2
      bbfi  = ebfi(ni)/r(i)
      bbb   = sqrt(bbm**2+bbfi**2)

      abb   = alog(bbb)
      ddd   = edd(ni)/r(i)**2
      ann   = alog(ddd)
       
      dll   = rshck1*edsm(ni)*bbb/bbm
      all   = alog(dll)

      algbb(nr) = abb
      algnn(nr) = ann
      algll(nr) = all

      do j=0,nmu
         do k=0,nw
            f(i,j,k) = 0.
         enddo
      enddo

      return 
      end

c ********************  end routine SP_refresh *********************
                                                  !!!! REWORK

c ********************  end routine GASDYN  *********************
c ============================ real calculation cycle starts here ==


      subroutine SP_simile(mm)                 !!!  REWORK !!!
      include 'coupler.h'
      common /SP_size /  nr,nmu,nw, dim1
      common /SP_radio/  nn,rmin,rshock,rmax,r(0:nRMax)
      common /SP_blast/  slamb,tblast,tblst1,rshck1,dlnt
      common /SP_smile/  mp,ni
     1       ,eta(0:6000),exx(0:6000),ezz(0:6000),efi(0:6000)
     2       ,evr(0:6000),evx(0:6000),evz(0:6000)
     3       ,edsm(0:6000),edfi(0:6000)
     4       ,ebm(0:6000),ebr(0:6000),ebt(0:6000),ebfi(0:6000)
     5       ,ebx(0:6000),ebz(0:6000),edd(0:6000)
     6       ,dbmds(0:6000),dbfids(0:6000)

      dimension axx(2500),ayy(2500),azz(2500)
      dimension tint(2500)

      dimension arr(2500),aro(2500)
      dimension vrr(2500),vtt(2500),vro(2500),vzz(2500),vfi(2500)
      dimension brr(2500),btt(2500),bro(2500),bzz(2500),bfi(2500)
      dimension add(2500),app(2500)
      
      real iCoupleHr             !Physical time in hours
      common/SP_couple/tCoupleHr

      include 'stdout.h'
c **************************************  concept: *****************
c     put in a similarity solution shifted so that eta=1 at i=nn   *
c     CONCEPT: rmin-rshock-rmax at start -- read  sim-solution     *
c              then place so that spans the desired region         *
c              Divide into 1001 !! equal dlnt regions              *
c **************************************  -------- *****************




c -------------------------------- beolvasas meglenne
      ausun = 200.
      hour  = 3600.
      cAUKm  = 1.5e8
      gauss = 1.e5
      dens0 = 1.e20
      pres0 = 1.e5
      bb0   = 100.
 
      
      do 30 i=1,iMax
         ax  = rx(i)
         ay  = ry(i)
         az  = rz(i)
        aro(i) = sqrt(ax**2 + ay**2)
        arr(i) = sqrt(ax**2 + ay**2 + az**2)
        axx(i) = ax
        ayy(i) = ay
        azz(i) = az

         ux  = vx(i)
         uy  = vy(i)
         uz  = vz(i)
        vrr(i) = (ux*ax +uy*ay + uz*az)/arr(i)
        vro(i) = (ux*ax +uy*ay)/aro(i)
        vtt(i) = (uz*aro(i) - vro(i)*az)/arr(i)
      vfi(i) = (uy*ax - ux*ay)/aro(i)
      vvv    = (ux**2 + uy**2 + uz**2)
        vzz(i) = uz       

         qx  = bx(i)
         qy  = by(i)
         qz  = bz(i)
      bbb    = sqrt(qx**2 + qy**2 + qz**2)
        brr(i) = (qx*ax +qy*ay + qz*az)/arr(i)
        bro(i) = (qx*ax +qy*ay)/aro(i)
        btt(i) = (qz*aro(i) - bro(i)*az)/arr(i)
      bfi(i) = (qy*ax - qx*ay)/aro(i)
        bzz(i) = qz       

        add(i) = dens(i)
        app(i) = pres(i)
30    continue

      write(iStdout,*) prefix,
     1       'Coupling time in hour         : ',tCoupleHr
      write(iStdout,*) prefix,
     1        'number of grids   : ',iMax
      write(iStdout,*) prefix,
     1        '1-st r-ro-z-      : ',arr(1),aro(1),azz(1)
      write(iStdout,*) prefix,
     1        '1-st br-bt-bfi    : ',brr(1),btt(1),bfi(1)
      write(iStdout,*) prefix,
     1        'last r-ro-z-roa   : ',arr(iMax),aro(iMax),azz(iMax)
      write(iStdout,*) prefix,
     1        'last br-bt-bfi    : ',brr(iMax),btt(iMax),bfi(iMax)
      write(iStdout,*) prefix      

c ==========================================  identify shock

      call SP_get_ishock
      ish=iShock
      call SP_get_rshock
      rsh=rshck1

      if(DoWriteAll)write(iStdout,*) prefix,
     1      'shock identified ish : ',ish,rsh
      if(DoWriteAll)write(iStdout,*) prefix,
     1       '     at distance r  : ',arr(ish),rsh
      if(DoWriteAll)write(iStdout,*) prefix      

c ====================================  identify end-points

       rmin1 = rmin*(rsh/rshock) 
       imin1 = 0
       rmax2 = rmax*(rsh/rshock)
       imax1 = 0
      do 50 i=1,iMax
      if (arr(i).lt.rmin1)  imin1 = i
      if (arr(i).lt.rmax2)  imax1 = i
50    continue

c      imin1,imax1 largest smaller (or 0)
      
c ================================  Now comes the ansatz
     
      tblast = tCoupleHr             
        nn = nr
c ---------------  compute dlnt=lnt/nn  -- needed to move across

      mp = mm
      ni = mp*(nn+1)
c ---------------------------------------------   inside :

      eta1  = rmin/rshock
      jj = mp
      if (imin1.eq.0) then
      if(DoWriteAll)write(iStdout,*) prefix,
     1        'imin1 = 0 -- not yet prepared ' 
      if(DoWriteAll)write(iStdout,*) prefix,
     1        'rshock/rmin-s desired : ',rmin,rshock 
      if(DoWriteAll)write(iStdout,*) prefix,
     1        'rshock/rmin-s read    : ',arr(1),rsh 
      call CON_stop(' ')
      endif
      i1 = imin1
      ii = i1+1
    
      eta(jj) = eta1
      exx(jj) = eta1*aro(ii)/arr(ii)
      ezz(jj) = eta1*azz(ii)/arr(ii)
     
      evr(jj) = vrr(ii)/rmin1*tblast
      evx(jj) = vro(ii)/rmin1*tblast
      evz(jj) = vzz(ii)/rmin1*tblast

      ebx(jj) = bro(ii)*arr(ii)**2
      ebz(jj) = bzz(ii)*arr(ii)**2
      ebr(jj) = brr(ii)*arr(ii)**2
      ebm(jj) = sqrt(ebx(jj)**2+ebz(jj)**2)

      edd(jj) = add(ii)*arr(ii)**2

                                           !!!  integral 1  !!!
       vv = vrr(ii)*tblast/slamb
      ss1 = alog((arr(ii)-vv)/(rmin1-vv))/slamb 

      if(DoWriteAll)write(iStdout,*) prefix,
     1     'AT R-MIN :        ',rmin
      if(DoWriteAll)write(iStdout,*) prefix,
     1     'Scaled to RMIN1   ',rmin1
      if(DoWriteAll)write(iStdout,*) prefix,
     1      'Found after       ',imin1,arr(imin1)
      if(DoWriteAll)write(iStdout,*) prefix,
     1       'j-eta-exx-ezz:    ',jj,eta(jj),exx(jj),ezz(jj)
      if(DoWriteAll)write(iStdout,*) prefix,
     1      'j-evr-evx-evz:    ',jj,evr(jj),evx(jj),evz(jj)
      if(DoWriteAll)write(iStdout,*) prefix,
     1       'j-ebx-ebz    :    ',jj,ebx(jj),ebz(jj)
      if(DoWriteAll)write(iStdout,*) prefix,
     1       'j-ebr-ebm    :    ',jj,ebr(jj),ebm(jj)
      if(DoWriteAll)write(iStdout,*) prefix,
     1        'j-edd        :    ',jj,edd(jj)
      if(DoWriteAll)write(iStdout,*) prefix      

c -----------------------------------------------  outside :

      eta2  = rmax/rshock
      eta(nn) = eta2
      jj = ni

      ii = imax1

      eta(jj) = eta2
      exx(jj) = eta2*aro(ii)/arr(ii)
      ezz(jj) = eta2*azz(ii)/arr(ii)
     
      evr(jj) = vrr(ii)/arr(ii)*tblst
      evx(jj) = vro(ii)/arr(ii)*tblst
      evz(jj) = vzz(ii)/arr(ii)*tblst

      ebx(jj) = bro(ii)*arr(ii)**2
      ebz(jj) = bzz(ii)*arr(ii)**2
      ebr(jj) = brr(ii)*arr(ii)**2
      ebm(jj) = sqrt(ebx(jj)**2+ebz(jj)**2)

      edd(jj) = add(ii)*arr(ii)**2
      ebfi(jj)= bfi(ii)*arr(ii)

                                           !!!  integral 2  !!!
       vv = vrr(ii)*tblast/slamb
      ss2 = alog((rmax2-vv)/(arr(ii)-vv))/slamb 

      if(DoWriteAll)write(iStdout,*) prefix,
     1       'AT R-MAX :        ',rmax
      if(DoWriteAll)write(iStdout,*) prefix,
     1      'Scaled to RMIN1   ',rmax2
      if(DoWriteAll)write(iStdout,*) prefix,
     1      'Found after       ',imax1,arr(imax1)
      if(DoWriteAll)write(iStdout,*) prefix,
     1      'j-eta-exx-ezz:    ',jj,eta(jj),exx(jj),ezz(jj)
      if(DoWriteAll)write(iStdout,*) prefix,
     1      'j-evr-evx-evz:    ',jj,evr(jj),evx(jj),evz(jj)
      if(DoWriteAll)write(iStdout,*) prefix,
     1      'j-ebx-ebz    :    ',jj,ebx(jj),ebz(jj)
      if(DoWriteAll)write(iStdout,*) prefix,
     1      'j-ebr-ebm    :    ',jj,ebr(jj),ebm(jj)
      if(DoWriteAll)write(iStdout,*) prefix,
     1      'j-edd        :    ',jj,edd(jj)
      if(DoWriteAll)write(iStdout,*) prefix      

c ---------------------------------------------  DIVIDE  
      sum = ss2
      tint(imax1) = ss2

      if(DoWriteAll)write(iStdout,*) prefix,
     1      'imax1=',imax1

CCC      do 110 ii=imax1-1,1,-1
      do 110 ii=imax1-1,2,-1
      vv = vrr(ii)*tblast/slamb
      ss = alog((arr(ii)-vv)/(arr(ii-1)-vv))/slamb 
      sum = sum + ss
      tint(ii) = sum
110   continue
      sum = tint(imin1+1) + ss1
      dlnt = sum/float(nn*mp)
       ede = edd(ni)/ebm(ni)*rmax/rshock*(slamb-evr(ni))*dlnt
       ede = float(mp)*ede
       edsm(ni) = ede*ebm(ni)/edd(ni)

      if(DoWriteAll)write(iStdout,*) prefix,
     1      'Selected s-lambda  : ',sum
      if(DoWriteAll)write(iStdout,*) prefix,
     1      'Integral r1-r2done : ',sum
      if(DoWriteAll)write(iStdout,*) prefix,
     1      'Dlnt - resulting   : ',dlnt
      if(DoWriteAll)write(iStdout,*) prefix,
     1      'ede                : ',ede       
      if(DoWriteAll)write(iStdout,*) prefix,
     1      'edsm(ni)           : ',edsm(ni)
      if(DoWriteAll)write(iStdout,*) prefix,
     1      '=============================================='
CCC      write(iStdout,*) prefix, 'PRESS ANY NUMBER'
CCC      read(*,*)  lull
      write(iStdout,*) prefix 

      jj = 0
      ii = imax1
      tio= 0.
       r2= rmax2
120   jj = jj + 1
       j = ni-jj
      tji = dlnt*float(jj)
121   continue
      tiu = tint(ii)
      if (tiu.lt.tji) then
        tio = tiu
           r2 = arr(ii)
          ii = ii-1
        go to 121
        endif

      if (j.ne.mp) then
           r1  = arr(ii) 
          ti  = (tiu-tio)*slamb
          tj  = (tji-tio)*slamb
          rj  = r2 - (r2-r1)*(1.-exp(-tj))/(1.-exp(-ti))
         eta(j) = rj/rsh 

       exx(j) = eta(j)*aro(ii)/arr(ii)
       ezz(j) = eta(j)*azz(ii)/arr(ii)

         evr(j) = vrr(ii)/rj*tblast
         evx(j) = vro(ii)/rj*tblast
         evz(j) = vzz(ii)/rj*tblast

         ebx(j) = bro(ii)*arr(ii)**2
         ebz(j) = bzz(ii)*arr(ii)**2
         ebr(j) = brr(ii)*arr(ii)**2
         ebm(j) = sqrt(ebx(j)**2+ebz(j)**2)
      
       edd(j) = add(ii)*arr(ii)**2
      endif
         edsm(j)= ede*ebm(j)/edd(j)*exp(-float(jj)*slamb*dlnt)

      if (j.gt.0) go to 120

                                     !!! Normalize - Adjust N,Bfi
      eddmx = edd(ni)
      ebfimx = ebfi(ni)
      sinmx  = exx(ni)/eta(ni)
      omega  = ebfi(ni)/ebr(ni)/sinmx*vrr(imax1)
      
      if(DoWriteAll)write(iStdout,*) prefix,
     1      'eddmx  : ',eddmx
      if(DoWriteAll)write(iStdout,*) prefix,
     1      'ebfimx : ',ebfimx
      if(DoWriteAll)write(iStdout,*) prefix,
     1      'sinmx  : ',sinmx

      edd(ni) = 1.

      do 130 j=ni-1,0,-1
      jj = ni-j
        edd(j) = edd(j)/eddmx
       ebfi(j) = ebfimx/sinmx*edd(j)*exx(j)/eta(j)
        timj = tblast*exp(float(j)*dlnt)
        efi(j) = omega*timj
130   continue

c ------------------------------------  now differences:
      do 210 jj=0,ni
       j1 = jj-mp
       j2 = jj+mp
      if (j1.lt.0) then  
       dex = exx(j2)-exx(jj)
       dez = ezz(j2)-ezz(jj)
       dbm = alog(ebm(j2)/ebm(jj))
       dbfi= alog(ebfi(j2)/ebfi(jj))
       go to 201
       endif
      if (j2.gt.ni) then  
       dex = exx(jj)-exx(j1)
       dez = ezz(jj)-ezz(j1)
       dbm = alog(ebm(jj)/ebm(j1))
       dbfi= alog(ebfi(jj)/ebfi(j1))
      else     
       dex = (exx(j2)-exx(j1))/2.
       dez = (ezz(j2)-ezz(j1))/2.
       dbm = alog(ebm(j2)/ebm(j1))/2.
       dbfi= alog(ebfi(j2)/ebfi(j1))/2.
      endif

201   continue
ccc   edsm(jj) = sqrt(dex**2 + dez**2)
        dbmds(jj) = dbm/edsm(jj) 
       dbfids(jj) = dbfi/edsm(jj) 
      
210   continue

c --------------------------------------  PLAY w fi - dfi

      efi(ni) = 0.
      do 310 jj=0,mp-1
        j = ni-jj
         efi(j) = 0.
310   continue
      do 320 j=ni,mp,-1
       j1= j-mp
       edfi(j) = ebfi(j)/ebm(j)*edsm(j) 
       efi(j1) = efi(j) + edfi(j)
320   continue
     
c --------------------------------- finally rescale TBLAST

      tblast = tblast*(rshock/rsh)**(1./slamb)

c ???????????????????????????????????????????????????

      call CON_io_unit_new(io)
      open(io,file='./SP/ujcme-eta',status='unknown')
      write(io,*) lineno,mm,ni,slamb,thr,tblast
      do 900 j=0,ni 
      write(io,911) j,eta(j),exx(j),ezz(j)
      write(io,912)   evr(j),evx(j),evz(j)
      write(io,912)   ebr(j),ebx(j),ebz(j)
      write(io,912)   ebm(j),ebfi(j),edd(j)
      write(io,913)   edsm(j)
900   continue
911   format(i8,3f12.6)
912   format(8x,3f12.6)
913   format(8x, f12.6)
      close(io)

      do 400 k=1,mp
      arc = 0.
      do j=k,ni,mp
         arc = arc + edsm(j)
      enddo
      if(DoWriteAll)write(iStdout,*) prefix,
     1      'k - meridian arclength in eta: ',k,arc
400   continue
c ???????????????????????????????????????????????????
    
      return
      end

c **************************  end of SP_simile***************


