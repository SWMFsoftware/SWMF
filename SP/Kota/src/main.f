      program sep   

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

      common /size  / nr,nmu,nw, dim1        
      common /suly  / wghtl,wghtmu,wghtw
      common /times/  time,tmax,dlnt0,dt1,dta,kfriss,kacc
      common /blast/  slamb,tblast,tblst1,rshck1,dlnt
      common /radio / nn,rmin,rshock,rmax,r(0:1000)
      common /scphys/ wind,omega,xscatt1
      common /gazdi/  ggamma,bbrad,vvmin,vvmax,ddmin,ddmax,
     1                ccmin,ccmax,aamin,aamax,bbmin,bbmax
c ----------------------------------------------------------
      jnext = 2
       write(*,*) 'first j-stop at ',jnext
      pi = 2.*asin(1.)
      call admit
      call cool   
      call pangle  
c ----------------------------------------------------------
c     initial conditions
      jstep = 0
      write(*,*) 'want to start from negative j-step : ',jstep
      jsep = 1
      write(*,*) 'Particles be calculated from j-sep : ',jsep
      mm = kfriss
      call simile(mm)
      call initial

      kstep = 0 
      istep = jstep*kfriss+kstep
      tblst1= tblast*exp(float(istep)*dlnt)
      rshck1= rshock*(tblst1/tblast)**slamb
       time = tblst1
        dt1 = tblst1*dlnt
      dta = dt1/float(kacc)

c ----------------------------------------------------------
      call peeloff                        
      call opentime
c ----------------------------------------------------------
      ! Time step loop
11    continue

      kstep = 0
      ! Grid refresh loop
21    kstep = kstep+1
      istep = jstep*kfriss+kstep
      tblst1= tblast*exp(float(istep)*dlnt)
      rshck1= rshock*(tblst1/tblast)**slamb
       time = tblst1
        dt1 = tblst1*dlnt
      dta = dt1/float(kacc)
      write(*,900) jstep,kstep,istep,time
900   format(5x,'step:',3i10,5x,'t[hour]: ',f12.6)
c ------------------------------------------ do dynamical step here
      call helios(time,dt1,jstep,kstep)
ccc   if (mod(jstep,200).eq.0.and.kstep.eq.1) read(*,*) lull

      if (jstep.lt.jsep) go to 30

      ! Loop for acceleration of particles on a static grid
      do 31 kk=1,kacc

      call source(dta)
      call deltal(dta)
      call lcycle(dta,wghtl)
      call update

      call subsub(dta)
31    continue
30    continue
c ------------------------------------------    dynamical step done

      if (kstep.lt.kfriss) go to 21
         jstep = jstep+1
      if (jstep.gt.0) call timevar(jstep,time)
      write(*,*) jstep,'  j-step done -- time: ',time
      if (jstep.eq.jnext) then
          if (time.le.1.) io = 11
          if (time.gt.1..and.time.le.2.)   io=12
          if (time.gt.2..and.time.le.6.)   io=13
          if (time.gt.6..and.time.le.12.)  io=14
          if (time.gt.12..and.time.le.48.) io=15
          if (time.gt.48.) io=16
         call alla(io,jstep,time)
       call csilla(io,istep,time)
       write(*,*) 'next j-stop ??'
c??      read(*,*)  jnext
       jnext = jnext+24
         write(*,*) jnext,'   auto - selected'
      end if

       call refresh
       if (time.lt.tmax) go to 11

      if (time.lt.tmax) go to 11
      ! End of main time step loop

      ! Finalize
      call closetime

      stop 'itt a vege'
      end

c ********************* end routine MASTER  *********************
