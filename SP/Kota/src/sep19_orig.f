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
!      common /scphys/ wind,omega,xscatt1 !'wind' is undefined
      common/scphys/ omega,xscatt1
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
11    continue
      kstep = 0
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
      call closetime
      stop 'itt a vege'
      end

c ********************* end routine MASTER  *********************
      
      subroutine admit
      common /size  / nr,nmu,nw, dim1
      common /suly /  wghtl,wghtmu,wghtw
      common /times/  time,tmax,dlnt0,dlnt1,dta,kfriss,kacc
      common /blast/  slamb,tblast,tblst1,rshck1,dlnt
      common /radio / nn,rmin,rshock,rmax,r(0:1000)
      common /gazdi/  ggamma,bbrad,vvmin,vvmax,ddmin,ddmax,
     1                ccmin,ccmax,aamin,aamax,bbmin,bbmax
!      common /inphys/ wind0,period0,xlambda0
!      commented out: is not a common block, variables are 
!      undefined, their use is not proper 
!      I.Sokolov<igorsok@umich.edu>
      common /partid/ iz,massa,ekpp,xlmbda0
!      common /scphys/ wind,omega,xscatt1 !'wind' is undefined
      common/scphys/ omega,xscatt1
      common /impuls/ pmin,pmax,ppin,dlnp,pp(0:600)
      common /energy/ emin,emax,eein,ee(0:600)
      common /speed / wmin,wmax,wwin,ww(0:600)
      common /quelle/ kinj,einj,pinj,qqp(0:600),qqx(0:1000)
      common /scatti/ qex,cmu(40),scmu(40),wsc(0:600),xsc(0:1000)
      common /obsrad/ krmax,krobs(5),robs(5)
      common /obserg/ kemax,keobs(5),eobs(5)
      common /convrt/ aukm,hour,valf
      data  pmass,clight,xkm  / 938., 3.e10, 1.e5 /
c ----------------------------------- scales:
      pi = 2.*asin(1.)
      aukm = 1.5e8
      hour = 3600.
      valf = 21.8
	e0 = pmass
c ------------------------------------ input data:
      in = 3
      open(in,file = 'violet.in')
      read(in,*) nr,nmu,nw, wghtl,wghtmu,wghtw
      write(*,*) nr,nmu,nw, wghtl,wghtmu,wghtw
      read(in,*) rmin,rmax,rshock
      write(*,*) rmin,rmax,rshock
      read(in,*) emin,emax,einj
      write(*,*) emin,emax,einj
      read(in,*) tmax,slamb,kfriss,kacc
      write(*,*) tmax,slamb,kfriss,kacc
      read(in,*) swind,fwind,bbrad,period
      write(*,*) swind,fwind,bbrad,period
      read(in,*) iz,massa,ekpp,xlmbda0,qex
      write(*,*) iz,massa,ekpp,xlmbda0,qex
      read(in,*) krmax,kemax
      Write(*,*) krmax,kemax
      read(in,*) (robs(k),k=1,krmax)
      write(*,*) (robs(k),k=1,krmax)
      read(in,*) (keobs(k),k=1,kemax)
      write(*,*) (keobs(k),k=1,kemax)
      write(*,*) 'beolvasas megvolt '
c --- close(in)
      dim = 3
      dim1= dim-1
c ------------------------------------ gasdynamics
      nn = nr
      mmdim = 2 
      vvmin = swind
      vvmax = fwind
      vcmin =  0.
      vcmax =  0.
      bbrnt =  bbrad
      bbrad =  bbrnt*valf*hour/aukm
      vamin =  0.
      vamax = vamin*sqrt(vvmin/vvmax)
      ddmin = 10. 
      ddmax = ddmin*vvmin/vvmax
      vvmin = vvmin*hour/aukm
      vvmax = vvmax*hour/aukm
      ccmin = vcmin*hour/aukm
      ccmax = vcmax*hour/aukm
      aamin = vamin*hour/aukm
      aamax = vamax*hour/aukm
      bbmin = aamin*sqrt(ddmin)
      bbmax = aamax*sqrt(ddmax)
c ------------------------------------ times:     
!     dt1 = dt0/kfriss  !dt0 is undefined at the moment
!     dta = dt1/kacc    !Commented out, I.Sokolov<igorsok@umich.edu>
c ------------------------------------ calculation
!       wind = wind0*hour/aukm!wind0 is undefined at the moment
!      swind = swind0*hour/aukm!swind0 is undefined at the moment
      if (period.gt.1.e3) then
	  omega = 0.
	  else
	  omega = 2.*pi/period/24.
	  end if

      xlambda0 = xlmbda0
      if (xlambda0.gt.1.e3) then
	  xscatt1 = 0.
	  else
	  xscatt1 = ((3-qex)*(1-qex)/3.)/xlambda0
	  end if
      return 
      end

c  *****************  end ADMIT   **************************
      
      subroutine cool   
      common /size  / nr,nmu,nw, dim1
      common /partid/ iz,massa,ekpp,xlmbda0
      common /impuls/ pmin,pmax,ppin,dlnp,pp(0:600)
      common /energy/ emin,emax,eein,ee(0:600)
      common /speed / wmin,wmax,wwin,ww(0:600)
      common /scatti/ qex,cmu(40),scmu(40),wsc(0:600),xsc(0:1000)
      common /quelle/ kinj,einj,pinj,qqp(0:600),qqx(0:1000)
      data  pmass,clight,xkm  / 938., 3.e10, 1.e5 /
c ----------------------------------- scales:
      pi = 2.*asin(1.)
      aukm = 1.5e8
      hour = 3600.
      valf = 21.8
	e0 = pmass
c ----------------------------------- scales:

	gv = 1000.
        e0 = pmass
      amass= float(massa)*pmass
	qz = float(iz)
      expo =-ekpp

      pmin = sqrt(emin*(emin+2.*e0))
      pmax = sqrt(emax*(emax+2.*e0))
      pinj = sqrt(einj*(einj+2.*e0))
      dlnp = alog(pmax/pmin)/float(nw)

      eein = einj
      ppin = pinj
      btin = ppin/(eein+e0)
      wwin = btin*clight*hour/xkm/aukm

      do 11 k = 0,nw
      pp(k) = pmin*exp(float(k)*dlnp)
       etot = sqrt(pp(k)**2+e0**2)
       beta = pp(k)/etot
      ee(k) = pp(k)**2/(etot+e0)
      ww(k) = beta*clight*hour/xkm/aukm
      qqp(k)= 0.
       rg   = float(massa)*pp(k)/abs(qz)/gv
      wsc(k)= ww(k)*rg**expo

      qqp(k)= exp(-ee(k)/eein)/ppin**3

11    continue

      qk  = alog(pinj/pmin)/dlnp 
      if (qk.ge.0) then
         kk  = qk
         kinj= kk
c        qqp(kk) = 1./pinj**2/dlnp
	 endif

      return
      end

c  *****************  end subroutine COOL     ******************

      subroutine peeloff
      common /size  / nr,nmu,nw, dim1
      common /impuls/ pmin,pmax,ppin,dlnp,pp(0:600)
      common /energy/ emin,emax,eein,ee(0:600)
      common /speed/  wmin,wmax,wwin,ww(0:600)
      common /scatti/ qex,cmu(40),scmu(40),wsc(0:600),xsc(0:1000)
      common /quelle/ kinj,einj,winj,qqw(0:600),qqx(0:1000)
      common /peel /  pex,fpeel(0:600),gpeel(0:600)
      common /solutn/ f(0:1000,0:40,0:600),df(0:1000,0:40,0:600)
     
      pex = 0.
      write(*,*) 'PEEL-OFF -- what exponent ??',pex
      do 30 k = 0,nw
      fact = (pp(k)/ppin)**pex
      fpeel(k) = 2./(1.+fact)
      gpeel(k) = pex*fact/(1.+fact)
      write(*,911) pex,k,ee(k),fact,gpeel(k),fpeel(k)
911   format(f10.2,i5,3f12.6,e14.4)

      qqw(k) = qqw(k)/fpeel(k)
      do 10 i=0,nr
      do 20 j=0,nmu
      f(i,j,k) = f(i,j,k)/fpeel(k)
20    continue
10    continue
30    continue

      return
      end

c  *****************  end subroutine PEELOFF  ******************

      subroutine pangle  
      common /size  / nr,nmu,nw, dim1
      common /pitch / mm,amu(0:40),sint(0:40),dmu
      common /scatti/ qex,cmu(40),scmu(40),wsc(0:600),xsc(0:1000)
c     linear grid in mu=cost
      mm = nmu
       m = mm/2
      dmu = 2./float(mm) 
      do 11 jj=0,m
      jm = mm-jj
      amu(jj) = 1. - float(jj)*dmu
      amu(jm) = - amu(jj)
      sint(jj) = sqrt(1.-amu(jj)**2)
      sint(jm) = sint(jj)
11    continue
       amu(0) = 1.
      sint(0) = 0.
       amu(m) = 0.
      sint(m) = 1.
      amu(mm) = -1
      sint(mm) = 0.
c  
      do 21 jj=1,m
      jm = mm-jj+1
        ccmu = (amu(jj)+amu(jj-1))/2.
        ssmu = 1.-ccmu**2
       absmu = abs(ccmu)
       scatty = ssmu*absmu**qex*3./(1.-qex)/(3.-qex)
       cmu(jj) =  ccmu
       scmu(jj) = scatty/dmu**2/2.   
       cmu(jm) = -ccmu
       scmu(jm) = scatty/dmu**2/2.
21    continue
c ---  scatty representing  (1/lambda)
      return
      end

c  *****************  end subroutine PITCH    ******************

      subroutine initial                                       
      common /size /  nr,nmu,nw, dim1
      common /radio / nn,rmin,rshock,rmax,r(0:1000)
      common /blast/  slamb,tblast,tblst1,rshck1,dlnt
      common /solutn/ f(0:1000,0:40,0:600),df(0:1000,0:40,0:600)
      common /plasma/ algbb(0:1000),algll(0:1000),algnn(0:1000),
     1                vr(0:1000)

      tblst1 = tblast
      rshck1 = rshock
       time  = tblst1

      do 10 i=0,nr
	 algbb(i) = 0.       
	 algll(i) = 0.      
	 algnn(i) = 0.       
	    vr(i) = 0.       

      do 30 k=0,nw
      do 20 j=0,nmu
       f(i,j,k) = 0.
      df(i,j,k) = 0.
20    continue
30    continue
10     continue
      return
      end

c **************************  end of INiT ****************

c ============================ real calculation cycle starts here ==

      subroutine refresh                     
      common /size  / nr,nmu,nw, dim1
      common /radio / nn,rmin,rshock,rmax,r(0:1000)
      common /times/  time,tmax,dt0,dt1,dta,kfriss,kacc
      common /blast/  slamb,tblast,tblst1,rshck1,dlnt
!      common /scphys/ wind,omega,xscatt1 !'wind' is undefined
      common/scphys/ omega,xscatt1
      common /elem /  zr(0:1000),zv(0:1000),zp(0:1000),zn(0:1000)
      common /magia/  qb(0:1000),zb(0:1000),zt(0:1000) 
      common /plasma/ algbb(0:1000),algll(0:1000),algnn(0:1000),
     1                vr(0:1000)
      common /solutn/ f(0:1000,0:40,0:600),df(0:1000,0:40,0:600)
      common /smile/  mp,ni
     1       ,eta(0:6000),exx(0:6000),ezz(0:6000),efi(0:6000)
     2       ,evr(0:6000),evx(0:6000),evz(0:6000)
     3       ,edsm(0:6000),edfi(0:6000)
     4       ,ebm(0:6000),ebr(0:6000),ebt(0:6000),ebfi(0:6000)
     5       ,ebx(0:6000),ebz(0:6000),edd(0:6000)
     6       ,dbmds(0:6000),dbfids(0:6000)

c **********************************  concept : *****************
c     shifts back to regions to be considered                   *
c ********************************** ---------- *****************
      
      do 11 i=0,nr-1
      i2=i+1

      algbb(i) = algbb(i2)
      algnn(i) = algnn(i2)
      algll(i) = algll(i2)
	 vr(i) = vr(i2)
      do 21 j=0,nmu
      do 31 k=0,nw
	 f(i,j,k) = f(i2,j,k)
31    continue
21    continue
11    continue
     
      i = nr
	  r(nr) = rshck1*eta(ni)
         vr(nr) = r(nr)/tblst1*evr(ni)
        
          bbr   = ebr(ni)/r(i)**2
          bbm   = ebm(ni)/r(i)**2
          bbfi  = ebfi(ni)/r(i)
          bbb   = sqrt(bbm**2+bbfi**2)
          abb   = alog(bbb)

	  ddd   = edd(ni)/r(i)**2
	  ann   = alog(ann)
       
	dll   = rshck1*edsm(ni)*bbb/bbm
	all   = alog(dll)

      algbb(nr) = abb
      algnn(nr) = ann
      algll(nr) = all

      do 121 j=0,nmu
      do 131 k=0,nw
	 f(i,j,k) = 0.
131   continue
121   continue

      return 
      end

c ********************  end routine REFRESH *********************
                                                  !!!! REWORK
      subroutine helios(t,dt,jj,kk)
      common /size / nr,nmu,nw, dim1
      common /gazdi/ ggamma,bbrad,vvmin,vvmax,ddmin,ddmax,
     1               ccmin,ccmax,aamin,aamax,bbmin,bbmax
      common /elem / zr(0:1000),zv(0:1000),zp(0:1000),zn(0:1000)
      common /magia/ qb(0:1000),qd(0:1000),zb(0:1000) 
      
      common /blast/  slamb,tblast,tblst1,rshck1,dlnt
      common /radio / nn,rmin,rshock,rmax,r(0:1000)
!      common /scphys/ wind,omega,xscatt1 !'wind' is undefined
      common/scphys/ omega,xscatt1
      common /spiral/ tll(0:1000), dbdl(0:1000),vl(0:1000)
      common /coeff / dvdt(0:1000),dldt(0:1000),dbdt(0:1000),
     1                dndt(0:1000)
      common /quelle/ kinj,einj,winj,qqw(0:600),qqx(0:1000)
      common /scatti/ qex,cmu(40),scmu(40),wsc(0:600),xsc(0:1000)
      common /plasma/ algbb(0:1000),algll(0:1000),algnn(0:1000),
     1                vr(0:1000)
      common /smile / mp,ni
     1       ,eta(0:6000),exx(0:6000),ezz(0:6000),efi(0:6000)
     2       ,evr(0:6000),evx(0:6000),evz(0:6000)
     3       ,edsm(0:6000),edfi(0:6000)
     4       ,ebm(0:6000),ebr(0:6000),ebt(0:6000),ebfi(0:6000)
     5       ,ebx(0:6000),ebz(0:6000),edd(0:6000)
     6       ,dbmds(0:6000),dbfids(0:6000)

c +++++++++++++++++++++++++++++++++++++++     concept:   +++++ 
c    Transfers between similarity solution and actual         +
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

      pi = 2.*asin(1.)

c -----------------------------------------------------
      rshck1 = rshock*(tblst1/tblast)**slamb 
      
      do 10 i=0,nr
      ii = mp*(i+1)-kk 

c -------------------------------   this one new for transport :

         r(i) = rshck1*eta(ii)
        vvr   = evr(ii)*r(i)/tblst1
      
        bbr   = ebr(ii)/r(i)**2
	bbm   = ebm(ii)/r(i)**2
	bbfi  = ebfi(ii)/r(i)
	bbb   = sqrt(bbm**2+bbfi**2)
        abb   = alog(bbb)

	dll   = rshck1*edsm(ii)*bbb/bbm
	all   = alog(dll)
       tll(i) = dll                           

	ddd   = edd(ii)/r(i)**2
	ann   = alog(ddd)

      dbdt(i) = (abb-algbb(i))/dt
      dndt(i) = (ann-algnn(i))/dt
      dldt(i) = (all-algll(i))/dt
      dvdt(i) = (vvr-vr(i))/dt
      
        vl(i) = 0.                            
c ??????????????????????????????????????????????????????
      if (i.lt.0) then
      write(*,*) 'I AM HERE IN HELIOS - kk,i,ii=  ', kk,i,ii
      write(*,*)   'R:       ',r(i)
      write(*,*)   'EBM-Fi : ',ebm(ii),ebfi(ii)
      write(*,*)   'BBM-Fi : ',bbm,bbfi,bbb
      write(*,*)   'BB     : ',abb,algbb(i),dbdt(i)
      write(*,*)   'LL     : ',all,algll(i),dldt(i)
      write(*,*)   'NN     : ',ann,algnn(i),dndt(i)
      read(*,*)    lull
      endif
c ??????????????????????????????????????????????????????

      algbb(i) = abb
      algnn(i) = ann
      algll(i) = all
         vr(i) = vvr                  

	qqx(i) = 0.        
         bbfac = 1./r(i)
	 bbfak = bbfac
	xsc(i) = xscatt1*bbfak
	if (eta(ii).lt.1.15) then
	   xsc(i) = 50.*xsc(i)
	   qqx(i) = 1/r(i)**2
	   endif
	if (eta(ii).lt.0.50) then
	   qqx(i) = 0.        
	   endif

10    continue
      
      do 20 i=0,nr
      if (i.eq.0) then
         dbdl(i) = (algbb(i+1)-algbb(i))/tll(i)
	 go to 20
      endif
      if (i.eq.nr) then
         dbdl(i) = (algbb(i)-algbb(i-1))/tll(i)
	 go to 20
      endif
         dbdl(i) = (algbb(i+1)-algbb(i-1))/tll(i)/2.
20    continue

c +++++++++++++++++++++++++++++++++++++++++++
      dbmax = 0.
      dbmin = 0.
      dlmax = 0.
      dlmin = 0.
      do 88 i=0,nr
      if (dbdt(i).gt.dbmax) then
	  dbmax = dbdt(i)
	  ibmax = i
	  endif
      if (dbdt(i).lt.dbmin) then
	  dbmin = dbdt(i)
	  ibmin = i
	  endif
      if (dldt(i).gt.dlmax) then
	  dlmax = dldt(i)
	  ilmax = i
	  endif
      if (dldt(i).lt.dlmin) then 
	  dlmin = dldt(i)
	  ilmin = i
	  endif
88    continue
      write(*,*) 'Min dldt: ',ilmin,dlmin
      write(*,*) 'Max dldt: ',ilmax,dlmax
      write(*,*) 'Min dbdt: ',ibmin,dbmin
      write(*,*) 'Max dbdt: ',ibmax,dbmax

      return
      end

c ********************  end routine GASDYN  *********************






      subroutine simile(mm)                 !!!  REWORK !!!
      
      common /size /  nr,nmu,nw, dim1
      common /radio/  nn,rmin,rshock,rmax,r(0:1000)
      common /blast/  slamb,tblast,tblst1,rshck1,dlnt
      common /smile/  mp,ni
     1       ,eta(0:6000),exx(0:6000),ezz(0:6000),efi(0:6000)
     2       ,evr(0:6000),evx(0:6000),evz(0:6000)
     3       ,edsm(0:6000),edfi(0:6000)
     4       ,ebm(0:6000),ebr(0:6000),ebt(0:6000),ebfi(0:6000)
     5       ,ebx(0:6000),ebz(0:6000),edd(0:6000)
     6       ,dbmds(0:6000),dbfids(0:6000)

      dimension rx(2500,20),ry(2500,20),rz(2500,20)
      dimension vx(2500,20),vy(2500,20),vz(2500,20)
      dimension bx(2500,20),by(2500,20),bz(2500,20)
      dimension dd(2500,20),pp(2500,20)
      dimension axx(2500),ayy(2500),azz(2500)
      dimension tint(2500)

      dimension arr(2500),aro(2500)
      dimension vrr(2500),vtt(2500),vro(2500),vzz(2500),vfi(2500)
      dimension brr(2500),btt(2500),bro(2500),bzz(2500),bfi(2500)
      dimension add(2500),app(2500)
      dimension thour(20),imax(20)

c **************************************  concept: *****************
c     put in a similarity solution shifted so that eta=1 at i=nn   *
c     CONCEPT: rmin-rshock-rmax at start -- read  sim-solution     *
c              then place so that spans the desired region         *
c              Divide into 1001 !! equal dlnt regions              *
c **************************************  -------- *****************

      in = 7
      write(*,*) 'Select line (1-2-3): '
      write(*,*) 'PRESS 1'
      read(*,*)  lineno
      write(*,*) 'lineno=',lineno
      if (lineno.eq.1) open(in,file = 'evolv1.cme' )
      if (lineno.eq.2) open(in,file = 'evolv2.cme' )
      if (lineno.eq.3) open(in,file = 'evolv3.cme' )

      read(in,*)
      read(in,*)
      read(in,*) lineno,kmax  !   kmax = no of snapshots
      read(in,*) 

      do 10 k=1,kmax 
      read(in,*)  kx,imx,tx
      thour(k) = tx
       imax(k) = imx

      do 20 i=1,imx 
         read(in,*)  i1,rx(i,k),ry(i,k),rz(i,k)
         read(in,*)  i2,vx(i,k),vy(i,k),vz(i,k)
         read(in,*)  i3,bx(i,k),by(i,k),bz(i,k)
         read(in,*)  i4,dd(i,k),pp(i,k)
20    continue

      read(in,*) 
10    continue
      close(in)

      write(*,*) 'beolvasas megvolt:'
      write(*,*) 'lineno,kmax :  ',lineno,kmax
      write(*,*)      

c -------------------------------- beolvasas meglenne
      ausun = 200.
      hour  = 3600.
      aukm  = 1.5e8
      gauss = 1.e5
      dens0 = 1.e20
      pres0 = 1.e5
      bb0   = 100.
 
      write(*,*) 'selected line: ',lineno
      write(*,*) 'Which snapshot to take (1-20): '
      write(*,*) 'PRESS 9'
      read(*,*)  kk    

      imx = imax(kk)
      thr = thour(kk)
      do 30 i=1,imx
	   ax  = rx(i,kk)
	   ay  = ry(i,kk)
	   az  = rz(i,kk)
        aro(i) = sqrt(ax**2 + ay**2)
        arr(i) = sqrt(ax**2 + ay**2 + az**2)
        axx(i) = ax
        ayy(i) = ay
        azz(i) = az

	   ux  = vx(i,kk)
	   uy  = vy(i,kk)
	   uz  = vz(i,kk)
        vrr(i) = (ux*ax +uy*ay + uz*az)/arr(i)
        vro(i) = (ux*ax +uy*ay)/aro(i)
        vtt(i) = (uz*aro(i) - vro(i)*az)/arr(i)
	vfi(i) = (uy*ax - ux*ay)/aro(i)
	vvv    = (ux**2 + uy**2 + uz**2)
        vzz(i) = uz       

	   qx  = bx(i,kk)
	   qy  = by(i,kk)
	   qz  = bz(i,kk)
	bbb    = sqrt(qx**2 + qy**2 + qz**2)
        brr(i) = (qx*ax +qy*ay + qz*az)/arr(i)
        bro(i) = (qx*ax +qy*ay)/aro(i)
        btt(i) = (qz*aro(i) - bro(i)*az)/arr(i)
	bfi(i) = (qy*ax - qx*ay)/aro(i)
        bzz(i) = qz       

        add(i) = dd(i,kk)
        app(i) = pp(i,kk)
30    continue

      write(*,*) 'consider snapshot : ',kk
      write(*,*) 't in hour         : ',thr
      write(*,*) 'number of grids   : ',imx
      write(*,*) '1-st r-ro-z-      : ',arr(1),aro(1),azz(1)
      write(*,*) '1-st br-bt-bfi    : ',brr(1),btt(1),bfi(1)
      write(*,*) 'last r-ro-z-roa   : ',arr(imx),aro(imx),azz(imx)
      write(*,*) 'last br-bt-bfi    : ',brr(imx),btt(imx),bfi(imx)
      write(*,*)      

c ==========================================  identify shock

	ish = 0
	dmax= 0.
	ddn = add(1)*arr(1)**2
      do 40 i=2,imx
        drr = arr(i)-arr(i-1) 
	ddo = ddn
	ddn = add(i)*arr(i)**2
        diff = (ddo-ddn)/(ddo+ddn)*arr(i)/drr
        if (diff.gt.dmax) then
	   ish = i-1
	  dmax = diff
	  endif
40    continue
           rsh = arr(ish)

      write(*,*) 'shock azonositva ish : ',ish,rsh
      write(*,*) '     tavolsagban r   : ',arr(ish),rsh
      write(*,*)      

c ====================================  identify end-points

       rmin1 = rmin*(rsh/rshock) 
       imin1 = 0
       rmax2 = rmax*(rsh/rshock)
       imax1 = 0
      do 50 i=1,imx
	if (arr(i).lt.rmin1)  imin1 = i
	if (arr(i).lt.rmax2)  imax1 = i
50    continue

c      imin1,imax1 largest smaller (or 0)
      
c ================================  Now comes the ansatz
     
      tblast = thr              
	  nn = nr
c ---------------  compute dlnt=lnt/nn  -- needed to move across

      mp = mm
      ni = mp*(nn+1)
c ---------------------------------------------   inside :

      eta1  = rmin/rshock
      jj = mp
      if (imin1.eq.0) then
	write(*,*)  'imin1 = 0 -- not yet prepared ' 
	write(*,*)  'rshock/rmin-s desired : ',rmin,rshock 
	write(*,*)  'rshock/rmin-s read    : ',arr(1),rsh 
	stop 
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

      write(*,*) 'AT R-MIN :        ',rmin
      write(*,*) 'Scaled to RMIN1   ',rmin1
      write(*,*) 'Found after       ',imin1,arr(imin1)
      write(*,*) 'j-eta-exx-ezz:    ',jj,eta(jj),exx(jj),ezz(jj)
      write(*,*) 'j-evr-evx-evz:    ',jj,evr(jj),evx(jj),evz(jj)
      write(*,*) 'j-ebx-ebz    :    ',jj,ebx(jj),ebz(jj)
      write(*,*) 'j-ebr-ebm    :    ',jj,ebr(jj),ebm(jj)
      write(*,*) 'j-edd        :    ',jj,edd(jj)
      write(*,*)      

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

      write(*,*) 'AT R-MAX :        ',rmax
      write(*,*) 'Scaled to RMIN1   ',rmax2
      write(*,*) 'Found after       ',imax1,arr(imax1)
      write(*,*) 'j-eta-exx-ezz:    ',jj,eta(jj),exx(jj),ezz(jj)
      write(*,*) 'j-evr-evx-evz:    ',jj,evr(jj),evx(jj),evz(jj)
      write(*,*) 'j-ebx-ebz    :    ',jj,ebx(jj),ebz(jj)
      write(*,*) 'j-ebr-ebm    :    ',jj,ebr(jj),ebm(jj)
      write(*,*) 'j-edd        :    ',jj,edd(jj)
      write(*,*)      

c ---------------------------------------------  DIVIDE  
      sum = ss2
      tint(imax1) = ss2

      do 110 ii=imax1-1,1,-1
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

      write(*,*) 'Selected s-lambda  : ',sum
      write(*,*) 'Integral r1-r2done : ',sum
      write(*,*) 'Dlnt - resulting   : ',dlnt
      write(*,*) 'ede                : ',ede       
      write(*,*) 'edsm(ni)           : ',edsm(ni)
      write(*,*) 'PRESS ANY NUMBER'
      read(*,*)  lull
      write(*,*) 

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
      
      write(*,*) 'eddmx  : ',eddmx
      write(*,*) 'ebfimx : ',ebfimx
      write(*,*) 'sinmx  : ',sinmx

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
      io = 11
      open(io,file='ujcme-eta')
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
      do 401 j=k,ni,mp
401   arc = arc + edsm(j)
      write(*,*) 'k - meridian arclength in eta: ',k,arc
400   continue
c ???????????????????????????????????????????????????
    
      return
      end

c **************************  end of SIMILE***************



      subroutine source(dt)          
      common /size  / nr,nmu,nw, dim1
      common /pitch / mm,amu(0:40),sint(0:40),dmu
      common /impuls/ pmin,pmax,ppin,dlnp,pp(0:600)
      common /solutn/ f(0:1000,0:40,0:600),df(0:1000,0:40,0:600)
      common /quelle/ kinj,einj,pinj,qqp(0:600),qqx(0:1000)
      do 31 k=0,nw
      do 11 i=1,nr
	 qq = qqp(k)*qqx(i)
      do 21 j=0,nmu
	 df(i,j,k) = qq*dt
21    continue
11    continue
31    continue

      return
      end
c  *****************  end subroutine SOURCE   ******************


      subroutine deltal(dt)                      !!! REWORK DLOGT!!!
      common /size  / nr,nmu,nw, dim1
      common /pitch / mm,amu(0:40),sint(0:40),dmu
      common /speed / wmin,wmax,wwin,ww(0:600)
      common /spiral/ tll(0:1000), dbdl(0:1000),vl(0:1000)
      common /coeff / dvdt(0:1000),dldt(0:1000),dbdt(0:1000),
     1                dndt(0:1000)
      common /solutn/ f(0:1000,0:40,0:600),df(0:1000,0:40,0:600)
      data  cau / 7.2/

      m = mm/2
      nr1 = nr-1

      do 30 k=0,nw
      do 20 jj=0,m-1

c -------------------------------- out direction(s):
      j = jj
	 df(0,j,k) = 0.
      vv =  ww(k)*amu(j)
      do 11 i=1,nr
	 ae = 1.+ vv*vl(i)/cau**2
         fact = vv*dt/tll(i)/ae
	 if (i.eq.1) then
	 df(i,j,k) = df(i,j,k) - fact*(f(i,j,k)-f(i-1,j,k)) 
	 else
	 df(i,j,k) = df(i,j,k)  
     1    - fact*(1.5*f(i,j,k)-2.*f(i-1,j,k)+0.5*f(i-2,j,k)) 
	 endif
11    continue
      
c --------------------------------  IN direction(s):
      j = nmu-jj
      vv = ww(k)*amu(j)
	 i=nr  
         fact = vv*dt/tll(i)
         df(i,j,k) = df(i,j,k) + fact*f(i,j,k)
      do 12 i=nr1,1,-1
	 ae = 1.+ vv*vl(i)/cau**2
         fact = vv*dt/tll(i)/ae
	 if (i.eq.nr1) then
c??      df(i,j,k) = df(i,j,k) + fact*(f(i,j,k)-f(i+1,j,k)) 
	 df(i,j,k) = df(i,j,k) 
     1    + fact*(1.5*f(i,j,k)-2.*f(i+1,j,k)) 
	 else
	 df(i,j,k) = df(i,j,k) 
     1    + fact*(1.5*f(i,j,k)-2.*f(i+1,j,k)+0.5*f(i+2,j,k)) 
	 end if
12    continue
c -------------------------------- end IN direction(s)
         df(0,j,k) = 0.
	 
20    continue
30    continue
      return
      end
c  *****************  end subroutine DELTA-L  ******************
                                                 !!! REWORK DLOGT!!!
      subroutine lcycle(dt,wght)
      common /size  / nr,nmu,nw, dim1
      common /pitch / mm,amu(0:40),sint(0:40),dmu
      common /speed / wmin,wmax,wwin,ww(0:600)
      common /spiral/ tll(0:1000), dbdl(0:1000),vl(0:1000)
      common /coeff / dvdt(0:1000),dldt(0:1000),dbdt(0:1000),
     1                dndt(0:1000)
      common /solutn/ f(0:1000,0:40,0:600),df(0:1000,0:40,0:600)
      data  cau / 7.2/

      m = mm/2
      nr1 = nr-1
      do 30 k=0,nw
      do 20 jj=0,m-1

      j = jj
      vv = ww(k)*amu(j)

	df(0,j,k) = 0.
      do 11 i=1,nr
	  ae = 1.+ vv*vl(i)/cau**2
	fact = wght*vv*dt/tll(i)/ae
	if (i.eq.1) then
	df(i,j,k) = (df(i,j,k)+fact*df(i-1,j,k))/(1.+fact)
	else
	df(i,j,k) = (df(i,j,k) 
     1             + fact*(2.*df(i-1,j,k)-.5*df(i-2,j,k)))
     2             /(1.+1.5*fact)
	endif
11    continue

      j = mm-jj
      vv = ww(k)*amu(j)

	i=nr
	  ae = 1.+ vv*vl(i)/cau**2
	fact = wght*vv*dt/tll(i)/ae
        df(nr,j,k) = df(nr,j,k)/(1.-fact)
      do 12 i=nr1,1,-1
	  ae = 1.+ vv*vl(i)/cau**2
	fact = wght*vv*dt/tll(i)/ae
	if (i.eq.nr1) then
c???    df(i,j,k) = (df(i,j,k)-fact*df(i+1,j,k))/(1.-fact)
	df(i,j,k) = (df(i,j,k) - fact*(2.*df(i+1,j,k))) 
     2             /(1.-1.5*fact)
	else
	df(i,j,k) = (df(i,j,k) 
     1             - fact*(2.*df(i+1,j,k)-0.5*df(i+2,j,k)))
     2             /(1.-1.5*fact)
        endif
12    continue
	df(0,j,k) = 0.

20    continue
30    continue

      return
      end
c  *****************  end subroutine L-CYCLE  ******************
                                                                    
      subroutine deltamu(i,dt)
      common /size  / nr,nmu,nw, dim1
      common /pitch / mm,amu(0:40),sint(0:40),dmu
      common /speed / wmin,wmax,wwin,ww(0:600)
      common /spiral/ tll(0:1000), dbdl(0:1000),vl(0:1000)
      common /coeff / dvdt(0:1000),dldt(0:1000),dbdt(0:1000),
     1                dndt(0:1000)
      common /scatti/ qex,cmu(40),scmu(40),wsc(0:600),xsc(0:1000)
      common /solutn/ f(0:1000,0:40,0:600),df(0:1000,0:40,0:600)
      data  cau / 7.2/
      nr1 = nr-1
       mm = nmu
      mm1 = mm-1
ccc   do 11 i=1,nr
      do 31 k=0,nw
      
      fact = xsc(i)*wsc(k)*dt
      df(i,0,k) = df(i,0,k) + 2.*fact*scmu(1)*(f(i,1,k)-f(i,0,k))
      df(i,mm,k)= df(i,mm,k)+ 2.*fact*scmu(mm)*(f(i,mm1,k)-f(i,mm,k))

      do 21 j=1,mm1
      ae = 1.+ ww(k)*amu(j)*vl(i)/cau**2
      dte= dt/ae
       
      foc =-sint(j)**2*(0.5*ww(k)*dbdl(i)+dvdt(i)/ww(k) +
     1                  amu(j)*(dldt(i) + 0.5*dbdt(i)))
      foc = dte/dmu*foc

      if (foc.ge.0.) then
	 if (j.eq.mm1) then
	 df(i,j,k) = df(i,j,k) - 
     1           foc*(f(i,j,k)-f(i,j+1,k)) 
	 else
	 df(i,j,k) = df(i,j,k) - 
     1           foc*(1.5*f(i,j,k)-2.*f(i,j+1,k)+.5*f(i,j+2,k)) 
	 end if
      else
	 if (j.eq.1) then
	 df(i,j,k) = df(i,j,k) + 
     1           foc*(f(i,j,k)-f(i,j-1,k)) 
	 else
	 df(i,j,k) = df(i,j,k) + 
     1           foc*(1.5*f(i,j,k)-2.*f(i,j-1,k)+.5*f(i,j-2,k)) 
	 end if
      end if

c --- add scattering :

	sc1 = fact/ae*scmu(j)
	sc2 = fact/ae*scmu(j+1)
	df(i,j,k) = df(i,j,k) +
     1           sc1*(f(i,j-1,k)-f(i,j,k))+sc2*(f(i,j+1,k)-f(i,j,k))

21    continue
31    continue
11    continue
      return
      end
c  *****************  end subroutine DELTA-MU  ******************
                                                 !!! REWORK DLOGT!!!
      subroutine mucycle(i,dt,wght)
      common /size  / nr,nmu,nw, dim1
      common /pitch / mm,amu(0:40),sint(0:40),dmu
      common /speed / wmin,wmax,wwin,ww(0:600)
      common /scatti/ qex,cmu(40),scmu(40),wsc(0:600),xsc(0:1000)
      common /spiral/ tll(0:1000), dbdl(0:1000),vl(0:1000)
      common /coeff / dvdt(0:1000),dldt(0:1000),dbdt(0:1000),
     1                dndt(0:1000)
      common /solutn/ f(0:1000,0:40,0:600),df(0:1000,0:40,0:600)
      dimension a2(0:1000),a1(0:1000),bb(0:1000),
     1          c1(0:1000),c2(0:1000),xy(0:1000)
      data  cau / 7.2/
      nr1 = nr-1
       mm = nmu
      mm1 = mm-1
ccc   do 11 i=1,nr
      do 31 k=0,nw
	 fact = wght*wsc(k)*xsc(i)*dt
c --  mu = 1  (or j=0)
	 ae = 1.+ww(k)*vl(i)/cau**2
        sc1 = fact*scmu(1)/ae
	a2(0) = 0.
	a1(0) = 0.
	bb(0) = 1. + 2.*sc1
	c1(0) = -2.*sc1
	c2(0) = 0.
	xy(0) = df(i,0,k)
c --  mu =-1  (or j=nmu)
	 ae = 1.-ww(k)*vl(i)/cau**2
        sc1 = fact*scmu(mm)/ae
	a2(nmu) = 0.
	a1(nmu) = -2.*sc1
	bb(nmu) = 1. + 2.*sc1
	c1(nmu) = 0.
	c2(nmu) = 0.
	xy(nmu) = df(i,nmu,k)
      do 21 j=1,mm1
	 ae = 1.+ww(k)*vl(i)*amu(j)/cau**2
c --- focusing :
       foc =-sint(j)**2*(0.5*ww(k)*dbdl(i)+dvdt(i)/ww(k) +
     1                   amu(j)*(dldt(i) + 0.5*dbdt(i) ))
       foc = wght*dt/dmu*foc/ae
      if (foc.ge.0.) then
	 a2(j) = 0.
	 a1(j) = 0.
	 if (j.eq.mm1) then
	 bb(j) = 1. + foc
	 c1(j) = -foc
	 c2(j) =  0.         
	 else
	 bb(j) = 1. + 1.5*foc
	 c1(j) = -2.*foc
	 c2(j) =  0.5*foc
	 end if
      else
	 if (j.eq.1) then
	 a2(j) = 0.
	 a1(j) = foc
	 bb(j) = 1. - foc
	 else
	 a2(j) = -0.5*foc
	 a1(j) =  2.0*foc
	 bb(j) =  1. - 1.5*foc
	 end if
	 c1(j) = 0.
	 c2(j) = 0.
      end if
	 xy(j) = df(i,j,k)
c --- add scattering :
	sc1 = fact*scmu(j)/ae
	sc2 = fact*scmu(j+1)/ae
	a1(j) = a1(j) -sc1 
	bb(j) = bb(j) + sc1 + sc2
	c1(j) = c1(j) -sc2
21    continue
c --- now solve and fill :
      call solve5(mm,a2,a1,bb,c1,c2,xy)
      do 22 j=0,mm 
      df(i,j,k) = xy(j)
22    continue
31    continue
11    continue
      return
      end
c  *****************  end subroutine MU-CYCLE  ******************
                                                 !!! REWORK DLOGT!!!
      subroutine deltap(i,dt)
      common /size  / nr,nmu,nw, dim1
      common /pitch / mm,amu(0:40),sint(0:40),dmu
      common /impuls/ pmin,pmax,ppin,dlnp,pp(0:600)
      common /speed / wmin,wmax,wwin,ww(0:600)
      common /spiral/ tll(0:1000), dbdl(0:1000),vl(0:1000)
      common /coeff / dvdt(0:1000),dldt(0:1000),dbdt(0:1000),
     1                dndt(0:1000)
      common /solutn/ f(0:1000,0:40,0:600),df(0:1000,0:40,0:600)
      common /peel /  pex,fpeel(0:600),gpeel(0:600)
      data  cau / 7.2/
      nr1 = nr-1
       mm = nmu
      nw1 = nw-1
ccc   do 11 i=1,nr  
      do 21 j=0,nmu
       acc0= -amu(j)**2*dldt(i) + 0.5*sint(j)**2*dbdt(i)
       acc1= -amu(j)*dvdt(i)
      do 31 k=0,nw
	ae = 1.+vl(i)*ww(k)*amu(j)/cau**2
       acc = dt*(acc0 + acc1/ww(k))/dlnp/ae
      if (acc.ge.0.) then
	 if (k.eq.0) then
         df(i,j,k) = df(i,j,k) - acc*f(i,j,k)
	 else
	 if (k.eq.1) then
	 df(i,j,k) = df(i,j,k) - acc*(f(i,j,k)-f(i,j,k-1)) 
	 else
	 df(i,j,k) = df(i,j,k) - acc*(f(i,j,k)-f(i,j,k-1)) 
	 endif
	 endif
       else
	 if (k.eq.nw) then
	 df(i,j,k) = df(i,j,k) + acc*f(i,j,k) 
	 else
	 if (k.eq.nw1) then
	 df(i,j,k) = df(i,j,k) + acc*(f(i,j,k)-f(i,j,k+1)) 
	 else
	 df(i,j,k) = df(i,j,k) + acc*(f(i,j,k)-f(i,j,k+1)) 
	 end if
	 end if
      end if
	 df(i,j,k) = df(i,j,k) + acc*dlnp*gpeel(k)*f(i,j,k)
31    continue
21    continue
11    continue
      return
      end
c  *****************  end subroutine DELTA-P   ******************
                                                 !!! REWORK DLOGT!!!
      subroutine pcycle(i,dt,wght)
      common /size  / nr,nmu,nw, dim1
      common /pitch / mm,amu(0:40),sint(0:40),dmu
      common /impuls/ pmin,pmax,ppin,dlnp,pp(0:600)
      common /speed / wmin,wmax,wwin,ww(0:600)
      common /spiral/ tll(0:1000), dbdl(0:1000),vl(0:1000)
      common /coeff / dvdt(0:1000),dldt(0:1000),dbdt(0:1000),
     1                dndt(0:1000)
      common /solutn/ f(0:1000,0:40,0:600),df(0:1000,0:40,0:600)
      common /peel /  pex,fpeel(0:600),gpeel(0:600)
      dimension a2(0:1000),a1(0:1000),bb(0:1000),
     1          c1(0:1000),c2(0:1000),xy(0:1000)
      data  cau / 7.2/
      nr1 = nr-1
      nw1 = nw-1
ccc   do 10 i=1,nr 
      do 20 j=0,mm
       acc0= - amu(j)**2*dldt(i)+0.5*sint(j)**2*dbdt(i)
       acc1= - amu(j)*dvdt(i)

      do 31 k=0,nw 
	ae = 1.+vl(i)*ww(k)*amu(j)/cau**2
       acc = acc0+acc1/ww(k)
       acc = wght*dt*acc/ae
       pac = acc*gpeel(k)
       acc = acc/dlnp
       pac = 0.

       if (acc.gt.0.) then
	 if (k.eq.0) then
	   bb(k) = 1.+acc - pac
	 else
	   a2(k) = 0.0    
	   a1(k) =-acc
	   bb(k) = 1.+acc - pac
	 endif
	 c1(k) = 0.
	 c2(k) = 0.
       else
	 a2(k) = 0.
	 a1(k) = 0.
	 if (k.eq.nw) then
	   bb(k) = 1.-acc - pac
	 else
	   bb(k) = 1.-acc - pac
	   c1(k) = acc
	   c2(k) = 0.0            
	 endif
       endif
	 xy(k) = df(i,j,k)
31    continue

c --- solve and fill :
      call solve5(nw,a2,a1,bb,c1,c2,xy)
      do 32 k = 0,nw
      df(i,j,k) = xy(k)
32    continue
20    continue
10    continue
      return
      end
c  *****************  end subroutine P-CYCLE   ******************
      subroutine subsub(dt)
      common /size  / nr,nmu,nw, dim1
      common /suly  / wghtl,wghtmu,wghtw
      common /scatti/ qex,cmu(40),scmu(40),wsc(0:600),xsc(0:1000)
      common /spiral/ tll(0:1000), dbdl(0:1000),vl(0:1000)
      common /coeff / dvdt(0:1000),dldt(0:1000),dbdt(0:1000),
     1                dndt(0:1000)

      do 10 i=1,nr
      ksub = 1
      acc1 = dldt(i)
      acc2 = dbdt(i)

      dts = dt/float(ksub)
      do 23 kk = 1,ksub

      call deltamu(i,dts)
      call mucycle(i,dts,wghtmu)

      call deltap(i,dts)
      call pcycle(i,dts,wghtw)
      call updatsub(i)

23    continue
10    continue

      return
      end

c  *****************  end subroutine SUB-SUB   ******************
      subroutine update         
      common /size  / nr,nmu,nw, dim1
      common /solutn/ f(0:1000,0:40,0:600),df(0:1000,0:40,0:600)
      do 20 j=0,nmu
      do 30 k=0,nw
      do 10 i=1,nr
       f(i,j,k) = f(i,j,k) + df(i,j,k)
      df(i,j,k) = 0.
10    continue
       f(0,j,k) = 0.
      df(0,j,k) = 0.
30    continue
20    continue
      return 
      end
c  *****************  end subroutine UPDATE    ******************
      subroutine updatsub(i)    
      common /size  / nr,nmu,nw, dim1
      common /solutn/ f(0:1000,0:40,0:600),df(0:1000,0:40,0:600)
      do 20 j=0,nmu
      do 30 k=0,nw
       f(i,j,k) = f(i,j,k) + df(i,j,k)
      df(i,j,k) = 0.
30    continue
20    continue
      return 
      end
c  *****************  end subroutine UPDATE    ******************

c ==================== SOLVER ROUTINES: =========================== 

      subroutine solve5(nx,a2,a1,bx,c1,c2,xy)
      dimension a2(0:1000),a1(0:1000),bx(0:1000),
     1          c1(0:1000),c2(0:1000),xy(0:1000)
c --  solve penta-diagonal
      c1(0) = c1(0)/bx(0)
      c2(0) = c2(0)/bx(0)
      xy(0) = xy(0)/bx(0)
      aax = a1(1)
      bbx = bx(1)-aax*c1(0)
      c1(1) = (c1(1)-aax*c2(0))/bbx
      c2(1) =  c2(1)/bbx
      xy(1) = (xy(1)-aax*xy(0))/bbx
      do 11 i=2,nx 
      aa2 = a2(i)
      aa1 = a1(i)-aa2*c1(i-2)
      bbx = bx(i)-aa2*c2(i-2)-aa1*c1(i-1)
      c1(i) = (c1(i)-aa1*c2(i-1))/bbx
      c2(i) =  c2(i)/bbx
      xy(i) = (xy(i)-aa1*xy(i-1)-aa2*xy(i-2))/bbx
11    continue
 
      nx1 = nx-1
      xy(nx1) = xy(nx1)-c1(nx1)*xy(nx)
      do 21 i=nx-2,0,-1
      xy(i) = xy(i)-c1(i)*xy(i+1)-c2(i)*xy(i+2)
21    continue 
     
      return
      end
c  *****************  end subroutine SOLVE-5  ******************
      subroutine solvaa(nx,a2,a1,bx,xy)
      dimension a2(0:1000),a1(0:1000),bx(0:1000),xy(0:1000)
c --  solve lower diagonal
      xy(0) = xy(0)/bx(0)
      xy(1) = (xy(1)-a1(1)*xy(0))/bx(1)
      do 11 i=2,nx
      xy(i) = (xy(i)-a1(i)*xy(i-1)-a2(i)*xy(i-2))/bx(i)
11    continue 
     
      return
      end
c  *****************  end subroutine SOLV-AA  ******************
      subroutine solvcc(nx,bx,c1,c2,xy)
      dimension bx(0:1000),c1(0:1000),c2(0:1000),xy(0:1000)
c --  solve upper diagonal
      nx1 = nx-1
      xy(nx) = xy(nx)/bx(nx)
      xy(nx1) = (xy(nx1)-c1(nx1)*xy(nx))/bx(nx1)
      do 21 i=nx-2,0,-1
      xy(i) = (xy(i)-c1(i)*xy(i+1)-c2(i)*xy(i+2))/bx(i)
21    continue 
     
      return
      end

c  *****************  end subroutine SOLV-CC  ******************

c ========================  OUTPUT ROUTINES to be revised ========

      subroutine alla(io,jst,t)      
      common /size  / nr,nmu,nw, dim1
      common /blast/  slamb,tblast,tblst1,rshck1,dlnt
      common /gazdi/  ggamma,bbrad,vvmin,vvmax,ddmin,ddmax,
     1                ccmin,ccmax,aamin,aamax,bbmin,bbmax
      common /radio / nn,rmin,rshock,rmax,r(0:1000)
      common /plasma/ algbb(0:1000),algll(0:1000),algnn(0:1000),
     1                vr(0:1000)
      common /spiral/ tll(0:1000), dbdl(0:1000),vl(0:1000)
      common /coeff / dvdt(0:1000),dldt(0:1000),dbdt(0:1000),
     1                dndt(0:1000)
      common /scatti/ qex,cmu(40),scmu(40),wsc(0:600),xsc(0:1000)
      
      if (io.lt.10) stop  'io kicsi'
      if (io.gt.16) stop  'io nagy '
      if (io.eq.10) open(io,file = 'alla0.out')
      if (io.eq.11) open(io,file = 'alla1.out')
      if (io.eq.12) open(io,file = 'alla2.out')
      if (io.eq.13) open(io,file = 'alla3.out')
      if (io.eq.14) open(io,file = 'alla4.out')
      if (io.eq.15) open(io,file = 'alla5.out')
      if (io.eq.16) open(io,file = 'alla6.out')

      nn=nr
      write(io,901) jst,t
      write(io,902) jst,t,rshck1,slamb
      write(io,*)
      write(io,910)      
      write(io,*)
      do 10 i=0,nn
      write(io,911) i,r(i),vr(i),algbb(i),algll(i),algnn(i)
10    continue
      write(io,*)
      write(io,920)
      write(io,*)
      do 20 i=0,nn
      write(io,921) i,r(i),dvdt(i),dldt(i),dbdt(i),dndt(i)
20    continue
      write(io,*)
      write(io,930)
      write(io,*)
      do 30 i=0,nn
      write(io,931) i,r(i),tll(i),dbdl(i),xsc(i)
30    continue
      close(io)
901   format(5x,i5,' th step -- time[hour] = ',f10.3)
902   format(i6,3f12.6)
910   format(6x,'Radius - Vradial - logB - logL - LogN :')
911   format(i6,5f12.6)
920   format(6x,'Radius - DVDt - DLDt - DBDt - DNDt :')
921   format(i6,5f12.6)
930   format(6x,'Radius - TLL - DBDL - XSC(1/lambda): ')
931   format(i6,4f12.6)
      return
      end

c ****************************    end ALL output ******************

      subroutine csilla(io,ist,t)
      common /size  / nr,nmu,nw, dim1
      common /radio / nn,rmin,rshock,rmax,r(0:1000)
      common /azimut/ fi(0:1000)
      common /impuls/ pmin,pmax,ppin,dlnp,pp(0:600)
      common /energy/ emin,emax,eein,ee(0:600)
      common /speed / wmin,wmax,wwin,ww(0:600)
      common /pitch / mm,amu(0:40),sint(0:40),dmu
      common /solutn/ f(0:1000,0:40,0:600),df(0:1000,0:40,0:600)
      common /peel /  pex,fpeel(0:600),gpeel(0:600)
      common /obsrad/ krmax,krobs(5),robs(5)
      common /obserg/ kemax,keobs(5),eobs(5)
      dimension ff(0:1000,0:600),s(0:1000,0:600),ffx(4)
      dimension ract(10)
      if (io.lt.10) stop  'io kicsi'
      if (io.gt.16) stop  'io nagy '
      if (io.eq.10) open(io,file = 'csilla0.out')
      if (io.eq.11) open(io,file = 'csilla1.out')
      if (io.eq.12) open(io,file = 'csilla2.out')
      if (io.eq.13) open(io,file = 'csilla3.out')
      if (io.eq.14) open(io,file = 'csilla4.out')
      if (io.eq.15) open(io,file = 'csilla5.out')
      if (io.eq.16) open(io,file = 'csilla6.out')
      do 11 i=0,nr
      do 13 k=0,nw
      aa = 0.5*(f(i,0,k)+f(i,nmu,k))
      ss = 0.5*(f(i,0,k)-f(i,nmu,k))
      do 12 j=1,nmu-1
      aa = aa + f(i,j,k)
      ss = ss + f(i,j,k)*amu(j)
12    continue
      ff(i,k) = aa/float(nmu)
       s(i,k) = ss/float(nmu)
13    continue
11    continue
      do 21 kk=1,krmax
      iobs = krobs(kk)
      ract(kk) = r(iobs)
21    continue
      do 22 kk=1,kemax
      kobs = keobs(kk)
      eobs(kk) = ee(kobs)
22    continue
c --  energy spectra:
      write(io,*) 
      write(io,101) ist,t 
      write(io,102) krmax,(ract(kk),kk=1,krmax)
      write(io,*) 
      do 51 k = 0,nw
      do 52 kk = 1,krmax
	 ii = krobs(kk)
	 ffx(kk) = ff(ii,k)*fpeel(k)*pp(k)**2
52    continue
	 write(io,111) k,ee(k),(ffx(kk),kk=1,krmax)
51    continue
         write(io,*) 
c --  radial dependence :
      write(io,*) 
      write(io,201) ist,t
      write(io,202) kemax,(eobs(kk),kk=1,kemax)
      write(io,*) 
      do 61 i = 0,nr
      do 62 kk = 1,kemax
	 ke  = keobs(kk)
	 ffx(kk) = ff(i,ke)*fpeel(ke)*pp(ke)**2
62    continue
	 write(io,211) i,r(i),fi(i),(ffx(kk),kk=1,kemax)
61    continue
         write(io,*) 
c --  pitch - angles :
      write(io,*) 
      write(io,300) ist,t 
      do 71 kke = 1,kemax
      ke = keobs(kke)
      write(io,*) 
      write(io,301) eobs(kke)
      write(io,*) 
      do 72 j = 0,nmu
      do 73 kkr=1,krmax
      ii = krobs(kkr)
      ffx(kkr) =  f(ii,j,ke)/ff(ii,ke)
73    continue
        write(io,311) j,amu(j),(ffx(kk),kk=1,krmax)
72    continue
        write(io,*) 
71    continue
      close(io)
c --- report:
      ixx = io-10
      write(*,*) ixx, ' - output done at step/time[hour] :  ',ist,t
      return
c --  formats:
101   format(2x,'Energy spectra - after step/hour: ',i8,f10.2)
102   format(2x,'at ',i3,'  radii [au] : ',4f10.2)
111   format(i4,f12.6,4e14.4)
201   format(2x,'Radial dependence - after step/hour : ',i8,f10.2)
202   format(2x,'at ',i3,'  energy [MeV/n] : ',5f8.3)
211   format(2x,i4,2f8.3,5e14.4)
300   format(2x,'Pitch-angle distribitions after step/hour : ',i8,f10.2)
301   format(2x,'Energy [MeV/n} : ',f10.3)
311   format(2x,i3,f10.3,4f14.6)
      end

c  *****************  end subroutine CSILLA   ******************

      subroutine opentime 
      common /size  / nr,nmu,nw, dim1
      common /radio / nn,rmin,rshock,rmax,r(0:1000)
      common /impuls/ pmin,pmax,ppin,dlnp,pp(0:600)
      common /energy/ emin,emax,eein,ee(0:600)
      common /speed / wmin,wmax,wwin,ww(0:600)
      common /pitch / mm,amu(0:40),sint(0:40),dmu
      common /solutn/ f(0:1000,0:40,0:600),df(0:1000,0:40,0:600)
      common /obsrad/ krmax,krobs(5),robs(5)
      common /obserg/ kemax,keobs(5),eobs(5)
      do 20 ke=1,kemax
      kk = keobs(ke)
20    eobs(ke) = ee(kk)
      do 11 kr=1,krmax
      io = 20+kr
      ip = 30+kr
      if (kr.eq.1) open(io,file='timevar.r1')
      if (kr.eq.2) open(io,file='timevar.r2')
      if (kr.eq.3) open(io,file='timevar.r3')
      if (kr.eq.4) open(io,file='timevar.r4')
      if (kr.eq.5) open(io,file='timevar.r5')
      write(io,*)
      write(io,711) robs(kr) 
      write(io,712) (eobs(ke),ke=1,kemax) 
      write(io,*)
      if (kr.eq.1) open(ip,file='plasma.r1')
      if (kr.eq.2) open(ip,file='plasma.r2')
      if (kr.eq.3) open(ip,file='plasma.r3')
      if (kr.eq.4) open(ip,file='plasma.r4')
      if (kr.eq.5) open(ip,file='plasma.r5')
      write(ip,*)
      write(ip,711) robs(kr) 
      write(ip,722)  
      write(ip,*)
11    continue
711   format(5x,'Time-variation at radius[AU]: ',f10.2)
712   format(5x,'Energies [MeV]:',5f10.3) 
722   format(5x,'Plasma:') 
      return
      end

c  *****************  end subroutine OPENTIME  ******************

      subroutine closetime 
      common /obsrad/ krmax,krobs(5),robs(5)
      common /obserg/ kemax,keobs(5),eobs(5)
      do 11 kr=1,krmax
      io = 20+kr
      ip = 30+kr
      write(io,*)
      write(io,*) 'the end'
      close(io)
      write(ip,*)
      write(ip,*) 'the end'
      close(ip)
11    continue 
      return
      end

c  *****************  end subroutine OPENTIME  ******************

      subroutine timevar(jst,t)
      common /size  / nr,nmu,nw, dim1
      common /radio / nn,rmin,rshock,rmax,r(0:1000)
      common /impuls/ pmin,pmax,ppin,dlnp,pp(0:600)
      common /energy/ emin,emax,eein,ee(0:600)
      common /speed / wmin,wmax,dlnw,ww(0:600)
      common /pitch / mm,amu(0:40),sint(0:40),dmu
      common /solutn/ f(0:1000,0:40,0:600),df(0:1000,0:40,0:600)
      common /peel /  pex,fpeel(0:600),gpeel(0:600)
      common /obsrad/ krmax,krobs(5),robs(5)
      common /obserg/ kemax,keobs(5),eobs(5)
      common /plasma/ algbb(0:1000),algll(0:1000),algnn(0:1000),
     1                vvr(0:1000)
      common /coeff / dvdt(0:1000),dldt(0:1000),dbdt(0:1000),
     1                dndt(0:1000)
      common /convrt/ aukm,hour,valf
      dimension ffe(10)

      do 10 kr=1,krmax
      io = 20+kr
      ip = 30+kr
      rr = robs(kr)
      ff0 = 1.0

      i=0
20    if (r(i).lt.rr) then
      if (i.eq.nr) stop 'miako tyukanyo'
      i=i+1
      go to 20
      endif
      i1 = i-1
      i2 = i
      fr1 = (rr-r(i1))/(r(i2)-r(i1))
      fr2 = 1.-fr1
      krobs(kr) = i1

      do 30 ke=1,kemax
      kk = keobs(ke)
      ss1 = 0.5*(f(i1,0,kk)+f(i1,nmu,kk))
      ss2 = 0.5*(f(i2,0,kk)+f(i2,nmu,kk))
      do 40 jj=1,nmu-1
      ss1 = ss1+f(i1,jj,kk)
      ss2 = ss2+f(i2,jj,kk)
40    continue
      ss1 = ss1/float(nmu)
      ss2 = ss2/float(nmu)
      ffe(ke)=ff0*(fr2*ss1+fr1*ss2)*fpeel(kk)*pp(kk)**2
30    continue
      
      vv =  fr2*vvr(i1) + fr1*vvr(i2)
      vv = aukm/hour*vv
      abb= fr2*algbb(i1)+fr1*algbb(i2)
      all= fr2*algll(i1)+fr1*algll(i2)
      add= abb-all
      ann= fr2*algnn(i1)+fr1*algnn(i2)
      write( *,788) jst,t,vv,add,all,abb,ann
      write(ip,788) jst,t,vv,add,all,abb,ann

      write(io,777) jst,t,(ffe(ke),ke=1,kemax)
      write( *,777) jst,t,(ffe(ke),ke=1,kemax)
10    continue
      return
777   format(i5,f8.2,5e14.4)
788   format(i5,f8.2,5f14.6)
      end

c  *****************  end subroutine TIMEVAR  ******************

