!     All routines here include the following 'param.h' header
!     integer::nRMax,nMuMax,nPMax
!     parameter(nRMax=1000,nMuMax=40,nPMax=600)

!BOP
!ROUTINE: SP_cool - set functions of the energy coordinate
!INTERFACE:
      subroutine SP_cool
!DESCRIPTION:
!Subroutine SP_cool sets pmin,pmax,pinj,dlnp,
!arrays pp,Speed_I,ee - momentum,speed and energy as a 
!function of the energetic coordinate.
!Also sets the scattering length dependence on energy
!EOP 
      include 'param.h'   
      common /SP_size  / nr,nmu,nw, dim1
      common /SP_partid/ iz,massa,ekpp,xlmbda0
      common /SP_impuls/ pmin,pmax,ppin,dlnp,Momentum_I(0:nPMax)
      common /SP_energy/ emin,emax,eein,ee(0:nPMax)
      common /SP_speed / wmin,wmax,wwin,Speed_I(0:nPMax)
      common /SP_scatti/ qex,Scattering_I(nMuMax),wsc(0:nPMax),
     1     xsc(0:nRMax)
      common /SP_quelle/ kinj,einj,pinj,qqp(0:nPMax),qqx(0:nRMax)
      common /SP_convrt/ cAUKm,hour,valf
      include 'stdout.h'
      data  pmass,clight,xkm  / 938., 3.e10, 1.e5 /
c ----------------------------------- scales:
      pi = 2.*asin(1.)
      cAUKm = 1.5e8
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
      wwin = btin*clight*hour/xkm/cAUKm

      do 11 k = 0,nw
      Momentum_I(k) = pmin*exp(float(k)*dlnp)
       etot = sqrt(Momentum_I(k)**2+e0**2)
       beta = Momentum_I(k)/etot
      ee(k) = Momentum_I(k)**2/(etot+e0)
      Speed_I(k) = beta*clight*hour/xkm/cAUKm
      qqp(k)= 0.
       rg   = float(massa)*Momentum_I(k)/abs(qz)/gv
      wsc(k)= Speed_I(k)*rg**expo

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

c  *****************  end subroutine SP_cool     ******************
!BOP
!ROUTINE: SP_peeloff - multiples distribution function by p^qex if desired
!INTERFACE:
      subroutine SP_peeloff
!EOP
      include 'param.h'   
      common /SP_size  / nr,nmu,nw, dim1
      common /SP_impuls/ pmin,pmax,ppin,dlnp,Momentum_I(0:nPMax)
      common /SP_energy/ emin,emax,eein,ee(0:nPMax)
      common /SP_speed/  wmin,wmax,wwin,Speed_I(0:nPMax)
      common /SP_scatti/ qex,Scattering_I(nMuMax),wsc(0:nPMax),
     1     xsc(0:nRMax)
      common /SP_quelle/ kinj,einj,winj,qqw(0:nPMax),qqx(0:nRMax)
      common /SP_peel /  pex,fpeel(0:nPMax),gpeel(0:nPMax)
      common /SP_solutn/ f(0:nRMax,0:nMuMax,0:nPMax),
     1     df(0:nRMax,0:nMuMax,0:nPMax)
      include 'stdout.h'
      pex = 0.
      if(DoWriteAll)write(iStdout,*) prefix,
     1      'PEEL-OFF -- what exponent ??',pex
      do 30 k = 0,nw
      fact = (Momentum_I(k)/ppin)**pex
      fpeel(k) = 2./(1.+fact)
      gpeel(k) = pex*fact/(1.+fact)
      if(DoWriteAll)write(iStdout,911) prefix,
     1        pex,k,ee(k),fact,gpeel(k),fpeel(k)
911   format(a,f10.2,i5,3f12.6,e14.4)
      qqw(k) = qqw(k)/fpeel(k)
      do 10 i=0,nr
      do 20 j=0,nmu
      f(i,j,k) = f(i,j,k)/fpeel(k)
20    continue
10    continue
30    continue

      return
      end

c  *****************  end subroutine SP_peeloff  ******************
!BOP
!ROUTINE: SP_pangle - sets \mu, 1-\mu^2, \mu-dependence for scattering
!INTERFACE:
      subroutine SP_pangle
!EOP
      include 'param.h'    
      common /SP_size  / nr,nmu,nw, dim1
      common /SP_pitch / mm,Mu_I(0:nMuMax),SinMu2_I(0:nMuMax),dmu
      common /SP_scatti/ qex,Scattering_I(nMuMax),wsc(0:nPMax),
     1     xsc(0:nRMax)
      include 'stdout.h'
c     linear grid in mu=cost
      mm = nmu
       m = mm/2

      !Uniformly spaced grid, -1<=\mu<=1, index value 0 is for \mu=1
      !index nMu is for \mu=-1, d\mu is positive, 
      dmu = 2./float(mm) 
      do 11 jj=0,m
      jm = mm-jj
      Mu_I(jj) = 1. - float(jj)*dmu
      Mu_I(jm) = - Mu_I(jj)
      SinMu2_I(jj) = 1.-Mu_I(jj)**2
      SinMu2_I(jm) = SinMu2_I(jj)
11    continue
       Mu_I(0) = 1.
      SinMu2_I(0) = 0.
       Mu_I(m) = 0.
      SinMu2_I(m) = 1.
      Mu_I(mm) = -1
      SinMu2_I(mm) = 0.
c  
      do 21 jj=1,m
      jm = mm-jj+1
        ccmu = (Mu_I(jj)+Mu_I(jj-1))/2.
        ssmu = 1.-ccmu**2
       absmu = abs(ccmu)
       scatty = ssmu*absmu**qex*3./(1.-qex)/(3.-qex)
       Scattering_I(jj) = scatty/dmu**2/2.   
       Scattering_I(jm) = scatty/dmu**2/2.
21    continue
      !Diffusion coefficient is as follows:
      !D_{\mu\mu}=xscat(r)*fact(energy)*(1-\mu^2)|\mu|^q
      !q can be related to the spectral index of turbulence
c ---  scatty representing  (1/lambda)
      return
      end

c  *****************  end subroutine SP_pangle    ******************
!BOP
!ROUTINE: SP_initial - nullify distribution function and related arrays 
!INTERFACE:
      subroutine SP_initial
!EOP                   
      include 'coupler.h'
      common /SP_size /  nr,nmu,nw, dim1
      common /SP_solutn/ f(0:nRMax,0:nMuMax,0:nPMax),
     1     df(0:nRMax,0:nMuMax,0:nPMax)
      include 'stdout.h'
      integer iFile

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
10    continue
      if(DoRestart)then
         call CON_io_unit_new(iFile)
         open(iFile,file='./SP/restartIN/f',
     1        status='old',form='unformatted')
         read(iFile)f
         close(iFile)
      end if
      return
      end
!BOP
!ROUTINE: SP_save_f - saves the distribution function 
!INTERFACE:
      subroutine SP_save_f
!EOP                   
      include 'coupler.h'
      common /SP_size /  nr,nmu,nw, dim1
      common /SP_solutn/ f(0:nRMax,0:nMuMax,0:nPMax),
     1     df(0:nRMax,0:nMuMax,0:nPMax)
      include 'stdout.h'
      integer iFile
         call CON_io_unit_new(iFile)
         open(iFile,file='./SP/restartOUT/f',
     1          status='replace',FORM='unformatted')
         write(iFile)f
         close(iFile)
      return
      end

c **************************  end of INiT ****************
!BOP
!ROUTINE: SP_helios - sets multipliers in the kinetic equation coming from MHD
!INTERFACE:
      subroutine SP_helios(t,dt,jj,kk)
!DESCRIPTION:
!Uses the values obtained from MHD and calculates the Lagrangian time 
!derivatives which are multipliers in the kinetic equation.
!Sets actual values for spatial dependance of the scattering length
!Input parameters: t - is not used, dt - time step, jj - is not used
!is not used in the framework
!EOP
      include 'coupler.h'
      common /SP_size / nr,nmu,nw, dim1
      common /SP_gazdi/ ggamma,bbrad,vvmin,vvmax,ddmin,ddmax,
     1               ccmin,ccmax,aamin,aamax,bbmin,bbmax
      common /SP_elem / zr(0:nRMax),zv(0:nRMax),zp(0:nRMax),zn(0:nRMax)
      common /SP_magia/ qb(0:nRMax),qd(0:nRMax),zb(0:nRMax) 
      
      common /SP_blast/  slamb,tblast,tblst1,rshck1,dlnt
      common /SP_radio / nn,rmin,rshock,rmax,r(0:nRMax)
!      common /scphys/ wind,omega,xscatt1 !'wind' is undefined
      common/SP_scphys/ omega,xscatt1
      common /SP_spiral/ tll(0:nRMax), dbdl(0:nRMax),vl(0:nRMax)
      common /SP_coeff / dvdt(0:nRMax),dldt(0:nRMax),dbdt(0:nRMax),
     1                dndt(0:nRMax)
      common /SP_quelle/ kinj,einj,winj,qqw(0:nPMax),qqx(0:nRMax)
      common /SP_scatti/ qex,Scattering_I(nMuMax),wsc(0:nPMax),
     1     xsc(0:nRMax)
      common /SP_smile / mp,ni
     1       ,eta(0:6000),exx(0:6000),ezz(0:6000),efi(0:6000)
     2       ,evr(0:6000),evx(0:6000),evz(0:6000)
     3       ,edsm(0:6000),edfi(0:6000)
     4       ,ebm(0:6000),ebr(0:6000),ebt(0:6000),ebfi(0:6000)
     5       ,ebx(0:6000),ebz(0:6000),edd(0:6000)
     6       ,dbmds(0:6000),dbfids(0:6000)
      logical UseSelfSimilarity,UseRefresh
      common/SP_log/UseSelfSimilarity,UseRefresh
      include 'stdout.h'
c +++++++++++++++++++++++++++++++++++++++     concept:   +++++ 
c    Transfers between similarity solution and actual, 
!    as long as the self-similar solution is used.
!    Otherwise needs the data through the same set of Lagrangian 
!    points to calculate the Lagrangian time derivatives.
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

      pi = 2.*asin(1.)

c -----------------------------------------------------
      if(UseSelfSimilarity)
     1     rshck1 = rshock*(tblst1/tblast)**slamb 
      
      do 10 i=0,nr
         if(UseSelfSimilarity)then
            !     The lagrangian derivatives and 
            !  scattering length is taken from self-similar soln   
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
            dldt(i) = (all-algll(i))/dt
            dvdt(i) = (vvr-vr(i))/dt
         else
            !The dirvatives and diffusion coefficient are
            !taken from the updated coupling variables and
            !from previous ones
            !NOTE: arrays rx,ry,...pp start from the index value 1
            !tll,d*dt d*ds - from index value 0
            r(i)=sqrt(rx(i+1)**2
     1              +ry(i+1)**2
     2              +rz(i+1)**2)
            if(i.eq.0)then
               tll(i)=sqrt((rx(2)-rx(1))**2
     1              +(ry(2)-ry(1))**2
     2              +(rz(2)-rz(1))**2)
            elseif(i.eq.nr)then
               tll(i)=sqrt((rx(nr)-rx(nr+1))**2
     1              +(ry(nr)-ry(nr+1))**2
     2              +(rz(nr)-rz(nr+1))**2)
            else
               tll(i)=sqrt((rx(i+2)-rx(i))**2
     1              +(ry(i+2)-ry(i))**2
     2              +(rz(i+2)-rz(i))**2)/2.0
            end if
            abb=sqrt(bx(i+1)**2
     1              +by(i+1)**2
     2              +bz(i+1)**2)
            dvdt(i)=(bx(i+1)*(vx(i+1)-vxOld(i+1))
     1              +by(i+1)*(vy(i+1)-vyOld(i+1))
     2              +bz(i+1)*(vz(i+1)-vzOld(i+1)))/(abb*dt)
            abb=alog(abb)
            ann=alog(dens(i+1))
            !Each used value of vx,vy,vz should be saved:
            vxOld(i+1)=vx(i+1)
            vyOld(i+1)=vy(i+1)
            vzOld(i+1)=vz(i+1)
            !vr is used in the output files:
            vr(i) =       vx(i+1)*rx(i+1)
            vr(i) = vr(i)+vy(i+1)*ry(i+1) 
            vr(i) = vr(i)+vz(i+1)*rz(i+1)
            vr(i) = vr(i)/r(i)
         end if
            dbdt(i) = (abb-algbb(i))/dt
            dndt(i) = (ann-algnn(i))/dt
            
            
            vl(i) = 0.                            
c ??????????????????????????????????????????????????????
      if (i.lt.0) then
      write(iStdout,*) prefix,
     1         'I AM HERE IN SP_helios - kk,i,ii=  ', kk,i,ii
      write(iStdout,*) prefix,
     1           'R:       ',r(i)
      write(iStdout,*) prefix,
     1          'EBM-Fi : ',ebm(ii),ebfi(ii)
      write(iStdout,*) prefix,
     1                     'BBM-Fi : ',bbm,bbfi,bbb
      write(iStdout,*) prefix,
     1          'BB     : ',abb,algbb(i),dbdt(i)
      write(iStdout,*) prefix,
     1          'LL     : ',all,algll(i),dldt(i)
      write(iStdout,*) prefix,
     1         'NN     : ',ann,algnn(i),dndt(i)
      read(*,*)    lull
      endif
c ??????????????????????????????????????????????????????

      algbb(i) = abb
      algnn(i) = ann

      if(UseSelfSimilarity)then
         algll(i) = all
         vr(i) = vvr                  
      else
         dldt(i)=dbdt(i)-dndt(i)
      end if

      qqx(i) = 0.        
         bbfac = 1./r(i)
       bbfak = bbfac
      xsc(i) = xscatt1*bbfak
      if (r(i).lt.1.15*rshck1) then
         xsc(i) = 50.*xsc(i)
         qqx(i) = 1/r(i)**2
      endif
      if (r(i).lt.0.50*rshck1) then
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
      ibmax = 1
      ibmin = 1
      ilmax = 1
      ilmin = 1
      dbmax = dbdt(1)
      dbmin = dbdt(1)
      dlmax = dldt(1)
      dlmin = dldt(1)
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
      if(DoWriteAll)write(iStdout,*) prefix,
     1         'Min dldt: ',ilmin,dlmin
      if(DoWriteAll)write(iStdout,*) prefix,
     1          'Max dldt: ',ilmax,dlmax
      if(DoWriteAll)write(iStdout,*) prefix,
     1           'Min dbdt: ',ibmin,dbmin
      if(DoWriteAll)write(iStdout,*) prefix,
     1           'Max dbdt: ',ibmax,dbmax
      return
      end
!===============================================================
!BOP
!ROUTINE: SP_source - calculates the SP_source of injected particles
!INTERFACE: 
      subroutine SP_source(dt)
!EOP
      include 'param.h'
      common /SP_size  / nr,nmu,nw, dim1
      common /SP_pitch / mm,Mu_I(0:nMuMax),SinMu2_I(0:nMuMax),dmu
      common /SP_impuls/ pmin,pmax,ppin,dlnp,Momentum_I(0:nPMax)
      common /SP_solutn/ f(0:nRMax,0:nMuMax,0:nPMax),
     1     df(0:nRMax,0:nMuMax,0:nPMax)
      common /SP_quelle/ kinj,einj,pinj,qqp(0:nPMax),qqx(0:nRMax)
      include 'stdout.h'
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
c  *****************  end subroutine SP_source   ******************

!BOP
!ROUTINE: SP_deltal - SP_update the input from spatial derivative of f
!INTERFACE:
      subroutine SP_deltal(dt)                      
!DESCIPTION
!Advance the solution of the following equation
!\partial f/\partial t+ v\mu\partial f/\partial s =0
!through a time step
!EOP
!!! REWORK DLOGT!!!
      include 'param.h'
      common /SP_size  / nr,nmu,nw, dim1
      common /SP_pitch / mm,Mu_I(0:nMuMax),SinMu2_I(0:nMuMax),dmu
      common /SP_speed / wmin,wmax,wwin,Speed_I(0:nPMax)
      common /SP_spiral/ tll(0:nRMax), dbdl(0:nRMax),vl(0:nRMax)
      common /SP_coeff / dvdt(0:nRMax),dldt(0:nRMax),dbdt(0:nRMax),
     1                dndt(0:nRMax)
      common /SP_solutn/ f(0:nRMax,0:nMuMax,0:nPMax),
     1     df(0:nRMax,0:nMuMax,0:nPMax)
      include 'stdout.h'
      data  cau / 7.2/

      m = mm/2
      nr1 = nr-1

      do 30 k=0,nw
      do 20 jj=0,m-1

c -------------------------------- out direction(s):
      j = jj
       df(0,j,k) = 0.
      vv =  Speed_I(k)*Mu_I(j)
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
      vv = Speed_I(k)*Mu_I(j)
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
c  *****************  end subroutine SP_deltal  ******************
                                                 !!! REWORK DLOGT!!!
!BOP
!ROUTINE: SP_lcycle - corrects the increment from df/ds using implicit scheme
!INTERFACE:
      subroutine SP_lcycle(dt,wght)
!DESCRIPTION:
!Implicit scheme with two-diagonal matrix (oly one lower diagonal or
!only one higher diagonal, depending on the \mu sign) is solved directly
!EOP
      include 'param.h'
      common /SP_size  / nr,nmu,nw, dim1
      common /SP_pitch / mm,Mu_I(0:nMuMax),SinMu_I(0:nMuMax),dmu
      common /SP_speed / wmin,wmax,wwin,Speed_I(0:nPMax)
      common /SP_spiral/ tll(0:nRMax), dbdl(0:nRMax),vl(0:nRMax)
      common /SP_coeff / dvdt(0:nRMax),dldt(0:nRMax),dbdt(0:nRMax),
     1                dndt(0:nRMax)
      common /SP_solutn/ f(0:nRMax,0:nMuMax,0:nPMax),
     1     df(0:nRMax,0:nMuMax,0:nPMax)
      include 'stdout.h'
      data  cau / 7.2/

      m = mm/2
      nr1 = nr-1
      do 30 k=0,nw
      do 20 jj=0,m-1

      j = jj
      vv = Speed_I(k)*Mu_I(j)

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
      vv = Speed_I(k)*Mu_I(j)

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
!BOP
!ROUTINE: SP_deltamu - SP_update the input from pitch angle derivative of f
!INTERFACE:                      
      subroutine SP_deltamu(i,dt)
!DESCIPTION
!Advances the solution of the following equation
!\end{verbatim}
!\begin{equation}
!\partial f/\partial t+ (1-\mu^2)\left(-\frac12v\frac{\partial ln B}
!{\partial s}-\mu(frac12\frac{d ln B}{d t}+\frac{d ln \delta s}{dt}-
!{\bf b}\cdot \frac{d{\bf u}}{dt}\left)
!{\partial f/\partial \mu =\frac{\partial}{\partial \mu}
!\left(D_{\mu\mu}\frac{\partial f}{\partial \mu}\right)
!\end{equation}
!\begin{verbatim}
!throught one time step
!EOP  
      include 'param.h'
      common /SP_size  / nr,nmu,nw, dim1
      common /SP_pitch / mm,Mu_I(0:nMuMax),SinMu2_I(0:nMuMax),dmu
      common /SP_speed / wmin,wmax,wwin,Speed_I(0:nPMax)
      common /SP_spiral/ tll(0:nRMax), dbdl(0:nRMax),vl(0:nRMax)
      common /SP_coeff / dvdt(0:nRMax),dldt(0:nRMax),dbdt(0:nRMax),
     1                dndt(0:nRMax)
      common /SP_scatti/ qex,Scattering_I(nMuMax),wsc(0:nPMax),
     1     xsc(0:nRMax)
      common /SP_solutn/ f(0:nRMax,0:nMuMax,0:nPMax),
     1     df(0:nRMax,0:nMuMax,0:nPMax)
      include 'stdout.h'
      data  cau / 7.2/
      nr1 = nr-1
       mm = nmu
      mm1 = mm-1
ccc   do 11 i=1,nr
      do 31 k=0,nw
      
      fact = xsc(i)*wsc(k)*dt
      df(i,0,k) = df(i,0,k) + 2.*fact*Scattering_I(1)*(f(i,1,k)-f(i,0,k))
      df(i,mm,k)= df(i,mm,k)+ 2.*fact*Scattering_I(mm)*(f(i,mm1,k)-f(i,mm,k))

      do 21 j=1,mm1
      ae = 1.+ Speed_I(k)*Mu_I(j)*vl(i)/cau**2
      dte= dt/ae
      
      foc =-SinMu2_I(j)*(0.5*Speed_I(k)*dbdl(i)+dvdt(i)/Speed_I(k) +
     1                  Mu_I(j)*(dldt(i) + 0.5*dbdt(i)))
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

      sc1 = fact/ae*Scattering_I(j)
      sc2 = fact/ae*Scattering_I(j+1)
      df(i,j,k) = df(i,j,k) +
     1           sc1*(f(i,j-1,k)-f(i,j,k))+sc2*(f(i,j+1,k)-f(i,j,k))

21    continue
31    continue
11    continue
      return
      end
c  *****************  end subroutine DELTA-MU  ******************
                                                 !!! REWORK DLOGT!!!
!BOP
!ROUTINE: SP_mucycle - corrects the increment from df/d\mu using implicit scheme
!INTERFACE:
      subroutine SP_mucycle(i,dt,wght)
!EOP
      include 'param.h'
      common /SP_size  / nr,nmu,nw, dim1
      common /SP_pitch / mm,Mu_I(0:nMuMax),SinMu2_I(0:nMuMax),dmu
      common /SP_speed / wmin,wmax,wwin,Speed_I(0:nPMax)
      common /SP_scatti/ qex,Scattering_I(nMuMax),wsc(0:nPMax),
     1     xsc(0:nRMax)
      common /SP_spiral/ tll(0:nRMax), dbdl(0:nRMax),vl(0:nRMax)
      common /SP_coeff / dvdt(0:nRMax),dldt(0:nRMax),dbdt(0:nRMax),
     1                dndt(0:nRMax)
      common /SP_solutn/ f(0:nRMax,0:nMuMax,0:nPMax),
     1     df(0:nRMax,0:nMuMax,0:nPMax)
      dimension a2(0:nRMax),a1(0:nRMax),bb(0:nRMax),
     1          c1(0:nRMax),c2(0:nRMax),xy(0:nRMax)
      include 'stdout.h'
      data  cau / 7.2/
      nr1 = nr-1
       mm = nmu
      mm1 = mm-1
      do 31 k=0,nw
       fact = wght*wsc(k)*xsc(i)*dt
c --  mu = 1  (or j=0)
       ae = 1.+Speed_I(k)*vl(i)/cau**2
        sc1 = fact*Scattering_I(1)/ae
      a2(0) = 0.
      a1(0) = 0.
      bb(0) = 1. + 2.*sc1
      c1(0) = -2.*sc1
      c2(0) = 0.
      xy(0) = df(i,0,k)
c --  mu =-1  (or j=nmu)
       ae = 1.-Speed_I(k)*vl(i)/cau**2
        sc1 = fact*Scattering_I(mm)/ae
      a2(nmu) = 0.
      a1(nmu) = -2.*sc1
      bb(nmu) = 1. + 2.*sc1
      c1(nmu) = 0.
      c2(nmu) = 0.
      xy(nmu) = df(i,nmu,k)
      do 21 j=1,mm1
       ae = 1.+Speed_I(k)*vl(i)*Mu_I(j)/cau**2
c --- focusing :
       foc =-SinMu2_I(j)*(0.5*Speed_I(k)*dbdl(i)+dvdt(i)/Speed_I(k) +
     1                   Mu_I(j)*(dldt(i) + 0.5*dbdt(i) ))
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
      sc1 = fact*Scattering_I(j)/ae
      sc2 = fact*Scattering_I(j+1)/ae
      a1(j) = a1(j) -sc1 
      bb(j) = bb(j) + sc1 + sc2
      c1(j) = c1(j) -sc2
21    continue
c --- now solve and fill :
      call SP_solve5(mm,a2,a1,bb,c1,c2,xy)
      do 22 j=0,mm 
      df(i,j,k) = xy(j)
22    continue
31    continue
11    continue
      return
      end
c  *****************  end subroutine MU-CYCLE  ******************
                                                 !!! REWORK DLOGT!!!
                      
!BOP
!ROUTINE: SP_deltap - SP_update the input from the momentum derivative of f
!INTERFACE:
      subroutine SP_deltap(i,dt)
!DESCIPTION
!Updates the solution of the following equation
!\end{verbatim}
!\begin{equation}
!d f/d t+ \left((1-\mu^2)frac12\frac{d ln B}{d t}-
!-\mu^2\frac{d ln \delta s}{dt}-
!\frac{\mu{\bf b}}v\cdot \frac{d{\bf u}}{dt}\left)
!p{\partial f/\partial p = 0
!\end{equation}
!\begin{verbatim}
!throught one time step
!EOP  
      include 'param.h'
      common /SP_size  / nr,nmu,nw, dim1
      common /SP_pitch / mm,Mu_I(0:nMuMax),SinMu2_I(0:nMuMax),dmu
      common /SP_impuls/ pmin,pmax,ppin,dlnp,Momentum_I(0:nPMax)
      common /SP_speed / wmin,wmax,wwin,Speed_I(0:nPMax)
      common /SP_spiral/ tll(0:nRMax), dbdl(0:nRMax),vl(0:nRMax)
      common /SP_coeff / dvdt(0:nRMax),dldt(0:nRMax),dbdt(0:nRMax),
     1                dndt(0:nRMax)
      common /SP_solutn/ f(0:nRMax,0:nMuMax,0:nPMax),
     1     df(0:nRMax,0:nMuMax,0:nPMax)
      common /SP_peel /  pex,fpeel(0:nPMax),gpeel(0:nPMax)
      data  cau / 7.2/
      nr1 = nr-1
       mm = nmu
      nw1 = nw-1
ccc   do 11 i=1,nr  
      do 21 j=0,nmu
       acc0= -Mu_I(j)**2*dldt(i) + 0.5*SinMu2_I(j)*dbdt(i)
       acc1= -Mu_I(j)*dvdt(i)
      do 31 k=0,nw   
      ae = 1.+vl(i)*Speed_I(k)*Mu_I(j)/cau**2   
       acc = dt*(acc0 + acc1/Speed_I(k))/dlnp/ae
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
!BOP
!ROUTINE: SP_pcycle - corrects the increment from df/dp using implicit scheme
!INTERFACE:
      subroutine SP_pcycle(i,dt,wght)
!EOP
      include 'param.h'
      common /SP_size  / nr,nmu,nw, dim1
      common /SP_pitch / mm,Mu_I(0:nMuMax),SinMu2_I(0:nMuMax),dmu
      common /SP_impuls/ pmin,pmax,ppin,dlnp,Momentum_I(0:nPMax)
      common /SP_speed / wmin,wmax,wwin,Speed_I(0:nPMax)
      common /SP_spiral/ tll(0:nRMax), dbdl(0:nRMax),vl(0:nRMax)
      common /SP_coeff / dvdt(0:nRMax),dldt(0:nRMax),dbdt(0:nRMax),
     1                dndt(0:nRMax)
      common /SP_solutn/ f(0:nRMax,0:nMuMax,0:nPMax),
     1     df(0:nRMax,0:nMuMax,0:nPMax)
      common /SP_peel /  pex,fpeel(0:nPMax),gpeel(0:nPMax)
      dimension a2(0:nRMax),a1(0:nRMax),bb(0:nRMax),
     1          c1(0:nRMax),c2(0:nRMax),xy(0:nRMax)
      include 'stdout.h'
      data  cau / 7.2/
      nr1 = nr-1
      nw1 = nw-1
ccc   do 10 i=1,nr 
      do 20 j=0,mm
       acc0= - Mu_I(j)**2*dldt(i)+0.5*SinMu2_I(j)**2*dbdt(i)
       acc1= - Mu_I(j)*dvdt(i)

      do 31 k=0,nw
      ae = 1.+vl(i)*Speed_I(k)*Mu_I(j)/cau**2
       acc = acc0+acc1/Speed_I(k)
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
      call SP_solve5(nw,a2,a1,bb,c1,c2,xy)
      do 32 k = 0,nw
      df(i,j,k) = xy(k)
32    continue
20    continue
10    continue
      return
      end
c  *****************  end subroutine P-CYCLE   ******************
!BOP
!ROUTINE: SP_subsub - collect the inputs from \mu and momentum derivatives
      subroutine SP_subsub(dt)
!EOP
      include 'param.h'
      common /SP_size  / nr,nmu,nw, dim1
      common /SP_suly  / wghtl,wghtmu,wghtw
      common /SP_scatti/ qex,Scattering_I(nMuMax),wsc(0:nPMax),
     1     xsc(0:nRMax)
      common /SP_spiral/ tll(0:nRMax), dbdl(0:nRMax),vl(0:nRMax)
      common /SP_coeff / dvdt(0:nRMax),dldt(0:nRMax),dbdt(0:nRMax),
     1                dndt(0:nRMax)
      include 'stdout.h'
      do 10 i=1,nr
      ksub = 1
      acc1 = dldt(i) !To be used to calculate ksub - under development
      acc2 = dbdt(i)

      dts = dt/float(ksub)
      do 23 kk = 1,ksub

      call SP_deltamu(i,dts)
      call SP_mucycle(i,dts,wghtmu)

      call SP_deltap(i,dts)
      call SP_pcycle(i,dts,wghtw)
      call SP_updatsub(i)

23    continue
10    continue

      return
      end

c  *****************  end subroutine SUB-SUB   ******************
!BOP
!ROUTINE: SP_update - adds df to f
!INTERFACE:
      subroutine SP_update        
!EOP
      include 'param.h'
      common /SP_size  / nr,nmu,nw, dim1
      common /SP_solutn/ f(0:nRMax,0:nMuMax,0:nPMax),
     1     df(0:nRMax,0:nMuMax,0:nPMax)
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
c  *****************  end subroutine SP_update    ******************
!BOP
!ROUTINE: SP_update - adds df to f for a given spatial point
!INTERFACE:
      subroutine SP_updatsub(i)    
!EOP
      include 'param.h'
      common /SP_size  / nr,nmu,nw, dim1
      common /SP_solutn/ f(0:nRMax,0:nMuMax,0:nPMax),
     1     df(0:nRMax,0:nMuMax,0:nPMax)
      do 20 j=0,nmu
      do 30 k=0,nw
       f(i,j,k) = f(i,j,k) + df(i,j,k)
      df(i,j,k) = 0.
30    continue
20    continue
      return 
      end
c  *****************  end subroutine SP_update    ******************

c ==================== SOLVER ROUTINES: =========================== 
!BOP
!ROUTINE: SP_solve5 - solves A\cdot x=b equation, with 5-diagonal matrix A
!INTERFACE
      subroutine SP_solve5(nx,a2,a1,bx,c1,c2,xy)
!EOP
      include 'param.h'
      dimension a2(0:nRMax),a1(0:nRMax),bx(0:nRMax),
     1          c1(0:nRMax),c2(0:nRMax),xy(0:nRMax)
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
      subroutine SP_solvaa(nx,a2,a1,bx,xy)
      include 'param.h'
      dimension a2(0:nRMax),a1(0:nRMax),bx(0:nRMax),xy(0:nRMax)
c --  solve lower diagonal
      xy(0) = xy(0)/bx(0)
      xy(1) = (xy(1)-a1(1)*xy(0))/bx(1)
      do 11 i=2,nx
      xy(i) = (xy(i)-a1(i)*xy(i-1)-a2(i)*xy(i-2))/bx(i)
11    continue 
     
      return
      end
c  *****************  end subroutine SOLV-AA  ******************
      subroutine SP_solvcc(nx,bx,c1,c2,xy)
      include 'param.h'
      dimension bx(0:nRMax),c1(0:nRMax),c2(0:nRMax),xy(0:nRMax)
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

      subroutine SP_alla(io,jst,t)
      include 'coupler.h'
      common /SP_size  / nr,nmu,nw, dim1
      common /SP_blast/  slamb,tblast,tblst1,rshck1,dlnt
      common /SP_gazdi/  ggamma,bbrad,vvmin,vvmax,ddmin,ddmax,
     1                ccmin,ccmax,aamin,aamax,bbmin,bbmax
      common /SP_radio / nn,rmin,rshock,rmax,r(0:nRMax)
      common /SP_spiral/ tll(0:nRMax), dbdl(0:nRMax),vl(0:nRMax)
      common /SP_coeff / dvdt(0:nRMax),dldt(0:nRMax),dbdt(0:nRMax),
     1                dndt(0:nRMax)
      common /SP_scatti/ qex,Scattering_I(nMuMax),
     1     wsc(0:nPMax),xsc(0:nRMax)
      integer  iFile
      include 'stdout.h'

      call CON_io_unit_new(iFile)
      if (io.lt.10) then
         write(iStdOut,*)prefix, 'io kicsi'
         call CON_stop(' ')
      end if 
      if (io.gt.16) then
         write(iStdOut,*)prefix, 'io nagy '
         call CON_stop(' ')
      end if
      if (io.eq.10) open(iFile,file = './SP/alla0.out',status='unknown')
      if (io.eq.11) open(iFile,file = './SP/alla1.out',status='unknown')
      if (io.eq.12) open(iFile,file = './SP/alla2.out',status='unknown')
      if (io.eq.13) open(iFile,file = './SP/alla3.out',status='unknown')
      if (io.eq.14) open(iFile,file = './SP/alla4.out',status='unknown')
      if (io.eq.15) open(iFile,file = './SP/alla5.out',status='unknown')
      if (io.eq.16) open(iFile,file = './SP/alla6.out',status='unknown')

      nn=nr
      write(iFile,901) jst,t
      write(iFile,902) jst,t,rshck1,slamb
      write(iFile,*)
      write(iFile,910)      
      write(iFile,*)
      do 10 i=0,nn
      write(iFile,911) i,r(i),vr(i),algbb(i),algll(i),algnn(i)
10    continue
      write(iFile,*)
      write(iFile,920)
      write(iFile,*)
      do 20 i=0,nn
      write(iFile,921) i,r(i),dvdt(i),dldt(i),dbdt(i),dndt(i)
20    continue
      write(iFile,*)
      write(iFile,930)
      write(iFile,*)
      do 30 i=0,nn
      write(iFile,931) i,r(i),tll(i),dbdl(i),xsc(i)
30    continue
      close(iFile)
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

      subroutine SP_csilla(io,ist,t)
      include 'param.h'
      common /SP_size  / nr,nmu,nw, dim1
      common /SP_radio / nn,rmin,rshock,rmax,r(0:nRMax)
      common /SP_azimut/ fi(0:nRMax)
      common /SP_impuls/ pmin,pmax,ppin,dlnp,Momentum_I(0:nPMax)
      common /SP_energy/ emin,emax,eein,ee(0:nPMax)
      common /SP_speed / wmin,wmax,wwin,Speed_I(0:nPMax)
      common /SP_pitch / mm,Mu_I(0:nMuMax),SinMu2_I(0:nMuMax),dmu
      common /SP_solutn/ f(0:nRMax,0:nMuMax,0:nPMax),
     1     df(0:nRMax,0:nMuMax,0:nPMax)
      common /SP_peel /  pex,fpeel(0:nPMax),gpeel(0:nPMax)
      common /SP_obsrad/ krmax,krobs(5),robs(5)
      common /SP_obserg/ kemax,keobs(5),eobs(5)
      dimension ff(0:nRMax,0:nPMax),s(0:nRMax,0:nPMax),ffx(5)
      dimension ract(10)
      integer  iFile
      include 'stdout.h'
      call CON_io_unit_new(iFile)
       if (io.lt.10) then
         write(iStdOut,*)prefix, 'io kicsi'
         call CON_stop(' ')
      end if 
      if (io.gt.16) then
         write(iStdOut,*)prefix, 'io nagy '
         call CON_stop(' ')
      end if
      if (io.eq.10) open(iFile,file = 'SP/csilla0.out',status='unknown')
      if (io.eq.11) open(iFile,file = 'SP/csilla1.out',status='unknown')
      if (io.eq.12) open(iFile,file = 'SP/csilla2.out',status='unknown')
      if (io.eq.13) open(iFile,file = 'SP/csilla3.out',status='unknown')
      if (io.eq.14) open(iFile,file = 'SP/csilla4.out',status='unknown')
      if (io.eq.15) open(iFile,file = 'SP/csilla5.out',status='unknown')
      if (io.eq.16) open(iFile,file = 'SP/csilla6.out',status='unknown')
      do 11 i=0,nr
      do 13 k=0,nw
      aa = 0.5*(f(i,0,k)+f(i,nmu,k))
      ss = 0.5*(f(i,0,k)-f(i,nmu,k))
      do 12 j=1,nmu-1
      aa = aa + f(i,j,k)
      ss = ss + f(i,j,k)*Mu_I(j)
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
      write(iFile,*) 
      write(iFile,101) ist,t 
      write(iFile,102) krmax,(ract(kk),kk=1,krmax)
      write(iFile,*) 
      do 51 k = 0,nw
      do 52 kk = 1,krmax
       ii = krobs(kk)
       ffx(kk) = ff(ii,k)*fpeel(k)*Momentum_I(k)**2
52    continue
       write(iFile,111) k,ee(k),(ffx(kk),kk=1,krmax)
51    continue
         write(iFile,*) 
c --  radial dependence :
      write(iFile,*) 
      write(iFile,201) ist,t
      write(iFile,202) kemax,(eobs(kk),kk=1,kemax)
      write(iFile,*) 
      do 61 i = 0,nr
         do 62 kk = 1,kemax
            ke  = keobs(kk)
            ffx(kk) = ff(i,ke)*fpeel(ke)*Momentum_I(ke)**2
 62      continue
         write(iFile,211) i,r(i),fi(i),(ffx(kk),kk=1,kemax)
 61   continue
         write(iFile,*) 
c --  pitch - angles :
      write(iFile,*) 
      write(iFile,300) ist,t 
      do 71 kke = 1,kemax
      ke = keobs(kke)
      write(iFile,*) 
      write(iFile,301) eobs(kke)
      write(iFile,*) 
      do 72 j = 0,nmu
      do 73 kkr=1,krmax
      ii = krobs(kkr)

      if(ff(ii,ke) > 0.)then
         ffx(kkr) =  f(ii,j,ke)/ff(ii,ke)
      else
         ffx(kkr) = -1.0
      endif

73    continue
        write(iFile,311) j,Mu_I(j),(ffx(kk),kk=1,krmax)
72    continue
        write(iFile,*) 
71    continue
      close(iFile)
c --- report:
      ixx = io-10
      write(iStdout,*) prefix,
     1      ixx, ' - output done at step/time[hour] :  ',ist,t
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

c  *****************  end subroutine SP_csilla   ******************

      subroutine SP_opentime
      include 'param.h'
      common /SP_size  / nr,nmu,nw, dim1
      common /SP_radio / nn,rmin,rshock,rmax,r(0:nRMax)
      common /SP_impuls/ pmin,pmax,ppin,dlnp,Momentum_I(0:nPMax)
      common /SP_energy/ emin,emax,eein,ee(0:nPMax)
      common /SP_speed / wmin,wmax,wwin,Speed_I(0:nPMax)
      common /SP_pitch / mm,Mu_I(0:nMuMax),SinMu2_I(0:nMuMax),dmu
      common /SP_solutn/ f(0:nRMax,0:nMuMax,0:nPMax),
     1     df(0:nRMax,0:nMuMax,0:nPMax)
      common /SP_obsrad/ krmax,krobs(5),robs(5)
      common /SP_obserg/ kemax,keobs(5),eobs(5)
      common /SP_iFile/  io(5),ip(5)
      include 'stdout.h'
      
      do ke=1,kemax
         kk = keobs(ke)
         eobs(ke) = ee(kk)
         write(iStdout,*) prefix,'The #',ke,' level of energy =',eobs(ke)
      enddo
  
      do 11 kr=1,krmax
      
      call CON_io_unit_new(io(kr))
      if (kr.eq.1) open(io(kr),file='./SP/timevar.r1',status='unknown')
      if (kr.eq.2) open(io(kr),file='./SP/timevar.r2',status='unknown')
      if (kr.eq.3) open(io(kr),file='./SP/timevar.r3',status='unknown')
      if (kr.eq.4) open(io(kr),file='./SP/timevar.r4',status='unknown')
      if (kr.eq.5) open(io(kr),file='./SP/timevar.r5',status='unknown')
      write(io(kr),*)
      write(io(kr),711) robs(kr) 
      write(io(kr),712) (eobs(ke),ke=1,kemax) 
      write(io(kr),*)
      call CON_io_unit_new(ip(kr))
      if (kr.eq.1) open(ip(kr),file='./SP/plasma.r1',status='unknown')
      if (kr.eq.2) open(ip(kr),file='./SP/plasma.r2',status='unknown')
      if (kr.eq.3) open(ip(kr),file='./SP/plasma.r3',status='unknown')
      if (kr.eq.4) open(ip(kr),file='./SP/plasma.r4',status='unknown')
      if (kr.eq.5) open(ip(kr),file='./SP/plasma.r5',status='unknown')
      write(ip(kr),*)
      write(ip(kr),711) robs(kr) 
      write(ip(kr),722)  
      write(ip(kr),*)
11    continue
711   format(5x,'Time-variation at radius[AU]: ',f10.2)
712   format(5x,'Energies [MeV]:',5f10.3) 
722   format(5x,'Plasma:') 
      return
      end

c  *****************  end subroutine SP_opentime  ******************

      subroutine SP_closetime 
      common /SP_obsrad/ krmax,krobs(5),robs(5)
      common /SP_obserg/ kemax,keobs(5),eobs(5)
      common /SP_iFile/  io(5),ip(5)
      do 11 kr=1,krmax
      write(io(kr),*)
      write(io(kr),*) 'the end'
      close(io(kr))
      write(ip(kr),*)
      write(ip(kr),*) 'the end'
      close(ip(kr))
11    continue 
      return
      end

c  *****************  end subroutine SP_opentime  ******************

      subroutine SP_timevar(jst,t)
      include 'coupler.h'
      common /SP_size  / nr,nmu,nw, dim1
      common /SP_radio / nn,rmin,rshock,rmax,r(0:nRMax)
      common /SP_impuls/ pmin,pmax,ppin,dlnp,Momentum_I(0:nPMax)
      common /SP_energy/ emin,emax,eein,ee(0:nPMax)
      common /SP_speed / wmin,wmax,dlnw,Speed_I(0:nPMax)
      common /SP_pitch / mm,Mu_I(0:nMuMax),SinMu2_I(0:nMuMax),dmu
      common /SP_solutn/ f(0:nRMax,0:nMuMax,0:nPMax),
     1     df(0:nRMax,0:nMuMax,0:nPMax)
      common /SP_peel /  pex,fpeel(0:nPMax),gpeel(0:nPMax)
      common /SP_obsrad/ krmax,krobs(5),robs(5)
      common /SP_obserg/ kemax,keobs(5),eobs(5)
      common /SP_coeff / dvdt(0:nRMax),dldt(0:nRMax),dbdt(0:nRMax),
     1                dndt(0:nRMax)
      common /SP_convrt/ cAUKm,hour,valf
      dimension ffe(10)
      common /SP_iFile/  io(5),ip(5)
      include 'stdout.h'
      do 10 kr=1,krmax
      rr = robs(kr)
      ff0 = 1.0

      i=0
20    if (r(i).lt.rr) then
      if (i.eq.nr) then
         write(iStdOut,*)prefix,  'miako tyukanyo'
         call CON_stop(' ')
      end if 
      if (i.eq.nr-1) then
         robs(kr)=r(i);rr=robs(kr)   
         write(iStdOut,*)prefix,'Reset robs(',kr,')=',rr
      end if
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
      ffe(ke)=ff0*(fr2*ss1+fr1*ss2)*fpeel(kk)*Momentum_I(kk)**2
30    continue
      
      vv =  fr2*vr(i1) + fr1*vr(i2)
      vv = cAUKm/hour*vv

      abb= fr2*algbb(i1)+fr1*algbb(i2)
      all= fr2*algll(i1)+fr1*algll(i2)
      add= abb-all
      ann= fr2*algnn(i1)+fr1*algnn(i2)
      if(DoWriteAll)write(iStdout,7881)prefix,
     1      jst,t,vv,add,all,abb,ann
      write(ip(kr),788) jst,t,vv,add,all,abb,ann

      write(io(kr),777) jst,t,(ffe(ke),ke=1,kemax)
       if(DoWriteAll)write( iStdout,7771) prefix,
     1      jst,t,(ffe(ke),ke=1,kemax)
10    continue
      return
777   format(i5,f8.2,5e14.4)
788   format(i5,f8.2,5f14.6)
7771  format(a,i5,f8.2,5e14.4)
7881  format(a,i5,f8.2,5f14.6)
      end

c  *****************  end subroutine SP_timevar  ******************

