      program sep   
      implicit none
      integer iLine         !The # of the magnetic filed line to process
      integer iSnapshot     !The # of coupling 
      real TimeToRead,tSimulation
      common/SP_contime/TimeToRead,tSimulation
      real tFinal
      real tCoupleHr
      common /SP_couple/tCoupleHr
      include 'stdout.h'
      logical UseSelfSimilarity,UseRefresh
      common/SP_log/UseSelfSimilarity,UseRefresh

      call set_stdout(' :  ')      
      UseSelfSimilarity=.true.
      call SP_admit 
      call sp_init


      if(DoWriteAll)write(iStdout,*) prefix, 
     1     'Select line (1-2-3): '
CCC   if(DoWriteAll)write(iStdout,*) prefix, 'PRESS 1'
CCC   read(*,*)  lineno

 
      iLine=1
      if(DoWriteAll)write(iStdout,*) prefix, 
     1     'lineno=',iLine
CCC   if(DoWriteAll)write(iStdout,*) prefix, 'PRESS 9'
CCC   read(*,*)  iSnapshot
      iSnapshot=9            


      call sp_get_from_ih(iLine,iSnapshot)
      
!     Physical time to start the SP_update
      tSimulation=tCoupleHr*3.60d3 
!     The start time is read from the snapshot
      if(DoWriteAll)write(iStdout,*)'tSimulation=',tSimulation
      
!     Final value of physical time is not defined
      tFinal=1.0e8
      
      call SP_MASTER(tSimulation,tFinal)
      write(iStdout,*)'tSimulation=',tSimulation 
      
      call SP_closetime  !Finalize
      stop
      end 
!=============================================================!
!BOP
!ROUTINE: SP_admit - reads the input parameters
!INTERFACE:
      subroutine SP_admit
!DESCRIPTION:
!Read input parameters from list
!in case of self similar solution, sets its parameters
!sets the scattering length
!EOP
      include 'coupler.h'
      common /SP_size  / nr,nmu,nw, dim1
      common /SP_suly /  wghtl,wghtmu,wghtw
      common /SP_times/  time,tmax,dlnt0,dlnt1,dta,kfriss,kacc
      common /SP_blast/  slamb,tblast,tblst1,rshck1,dlnt
      common /SP_radio / nn,rmin,rshock,rmax,r(0:nRMax)
      common /SP_gazdi/  ggamma,bbrad,vvmin,vvmax,ddmin,ddmax,
     1                ccmin,ccmax,aamin,aamax,bbmin,bbmax

!      common /inphys/ wind0,period0,xlambda0
!      commented out: is not a common block, variables are 
!      undefined, their use is not proper 
!      I.Sokolov<igorsok@umich.edu>

      common /SP_partid/ iz,massa,ekpp,xlmbda0
!      common /scphys/ wind,omega,xscatt1 !'wind' is undefined
      common/SP_scphys/ omega,xscatt1
      common /SP_impuls/ pmin,pmax,ppin,dlnp,pp(0:nPMax)
      common /SP_energy/ emin,emax,eein,ee(0:nPMax)
      common /SP_speed / wmin,wmax,wwin,ww(0:nPMax)
      common /SP_quelle/ kinj,einj,pinj,qqp(0:nPMax),qqx(0:nRMax)
      common /SP_scatti/ qex,cmu(nMuMax),scmu(nMuMax),wsc(0:nPMax),
     1     xsc(0:nRMax)
      common /SP_obsrad/ krmax,krobs(5),robs(5)
      common /SP_obserg/ kemax,keobs(5),eobs(5)
      common /SP_convrt/ cAUKm,hour,valf
      logical UseSelfSimilarity,UseRefresh
      common/SP_log/UseSelfSimilarity,UseRefresh
      include 'stdout.h'
      integer iFile
      data  pmass,clight,xkm  / 938., 3.e10, 1.e5 / 

c ----------------------------------- scales:
      pi = 2.*asin(1.)
      cAUKm = 1.5e8
      hour = 3600.
      valf = 21.8
      e0 = pmass
c ------------------------------------ input data:
      call CON_io_unit_new(iFile)
      open(iFile,FILE='violet.in')
      read(iFile,*) nr,nmu,nw, wghtl,wghtmu,wghtw
      write(iStdout,*) prefix,
     1      'nr=',nr,'  nmu=',nmu,'  nw=',nw, 
     2 '  wght1=',wghtl,'  wghtmu=',wghtmu,
     3 '  wghtw=',wghtw
      read(iFile,*) rmin,rmax,rshock
      write(iStdout,*) prefix,
     1    'rmin=',rmin,'  rmax=',rmax,'  rshock=',rshock
      read(iFile,*) emin,emax,einj
      write(iStdout,*) prefix,
     1 'emin=',emin,'  emax=',emax,'  einj=',einj
      read(iFile,*) 
     1 tmax,slamb,kfriss,kacc
      UseRefresh=UseSelfSimilarity
      write(iStdout,*) prefix,
     1 'tmax=',tmax,'  slamb=',slamb,
     2  '  kfriss=',kfriss,'  kacc=',kacc
      read(iFile,*) swind,fwind,bbrad,period
      write(iStdout,*) prefix,
     1 'swind=',swind,'  fwind=',fwind,
     2 '  bbrad=',bbrad,'  period=',period
      read(iFile,*) iz,massa,ekpp,xlmbda0,qex
      write(iStdout,*) prefix,
     1 'iz=',iz,'  massa=',massa,'  ekpp=',
     2 ekpp,'  xlambda0=',xlmbda0,'  qex=',qex
      read(iFile,*) krmax,kemax
      write(iStdout,*) prefix,
     1 'krmax=',krmax,'  kemax=',kemax
      read(iFile,*) (robs(k),k=1,krmax)
      write(iStdout,*) prefix,
     1       ('  robs(',k,')=',robs(k),k=1,krmax)
      read(iFile,*) (keobs(k),k=1,kemax)
      close(iFile)
      write(iStdout,*) prefix, 
     1      ('  keobs(',k,')=',keobs(k),k=1,kemax)
      write(iStdout,*) prefix, 
     1       'particle energy in megavolts'
      if(UseSelfSimilarity)then
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
         bbrad =  bbrnt*valf*hour/cAUKm
         vamin =  0.
         vamax = vamin*sqrt(vvmin/vvmax)
         ddmin = 10. 
         ddmax = ddmin*vvmin/vvmax
         vvmin = vvmin*hour/cAUKm
         vvmax = vvmax*hour/cAUKm
         ccmin = vcmin*hour/cAUKm
         ccmax = vcmax*hour/cAUKm
         aamin = vamin*hour/cAUKm
         aamax = vamax*hour/cAUKm
         bbmin = aamin*sqrt(ddmin)
         bbmax = aamax*sqrt(ddmax)
c ------------------------------------ times:     
!     dt1 = dt0/kfriss  !dt0 is undefined at the moment
!     dta = dt1/kacc    !Commented out, I.Sokolov<igorsok@umich.edu>
      end if
!       wind = wind0*hour/aukm!wind0 is undefined at the moment
!      swind = swind0*hour/aukm!swind0 is undefined at the moment
c ------------------------------------ calculation
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
      rTransient=1.0
      return 
      end

c  *****************  end SP_admit   **************************
!===========================================================
      subroutine sp_get_from_ih(iLine,iSnapshot)
      implicit none
      include 'stdout.h'
      include 'coupler.h'
      real tCoupleHr
      common /SP_couple/tCoupleHr
      integer iLine
      integer iFile
      integer iLoop
      integer iLoopSnapshot
      integer iSnapshot
      integer iMisc
C ======= START COMMUNICATION WITH IH COMPONENT HERE =====================

      call CON_io_unit_new(iFile)

      if (iLine.eq.1) open(iFile,file = 'evolv1.cme' )
      if (iLine.eq.2) open(iFile,file = 'evolv2.cme' )
      if (iLine.eq.3) open(iFile,file = 'evolv3.cme' )

      read(iFile,*)
      read(iFile,*)
      read(iFile,*) iLine,iMisc  
      read(iFile,*) 
      write(*,*)DoWriteAll
      if(DoWriteAll)write(iStdout,*) prefix,
     1           'energy in megavolts'
      if(DoWriteAll)write(iStdout,*) prefix,
     1         'lineno,kmax :  ',iLine,iMisc
      if(DoWriteAll)write(iStdout,*) prefix   
      if(DoWriteAll)write(iStdout,*) prefix,
     1           'selected line: ',iLine
      if(DoWriteAll)write(iStdout,*) prefix,
     1       'Which snapshot to take (1-20): '
   
      do  iLoopSnapshot=1,iSnapshot
         read(iFile,*)  iMisc,iMax,tCoupleHr
      

         do iLoop=1,iMax 
            read(iFile,*)  iMisc,rx(iLoop),ry(iLoop),rz(iLoop)
            read(iFile,*)  iMisc,vx(iLoop),vy(iLoop),vz(iLoop)
            read(iFile,*)  iMisc,bx(iLoop),by(iLoop),bz(iLoop)
            read(iFile,*)  iMisc,dens(iLoop),pres(iLoop)
         end do

         read(iFile,*) 
      end do
      close(iFile)
      if(DoWriteAll)write(iStdout,*) prefix,
     1      'consider snapshot : ',iSnapshot
      end subroutine sp_get_from_ih
C ======= END COMMUNICATION WITH IH COMPONENT HERE =====================

