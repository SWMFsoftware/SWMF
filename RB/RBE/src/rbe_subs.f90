!******************************************************************************
!
!                               rbe_swmf.f
! 
!  This code is a modification of rbe_v02.f.  The code is created to be 
!  integrated into the Space Weather Modeling Framework.
!
!  Created on 21 March 2006 by Mei-Ching Fok, Code 612.2, NASA GSFC.
!
!  Radiation Belt-Ring Current Forecasting Model, Version 02.
!
!  Contact: Mei-Ching Fok at mei-ching.h.fok@nasa.gov, 301-286-1083.
!
!  This code is designed to model radiation belt electrons and ions with 
!  energy ranges 10 keV - 4 MeV for e-, 10 keV - 1 MeV for H+.
!
!  A program calculates temporal evolutions of radiation belt particle fluxes,
!  considering drift, losscone loss, radial, pitch-angle and energy diffusions.
!
!  For particles in the loss cone, they are assumed to decay with a
!  lifetime of half of the bounce period. 
!
!  Model limit:
!       r:    r < rb        rb is set at 12 RE
!      LT:    0 - 24 
!       M:    corresponding to 10 keV - 4 MeV for e-, 10 keV - 1 MeV for H+.
!       K:    corresponding to sine of equatorial pitch angle from 0 - 1.
!
!  Boundary conditions:
!    psd(xlati(0))=0, psd(xlati(ib+1))=fb
!    psd(w(iw1-1))=psd(w(iw1)), psd(w(iw2+1))=psd(w(iw2))   for given K, w is M 
!    psd(si(0))=psd(si(1)), psd(si(ik+1))=psd(si(ik))  for given energy, si is K
!  
!  Magnetic field model: Tsyganenko 96, 04 model (t96_01.f or t04_s.f) or MHD.
!
!  Electric field model: Weimer 2k model (w2k.f, w2k.dat) or MHD.
!
!  Plasmasphere model: Dan Ober's model (pbo_2.f)
!
!  To compile: 
!     ifortd -o rbe_swmf.out rbe_swmf.f t96_01.f t04_s.f geopack.f w2k.f pbo_2.f
!     alias ifortd  'ifort -r8 -132 -save'
! 
!  Input files: rbe_swmf.dat
!               storm.SWIMF
!               storm.symH 
!               w2k.dat
!               rbe_*.fin for initial run
!               outname_*_c.f2 for continuous run
!
!  Output files: outname_*.fls (* can be e or h for electrons or protons)
!******************************************************************************

! 13 September 2006, Mei-Ching Fok note:
! This file was modified from rbe_sub_old.f90:
! (1) In subroutine fieldpara, integration along field line uses Taylor 
!     expansion method.

!=============================================================================

subroutine rbe_init
  ! Initial setup for the rbe model

  use rbe_constant
  use rbe_cread1
  use rbe_cread2
  use rbe_cgrid
  use rbe_cfield
  use rbe_ccepara
  use rbe_convect
  use rbe_cinitial
  use rbe_cboundary
  use rbe_plasmasphere
  use rbe_time
  use ModRbTime
  use ModTimeConvert, ONLY: time_int_to_real,time_real_to_int

  !\
  ! Set the Time parameters
  !/
   if (dt > dtmax) dt=dtmax
   t=tstart-trans

  iYear = iStartTime_I(1)
!  iDOY  = get_doy(iStartTime_I(1),iStartTime_I(2),iStartTime_I(3))
!  IYD=mod(iYear,100)*1000+iDOY
  call time_int_to_real(iStartTime_I,CurrentTime)
  StartTime=CurrentTime

  !\  
  !set outname
  !/

  !set outname
  if (UseSeparatePlotFiles) then
     CurrentTime = StartTime+t
     call time_real_to_int(CurrentTime,iCurrentTime_I)
     write(outnameSep,"(i4.4,i2.2,i2.2,a,i2.2,i2.2)") & 
          iCurrentTime_I(1),iCurrentTime_I(2),iCurrentTime_I(3),'_',&
          iCurrentTime_I(4),iCurrentTime_I(5)
     write(outnameSepOrig,"(i4.4,i2.2,i2.2,a,i2.2,i2.2)") & 
          iStartTime_I(1),iStartTime_I(2),iStartTime_I(3),'_',&
          iStartTime_I(4),iStartTime_I(5)
  endif

  
  ! latitudes from Meredith as a function of L, MLT and Kp.
  call readChorusIntensity
  ! Read Horne's chorus pitch-angle and energy diffusion coeff 
  call readChorusDiffCoef
  
  if (IsStandAlone) then
     call timing_start('rbe_grids')
     call grids(re,rc,xme,xmp,q,c,js)
     call timing_stop('rbe_grids')
  endif
  call timing_start('rbe_field')
  call fieldpara(t,dt,c,q,rc,re,xlati,xmlt,&
       phi,w,si,xmass,xme,xmp)
  call timing_stop('rbe_field')
  call timing_start('rbe_initial')
  call initial(itype,ekev,xjac,ro,gride,c,xmass,d4,js,irm,&
       iba,init,il,ie)
  call timing_stop('rbe_initial')
  if (js.eq.2) call cepara(dt,ekev,Hdens,v,irm,iw1,iw2)

  call timing_start('rbe_convection')
  call convection(t,tstart,xlati,phi,&
       rc,xnsw0,vsw0,Bx0,By0,Bz0)
  call timing_stop('rbe_convection')
  call timing_start('rbe_vdrift')
  call Vdrift(re,rc,xme,dphi,xlati,ekev,potent,js,irm,iw1,iw2)
  call timing_stop('rbe_vdrift')
  call timing_start('rbe_boundary')
  call boundary(t,tstart,f2,v,xjac,xmass,p,xktd,xnd,&
     vswb0,xnswb0,itype,ibset,irm,irm0,iba)
  call timing_stop('rbe_boundary')

  ! Setup for the plasmasphere model
  par(1)=-1.0              ! disable the model's internal E field model
  par(2)=-1.0              !

  if (iplsp.eq.1) then
     call timing_start('rbe_initmain')
     call initmain(colat(ir),colat(1)) 
     call timing_stop('rbe_initmain')
     call timing_start('rbe_set_flux_tube')
     call setfluxtubevol(colat,ir,xmltd,ip,volume)
     call timing_stop('rbe_set_flux_tube')
     call timing_start('rbe_set_xy_grid')
     call setxygrid(colat,ir,xmltd,ip,xo,yo,gridoc)
     call timing_stop('rbe_set_xy_grid')
     call timing_start('rbe_init_density')
     call initdensity(itype)  ! saturated plasmasphere when itype=1
     call timing_stop('rbe_init_density')
     if (itype.eq.1) then     ! initial warm up of the plasmasphere
        call timing_start('rbe_set_pot')
        call setpot(colat,ir,xmltd,ip,potent)
        call timing_stop('rbe_set_pot')
        delt=86400.              
        call timing_start('rbe_plasmasphere')
        call plasmasphere(delt,par)
        call timing_stop('rbe_plasmasphere')
     endif
     call timing_start('rbe_get_density')
     call getdensity(colat,ir,xmltd,ip,density)
     call timing_stop('rbe_get_density')
  endif

end subroutine rbe_init
!============================================================================
subroutine rbe_run

  use rbe_constant
  use rbe_cread1
  use rbe_cread2
  use rbe_cgrid
  use rbe_cfield
  use rbe_ccepara
  use rbe_convect
  use rbe_cinitial
  use rbe_cboundary
  use rbe_plasmasphere
  use rbe_time
  use rbe_cvdrift
  use ModChorusDiffCoef
  use ModChorusIntensity
  use ModWpower
  use ModIoUnit, ONLY: UnitTmp_
  use ModRbTime
  use ModTimeConvert, ONLY: time_real_to_int
  use ModPrerunField, ONLY: UsePrerun, read_prerun, read_prerun_IE

  ! print initial fluxes
  if (t.eq.tstart .and. itype.eq.1) call rbe_save_result(.false.,.true.)

  if (t.gt.(tstart-trans).and.ires.eq.1.and.mod(t,tf).eq.0.) then
     call fieldpara(t,dt,c,q,rc,re,xlati,&
          xmlt,phi,w,si,xmass,xme,xmp)
     call boundary(t,tstart,f2,v,xjac,xmass,p,xktd,xnd,&
     vswb0,xnswb0,itype,ibset,irm,irm0,iba)
     call E_change(f2,d4,ekev,elb,eub,e_l,ecbf,iba,iw1,iw2)
     if (iplsp.eq.1) then
        call setfluxtubevol(colat,ir,xmltd,ip+1,volume)
        call setxygrid(colat,ir,xmltd,ip+1,xo,yo,gridoc)
     endif
     if (js.eq.2) call cepara(dt,ekev,Hdens,v,irm,iw1,iw2)
  endif

  if (t.gt.(tstart-trans)) then
     call convection(t,tstart,xlati,phi,&
          rc,xnsw0,vsw0,Bx0,By0,Bz0)
     call Vdrift(re,rc,xme,dphi,xlati,ekev,potent,js,irm,iw1,iw2)
     call boundary(t,tstart,f2,v,xjac,xmass,p,xktd,xnd,&
     vswb0,xnswb0,itype,ibset,irm,irm0,iba)
  endif

  if (t.ge.tstart) then
     if (idfe.eq.1.or.idfa.eq.1) call Wpower(t,density,tKp,xKph, &
          chorusI,wLshell,wmlt,Bo,ro,xmlto,irm,nKp)
     if (idfe.eq.1) call diffusee_E(f2,dt,y,bm,ekev,ro,w,xjac, &
          CHpower,ompe,cLshell,ompea,ckeV,cPA,cDEE,iw1,iw2,iba)
     if (idfa.eq.1) call diffusea_a(f2,xjac,ro,ekev,y,tya,&
          dt,CHpower,ompe,cLshell,ompea,ckeV,cPA,cDaa,iba,iw1,iw2)
     call E_change(f2,d4,ekev,elb,eub,e_l,ecce,iba,iw1,iw2)
  endif
  
  call drift(t,dt,f2,vl,vp,ro,rb,fb,dlati,dphi,iba,iw1,iw2,irm)
  call E_change(f2,d4,ekev,elb,eub,e_l,ecdt,iba,iw1,iw2)
  call losscone(f2,tcone,rmir,rc,iba,iw1,iw2)
  call E_change(f2,d4,ekev,elb,eub,e_l,eclc,iba,iw1,iw2)

  if (js.eq.2) then
     call charexchange(f2,achar,iba,iw1,iw2)  
     call E_change(f2,d4,ekev,elb,eub,e_l,ecce,iba,iw1,iw2)
  endif

  t=t+dt
  !update field if using prerun MHD field
  if (usePrerun) call read_prerun(t)
  if (usePrerun .and. iConvect==2) call read_prerun_IE(t)
  !update outname 
  CurrentTime = StartTime+t
  call time_real_to_int(CurrentTime,iCurrentTime_I)
  if (UseSeparatePlotFiles) then
     write(outnameSep,"(i4.4,i2.2,i2.2,a,i2.2,i2.2)") & 
          iCurrentTime_I(1),iCurrentTime_I(2),iCurrentTime_I(3),'_',&
          iCurrentTime_I(4),iCurrentTime_I(5)  
  endif


  if (ires.eq.1.and.mod(t,tf).eq.0.) then
     call fieldpara(t,dt,c,q,rc,re,xlati,&
          xmlt,phi,w,si,xmass,xme,xmp)
     call boundary(t,tstart,f2,v,xjac,xmass,p,xktd,xnd,&
     vswb0,xnswb0,itype,ibset,irm,irm0,iba)
     call E_change(f2,d4,ekev,elb,eub,e_l,ecbf,iba,iw1,iw2)
     if (iplsp.eq.1) then
        call setfluxtubevol(colat,ir,xmltd,ip+1,volume)
        call setxygrid(colat,ir,xmltd,ip+1,xo,yo,gridoc)
     endif
     if (js.eq.2) call cepara(dt,ekev,Hdens,v,irm,iw1,iw2)
  endif

  call convection(t,tstart,xlati,phi,&
       rc,xnsw0,vsw0,Bx0,By0,Bz0)
  call Vdrift(re,rc,xme,dphi,xlati,ekev,potent,js,irm,iw1,iw2)

  call boundary(t,tstart,f2,v,xjac,xmass,p,xktd,xnd,&
     vswb0,xnswb0,itype,ibset,irm,irm0,iba)

  if (js.eq.2) then
     call charexchange(f2,achar,iba,iw1,iw2)  
     call E_change(f2,d4,ekev,elb,eub,e_l,ecce,iba,iw1,iw2)
  endif

  call losscone(f2,tcone,rmir,rc,iba,iw1,iw2)
  call E_change(f2,d4,ekev,elb,eub,e_l,eclc,iba,iw1,iw2)

  call drift(t,dt,f2,vl,vp,ro,rb,fb,dlati,dphi,iba,iw1,iw2,irm)
  call E_change(f2,d4,ekev,elb,eub,e_l,ecdt,irm,iw1,iw2)

  if (t.ge.tstart) then
     if (idfe.eq.1.or.idfa.eq.1) call Wpower(t,density,tKp,xKph, &
          chorusI,wLshell,wmlt,Bo,ro,xmlto,irm,nKp)
     if (idfa.eq.1) call diffusea_a(f2,xjac,ro,ekev,y,tya,&
          dt,CHpower,ompe,cLshell,ompea,ckeV,cPA,cDaa,iba,iw1,iw2)
     if (idfe.eq.1) call diffusee_E(f2,dt,y,bm,ekev,ro,w,xjac, &
          CHpower,ompe,cLshell,ompea,ckeV,cPA,cDEE,iw1,iw2,iba)
     call E_change(f2,d4,ekev,elb,eub,e_l,ecce,iba,iw1,iw2)
  endif

  t=t+dt
  !update outname 
  CurrentTime = StartTime+t
  call time_real_to_int(CurrentTime,iCurrentTime_I)
  if (UseSeparatePlotFiles) then
     write(outnameSep,"(i4.4,i2.2,i2.2,a,i2.2,i2.2)") & 
          iCurrentTime_I(1),iCurrentTime_I(2),iCurrentTime_I(3),'_',&
          iCurrentTime_I(4),iCurrentTime_I(5)  
  endif
  ! update the plasmasphere density
  if (iplsp.eq.1) then
     call setpot(colat,ir,xmltd,ip+1,potent)
     delt=2.*dt
     call plasmasphere(delt,par)
     call getdensity(colat,ir,xmltd,ip+1,density)
  endif

  !  Print results
  if (t.gt.tstart.and.mod(t,tint).eq.0.) &
       call rbe_save_result(.true. .and. IsStandAlone, .true.)

  open(unit=UnitTmp_,file='RB/rbe_swmf.log')
  write(UnitTmp_,'(a8)') outname
  write(UnitTmp_,*) 'istep, t(hour)   ',istep,t/3600. 
  close(UnitTmp_)

end subroutine rbe_run

!============================================================================

subroutine rbe_save_result(DoSaveRestart, DoSavePlot)

  use rbe_time
  use rbe_cgrid
  use rbe_cread1
  use rbe_cread2
  use rbe_cfield
  use rbe_cinitial
  use rbe_convect
  use rbe_ccepara
  use rbe_cboundary
  use rbe_plasmasphere

  implicit none

  logical, intent(in) :: DoSaveRestart ! logical sets whether to save restart
  logical, intent(in) :: DoSavePlot    ! logical sets whether to save plot

  call p_result(t,tstart,f2,rc,xlati,ekev,y,p,ro,xmlto,xmlt,&
       xjac,gride,gridp,gridy,bo,xnsw,vsw,Bx,By,Bz,vswb,&
       xnswb,parmod,ecbf,ecdt,eclc,ecce,density,iprint,ntime,irm,&
       iplsp,iw1,iw2,itype,DoSaveRestart, DoSavePlot)

end subroutine rbe_save_result

! ***************************************************************************
!                             readInputData
!      read parameters: dt, tmax, species, storm,IMF and solar wind data ...,
!      and open files to read and write.
! ***************************************************************************
subroutine readInputData

  use rbe_time, ONLY: dt
  use rbe_constant
  use rbe_cread1
  use rbe_cread2
  use rbe_io_unit
  use ModIoUnit, ONLY: UnitTmp_, io_unit_new

  real:: bxw1(nswmax),byw1(nswmax),bzw1(nswmax),xnsw1(nswmax),vsw1(nswmax)
  integer :: isymH(60),Kp8(8)
  character (len=80):: header
  character  pmz(8)*1
  !---------------------------------------------------------------------------

  iprint=2                  ! 1=print result @ tmax, 2=print every tint
 
  if (IsStandAlone) then
     ntime=int((tmax-tstart)/tint)+1
  else
     ntime=2
  endif

  if (iprint.eq.1) then
     tint=tmax     
     ntime=1
  endif
  tf1=400.                       ! time range of smoothing IMF and SW data
  tf=300.                        ! time resolution in changing B
  hlosscone=100.                 ! loss cone altitude in km

  rc=(re+abs(hlosscone)*1000.)/re  ! losscone in Re
  
  if (mod(tf,dt).ne.0.) then
     write(*,*) 'RBE ERROR: mod(tf,dt).ne.0.'
     call CON_stop('RBE ERROR')
  endif

  if (itype.eq.2) trans=0.
  nstep=ifix((tmax-tstart)/dt/2.+0.5)
  nstept=ifix(trans/dt/2.+0.5)

  if (js.eq.1) st2='_e'
  if (js.eq.2) st2='_h'

  !.....Open files in case of original run and continuous run 
  iUnit1 = io_unit_new()
  open(unit=iUnit1,file='RB/rbe'//st2//'.fin',status='old')
  read(iUnit1,*) init        ! 1=NASA RB model
  read(iUnit1,*) il,ie
  if (itype.eq.2) then
     iUnit2 = io_unit_new()
     open(unit=iUnit2,file='RB/restartIN/'//outname//st2//'_c.f2',status='old',form='unformatted')
     init=0              ! continuous run
  endif

  !.....open file to read the SW-IMF data   
  if (.not. UseGm)then
     open(unit=UnitTmp_,file='RB/'//storm//'.SWIMF',status='old')
     read(UnitTmp_,*) iyear,iday,ihour        ! time corresponding to t=0
     read(UnitTmp_,*) swlag    ! time in sec for sw travel from s/c to subsolar point
     
     read(UnitTmp_,*) nsw
     read(UnitTmp_,'(a80)') header
     j=1
     do i=1,nsw
        read(UnitTmp_,*) idy,month,iyr,ihr,minu,sec,xnsw1(i),vsw1(i)        
        call modd_dayno(iyr,month,idy,iday1,j)                     
        tsw(i)=swlag+(iday1-iday)*86400.+(ihr-ihour)*3600.+minu*60.+sec 
     enddo
     do i=1,nsw                   ! smooth solar wind data
        tti=tsw(i)-tf1
        ttf=tsw(i)+tf1
        call locate1(tsw,nsw,tti,j1)
        call locate1(tsw,nsw,ttf,j_2)
        j2=j_2+1
        if (j1.eq.0) j1=1
        if (j2.gt.nsw) j2=nsw
        xnswa(i)=0.
        vswa(i)=0.
        do j=j1,j2
           xnswa(i)=xnswa(i)+xnsw1(j)/(j2-j1+1)       
           vswa(i)=vswa(i)+vsw1(j)/(j2-j1+1)
        enddo
     enddo
     
     read(UnitTmp_,*) nimf
     if (nsw.gt.nswmax.or.nimf.gt.nswmax) then
        write(6,*) 'RBE Error: nsw.gt.nswmax.or.nimf.gt.nswmax'
        call CON_stop('RBE ERROR')
     endif
     read(UnitTmp_,'(a80)') header
     j=1
     do i=1,nimf
        read(UnitTmp_,*) idy,month,iyr,ihr,minu,sec,bxw1(i),byw1(i),bzw1(i)
        call modd_dayno(iyr,month,idy,iday1,j)                     
        timf(i)=swlag+(iday1-iday)*86400.+(ihr-ihour)*3600.+minu*60.+sec 
     enddo
     do i=1,nimf                  ! smooth IMF data
        tti=timf(i)-tf1
        ttf=timf(i)+tf1
        call locate1(timf,nimf,tti,j1)
        call locate1(timf,nimf,ttf,j_2)
        j2=j_2+1
        if (j1.eq.0) j1=1
        if (j2.gt.nimf) j2=nimf
        bxw(i)=0.
        byw(i)=0.
        bzw(i)=0.
        do j=j1,j2
           bxw(i)=bxw(i)+bxw1(j)/(j2-j1+1)
           byw(i)=byw(i)+byw1(j)/(j2-j1+1)
           bzw(i)=bzw(i)+bzw1(j)/(j2-j1+1)
        enddo
     enddo
     close(UnitTmp_)
     
     ! Read Dst (symH) and Kp data
     open(unit=UnitTmp_,file='RB/'//storm//'.symHKp',status='old')
     read(UnitTmp_,*) nhour
     ndst=nhour*60                      ! 1-minute resolution         
     if (ndst.gt.ndstmax) then
        print *,'Error: ndst.gt.ndstmax'
        call CON_stop('RBE ERROR')
     endif
     j=1
     do i=1,nhour
        read(UnitTmp_,'(14x,2i2,1x,i2,13x,60i6)') month,idy,ihr,isymH
        call modd_dayno(iyr,month,idy,iday1,j)
        tdst0=(iday1-iday)*86400.+(ihr-ihour)*3600.
        do k=1,60
           m=(i-1)*60+k
           tdst(m)=tdst0+k*60.+30.
           dsth(m)=float(isymH(k))
        enddo
     enddo
     ! read Kp data
     read(UnitTmp_,*) nday
     nKp=nday*8
     j=1
     do i=1,nday
        read(UnitTmp_,'(i4,2i2,1x,8(i1,a1))') iyr,month,idy,Kp8(1),pmz(1), &
             Kp8(2),pmz(2),Kp8(3),pmz(3),Kp8(4),pmz(4),Kp8(5),pmz(5),&
             Kp8(6),pmz(6),Kp8(7),pmz(7),Kp8(8),pmz(8)
        call modd_dayno(iyr,month,idy,iday1,j)
        tKp0=(iday1-iday)*86400.
        do k=1,8
           m=(i-1)*8+k
           tKp(m)=tKp0+k*10800.-5400.
           dKp=0.
           if (pmz(k).eq.'-') dKp=-0.33
           if (pmz(k).eq.'+') dKp=0.33
           xKph(m)=float(Kp8(k))+dKp
        enddo
     enddo
     close(UnitTmp_)
  endif
end subroutine readInputData
!*****************************************************************************
!                   readChorusIntensity
!  Routine reads lower band chorus (0.1fce < f < 0.5fce) wave intensity at low
!  latitudes from Meredith as a function of L, MLT and Kp.
!*****************************************************************************
      subroutine readChorusIntensity
      use rbe_grid
      use ModIoUnit, ONLY: UnitTmp_
      use ModChorusIntensity
      parameter (irw1=70)
      real wLshell1(irw1)
      character header*80
      !------------------------------------------------------------------------
      
      nav=irw1/irw
      
      do i=1,irw1
         wLshell1(i)=(i-1)*0.1+1.05       ! L = 1.05 - 7.95
      enddo
      do j=1,ipw
         wmlt(j)=(j-1)*1.+0.5            ! mlt = 0.5 - 23.5
      enddo

 ! Calculate wLshell by averaging over wLshell
      do i=1,irw
         i1=1+(i-1)*nav
         i2=i1+nav-1
         wLshell(i)=sum(wLshell1(i1:i2))/nav
      enddo

! Read and average the wave intensity data
      open(unit=UnitTmp_,file='RB/B_wave_eq.dat',status='old')
      read(UnitTmp_,'(a80)') header
      read(UnitTmp_,'(a80)') header
      do k=1,3       ! k=1: Kp < 2,  k=2: 2 <= Kp < 4,  k=3: Kp >= 4
         do n=1,9
            read(UnitTmp_,'(a80)') header
         enddo
         do j=1,ipw
            do i=1,irw
               rB_wave=0.
               do ii=1,nav
                  read(UnitTmp_,*) rmlt,rl,rB_wave1
                  if (rB_wave1.eq.99999.0) rB_wave1=0.0
                  rB_wave=rB_wave+rB_wave1
               enddo
               chorusI(i,j,k)=rB_wave/nav     ! wave intensity in pT^2
            enddo
         enddo
      enddo
      close(UnitTmp_)
 
      return
    end subroutine readChorusIntensity


!*****************************************************************************
!                           readChorusDiffCoef
!  Routine reads Horne's chorus pitch-angle and energy diffusion coeff.
!*****************************************************************************
      subroutine readChorusDiffCoef
      use rbe_grid
      use ModIoUnit, ONLY: UnitTmp_
      use ModChorusDiffCoef
      real xx(ipa),yya(ipa),yye(ipa)
      character header*111,lhead(irc)*5,dhead(ipe)*5,ehead(iwc)*6
 
!     pitch angles from L shell = 6.5
      cPA(1:ipa)=(/ 0.,2.58258,3.55388,4.52518,5.49649,6.46779,7.43910,&
                  8.41040,9.38171,10.35301,11.32432,12.29562,13.26693,&
                  14.23823,15.20954,16.18084,17.15215,18.12345,19.09476,&
                  20.06606,21.03736,22.00867,22.97997,23.95128,24.92258,&
                  25.89389,26.86519,27.83650,28.80780,29.77911,30.75041,&
                  31.72172,32.69302,33.66433,34.63563,35.60694,36.57824,&
                  37.54955,38.52085,39.49215,40.46346,41.43476,42.40607,&
                  43.37737,44.34868,45.31998,46.29129,47.26259,48.23390,&
                  49.20520,50.17651,51.14781,52.11912,53.09042,54.06173,&
                  55.03303,56.00433,56.97564,57.94694,58.91825,59.88955,&
                  60.86086,61.83216,62.80347,63.77477,64.74608,65.71738,&
                  66.68869,67.65999,68.63130,69.60260,70.57391,71.54521,&
                  72.51652,73.48782,74.45912,75.43043,76.40173,77.37304,&
                  78.34434,79.31565,80.28695,81.25826,82.22956,83.20087,&
                  84.17217,85.14348,86.11478,87.08609,88.05739,89.02870 /)
 
      cLshell(1:irc)=(/2.5,3.5,4.5,5.5,6.5/)            ! L shell
      ompea(1:ipe)=(/1.5,2.5,5.,7.5,10./)               ! ratio of fpe/fpc
      ckeV(1:iwc)=(/10.,30.,100.,300.,1000.,3000./)     ! electron energy in keV

! File heads
      do i=1,irc
         write(lhead(i),'("_l",f3.1)') cLshell(i)
      enddo
      do i=1,ipe
         if (ompea(i).ne.5..and.ompea(i).ne.10.) &
             write(dhead(i),'("_d",f3.1)') ompea(i)
         if (ompea(i).eq.5.) write(dhead(i),'("_d",i1)') ifix(ompea(i))
         if (ompea(i).eq.10.) write(dhead(i),'("_d",i2)') ifix(ompea(i))
      enddo
      do i=1,iwc
         if (ckeV(i).lt.100.) write(ehead(i),'("_e",i2)') ifix(ckeV(i))
         if (ckeV(i).ge.100..and.ckeV(i).lt.1000.) &
             write(ehead(i),'("_e",i3)') ifix(ckeV(i))
         if (ckeV(i).ge.1000.) write(ehead(i),'("_e",i4)') ifix(ckeV(i))
      enddo
 
! Read Daa, DEE and map to the cPA grid.
      do i=1,irc
         do j=1,ipe 
            do k=1,iwc 

               ! read daa, dee
               open(unit=UnitTmp_,file='RB/Horne_chorus/fb'//lhead(i)//trim(dhead(j))//&
                    trim(ehead(k))//'.out',status='old')
               do n=1,24
                  read(UnitTmp_,'(a111)') header
               enddo
               do m=1,ipa
                  read(UnitTmp_,*) pa,daa,dae,dee,dap,dpp,rmirror
                  if (i.eq.irc) then       ! at L=6.5, no interpolation 
                     cDaa(i,j,k,m)=daa
                     cDEE(i,j,k,m)=dee
                  else
                     xx(m) =pa
                     yya(m)=daa
                     yye(m)=dee
                  endif
               enddo
               close(UnitTmp_)

               ! Interpolation of Daa, Dee for other L-shells
               if (i.lt.irc) then
                  do m=1,ipa
                     call lintp(xx,yya,ipa,cPA(m),ym)
                      cDaa(i,j,k,m)=ym   ! Daa for wave amplitude of 100 pT
                     call lintp(xx,yye,ipa,cPA(m),ym)
                     cDEE(i,j,k,m)=ym   ! DEE for wave amplitude of 100 pT
                  enddo
               endif

            enddo
         enddo
      enddo
 
      return
    end subroutine readChorusDiffCoef

!***********************************************************************
!                        grids
!                  set up all the grids
!***********************************************************************
subroutine grids(re,rc,xme,xmp,q,c,js)

  use rbe_cgrid

  real xmass1(ns),gride_e12(12),gride_e9(9),gride_i(12),xlat_data(0:ir+1),&
       y_data(12)

  data xmass1/5.4462e-4,1./  ! mass no. of e-, H+
  data y_data/0.010021,0.030708,0.062026,0.086108,0.16073,0.27682,&
       0.43083,0.60149,0.75379,0.86379,0.94890,0.98827/

  !.....Energy grid (keV)  of radiation-belt electrons, logarithmically spaced
  data (gride_e12(k),k=1,12)/10.000,17.241,29.724,51.245,88.349,152.32,&
       262.61,452.75,780.56,1345.7,2320.1,4000.0/
  data (gride_e9(k),k=1,9)/10.000,21.147,44.721,94.574,200.00,422.95,&
       894.43,1891.5,4000.0/
  !.....Energy grid (keV) of radiation-belt ions, logarithmically spaced
  data (gride_i(k),k=1,12)/10.000,15.199,23.101,35.112,53.367,81.113,&
       123.28,187.38,284.80,432.88,657.93,1000.0/
  data xlat_data/9.8403,11.809,13.777,15.742,17.705,19.665,&
       21.622,23.576,25.527,27.473,29.414,31.350,&
       33.279,35.200,37.112,39.012,40.897,42.763,&
       44.604,46.409,48.163,49.837,51.382,52.799,&
       54.100,55.293,56.387,57.392,58.313,59.158,&
       59.933,60.645,61.297,61.896,62.445,62.949,&
       63.412,63.862,64.312,64.762,65.212,65.662,&
       66.112,66.562,67.012,67.462,67.912,68.362,&
       68.812,69.262,69.712,70.162,70.612/

  pi=acos(-1.)

  !.....Set up latitude grid at the ionosphere. Dipole field is
  !     assumed in the ionosphere
  do i=1,ir
     xlati(i)=xlat_data(i)*pi/180.  
     dlati(i)=0.5*(xlat_data(i+1)-xlat_data(i-1))*pi/180
     xlati_deg=xlati(i)*180./pi
     colat(i)=90.-xlati_deg    ! colatitude in deg
  enddo

  dphi=2.*pi/ip      ! ip have to be factor of nlt
  do j=1,ip
     phi(j)=(j-1)*dphi          ! magnetic local time in radian
     xmlt(j)=phi(j)*12./pi      ! magnetic local time in hour
     xmltd(j)=xmlt(j)*15.       ! magnetic local time in degree
  end do

  do i=1,ns
     xmass(i)=xmp*xmass1(i)     ! mass of each species (kg)
  end do

  ! Set up gride and gridp
  Eo=511.*xmass(js)/xmass(1)    ! rest energy in keV
  do k=1,je
     if (js.eq.1) then
        if (je.eq.12) gride(k)=gride_e12(k)
        if (je.eq.9) gride(k)=gride_e9(k)
     endif
     if (js.eq.2) gride(k)=gride_i(k)
     gridp(k)=sqrt(gride(k)*(gride(k)+2.*Eo))*1.6e-16/c
  enddo

  !.....calculate magnetic moment (w) in joule/tesla, rw depends on iw
  b2=xme/(2.*re)**3          ! equatorial B at L = 2
  w(1)=gride(1)*q/b2         ! min. M such that ekev(2,1,1,ik) ~ gride(1)
  if (js.eq.1) rw=2.20       ! rw such that ekev(iba(1),1,iw,1) ~ gride(je)
  if (js.eq.2) rw=1.90       ! rw such that ekev(iba(1),1,iw,1) ~ gride(je)
  if (iw.eq.15) rw=rw*rw
  rw1=(rw-1.)/sqrt(rw)
  w(0)=w(1)/rw
  do k=1,iw      ! This setup makes w(k+0.5)=sqrt(w(k)*w(k+1))
     w(k)=w(k-1)*rw       
     dw(k)=w(k)*rw1                                                  
  enddo                                 !  |___.___|____.____|______.______|
  w(iw+1)=w(iw)*rw                      !             w(k)   <   dw(k+1)   >

  !.....Grids in K
  si(1)=58.
  rsi=1.45
  if (ik.eq.15) rsi=rsi*rsi
  rs1=(rsi-1.)/sqrt(rsi) ! in following sutup: si(m+0.5)=sqrt(si(m)*si(m+1))
  si(0)=si(1)/rsi
  do m=1,ik
     si(m)=si(m-1)*rsi                    
     ds(m)=si(m)*rs1                    !  |___.___|____.____|______.______|
  enddo                                 !            si(m)   <   ds(m+1)   >
  si(ik+1)=si(ik)*rsi

  !.....calculate dlati*dphi*dm*ds
  do i=1,ir
     d2=dlati(i)*dphi
     do k=1,iw
        do m=1,ik
           d4(i,k,m)=d2*dw(k)*ds(m)
        enddo
     enddo
  end do

  !.....Calculate the jacobian jac
  xjac1=4.*sqrt(2.)*pi*xmass(js)*xme/rc/re  ! factor of jacobian
  sqrtm=sqrt(xmass(js))
  do i=1,ir
     xjac2=sin(2.*xlati(i))
     do k=1,iw
        xjac(i,k)=xjac1*xjac2*sqrt(w(k))*sqrtm
     end do
  end do

  !.....Set up y grids of output f and flux
  im=12/ig
  do m=1,ig
     m1=m*im
     gridy(m)=y_data(m1)
  enddo
end subroutine grids


!***********************************************************************
!                             fieldpara
! Routine calculates kinetic energy, velocity, y, latitude and altitude
! at mirror point, etc, for given magnetic moment, K and position for a
! given magnetic field configuration.             
!***********************************************************************
subroutine fieldpara(t,dt,c,q,rc,re,xlati,xmlt,phi,w,si,&
     xmass,xme,xmp)
  use rbe_grid
  use rbe_cfield
  use rbe_cread2, ONLY: tf,dsth,tdst,byw,bzw,timf,xnswa,vswa,tsw,&
       ndst,nimf,nsw,js,iyear,iday,imod, UseEllipse
  use ModNumConst, ONLY: pi => cPi
  common/geopack/aa(10),sps,cps,bb(3),ps,cc(11),kk(2),dd(8)

  external tsyndipoleSM,MHD_B
  parameter (np=1000,nd=3)
  real xlati(ir),phi(ip),w(0:iw+1),si(0:ik+1),xmass(ns),&
       si3(np),bm1(np),rm(np),rs(np),dss(np),&
       h3(np),bs(np),bba(np),&
       x1(ir),xmlt(ip),bme(0:ir,ip,ik),xli(0:ir),&
       ra(np),dssa(np),tya3(np)
  ! coeff and integrals in Taylor expansion
  parameter (nTaylor=10)
  real a_I(0:nTaylor),b_I(0:nTaylor),sumBn(0:nTaylor),sumhBn(0:nTaylor), &
       BnI(0:nTaylor,np),hBnI(0:nTaylor,np)
  integer :: iLatTest = -1, iLonTest=-1

  integer :: imax
  real :: R_12,R_24,xmltr,xBoundary(ip),MajorAxis,MinorAxis,&
       MajorAxis2,MinorAxis2, sin2, Req2, xo1,xc, xCenter,rell2
  !--------------------------------------------------------------------------

  a_I(0)=1.
  b_I(0)=1.
  do iTaylor=1,nTaylor
     a_I(iTaylor)=a_I(iTaylor-1)*(2.*iTaylor-3.)/(2.*iTaylor)
     b_I(iTaylor)=b_I(iTaylor-1)*(2.*iTaylor-1.)/(2.*iTaylor)
  enddo

  rb=10.               ! nightside outer boundary in RE
  iopt=1               ! dummy parameter in t96_01 and t04_s 
  rlim=2.*rb
  xmltlim=2.           ! limit of field line warping in hour
  dre=0.06             ! r interval below the surface of the Earth
  n5=16                ! no. of point below the surface of the Earth

  ! Save irm0
  irm0(1:ip)=irm(1:ip)
  
  if (imod <= 2) then
     !  Determine parmod
     call TsyParmod(t,tf,tsw,xnswa,vswa,nsw,tdst,dsth,ndst,timf,byw,bzw,&
          nimf,xmp,imod,parmod)
     
     !  Call recalc to calculate the dipole tilt
     isec=mod(ifix(t),60)
     min1=ifix(t)/60
     minu=mod(min1,60)
     ihour1=ifix(t)/3600
     ihour=mod(ihour1,24)
     iday1=iday+ifix(t)/86400
     if (imod.le.2) then     
        ps=0.                ! force ps = 0 when using Tsy models
        cps=cos(ps)
        sps=sin(ps)
     else
        call recalc(iyear,iday1,ihour,min,isec)
     endif
  endif
  !  Start field line tracing.  
  call timing_start('rbe_trace')
  LONGITUDE: do j=1,ip
     irm(j)=ir
     LATITUDE: do i=1,ir
        iout=0
        xlati1=xlati(i)
        xli(i)=rc/cos(xlati1)/cos(xlati1)
        phi1=phi(j)+pi                  ! +x corresponing to noon
        
        if (imod.le.2) call tsy_trace(i,rlim,re,rc,xlati1,phi1,t,ps,parmod,&
             imod,np,npf1,dssa,bba,volume1,ro1,xmlt1,bo1,ra)
        if (imod.eq.3) call MHD_trace(xlati1,phi(j),re,i,j,np, &
             npf1,dssa,bba,volume1,ro1,xmlt1,bo1,ra)

        
        if (i==iLatTest .and. j==iLonTest) then
           write(*,*) npf1,xlati1*180.0/3.14,xmlt1,ro1
           call RB_plot_fieldline(npf1,i,j,dssa,ra,bba) 
        endif

        volume(i,j)=volume1
        ro(i,j)=ro1
        if (i.gt.1) then
           if (ro(i,j).lt.ro(i-1,j)) ro(i,j)=ro(i-1,j)
        endif
        xmlto(i,j)=xmlt1
        bo(i,j)=bo1
        phi1=xmlt1*pi/12.+pi         ! phi1=0 corresponing to noon
        
        xo(i,j)=ro1*cos(phi1)
        yo(i,j)=ro1*sin(phi1)
        gridoc(i,j)=1.
        if (npf1.eq.0) gridoc(i,j)=0.
        
        !xmlto(i,j)=xmlt_1
        !if (j.gt.1.and.xmlto(i,j).lt.0.) xmlto(i,j)=xmlto(i,j)+24.
        !ro(i,j)=ra(ieq)
        !bo(i,j)=bba(ieq)
        !xo(i,j)=xa(ieq)
        !yo(i,j)=ya(ieq)
        !if (iout.ge.1) gridoc(i,j)=0.
        !if (iout.eq.0) gridoc(i,j)=1.
        
        
        if (npf1.eq.0) then             ! open field line
           irm(j)=i-1
           exit LATITUDE                   ! do next j                 
        endif
        
        if (xmlt1.gt.xmltlim.and.xmlt1.lt.(24.-xmltlim).and.&
             abs(xmlt1-xmlt(j)).gt.xmltlim) then   ! big warping in mlt
           irm(j)=i-1
           exit LATITUDE
        endif
        
        dss2=dssa(npf1)/2.      ! find the middle point
        call locate1(dssa,npf1,dss2,im)
        im1=im
        if ((dssa(im+1)-dss2).lt.(dss2-dssa(im))) im1=im+1

        npf=n5           ! make sure B decreases to bba(im1) and rises
        dssm=0.
        do m=1,npf1
           if (m.lt.npf1) dssm=dssm+(dssa(m+1)-dssa(m))
           igood=-1
           if (m.eq.1.or.m.eq.im1) igood=1
           if (m.gt.1.and.m.lt.im1.and.bba(m).gt.bba(im1).and.&
                bba(m).lt.bm1(npf)) igood=1     ! B should be decreasing
           if (m.gt.im1.and.bba(m).gt.bba(im1).and.&
                bba(m).gt.bm1(npf)) igood=1     ! B should be increasing
           if (igood.eq.1) then
              npf=npf+1
              bm1(npf)=bba(m)
              rm(npf)=ra(m)
              if (m.lt.npf1) dss(npf)=dssm
              dssm=0.                      ! reset dssm
              if (m.eq.im1) im2=npf        ! new im1
           endif
        enddo

        do m=1,n5               ! Add n5 points below rc
           rm1=rc-m*dre
           rme=rm1*re  
           rm(n5-m+1)=rm1
           rm(npf+m)=rm1
           bm1(n5-m+1)=sqrt(4.-3.*rm1/xli(i))*xme/rme**3 !assume dipole
           bm1(npf+m)=bm1(n5-m+1)
        enddo

        npf=npf+n5
        n=npf-1  ! new no. of intervals from N to S hemisphere
        do m=1,n
           rs(m)=0.5*(rm(m+1)+rm(m))
           bs(m)=0.5*(bm1(m+1)+bm1(m))
        enddo
        do m=1,n
           if (m.le.n5.or.m.ge.(npf-n5)) then
              rl=rs(m)/xli(i)
              cost=sqrt(1.-rl*rl)
              dss(m)=dre*sqrt(3.*cost*cost+1.)/2./cost
           endif
        enddo

        !                                <--dss(m)-->
        !      rs(1)                         rs(m)                     rs(n)
        !      bs(1)                         bs(m)                     bs(n)
        !    |-------|-----|......|------|----------|----|.......|---|-------|
        !  rm(1)                       rm(m)                               rm(n+1)
        !  bm(1)                       bm(m)                               bm(n+1)

        
        ! Field line integration using Taylor expansion method
        call timing_start('rbe_taylor')
        sumBn(0:nTaylor)=0.
        sumhBn(0:nTaylor)=0.
        do m=im2-1,1,-1 !im2 = middle of field line
           ! search for the southern conjugate point
           n8=npf
           SEARCH: do ii=m,n
              if (bm1(ii+1).ge.bm1(m)) then
                 n8=ii+1
                 EXIT SEARCH
              endif
           enddo SEARCH
           n7=n8-1
           
                      
           ! field line integration at the northern hemisphere
           bs_n=1.
           do iTaylor=0,nTaylor
              bsndss=bs_n*dss(m)
              ! Sum_m=1^im2-1 (ds*B^n)
              sumBn(iTaylor)=sumBn(iTaylor)+bsndss
              sumhBn(iTaylor)=sumhBn(iTaylor)+hden(rs(m))*bsndss
              bs_n=bs_n*bs(m)
           enddo

           ! field line integration at the southern hemisphere
           m0=m+1
           if (m.lt.(im2-1)) m0=n70
           do mir=m0,n7-1
              bs_n=1.
              do iTaylor=0,nTaylor
                 bsndss=bs_n*dss(mir)
                 ! Sum_m=im2^n7-1 (ds*B^n)
                 sumBn(iTaylor)=sumBn(iTaylor)+bsndss
                 sumhBn(iTaylor)=sumhBn(iTaylor)+hden(rs(mir))*bsndss
                 bs_n=bs_n*bs(mir)
              enddo
           enddo

           ! add the partial segment near the southern conjugate point
           dssp=dss(n7)*(bm1(m)-bm1(n7))/(bm1(n8)-bm1(n7)) ! partial ds
           bsp =0.5*(bm1(m)+bm1(n7))
           bs_n=1.
           do iTaylor=0,nTaylor
              bsndss=bs_n*dssp
              ! Sum_m=1^n7-1 (ds*B^n) + correction
              BnI(iTaylor,m)=sumBn(iTaylor)+bsndss
              hBnI(iTaylor,m)=sumhBn(iTaylor)+hden(rs(n7))*bsndss
              bs_n=bs_n*bsp
           enddo

           n70=n7
        enddo

        ! Set up arrarys: si3, tya3, h3
        do m=1,im2-1
           ss=0.
           ss1=0.
           ss2=0.
           bm_n=1.
           do iTaylor=0,nTaylor
              BnIbm_n=BnI(iTaylor,m)/bm_n
              
              ! ss = int_s(1)^s(m) sqrt(B(m)-B(s)) ds
              ss=ss+a_I(iTaylor)*BnIbm_n
              ! ss1 = int_s(1)^s(m) 1/sqrt(B(m)-B(s)) ds              
              ss1=ss1+b_I(iTaylor)*BnIbm_n
              ! ss2 = int_s(1)^s(m) n_H(s)/sqrt(B(m)-B(s)) ds              
              ss2=ss2+b_I(iTaylor)*hBnI(iTaylor,m)/bm_n
              bm_n=bm_n*bm1(m)
           enddo
           si3(m)=re*sqrt(bm1(m))*ss
           tya3(m)=ss1/ro1/2.
           h3(m)=ss2/ss1
        enddo

        call timing_stop('rbe_taylor')

        si3(im2)=0.               ! equatorially mirroring
        tya3(im2)=tya3(im2-1)
        h3(im2)=hden(rm(im2))

        ! Calculate y, rmir (dist. of mirror point), T(y), bounced average [H]
        do m=0,ik+1
           sim=si(m)                 ! get Bm @ given K & location
           call lintp(si3,bm1,im2,sim,bmmx)
           if (m.ge.1.and.m.le.ik) bme(i,j,m)=bmmx
           if (i.ge.1) then
              if (m.ge.1.and.m.le.ik) bm(i,j,m)=bmmx
              y(i,j,m)=sqrt(bo(i,j)/bmmx)
              if (y(i,j,m).gt.1.) y(i,j,m)=1.
              call lintp(si3,rm,im2,sim,rmm)
              if (m.ge.1.and.m.le.ik) rmir(i,j,m)=rmm
              call lintp(si3,tya3,im2,sim,tya33)
              tya(i,j,m)=tya33
              call lintp(si3,h3,im2,sim,h33)
              if (m.ge.1.and.m.le.ik) Hdens(i,j,m)=h33  ! bounce-ave [H]
           endif
        enddo

     enddo LATITUDE                              ! end of i loop
  enddo LONGITUDE                              ! end of j loop

  ! Fill in volumes and coordinates for open (?) field lines
  do j = 1, ip
     do i=irm(j)+1,ir
        volume(i,j)=volume(irm(j),j)
        xo(i,j)=xo(irm(j),j)
        yo(i,j)=yo(irm(j),j)
     enddo
  end do

  call timing_stop('rbe_trace')

  ! Peridic boundary condition
  !do i=1,ir        
  !   volume(i,ip+1)=volume(i,1)
  !   xo(i,ip+1)=xo(i,1)
  !   yo(i,ip+1)=yo(i,1)
  !   gridoc(i,ip+1)=gridoc(i,1)
  !enddo

  !.....Calculate p, E, v, et (mc^2) at grid boundaries, and lifetime for
  !     loss cone particles
  c2mo=c*c*xmass(js)
  c4mo2=c2mo*c2mo
  do j=1,ip
     do i=1,irm(j)
        if (i.ge.1) ro2=2.*ro(i,j)*re
        do m=1,ik
           pp1=sqrt(2.*xmass(js)*bme(i,j,m))
           if (i.ge.1) tcone1=ro2*tya(i,j,m)

           do k=1,iw
              pijkm=pp1*sqrt(w(k))
              pc=pijkm*c
              c2m=sqrt(pc*pc+c4mo2)
              e=c2m-c2mo                 ! E in J
              ekev(i,j,k,m)=e/1000./q    ! E in keV
              if (i.ge.1) then
                 gamma(i,j,k,m)=c2m/c2mo
                 p(i,j,k,m)=pijkm  
                 v(i,j,k,m)=pc*c/c2m
!                 write(*,*)'!!!!tcone1,v(i,j,k,m)',tcone1,v(i,j,k,m)
                 tcone2=tcone1/v(i,j,k,m)      ! Tbounce/2
                 x=dt/tcone2
                 tcone(i,j,k,m)=0. 
                 if (x.le.80.) tcone(i,j,k,m)=exp(-x)
              endif
           enddo

        enddo
     enddo
  enddo

  ! Reduce irm by 1 for convenience in calculating vl at i=irm+0.5
  do j=1,ip
     irm(j)=irm(j)-1
  enddo

!  ! Find iba
!  do j=1,ip
!     do i=1,irm(j)
!        x1(i)=ro(i,j)
!     enddo
!     call locate1(x1,irm(j),rb,ib)
!     iba(j)=ib
!  enddo
!
!    ! Find iba
  if (UseEllipse) then
     R_24=rb                 ! boundary distance at midnight
     do j=1,ip
        imax=irm(j)
        xmltr=xmlto(imax,j)*pi/12.
        xBoundary(j)=-ro(imax,j)*cos(xmltr)
     enddo
     R_12=0.95*maxval(xBoundary)    ! boundary distance at noon
     MajorAxis=0.5*(R_12+R_24)      ! major axis
     if (R_12 < R_24) then
        MinorAxis = R_12       ! minor axis
     else
        MinorAxis = R_24       ! minor axis
     endif
     xCenter=0.5*(R_12-R_24)        ! center on x-axis
     MajorAxis2=MajorAxis*MajorAxis
     MinorAxis2=MinorAxis*MinorAxis
     do j=1,ip
        find_ib: do i=irm(j),1,-1
           xmltr=xmlto(i,j)*pi/12.
           sin2=sin(xmltr)*sin(xmltr)
           Req2=ro(i,j)*ro(i,j)
           xo1=-ro(i,j)*cos(xmltr)
           xc=xo1-xCenter
           rell2=MinorAxis2*(1.-xc*xc/MajorAxis2)/sin2 ! r^2 onellipse at x=xc
           if (Req2.le.rell2) then
              iba(j)=i
              exit find_ib
           endif
        enddo find_ib
     enddo
  else
     !use circle
     do j=1,ip
        do i=1,irm(j)
           x1(i)=ro(i,j)
        enddo
        call locate1(x1,irm(j),rb,ib)
        iba(j)=ib
     enddo
  endif


end subroutine fieldpara


!*******************************************************************************
!                             TsyParmod
!  Rountine calculates the parmod in Tsyganenko model.
!*******************************************************************************
subroutine TsyParmod(t,tf,tsw,xnswa,vswa,nsw,tdst,dsth,ndst,timf,&
     byw,bzw,nimf,xmp,imod,parmod)
  use rbe_cread2, ONLY: ismo
  real tsw(nsw),xnswa(nsw),vswa(nsw),tdst(ndst),dsth(ndst),timf(nimf),&
       byw(nimf),bzw(nimf),parmod(10),w04(6),rr(6),xlamb(6),beta(6),gamm(6)

  ! Parameters for T04_S model
  data rr/0.39,0.7,0.031,0.58,1.15,0.88/     ! relaxation rate in hour^-1
  data xlamb/0.39,0.46,0.39,0.42,0.41,1.29/
  data beta/0.8,0.18,2.32,1.25,1.6,2.4/
  data gamm/0.87,0.67,1.32,1.29,0.69,0.53/

  parmod(1:10)=0.             ! initial values

  tsmo=1.7                    ! running window time normalized by tf
  if (ismo.eq.0) tsmo=0.
  tti=t-tsmo*tf
  ttf=t+(1.+tsmo)*tf

  !  parmod(1): solar wind pressure in nPa
  call locate1(tsw,nsw,tti,j1)
  call locate1(tsw,nsw,ttf,j_2)
  j2=j_2+1
  if (j1.eq.0) j1=1
  if (j2.gt.nsw) j2=nsw
  xnsw=0.
  vsw=0.
  do j=j1,j2
     xnsw=xnsw+xnswa(j)/(j2-j1+1)        ! find average SW ./. tti and ttf
     vsw=vsw+vswa(j)/(j2-j1+1)
  enddo
  parmod(1)=xmp*xnsw*1.e6*vsw*vsw*1.e6/1.e-9    ! Pdyn in nPa
  if (parmod(1).lt.2.0) parmod(1)=2.0  ! set min parmod(1) to 2.0

  !  parmod(2): Dst
  call locate1(tdst,ndst,tti,j1)
  call locate1(tdst,ndst,ttf,j_2)
  j2=j_2+1
  if (j1.eq.0) j1=1
  if (j2.gt.ndst) j2=ndst
  dst=0.
  do j=j1,j2
     dst=dst+dsth(j)/(j2-j1+1)        ! find average dst ./. tti and ttf
  enddo
  parmod(2)=dst

  !  parmod(3:4): IMF By, Bz in nT
  call locate1(timf,nimf,tti,j1)
  call locate1(timf,nimf,ttf,j_2)
  j2=j_2+1
  if (j1.eq.0) j1=1
  if (j2.gt.nimf) j2=nimf
  byimf=0.
  bzimf=0.
  do j=j1,j2
     byimf=byimf+byw(j)/(j2-j1+1)    ! find average IMF ./. tti and ttf
     bzimf=bzimf+bzw(j)/(j2-j1+1)
  enddo
  parmod(3)=byimf
  parmod(4)=bzimf

  !  Limit the values of parmod(1:4) in t96 model (imod=1)
  if (imod.eq.1) then
     if (parmod(1).gt.10.) parmod(1)=10.              
     if (parmod(2).lt.-100.) parmod(2)=-100.
     if (parmod(2).gt.20.) parmod(2)=20.
     if (parmod(3).lt.-10.) parmod(3)=-10.
     if (parmod(3).gt.10.) parmod(3)=10.
     if (parmod(4).lt.-10.) parmod(4)=-10.
     if (parmod(4).gt.10.) parmod(4)=10.
  endif

  !  parmod(5:10) for t04_s: w04(1:6) defined in Tsyganenko and Sitnov, 2005
  if (imod.eq.2) then
     tti=t-100.*3600.              ! 100 hours before t
     call locate1(tsw,nsw,tti,j1)    
     if (j1.eq.0) j1=1
     call locate1(tsw,nsw,t,j2)
     if (j2.eq.0) j2=1
     w04(1:6)=0.
     do j=j1,j2      ! average over preceding hours
        tk=tsw(j)
        tdiff=(tk-t)/3600.   ! time difference in hour
        if (tk.lt.timf(1)) tk=timf(1)
        if (tk.gt.timf(nimf)) tk=timf(nimf)        
        call lintp(timf,bzw,nimf,tk,bz1)
        if (bz1.lt.0.) Bs1=-bz1
        if (bz1.ge.0.) goto 1       ! +ve Bz, no contribution to w04
        xnsw_n=xnswa(j)/5.          ! normalized sw density
        vsw_n=vswa(j)/400.          ! normalized sw velocity
        Bs_n=Bs1/5.                 ! normalized Bs
        do m=1,6
           ert=exp(rr(m)*tdiff)
           Sk=xnsw_n**xlamb(m)*vsw_n**beta(m)*Bs_n**gamm(m)
           w04(m)=w04(m)+Sk*ert
        enddo
1       continue
     enddo
     del_t=(tsw(j2)-tsw(j1))/(j2-j1+1)/3600.      ! delta t in hour
     if (del_t.le.0.) del_t=1./12.                
     do m=1,6
        w04(m)=w04(m)*rr(m)*del_t
        parmod(m+4)=w04(m)
     enddo
  endif        ! end of if (imod.eq.2) 

  ! Set limit to parmod for t04_s
  if (imod.eq.2) then
     if (parmod(2).lt.-300.) parmod(2)=-300.     ! Dst
     if (parmod(8).gt.25.) parmod(8)=25.         ! partial ring current
     if (parmod(10).gt.100) parmod(10)=100.      ! region 2 current
  endif
end subroutine TsyParmod


!***********************************************************************
!                         cepara
! routine calculates the cross-section of charge exchange of energetic
! H+ with the cold neutral hydrogen and the charge exchange
! rate
!***********************************************************************
subroutine cepara(dt,ekev,Hdens,v,irm,iw1,iw2)
  use rbe_ccepara

  real ekev(ir,ip,iw,ik),v(ir,ip,iw,ik),Hdens(ir,ip,ik)
  integer iw1(ik),iw2(ik),irm(ip)

  !.....Chamberlain model of [H] fitted by an exponential function (hden) 
  !     is used. The fit matches Rairden et al. [1986] for radial distance 
  !     from 1.08 to 12 Re.
  !     Calculate charge exchange crosss-section of ring species js with H
  !     and then the charge exchange decay rate achar

  do j=1,ip
     do i=1,irm(j)
        do m=1,ik
           do k=iw1(m),iw2(m)
              x=log10(ekev(i,j,k,m))
              if(x.lt.-2.) x=-2.
              d=-18.767-0.11017*x-3.8173e-2*x**2-0.1232*x**3-5.0488e-2*x**4
              sigma=10.**d           ! cross section of h+ in m2
              alpha=v(i,j,k,m)*sigma*Hdens(i,j,m)*dt
              achar(i,j,k,m)=exp(-alpha) ! charge. exchange decay rate
           enddo
        enddo
     enddo
  enddo

end subroutine cepara


!!****************************************************************************
!!                            diffusee
!!  Routine calculates the change of electron distributions due to
!!  diffusion in E
!!
!!  The expression for DEE is from Summers and Ma, 2000 (JGR, p2625).
!!****************************************************************************
!subroutine diffusee(f2,dt,ekev,bo,ro,w,dw,density,bm,q,xmp,&
!     xmass,c,iw1,iw2,irm)
!  use rbe_grid
!  real bo(ir,ip),density(ir,ip+1),w(0:iw+1),dw(iw),xmass(ns),&
!       bm(ir,ip,ik),f2(ir,ip,iw,ik),fl(0:iw+1),&
!       um(iw),up(iw),um_1(iw),up1(iw),&
!       fr(iw),a1d(iw),b1d(iw),c1d(iw),ekev(0:ir,ip,iw,ik),ro(ir,ip)
!  integer iw1(ik),iw2(ik),irm(ip)
!
!  pi=acos(-1.)
!  e_mass=xmass(1)          ! electron mass in kg
!  r_mass=xmass(2)/xmass(1) ! mass ratio of proton to electron
!  epsilon0=8.8542e-12      ! permittivity of free space
!  xmuo=4.*pi*1.e-7         ! permeability of free space
!  Eo=511.                  ! electron rest energy in keV
!  xlam=0.5                 ! implicitness in solving diffusion equation
!  alam=1.-xlam
!  waveq=5./3.              ! turbulence spectral index
!  deltB=10.0e-11           ! wave amplitude in Tesla
!  r_up=8.                  ! r bound of diffusion in E
!  den_min=1.e6             ! min density 1/m^3
!
!  do j=1,ip
!     do i=1,irm(j)
!
!        den_loc=density(i,j)        ! the value of the density in SI unit
!        if (den_loc.gt.den_min.and.ro(i,j).le.r_up) then 
!           b0_loc=bo(i,j)
!           dbbe=b0_loc*deltB/e_mass      ! D0 in s^-1, from Eq.14 of Summers
!           D0=pi*(waveq-1.)/8.*q*dbbe*dbbe*(epsilon0/den_loc)**1.5/sqrt(xmp)
!           Va=b0_loc/sqrt(xmuo*den_loc*xmass(2))   ! Alfven speed
!           beta_a=Va/c
!           Ec=sqrt(1.+r_mass*beta_a*beta_a)-1.   ! critical energy, Eq.20
!           Ec_keV=Ec*Eo   ! critical energy in keV, resonance only at E > Ec
!           do m=1,ik
!              k1=iw2(m)
!              do k=iw1(m),iw2(m)
!                 if (Ec_keV.le.ekev(i,j,k,m)) then
!                    k1=k
!                    goto 1
!                 endif
!              enddo
!1             if (k1.eq.iw2(m)) then
!                 write(6,*) 'Error: k1.eq.iw2(m)'
!                 call CON_stop('RBE ERROR')
!              endif
!              k2=iw2(m)
!              iww=k2-k1+1
!              bm_loc=bm(i,j,m)           ! magnetic field at mirror point
!              Wo=Eo*1.6e-16/bm_loc       ! magnetic moment in Joule/Tesla
!              Wo2=Wo*Wo
!
!              ! Find DMM, up and um
!              u_mx=0.
!              do k=k1,k2
!                 factor_1=dw(k)*(w(k)-w(k-1))/Wo2   ! normalized fact.
!                 factor1=dw(k)*(w(k+1)-w(k))/Wo2    !
!                 WM_1=sqrt(w(k)*w(k-1))             ! M@lower grid
!                 WM1=sqrt(w(k)*w(k+1))              ! M@upper grid
!                 gjacobiM_1=sqrt(WM_1)
!                 gjacobiM1=sqrt(WM1)
!
!                 if (k.gt.1) e_1=sqrt(ekev(i,j,k-1,m)*ekev(i,j,k,m))
!                 if (k.eq.1) e_1=ekev(i,j,k,m)*sqrt(ekev(i,j,k,m)/&
!                      ekev(i,j,k+1,m))
!                 if (k.lt.iw) e1=sqrt(ekev(i,j,k,m)*ekev(i,j,k+1,m))
!                 if (k.eq.iw) e1=ekev(i,j,k,m)*sqrt(ekev(i,j,k,m)/&
!                      ekev(i,j,k-1,m))
!                 e_em=e_1/Eo             ! normalized kinetic energy
!                 e_ep=e1/Eo              ! normalized kinetic energy
!                 DEEm=D0*sqrt(e_em*(e_em+2.))/(e_em+1.)  ! Eq.6 from Summers
!                 DEEp=D0*sqrt(e_ep*(e_ep+2.))/(e_ep+1.)  ! Eq.6 from Summers
!                 DMMm=DEEm*(e_em+1.)*(e_em+1.)           ! normalized DMM
!                 DMMp=DEEp*(e_ep+1.)*(e_ep+1.)           ! normalized DMM
!                 um(k)=dt*DMMm*gjacobiM_1/factor_1/sqrt(w(k))
!                 um_1(k)=dt*DMMm*gjacobiM_1/factor_1/sqrt(w(k-1))
!                 up(k)=dt*DMMp*gjacobiM1/factor1/sqrt(w(k))
!                 up1(k)=dt*DMMp*gjacobiM1/factor1/sqrt(w(k+1))
!                 ump_mx=max(abs(up1(k)),abs(um_1(k)))
!                 if (ump_mx.gt.u_mx) u_mx=ump_mx
!              enddo             ! end k loop
!
!              ! reduce time step size if up or um is too large
!              irun=ifix(u_mx)+1
!              do k=k1,k2
!                 um(k)=um(k)/irun
!                 um_1(k)=um_1(k)/irun
!                 up1(k)=up1(k)/irun
!                 up(k)=up(k)/irun
!                 a1d(k)=-xlam*um_1(k)
!                 b1d(k)=1.+xlam*(um(k)+up(k))
!                 c1d(k)=-xlam*up1(k)
!              enddo
!              ! Start diffusion in M
!              do k=k1,k2
!                 fl(k)=f2(i,j,k,m)
!              enddo
!              if (k1.gt.iw1(m)) fl(k1-1)=f2(i,j,k1-1,m)  ! b. c. at low E
!              do nrun=1,irun
!                 if (k1.eq.iw1(m)) fl(k1-1)=fl(k1)*sqrt(w(k1-1)/w(k1)) 
!                 fl(k2+1)=fl(k2)*sqrt(w(k2+1)/w(k2))     ! b. c. at high E
!                 fr(k1)=alam*um_1(k1)*fl(k1-1)+(1.-alam*(up(k1)+um(k1)))*&
!                      fl(k1)+alam*up1(k1)*fl(k1+1)+xlam*um_1(k1)*fl(k1-1)
!                 do k=k1+1,k2-1    ! calculate the RHS of matrix equation
!                    fr(k)=alam*um_1(k)*fl(k-1)+(1.-alam*(up(k)+um(k)))*&
!                         fl(k)+alam*up1(k)*fl(k+1)
!                 enddo
!                 fr(k2)=alam*um_1(k2)*fl(k2-1)+(1.-alam*(up(k2)+um(k2)))*&
!                      fl(k2)+alam*up1(k2)*fl(k2+1)+xlam*up1(k2)*fl(k2+1)
!                 call tridag(a1d(k1),b1d(k1),c1d(k1),fr(k1),fl(k1),iww,ier) 
!              enddo
!
!              ! get back f2
!              do k=k1,k2
!                 f2(i,j,k,m)=fl(k)
!              enddo
!
!           enddo       ! end of m=1,ik
!        endif          ! end of if (den_loc.gt.den_min.and.ro(i,j).le.r_up)
!
!     enddo       ! end of i loop
!  enddo          ! end of j loop
!
!end subroutine diffusee


!!****************************************************************************
!!                             diffuesa_a
!! Routine solves pitch-angle diffusion in ao, the equatorial pitch angle.
!!****************************************************************************
!subroutine diffusea_a(f2,xjac,ro,ekev,y,tya,tanA2,dt,irm,iw1,iw2)
!  use rbe_grid
!  parameter (ie=40) 
!  real f2(ir,ip,iw,ik),xjac(ir,iw),ro(ir,ip),ekev(0:ir,ip,iw,ik),&
!       y(ir,ip,0:ik+1),tya(ir,ip,0:ik+1),tanA2(ir,ip,0:ik+1),ein_log(ie),&
!       f1d(iw),e1d(iw),ein(ie),ao(0:ik+1),dao(ik),Gjac(0:ik+1),DD(0:ik+1),&
!       um(ik),um_1(ik),up(ik),up1(ik),a1d(ik),b1d(ik),c1d(ik),f0(0:ik+1),&
!       f2d(ie,ik),fcon(ie,ik),df(ie,ik),df1(ie),fr(ik),ekevlog(iw,ik) 
!  integer irm(ip),iw1(ik),iw2(ik)
!
!  u_max=100.         ! maximum value of mu 
!  u_max_log=log10(u_max)   
!  xl1=2.             ! lower limit of diffusion region
!  xl2=8.             ! upper limit of diffusion region
!  daa=1.0e-3         ! pitch-angle diffusion coeff. in s^-1
!
!  do j=1,ip
!     do i=1,irm(j)
!        if (ro(i,j).le.xl2.and.ro(i,j).ge.xl1) then
!
!           ! Set up the energy grid, ein
!           emin=1.e20
!           do m=1,ik
!              emin=min(emin,ekev(i,j,iw1(m),m))
!           enddo
!           emax=0.
!           do m=1,ik
!              emax=max(emax,ekev(i,j,iw2(m),m))
!           enddo
!           ein(1)=emin
!           ein(ie)=emax
!           ein_log(1)=log10(emin)
!           ein_log(ie)=log10(emax)
!           x0=ein_log(1)
!           x2=ein_log(ie)
!           x=(x2-x0)/(ie-1)
!           do k=2,ie-1
!              ein_log(k)=x0+(k-1)*x
!              ein(k)=10.**ein_log(k)
!           enddo
!
!           ! Map psd to ein grid, f2d 
!           do m=1,ik
!              do k=1,iw
!                 f1d(k)=-50.                      ! f1d is log(psd)
!                 if (f2(i,j,k,m).gt.0.) f1d(k)=log10(f2(i,j,k,m)/xjac(i,k))
!                 ekevlog(k,m)=log10(ekev(i,j,k,m))     
!                 e1d(k)=ekevlog(k,m)              ! e1d is log(ekev)
!              enddo
!              do k=1,iw
!                 if (k.lt.iw1(m)) f1d(k)=f1d(iw1(m))
!                 if (k.gt.iw2(m)) f1d(k)=f1d(iw2(m))
!              enddo
!              do k=1,ie
!                 if (ein_log(k).ge.e1d(1).and.ein_log(k).le.e1d(iw)) then 
!                    call lintp(e1d,f1d,iw,ein_log(k),x)
!                    f2d(k,m)=10.**x      ! f2d is psd
!                 endif
!              enddo
!           enddo
!
!           ! calcuate ao, dao, and Gjac
!           do m=0,ik+1
!              ao(m)=asin(y(i,j,m))
!              Gjac(m)=tya(i,j,m)*y(i,j,m)*sqrt(1.-y(i,j,m)*y(i,j,m))
!           enddo
!           do m=1,ik
!              dao(m)=0.5*(ao(m+1)-ao(m-1))
!           enddo
!
!           ! Calculate conservative psd in (ein,ao), fcon
!           do m=1,ik
!              do k=1,ie
!                 if (ein_log(k).lt.ekevlog(1,m)) f2d(k,m)=f2d(k,m-1)
!              enddo
!           enddo
!           do m=ik,1,-1
!              do k=1,ie
!                 if (ein_log(k).gt.ekevlog(iw,m)) f2d(k,m)=f2d(k,m+1)
!                 fcon(k,m)=f2d(k,m)*Gjac(m)
!              enddo
!           enddo
!
!           do k=1,ie      ! *****
!
!              ! calculate DD, Daoao*Gjac
!              do m=0,ik+1
!                 Daoao=tan(ao(m))*tan(ao(m))*daa*tanA2(i,j,m)
!                 DD(m)=Daoao*Gjac(m)
!              enddo
!
!              ! calculate up and um
!              u_mx=0.
!              do m=1,ik
!                 factor_1=dao(m)*(ao(m)-ao(m-1))
!                 factor1=dao(m)*(ao(m+1)-ao(m))
!                 DDm=0.5*(DD(m)+DD(m-1))
!                 DDp=0.5*(DD(m)+DD(m+1))
!                 um(m)=dt*DDm/factor_1/Gjac(m)
!                 um_1(m)=dt*DDm/factor_1/Gjac(m-1)
!                 up(m)=dt*DDp/factor1/Gjac(m)
!                 up1(m)=dt*DDp/factor1/Gjac(m+1)
!                 ump_mx=max(abs(up1(m)),abs(um_1(m)))
!                 if (ump_mx.gt.u_mx) u_mx=ump_mx
!              enddo
!
!              ! reduce time step if u_mx > u_max
!              irun=ifix(u_mx/u_max)+1
!              do m=1,ik
!                 um(m)=um(m)/irun
!                 um_1(m)=um_1(m)/irun
!                 up(m)=up(m)/irun
!                 up1(m)=up1(m)/irun
!                 a1d(m)=-xlam*um_1(m)
!                 b1d(m)=1.+xlam*(um(m)+up(m))
!                 c1d(m)=-xlam*up1(m)
!              enddo
!
!              ! determine the implicitness, xlam
!              if (u_mx.le.1.) then
!                 xlam=0.5              ! Crank-Nicolson
!              else
!                 u_mx_log=log10(u_mx)
!                 xlam=0.5*u_mx_log/u_max_log+0.5
!              endif
!              alam=1.-xlam
!
!              ! start diffusion in ao
!              do m=1,ik
!                 f0(m)=fcon(k,m)
!              enddo
!              do n=1,irun
!                 f0(0)=f0(1)*Gjac(0)/Gjac(1)             ! psd(0)=psd(1)
!                 f0(ik+1)=f0(ik)*Gjac(ik+1)/Gjac(ik)     ! psd(ik+1)=psd(ik)
!                 fr(1)=alam*um_1(1)*f0(0)+(1.-alam*(up(1)+um(1)))*f0(1)&
!                      +alam*up1(1)*f0(2)+xlam*um_1(1)*f0(0)
!                 do m=2,ik-1    ! calculate the RHS of matrix equation
!                    fr(m)=alam*um_1(m)*f0(m-1)+(1.-alam*(up(m)+um(m)))*f0(m)&
!                         +alam*up1(m)*f0(m+1)
!                 enddo
!                 fr(ik)=alam*um_1(ik)*f0(ik-1)+(1.-alam*(up(ik)+um(ik)))*f0(ik)&
!                      +alam*up1(ik)*f0(ik+1)+xlam*up1(ik)*f0(ik+1)
!                 call tridag(a1d,b1d,c1d,fr,f0(1),ik,ier)
!              enddo
!              do m=1,ik
!                 df(k,m)=f0(m)/Gjac(m)-f2d(k,m)    ! df is differential psd
!              enddo
!
!           enddo         ! end of do k=1,ie      ! *****
!
!           ! map psd back to M grid
!           do m=1,ik
!              do k=1,ie
!                 df1(k)=df(k,m)             
!              enddo
!              do k=iw1(m),iw2(m)
!                 call lintp(ein_log,df1,ie,ekevlog(k,m),dpsd)
!                 f2(i,j,k,m)=f2(i,j,k,m)+xjac(i,k)*dpsd
!                 if (f2(i,j,k,m).lt.0.) f2(i,j,k,m)=0.
!              enddo
!           enddo
!
!        endif       ! end of if (ro(i,j).le.xl2.and.ro(i,j).ge.xl1)
!     enddo
!  enddo
!
!end subroutine diffusea_a


!****************************************************************************
!                                convection
!  Routine calculates the velocities of convection
!****************************************************************************
subroutine convection(t,tstart,xlati,phi,&
     rc,xnsw0,vsw0,Bx0,By0,Bz0)

  use EIE_ModWeimer, ONLY: setmodel00, boundarylat00, epotval00
  use ModNumConst, ONLY: pi => cPi
  use rbe_convect
  use rbe_cread2,ONLY:nimf,timf,&
       bxw,byw,bzw,nsw,tsw,xnswa,vswa,iconvect,itype,UseGm,UseIe
  
  real xlati(ir),phi(ip)
  logical UseAL

  COMMON /GEOPACK/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,&
       CPS,SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,A33,DS3,&
       K,IY,CGST,SGST,BA(6)

  UseAL=.false.
  ALindex=10.          ! arbitrary value 

  !  Setup for Weimer's electric field model if iconvect = 1
  if (iconvect.eq.1) then   
     Tilt=psi*180./pi          ! dipole tilt angle in degree
     if (t.eq.tstart.and.itype.eq.2) then
        xnsw=xnsw0
        vsw=vsw0
        bx=Bx0
        by=By0
        bz=Bz0
     else
        if (UseGm)then
           xnsw=xnswa(1)
           vsw=vswa(1)
           bx=bxw(1)
           by=byw(1)
           bz=bzw(1)
           t1=t
        else
           t1=t
           if (t1.lt.tsw(1)) t1=tsw(1)
           if (t1.gt.tsw(nsw)) t1=tsw(nsw)
           call lintp(tsw,xnswa(1:nsw),nsw,t1,xnsw)
           call lintp(tsw,vswa(1:nsw),nsw,t1,vsw)
           t1=t
           if (t1.lt.timf(1)) t1=timf(1)
           if (t1.gt.timf(nimf)) t1=timf(nimf)
           call lintp(timf,bxw(1:nimf),nimf,t1,bx)
           call lintp(timf,byw(1:nimf),nimf,t1,by)
           call lintp(timf,bzw(1:nimf),nimf,t1,bz)
        endif
     endif
     angle=atan2(by,bz)*180./pi       ! degrees from northward toward +Y
     Bt=sqrt(by*by+bz*bz)             ! Magnitude of IMF in Y-Z plane in nT
     call SetModel00(angle,Bt,Tilt,vsw,xnsw,ALindex,UseAL)
  endif

  !  Find potential (in Volt) at the ionosphere
  if (.not. UseIE) then
     do i=1,ir
        gLAT=acos(cos(xlati(i))/sqrt(rc))*180./pi !invariant lat. in degree
        do j=1,ip
           gMLT=phi(j)*12./pi    ! convert mlt from radians to hour 0-24.
           BnLat=BoundaryLat00(gMLT)    ! boundary latitude 
           if (iconvect.eq.1) then
              potent(i,j)=0.0
              if (gLAT.gt.BnLat) potent(i,j)=EpotVal00(gLAT,gMLT)*1000.  ! Volt
           else
              !              potent(i,j)=MHD_pot(gLAT,gMLT)
           endif
        enddo
     enddo
  end if
!  potent(0:ir+1,ip+1)=potent(0:ir+1,1)      ! periodic boundary condition

end subroutine convection


!*****************************************************************************
!                               Vdrift
!  Routine calculates drift velocity at grid faces: magnetic drifts +
!  electric drift + corotation.
! 
!  Note: This routine was created in February 2006. The algorithm was 
!        suggested by Daniel Spicer. Combining Vdrift with driftlp greatly
!        reduces the numerical diffusion for high-energy particles.
!*****************************************************************************
subroutine Vdrift(re,rc,xme,dphi,xlati,ekev,potent,js,irm,iw1,iw2)
  use rbe_cvdrift
  use ModNumConst, ONLY: pi => cPi
  
  real xlati(ir),ekev(ir,ip,iw,ik),potent(ir,ip),ham(ir,ip), kfactor
  integer iw1(ik),iw2(ik),irm(ip)

  !--------------------------------------------------------------------------

  dphi2=dphi*2.
  kfactor=xme/rc/re
  cor=2.*pi/86400.                  ! corotation speed in rad/s  
  if (js.eq.1) icharge=-1           ! electrons 
  if (js.eq.2) icharge=1            ! ions 

  do m=1,ik
     do k=iw1(m),iw2(m)                 

        ! ham: Hamiltonian/q
        do j=1,ip
           do i=1,irm(j)+1
              ham(i,j)=icharge*ekev(i,j,k,m)*1000.+potent(i,j)
           enddo
        enddo
        
        ! calculate drift velocities vl and vp
        do i=1,ir-1
           i_1=i-1
           if (i_1.lt.1) i_1=1
           ksai=kfactor*sin(2.*xlati(i))
           xlati1=0.5*(xlati(i)+xlati(i+1))
           ksai1=kfactor*sin(2.*xlati1)         ! ksai at i+0.5
           dlat2=xlati(i+1)-xlati(i_1)
           do j=1,ip
              j0=j-1
              if (j0.lt.1) j0=j0+ip
              j2=j+1
              if (j2.gt.ip) j2=j2-ip
              
              ! calculate vl
              if (irm(j0).ge.i.and.irm(j2).ge.i) then
                 sf0=0.5*(ham(i+1,j0)+ham(i,j0))
                 sf2=0.5*(ham(i+1,j2)+ham(i,j2))
                 vl(i,j,k,m)=-(sf2-sf0)/dphi2/ksai1     ! vl at (i+0.5,j)
              else
                 vl(i,j,k,m)=vl(i-1,j,k,m)
              endif
              
              ! calculate vp
              if (irm(j2).ge.i.and.irm(j).ge.i) then
                 sf0=0.5*(ham(i_1,j2)+ham(i_1,j))
                 sf2=0.5*(ham(i+1,j2)+ham(i+1,j))
                 vp(i,j,k,m)=cor+(sf2-sf0)/dlat2/ksai   ! vp@(i,j+0.5)
              else
                 vp(i,j,k,m)=vp(i-1,j,k,m)  
              endif
           end do      ! end of j loop
        end do         ! end of i loop
     end do            ! end of k loop
  end do               ! end of m loop            

  
end subroutine Vdrift


!******************************************************************************
!                           initial
! Initially set up the distribution function (f2) 
!******************************************************************************
subroutine initial(itype,ekev,xjac,ro,gride,c,xmass,d4,js,irm,&
     iba,init,il,ie)

  use rbe_cinitial
  use rbe_io_unit, ONLY: iUnit1, iUnit2
  use ModIoUnit, ONLY: UnitTmp_
  use ModNumConst, ONLY: pi => cPi
  use rbe_cread1

  real:: ekev(ir,ip,iw,ik),xjac(ir,iw),roi(il),ro(ir,ip),xmass(ns),&
       ei_MeV(ie),gride(je),ei(ie),eilog(ie),fi(il,ie),psdi(il,ie),&
       d4(ir,iw,ik)
  integer:: iba(ip),irm(ip)
  !----------------------------------------------------------------------------

  e_rest=xmass(js)*c*c/1.6e-16       ! rest energy in keV

  do i=1,ir
     do j=1,ip
        do k=1,iw
           do m=1,ik
              f2(i,j,k,m)=0.          ! set to zero initiallly
           enddo
        enddo
     enddo
  enddo

  ! Find iw1(m), iw2(m)
  e_min=0.5*gride(1)
  e_max=3.*gride(je)
  if (itype.eq.1) then
     do m=1,ik
        iw1(m)=1
        do k=iw,1,-1
           if (ekev(1,1,k,m).lt.e_min) then
              iw1(m)=k
              goto 1
           endif
        enddo
1       iw2(m)=iw
        do k=1,iw
           if (ekev(iba(1),1,k,m).gt.e_max) then
              iw2(m)=k
              goto 2
           endif
        enddo
2    enddo
  else
     read(iUnit2) f2             
     read(iUnit2) iw1
     read(iUnit2) iw2
     read(iUnit2) xnsw0,vsw0,Bx0,By0,Bz0,vswb0,xnswb0
  endif

  ! Set up f2 and f at t=tstart for initial run
  if (init.gt.0) then     ! initial flux from AE8MAX (e-) or AP8MAX+CCE (H+)
     read(iUnit1,*) (roi(i),i=1,il)    ! radial distance at which input is taken
     read(iUnit1,*) (ei_MeV(k),k=1,ie)  ! in MeV
     do k=1,ie
        ei(k)=ei_MeV(k)*1000.        ! convert ei to keV
        eilog(k)=log10(ei(k))
        eir=ei(k)/e_rest             ! normalized energy in rest energy
        pp=xmass(js)*c*sqrt(eir*(eir+2.))
        read(iUnit1,*) (fi(i,k),i=1,il)      ! average flux in /cm2secMeV
        do i=1,il
           flux=fi(i,k)/4./pi/1000.                  ! flux in /cm2secsrkeV
           psdi(i,k)=log10(flux/(1.6e19*pp)/pp)      ! log(psd)
        enddo
     enddo
     ! get psd at the grids by linear interpol. Then find conser. density
     do j=1,ip
        do i=1,irm(j)
           roii=ro(i,j)
           if (roii.lt.roi(1)) roii=roi(1)    ! flat dist. at low L
           if (roii.gt.roi(il)) roii=roi(il)  ! flat dist. @ high L
           do m=1,ik
              do k=1,iw2(m)
                 e1=log10(ekev(i,j,k,m))
                 if (e1.lt.eilog(1)) e1=eilog(1)    ! flat dist. at low E
                 if (e1.gt.eilog(ie)) e1=eilog(ie)  ! flat dist. at high E
                 call lintp2(roi,eilog,psdi,il,ie,roii,e1,x)
                 psd=10.**x
                 f2(i,j,k,m)=psd*xjac(i,k)*1.e20*1.e19       ! conser f2
              enddo
           enddo
        enddo
     enddo
  endif                     ! end of if (init.gt.0)
  close(iUnit1)                 ! close the file for initial condition

  ! Setup elb, eub, e_l, ecbf, ecdt, eclc, ecce
  elb=e_min             ! keV
  eub=e_max             ! keV
  do i=1,ir
     e_l(i)=0.
     ecbf(i)=0.
     ecdt(i)=0.
     eclc(i)=0.
     ecce(i)=0.
     do j=1,ip
        do m=1,ik
           do k=iw1(m),iw2(m)
              ekev1=ekev(i,j,k,m)
              if (ekev1.ge.elb.and.ekev1.le.eub) &
                   e_l(i)=e_l(i)+f2(i,j,k,m)*ekev1*d4(i,k,m)
           enddo
        enddo
     enddo
  enddo
  if (itype.eq.1) then
     if(UseSeparatePlotFiles) then
        open(unit=UnitTmp_,file='RB/'//outnameSepOrig//st2//'.ec')
     else
        open(unit=UnitTmp_,file='RB/'//outname//st2//'.ec')
     end if
     write(UnitTmp_,*) elb,eub,'      ! elb,eub'
     close(UnitTmp_)
  else
     read(iUnit2) ecbf
     read(iUnit2) ecdt
     read(iUnit2) eclc
     read(iUnit2) ecce
     close(iUnit2)
  endif

end subroutine initial


!***************************************************************************
!                              boundary
!  Rountine finds the instantaneous boundary conditions on the dayside and
!  nightside side.
!***************************************************************************
subroutine boundary(t,tstart,f2,v,xjac,xmass,p,xktd,xnd,&
     vswb0,xnswb0,itype,ibset,irm,irm0,iba)
  
  use rbe_cboundary
  use rbe_cread2,  ONLY:js,tsw,xnswa,vswa,nsw,UseGm,tint, UseMhdBoundary
  use ModIoUnit,   ONLY: UnitTmp_
  use ModNumConst, ONLY: pi => cPi
  use rbe_cread1
  use ModGmRb,     ONLY: StateIntegral_IIV,AveDens_,AveP_
  real v(ir,ip,iw,ik),xjac(ir,iw),p(ir,ip,iw,ik),&
       xmass(ns),f2(ir,ip,iw,ik)
  integer iba(ip),irm(ip),irm0(ip)

  ibset=4                      ! Kappa distribution on the nightside

  if (js.eq.1) xkappa=3.        ! electrons
  if (js.eq.2) xkappa=4.        ! ions
  xk1=xkappa+1.
  xk2=xkappa-0.5

  !  Calculate temperature and density at boundary using vswb and xnswb
  if (t.eq.tstart.and.itype.eq.2) then
     vswb=vswb0
     xnswb=xnswb0
  else
     if (UseGM) then
        xnswb = xnswa(1)
        vswb  = vswa(1)
     else
        t1=t-2.*3600.                             ! 2 hours lag from SW to PS
        if (t1.lt.tsw(1)) t1=tsw(1)
        if (t1.gt.tsw(nsw)) t1=tsw(nsw)
        call lintp(tsw,xnswa(1:nsw),nsw,t1,xnswb)
        call lintp(tsw,vswa(1:nsw),nsw,t1,vswb)       ! vsw in km/s
     endif
  endif
  xnn=(0.02*xnswb+0.316)*sqrt(xmass(js)/xmass(2))*1.e6      ! xnn in m^-3
  if (js.eq.1) xktn=0.016*vswb-2.4          ! e- kT or Eo in keV
  if (js.eq.2) xktn=0.012*vswb-1.8          ! H+ kT or Eo in keV
  if (xktn.le.0.) then
     write(*,*) ' xktn.le.0. '
     call CON_stop('RBE ERROR')
  endif

  !  Assume a Maxwellian at nightside when ibset=3 and Kappa when ibset=4
  if (ibset.eq.3) factorn=xnn/(2.*pi*xmass(js)*xktn*1.6e-16)**1.5  
  if (ibset.eq.4) factorn=xnn*exp(gammln(xk1))/exp(gammln(xk2))/&
       (2.*pi*xkappa*xmass(js)*xktn*1.6e-16)**1.5

  !  Assume a Maxwellian at the dayside magnetopause
  xktd=0.3                  ! temperature in keV at dayside magnetopause
  xnd=xnswb           ! density in cm-3. Summer 2006: it was set to 0 by mistake
  xktd1=xktd*1000.          ! kT in eV  
  xnd1=xnd*1.e6             ! n in m^-3
  factord=xnd1/(2.*pi*xmass(js)*xktd1*1.6e-19)**1.5   

  !  Calculate psd at boundary. Use nightside condition for all local times
  chmass=1.6e-19/xmass(js)
  do j=1,ip
     ib=iba(j)
     ib1=ib+1
     if (ib1.gt.irm(j)) ib1=irm(j)
     
     ! If using MHD values for boundary, set factors
     if (UseMhdBoundary) then
         xktn = &
             StateIntegral_IIV(irm(j),j,AveP_)&
             / StateIntegral_IIV(irm(j),j,AveDens_) *6.2415e15 !J-->KeV 
        xnn  = StateIntegral_IIV(irm(j),j,AveDens_)
        !  Assume a Maxwellian at nightside when ibset=3 and Kappa when ibset=4
        if (ibset.eq.3) factorn=xnn/(2.*pi*xmass(js)*xktn*1.6e-16)**1.5  
        if (ibset.eq.4) factorn=xnn*exp(gammln(xk1))/exp(gammln(xk2))/&
             (2.*pi*xkappa*xmass(js)*xktn*1.6e-16)**1.5
     endif
     
     do m=1,ik
        do k=1,iw
           if (ibset.eq.3) then               ! Maxwellian
              v2=v(ib1,j,k,m)*v(ib1,j,k,m)
              fbb=factorn*exp(-v2/2./xktn/1000./chmass)   ! xktn in keV 
              fb(j,k,m)=fbb*xjac(ib1,k)
           elseif (ibset.eq.4) then           ! Kappa
              ekb=0.5*p(ib1,j,k,m)/xmass(js)*p(ib1,j,k,m)/1.6e-16
              fbb=factorn/(1.+ekb/xktn/xkappa)**xk1
              fb(j,k,m)=fbb*xjac(ib1,k)
           endif                           ! end of if(ibset.eq.*)
           ! fill f2 at boundaries
           if (t.gt.tstart) then
              do i=irm0(j)+1,ib
                 f2(i,j,k,m)=fbb*xjac(i,k) ! psd=fb between irm0 and ib
                 if (m.eq.1.and.k.eq.1) write(*,*) 't,i,j,irm0(j),ib ', &
                      t,i,j,irm0(j),ib
              enddo
           endif
           
           do i=ib+1,irm(j)
              f2(i,j,k,m)=fbb*xjac(i,k)    ! psd=fb between ib and irm
           enddo
           do i=irm(j)+1,ir
              f2(i,j,k,m)=0.
           enddo
        enddo                                 ! end of k loop
     enddo                                    ! end of m loop
  enddo                                       ! end of j loop
  
  !  Write boundary condition in file outname_st2.bc
  xnn_cm3=xnn/1.e6
  if (t.ge.tstart) then
  thour=t/3600.
     if (itype.eq.1.and.t.eq.tstart) then
        if(UseSeparatePlotFiles) then
           open(unit=UnitTmp_,file='RB/'//outnameSepOrig//st2//'.bc')
        else
           open(unit=UnitTmp_,file='RB/'//outname//st2//'.bc')
        end if
        write(UnitTmp_,*)&
             '                   nightside BC           dayside BC'
        write(UnitTmp_,*) &
             '     t(hour)    n(cm^-3)    kT(keV)   n(cm^-3)    kT(keV)'
        write(UnitTmp_,'(f12.2,2(f11.4,f11.3))') thour,xnn_cm3,xktn,xnd,xktd
        close(UnitTmp_)
     elseif (mod(t,tint).eq.0.) then
        if(UseSeparatePlotFiles) then
           open(unit=UnitTmp_,file='RB/'//outnameSepOrig//st2//'.bc',&
                status='old',position='append')
        else
           open(unit=UnitTmp_,file='RB/'//outname//st2//'.bc',&
                status='old',position='append')
        end if
        write(UnitTmp_,'(f12.2,2(f11.4,f11.3))') thour,xnn_cm3,xktn,xnd,xktd
        close(UnitTmp_)
     endif
  endif
  
  ! Reset irm0
  irm0(1:ip)=irm(1:ip)
end subroutine boundary


! **************************************************************************
!                           p_result
!  Routine prints output fluxes and informations needed for continuous run.
!***************************************************************************
subroutine p_result(t,tstart,f2,rc,xlati,ekev,y,p,ro,xmlto,xmlt,&
     xjac,gride,gridp,gridy,bo,xnsw,vsw,Bx,By,Bz,&
     vswb,xnswb,parmod,ecbf,ecdt,eclc,ecce,density,iprint,&
     ntime,irm,iplsp,iw1,iw2,itype, DoSaveRestart, DoSavePlot)

  use rbe_grid
  use rbe_cread1,ONLY: UseSeparatePlotFiles
  use rbe_cread2,ONLY: js,storm,DoSaveIe
  use ModIoUnit, ONLY: UnitTmp_
  use ModNumConst, ONLY: pi => cPi
  use ModRbSat,    ONLY: write_rb_sat, nRbSats, DoWriteSats
  use rbe_cread1
  use ModWriteTec, ONLY: write_tec, DoWriteTec
  real xlati(ir),ekev(ir,ip,iw,ik),y(ir,ip,0:ik+1),bo(ir,ip),&
       xjac(ir,iw),gride(je),gridp(je),gridy(ig),f2(ir,ip,iw,ik),ro(ir,ip),&
       xmlto(ir,ip),f(ir,ip,iw,ik),xlati1(ir),p(ir,ip,iw,ik),xmlt(ip),&
       flx(ir,ip,je,ig),ecbf(ir),ecdt(ir),eclc(ir),ecce(ir),&
       psd(ir,ip,iw,ik),ebound(je+1),density(ir,ip),parmod(10)
  integer iw1(ik),iw2(ik),irm(ip),iSat
  logical, intent(in) :: DoSaveRestart, DoSavePlot
  hour=t/3600.
  do i=1,ir
     xlati1(i)=xlati(i)*180./pi   ! lat. at ionosphere in degree
  enddo

  ! Find Ebound
  do k=2,je
     Ebound(k)=sqrt(gride(k-1)*gride(k))
  enddo
  Ebound(1)=gride(1)*gride(1)/Ebound(2)
  Ebound(je+1)=gride(je)*gride(je)/Ebound(je)

  ! Open the files to write fluxes and psd                  
  iwh=ifix(0.5*(iw+1))
  ikh=ifix(0.5*(ik+1))
  if (DoSavePlot) then
     if (UseSeparatePlotFiles) then
        open(unit=UnitTmp_,file='RB/plots/'//outnameSep//st2//'.fls',&
             status='unknown')
        write(UnitTmp_,'(f10.5,5i6,"         ! rc(Re),ir,ip,je,ig,ntime")')&
             rc,ir,ip,je,ig,ntime
        write(UnitTmp_,'(6f9.3)') (gride(k),k=1,je)
        !     write(UnitTmp_,'(7f9.3)') (Ebound(k),k=1,je+1)
        write(UnitTmp_,'(6f9.5)') (gridy(m),m=1,ig)
        write(UnitTmp_,'(10f8.3)') (xlati1(i),i=1,ir)
        if (iprint.eq.1) then
           close(UnitTmp_)
           return
        endif
     else
        if (t.eq.tstart) then
           open(unit=UnitTmp_,file='RB/plots/'//outname//st2//'.fls',status='unknown')
           !        open(unit=13,file=outname//st2//'.psd',status='unknown')
           write(UnitTmp_,'(f10.5,5i6,"         ! rc(Re),ir,ip,je,ig,ntime")')&
                rc,ir,ip,je,ig,ntime
           write(UnitTmp_,'(6f9.3)') (gride(k),k=1,je)
           !     write(UnitTmp_,'(7f9.3)') (Ebound(k),k=1,je+1)
           write(UnitTmp_,'(6f9.5)') (gridy(m),m=1,ig)
           write(UnitTmp_,'(10f8.3)') (xlati1(i),i=1,ir)
           !        write(13,'(f10.5,5i6,"         ! rc(Re),iwh,ikh,ir,ip,ntime")')
           !    *         rc,iwh,ikh,ir,ip,ntime
           !        write(13,'(1p,7e11.3)') (w(k),k=1,iw,2)
           !        write(13,'(1p,7e11.3)') (si(m),m=1,ik,2)
           !        write(13,'(10f8.3)') (xlati1(i),i=1,ir)
           if (iprint.eq.1) then
              close(UnitTmp_)
              !           close(13)
              return
           endif
        else                                                  ! in pbo_2.f
           open(unit=UnitTmp_,file='RB/plots/'//outname//st2//'.fls',status='old',position='append')
           !        open(unit=13,file=outname//st2//'.psd',status='old',position='append')
        endif
     endif
  endif
  ! Convert f2 to f (differential flux)
  f(:,:,:,:)=0.0
  do i=1,ir
     do j=1,ip
        do m=1,ik
           do k=iw1(m),iw2(m)
              psd(i,j,k,m)=f2(i,j,k,m)/xjac(i,k)/1.e20/1.e19 !mug^-3cm^-6s^3
              f(i,j,k,m)=0.                          ! differential flux
              if (i.le.irm(j))&
                   f(i,j,k,m)=psd(i,j,k,m)*(1.6e19*p(i,j,k,m))*p(i,j,k,m)
           enddo
        enddo
     enddo
  enddo

  ! Set values beyond irm
  do j=1,ip
     do i=irm(j)+1,ir
        ro(i,j)=ro(irm(j),j) 
        bo(i,j)=bo(irm(j),j)
        xmlto(i,j)=xmlto(irm(j),j)
        do m=0,ik+1
           y(i,j,m)=y(irm(j),j,m)
        enddo
     enddo
  enddo

  ! Calculate and write fluxes at fixed E and y grides. 
  call fluxes(f,y,p,gridp,ekev,gride,gridy,irm,iw1,iw2,flx) 
  if (DoSavePlot) write(UnitTmp_,'(f7.2,10f9.2,"   hour,parmod(1:10)")') hour,parmod         
  !     write(13,'(1x,7hhour = ,f6.2,10f9.2,"   parmod(1:10)")') hour,parmod
  do i=1,ir             ! Write fluxes @ fixed E & y grids
     do j=1,ip
        if (i.gt.irm(j)) then
           ro(i,j)=ro(irm(j),j)
           xmlto(i,j)=xmlto(irm(j),j)
           bo(i,j)=bo(irm(j),j)
           density(i,j)=density(irm(j),j)
           flx(i,j,1:je,1:ig)=flx(irm(j),j,1:je,1:ig)
        endif
        if (DoSavePlot) write(UnitTmp_,'(f7.2,f6.1,2f8.3,1pe11.3,0p,i5,1pe11.3)')&
             xlati1(i),xmlt(j),ro(i,j),xmlto(i,j),bo(i,j),irm(j),density(i,j)
        !           write(13,'(f7.2,f6.1,2f8.3,1pe11.3,0p,i5)')
        !    *                xlati1(i),xmlt(j),ro(i,j),xmlto(i,j),bo(i,j),irm(j)
        if (DoSavePlot) then
           do k=1,je
              write(UnitTmp_,'(1p,12e11.3)') (flx(i,j,k,m),m=1,ig)
           enddo
        endif
        !           do k=1,iw,2
        !              write(13,'(1p,12e11.3)') (psd(i,j,k,m),m=1,ik,2)
        !           enddo
     enddo
  enddo
  if (DoSavePlot) close(UnitTmp_)
  if (DoWriteTec .and. DoSavePlot) then
     call write_tec(t,flx,ebound)
  endif
  !Write out any sats that are being tracked
  if(DoWriteSats .and. DoSavePlot) then
     do iSat=1,nRbSats
        call write_rb_sat(iSat,ir,ip,je,ig,flx)
     enddo
  endif
  !     close(13)

  ! Write energy changes from various processes
  if (DoSavePlot) then
     if(UseSeparatePlotFiles) then
        open(unit=UnitTmp_,file='RB/'//outnameSepOrig//st2//'.ec',status='old',position='append')
     else
        open(unit=UnitTmp_,file='RB/'//outname//st2//'.ec',status='old',position='append')
     endif
     write(UnitTmp_,*) hour,'     ! hour'
     write(UnitTmp_,*) '   i  ro(i,1)       ecbf          ecdt          ecce',&
          '          eclc'
     do i=1,ir
        ro1=ro(i,1)
        write(UnitTmp_,'(i4,f9.2,1p,4e14.4)') i,ro1,ecbf(i),ecdt(i),ecce(i),eclc(i)
     enddo
     write(UnitTmp_,'(8x,"total",1p,4e14.4)') &
          sum(ecbf),sum(ecdt),sum(ecce),sum(eclc)
     close(UnitTmp_)
  endif

  ! Open files to write all the information for continous run        
  if (DoSaveRestart) then
     if(UseSeparatePlotFiles)then
        open(unit=UnitTmp_,file='RB/restartOUT/'//outnameSep//st2//'_c.f2',form='unformatted')
     else
        open(unit=UnitTmp_,file='RB/restartOUT/'//outname//st2//'_c.f2',form='unformatted')
     end if
     write(UnitTmp_) f2   
     write(UnitTmp_) iw1
     write(UnitTmp_) iw2
     write(UnitTmp_) xnsw,vsw,Bx,By,Bz,vswb,xnswb
     write(UnitTmp_) ecbf
     write(UnitTmp_) ecdt
     write(UnitTmp_) eclc
     write(UnitTmp_) ecce
     close(UnitTmp_)
     if (iplsp.eq.1) call saveplasmasphere(t,tstart,itype)
     
     ! Write the restart.H file to be included at restart
     open(unit=UnitTmp_,file='RB/restartOUT/restart.H')
     
     write(UnitTmp_,'(a)') '#TIMESIMULATION'
     write(UnitTmp_,'(es15.8,a25)') t,'tSimulation'
     write(UnitTmp_,'(a)') '#RESTART'
     write(UnitTmp_,'(a,a25)') 'T', 'IsRestart'
     write(UnitTmp_,'(a)') '#INPUTDATA'
     write(UnitTmp_,'(a,a25)') storm, 'NameStorm'
     write(UnitTmp_,'(a)') '#SPECIES'
     if (js == 1 ) write(UnitTmp_,'(a,a32)') 'e','NameSpecies'
     if (js == 2 ) write(UnitTmp_,'(a,a32)') 'H','NameSpecies'
     
     close(UnitTmp_)
  endif
  ! Write to log file
  if (t.eq.tstart .and. DoSavePlot) write(*,'(a8)') outname
  if (DoSavePlot) write(*,*) 't(hour)   ',t/3600. 
  
  ! Write the potential values
  if (DoSaveIe .and. DoSavePlot) call RB_plot_potential
end subroutine p_result


!***********************************************************************
!                                fluxes
! Routine calculates the fluxes at fixed energy and pitch angle grids
!***********************************************************************
subroutine fluxes(f,y,p,gridp,ekev,gride,gridy,irm,iw1,iw2,flx)
  use rbe_grid
  real aloge(je),y(ir,ip,0:ik+1),y1(ik),ekev(ir,ip,iw,ik),e1(iw),h1(iw),&
       f2d(je,ik),xj(ik),f(ir,ip,iw,ik),gride(je),gridp(je),gridy(ig),&
       flx(ir,ip,je,ig),p(ir,ip,iw,ik)
  integer iw1(ik),iw2(ik),irm(ip)

  do k=1,je
     aloge(k)=log10(gride(k))
  end do

  ! Calculate fluxes at certain energy and pitch-angle grids 
  do j=1,ip
     do i=1,irm(j)

        do m=1,ik
           k1=iw1(m)
           k2=iw2(m)
           y1(m)=y(i,j,m)
           do k=1,iw
              e1(k)=log10(ekev(i,j,k,m))
              x=f(i,j,k,m)        
              if (x <= 1.e-50) then
                 h1(k)=-50.
              else
                 h1(k)=log10(x)
              endif
           end do
           do k=1,je
              if (aloge(k).le.e1(k1)) then
                 f2d(k,m)=h1(k1)+2.*log10(gridp(k)/p(i,j,k1,m))  ! flat psd
              elseif (aloge(k).ge.e1(k2)) then
                 f2d(k,m)=h1(k2)+2.*log10(gridp(k)/p(i,j,k2,m))  ! flat psd
              else
                 call lintp(e1,h1,iw,aloge(k),x)
                 f2d(k,m)=x        ! log(flux) at fixed E grid
              endif
           end do
        end do
        do k=1,je
           do m=1,ik
              xj(m)=f2d(k,m)
           end do
           do m=1,ig
              if (gridy(m).le.y1(ik)) then
                 x=xj(ik)   ! flat dist. when beyond grid 
              elseif (gridy(m).ge.y1(1)) then
                 x=xj(1)    ! flat dist. when beyond grid 
              else
                 call lintp(y1,xj,ik,gridy(m),x)
              endif
              flx(i,j,k,m)=10.**x       !flux @ fixed E & y grids   
           enddo
        enddo

     enddo
  enddo

end subroutine fluxes


!***********************************************************************
!                             losscone
! Routine removes losscone particles with a lifetime of half the bounce
! period
!***********************************************************************
subroutine losscone(f2,tcone,rmir,rc,iba,iw1,iw2)
  use rbe_grid
  real f2(ir,ip,iw,ik),tcone(ir,ip,iw,ik),rmir(ir,ip,ik)
  integer iw1(ik),iw2(ik),iba(ip)

  do j=1,ip
     do i=1,iba(j)
        do m=1,ik
           do k=iw1(m),iw2(m)

              if (rmir(i,j,m).lt.rc) then
                 if (tcone(i,j,k,m).eq.0.) then
                    f2(i,j,k,m)=0.
                 else 
                    f2(i,j,k,m)=f2(i,j,k,m)*tcone(i,j,k,m)
                 endif
              endif

           enddo
        enddo
     enddo
  enddo

end subroutine losscone


!***********************************************************************
!                            drift
!  Routine calculate the change of distribution function due to drift.
!  Time step is reduced if Courant number is greater 1.
!***********************************************************************
subroutine drift(t,dt,f2,vl,vp,ro,rb,fb,dlati,dphi,iba,&
     iw1,iw2,irm)
  use rbe_grid
  use rbe_cread2, ONLY: UseSplitting

  real vl(ir,ip,iw,ik),vp(ir,ip,iw,ik),&
       f2(ir,ip,iw,ik),ro(ir,ip),fb(ip,iw,ik),dlati(ir),cl(ir,ip),cp(ir,ip)
  integer iba(ip),iw1(ik),iw2(ik),irm(ip)

  logical :: IsOdd
  !-----------------------------------------------------------------------
  do m=1,ik
     do k=iw1(m),iw2(m)      

        ! calculate maximum Courant number and nrun
        cmax=0.
        do j=1,ip
           do i=1,iba(j)
              cl1=dt/dlati(i)*vl(i,j,k,m)
              cp1=0.
              cp1=dt/dphi*vp(i,j,k,m)
              cmx=max(abs(cl1),abs(cp1))
              cmax=max(cmx,cmax)
           enddo
        enddo

        nrun=ifix(cmax/0.35)+1     ! Courant number can't be greater 0.35
        dt1=dt/nrun                ! new dt

        ! new Courant numbers
        do j=1,ip
           do i=1,ir               ! summer 2006: change from do i=1,iba(j) to
              cl(i,j)=dt1/dlati(i)*vl(i,j,k,m)   ! make sure all "velocities"
              cp(i,j)=dt1/dphi*vp(i,j,k,m)       ! are reduced.
           enddo
        enddo

        ! run driftl and driftp
        IsOdd = .true.
        do n=1,nrun
           if(UseSplitting)then
              if(IsOdd)then
                 call driftl(t,dt1,f2,k,m,vl,ro,rb,fb,dlati,iba,irm)
                 call driftp(t,dt1,f2,k,m,vp,fb,dphi,iba,irm)
              else
                 call driftp(t,dt1,f2,k,m,vp,fb,dphi,iba,irm)
                 call driftl(t,dt1,f2,k,m,vl,ro,rb,fb,dlati,iba,irm)
              end if
              IsOdd = .not. IsOdd
           else
              call driftlp(t,dt1,f2,k,m,vl,cl,cp,ro,rb,fb,dlati,&
                   iba,irm,n)
           end if

           ! Check for NaN-s: not (NaN < 2) and not (NaN > 1) is true
           do j=1,ip
              do i=1,iba(j)
                 if(.not. f2(i,j,k,m) < 2.0 .and. .not. f2(i,j,k,m) > 1.0)then
                    write(*,*)'i,j,k,m,n=',i,j,k,m,n
                    call CON_stop('RBE ERROR: NaN found in f2')
                 end if
              end do
           end do

        enddo
     enddo
  enddo

end subroutine drift


!*******************************************************************************
!                             driftl
!  Routine calculate the change of distribution function due to
!  drift in ionospheric latitude (radial drift at the eqautor).
!  Fluxes at boundary are assumed to be local time symmetric.
!*******************************************************************************
subroutine driftl(t,dt1,f2,k,m,vl,ro,rb,fb,dlati,iba,irm)
  use rbe_grid
  real f(ir),fa(0:ir),c(ir),vl(ir,ip,iw,ik),&
       f2(ir,ip,iw,ik),ro(ir,ip),fb(ip,iw,ik),ro1(ir),dlati(ir),al(0:ir)
  integer iba(ip),irm(ip)

  ibc=1                     ! advective boundary condition

  do j=1,ip
     ib=iba(j)
     do i=1,irm(j)
        ro1(i)=ro(i,j)
     enddo
     fb1=fb(j,k,m)
     do i=0,ib
        if (i==0) then
           al(i)=vl(i+1,j,k,m)
        else
           al(i)=vl(i,j,k,m)
        endif
        if (i.ge.1) c(i)=dt1/dlati(i)*al(i)             ! Courant number
     enddo

     do i=1,ib
        f(i)=f2(i,j,k,m)
     end do
     fb0=f(1)              ! psd at inner bound

     fa(0)=f(1)     
     call FLS_2or(ibc,fb0,fb1,ib,c,f,fa(1))         ! calculate fa(1:ib)
     !        call FCT_2or(ibc,fb0,fb1,ib,c,f,fa(1))         ! calculate fa(1:ib)

     do i=1,ib
        f2(i,j,k,m)=f(i)+dt1/dlati(i)*(al(i-1)*fa(i-1)-al(i)*fa(i))
        if (f2(i,j,k,m).lt.0.) then
           if (f2(i,j,k,m).gt.-1.e-50) then
              f2(i,j,k,m)=0.
           else
              write(6,*)' f2 < 0 in driftl ',i,j,k,m
              write(6,*) t,f2(i,j,k,m)
              write(6,*) ib,fb(j,k,m)
              do ii=1,ir
                 write(6,'(1p,5e13.3)') f(ii),fa(ii),c(ii),vl(ii,j,k,m)
              enddo
              call CON_stop('RBE ERROR')
           endif
        endif
     end do
  enddo

end subroutine driftl


!***********************************************************************
!                             driftp
!  Routine calculate the change of distribution function due to
!  azimuthal drift.
!***********************************************************************
subroutine driftp(t,dt1,f2,k,m,vp,fb,dphi,iba,irm)
  use rbe_grid
  real fa(ip),f(ip),c(ip),vp(ir,ip,iw,ik),f2(ir,ip,iw,ik),fb(ip,iw,ik)
  integer iba(ip),irm(ip)

  ibc=2                     ! periodic boundary condition
  fb0=0.                    ! dummy for periodic boundary condition
  fb1=0.                    ! dummy for periodic boundary condition

  do i=1,ir         

     do j=1,ip
        j1=j+1
        if (j1.gt.ip) j1=j1-ip
        c(j)=dt1/dphi*vp(i,j,k,m)               ! Courant number
     enddo

     do j=1,ip
        f(j)=f2(i,j,k,m)
        if (i.gt.iba(j)) f(j)=fb(j,k,m)
     end do

     call FLS_2or(ibc,fb0,fb1,ip,c,f,fa)         ! calculate fa
     !        call FCT_2or(ibc,fb0,fb1,ip,c,f,fa)         ! calculate fa

     do j=1,ip
        j0=j-1
        if (j0.lt.1) j0=j0+ip
        f2(i,j,k,m)=f(j)           ! initial value
        if (i.le.irm(j)) f2(i,j,k,m)=f(j)-c(j)*fa(j)+c(j0)*fa(j0)
        if (f2(i,j,k,m).lt.0.) then
           if (f2(i,j,k,m).gt.-1.e-50) then
              f2(i,j,k,m)=0.
           else
              write(6,*)' f2 < 0 in driftp ',i,j,k,m 
              write(6,*) t,f2(i,j,k,m)
              do jj=1,ip
                 write(6,'(1p,5e13.3)') f(jj),fa(jj),c(jj),vp(i,jj,k,m)
              enddo
              call CON_stop('RBE ERROR')
           endif
        endif
     enddo

  enddo

end subroutine driftp


!*******************************************************************************
!                             driftlp
!  Routine calculate the change of distribution function due to
!  drift in ionospheric latitude (radial drift at the eqautor) and local time.
!  Fluxes at boundary are assumed to be local time symmetric.
!*******************************************************************************
subroutine driftlp(t,dt1,f2,k,m,vl,cl,cp,ro,rb,fb,dlati,iba,irm,n)
  use rbe_grid
  real f(ir,ip),fal(ir,ip),fap(ir,ip),cl(ir,ip),cp(ir,ip),fupl(ir,ip),&
       vl(ir,ip,iw,ik),fupp(ir,ip),&
       f2(ir,ip,iw,ik),ro(ir,ip),fb(ip,iw,ik),dlati(ir)
  integer iba(ip),irm(ip)

  f(1:ir,1:ip)=f2(1:ir,1:ip,k,m)         ! initial f2

  ! Calculate fluxes at grid faces: fal, fap
  call interFlux(cl,cp,f,fal,fap,fupl,fupp)


  ! Update f2
  iup=0
1  do j=1,ip
     j_1=j-1
     if (j_1.lt.1) j_1=j_1+ip
     do i=2,iba(j)
        f2(i,j,k,m)=f(i,j)+dt1/dlati(i)*&
             (vl(i-1,j,k,m)*fal(i-1,j)-vl(i,j,k,m)*fal(i,j))+&
             cp(i,j_1)*fap(i,j_1)-cp(i,j)*fap(i,j)
        
        if (f2(i,j,k,m).lt.0..and.f2(i,j,k,m).gt.-1.e-20) f2(i,j,k,m)=0.
        if (f2(i,j,k,m).lt.0.) then
           if (iup.eq.0) then
              iup=1         ! use upwind scheme when high res. method fails
              fal(1:ir,1:ip)=fupl(1:ir,1:ip)
              fap(1:ir,1:ip)=fupp(1:ir,1:ip)
              goto 1
           else
              write(*,'(a,4i4)') 'f2 < 0 in driftlp, i,j,k,m ',i,j,k,m
              write(*,'(a,f12.0,1p,2e15.4)') 't,f2,f ',t,f2(i,j,k,m),f(i,j)
              write(*,'(a,4i6,1p,e15.4)')&
                   'irm(j),irm(j_1),iba(j),iba(j_1),fb(j,k,m) ',&
                   irm(j),irm(j_1),iba(j),iba(j_1),fb(j,k,m)
              write(*,*) 'dt1,dlati(i) ',dt1,dlati(i)
              write(*,'(a,1p,4e12.3)') 'vl(i-1),fal(i-1),vl(i),fal(i) ',&
                   vl(i-1,j,k,m),fal(i-1,j),vl(i,j,k,m),fal(i,j)
              write(*,'(a,1p,4e12.3)') 'cp(j_1),fap(j_1),cp(j),fap(j) ',&
                   cp(i,j_1),fap(i,j_1),cp(i,j),fap(i,j)
              call CON_stop('RBE ERROR')
           endif
        endif
     enddo
  enddo

end subroutine driftlp


!****************************************************************************
!                         charexchange
!  routine calculates the decay of ion distributions due to charge exchange.
!****************************************************************************
subroutine charexchange(f2,achar,iba,iw1,iw2)
  use rbe_grid
  real f2(ir,ip,iw,ik),achar(ir,ip,iw,ik)
  integer iw1(ik),iw2(ik),iba(ip)

  do j=1,ip
     do i=1,iba(j)
        do m=1,ik
           do k=iw1(m),iw2(m)
              f2(i,j,k,m)=f2(i,j,k,m)*achar(i,j,k,m)
           enddo
        enddo
     enddo
  enddo

end subroutine charexchange


!****************************************************************************
!                            E_change
! Routine calculates the energy change due to a given process
!****************************************************************************
subroutine E_change(f2,d4,ekev,elb,eub,e_l,ecxx,iba,iw1,iw2)
  use rbe_grid
  real f2(ir,ip,iw,ik),d4(ir,iw,ik),ekev(ir,ip,iw,ik),e_l(ir),ecxx(ir),&
       e0(ir)
  integer iw1(ik),iw2(ik),iba(ip)

  do i=1,ir
     ! Update e_l
     e0(i)=e_l(i)
     e_l(i)=0.
     do j=1,ip
        if (i.le.iba(j)) then
           do m=1,ik
              do k=iw1(m),iw2(m)
                 ekev1=ekev(i,j,k,m)
                 if (ekev1.ge.elb.and.ekev1.le.eub)&
                      e_l(i)=e_l(i)+f2(i,j,k,m)*ekev1*d4(i,k,m)
              enddo
           enddo
        endif
     enddo

     ! update ecxx
     dee=e_l(i)-e0(i)
     ecxx(i)=ecxx(i)+dee
  enddo

end subroutine E_change


!***********************************************************************
!                            FLS_2or
!  Routine calculates the inter-flux, f_i+0.5, using 2nd order flux
!  limited scheme with super-bee flux limiter method
!***********************************************************************
subroutine FLS_2or(ibc,fb0,fb1,ipt,c,f,fa)

  use rbe_cread2, ONLY: UseMcLimiter, BetaLimiter

  implicit none

  integer, intent(in) :: ibc, ipt
  real :: fb0, fb1, c(ipt),f(ipt),fa(ipt)

  real, allocatable :: f_new(:)

  integer :: i

  real :: x, xsign, xlimiter, fup, r, corr
  !----------------------------------------------------------------------
  allocate(f_new(0:ipt+2))

  f_new(1:ipt)=f(1:ipt)

  ! Set up boundary condition
  if (ibc.eq.1) then                  ! advective boundary condition
     f_new(0)=fb0
     f_new(ipt+1)=fb1
     f_new(ipt+2)=fb1
  else                                ! periodic boundary condition
     f_new(0)=f_new(ipt)
     f_new(ipt+1)=f_new(1)
     f_new(ipt+2)=f_new(2)
  endif

  ! find fa
  do i=1,ipt
     xsign=sign(1.,c(i))
     fup=0.5*(1.+xsign)*f_new(i)+0.5*(1.-xsign)*f_new(i+1)
     x=f_new(i+1)-f_new(i)
     if (abs(x).le.1.e-27) fa(i)=fup
     if (abs(x).gt.1.e-27) then
        if (xsign.eq.1.) r=(f_new(i)-f_new(i-1))/x
        if (xsign.eq.-1.) r=(f_new(i+2)-f_new(i+1))/x
        if (r.le.0.) fa(i)=fup
        if (r.gt.0.) then
           if(UseMcLimiter)then
              ! MC limiter with beta parameter
              xlimiter = min(BetaLimiter*r, BetaLimiter, 0.5*(1+r))
           else
              ! Superbee limiter
              xlimiter=max(min(2.*r,1.),min(r,2.))
           end if
           corr=-0.5*(c(i)-xsign)*x      
           fa(i)=fup+xlimiter*corr
        endif
     endif
  enddo

  deallocate(f_new)

end subroutine FLS_2or


!*******************************************************************************
!                            interFlux
!  Routine calculates the inter-flux, fal(i+0.5,j) and fap(i,j+0.5), using 
!  2nd order flux limited scheme with super-bee flux limiter method
!*******************************************************************************
subroutine interFlux(cl,cp,f,fal,fap,fupl,fupp)
  use rbe_grid
  use rbe_cread2, ONLY: UseMcLimiter, BetaLimiter, UseCentralDiff
  real cl(ir,ip),cp(ir,ip),f(ir,ip),fal(ir,ip),fap(ir,ip),fupl(ir,ip),&
       fwbc(0:ir+2,ip),fupp(ir,ip)
  integer iba(ip)
  
  fwbc(1:ir,1:ip)=f(1:ir,1:ip)        ! fwbc is f with boundary condition

  ! Set up boundary condition
  do j=1,ip
     fwbc(0,j)=f(1,j)
     fwbc(ir+1:ir+2,j)=f(ir,j)
  enddo

  ! find fa*
  do j=1,ip
     j_1=j-1
     j1=j+1
     j2=j+2
     if (j_1.lt.1) j_1=j_1+ip
     if (j1.gt.ip) j1=j1-ip
     if (j2.gt.ip) j2=j2-ip
     do i=1,ir-1
        ! find fal
        xsign=sign(1.,cl(i,j))
        fupl(i,j)=0.5*(1.+xsign)*fwbc(i,j)+0.5*(1.-xsign)*fwbc(i+1,j) ! upwind
        if (UseCentralDiff) then
           flw=0.5*fwbc(i,j)+0.5*fwbc(i+1,j)   ! Central diff
        else
           flw=0.5*(1.+cl(i,j))*fwbc(i,j)+0.5*(1.-cl(i,j))*fwbc(i+1,j)   ! LW
        end if
        x=fwbc(i+1,j)-fwbc(i,j)
        if (abs(x).le.1.e-27) fal(i,j)=fupl(i,j)
        if (abs(x).gt.1.e-27) then
           if (xsign.eq.1.) r=(fwbc(i,j)-fwbc(i-1,j))/x
           if (xsign.eq.-1.) r=(fwbc(i+2,j)-fwbc(i+1,j))/x
           if (r.le.0.) fal(i,j)=fupl(i,j)
           if (r.gt.0.) then
              xlimiter=max(min(2.*r,1.),min(r,2.))
              corr=flw-fupl(i,j)
              fal(i,j)=fupl(i,j)+xlimiter*corr
           endif
        endif
        ! find fap
        xsign=sign(1.,cp(i,j))
        fupp(i,j)=0.5*(1.+xsign)*fwbc(i,j)+0.5*(1.-xsign)*fwbc(i,j1)   ! upwind
        if (UseCentralDiff) then
           flw=0.5*fwbc(i,j)+0.5*fwbc(i,j1)   ! Central diff
        else
           flw=0.5*(1.+cp(i,j))*fwbc(i,j)+0.5*(1.-cp(i,j))*fwbc(i,j1)   ! LW
        end if
!        flw=0.5*(1.+cp(i,j))*fwbc(i,j)+0.5*(1.-cp(i,j))*fwbc(i,j1)   ! LW
        x=fwbc(i,j1)-fwbc(i,j)
        if (abs(x).le.1.e-27) fap(i,j)=fupp(i,j)
        if (abs(x).gt.1.e-27) then
           if (xsign.eq.1.) r=(fwbc(i,j)-fwbc(i,j_1))/x
           if (xsign.eq.-1.) r=(fwbc(i,j2)-fwbc(i,j1))/x
           if (r.le.0.) fap(i,j)=fupp(i,j)
           if (r.gt.0.) then
              if(UseMcLimiter)then
                 xlimiter = min(BetaLimiter*r, BetaLimiter, 0.5*(1+r))
              else
                 xlimiter = max(min(2.*r,1.),min(r,2.))
              end if
              !xlimiter=max(min(2.*r,1.),min(r,2.))
              corr=flw-fupp(i,j)
              fap(i,j)=fupp(i,j)+xlimiter*corr
           endif
        endif
     enddo              ! end of do i=1,ir-1
  enddo                 ! end of do j=1,ip

end subroutine InterFlux

!*******************************************************************************
  subroutine Wpower(t,density,tKp,xKph,chorusI,wLshell,wmlt,Bo,ro,xmlto,irm,nKp)
!*******************************************************************************
! Routine determines Chorus wave power (pT^2) and fpe/fce as a function of
! location
!
! Input: t,density,tKp,xKph,chorusI,wLshell,wmlt,Bo,ro,xmlto,irm,nKp
! Output: CHpower,ompe (through common block cWpower)

  use rbe_grid
  use ModWpower
  implicit none
  integer irm(ip),nKp,Kpl,i,j
  real t,tKp(nKp),xKph(nKp),chorusI(irw,ipw,3),wLshell(irw),wmlt(ipw),ro(ir,ip)
  real density(ir,ip),Bo(ir,ip),xmlto(ir,ip)
  real r_lo,r_up,xKp,chorusI2(irw,ipw+1),wmlt2(ipw+1),ro1,xmlt1,e_mass,epsilon0

  e_mass=9.11e-31          ! electron mass in kg
  epsilon0=8.8542e-12      ! permittivity of free space
  r_lo=1.5                 ! r lower boundary of Chorus diffusion
  r_up=6.5                 ! r upper boundary of Chorus diffusion

! Determine Kp level in Chorus wave amplitude data provided by Meredith
  call lintp(tKp,xKph,nKp,t,xKp)
  if (xKp.lt.2.) Kpl=1
  if (xKp.ge.2..and.xKp.lt.4.) Kpl=2
  if (xKp.ge.4.) Kpl=3
  chorusI2(1:irw,1:ipw)=chorusI(1:irw,1:ipw,Kpl)
  chorusI2(1:irw,ipw+1)=chorusI2(1:irw,1)
  wmlt2(1:ipw)=wmlt(1:ipw)
  wmlt2(ipw+1)=wmlt(1)+24.

! Determine the Chorus wave power, CHpower (in pT^2), and ompe (fpe/fce)
  CHpower=0.
  ompe=9999.         ! arbitrary big number
  do j=1,ip
     do i=1,irm(j)
        ompe(i,j)=sqrt(density(i,j)*e_mass/epsilon0)/bo(i,j)    ! fpe/fce
        ro1=ro(i,j)
        xmlt1=xmlto(i,j)
        if (xmlt1.lt.wmlt(1)) xmlt1=xmlt1+24.
        if (ro1.ge.r_lo.and.ro1.le.r_up) call lintp2(wLshell,wmlt2, &
           chorusI2,irw,ipw+1,ro1,xmlt1,CHpower(i,j))
     enddo
  enddo

  end subroutine Wpower


!****************************************************************************
!                            diffusee_E
!  Routine calculates the change of electron distributions due to
!  diffusion in E.
!
!  Chorus wave amplitude data were provided by Meredith and electron-Chorus
!  DEE is from Horne.
!****************************************************************************
      subroutine diffusee_E(f2,dt,y,bm,ekev,ro,w,xjac,CHpower,ompe,cLshell, &
                            ompea,ckeV,cPA,cDEE,iw1,iw2,iba)
      use rbe_grid
      implicit none
      real dt,pi,Eo,f2(ir,ip,iw,ik),fl(iw), &
           xlam,alam,Wpower,Wpower0,r_wave,um(iw),up(iw), &
           ompe1,ro1,PA1,cLshell(irc),ompea(ipe),CHpower(ir,ip),ompe(ir,ip),&
           fr(iw),a1d(iw),b1d(iw),c1d(iw),ekev(ir,ip,iw,ik),ro(ir,ip),&
           y(ir,ip,ik),ckeV(iwc),cPA(ipa),cDEE(irc,ipe,iwc,ipa),factor_1, &
           factor1,gjac_1,gjac1,gjac(iw),E_1(iw),E1(iw),Em,Ep, &
           DEEm,DEEp,DDm,DDp,u_mx,ump_mx,Enor(iw),dEn(iw),w(0:iw+1), &
           WM_1,WM1,xjac(ir,iw),bm(ir,ip,ik),Wo
      integer i,j,k,m,k1,k2,iww,iw1(ik),iw2(ik),iba(ip),irun,nrun,ier

      pi=acos(-1.)
      Eo=511.                  ! electron rest energy in keV
      xlam=0.5                 ! implicitness in solving diffusion equation
      alam=1.-xlam
      Wpower0=10000.           ! Horne's coeff based on wave power of 10000 pT^2
      r_wave=3.         ! only consider wave-particle interaction at ro > r_wave

      do j=1,ip
         do i=1,iba(j)
            ro1=ro(i,j)
            Wpower=CHpower(i,j)
            ompe1=ompe(i,j)                         ! fpe/fce
            if (ompe1.lt.ompea(1)) ompe1=ompea(1)   ! min value of fpe/fce

            if (ro1.gt.r_wave.and.ompe1.le.ompea(ipe).and.Wpower.gt.0.) then
               do m=1,ik
                  PA1=asin(y(i,j,m))*180./pi       ! pitch angle in degree
                  Wo=Eo*1.6e-16/bm(i,j,m)   ! normalization factor of mag moment
                  do k=1,iw
                     Enor(k)=ekev(i,j,k,m)/Eo
                     gjac(k)=(Enor(k)+1.)*sqrt(Enor(k)*(Enor(k)+2.))
                     WM_1=sqrt(w(k)*w(k-1))             ! M@lower grid
                     WM1=sqrt(w(k)*w(k+1))              ! M@upper grid
                     E_1(k)=sqrt(2.*WM_1/Wo+1.)-1. ! normalized kinetic energy
                     E1(k)=sqrt(2.*WM1/Wo+1.)-1.   ! normalized kinetic energy
                     dEn(k)=E1(k)-E_1(k)
                  enddo

                  ! determine k1 and k2, corresponding to ckeV(1) and ckeV(iwc)
                  k1=iw2(m)
                  findk1: do k=iw1(m),iw2(m)
                     if (ckeV(1).le.ekev(i,j,k,m)) then
                        k1=k
                        exit findk1
                     endif
                  enddo findk1
                  k2=iw1(m)
                  findk2: do k=iw2(m),iw1(m),-1
                     if (ckeV(iwc).ge.ekev(i,j,k,m)) then
                        k2=k
                        exit findk2
                     endif
                  enddo findk2
                  iww=k2-k1+1
                  if (k1.le.iw1(m).or.k2.ge.iw2(m)) then
                     write(*,*) ' Error: k1.le.iw1(m).or.k2.ge.iw2(m)'
                     stop
                  endif

                  u_mx=0.
                  do k=k1,k2
                     factor_1=dEn(k)*(Enor(k)-Enor(k-1))   ! normalized factor
                     factor1=dEn(k)*(Enor(k+1)-Enor(k))    !
                     gjac_1=(E_1(k)+1.)*sqrt(E_1(k)*(E_1(k)+2.))
                     gjac1=(E1(k)+1.)*sqrt(E1(k)*(E1(k)+2.))

                     ! find DEE/E^2 by interpolation
                     Em=E_1(k)*Eo             ! Em in keV
                     Ep=E1(k)*Eo              ! Ep in keV
                     call lintp4(cLshell,ompea,ckeV,cPA,cDEE,irc,ipe,iwc,ipa,&
                                 ro1,ompe1,Em,PA1,DEEm)
                     call lintp4(cLshell,ompea,ckeV,cPA,cDEE,irc,ipe,iwc,ipa,&
                                 ro1,ompe1,Ep,PA1,DEEp)
                     DEEm=DEEm*Wpower/Wpower0    ! these are actually DEE/E^2
                     DEEp=DEEp*Wpower/Wpower0    !
                     DDm=DEEm*E_1(k)*E_1(k)*gjac_1
                     DDp=DEEp*E1(k)*E1(k)*gjac1

                     um(k)=dt*DDm/factor_1/gjac(k)
                     up(k)=dt*DDp/factor1/gjac(k)
                     ump_mx=max(abs(up(k)),abs(um(k)))
                     if (ump_mx.gt.u_mx) u_mx=ump_mx
                  enddo             ! end k loop

                  ! reduce time step size if up or um is too large
                  irun=ifix(u_mx)+1
                  do k=k1,k2
                     um(k)=um(k)/irun
                     up(k)=up(k)/irun
                     a1d(k)=-xlam*um(k)
                     b1d(k)=1.+xlam*(um(k)+up(k))
                     c1d(k)=-xlam*up(k)
                  enddo
                  ! Start diffusion in E
                  do k=k1-1,k2+1
                     fl(k)=f2(i,j,k,m)/xjac(i,k)   ! fl is psd
                  enddo
                  do nrun=1,irun
                     fr(k1)=alam*um(k1)*fl(k1-1)+(1.-alam*(up(k1)+um(k1)))*&
                           fl(k1)+alam*up(k1)*fl(k1+1)+xlam*um(k1)*fl(k1-1)
                     do k=k1+1,k2-1    ! calculate the RHS of matrix equation
                        fr(k)=alam*um(k)*fl(k-1)+(1.-alam*(up(k)+um(k)))*&
                              fl(k)+alam*up(k)*fl(k+1)
                     enddo
                     fr(k2)=alam*um(k2)*fl(k2-1)+(1.-alam*(up(k2)+um(k2)))*&
                           fl(k2)+alam*up(k2)*fl(k2+1)+xlam*up(k2)*fl(k2+1)
                     call tridag(a1d(k1),b1d(k1),c1d(k1),fr(k1),fl(k1),iww,ier)
                  enddo

                  ! get back f2
                  do k=k1,k2
                     f2(i,j,k,m)=fl(k)*xjac(i,k)
                  enddo

               enddo       ! end of m=1,ik
            endif          ! end of if (ompe1.le.ompea(ipe).and.Wpower.gt.0.)

         enddo       ! end of i loop
      enddo          ! end of j loop

      end subroutine diffusee_E


!****************************************************************************
!                             diffusea_a
!  Routine solves pitch-angle diffusion in ao, the equatorial pitch angle.
!
!  Chorus wave amplitude data were provided by Meredith and electron-Chorus
!  Daa is from Horne.
!****************************************************************************
      subroutine diffusea_a(f2,xjac,ro,ekev,y,tya,dt,CHpower,ompe, &
                            cLshell,ompea,ckeV,cPA,cDaa,iba,iw1,iw2)
      use rbe_grid
      implicit none
      integer,parameter :: ie=40
      real f2(ir,ip,iw,ik),xjac(ir,iw),ro(ir,ip),ekev(ir,ip,iw,ik),&
           y(ir,ip,0:ik+1),tya(ir,ip,0:ik+1),ein_log(ie),&
           f1d(iw),e1d(iw),ein(ie),ao(0:ik+1),dao(ik),Gjac(0:ik+1),DD(0:ik+1),&
           um(ik),up(ik),a1d(ik),b1d(ik),c1d(ik),f0(0:ik+1),&
           f2d(ie,ik),df(ie,ik),df1(ie),fr(ik),ekevlog(iw,ik), &
           cLshell(irc),ompea(ipe),CHpower(ir,ip),ompe(ir,ip),ckeV(iwc), &
           cPA(ipa),cDaa(irc,ipe,iwc,ipa)
      real u_max,u_max_log,Wpower0,ompe1,Wpower,ro1,emin,emax,x0,x2,x,Daoao,DDp
      real u_mx,u_mx_log,factor1,factor_1,DDm,DDo,ump_mx,xlam,alam,dpsd,dt,ao_d
      real pi,r_wave
      integer iba(ip),iw1(ik),iw2(ik),i,j,m,k,irun,n,ier

  pi=acos(-1.)
  u_max=100.               ! maximum value of mu
  u_max_log=log10(u_max)
  Wpower0=10000.           ! Horne's coeff based on wave power of 10000 pT^2
  r_wave=3.             ! only consider wave-particle interaction at ro > r_wave

  do j=1,ip
     do i=1,iba(j)
        ompe1=ompe(i,j)                         ! fpe/fce
        if (ompe1.lt.ompea(1)) ompe1=ompea(1)   ! min value of fpe/fce
        Wpower=CHpower(i,j)
        ro1=ro(i,j)

        if (ro1.gt.r_wave.and.ompe1.le.ompea(ipe).and.Wpower.gt.0.) then
            ! Set up the energy grid, ein
            emin=1.e20
            do m=1,ik
               emin=min(emin,ekev(i,j,iw1(m),m))
            enddo
            emax=0.
            do m=1,ik
               emax=max(emax,ekev(i,j,iw2(m),m))
            enddo
            ein(1)=emin
            ein(ie)=emax
            ein_log(1)=log10(emin)
            ein_log(ie)=log10(emax)
            x0=ein_log(1)
            x2=ein_log(ie)
            x=(x2-x0)/(ie-1)
            do k=2,ie-1
               ein_log(k)=x0+(k-1)*x
               ein(k)=10.**ein_log(k)
            enddo

            ! Map psd to ein grid, f2d
            do m=1,ik
               do k=1,iw
                  f1d(k)=-50.                      ! f1d is log(psd)
                  if (f2(i,j,k,m).gt.0.) f1d(k)=log10(f2(i,j,k,m)/xjac(i,k))
                  ekevlog(k,m)=log10(ekev(i,j,k,m))
                  e1d(k)=ekevlog(k,m)              ! e1d is log(ekev)
               enddo
               do k=1,iw
                  if (k.lt.iw1(m)) f1d(k)=f1d(iw1(m))
                  if (k.gt.iw2(m)) f1d(k)=f1d(iw2(m))
               enddo

               do k=1,ie
                  if (ein_log(k).ge.e1d(1).and.ein_log(k).le.e1d(iw)) then
                     call lintp(e1d,f1d,iw,ein_log(k),x)
                     f2d(k,m)=10.**x      ! f2d is psd
                  endif
               enddo
            enddo
            ! Calculate f2d when ein is beyond the range of ekev
            do m=1,ik
               do k=1,ie
                  if (ein_log(k).lt.ekevlog(1,m)) f2d(k,m)=f2d(k,m-1)
               enddo
            enddo
            do m=ik,1,-1
               do k=1,ie
                  if (ein_log(k).gt.ekevlog(iw,m)) f2d(k,m)=f2d(k,m+1)
               enddo
            enddo

            ! calcuate ao, dao, and Gjac
            do m=0,ik+1
               ao(m)=asin(y(i,j,m))
               Gjac(m)=tya(i,j,m)*y(i,j,m)*sqrt(1.-y(i,j,m)*y(i,j,m))
            enddo
            do m=1,ik
               dao(m)=0.5*(ao(m+1)-ao(m-1))
            enddo

            df=0.
            do k=1,ie      ! *****
             if (ein(k).ge.ckeV(1).and.ein(k).le.ckeV(iwc)) then

               ! calculate DD, Daoao*Gjac
               do m=0,ik+1
                  ao_d=ao(m)*180./pi
                  call lintp4(cLshell,ompea,ckeV,cPA,cDaa,irc,ipe,iwc,ipa, &
                              ro1,ompe1,ein(k),ao_d,Daoao)
                  Daoao=Daoao*Wpower/Wpower0     ! scale daa with Wpower/Wpower0
                  DD(m)=Daoao*Gjac(m)
               enddo

               ! calculate up and um
               u_mx=0.
               do m=1,ik
                  factor_1=dao(m)*(ao(m)-ao(m-1))
                  factor1=dao(m)*(ao(m+1)-ao(m))
                  DDm=0.5*(DD(m)+DD(m-1))
                  DDp=0.5*(DD(m)+DD(m+1))
                  um(m)=dt*DDm/factor_1/Gjac(m)
                  up(m)=dt*DDp/factor1/Gjac(m)
                  ump_mx=max(abs(up(m)),abs(um(m)))
                  if (ump_mx.gt.u_mx) u_mx=ump_mx
               enddo

               ! reduce time step if u_mx > u_max
               irun=ifix(u_mx/u_max)+1
               do m=1,ik
                  um(m)=um(m)/irun
                  up(m)=up(m)/irun
               enddo

               ! determine the implicitness, xlam
               if (u_mx.le.1.) then
                  xlam=0.5              ! Crank-Nicolson
               else
                  u_mx_log=log10(u_mx)
                  xlam=0.5*u_mx_log/u_max_log+0.5
               endif
               alam=1.-xlam

               ! calculate a1d, b1d, c1d
               do m=1,ik
                  a1d(m)=-xlam*um(m)
                  b1d(m)=1.+xlam*(um(m)+up(m))
                  c1d(m)=-xlam*up(m)
               enddo

               ! start diffusion in ao
               f0(1:ik)=f2d(k,1:ik)
               do n=1,irun
                  f0(0)=f0(1)
                  f0(ik+1)=f0(ik)
                  fr(1)=alam*um(1)*f0(0)+(1.-alam*(up(1)+um(1)))*f0(1)+ &
                        alam*up(1)*f0(2)+xlam*um(1)*f0(0)
                  do m=2,ik-1    ! calculate the RHS of matrix equation
                     fr(m)=alam*um(m)*f0(m-1)+(1.-alam*(up(m)+um(m)))*f0(m)+ &
                           alam*up(m)*f0(m+1)
                  enddo
                  fr(ik)=alam*um(ik)*f0(ik-1)+(1.-alam*(up(ik)+um(ik)))*f0(ik) &
                         +alam*up(ik)*f0(ik+1)+xlam*up(ik)*f0(ik+1)
                  call tridag(a1d,b1d,c1d,fr,f0(1:ik),ik,ier)
               enddo
               do m=1,ik
                  df(k,m)=f0(m)-f2d(k,m)    ! df is differential psd
               enddo

             endif
            enddo         ! end of do k=1,ie      ! *****

            ! map psd back to M grid
            do m=1,ik
               do k=1,ie
                  df1(k)=df(k,m)
               enddo
               do k=iw1(m),iw2(m)
                  call lintp(ein_log,df1,ie,ekevlog(k,m),dpsd)
                  f2(i,j,k,m)=f2(i,j,k,m)+xjac(i,k)*dpsd
                  if (f2(i,j,k,m).lt.0.) f2(i,j,k,m)=0.
               enddo
            enddo
        endif   

     enddo
  enddo

  end subroutine diffusea_a

!!***********************************************************************
!!                            FCT_2or
!!  Routine calculates the inter-flux, f_i+0.5, using 2nd order flux
!!  corrected transport scheme with super-bee flux limiter method
!!***********************************************************************
!subroutine FCT_2or(ibc,fb0,fb1,ipt,c,f,fa)
!  real c(ipt),f(ipt),fa(ipt),f_new(0:ipt+3),fup(0:ipt+2)
!
!  f_new(1:ipt)=f(1:ipt)
!
!  ! Set up boundary condition
!  if (ibc.eq.1) then                  ! advective boundary condition
!     f_new(0)=fb0
!     f_new(ipt+1)=fb1
!     f_new(ipt+2)=fb1
!     f_new(ipt+3)=fb1
!  else                                ! periodic boundary condition
!     f_new(0)=f_new(ipt)
!     f_new(ipt+1)=f_new(1)
!     f_new(ipt+2)=f_new(2)
!     f_new(ipt+3)=f_new(3)
!  endif
!
!  ! find fup
!  do i=0,ipt+2
!     xsign=sign(1.,c(i))
!     fup(i)=0.5*(1.+xsign)*f_new(i)+0.5*(1.-xsign)*f_new(i+1)
!  enddo
!
!  ! find fa
!  do i=1,ipt
!     xsign=sign(1.,c(i))
!     x=f_new(i+1)-f_new(i)
!     if (abs(x).le.1.e-27) fa(i)=fup(i)
!     if (abs(x).gt.1.e-27) then
!        if (xsign.eq.1.) r=(fup(i)-fup(i-1))/(fup(i+1)-fup(i))
!        if (xsign.eq.-1.) r=(fup(i+2)-fup(i+1))/(fup(i+1)-fup(i))
!        if (r.le.0.) fa(i)=fup(i)
!        if (r.gt.0.) then
!           xlimiter=max(min(2.*r,1.),min(r,2.))
!           corr=-0.5*(c(i)-xsign)*x
!           fa(i)=fup(i)+xlimiter*corr
!        endif
!     endif
!  enddo
!
!end subroutine FCT_2or


!***********************************************************************
!                          closed
! Routine performs numerical integration using closed form 
!
!  S = y1dx1 + y2dx2 + .. + yidxi + .. + yn-1dxn-1 + yndxn
!                 
!  where yi is the value at the middle of dxi
!***********************************************************************
subroutine closed(n,y,dx,s)
  real y(n),dx(n)

  s=0.           
  do i=1,n     
     s=s+y(i)*dx(i) 
  enddo

end subroutine closed

!-----------------------------------------------------------------------
function hden(r)
  !-----------------------------------------------------------------------
  !  function calculates the H density at radial distance r (in RE)

  ! 1.08 = exobase in Rairden et al. [1986]
  real, parameter :: rexob=1.08           

  ar=max(0.0, log(r/rexob))   ! rexob = exobase in Rairden et al. [1986]
  hden=1.e6*exp(10.692-4.4431*ar**0.715831)      ! H density in m^-3

end function hden


!-----------------------------------------------------------------------
subroutine lintp(xx,yy,n,x,y)
  !-----------------------------------------------------------------------
  !  Routine does 1-D interpolation.  xx must be increasing or decreasing
  !  monotonically.  x is between xx(j) and xx(j+1)

  real xx(n),yy(n)
  ier = 0

  !  Make sure xx is increasing or decreasing monotonically
  do i=2,n
     if (xx(n).gt.xx(1).and.xx(i).lt.xx(i-1)) then
        write(*,*) ' lintp: xx is not increasing monotonically '
        write(*,*) n,(xx(j),j=1,n)
        call CON_stop('RBE ERROR')
     endif
     if (xx(n).lt.xx(1).and.xx(i).gt.xx(i-1)) then
        write(*,*) ' lintp: xx is not decreasing monotonically '
        write(*,*) 'i,xx(i),xx(i-1) ',i,xx(i-1),xx(i)
        write(*,*) n,(xx(j),j=1,n)
        call CON_stop('RBE ERROR')
     endif
  enddo

  !  Set ier=1 if out of range
  if (xx(n).gt.xx(1)) then
     if (x.lt.xx(1).or.x.gt.xx(n)) ier=1
  else
     if (x.gt.xx(1).or.x.lt.xx(n)) ier=1
  endif
  if (ier.eq.1) then
     write(*,*) ' Error: ier.eq.1'
     write(*,*) ' x  ',x
     write(*,*) ' xx  ',xx
     call CON_stop('RBE ERROR')
  endif
  !
  !    initialize lower and upper values
  !
  jl=1
  ju=n
  !
  !    if not dne compute a midpoint
  !
10 if(ju-jl.gt.1)then
     jm=(ju+jl)/2
     !
     !    now replace lower or upper limit
     !
     if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
        jl=jm
     else
        ju=jm
     endif
     !
     !    try again
     !
     go to 10
  endif
  !
  !    this is j
  !
  j=jl      ! if x.le.xx(1) then j=1
  !                 if x.gt.x(j).and.x.le.x(j+1) then j=j
  !                 if x.gt.x(n) then j=n-1
  d=xx(j+1)-xx(j)
  y=(yy(j)*(xx(j+1)-x)+yy(j+1)*(x-xx(j)))/d

end subroutine lintp
!
!
!------------------------------------------------------------------------------
subroutine lintp2(x,y,v,nx,ny,x1,y1,v1)
  !----------------------------------------------------------------------------
  !  Routine does 2-D interpolation.  x and y must be increasing or decreasing
  !  monotonically
  !
  real x(nx),y(ny),v(nx,ny)

  call locate1(x,nx,x1,i)
  if (i.gt.(nx-1)) i=nx-1      ! extrapolation if out of range
  if (i.lt.1) i=1              ! extrapolation if out of range
  i1=i+1
  a=(x1-x(i))/(x(i1)-x(i))

  call locate1(y,ny,y1,j)
  if (j.gt.(ny-1)) j=ny-1      ! extrapolation if out of range
  if (j.lt.1) j=1              ! extrapolation if out of range
  j1=j+1
  b=(y1-y(j))/(y(j1)-y(j))

  q00=(1.-a)*(1.-b)
  q01=(1.-a)*b
  q10=a*(1.-b)
  q11=a*b
  v1=q00*v(i,j)+q01*v(i,j1)+q10*v(i1,j)+q11*v(i1,j1)

end subroutine lintp2


!-----------------------------------------------------------------------------
subroutine lintp2a(x,y,v,nx,ny,x1,y1,v1)
  !----------------------------------------------------------------------------
  !  This sub program takes 2-d interplation. x is 2-D and y is 1-D.
  !
  implicit real*8 (a-h,o-z)
  real*8 x(nx,ny),y(ny),v(nx,ny),x1d(1000)   ! max(nx)=1000

  call locate1(y,ny,y1,j)
  j1=j+1
  if (j.eq.0.or.j1.gt.ny) then
     b=1.
     if (j.eq.0) j=j1
     if (j1.gt.ny) j1=j
  else
     b=(y1-y(j))/(y(j+1)-y(j))
  endif

  do ix=1,nx
     x1d(ix)=x(ix,j)
  enddo
  call locate1(x1d,nx,x1,i)
  i1=i+1
  if (i.eq.0.or.i1.gt.nx) then
     a=1.
     if (i.eq.0) i=i1
     if (i1.gt.nx) i1=i
  else
     a=(x1-x1d(i))/(x1d(i+1)-x1d(i))
  endif

  do ix=1,nx
     x1d(ix)=x(ix,j1)
  enddo
  call locate1(x1d,nx,x1,i2)
  i3=i2+1
  if (i2.eq.0.or.i3.gt.nx) then
     a1=1.
     if (i2.eq.0) i2=i3
     if (i3.gt.nx) i3=i2
  else
     a1=(x1-x1d(i2))/(x1d(i2+1)-x1d(i2))
  endif

  q00=(1.-a)*(1.-b)
  q01=(1.-a1)*b
  q10=a*(1.-b)
  q11=a1*b
  v1=q00*v(i,j)+q01*v(i2,j1)+q10*v(i1,j)+q11*v(i3,j1)

end subroutine lintp2a

!-------------------------------------------------------------------------------
subroutine lintp4(x,y,z,u,v,nx,ny,nz,nu,x1,y1,z1,u1,v1)

  !  This sub program takes 4-d interplation. 
  !
  real x(nx),y(ny),z(nz),u(nu),v(nx,ny,nz,nu)
  
  call locate1(x,nx,x1,i)
  call locate1(y,ny,y1,j)
  call locate1(z,nz,z1,k)
  call locate1(u,nu,u1,l)
  if (i.gt.(nx-1).or.i.lt.1.or.j.gt.(ny-1).or.j.lt.1.or.k.gt.(nz-1).or.&
       k.lt.1.or.l.gt.(nu-1).or.l.lt.1) then
     v1=0.                         ! v1=0 if out of range
     goto 1
  endif
  
  i1=i+1
  j1=j+1
  k1=k+1
  l1=l+1
  a1=(x1-x(i))/(x(i1)-x(i))
  b1=(y1-y(j))/(y(j1)-y(j))
  c1=(z1-z(k))/(z(k1)-z(k))
  d1=(u1-u(l))/(u(l1)-u(l))
  a=1.-a1
  b=1.-b1
  c=1.-c1
  d=1.-d1
  
  q0000=a*b*c*d
  q0001=a*b*c*d1
  q0010=a*b*c1*d
  q0011=a*b*c1*d1
  
  q0100=a*b1*c*d
  q0101=a*b1*c*d1
  q0110=a*b1*c1*d
  q0111=a*b1*c1*d1
  
  q1000=a1*b*c*d
  q1001=a1*b*c*d1
  q1010=a1*b*c1*d
  q1011=a1*b*c1*d1
  
  q1100=a1*b1*c*d
  q1101=a1*b1*c*d1
  q1110=a1*b1*c1*d
  q1111=a1*b1*c1*d1
  
  v1=q0000*v(i,j,k,l)+q0001*v(i,j,k,l1)+q0010*v(i,j,k1,l)+&
       q0011*v(i,j,k1,l1)+q0100*v(i,j1,k,l)+q0101*v(i,j1,k,l1)+&
       q0110*v(i,j1,k1,l)+q0111*v(i,j1,k1,l1)+q1000*v(i1,j,k,l)+&
       q1001*v(i1,j,k,l1)+q1010*v(i1,j,k1,l)+q1011*v(i1,j,k1,l1)+&
       q1100*v(i1,j1,k,l)+q1101*v(i1,j1,k,l1)+q1110*v(i1,j1,k1,l)+&
       q1111*v(i1,j1,k1,l1)
  
1 return
end subroutine lintp4


!--------------------------------------------------------------------------
subroutine locate1(xx,n,x,j)
  !--------------------------------------------------------------------------
  !  Routine return a value of j such that x is between xx(j) and xx(j+1).
  !  xx must be increasing or decreasing monotonically.
  !  If xx is increasing:
  !     If x=xx(m), j=m-1 so if x=xx(1), j=0  and if x=xx(n), j=n-1
  !     If x < xx(1), j=0  and if x > xx(n), j=n
  !  If xx is decreasing:
  !     If x=xx(m), j=m so if x=xx(1), j=1  and if x=xx(n), j=n
  !     If x > xx(1), j=0  and if x < xx(n), j=n

  real xx(n)

  !  Make sure xx is increasing or decreasing monotonically
  do i=2,n
     if (xx(n).gt.xx(1).and.xx(i).lt.xx(i-1)) then
        write(*,*) ' locate1: xx is not increasing monotonically '
        write(*,*) n, (xx(j),j=1,n)
        call CON_stop('RBE ERROR')
     endif
     if (xx(n).lt.xx(1).and.xx(i).gt.xx(i-1)) then
        write(*,*) ' locate1: xx is not decreasing monotonically '
        write(*,*) ' n, xx  ',n,xx
        call CON_stop('RBE ERROR')
     endif
  enddo

  jl=0
  ju=n+1
10 if(ju-jl.gt.1)then
     jm=(ju+jl)/2
     if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
        jl=jm
     else
        ju=jm
     endif
     go to 10
  endif
  j=jl

end subroutine locate1


!--------------------------------------------------------------------------
subroutine ilocate(ixx,n,ix,j)
  !--------------------------------------------------------------------------
  !  Routine modified from locate.  ilocate find the location of an integer
  !  in an integer array.

  integer ixx(n)

  jl=0
  ju=n+1
10 if(ju-jl.gt.1)then
     jm=(ju+jl)/2
     if((ixx(n).gt.ixx(1)).eqv.(ix.gt.ixx(jm)))then
        jl=jm
     else
        ju=jm
     endif
     go to 10
  endif
  j=jl

end subroutine ilocate

!--------------------------------------------------------------------------
subroutine indexx(n,arrin,indx)
  
  !  Routine arranges an array in ascending order.
  
  real arrin(n)
  integer indx(n)
  do  j=1,n
     indx(j)=j
  end do
  l=n/2+1
  ir=n
10 continue
  if(l.gt.1)then
     l=l-1
     indxt=indx(l)
     q=arrin(indxt)
  else
     indxt=indx(ir)
     q=arrin(indxt)
     indx(ir)=indx(1)
     ir=ir-1
     if(ir.eq.1)then
        indx(1)=indxt
        return
     endif
  endif
  i=l
  j=l+l
20 if(j.le.ir)then
     if(j.lt.ir)then
        if(arrin(indx(j)).lt.arrin(indx(j+1)))j=j+1
     endif
     if(q.lt.arrin(indx(j)))then
        indx(i)=indx(j)
        i=j
        j=j+j
     else
        j=ir+1
     endif
     go to 20
  endif
  indx(i)=indxt
  go to 10
  
end subroutine indexx

!-----------------------------------------------------------------------------
subroutine tridag(a,b,c,r,u,n,ier)
  
  !  Routine solves the tri-diagonal matrix equation.
  
  parameter (nmax=100)
  real gam(nmax),a(n),b(n),c(n),r(n),u(n)
  !
  !    problem can be simplified to n-1
  !
  if(b(1).eq.0.)then
     ier = 1
     return
  endif
  ier = 0
  bet=b(1)
  u(1)=r(1)/bet
  !
  !    decomposition and forward substitution
  !
  do j=2,n
     gam(j)=c(j-1)/bet
     bet=b(j)-a(j)*gam(j)
     !
     !    algotithm fails
     !
     if(bet.eq.0.)then
        ier = 2
        return
     endif
     u(j)=(r(j)-a(j)*u(j-1))/bet
  end do
  !
  !    back substitution
  !
  do j=n-1,1,-1
     u(j)=u(j)-gam(j+1)*u(j+1)
  end do

end subroutine tridag

!-----------------------------------------------------------------------------
subroutine modd_dayno(iyy,imo,idy,iday,j)
  !----------------------------------------------------------------------------
  !  Routine finds day number in a year for given month and day of the month, 
  !  and vice versa
  !
  !   imo: month
  !   idy: day number in the month
  !  iday: day number in the year
  !
  !  When j>0, find day number in a year for given month and day of the month.
  !  When j<0, find month and day of the month for given day number in a year.

  parameter (nm=12)
  integer imv(nm),imv_r(nm),imv_l(nm)
  data imv_r/0,31,59,90,120,151,181,212,243,273,304,334/ !days in each month
  data imv_l/0,31,60,91,121,152,182,213,244,274,305,335/  ! leap year

  !  Determine regular or leap year
  leap=0
  if (mod(iyy,4).eq.0) leap=1
  do i=1,nm
     if (leap.eq.0) then
        imv(i)=imv_r(i)
        iday_max=365
     else
        imv(i)=imv_l(i)
        iday_max=366
     endif
  enddo

  !  Find iday when j>0 and imo,idy when j<0
  if (j.ge.0) then
     iday=idy+imv(imo)
  else
     if (iday.gt.iday_max) then       ! year boundary
        iday=iday-iday_max
        iyy=iyy+1
     endif
     call ilocate(imv,nm,iday,imo)
     idy=iday-imv(imo)
  endif

end subroutine modd_dayno

!!===============================================================================
!!
!SUBROUTINE TRACE1(XIsm,YIsm,ZIsm,DIR,RLIM,R0,IOPT,PARMOD,EXNAME,INNAME,&
!     XFsm,YFsm,ZFsm,XXsm,YYsm,ZZsm,L,iout)
!  
!  !   A modification of TRACE in geopack_2005.f. This routine input sm xi, yi, zi, output 
!  !   iout, the tracing stops at the magnetic equator, and output fieldline points in
!  !   sm coordinates.
!  !
!  !  TRACES A FIELD LINE FROM AN ARBITRARY POINT OF SPACE TO THE EARTH'S
!  !  SURFACE OR TO A MODEL LIMITING BOUNDARY.
!  !
!  !  THE HIGHEST ORDER OF SPHERICAL HARMONICS IN THE MAIN FIELD EXPANSION USED
!  !  IN THE MAPPING IS CALCULATED AUTOMATICALLY. IF INNAME=IGRF_GSM, THEN AN IGRF MODEL
!  !  FIELD WILL BE USED, AND IF INNAME=DIP, A PURE DIPOLE FIELD WILL BE USED.
!  
!  !  IN ANY CASE, BEFORE CALLING TRACE, ONE SHOULD INVOKE RECALC, TO CALCULATE CORRECT
!  !  VALUES OF THE IGRF COEFFICIENTS AND ALL QUANTITIES NEEDED FOR TRANSFORMATIONS
!  !  BETWEEN COORDINATE SYSTEMS INVOLVED IN THIS CALCULATIONS.
!  !
!  !  ALTERNATIVELY, THE SUBROUTINE RECALC CAN BE INVOKED WITH THE DESIRED VALUES OF
!  !  IYEAR AND IDAY (TO SPECIFY THE DIPOLE MOMENT), WHILE THE VALUES OF THE DIPOLE
!  !  TILT ANGLE PSI (IN RADIANS) AND ITS SINE (SPS) AND COSINE (CPS) CAN BE EXPLICITLY
!  !  SPECIFIED AND FORWARDED TO THE COMMON BLOCK GEOPACK1 (11th, 12th, AND 16th ELEMENTS, RESP.)
!  !
!  !------------- INPUT PARAMETERS:
!  !
!  !   XI,YI,ZI - GSM COORDS OF INITIAL POINT (IN EARTH RADII, 1 RE = 6371.2 km),
!  !
!  !   DIR - SIGN OF THE TRACING DIRECTION: IF DIR=1.0 THEN WE MOVE ANTIPARALLEL TO THE
!  !     FIELD VECTOR (E.G. FROM NORTHERN TO SOUTHERN CONJUGATE POINT),
!  !     AND IF DIR=-1.0 THEN THE TRACING GOES IN THE OPPOSITE DIRECTION.
!  !
!  !   R0 -  RADIUS OF A SPHERE (IN RE) FOR WHICH THE FIELD LINE ENDPOINT COORDINATES
!  !     XF,YF,ZF  SHOULD BE CALCULATED.
!  !
!  !   RLIM - UPPER LIMIT OF THE GEOCENTRIC DISTANCE, WHERE THE TRACING IS TERMINATED.
!  !
!  !   IOPT - A MODEL INDEX; CAN BE USED FOR SPECIFYING AN OPTION OF THE EXTERNAL FIELD
!  !       MODEL (E.G., INTERVAL OF THE KP-INDEX). ALTERNATIVELY, ONE CAN USE THE ARRAY
!  !       PARMOD FOR THAT PURPOSE (SEE BELOW); IN THAT CASE IOPT IS JUST A DUMMY PARAMETER.
!  !
!  !   PARMOD -  A 10-ELEMENT ARRAY CONTAINING MODEL PARAMETERS, NEEDED FOR A UNIQUE
!  !      SPECIFICATION OF THE EXTERNAL FIELD. THE CONCRETE MEANING OF THE COMPONENTS
!  !      OF PARMOD DEPENDS ON A SPECIFIC VERSION OF THE EXTERNAL FIELD MODEL.
!  !
!  !   EXNAME - NAME OF A SUBROUTINE PROVIDING COMPONENTS OF THE EXTERNAL MAGNETIC FIELD
!  !    (E.G., T96_01).
!  !   INNAME - NAME OF A SUBROUTINE PROVIDING COMPONENTS OF THE INTERNAL MAGNETIC FIELD
!  !    (EITHER DIP OR IGRF_GSM).
!  !
!  !-------------- OUTPUT PARAMETERS:
!  !
!  !   XF,YF,ZF - GSM COORDS OF THE LAST CALCULATED POINT OF A FIELD LINE
!  !   XX,YY,ZZ - ARRAYS, CONTAINING COORDS OF FIELD LINE POINTS. HERE THEIR MAXIMAL LENGTH WAS
!  !      ASSUMED EQUAL TO 999.
!  !   L - ACTUAL NUMBER OF THE CALCULATED FIELD LINE POINTS. IF L EXCEEDS 999, TRACING
!  !     TERMINATES, AND A WARNING IS DISPLAYED.
!  !
!  !
!  !     LAST MODIFICATION:  MARCH 31, 2003.
!  !
!  !     AUTHOR:  N. A. TSYGANENKO
!  !
!  parameter (i_one=1,m_one=-1)
!  DIMENSION XX(1000),YY(1000),ZZ(1000), PARMOD(10)
!  DIMENSION XXsm(1000),YYsm(1000),ZZsm(1000)
!  COMMON /GEOPACK1/ AA(26),DD,BB(8)
!  EXTERNAL EXNAME,INNAME
!  !
!  iout=0
!  ERR=0.0001
!  L=0
!  !     DS=0.5*DIR
!  DS=0.3*DIR
!  Zrsm=Zism
!  call smgsm(Xism,Yism,Zism,XI,YI,ZI,i_one)
!  X=XI
!  Y=YI
!  Z=ZI
!  DD=DIR
!  AL=0.
!  !
!  !  here we call RHAND just to find out the sign of the radial component of the field
!  !   vector, and to determine the initial direction of the tracing (i.e., either away
!  !   or towards Earth):
!  !
!  CALL RHAND (X,Y,Z,R1,R2,R3,IOPT,PARMOD,EXNAME,INNAME)
!  AD=0.01
!  IF (X*R1+Y*R2+Z*R3.LT.0.) AD=-0.01
!  !
!  !     |AD|=0.01 and its sign follows the rule:
!  ! (1) if DIR=1 (tracing antiparallel to B vector) then the sign of AD is the same as of Br
!  ! (2) if DIR=-1 (tracing parallel to B vector) then the sign of AD is opposite to that of Br
!  !     AD is defined in order to initialize the value of RR (radial distance at previous step):
!  
!  RR=SQRT(X**2+Y**2+Z**2)+AD
!1 L=L+1
!  IF(L.GT.999) GOTO 7
!  XX(L)=X
!  YY(L)=Y
!  ZZ(L)=Z
!  call smgsm(Xsm,Ysm,Zsm,X,Y,Z,m_one)
!  XXsm(L)=Xsm
!  YYsm(L)=Ysm
!  ZZsm(L)=Zsm
!  RYZ=Y**2+Z**2
!  R2=X**2+RYZ
!  R=SQRT(R2)
!  
!  !  check if the line hit the outer tracing boundary; if yes, then terminate
!  !   the tracing (label 8):
!  
!  !     IF (R.GT.RLIM.OR.RYZ.GT.1600.D0.OR.X.GT.20.D0) GOTO 8
!  IF (R.GT.RLIM.OR.RYZ.GT.1600.D0.OR.X.GT.20.D0) then
!     iout=1
!     GOTO 8
!  endif
!  Zsign=Zsm*Zrsm
!  if (zsign.lt.0.) goto 99            ! cross magnetic equator
!  !
!  !  check whether or not the inner tracing boundary was crossed from outside,
!  !  if yes, then calculate the footpoint position by interpolation (go to label 6):
!  !
!  IF (R.LT.R0.AND.RR.GT.R) GOTO 6
!  
!  !  check if (i) we are moving outward, or (ii) we are still sufficiently
!  !    far from Earth (beyond R=5Re); if yes, proceed further:
!  !
!  IF (R.GE.RR.OR.R.GT.5.) GOTO 5
!  
!  !  now we moved closer inward (between R=3 and R=5); go to 3 and begin logging
!  !  previous values of X,Y,Z, to be used in the interpolation (after having
!  !  crossed the inner tracing boundary):
!  
!  IF (R.GE.3.) GOTO 3
!  !
!  !  we entered inside the sphere R=3: to avoid too large steps (and hence inaccurate
!  !  interpolated position of the footpoint), enforce the progressively smaller
!  !  stepsize values as we approach the inner boundary R=R0:
!  !
!  FC=0.2
!  IF(R-R0.LT.0.05) FC=0.05
!  AL=FC*(R-R0+0.2)
!  DS=DIR*AL
!  GOTO 4
!3 DS=DIR
!4 XR=X
!  YR=Y
!  ZR=Z
!5 RR=R
!  Xrsm=Xsm
!  Yrsm=Ysm
!  Zrsm=Zsm
!  CALL STEP (X,Y,Z,DS,ERR,IOPT,PARMOD,EXNAME,INNAME)
!  GOTO 1
!99 R1=-Zsm/(Zrsm-Zsm)
!  Xsm=Xsm-(Xsm-Xrsm)*R1
!  Ysm=Ysm-(Ysm-Yrsm)*R1
!  Zsm=Zsm-(Zsm-Zrsm)*R1
!  call smgsm(Xsm,Ysm,Zsm,X,Y,Z,i_one)
!  goto 8
!  
!  !
!  !  find the footpoint position by interpolating between the current and previous
!  !   field line points:
!  !
!6 R1=(R0-R)/(RR-R)
!  X=X-(X-XR)*R1
!  Y=Y-(Y-YR)*R1
!  Z=Z-(Z-ZR)*R1
!  GOTO 8
!7 WRITE (*,10)
!  L=999
!8 XF=X
!  YF=Y
!  ZF=Z
!  XFsm=Xsm
!  YFsm=Ysm
!  ZFsm=Zsm
!  do i=L,1000
!     XXsm(I)=XFsm
!     YYsm(I)=YFsm
!     ZZsm(I)=ZFsm
!  enddo
!  RETURN
!10 FORMAT(//,1X,'**** COMPUTATIONS IN THE SUBROUTINE TRACE ARE',&
!        ' TERMINATED: THE CURRENT NUMBER OF POINTS EXCEEDED 1000 ****'//)
!END SUBROUTINE TRACE1  
    

!-----------------------------------------------------------------------------
function gammln(xx)
  !---------------------------------------------------------------------------
  !    This is the nature log of gamma function

  !    double precision for recurrences
  !
  !     double precision x,cof,stp,half,one,fpf,tmp,ser
  real cof(6)
  data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,&
       -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
  data half,one,fpf/0.5d0,1.0d0,5.5d0/
  x=xx-one
  tmp=x+fpf
  tmp=(x+half)*log(tmp)-tmp
  ser=one
  do j=1,6
     x=x+one
     ser=ser+cof(j)/x
  end do
  gammln=tmp+log(stp*ser)

end function gammln
!--------------------------------------------------------------------------
function derivative_3pt(x0,x1,x2,f0,f1,f2,xj)
  !--------------------------------------------------------------------------
  ! Function calculates derivative df/dx at xj using 3-point formula 
  ! (R. Burden, and J. Faires, Numerical Analysis, Prindle, 
  ! Weber & Schmidt, 1985).

 der0=f0*(2.*xj-x1-x2)/((x0-x1)*(x0-x2))
 der1=f1*(2.*xj-x0-x2)/((x1-x0)*(x1-x2))
 der2=f2*(2.*xj-x0-x1)/((x2-x0)*(x2-x1))
 derivative_3pt=der0+der1+der2

end function derivative_3pt


!*****************************************************************************
subroutine tsyndipoleSM(imod,iopt,parmod,ps,t,xsm,ysm,zsm,bxsm,bysm,bzsm)
  !****************************************************************************
  !  Routine calculate the total (T96 or T04 external and dipole fields) field
  !  in SM.

  use EGM_ModTsyganenko, ONLY: t96_01, t04_s, dipole, dipole_t96

  parameter (i_1=1,m_1=-1)
  real parmod(10)

  call smgsm(xsm,ysm,zsm,xgsm,ygsm,zgsm,i_1)
  if (imod.eq.1) then 
     call t96_01(iopt,parmod,ps,xgsm,ygsm,zgsm,bxext,byext,bzext)
     call dipole_t96(ps,xgsm,ygsm,zgsm,bxint,byint,bzint)
  endif
  if (imod.eq.2) then
     call t04_s(iopt,parmod,ps,xgsm,ygsm,zgsm,bxext,byext,bzext)
     call dipole(ps,xgsm,ygsm,zgsm,bxint,byint,bzint)
  endif
  bx=bxint+bxext           ! gsm bx
  by=byint+byext           ! gsm by
  bz=bzint+bzext           ! gsm bz
  call smgsm(bxsm,bysm,bzsm,bx,by,bz,m_1)

end subroutine tsyndipoleSM
!------------------------------------------------------------------------------
subroutine rk4(Bfield,imod,iopt,parmod,ps,t,t0,h,x0,xend,xwrk,nd,f,tend)
  !----------------------------------------------------------------------------
  !   *** FOURTH-ORDER RUNGE-KUTTA ***
  !   Solve xend = x0 + fn*h                ! fn is the unit vector of f

  real x0(nd),xend(nd),x00(3),xwrk(4,nd),f(nd),parmod(10)
  external Bfield

  call Bfield(imod,iopt,parmod,ps,t,x0(1),x0(2),x0(3),f(1),f(2),f(3))

  fmag=sqrt(f(1)*f(1)+f(2)*f(2)+f(3)*f(3))
  do i=1,nd
     x00(i)=x0(i)
     xwrk(1,i)=h*f(i)/fmag
     xend(i)=x00(i)+xwrk(1,i)/2.
     x0(i)=xend(i)
  enddo
  call Bfield(imod,iopt,parmod,ps,t,x0(1),x0(2),x0(3),f(1),f(2),f(3))
  fmag=sqrt(f(1)*f(1)+f(2)*f(2)+f(3)*f(3))
  do i=1,nd
     xwrk(2,i)=h*f(i)/fmag
     xend(i)=x00(i)+xwrk(2,i)/2.
     x0(i)=xend(i)
  enddo
  call Bfield(imod,iopt,parmod,ps,t,x0(1),x0(2),x0(3),f(1),f(2),f(3))
  fmag=sqrt(f(1)*f(1)+f(2)*f(2)+f(3)*f(3))
  do i=1,nd
     xwrk(3,i)=h*f(i)/fmag
     xend(i)=x00(i)+xwrk(3,i)
     x0(i)=xend(i)
  enddo
  call Bfield(imod,iopt,parmod,ps,t,x0(1),x0(2),x0(3),f(1),f(2),f(3))
  fmag=sqrt(f(1)*f(1)+f(2)*f(2)+f(3)*f(3))
  do i=1,nd
     xwrk(4,i)=h*f(i)/fmag
  enddo
  do i=1,nd
     xend(i)=x00(i)+(xwrk(1,i)+2.*xwrk(2,i)&
          +2.*xwrk(3,i)+xwrk(4,i))/6.
  enddo
  tend=t0+h

end subroutine rk4

!!=============================================================================                                                                                             
!integer function get_doy(year, mon, day) result(doy)
!
!  implicit none
!
!  integer :: i
!  integer, dimension(1:12) :: dayofmon
!  integer :: year, mon, day
!
!  dayofmon(1) = 31
!  dayofmon(2) = 28
!  dayofmon(3) = 31
!  dayofmon(4) = 30
!  dayofmon(5) = 31
!  dayofmon(6) = 30
!  dayofmon(7) = 31
!  dayofmon(8) = 31
!  dayofmon(9) = 30
!  dayofmon(10) = 31
!  dayofmon(11) = 30
!  dayofmon(12) = 31
!
!  if (mod(year,4).eq.0) dayofmon(2) = dayofmon(1) + 1
!  doy = 0
!  do i = 1, mon-1
!     doy = doy + dayofmon(i)
!  enddo
!  doy = doy + day
!
!end function get_doy
