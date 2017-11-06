 
subroutine cimi_run(delta_t)
  use ModConst,       ONLY: cLightSpeed, cElectronCharge
  use ModCimiInitialize, ONLY: xmm,xk,dphi,dmm,dk,delE,dmu, xjac
  use ModCimi,        ONLY: f2,dt, Time, phot, Ppar_IC, Pressure_IC, &
                            PressurePar_IC,FAC_C, Bmin_C, &
                            OpDrift_, OpBfield_, OpChargeEx_, &
                            OpWaves_, OpStrongDiff_, OpLossCone_, &
                            rbsumLocal, rbsumGlobal, &
                            rcsumLocal,rcsumGlobal,&
                            driftin, driftout, IsStandAlone,&
                            preP,preF,Eje1,UseStrongDiff,&
                            eChangeOperator_VICI,nOperator,&
                            eChangeLocal,eChangeGlobal
  use ModCimiPlanet,  ONLY: re_m, dipmom, Hiono, nspec, amu_I, &
                            dFactor_I,tFactor_I
  use ModCimiTrace,  ONLY: &
       fieldpara, brad=>ro, ftv=>volume, xo,yo,rb,irm,&
       ekev,iba,bo,pp,Have, sinA, vel, alscone, iw2,xmlto, bm
  use ModGmCimi,      ONLY: Den_IC,UseGm
  use ModIeCimi,      ONLY: UseWeimer, pot
  use ModCimiPlot,    ONLY: Cimi_plot, Cimi_plot_fls, &
       			    Cimi_plot_psd,Cimi_plot_vl,Cimi_plot_vp, &
                            cimi_plot_log, cimi_plot_precip, Cimi_plot_Lstar,&
                            DtOutput, DtLogOut,DoSavePlot, &
                            DoSaveFlux, DoSavePSD, DoSaveDrifts, DoSaveLog
  use ModCimiRestart, ONLY: IsRestart
  use ModImTime
  use ModTimeConvert, ONLY: time_real_to_int
  use ModImSat,       ONLY: nImSats, write_im_sat, DoWriteSats, DtSatOut
  use ModCimiGrid,    ONLY: iProc,nProc,iComm,nLonPar,nLonPar_P,nLonBefore_P, &
                            MinLonPar,MaxLonPar,nt,np,neng,npit,nm,nk,dlat,&
                            energy,phi,sinao,xlat,xmlt
  use ModCimiBoundary,ONLY: cimi_set_boundary_mhd, cimi_set_boundary_empirical
  use ModMpi
  use ModWaveDiff,    ONLY: UseWaveDiffusion, ReadDiffCoef, WavePower, & 
                            diffuse_aa, diffuse_EE, diffuse_aE, DiffStartT, &
                            testDiff_aa, testDiff_EE, testDiff_aE, &
                            TimeAeIndex_I, AeIndex_I, interpolate_ae
  use ModCoupleSami,  ONLY: DoCoupleSami
  use DensityTemp,    ONLY: density, simple_plasmasphere
  use ModIndicesInterfaces
  use ModLstar,       ONLY: calc_Lstar1,calc_Lstar2 
  implicit none

  !regular variables
  integer n,nstep
  integer,save :: ib0(nt)
  real delta_t
  real flux(nspec,np,nt,neng,npit),&
       vlEa(nspec,np,nt,neng,npit),vpEa(nspec,np,nt,neng,npit)
  real achar(nspec,np,nt,nm,nk)
  real :: vl(nspec,0:np,nt,nm,nk)=0.0, vp(nspec,np,nt,nm,nk)=0.0, &
       psd(nspec,np,nt,nm,nk), fb(nspec,nt,nm,nk), rc
  integer iLat, iLon, iSpecies, iSat, iOperator
  logical, save :: IsFirstCall =.true.
  real  AE_temp,Kp_temp
  real Lstar_C(np,nt),Lstar_max,Lstarm(np,nt,nk),Lstar_maxm(nk)

  !Vars for mpi passing
  integer ::iSendCount,iM,iK,iLon1,iError,iEnergy,iPit,iRecvLower,iRecvUpper,iPe
  integer,allocatable :: iRecieveCount_P(:),iDisplacement_P(:)
  real :: BufferSend_C(np,nt),BufferRecv_C(np,nt)
  integer :: BufferSend_I(nt),BufferRecv_I(nt)
  integer,allocatable :: iBufferSend_I(:),iBufferRecv_I(:)
  integer :: iStatus_I(MPI_STATUS_SIZE)

  real,   allocatable :: ekevSEND_IIII(:,:,:,:),ekevRECV_IIII(:,:,:,:)
  real,   allocatable :: sinaSEND_III(:,:,:),sinaRECV_III(:,:,:)
  real,   allocatable :: F2SEND_IIIII(:,:,:,:,:),f2RECV_IIIII(:,:,:,:,:)

  integer :: tmp_I(6)
  !----------------------------------------------------------------------------

  if (dt==0) then
     nstep = 0
     dt = 0.0
  else
     nstep=nint(delta_t/dt)
     dt=delta_t/nstep         ! new dt
  endif

  ! Update CurrentTime and iCurrentTime_I
  CurrentTime = StartTime+Time
  call time_real_to_int(CurrentTime,iCurrentTime_I)
  
  ! do field line integration and determine vel, ekev, momentum (pp), etc.
  rc=(re_m+Hiono*1000.)/re_m        ! ionosphere distance in RE`

  call timing_start('cimi_fieldpara')
  call fieldpara(Time,dt,cLightSpeed,cElectronCharge,rc,re_m,xlat,xmlt,phi,xk,&
                 dipmom)
  call timing_stop('cimi_fieldpara')

  ! get Bmin, needs to be passed to GM for anisotropic pressure coupling
  Bmin_C = bo

  ! set the boundary and temperature and density (also sets the interior 
  ! density and temperature for I.C. but that is only if initial_f2 is called)
  if(.not.UseGm) then
     call cimi_set_boundary_empirical
  else
     call cimi_set_boundary_mhd
  endif

  ! Bcast DoWriteSats on firstcall
  if (IsFirstCall .and. nProc > 1) then
     call MPI_bcast(DoWriteSats,1,MPI_LOGICAL,0,iComm,iError)
  endif

  !  read wave models 
  if (IsFirstCall) then
     if (UseWaveDiffusion) call ReadDiffCoef(rc)
  endif
  
  ! setup initial distribution
  if (IsFirstCall .and. .not.IsRestart) then
     !set initial state when no restarting
     call initial_f2(nspec,np,nt,iba,amu_I,vel,xjac,ib0)
     IsFirstCall=.false.
  elseif(IsFirstCall .and. IsRestart) then
     ib0=iba
     IsFirstCall=.false.
  endif
  
  ! Calculate rbsumLocal and Global on first time and get energy
  ! contribution from Bfield change
  call sume_cimi(OpBfield_) 

  ! calculate boundary flux (fb) at the CIMI outer boundary at the equator
  call boundaryIM(nspec,neng,np,nt,nm,nk,iba,irm,amu_I,xjac,energy,vel,fb)

  if (Time.ge.DiffStartT .and. UseWaveDiffusion) then

     ! calculate wave power for the first time; the second time is in the loop
     call interpolate_ae(CurrentTime, AE_temp)

     ! Determines if the simple plasmasphere model needs to be used.
     if ( .not. DoCoupleSami ) then
        call timing_start('cimi_simp_psphere')
        call get_kp(CurrentTime, Kp_temp, iError)        
        call simple_plasmasphere(Kp_temp)
        call timing_stop('cimi_simp_psphere')
     end if

     call WavePower(Time,AE_temp,iba)
  end if

  ! calculate the ionospheric potential (if not using MHD potential)
  if (UseWeimer) then
     call timing_start('set_cimi_potential')
     call set_cimi_potential(CurrentTime,rc) 
     call timing_stop('set_cimi_potential')
  endif
  
  ! save the initial point (needs to be fixed and so not called)
  if (Time == 0.0 .and. nProc == 0 .and. DoSavePlot) then
     call timing_start('cimi_output')
     call cimi_output(np,nt,nm,nk,nspec,neng,npit,iba,ftv,f2,ekev, &
          sinA,energy,sinAo,delE,dmu,amu_I,xjac,pp,xmm,dmm,dk,xlat,dphi, &
          re_m,Hiono,vp,vL,flux,FAC_C,phot,Ppar_IC,Pressure_IC,PressurePar_IC, &
          vlEa,vpEa,psd)
     call timing_stop('cimi_output')
     
     call timing_start('calc_Lstar1')
     call calc_Lstar1(Lstar_C,Lstar_max,rc)
     call timing_stop('calc_Lstar1')

     call timing_start('calc_Lstar2')
     call calc_Lstar2(Lstarm,Lstar_maxm,rc)
     call timing_stop('calc_Lstar2')

     call timing_start('cimi_plot')
     call Cimi_plot(np,nt,xo,yo,Pressure_IC,PressurePar_IC,phot,Ppar_IC,Den_IC,&
          bo,ftv,pot,FAC_C,Time,dt,Lstar_C)
     call timing_stop('cimi_plot')

     call timing_start('cimi_plot_log')
     call cimi_plot_log(Time)
     call timing_stop('cimi_plot_log')

     if (DoSaveFlux) call Cimi_plot_fls(rc,flux,time,Lstar_C,Lstar_max)
     if (DoSavePSD) call Cimi_plot_psd(rc,psd,xmm,xk,time)
     if (DoSaveFlux.or.DoSavePSD) call &
        Cimi_plot_Lstar(rc,xk,time,Lstarm,Lstar_maxm)
     if (DoSaveDrifts) then
        call Cimi_plot_vl(rc,vlEa,time)
        call Cimi_plot_vp(rc,vpEa,time)
     endif
  endif
  

  ! calculate the drift velocity
  call timing_start('cimi_driftV')
  call driftV(nspec,np,nt,nm,nk,irm,re_m,Hiono,dipmom,dphi,xlat, &
       dlat,ekev,pot,vl,vp)
  call timing_stop('cimi_driftV')

  ! calculate the depreciation factor, achar, due to charge exchange loss
  call timing_start('cimi_ceparaIM')
  call ceparaIM(nspec,np,nt,nm,nk,irm,dt,vel,ekev,Have,achar)
  call timing_stop('cimi_ceparaIM')

  if (UseStrongDiff) then
     ! Calculate the strong diffusion lifetime for electrons
     call timing_start('cimi_StDiTime')
     call StDiTime(dt,vel,ftv,rc,re_m,dipmom,iba)
     call timing_stop('cimi_StDiTime')
  endif
  
  ! time loop
  do n=1,nstep

     call timing_start('cimi_driftIM')
     call driftIM(iw2,nspec,np,nt,nm,nk,dt,dlat,dphi,brad,rb,vl,vp, &
          fb,f2,driftin,driftout,ib0)
     call sume_cimi(OpDrift_)
     call timing_stop('cimi_driftIM')
     
     call timing_start('cimi_charexchange')
     call charexchangeIM(np,nt,nm,nk,nspec,iba,achar,f2)
     call sume_cimi(OpChargeEx_)
     call timing_stop('cimi_charexchange') 
     
     if (Time.ge.DiffStartT .and. UseWaveDiffusion) then

        call timing_start('cimi_WaveDiffusion')
        
        call timing_start('cimi_Diffuse_aa')
        call diffuse_aa(f2,dt,xjac,iba,iw2)
        call timing_stop('cimi_Diffuse_aa')

        call timing_start('cimi_Diffuse_EE')
        call diffuse_EE(f2,dt,xmm,xjac,iw2,iba)
        call timing_stop('cimi_Diffuse_EE')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!!!                                                  !!!
!!!     THIS NEEDS TO BE UNCOMMENTED ONCE            !!!
!!!     CROSS-DIFFUSION IS FIXED.  :::COLIN:::       !!!
!!!                                                  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
!!$        call timing_start('cimi_Diffuse_aE')
!!$        call diffuse_aE(f2,dt,xjac,iw2,iba,Time)
!!$        call timing_stop('cimi_Diffuse_aE')

        call sume_cimi(OpWaves_)
        call timing_stop('cimi_WaveDiffusion')
        
        call interpolate_ae(CurrentTime, AE_temp)

        if ( .not. DoCoupleSami ) then
           call timing_start('cimi_simp_psphere')
           call get_kp(CurrentTime, Kp_temp, iError)        
           call simple_plasmasphere(Kp_temp)
           call timing_stop('cimi_simp_psphere')
        end if

        call WavePower(Time,AE_temp,iba)
     endif

     if (UseStrongDiff) then
        call timing_start('cimi_StrongDiff')
        call StrongDiff(iba)
        call sume_cimi(OpStrongDiff_)
        call timing_stop('cimi_StrongDiff')
     endif
     
     call timing_start('cimi_lossconeIM')
     call lossconeIM(np,nt,nm,nk,nspec,iba,alscone,f2)
     call sume_cimi(OpLossCone_)
     call timing_stop('cimi_lossconeIM')
     
     Time = Time+dt
     ! Update CurrentTime and iCurrentTime_I
     CurrentTime = StartTime+Time
     call time_real_to_int(CurrentTime,iCurrentTime_I)
  enddo

  ! calculate precipitations accumulated over some time interval
  ! (DtLogOut in this case)
  if( (floor((Time+1.0e-5)/DtLogOut))/=&
       floor((Time+1.0e-5-delta_t)/DtLogOut)) then
     call timing_start('cimi_precip_calc')
     call cimi_precip_calc(rc,DtLogOut)
     call timing_stop('cimi_precip_calc')
  endif

  call timing_start('cimi_output')
  call cimi_output(np,nt,nm,nk,nspec,neng,npit,iba,ftv,f2,ekev, &
       sinA,energy,sinAo,delE,dmu,amu_I,xjac,pp,xmm,dmm,dk,xlat,dphi, &
       re_m,Hiono,vl,vp,flux,FAC_C,phot,Ppar_IC,Pressure_IC,PressurePar_IC,&
       vlEa,vpEa,psd)
  call timing_stop('cimi_output')
  
  ! When nProc >1 consolodate: phot, Ppar_IC, Pressure_IC,
  ! PressurePar_IC, fac and iba on iProc 0
  if (nProc>1) then    
     if (.not.allocated(iRecieveCount_P)) &
          allocate(iRecieveCount_P(nProc), iDisplacement_P(nProc))       
     !Gather to root
     iSendCount = np*nLonPar
     iRecieveCount_P=np*nLonPar_P
     iDisplacement_P = np*nLonBefore_P
     BufferSend_C(:,:) = FAC_C(:,:) 
     call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
          MPI_REAL, BufferRecv_C, iRecieveCount_P, iDisplacement_P, &
          MPI_REAL, 0, iComm, iError)
     if (iProc==0) FAC_C(:,:)=BufferRecv_C(:,:)

     do iSpecies=1,nspec
        BufferSend_C(:,:)=Pressure_IC(iSpecies,:,:)
        call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
             MPI_REAL, BufferRecv_C, iRecieveCount_P, iDisplacement_P, &
             MPI_REAL, 0, iComm, iError)
        if (iProc==0) Pressure_IC(iSpecies,:,:)=BufferRecv_C(:,:)

        BufferSend_C(:,:)=PressurePar_IC(iSpecies,:,:)
        call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
             MPI_REAL, BufferRecv_C, iRecieveCount_P, iDisplacement_P, &
             MPI_REAL, 0, iComm, iError)
        if (iProc==0) PressurePar_IC(iSpecies,:,:)=BufferRecv_C(:,:)

        BufferSend_C(:,:)=phot(iSpecies,:,:)
        call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
             MPI_REAL, BufferRecv_C, iRecieveCount_P, iDisplacement_P, &
             MPI_REAL, 0, iComm, iError)
        if (iProc==0) phot(iSpecies,:,:)=BufferRecv_C(:,:)

        BufferSend_C(:,:)=Ppar_IC(iSpecies,:,:)
        call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
             MPI_REAL, BufferRecv_C, iRecieveCount_P, iDisplacement_P, &
             MPI_REAL, 0, iComm, iError)
        if (iProc==0) Ppar_IC(iSpecies,:,:)=BufferRecv_C(:,:)

        BufferSend_C(:,:)=Den_IC(iSpecies,:,:)
        call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
             MPI_REAL, BufferRecv_C, iRecieveCount_P, iDisplacement_P, &
             MPI_REAL, 0, iComm, iError)
        if (iProc==0) Den_IC(iSpecies,:,:)=BufferRecv_C(:,:)
     enddo

     BufferSend_I(:) = iba(:)
     call MPI_GATHERV(BufferSend_I(MinLonPar:MaxLonPar), nLonPar, &
          MPI_INTEGER, BufferRecv_I, nLonPar_P, nLonBefore_P, &
          MPI_INTEGER, 0, iComm, iError)
     if (iProc==0) iba(:)=BufferRecv_I(:)

     BufferSend_C(:,:) = Bmin_C(:,:)
     call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
          MPI_REAL, BufferRecv_C, iRecieveCount_P, iDisplacement_P, &
          MPI_REAL, 0, iComm, iError)
     if (iProc==0) Bmin_C(:,:)=BufferRecv_C(:,:)
     
  endif
  
  ! On processor O, gather info and save plots
  ! When time to write output, consolodate xo,yo,flux,pot,ftv, bo, and irm 
  ! on iProc 0
  if (nProc>1 .and. DoSavePlot  .and.&
       (floor((Time+1.0e-5)/DtOutput))/=&
       floor((Time+1.0e-5-delta_t)/DtOutput)) then
     
!     call MPI_GATHERV(pot(:,MinLonPar:MaxLonPar), iSendCount, MPI_REAL, &
!          pot, iRecieveCount_P, iDisplacement_P, MPI_REAL, &
!          0, iComm, iError)

     BufferSend_C(:,:)=ftv(:,:)
     call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
          MPI_REAL, BufferRecv_C, iRecieveCount_P, iDisplacement_P, &
          MPI_REAL, 0, iComm, iError)
     if (iProc==0) ftv(:,:)=BufferRecv_C(:,:)
     BufferSend_C(:,:)=xo(:,:)
     call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
          MPI_REAL, BufferRecv_C, iRecieveCount_P, iDisplacement_P, &
          MPI_REAL, 0, iComm, iError)
     if (iProc==0) xo(:,:)=BufferRecv_C(:,:)
     BufferSend_C(:,:)=yo(:,:)
     call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
          MPI_REAL, BufferRecv_C, iRecieveCount_P, iDisplacement_P, &
          MPI_REAL, 0, iComm, iError)
     if (iProc==0) yo(:,:)=BufferRecv_C(:,:)
     BufferSend_C(:,:)=bo(:,:)
     call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
          MPI_REAL, BufferRecv_C, iRecieveCount_P, iDisplacement_P, &
          MPI_REAL, 0, iComm, iError)
     if (iProc==0) bo(:,:)=BufferRecv_C(:,:)
     BufferSend_C(:,:)=xmlto(:,:)
     call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
          MPI_REAL, BufferRecv_C, iRecieveCount_P, iDisplacement_P, &
          MPI_REAL, 0, iComm, iError)
     if (iProc==0) xmlto(:,:)=BufferRecv_C(:,:)
     BufferSend_C(:,:)=brad(:,:)
     call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
          MPI_REAL, BufferRecv_C, iRecieveCount_P, iDisplacement_P, &
          MPI_REAL, 0, iComm, iError)
     if (iProc==0) brad(:,:)=BufferRecv_C(:,:)
     if ( .not. DoCoupleSami ) then

        BufferSend_C(:,:) = density(:,:)
        call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
             MPI_REAL, BufferRecv_C, iRecieveCount_P, iDisplacement_P, &
             MPI_REAL, 0, iComm, iError)
        if (iProc==0) density(:,:) = BufferRecv_C(:,:)
        
     end if
     BufferSend_I(:)=irm(:)
     call MPI_GATHERV(BufferSend_I(MinLonPar:MaxLonPar), nLonPar, &
          MPI_INTEGER, BufferRecv_I, nLonPar_P, nLonBefore_P, &
          MPI_INTEGER, 0, iComm, iError)
     if (iProc==0) irm(:)=BufferRecv_I(:)
  elseif (nProc > 1 .and. DoWriteSats .and. &
       (floor((Time+1.0e-5)/DtSatOut))/=&
       floor((Time+1.0e-5-delta_t)/DtSatOut)) then
     BufferSend_C(:,:)=bo(:,:)
     call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
          MPI_REAL, BufferRecv_C, iRecieveCount_P, iDisplacement_P, &
          MPI_REAL, 0, iComm, iError)     
     if (iProc==0) bo(:,:)=BufferRecv_C(:,:)
  endif


  if (nProc>1 .and. ((DoSavePlot  .and.&
       (floor((Time+1.0e-5)/DtOutput))/=&
       floor((Time+1.0e-5-delta_t)/DtOutput)) .or.&
       ((floor((Time+1.0e-5)/DtSatOut))/=&
       floor((Time+1.0e-5-delta_t)/DtSatOut).and.DoWriteSats))) then
     do  iSpecies=1,nspec
        do iM=1,nm
           do iK=1,nk

              BufferSend_C(:,:)=psd(iSpecies,:,:,im,iK)
              call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar),iSendCount, &
                   MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P, &
                   MPI_REAL, 0, iComm, iError)
              if (iProc==0) psd(iSpecies,:,:,im,ik)=BufferRecv_C(:,:)

           end do
        end do
        
        do iEnergy=1,neng
           do iPit=1,nPit

              BufferSend_C(:,:)=flux(iSpecies,:,:,iEnergy,iPit)
              call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar),iSendCount, &
                   MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P, &
                   MPI_REAL, 0, iComm, iError)
              if (iProc==0) flux(iSpecies,:,:,iEnergy,iPit)=BufferRecv_C(:,:)

              BufferSend_C(:,:)=vlEa(iSpecies,:,:,iEnergy,iPit)
              call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar),iSendCount, &
                   MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P, &
                   MPI_REAL, 0, iComm, iError)
              if (iProc==0) vlEa(iSpecies,:,:,iEnergy,iPit)=BufferRecv_C(:,:)

              BufferSend_C(:,:)=vpEa(iSpecies,:,:,iEnergy,iPit)
              call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar),iSendCount, &
                   MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P, &
                   MPI_REAL, 0, iComm, iError)
              if (iProc==0) vpEa(iSpecies,:,:,iEnergy,iPit)=BufferRecv_C(:,:)

           enddo
        enddo
     enddo
  endif

  if (nProc>1 .and. (DoSaveLog .and. &
          (floor((Time+1.0e-5)/DtLogOut))/=&
          floor((Time+1.0e-5-delta_t)/DtLogOut))) then
      do  iSpecies=1,nspec
         do iEnergy=1,neng+2
            BufferSend_C(:,:)=preP(iSpecies,:,:,iEnergy)
            call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar),iSendCount, &
                 MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P, &
                 MPI_REAL, 0, iComm, iError)
            if (iProc==0) preP(iSpecies,:,:,iEnergy)=BufferRecv_C(:,:) 
            BufferSend_C(:,:)=preF(iSpecies,:,:,iEnergy)  
            call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar),iSendCount, &
                 MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P, &
                 MPI_REAL, 0, iComm, iError)
            if (iProc==0) preF(iSpecies,:,:,iEnergy)=BufferRecv_C(:,:)
         enddo ! Do loop over iEnergy
         
         BufferSend_C(:,:)=Eje1(iSpecies,:,:)
         call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar),iSendCount, &
                  MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P, &
                  MPI_REAL, 0, iComm, iError)
         if (iProc==0) Eje1(iSpecies,:,:)=BufferRecv_C(:,:)
      enddo !Do loop over Species
      
   endif

   ! Gather bm to root for multiple processes for Lstar2 calculation.
   if (nProc > 1 .and. &
        DoSavePlot .and. &
        (floor((Time+1.0e-5)/DtOutput))/=&
        floor((Time+1.0e-5-delta_t)/DtOutput)) then

      do iK = 1,nk
         BufferSend_C(:,:)=bm(:,:,iK)
         call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar),iSendCount, &
              MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P, &
              MPI_REAL, 0, iComm, iError)
         if (iProc==0) bm(:,:,iK)=BufferRecv_C(:,:)
      end do      
   endif
   
  !Gather to root
  !    iSendCount = np*nLonPar
  !    iRecieveCount_P=np*nLonPar_P
  !    iDisplacement_P = np*nLonBefore_P
  !    call MPI_GATHERV(ftv(:,MinLonPar:MaxLonPar), iSendCount, MPI_REAL, &
  !         ftv, iRecieveCount_P, iDisplacement_P, MPI_REAL, &
  !            0, iComm, iError)
  if (iProc == 0) then
     ! do main plotting
     if (DoSavePlot.and.&
          (floor((Time+1.0e-5)/DtOutput))/=&
          floor((Time+1.0e-5-delta_t)/DtOutput)) then
        call timing_start('calc_Lstar1')
        call calc_Lstar1(Lstar_C,Lstar_max,rc)
        call timing_stop('calc_Lstar1')
        
        call timing_start('calc_Lstar2')
        call calc_Lstar2(Lstarm,Lstar_maxm,rc)
        call timing_stop('calc_Lstar2')

        call timing_start('cimi_plot')
        call Cimi_plot(np,nt,xo,yo,Pressure_IC,PressurePar_IC,phot,&
             Ppar_IC,Den_IC,bo,ftv,pot,FAC_C,Time,dt,Lstar_C)
        call timing_stop('cimi_plot')

        if (DoSaveFlux) call Cimi_plot_fls(rc,flux,time,Lstar_C,Lstar_max)
        if (DoSavePSD) call Cimi_plot_psd(rc,psd,xmm,xk,time)
        if (DoSaveDrifts) then
           call Cimi_plot_vl(rc,vlEa,time)
           call Cimi_plot_vp(rc,vpEa,time)
        endif
     endif
  
     ! Write Sat Output
     if(DoWriteSats .and. DoSavePlot .and. &
          (floor((Time+1.0e-5)/DtSatOut))/=&
          floor((Time+1.0e-5-delta_t)/DtSatOut))then
        do iSat=1,nImSats
           call timing_start('cimi_write_im_sat')
           call write_im_sat(iSat,np,nt,neng,npit,flux)
           call timing_stop('cimi_write_im_sat')
        enddo
     endif

     ! Write Logfile
     if(DoSaveLog .and. &
          (floor((Time+1.0e-5)/DtLogOut))/=&
          floor((Time+1.0e-5-delta_t)/DtLogOut))then
        call timing_start('cimi_plot_log')
        call cimi_plot_log(Time)
        call timing_stop('cimi_plot_log')
     endif

     ! Write precipitation file
     if(DoSaveLog .and. &
          (floor((Time+1.0e-5)/DtLogOut))/=&
          floor((Time+1.0e-5-delta_t)/DtLogOut))then
        call timing_start('cimi_plot_precip')
        call cimi_plot_precip(rc,Time)
        call timing_stop('cimi_plot_precip')
     endif

  endif

end subroutine Cimi_run

!-----------------------------------------------------------------------------
subroutine cimi_init
  !---------------------------------------------------------------------------
  ! Routine does CIMI initialization: fill arrays
  !
  ! Input: np,nt,neng,npit,nspec,re_m,dipmom,Hiono
  ! Output: xlat,xmlt,energy,sinAo (through augments)
  !         xmm1,xk1,phi1,dlat1,dphi1,dmm1,dk1,delE1,dmu1,xjac,d4,amu (through 
  !         common block cinitialization

  use ModPlanetConst, ONLY: Earth_,DipoleStrengthPlanet_I,rPlanet_I
  use ModConst,       ONLY: cElectronCharge
  use ModNumConst,    ONLY: cDegToRad,cRadToDeg,cPi
  use ModCimiPlanet,  ONLY: re_m, dipmom, Hiono, amu_I, nspec
  use ModCimiInitialize
  use ModCimiRestart, ONLY: IsRestart, cimi_read_restart
  use ModImTime
  use ModCimiGrid,    ONLY: iProcLeft, iProcRight, iLonLeft, iLonRight, &
       d4Element_C, MinIonEnergy, MaxIonEnergy, iProc,nLonPar, nt, nProc,&
       nLonBefore_P, nLonPar_P, MinLonPar, MaxLonPar, nLonPar_P, &
       iLonMidnight, iProcMidnight, nLonPar_P, np, xlatr, sinao, npit, energy, &
       dLat,ebound, iComm, phi, xlat,xlat_data,xmlt 
  use ModTimeConvert, ONLY: time_int_to_real,time_real_to_int
  use ModMpi

  implicit none

  integer i,n,k,m,iPe, iError

  real aloga,eratio
  real sina0,sina1
  real rw,rw1,rsi,rs1
  real xjac1,sqrtm
  real d2,energy_ion(1:neng),energy_ele(1:neng)


  ! Set up proc distribution
  if (iProc < mod(nt,nProc))then
     nLonPar=(nt+nProc-1)/nProc
  else
     nLonPar=nt/nProc
  endif

  ! Figure out the min and max longitude range for each proc
  if (.not.allocated(nLonBefore_P)) allocate(nLonBefore_P(0:nProc-1))
  if (.not.allocated(nLonPar_P))    allocate(nLonPar_P(0:nProc-1))
  call MPI_allgather(nLonPar,1,MPI_INTEGER,nLonPar_P,1,MPI_INTEGER,iComm,iError)
  nLonBefore_P(0) = 0
  do iPe = 1, nProc - 1
     nLonBefore_P(iPe) = sum(nLonPar_P(0:iPe-1))
  end do

  if (iProc == 0) then
     MinLonPar=1
     MaxLonPar=nLonPar
  else
     MinLonPar=sum(nLonPar_P(0:iProc-1))+1
     MaxLonPar=minLonPar+nLonPar-1
  endif

  !Define neighbors and ghost cell indicies
  iProcLeft=iProc-1
  iLonLeft=MinLonPar-1
  if (iProcLeft < 0) then 
     iProcLeft=nProc-1
     iLonLeft = nt
  endif
  iProcRight=iProc+1
  iLonRight=MaxLonPar+1
  if (iProcRight == nProc) then
     iProcRight=0
     iLonRight=1
  endif

  !Define and iLonMidnightiProcMidnight, needed for setting iw2 in fieldpara
  iLonMidnight=nt/2
  if (nProc>1) then
     PROCLIST: do iPe=0,nProc-1
        if (nLonBefore_P(iPe)<iLonMidnight .and. &
             nLonBefore_P(iPe)+nLonPar_P(iPe)>=iLonMidnight)then
           iProcMidnight=iPe
        endif
     enddo PROCLIST
  else
     iProcMidnight=0
  endif
  ! Set start time

  call time_int_to_real(iStartTime_I,CurrentTime)
  StartTime=CurrentTime

  ! Define constants
  re_m = rPlanet_I(Earth_)                            ! earth's radius (m)
  dipmom=abs(DipoleStrengthPlanet_I(Earth_)*re_m**3)  ! earth's dipole moment
  
  ! CIMI xlat and xmlt grids
  do i=1,np
     xlat(i)=xlat_data(i)
     dlat(i)=0.5*(xlat_data(i+1)-xlat_data(i-1))*cDegToRad    ! dlat in radian
  enddo
  xlatr=xlat*cDegToRad  
  dphi=2.*cPi/nt
  do i=1,nt
     phi(i)=(i-1)*dphi
     xmlt(i)=mod(phi(i)*12.0/cPi + 12.0,24.0)   
  enddo

  ! CIMI output grids: energy, sinAo, delE1, dmu1

!!$  Replacing old grids in CIMI with those of Mei-Ching's standalone CIMI.
!!$  -Colin 07/23/2015.
  
!!$  energy_ion=(/1.0000,1.6795,2.8209,4.7378,7.9574,13.365, &
!!$       22.447,37.701,63.320,106.35,178.62,300.00/)
!!$  
!!$  energy_ele=10.*(/1.0000,1.6795,2.8209,4.7378,7.9574,13.365, &
!!$       22.447,37.701,63.320,106.35,178.62,300.00/)

!  delE=0.5243*energy_ion

!!$  sinAo=(/0.010021,0.030708,0.062026,0.086108,0.16073,0.27682, &
!!$       0.430830,0.601490,0.753790,0.863790,0.94890,0.98827/)
!!$  dmu=(/0.000207365,0.000868320,0.00167125,0.00489855,0.0165792,0.0404637, &
!!$       0.078819500,0.121098000,0.14729600,0.16555900,0.1738560,0.2486830/)

  energy_ion(1)=MinIonEnergy
  energy_ion(neng)=MaxIonEnergy
  aloga=log10(energy_ion(neng)/energy_ion(1))/(neng-1)
  eratio=10.**aloga
  do k=2,neng-1
     energy_ion(k)=energy_ion(k-1)*eratio
  enddo  
  energy_ele(1:neng)=10.*energy_ion(1:neng)

  sinAo=(/0.009417,0.019070,0.037105,0.069562,0.122536,0.199229, &
       0.296114,0.405087,0.521204,0.638785,0.750495,0.843570, &
       0.910858,0.952661,0.975754,0.988485,0.995792,0.998703/)

  do m=1,npit
     if (m.eq.1) sina0=0.
     if (m.gt.1) sina0=0.5*(sinAo(m)+sinAo(m-1))
     if (m.lt.npit) sina1=0.5*(sinAo(m)+sinAo(m+1))
     if (m.eq.npit) sina1=1.
     dmu(m)=sqrt(1.-sina0*sina0)-sqrt(1.-sina1*sina1)
  enddo
    
!  do k=2,neng
!     Ebound(k)=sqrt(energy(k-1)*energy(k))
!  enddo
!  Ebound(1) = energy(1)**2.0/Ebound(2)
!  Ebound(neng+1)=energy(neng)**2.0/Ebound(neng)


! setup energies:
  do n=1,nspec-1
  energy(n,1:neng) = energy_ion(1:neng)
  delE(n,1:neng)   = 0.5243*energy_ion(1:neng)
  enddo
  energy(nspec,1:neng)=energy_ele(1:neng)
  delE(nspec,1:neng)   = 0.5243*energy_ele(1:neng)

! From CIMI:
! Setup Ebound
  do n=1,nspec
     do k=1,neng-1
         Ebound(n,k)=sqrt(energy(n,k)*energy(n,k+1))
     enddo
     Ebound(n,neng)=energy(n,neng)*energy(n,neng)/Ebound(n,neng-1)
  enddo

  ! CIMI magnetic moment, xmm1
 do n=1,nspec
  xmm(n,1)=energy(n,1)*cElectronCharge/(dipmom/(2*re_m)**3.0)
!  dmm(n,1)=xmm(n,1)*2.              
  rw=1.55 
  rw1=(rw-1.)/sqrt(rw)
  xmm(n,0)=xmm(n,1)/rw                      
  !do i=2,nm                    
  !   dmm(n,i)=dmm(n,1)*rw**(i-1)           
  !   xmm(n,i)=xmm(n,i-1)+0.5*(dmm(n,i-1)+dmm(n,i))
  !enddo
  do i=1,nm            ! This setup makes xmm(k+0.5)=sqrt(xmm(k)*xmm(k+1))
         xmm(n,i)=xmm(n,i-1)*rw
         dmm(n,i)=xmm(n,i)*rw1
  enddo
      xmm(n,nm+1)=xmm(n,nm)*rw
 enddo

  ! CIMI K, xk
  rsi=1.47
  xk(1)=40.*rsi
  rs1=(rsi-1.)/sqrt(rsi) ! in following sutup: xk(i+0.5)=sqrt(xk(i)*xk(i+1))
  do i=1,nk
     if (i.gt.1) xk(i)=xk(i-1)*rsi
     dk(i)=xk(i)*rs1                 
  enddo
  
  xk(0)=40.
  xk(nk+1)=xk(nk)*rsi

  ! Calculate Jacobian, xjac
  do n=1,nspec 
     xjac1=4.*sqrt(2.)*cPi*(1.673e-27*amu_I(n))*dipmom/(re_m+Hiono*1000.)
     sqrtm=sqrt(1.673e-27*amu_I(n))
     do i=1,np
        do k=1,nm
           xjac(n,i,k)=xjac1*sin(2.*xlatr(i))*sqrt(xmm(n,k))*sqrtm
        enddo
     enddo
  enddo

  ! Calculate d4Element_C: dlat*dphi*dmm*dk
    do n=1,nspec
      do i=1,np
         d2=dlat(i)*dphi
         do k=1,nm
            do m=1,nk
               d4Element_C(n,i,k,m)=d2*dmm(n,k)*dk(m)
            enddo
         enddo
      enddo
    enddo

  if(IsRestart) then
     !set initial state when restarting
     call cimi_read_restart
  endif

end subroutine cimi_init

!-------------------------------------------------------------------------------
subroutine initial_f2(nspec,np,nt,iba,amu_I,vel,xjac,ib0)
  !-----------------------------------------------------------------------------
  ! Routine setup initial distribution.
  ! 
  ! Input: nspec,np,nt,iba,Den_IC,Temp_IC,amu,vel,xjac
  ! Output: ib0,f2,rbsum,xleb,xled,xlel,xlee,xles,driftin,driftout
  !         (through common block cinitial_f2)
  use ModIoUnit, ONLY: UnitTmp_
  use ModGmCimi, ONLY: Den_IC, Temp_IC, Temppar_IC, DoAnisoPressureGMCoupling
  use ModCimi,   ONLY: f2, nOperator, driftin, driftout, &
       eChangeOperator_VICI, echangeLocal, eChangeGlobal, &
       pChangeOperator_VICI, eTimeAccumult_ICI, pTimeAccumult_ICI, &
       rbsumLocal, rbsumGlobal, rcsumLocal, rcsumGlobal
  use ModCimiInitialize,   ONLY: IsEmptyInitial, IsDataInitial, IsRBSPData, &
       IsGmInitial
  use ModCimiGrid,ONLY: nm,nk,MinLonPar,MaxLonPar,iProc,nProc,iComm,d4Element_C,neng
  use ModCimiTrace, ONLY: sinA,ro, ekev,pp,iw2,irm
  use ModMpi
  use ModWaveDiff, ONLY:  testDiff_aa,testDiff_EE,testDiff_aE

  implicit none

  integer,parameter :: np1=51,nt1=48,nspec1=1  
  !integer,parameter :: nm=35,nk=28 ! dimension of CIMI magnetic moment and K
 
  integer nspec,np,nt,iba(nt),ib0(nt),n,j,i,k,m, iError
  real amu_I(nspec),vel(nspec,np,nt,nm,nk)
  real velperp2, velpar2
  real xjac(nspec,np,nm),pi,xmass,chmass,f21,vtchm
  real Tempperp_IC(nspec,np,nt)
  real xleb(nspec),xled(nspec),xlel(nspec),xlee(nspec),xles(nspec)

  ! Variables needed for data initialization 
  integer :: il, ie, iunit
  real, allocatable :: roi(:), ei(:), fi(:,:)
  real :: roii, e1,x, fluxi,psd2,etemp
  
  character(11) :: NameFile='quiet_x.fin'
  character(5) :: FilePrefix='xxxxx'
  pi=acos(-1.)

  ib0=iba
  f2=0.

  if (IsEmptyInitial) then
     ! Set initial f2 to a small number
     f2(:,:,:,:,:)=1.0e-40
  elseif(IsGmInitial) then
     ! Set initial f2 based on Maxwellian or bi-Maxwellian
     if(DoAnisoPressureGMCoupling) &
          Tempperp_IC(:,:,:) = (3*Temp_IC(:,:,:) - Temppar_IC(:,:,:))/2.
     do n=1,nspec
        xmass=amu_I(n)*1.673e-27
        chmass=1.6e-19/xmass
        do j=MinLonPar,MaxLonPar
           do i=1,iba(j)
              if(DoAnisoPressureGMCoupling)then
                 f21=Den_IC(n,i,j)/(2.*pi*xmass*Temppar_IC(n,i,j)*1.6e-19)**0.5 &
                      /(2.*pi*xmass*Tempperp_IC(n,i,j)*1.6e-19)
              else
                 f21=Den_IC(n,i,j)/(2.*pi*xmass*Temp_IC(n,i,j)*1.6e-19)**1.5
              end if
              do k=1,nm
                 do m=1,nk
                    if(DoAnisoPressureGMCoupling)then
                       velperp2 = (vel(n,i,j,k,m)*sinA(i,j,m))**2
                       velpar2 = vel(n,i,j,k,m)**2 - velperp2
                       vtchm = -velpar2/(2*Temppar_IC(n,i,j)*chmass) &
                            -velperp2/(2*Tempperp_IC(n,i,j)*chmass)
                    else                    
                       vtchm = -vel(n,i,j,k,m)**2/(2*Temp_IC(n,i,j)*chmass)
                    end if
                    f2(n,i,j,k,m)=xjac(n,i,k)*f21*exp(vtchm)
                 end do
              end do
           end do
        end do
     end do
  elseif(IsDataInitial) then
     do n=1,nspec
        if (IsRBSPData) then
           FilePrefix='RBSP'
        else
           FilePrefix='quiet'
        endif
        !set the file name, open it and read it
        SELECT CASE (n)
        CASE (1)
           if (n==nspec) then
              NameFile=TRIM(FilePrefix) // '_e.fin'
           else
              NameFile=TRIM(FilePrefix) // '_h.fin'
           endif
        CASE (2)
           if (n==nspec) then
              NameFile=TRIM(FilePrefix) // '_e.fin'
           else
              NameFile=TRIM(FilePrefix) // '_o.fin'
           endif
        CASE DEFAULT
           NameFile=TRIM(FilePrefix) // '_e.fin'
        END SELECT
        open(unit=UnitTmp_,file='IM/'//NameFile,status='old')
        read(UnitTmp_,*) il,ie
        allocate (roi(il),ei(ie),fi(il,ie))
        read(UnitTmp_,*) iunit   ! 1=flux in (cm2 s sr keV)^-1, 2=in (cm2 s MeV)^-1
        read(UnitTmp_,*) roi
        read(UnitTmp_,*) ei      ! ei in keV
        read(UnitTmp_,*) fi
        close(UnitTmp_)
        if(iunit.eq.2) fi(:,:)=fi(:,:)/4./pi/1000. !con.To(cm^2 s sr keV)^-1\

        ei(:)=log10(ei(:))                      ! take log of ei         
        fi(:,:)=log10(fi(:,:))                  ! take log of fi
        
        !interpolate data from quiet.fin files to CIMI grid
        do j=MinLonPar,MaxLonPar
           do i=1,irm(j)
              roii=ro(i,j)
              do m=1,nk
                 do k=1,iw2(n,m)
                    e1=log10(ekev(n,i,j,k,m)) 
                   ! if (e1.le.ei(ie)) then
                    if (e1.ge.ei(ie)) e1=ei(ie)    ! flat dist at high E
                       if (e1.lt.ei(1)) e1=ei(1)    ! flat dist. at low E
                       if (roii.lt.roi(1)) roii=roi(1) ! flat dist @ lowL
                       if (roii.gt.roi(il)) roii=roi(il) ! flat @ high L
                       call lintp2IM(roi,ei,fi,il,ie,roii,e1,x)
                       fluxi=10.**x          ! flux in (cm^2 s sr keV)^-1
                       psd2=fluxi/(1.6e19*pp(n,i,j,k,m))/pp(n,i,j,k,m)
                   !   if (testDiff_aa) psd2=psd2*sinA(i,j,m)*sinA(i,j,m)  !  
                       if (testDiff_EE.or.testDiff_aE)  psd2=1.
                       f2(n,i,j,k,m)=psd2*xjac(n,i,k)*1.e20*1.e19 
                   ! endif
                 enddo                            ! end of k loop
              enddo                               ! end of m loop
              
           enddo                                  ! end of i loop
        enddo                                     ! end of j loop
        deallocate (roi,ei,fi)
        !f2(:,1,:,:,:)=f2(:,2,:,:,:)
     enddo                                        ! end of n loop
  end if

! Setup variables for energy gain/loss from each process
  eChangeOperator_VICI(1:nspec,1:np,1:nt,1:neng+2,1:nOperator)=0.0
  pChangeOperator_VICI(1:nspec,1:np,1:nt,1:neng+2,1:nOperator)=0.0
  eTimeAccumult_ICI(1:nspec,1:np,1:nt,1:neng+2) = 0.0
  pTimeAccumult_ICI(1:nspec,1:np,1:nt,1:neng+2) = 0.0

  eChangeLocal(1:nspec, 1:nOperator) = 0.
  eChangeGlobal(1:nspec, 1:nOperator) = 0.

  ! Total energy for the simulation domain
  rbsumLocal(1:nspec) = 0. 
  rbsumGlobal(1:nspec) = 0.

  ! Total energy for L < 6.6 R_E
  rcsumLocal(1:nspec) = 0.
  rcsumGlobal(1:nspec) = 0.

  driftin(1:nspec)=0.      ! energy gain due injection
  driftout(1:nspec)=0.     ! energy loss due drift-out loss

end subroutine initial_f2


!-------------------------------------------------------------------------------
subroutine boundaryIM(nspec,neng,np,nt,nm,nk,iba,irm,amu_I,xjac,energy,vel,fb)
  !-----------------------------------------------------------------------------
  !  *Relativistic version of boundary PSD.*
  ! Routine setup the boundary distribution for the CIMI. Distribution
  ! at the boundary is assumed to be Maxwellian for ions and a kappa
  ! (kappa=3) distribution for electrons. Boundary temperature and
  ! density are from MHD (or empirical models if Weimer/Tsyganenko fields are
  ! used)
  !
  ! Input: nspec,np,nt,nm,nk,iba,irm,amu,xjac,Den_IC,Temp_IC,vel
  ! Output: fb
  ! 
  Use ModGmCimi, ONLY:  Temp_IC, Temppar_IC, DoAnisoPressureGMCoupling
  use ModCimi,        ONLY: MinLonPar,MaxLonPar, f2
  use ModCimiGrid, ONLY: MinLonPar,MaxLonPar
  use ModCimiTrace,  ONLY: sinA,ekev,iw2,pp
  use ModCimiBoundary,ONLY: BoundaryDens_IC,BoundaryTemp_IC,BoundaryTempPar_IC

  implicit none

  integer nspec,neng,np,nt,nm,nk,iba(nt),irm(nt),j,n,k,m,ib1
  real amu_I(nspec),energy(nspec,neng),xjac(nspec,np,nm)
  real vel(nspec,np,nt,nm,nk),fb(nspec,nt,nm,nk),pi,xmass,chmass,fb1,fb_temp,vtchm
  real velperp2, velpar2,fbb,fbb1,parE1,perE1,den1,temp32,y2,x2,Psq,Pper2,Ppar2
  real erpp,Vsq,Vper2,Vpar2
  real BoundaryTempperp_IC(nspec,nt)
  real kappa,kappa_plus_one,kappa_minus_half,ln_gamma_diff,gamma_ratio
  real ln_gamma
  real e_min

  kappa=3.
 ! zk1=zkappa+1.
  kappa_plus_one=kappa+1.
  kappa_minus_half=kappa-0.5
  ln_gamma_diff=ln_gamma(kappa_plus_one)-ln_gamma(kappa_minus_half)
  gamma_ratio=exp(ln_gamma_diff)

  pi=acos(-1.)

   !   defining Pperp and Ppar for both cases, anisotropic and isotropic.
   !   Advantages: use the same form of Maxwellian/kappa dist
   !   since in the case of isotropic PSD 
  
  if(DoAnisoPressureGMCoupling) then    !  
          BoundaryTempperp_IC(:,:) = &
             (3*BoundaryTemp_IC(:,:) - BoundaryTempPar_IC(:,:))/2
  else
         BoundaryTempperp_IC(:,:) =  BoundaryTemp_IC(:,:)
         BoundaryTempPar_IC(:,:) = BoundaryTemp_IC(:,:)
  endif

  fb(1:nspec,MinLonPar:MaxLonPar,1:nm,1:nk)=0.
  do n=1,nspec
     e_min=0.
     if (n==nspec) e_min=0.5*energy(n,1) 
     xmass=amu_I(n)*1.673e-27
     chmass=1.6e-19/xmass 
     do j=MinLonPar,MaxLonPar
         parE1=BoundaryTempPar_IC(n,j)/1000.      ! characteristic in keV
         perE1=BoundaryTempPerp_IC(n,j)/1000.     ! characteristic in keV
         den1 =BoundaryDens_IC(n,j)                     ! density in m-3
        temp32=sqrt(parE1)*perE1*1.6e-16**1.5
        ib1=iba(j)+1   ! difference bewteen cimi standalone where ib1=iba(j)
        if (ib1.gt.irm(j)) ib1=irm(j)
        if (n .eq. nspec) then    !    kappa, for electrons only
            fbb1=den1*gamma_ratio/temp32/(2.*pi*kappa*xmass)**1.5
        else
            fbb1=den1/temp32/(2.*pi*xmass)**1.5 ! Maxwellian
        endif       

  !      do k=1,nm   ! magnetic moment
          do m=1,nk   ! K invariant 
              y2=sinA(ib1,j,m)*sinA(ib1,j,m)
              x2=1.-y2
            do k=1,iw2(n,m)   ! magnetic moment            
              fbb=0.
              if ((ekev(n,ib1,j,k,m) .gt. e_min) .and. &
                   (BoundaryDens_IC(n,j).gt.0.)) then

               if (n .eq. nspec) then   ! KAPPA FOR ELECTRONS
                  Psq=pp(n,ib1,j,k,m)*pp(n,ib1,j,k,m)
                  Ppar2=Psq*x2
                  Pper2=Psq*y2
                  erpp=(Ppar2/parE1+Pper2/perE1)/2./xmass/1.6e-16
                  fbb=fbb1/(1.+erpp/kappa)**kappa_plus_one                    

                                 else ! MAXWELLIAN OTHERWISE:
                  Vsq=vel(n,ib1,j,k,m)*vel(n,ib1,j,k,m)
                  Vpar2=Vsq*x2
                  Vper2=Vsq*y2
                  erpp=(Vpar2/parE1+Vper2/perE1)/2./1000./chmass
                  fbb=0.
                  if (erpp.lt.500.) fbb=fbb1/exp(erpp)                      
               endif 
             endif      ! end of ekev(n,ib,j,k,m).gt.e_min.and.den1.gt.0


              fb(n,j,k,m)=fbb*xjac(n,ib1,k)
 !        if ((n.eq.nspec).and.(j.eq.1).and.(m.eq.10)) write(*,*) fbb,xjac(n,ib1,k), &
 !        ekev(n,ib1,j,k,m), vel(n,ib1,j,k,m),xmass,chmass, &
 !        Psq,Ppar2,Vsq,Vpar2, 'f,jac,e,vel,m,chm,Psq,Ppar,Vsq,Vpar'   !   NB Jan 20 2017
           enddo                ! end of m loop
        enddo                   ! end of k loop
     enddo                      ! end of j loop
  enddo                         ! end of n loop



end subroutine boundaryIM
   
!-----------------------------------------------------------------------------
function ln_gamma(xx)
!-----------------------------------------------------------------------------
!    ln(gamma(xx))
!    Added from Mei-Ching's stanadalone CIMI to calculate the natural
!    logarithm of the gamma function, which is needed for calculating
!    kappa distributions for the electrons.  -Colin, 07/25/2015.

  real cof(6)
  real ln_gamma
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
  enddo
  ln_gamma=tmp+log(stp*ser)
     
end function ln_gamma
 
!-------------------------------------------------------------------------------
subroutine ceparaIM(nspec,np,nt,nm,nk,irm,dt,vel,ekev,Have,achar)
  !-----------------------------------------------------------------------------
  ! Routine calculates the depreciation factor of H+, achar, due to charge
  ! exchange loss
  !
  ! Input: irm,nspec,np,nt,nm,nk,dt,vel,ekev,Have     ! Have: bounce-ave [H]
  ! Output: achar
  use ModCimiPlanet,  ONLY: a0_I,a1_I,a2_I,a3_I,a4_I
  use ModCimi,       ONLY: MinLonPar,MaxLonPar
  
  implicit none

  integer np,nt,nspec,nk,irm(nt),nm,i,j,k,m,n
  real vel(nspec,np,nt,nm,nk),ekev(nspec,np,nt,nm,nk),Have(np,nt,nk)
  real achar(nspec,np,nt,nm,nk),dt,Havedt,x,d,sigma,alpha

  do n=1,nspec-1
     do j=MinLonPar,MaxLonPar
        do i=1,irm(j)
           do m=1,nk
              Havedt=Have(i,j,m)*dt
              do k=1,nm
                 x=log10(ekev(n,i,j,k,m))
                 if (x.lt.-2.) x=-2.
                 d=a0_I(n)+a1_I(n)*x+a2_I(n)*x**2+a3_I(n)*x**3+a4_I(n)*x**4
                 sigma=10.**d        ! charge exchange cross section of H+ in m2
                 alpha=vel(n,i,j,k,m)*sigma*Havedt
                 achar(n,i,j,k,m)=exp(-alpha) ! charge. exchange decay rate
              enddo
           enddo
        enddo
     enddo
  enddo

end subroutine ceparaIM

!-------------------------------------------------------------------------------
subroutine set_cimi_potential(CurrentTime,rc)
  !-----------------------------------------------------------------------------
  ! Routine sets the ionospheric potentials from weimer when not using MHD input
  !
  use ModCimiGrid, ONLY: xlatr, xmlt, np, nt
  use ModIeCimi,   ONLY: pot
!  use ModTsyInput, ONLY: xnswa,vswa,bxw,byw,bzw
  use ModMpi 
  use EIE_ModWeimer, ONLY: setmodel00, boundarylat00, epotval00
  use ModNumConst, ONLY: cPi, cRadToDeg
  use ModIndicesInterfaces
  use ModCimiTrace, ONLY: UsePotential
  implicit none

  real, intent(in) :: CurrentTime,rc

  logical :: UseAL
  real :: ALindex
  real :: xnsw,vsw,bx,by,bz !solar wind values at current time
  real :: Tilt, angle, Bt, gLAT, gMLT, BnLat
  
  integer :: i, iError, j
  
  ! geopack common block stuff
  COMMON /GEOPACK/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,&
       CPS,SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,A33,DS3,&
       K,IY,CGST,SGST,BA(6)
  real::ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,&
       CPS,SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,A33,DS3,&
       CGST,SGST,BA
  integer::K,IY

  !----------------------------------------------------------------------------
  
  UseAL=.false.
  ALindex=10.          ! arbitrary value 

  !  Setup for Weimer's electric field model if iconvect = 1
  Tilt=psi*cRadToDeg          ! dipole tilt angle in degree
!  if (UseGm)then
!     ! When using GM, it is assumed that the solar wind will be passed 
!     xnsw=xnswa(1)
!     vsw=vswa(1)
!     bx=bxw(1)
!     by=byw(1)
!     bz=bzw(1)
!  else
     ! When GM is not used set the solar wind parameters form the imf.dat file
     ! defined in the PARAM.in 

     call get_IMF_Bz(CurrentTime, bz, iError)
          
     call get_IMF_By(CurrentTime, by, iError)
          
     call get_SW_V(CurrentTime, vsw, iError)
     
     call get_SW_N(CurrentTime, xnsw, iError)
     
     !only use the absolute value of sw velocity
     vsw=abs(vsw)
     
     if (iError /= 0) then
        call con_stop&
             ("IM_ERROR: Problem setting solar wind in set_cimi_potential")
     end if
!  endif
  angle=atan2(by,bz)*180./cPi       ! degrees from northward toward +Y
  Bt=sqrt(by*by+bz*bz)             ! Magnitude of IMF in Y-Z plane in nT
  call SetModel00(angle,Bt,Tilt,vsw,xnsw,ALindex,UseAL)
  

  !  Find potential (in Volt) at the ionosphere
  do i=1,np
     gLAT=acos(cos(xlatr(i))/sqrt(rc))*cRadToDeg !invariant lat. in degree
     do j=1,nt
        gMLT=xmlt(j)                 ! convert mlt from radians to hour 0-24.
        BnLat=BoundaryLat00(gMLT)    ! boundary latitude 
        if (UsePotential) then
           pot(i,j)=0.0
           if (gLAT.gt.BnLat) pot(i,j)=EpotVal00(gLAT,gMLT)*1000.  ! Volt
        else
           pot(i,j)=0.0
        endif
     enddo
  enddo
  
end subroutine set_cimi_potential
  


!-------------------------------------------------------------------------------
subroutine driftV(nspec,np,nt,nm,nk,irm,re_m,Hiono,dipmom,dphi,xlat, &
     dlat,ekev,pot,vl,vp)
  !-----------------------------------------------------------------------------
  ! Routine calculates the drift velocities
  !
  ! Input: re_m,Hiono,dipmom,dphi,xlat,dlat,ekev,pot,nspec,np,nt,nm,nk,irm
  ! Output: vl,vp
  use ModCimiGrid, ONLY: iProc,nProc,iComm,MinLonPar,MaxLonPar, &
       iProcLeft, iLonLeft, iProcRight, iLonRight
  use ModMpi
  use ModCimiTrace, ONLY: UseCorotation
  implicit none

  integer nspec,np,nt,nm,nk,irm(nt),n,i,ii,j,k,m,i0,i2,j0,j2,icharge
  real kfactor,xlat(np),xlatr(np),dlat(np),ekev(nspec,np,nt,nm,nk),pot(np,nt)
  real ksai,ksai1,xlat1,sf0,sf2,dlat2,re_m,Hiono,dipmom,dphi,pi,dphi2,cor
  real ham(np,nt),vl(nspec,0:np,nt,nm,nk),vp(nspec,np,nt,nm,nk)

  ! MPI status variable
  integer :: iStatus_I(MPI_STATUS_SIZE), iError

  pi=acos(-1.)
  dphi2=dphi*2.
  kfactor=dipmom/(re_m+Hiono*1000.)
  if (UseCorotation) then
     cor=2.*pi/86400.                        ! corotation speed in rad/s
  else
     cor=0.
  endif
  xlatr=xlat*pi/180.

  nloop: do n=1,nspec
     if (n < nspec) then
        icharge=1
     else
        icharge=-1
     endif

     mloop: do m=1,nk
        kloop: do k=1,nm  

           ! ham: Hamiltonian/q
           ham(1:np,1:nt)=icharge*ekev(n,1:np,1:nt,k,m)*1000.+pot(1:np,1:nt)

           ! When nProc>1 exchange ghost cell info for ham and irm
           if (nProc >1) then
              !send to neigboring Procs
              call MPI_send(ham(1:np,MaxLonPar),np,MPI_REAL,iProcRight,&
                   1,iComm,iError)
              call MPI_send(ham(1:np,MinLonPar),np,MPI_REAL,iProcLeft,&
                   2,iComm,iError)
              call MPI_send(irm(MaxLonPar),1,MPI_INTEGER,iProcRight,&
                   3,iComm,iError)
              call MPI_send(irm(MinLonPar),1,MPI_INTEGER,iProcLeft,&
                   4,iComm,iError)
              !recieve from neigboring Procs
              call MPI_recv(ham(1:np,iLonLeft),np,MPI_REAL,iProcLeft,&
                   1,iComm,iStatus_I,iError)
              call MPI_recv(ham(1:np,iLonRight),np,MPI_REAL,iProcRight,&
                   2,iComm,iStatus_I,iError)
              call MPI_recv(irm(iLonLeft),1,MPI_INTEGER,iProcLeft,&
                   3,iComm,iStatus_I,iError)
              call MPI_recv(irm(iLonRight),1,MPI_INTEGER,iProcRight,&
                   4,iComm,iStatus_I,iError)
              
           endif

           ! calculate drift velocities vl and vp
           iloop: do i=0,np
              ii=i
              if (i.eq.0) ii=1
              if (i.ge.1) ksai=kfactor*sin(2.*xlatr(i))
              if (i.lt.np) xlat1=0.5*(xlatr(ii)+xlatr(i+1))    ! xlat(i+0.5)
              ksai1=kfactor*sin(2.*xlat1)                   ! ksai at i+0.5
              jloop: do j=MinLonPar,MaxLonPar
                 j0=j-1
                 if (j0.lt.1) j0=j0+nt
                 j2=j+1
                 if (j2.gt.nt) j2=j2-nt

                 ! calculate vl 
                 if (irm(j0).gt.i.and.irm(j2).gt.i) then
                    sf0=0.5*ham(ii,j0)+0.5*ham(i+1,j0)
                    sf2=0.5*ham(ii,j2)+0.5*ham(i+1,j2)
                    vl(n,i,j,k,m)=-(sf2-sf0)/dphi2/ksai1   ! vl at (i+0.5,j)
                 else
                    vl(n,i,j,k,m)=vl(n,i-1,j,k,m)
                 endif

                 ! calculate vp
                 if (i.ge.1) then
                    if (irm(j2).gt.i) then
                       i0=i-1
                       if (i.eq.1) i0=1
                       i2=i+1
                       if (i.eq.np) i2=np
                       dlat2=xlatr(i2)-xlatr(i0)
                       sf0=0.5*(ham(i0,j2)+ham(i0,j))
                       sf2=0.5*(ham(i2,j2)+ham(i2,j))
                       vp(n,i,j,k,m)=cor+(sf2-sf0)/dlat2/ksai  ! vp at (i,j+0.5)
                    else
                       vp(n,i,j,k,m)=vp(n,i-1,j,k,m)
                    endif
                 endif

              enddo jloop
           enddo iloop
        enddo kloop
     enddo mloop
  enddo nloop

end subroutine driftV


!-------------------------------------------------------------------------------
subroutine driftIM(iw2,nspec,np,nt,nm,nk,dt,dlat,dphi,brad,rb,vl,vp, &
     fb,f2,driftin,driftout,ib0)
  !-----------------------------------------------------------------------------
  ! Routine updates f2 due to drift
  !
  ! Input: iw2,nspec,np,nt,nm,nk,iba,dt,dlat,dphi,brad,rb,vl,vp,fbi
  ! Input/Output: f2,ib0,driftin,driftout
  use ModCimiGrid, ONLY: MinLonPar, MaxLonPar
  use ModCimiTrace, ONLY: iba, ekev
  use ModCimiGrid, ONLY: iProc,nProc,iComm,MinLonPar,MaxLonPar, &
       iProcLeft, iLonLeft, iProcRight, iLonRight, d4Element_C    
  use ModMpi
  implicit none

  integer nk,nspec,np,nt,nm,ib0(nt)
  integer n,i,j,k,m,j1,j_1,ibaj,ib,ibo,nrun,nn
  integer iw2(nspec,nk)
  real dt,dlat(np),dphi,brad(np,nt),&
       vl(nspec,0:np,nt,nm,nk),vp(nspec,np,nt,nm,nk)
  real rb,fb(nspec,nt,nm,nk),f2(nspec,np,nt,nm,nk)
  real f2d(np,nt),cmax,cl1,cp1,cmx,dt1,fb0(nt),fb1(nt),fo_log,fb_log,f_log
  real slope,cl(np,nt),cp(np,nt),fal(0:np,nt),fap(np,nt),&
       fupl(0:np,nt),fupp(np,nt)
  real driftin(nspec),driftout(nspec),dEner,dPart,&
       dEnerLocal_C(nt),dPartLocal_C(nt)
  logical :: UseUpwind=.false.

  ! MPI status variable
  integer :: iStatus_I(MPI_STATUS_SIZE)
  integer :: iError
  real, allocatable :: cmax_P(:)
  
  if (.not.allocated(cmax_P) .and. nProc>1) allocate(cmax_P(nProc))
  
  ! When nProc>1 pass iba from neighboring procs
  if (nProc >1) then
     !send to neigboring Procs
     call MPI_send(iba(MinLonPar),1,MPI_INTEGER,iProcLeft,&
          1,iComm,iError)
     !recieve from neigboring Procs
     call MPI_recv(iba(iLonRight),1,MPI_INTEGER,iProcRight,&
          1,iComm,iStatus_I,iError)
  endif

  nloop: do n=1,nspec
     mloop: do m=1,nk
        kloop: do k=1,iw2(n,m)
           f2d(1:np,1:nt)=f2(n,1:np,1:nt,k,m)         ! initial f2
           ! find nrun and new dt (dt1)
           cmax=0.
           do j=MinLonPar,MaxLonPar
              j1=j+1
              if (j1.gt.nt) j1=j1-nt
              ibaj=max(iba(j),iba(j1))
              do i=1,ibaj
                 cl1=dt/dlat(i)*vl(n,i,j,k,m)
                 cp1=dt/dphi*vp(n,i,j,k,m)
                 cmx=max(abs(cl1),abs(cp1))
                 cmax=max(cmx,cmax)
              enddo
           enddo

           !get same cmax on all procs
           if (nProc > 1) then
              call MPI_allgather(cmax,1,MPI_REAL,cmax_P,1,MPI_REAL,iComm,iError)
              cmax=maxval(cmax_P)
           endif

           nrun=ifix(cmax/0.50)+1     ! nrun to limit the Courant number
           dt1=dt/nrun                ! new dt
           ! Setup boundary fluxes and Courant numbers
           do j=MinLonPar,MaxLonPar
              ib=iba(j)
              ibo=ib0(j)
              fb0(j)=f2d(1,j)                   ! psd at inner boundary
              fb1(j)=fb(n,j,k,m)              ! psd at outer boundary
              if (ib.gt.ibo) then             ! during dipolarization
                 fo_log=-50.
                 if (f2d(ibo,j).gt.1.e-50) fo_log=log10(f2d(ibo,j))
                 fb_log=-50.
                 if (fb1(j).gt.1.e-50) fb_log=log10(fb1(j))
                 slope=(fo_log-fb_log)/(brad(ibo,j)-rb)
                 do i=ibo+1,ib
                    f_log=fo_log+slope*(brad(i,j)-brad(ibo,j))
                    f2d(i,j)=10.**f_log
                 enddo
              endif
              do i=1,np
                 cl(i,j)=dt1/dlat(i)*vl(n,i,j,k,m)
                 cp(i,j)=dt1/dphi*vp(n,i,j,k,m)
              enddo
           enddo

           ! run drift nrun times
           do nn=1,nrun
              UseUpwind=.false.
              ! When nProc>1, pass fb0, fb1, and f2d
              !send to neigboring Procs
             if (nProc>1) then
                 !send f2d ghostcells
                 call MPI_send(f2d(1:np,MaxLonPar),np,MPI_REAL,iProcRight,&
                      3,iComm,iError)
                 call MPI_send(f2d(1:np,MinLonPar:MinLonPar+1),2*np,MPI_REAL,&
                      iProcLeft,4,iComm,iError)
                 !recieve f2d ghostcells from neigboring Procs
                 call MPI_recv(f2d(1:np,iLonLeft),np,MPI_REAL,iProcLeft,&
                      3,iComm,iStatus_I,iError)
                 call MPI_recv(f2d(1:np,iLonRight:iLonRight+1),2*np,MPI_REAL,&
                      iProcRight,4,iComm,iStatus_I,iError)

                 !send fb0 ghostcells
                 call MPI_send(fb0(MinLonPar:MinLonPar+1),2,MPI_REAL,&
                      iProcLeft,5,iComm,iError)
                 call MPI_send(fb0(MaxLonPar),1,MPI_REAL,iProcRight,&
                      6,iComm,iError)
                 !recieve fb0 from neigboring Procs
                 call MPI_recv(fb0(iLonRight:iLonRight+1),2,MPI_REAL,&
                      iProcRight,5,iComm,iStatus_I,iError)
                 call MPI_recv(fb0(iLonLeft),1,MPI_REAL,iProcLeft,&
                      6,iComm,iStatus_I,iError)

                 !send fb1 ghostcells
                 call MPI_send(fb1(MinLonPar:MinLonPar+1),2,MPI_REAL,&
                      iProcLeft,7,iComm,iError)
                 call MPI_send(fb1(MaxLonPar),1,MPI_REAL,iProcRight,&
                      8,iComm,iError)
                 !recieve fb1 from neigboring Procs
                 call MPI_recv(fb1(iLonRight:iLonRight+1),2,MPI_REAL,&
                      iProcRight,7,iComm,iStatus_I,iError)
                 call MPI_recv(fb1(iLonLeft),1,MPI_REAL,iProcLeft,&
                      8,iComm,iStatus_I,iError)
              endif
              call FLS_2D(np,nt,iba,fb0,fb1,cl,cp,f2d,fal,fap,fupl,fupp)
              fal(0,1:nt)=f2d(1,1:nt)
              ! When nProc>1 pass needed ghost cell info for fap,fupp and cp
              if (nProc>1) then
                 !send fap ghostcells
                 call MPI_send(fap(:,MaxLonPar),np,MPI_REAL,iProcRight,&
                      9,iComm,iError)
                 !recieve fap from neigboring Procs
                 call MPI_recv(fap(:,iLonLeft),np,MPI_REAL,iProcLeft,&
                      9,iComm,iStatus_I,iError)

                 !send fupp ghostcells
                 call MPI_send(fupp(:,MaxLonPar),np,MPI_REAL,iProcRight,&
                      10,iComm,iError)
                 !recieve fupp from neigboring Procs
                 call MPI_recv(fupp(:,iLonLeft),np,MPI_REAL,iProcLeft,&
                      10,iComm,iStatus_I,iError)

                 !send cp ghostcells
                 call MPI_send(cp(:,MaxLonPar),np,MPI_REAL,iProcRight,&
                      11,iComm,iError)
                 !recieve cp from neigboring Procs
                 call MPI_recv(cp(:,iLonLeft),np,MPI_REAL,iProcLeft,&
                      11,iComm,iStatus_I,iError)
              endif

              jloop: do j=MinLonPar,MaxLonPar
                 j_1=j-1
                 if (j_1.lt.1) j_1=j_1+nt
                 iloop: do i=1,iba(j)
                    f2d(i,j)=f2d(i,j)+dt1/dlat(i)* &
                         (vl(n,i-1,j,k,m)*fal(i-1,j)-vl(n,i,j,k,m)*fal(i,j))+ &
                         cp(i,j_1)*fap(i,j_1)-cp(i,j)*fap(i,j)
                    if (f2d(i,j).lt.0.) then
                       if (f2d(i,j).gt.-1.e-30) then
                          f2d(i,j)=0.
                       else
                          write(*,*)'IM WARNING: f2d < 0 in drift ',n,i,j,k,m
                          write(*,*)'IM WARNING: '//&
                               'Retrying step with upwind scheme'
                          UseUpwind=.true.
                          exit jloop
                       endif
                    endif
                 enddo iloop
                 ! Calculate gain or loss at the outer boundary
                 dPartLocal_C(j) = -dt1/dlat(iba(j)) &
                      *vl(n,iba(j),j,k,m)*fal(iba(j),j)*&
                      d4Element_C(n,iba(j),k,m)
                 dEnerLocal_C(j)=ekev(n,iba(j),j,k,m)*dPartLocal_C(j)
              enddo jloop

              ! When regular scheme fails, try again with upwind scheme before 
              ! returning an error
              if (UseUpwind) then
                 fupl(0,1:nt)=f2d(1,1:nt)
                 do j=MinLonPar,MaxLonPar
                    j_1=j-1
                    if (j_1.lt.1) j_1=j_1+nt
                    iLoopUpwind: do i=1,iba(j)
                       f2d(i,j)=f2d(i,j)+dt1/dlat(i)* &
                        (vl(n,i-1,j,k,m)*fupl(i-1,j)-vl(n,i,j,k,m)*fupl(i,j))+ &
                        cp(i,j_1)*fupp(i,j_1)-cp(i,j)*fupp(i,j)
                       if (f2d(i,j).lt.0.) then
                          if (f2d(i,j).gt.-1.e-30) then
                             f2d(i,j)=0.
                          else
                             !write(*,*)'IM ERROR: f2d < 0 in drift ',n,i,j,k,m
                             !call CON_STOP('CIMI dies in driftIM')
                             !write(*,*)'IM WARNING: f2d < 0 in drift ',n,i,j,k,m
                             !write(*,*)'IM WARNING: upwind scheme failed, setting f2d(i,j)=0.0'
                             !write(*,*)'IM WARNING: repeated failure may need to be examined'
                             ! should now have f2d(i,j)=0. but didnt before

                             write(*,*)'IM WARNING: f2d < 0 in drift ',n,i,j,k,m
                             write(*,*)'IM WARNING: '//&
                                  'upwind scheme failed, making iba(j)=i'
                             write(*,*)'IM WARNING: '//&
                                  'repeated failure may need to be examined'
                             f2d(i,j)=0.0
                             iba(j)=i
                             exit iLoopUpwind
                          endif
                       endif
                    enddo iLoopUpwind
                    ! Calculate gain or loss at the outer boundary
                    dPartLocal_C(j)=-dt1/dlat(iba(j))*&
                         vl(n,iba(j),j,k,m)*fupl(iba(j),j)*&
                         d4Element_C(n,iba(j),k,m)
                    dEnerLocal_C(j)=ekev(n,iba(j),j,k,m)*dPartLocal_C(j)
                 enddo
              endif
              !sum all dEner to root proc
              if(nProc>1) then
                 call MPI_REDUCE (sum(dPartLocal_C(MinLonPar:MaxLonPar)), &
                      dPart, 1, MPI_REAL, MPI_SUM, 0, iComm, iError)
                 call MPI_REDUCE (sum(dEnerLocal_C(MinLonPar:MaxLonPar)), &
                      dEner, 1, MPI_REAL, MPI_SUM, 0, iComm, iError)
              else
                 dPart=sum(dPartLocal_C)
                 dEner=sum(dEnerLocal_C)
              endif
              
              if(iProc==0) then
                 if (dPart.gt.0.) driftin(n)=driftin(n)+dEner
                 if (dPart.lt.0.) driftout(n)=driftout(n)+dEner
              else
                 driftin(n)=0.
                 driftout(n)=0.
              endif
           enddo          ! end of do nn=1,nrun
           f2(n,1:np,1:nt,k,m)=f2d(1:np,1:nt)
        enddo kloop
     enddo mloop
  enddo nloop


  ! Update ib0
  ib0(1:nt)=iba(1:nt)
end subroutine driftIM


!-------------------------------------------------------------------------------
subroutine charexchangeIM(np,nt,nm,nk,nspec,iba,achar,f2)
  !-----------------------------------------------------------------------------
  ! Routine updates f2 due to charge exchange loss
  !
  ! Input: np,nt,nm,nk,nspec,achar   ! charge exchange depreciation of H+ 
  ! Input/Output: f2
  use ModCimi,       ONLY: MinLonPar,MaxLonPar
  implicit none

  integer np,nt,nm,nk,nspec,iba(nt),n,i,j
  real achar(nspec,np,nt,nm,nk),f2(nspec,np,nt,nm,nk)

  do n=1,nspec-1             
     do j=MinLonPar,MaxLonPar
        do i=1,iba(j)
           f2(n,i,j,1:nm,1:nk)=f2(n,i,j,1:nm,1:nk)*achar(n,i,j,1:nm,1:nk)
        enddo
     enddo
  enddo
 
end subroutine charexchangeIM

!******************************************************************************
!                                StDiTime                                      
!  Routine calculate the strong diffusion lifetime for electrons.     
!*****************************************************************************
subroutine StDiTime(dt,vel,volume,rc,re_m,xme,iba)
  use ModCimi,       ONLY: SDtime
  use ModCimiGrid,   ONLY: np,nt,nm,nk, xlatr,MinLonPar,MaxLonPar
  use ModCimiPlanet, ONLY: nspec
  real vel(nspec,np,nt,nm,nk),volume(np,nt)
  integer iba(nt)
  
  
  eb=0.25                         ! fraction of back scatter e-   
  xmer3=xme/(rc*re_m)**3
  
  do j=MinLonPar,MaxLonPar
     do i=1,iba(j)
        !              xlat2=xlati(i)*xlati(i)!- from M.-Ch., Aug 1 2007  
        !              Bi=xmer3*sqrt(3.*xlat2+1.)     
        sinlat2=sin(xlatr(i))*sin(xlatr(i))
        Bi=xmer3*sqrt(3.*sinlat2+1.)      ! magnetic field at ionosphere 
        
        vBe=2.*volume(i,j)*Bi/(1.-eb)
        do k=1,nm
           do m=1,nk
              SDtime1=vBe/vel(nspec,i,j,k,m) !strong diff T,(gamma*mo/p = 1/v)
              SDtime(i,j,k,m)=exp(-dt/SDtime1)
           enddo
        enddo
     enddo
  enddo
  
  return
end subroutine StDiTime


!***********************************************************************   
!                            StrongDiff                                     
!  Routine calculate the change of electron psd (f2) by strong diffusion  
!***********************************************************************        
subroutine StrongDiff(iba)                               
  use ModCimi,       ONLY: SDtime,f2
  use ModCimiGrid,   ONLY: np,nt,nm,nk,MinLonPar,MaxLonPar
  use ModCimiPlanet, ONLY: nspec  
  implicit none
  integer iba(nt),i,j,k,m
  
  do j=MinLonPar,MaxLonPar
     do i=2,iba(j)
        do m=1,nk
           do k=1,nm
              f2(nspec,i,j,k,m)=f2(nspec,i,j,k,m)*SDtime(i,j,k,m)
           enddo
        enddo
     enddo
  enddo
  
  return
end subroutine StrongDiff



!-------------------------------------------------------------------------------
subroutine lossconeIM(np,nt,nm,nk,nspec,iba,alscone,f2)
  !-----------------------------------------------------------------------------
  ! Routine calculate the change of f2 due to lossconeIM loss
  ! 
  ! Input: np,nt,nm,nk,nspec,iba,alscone
  ! Input/Output: f2
  use ModCimi,       ONLY: MinLonPar,MaxLonPar
  implicit none

  integer np,nt,nm,nk,nspec,iba(nt),n,i,j,k,m
  real alscone(nspec,np,nt,nm,nk),f2(nspec,np,nt,nm,nk)

  do n=1,nspec
     do j=MinLonPar,MaxLonPar
        do i=1,iba(j)
           do k=1,nm
              do m=1,nk
                 if (alscone(n,i,j,k,m).lt.1.) &
                      f2(n,i,j,k,m)=f2(n,i,j,k,m)*alscone(n,i,j,k,m)
              enddo
           enddo
        enddo
     enddo
  enddo

end subroutine lossconeIM

!==============================================================================

subroutine sume(xle)
!-------------------------------------------------------------------------------
! Routine updates rbsum and xle
! 
! Input: f2,ekev,iba
! Input/Output: rbsum,xle
  use ModCimi,       ONLY: rbsumLocal
  use ModCimiTrace, ONLY: iba
  use ModCimiGrid,   ONLY: nProc,iProc,iComm
  use ModCimiPlanet, ONLY: nspec
  use ModMPI
  implicit none
  
  real, intent(inout):: xle(nspec)
  
  integer n,i,j,k,m,iError
  real rbsumLocal0,xleChange,xleChangeLocal

  do n=1,nspec
     rbsumLocal0=rbsumLocal(n)
     
     call calc_rbsumlocal(n)

     xleChangeLocal=rbsumLocal(n)-rbsumLocal0
     
     if (nProc >1) call MPI_REDUCE (xleChangeLocal, xleChange, 1, MPI_REAL, &
           MPI_SUM, 0, iComm, iError)

     if(iProc==0) then 
        xle(n)=xle(n)+xleChange
     else
        xle(n)=0.0
     endif
 enddo

end subroutine sume

!==============================================================================

subroutine calc_rbsumlocal(iSpecies)
  use ModCimi,       ONLY: f2,rbsumLocal
  use ModCimiGrid,   ONLY: np,nm,nk,MinLonPar,MaxLonPar,d4Element_C
  use ModCimiTrace, ONLY: iba, ekev
  implicit none

  integer, intent(in) :: iSpecies
  real    :: weight
  integer :: i,j,k,m
  !-----------------------------------------------------------------------------
  rbsumLocal(iSpecies)=0.
  do j=MinLonPar,MaxLonPar
     do i=1,iba(j)
        do k=1,nm
           do m=1,nk
              weight=f2(iSpecies,i,j,k,m)*d4Element_C(iSpecies,i,k,m)*ekev(iSpecies,i,j,k,m)
              rbsumLocal(iSpecies)=rbsumLocal(iSpecies)+weight        ! rbsum in keV
           enddo
        enddo
     enddo
  enddo

end subroutine calc_rbsumlocal

!==============================================================================

subroutine sume_cimi(OperatorName)
! Routine updates rbsum and xle
! 
! Input: f2,ekev,iba
! Input/Output: rbsum,xle

! in CIMI: eChangeOperator=xle(ns,ir,ip,nOperator,je+2) Here we create
!   	one array for all operators (plus one dimension)
!  
! eTimeAccumultv=esum(ns,ir,ip,je+2) No need for additional dimension
! 	because it's total energy
!  
!!!!!!!!!!!!!!!!
  use ModCimi,       ONLY: &
       f2,rbsum=>rbsumLocal,rbsumGlobal,rcsum=>rcsumLocal,rcsumGlobal,&
       xle=>eChangeOperator_VICI,ple=>pChangeOperator_VICI, &
       eChangeGlobal,eChangeLocal, &
       esum=>eTimeAccumult_ICI,psum=>pTimeAccumult_ICI
  use ModCimiTrace, ONLY: iba,irm,ekev,ro,iw2
  use ModCimiGrid,   ONLY: &
       nProc,iProc,iComm, MinLonPar,MaxLonPar,d4Element_C, &
       ip=>np,ir=>nt,im=>nm,ik=>nk,je=>neng,Energy,Ebound
  use ModCimiPlanet, ONLY: nspec
  use ModMPI
  use ModWaveDiff,    ONLY: testDiff_aE 

  implicit none

  integer OperatorName   ! defined in ModCimi
  real    :: weight,ekev1,weighte,dee,dpe
  integer n,i,j,k,m,iError,kk
  real gride1(0:je+1),e0(ip,ir,je+2),p0(ip,ir,je+2)

! Set up gride1(0) and gride1(je+1)
      gride1(0)=0.
      gride1(je+1)=1.e10     ! arbitrary large number

! Calculate esum, psum, etc.
    do n=1,nspec
       eChangeLocal(n, OperatorName) = 0.
       e0(1:ip,MinLonPar:MaxLonPar,1:je+2)=&
            esum(n,1:ip,MinLonPar:MaxLonPar,1:je+2)
       p0(1:ip,MinLonPar:MaxLonPar,1:je+2)=&
            psum(n,1:ip,MinLonPar:MaxLonPar,1:je+2)
       esum(n,1:ip,MinLonPar:MaxLonPar,1:je+2)=0.
       psum(n,1:ip,MinLonPar:MaxLonPar,1:je+2)=0.
       gride1(1:je)=Ebound(n,1:je)

       do j=MinLonPar,MaxLonPar
          do i=1,irm(j)
                do m=1,ik
                   do k=1,iw2(n,m)
                      ekev1=ekev(n,i,j,k,m)
                      weight=d4Element_C(n,i,k,m)*f2(n,i,j,k,m)
                      weighte=ekev1*weight
                      psum(n,i,j,je+2)=psum(n,i,j,je+2)+weight
                      esum(n,i,j,je+2)=esum(n,i,j,je+2)+weighte
                      rbsum(n)=rbsum(n)+weighte
                      if (ro(i,j).le.6.6) rcsum(n)=rcsum(n)+weighte
                      kkloop: do kk=1,je+1
                         if ( ekev1 .gt. gride1(kk-1) .and. &
                              ekev1 .le. gride1(kk) ) then
                            psum(n,i,j,kk)=psum(n,i,j,kk)+weight
                            esum(n,i,j,kk)=esum(n,i,j,kk)+weighte
                            exit kkloop
                         endif
                      enddo kkloop
                   enddo
                enddo

             do kk=1,je+2
                dee = esum(n,i,j,kk) - e0(i,j,kk)
                xle(n,i,j,kk,OperatorName) = &
                     xle(n,i,j,kk,OperatorName) + dee
                dpe=psum(n,i,j,kk)-p0(i,j,kk)
                ple(n,i,j,kk,OperatorName) = &
                     ple(n,i,j,kk,OperatorName) + dpe
             enddo

             eChangeLocal(n, OperatorName) = &
                  eChangeLocal(n, OperatorName) + &
                  xle(n,i,j,je+2,OperatorName)
             
          enddo                ! end of do i=1,iba(j)
       enddo                   ! end of do j=MinLonPar,MaxLonPar

       if (nProc > 1) then
          call MPI_REDUCE(&
               rbsum(n), rbsumGlobal(n), 1, &
               MPI_REAL, MPI_SUM, 0, iComm, iError)
          call MPI_REDUCE(&
               rcsum(n), rcsumGlobal(n), 1, &
               MPI_REAL, MPI_SUM, 0, iComm, iError)
          call MPI_REDUCE(&
               eChangeLocal(n, OperatorName),&
               eChangeGlobal(n, OperatorName), 1, &
               MPI_REAL, MPI_SUM, 0, iComm, iError)
       else 
          rbsumGlobal(n) = rbsum(n)
          rcsumGlobal(n) = rcsum(n)
          eChangeGlobal(n,OperatorName) = &
               eChangeLocal(n,OperatorName)
       endif

    enddo                      ! end of do n=1,nSpecies

    if (testDiff_aE) &
         write(*,*) 'tot particles, el: ',psum(nspec,30,5,je+2)

end subroutine sume_cimi


!==============================================================================

subroutine cimi_output(np,nt,nm,nk,nspec,neng,npit,iba,ftv,f2,ekev, &
     sinA,energy,sinAo,delE,dmu,amu_I,xjac,pp,xmm, &
     dmm,dk,xlat,dphi,re_m,Hiono,vl,vp, &
     flux,fac,phot,Ppar_IC,Pressure_IC,PressurePar_IC,vlEa,vpEa,psd)
  !-----------------------------------------------------------------------------
  ! Routine calculates CIMI output, flux, fac and phot from f2
  ! Routine also converts the particle drifts from (m,K) space to (E,a) space
  !
  ! Input: np,nt,nm,nk,nspec,neng,npit,iba,ftv,f2,ekev,sinA,energy,sinAo,xjac
  !        delE,dmu,amu_I,xjac,pp,xmm,dmm,dk,xlat,dphi,re_m,Hiono,vl,vp
  ! Output: flux,fac,phot,Ppar_IC,Den_IC,Temp_IC,vlEa,vpEa,psd
  Use ModGmCimi, ONLY: Den_IC, Temp_IC
  use ModConst,   ONLY: cProtonMass
  use ModNumConst,ONLY: cPi, cDegToRad
  use ModCimiGrid,ONLY: iProc,nProc,iComm,MinLonPar,MaxLonPar,&
       iProcLeft, iLonLeft, iProcRight, iLonRight
  use ModMpi
  use ModCimiTrace, ONLY: ro
  implicit none

  integer np,nt,nm,nk,nspec,neng,npit,iba(nt),i,j,k,m,n,j1,j_1
  real f2(nspec,np,nt,nm,nk),ekev(nspec,np,nt,nm,nk),sinA(np,nt,nk),re_m,Hiono,rion,vl(nspec,np,nt,nm,nk),vp(nspec,np,nt,nm,nk)
  real ftv(np,nt),ftv1,energy(nspec,neng),sinAo(npit),delE(nspec,neng),dmu(npit),aloge(nspec,neng)
  real flux2D(nm,nk),pp(nspec,np,nt,nm,nk),xjac(nspec,np,nm)
  real vl2D(nm,nk),vp2D(nm,nk)
  real sinA1D(nk),cosA2(nk),flx,ekev2D(nm,nk),flx_lo,pf(nspec),delEE(nspec,neng),pi,cosAo2(npit)
  real vl_lo,vp_lo
  real sina1,sina0,dcosa
  real amu_I(nspec),amu1,psd1,psd(nspec,np,nt,nm,nk),fave(nspec,np,nt,neng)
  real xmm(nspec,0:nm+1),dmm(nspec,nm),dk(nk),xlat(np),xlatr(np),dphi,eta(nspec,np,nt,nm,nk)
  real flux(nspec,np,nt,neng,npit),detadi,detadj,dwkdi,dwkdj
  real vlEa(nspec,np,nt,neng,npit),vpEa(nspec,np,nt,neng,npit)
  real fac(np,nt),phot(nspec,np,nt),Ppar_IC(nspec,np,nt)
  real Pressure_IC(nspec,np,nt), PressurePar_IC(nspec,np,nt)
  real Pressure0, Pressure1, PressurePar1, Coeff
  integer :: iStatus_I(MPI_STATUS_SIZE), iError
  logical, parameter :: DoCalcFac=.true.
  flux=0.
  fac=0.
  eta=0.
  phot=0.
  Ppar_IC = 0.
  PressurePar_IC = 0.

  ! Some constants for pressure, fac calculations
  rion=re_m+Hiono*1000.                      ! ionosphere distance in meter
  do n=1,nspec
     pf(n)=4.*cPi*1.e4/3.*sqrt(2.*cProtonMass*amu_I(n))*sqrt(1.6e-16)*1.e9  ! phot(nPa)
     delEE(n,1:neng) = delE(n,1:neng) * sqrt(energy(n,1:neng))
     aloge(n,1:neng) = log10(energy(n,1:neng)) 
  enddo
!  delEE=delE*sqrt(energy)
  xlatr=xlat*cDegToRad

  ! Calculate CIMI ion density (m^-3), Den_IC, and flux (cm^-2 s^-1 keV^-1 sr^-1)
  ! at fixed energy & pitch-angle grids 

 ! aloge=log10(energy)
  jloop1: do j=MinLonPar,MaxLonPar
     iloop1: do i=1,iba(j)
        ftv1=ftv(i,j)     ! ftv1: flux tube volume in m^3/Wb
        nloop: do n=1,nspec
           Pressure0=0.0
           Pressure1=0.0
           PressurePar1=0.0
           Den_IC(n,i,j)=0.0
           amu1=amu_I(n)**1.5
!!!! Calculate Den_IC, and 2D flux, fl2D(log), ekev2D(log) and sinA1D
           do m=1,nk
              sinA1D(m) = sinA(i,j,m)
              cosA2(m) = 1 - sinA1D(m)**2
           end do
           do m=1,nk
              if (m.eq.1) sina0=1.
              if (m.gt.1) sina0=0.5*(sinA1D(m)+sinA1D(m-1))
              if (m.eq.nk) sina1=0.
              if (m.lt.nk) sina1=0.5*(sinA1D(m)+sinA1D(m+1))
              dcosa=sqrt(1.-sina1*sina1)-sqrt(1.-sina0*sina0)
              do k=1,nm
                 !write(*,*) 'n,i,k,xjac(n,i,k)',n,i,k,xjac(n,i,k)
                 psd1=f2(n,i,j,k,m)/1.e20/1.e19/xjac(n,i,k)  ! mug^-3cm^-6s^3
                 flx=psd1*(1.6e19*pp(n,i,j,k,m))*pp(n,i,j,k,m)
                 flux2D(k,m)=-50.
                 if (flx.gt.1.e-50) flux2D(k,m)=log10(flx)
                 ekev2D(k,m)=log10(ekev(n,i,j,k,m))
                 eta(n,i,j,k,m)=amu1*1.209*psd1*sqrt(xmm(n,k))*dmm(n,k)*dk(m)
                 psd(n,i,j,k,m)=psd1

                 ! Stores drift velocities in (m,k)
                 ! NOTE:: vl is calculated at dlambda/dt [rad/s] in
                 ! ionosphere.  This variable still needs to be
                 ! converted to radial velocities [km/s].
                 vl2D(k,m)=vl(n,i,j,k,m)
                 vp2D(k,m)=vp(n,i,j,k,m)
                 
                 ! The old rho and p calculation based on RCM method:
                 !   Den_IC(n,i,j) = Den_IC(n,i,j)+eta(n,i,j,k,m)/ftv1
                 !   Pressure0     = eta(n,i,j,k,m)*ekev(i,j,k,m)/ftv1
                 ! might be incorrect, giving different results from the 
                 ! following calculation based on integration of flux.
                 
                 ! Number density comes from the integration of "psd":
                 ! n = int(psd*dp^3) = int(flx/p^2*4*pi*p^2*sinA*dpdA)
                 ! with M = p^2/(2*m0*Bm) --> dp = p/2M*dM
                 ! so n = 2*pi*int(flx*p/M*dcosAdM)
                 ! 
                 ! Total pressure and parallel pressure are from
                 !   P    = 4*pi/3*int(E*flx*p/M*dcosAdM)
                 !   Ppar = 4*pi*int(E*flx*p/M*(cosA)^2*dcosAdM)
                 
                 Den_IC(n,i,j) = Den_IC(n,i,j) & 
                      + flx*pp(n,i,j,k,m)/xmm(n,k)*dmm(n,k)*dcosa
                 Pressure0 = ekev(n,i,j,k,m)*flx*pp(n,i,j,k,m)/xmm(n,k)*dmm(n,k)*dcosa
                 Pressure1 = Pressure1 + Pressure0
                 PressurePar1 = PressurePar1 + 3.*Pressure0*cosA2(m)
              enddo
           enddo
           
           Den_IC(n,i,j) = Den_IC(n,i,j)*2*cPi/1.6e-20   ! density in m^-3
           !Coeff = 1.6e-16*2./3.*1.e9                   ! for the old p
           Coeff = 4.*cPi/3.*1.e4*1.e9
           Pressure_IC(n,i,j) = Pressure1*Coeff          ! pressure in nPa
           PressurePar_IC(n,i,j) = PressurePar1*Coeff

!!!! Map flux to fixed energy and pitch-angle grids (energy, sinAo)
           do k=1,neng
              do m=1,npit
!!$                 if ( aloge(n,k) < minval(ekev2D) &
!!$                      .or. aloge(n,k) > maxval(ekev2D) &
!!$                      .or. sinAo(m) < minval(sinA1D) &
!!$                      .or. sinAo(m) > maxval(sinA1D) ) then
!!$                    flx_lo=-50.0
!!$                    vl_lo=0.0
!!$                    vp_lo=0.0
!!$                 else
                 call lintp2aIM(ekev2D,sinA1D,flux2D,nm,nk,aloge(n,k),&
                      sinAo(m),flx_lo)

!!!! Map radial drift to fixed energy and pitch-angle grids (energy, sinAo)
                 call lintp2aIM(ekev2D,sinA1D,vl2D,nm,nk,aloge(n,k),&
                      sinAo(m),vl_lo)

!!!! Map poloidal drift to fixed energy and pitch-angle grids (energy, sinAo)
                 call lintp2aIM(ekev2D,sinA1D,vp2D,nm,nk,aloge(n,k),&
                      sinAo(m),vp_lo)
                 
                 flux(n,i,j,k,m)=10.**flx_lo
                 vlEa(n,i,j,k,m)=vl_lo
                 vpEa(n,i,j,k,m)=re_m*ro(i,j)*1e-3*vp_lo
              enddo
           enddo
        enddo nloop
     enddo iloop1
  enddo jloop1

  ! Calculate pressure of the 'hot' ring current, phot, and temperature, Temp_IC
  jloop2: do j=MinLonPar,MaxLonPar
     iloop2: do i=1,iba(j)
!!!! calculate pitch-angle averaged flux
        do n=1,nspec
           do k=1,neng
              fave(n,i,j,k)=0.
              do m=1,npit
                 fave(n,i,j,k)=fave(n,i,j,k)+flux(n,i,j,k,m)*dmu(m)
              enddo
           enddo
        enddo
!!!! calculate pressure and temperature
        do n=1,nspec
           do k=1,neng
              phot(n,i,j)=phot(n,i,j)+fave(n,i,j,k)*delEE(n,k)*pf(n) ! phot in nPa
           enddo
           Temp_IC(n,i,j)=0.
           if (Den_IC(n,i,j).gt.0.) &
                Temp_IC(n,i,j)=phot(n,i,j)*1.e-9/Den_IC(n,i,j)/1.6e-19   ! eV
        enddo
!!!! calculate parallel pressure
        cosAo2(1:npit) = 1-sinAo(1:npit)**2  !store 1-sinAo^2
        do n=1,nspec
           do k=1,neng
              do m=1,npit
                 Ppar_IC(n,i,j) = Ppar_IC(n,i,j) + flux(n,i,j,k,m) &
                      *cosAo2(m)*dmu(m)*delEE(n,k)*pf(n)*3.
              enddo
           enddo
        enddo
     enddo iloop2
  enddo jloop2


  if (DoCalcFac) then
     ! Calculate field aligned current, fac
     ! First get ghost cell info for eta and ekev when nProc>1
     if (nProc>1) then
        !send ekev ghostcells
        do k=1,nm
           do m=1,nk
             ! call MPI_send(ekev(:,MaxLonPar,k,m),np,MPI_REAL,iProcRight,&
             !      1,iComm,iError)
             ! call MPI_send(ekev(:,MinLonPar,k,m),np,MPI_REAL,iProcLeft,&
             !      2,iComm,iError)
              !recieve ekev from neigboring Procs
             ! call MPI_recv(ekev(:,iLonLeft,k,m),np,MPI_REAL,iProcLeft,&
             !      1,iComm,iStatus_I,iError)
             ! call MPI_recv(ekev(:,iLonRight,k,m),np,MPI_REAL,iProcRight,&
             !      2,iComm,iStatus_I,iError)
              do n=1,nspec
                call MPI_send(ekev(n,:,MaxLonPar,k,m),np,MPI_REAL,iProcRight,&
                   1,iComm,iError)
                call MPI_send(ekev(n,:,MinLonPar,k,m),np,MPI_REAL,iProcLeft,&
                   2,iComm,iError)
              !recieve ekev from neigboring Procs
              call MPI_recv(ekev(n,:,iLonLeft,k,m),np,MPI_REAL,iProcLeft,&
                   1,iComm,iStatus_I,iError)
              call MPI_recv(ekev(n,:,iLonRight,k,m),np,MPI_REAL,iProcRight,&
                   2,iComm,iStatus_I,iError)
                 !send eta ghostcells
                 call MPI_send(eta(n,:,MaxLonPar,k,m),np,MPI_REAL,iProcRight,&
                      1,iComm,iError)
                 call MPI_send(eta(n,:,MinLonPar,k,m),np,MPI_REAL,iProcLeft,&
                      2,iComm,iError)
                 !recieve eta from neigboring Procs
                 call MPI_recv(eta(n,:,iLonLeft,k,m),np,MPI_REAL,iProcLeft,&
                      1,iComm,iStatus_I,iError)
                 call MPI_recv(eta(n,:,iLonRight,k,m),np,MPI_REAL,iProcRight,&
                      2,iComm,iStatus_I,iError)
              enddo
           enddo
        enddo
     endif
     jloop3: do j=MinLonPar,MaxLonPar
        j1=j+1
        j_1=j-1
        if (j1.gt.nt) j1=j1-nt          !! periodic boundary condition
        if (j_1.lt.1) j_1=j_1+nt        !!
        iloop3: do i=2,iba(j)-1
           do k=1,nm
              do m=1,nk
                do n=1,nspec
                 dwkdi=(ekev(n,i+1,j,k,m)-ekev(n,i-1,j,k,m))/(xlatr(i+1)-xlatr(i-1))
                 dwkdj=(ekev(n,i,j1,k,m)-ekev(n,i,j_1,k,m))/(2.*dphi)
                    detadi=(eta(n,i+1,j,k,m)-eta(n,i-1,j,k,m))/(xlatr(i+1)-xlatr(i-1))
                    detadj=(eta(n,i,j1,k,m)-eta(n,i,j_1,k,m))/(2.*dphi)
                    fac(i,j)=fac(i,j)+(detadi*dwkdj-detadj*dwkdi)
                 enddo
              enddo
           enddo
           fac(i,j)=1.6e-16*fac(i,j)/cos(xlatr(i))/rion**2    ! fac in Amp/m^2
        enddo iloop3
     enddo jloop3
  else
     fac(:,:)=0.0
  endif
end subroutine cimi_output


subroutine cimi_precip_calc(rc,dsec)

  use ModCimi,       ONLY: &
       preF, preP, Eje1, &
       xlel=>eChangeOperator_VICI, plel=>pChangeOperator_VICI, &
       OpLossCone_, OpLossCone0_
  use ModCimiTrace, ONLY: iba
  use ModCimiGrid,   ONLY: &
       nProc,iProc,iComm,MinLonPar,MaxLonPar,nt,np,neng,xlatr,xmlt,dlat
  use ModCimiPlanet, ONLY: nspec,re_m
  use ModMPI
  use ModCimiInitialize, ONLY: dphi
  
  implicit none
  
  real :: rc,dsec,dlel,dplel,area,area1,Asec
  integer :: n,i,j,k
  
  preF(1:nspec,1:np,1:nt,1:neng+2)=0.
  preP(1:nspec,1:np,1:nt,1:neng+2)=0.
  Eje1(1:nspec,1:np,1:nt)=0.
  
  area1=rc*rc*re_m*re_m*dphi
  
  do n=1,nspec
     do j=MinLonPar,MaxLonPar
        do i=1,iba(j)
           area=area1*cos(xlatr(i))*dlat(i)            ! area in m^2
           Asec=area*dsec
           do k=1,neng+2
              dlel=xlel(n,i,j,k,OpLossCone_)-xlel(n,i,j,k,OpLossCone0_)
              dplel=plel(n,i,j,k,OpLossCone_)-plel(n,i,j,k,OpLossCone0_)
              if (dlel.lt.0..and.dplel.lt.0.) then
                 preF(n,i,j,k)=-dlel*1.6e-13/Asec     ! E flux in mW/m2
                 preP(n,i,j,k)=-dplel                 ! number of particles
                 
                 ! meanE for E>gride(je)
                 if (k.eq.neng+1) Eje1(n,i,j)=dlel/dplel       
                 
              endif
           enddo
        enddo
     enddo
  enddo

  ! Overwrites the OpLossCone0_ array with the current time information.
  xlel(:,:,:,:,OpLossCone0_) = xlel(:,:,:,:,OpLossCone_)
  plel(:,:,:,:,OpLossCone0_) = plel(:,:,:,:,OpLossCone_)
  
end subroutine cimi_precip_calc


!-------------------------------------------------------------------------------
subroutine FLS_2D(np,nt,iba,fb0,fb1,cl,cp,f2d,fal,fap,fupl,fupp)
!-------------------------------------------------------------------------------
  !  Routine calculates the inter-flux, fal(i+0.5,j) and fap(i,j+0.5), using
  !  2nd order flux limited scheme with super-bee flux limiter method
  !
  !  Input: np,nt,iba,fb0,fb1,cl,cp,f2d
  !  Output: fal,fap
  use ModCimi, ONLY: UseMcLimiter, BetaLimiter
  use ModCimi,       ONLY: MinLonPar,MaxLonPar
  implicit none

  integer np,nt,iba(nt),i,j,j_1,j1,j2,ib
  real cl(np,nt),cp(np,nt),f2d(np,nt),fal(0:np,nt),fap(np,nt),fwbc(0:np+2,nt)
  real fb0(nt),fb1(nt),x,fup,flw,xsign,corr,xlimiter,r

  real,intent(out) :: fupl(0:np,nt), fupp(np,nt)
  fwbc(1:np,1:nt)=f2d(1:np,1:nt)        ! fwbc is f2d with boundary condition

  ! Set up boundary condition
  fwbc(0,1:nt)=fb0(1:nt)
  do j=MinLonPar,MaxLonPar
     ib=iba(j)
     fwbc(ib+1:np+2,j)=fb1(j)
  enddo

  ! find fal and fap
  jloop: do j=MinLonPar,MaxLonPar
     j_1=j-1
     j1=j+1
     j2=j+2
     if (j_1.lt.1) j_1=j_1+nt
     if (j1.gt.nt) j1=j1-nt
     if (j2.gt.nt) j2=j2-nt
     iloop: do i=1,np
        ! find fal
        xsign=sign(1.,cl(i,j))
        fupl(i,j)=0.5*(1.+xsign)*fwbc(i,j)+0.5*(1.-xsign)*fwbc(i+1,j) ! upwind
        flw=0.5*(1.+cl(i,j))*fwbc(i,j)+0.5*(1.-cl(i,j))*fwbc(i+1,j)   ! LW
        x=fwbc(i+1,j)-fwbc(i,j)
        if (abs(x).le.1.e-27) fal(i,j)=fupl(i,j)
        if (abs(x).gt.1.e-27) then
           if (xsign.eq.1.) r=(fwbc(i,j)-fwbc(i-1,j))/x
           if (xsign.eq.-1.) r=(fwbc(i+2,j)-fwbc(i+1,j))/x
           if (r.le.0.) fal(i,j)=fupl(i,j)
           if (r.gt.0.) then
              if(UseMcLimiter)then
                 xlimiter = min(BetaLimiter*r, BetaLimiter, 0.5*(1+r))
              else
                 xlimiter = max(min(2.*r,1.),min(r,2.))
              end if
              corr=flw-fupl(i,j)
              fal(i,j)=fupl(i,j)+xlimiter*corr
           endif
        endif
        ! find fap
        xsign=sign(1.,cp(i,j))
        fupp(i,j)=0.5*(1.+xsign)*fwbc(i,j)+0.5*(1.-xsign)*fwbc(i,j1) ! upwind
        flw=0.5*(1.+cp(i,j))*fwbc(i,j)+0.5*(1.-cp(i,j))*fwbc(i,j1)   ! LW
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
              corr=flw-fupp(i,j)
              fap(i,j)=fupp(i,j)+xlimiter*corr
           endif
        endif
     enddo iloop
  enddo jloop

end subroutine FLS_2D

! OLD LINTP
!!-----------------------------------------------------------------------
!subroutine lintp(xx,yy,n,x,y,ier)
!  !-----------------------------------------------------------------------
!  !  Routine does 1-D interpolation.  xx must be increasing or decreasing
!  !  monotonically.  x is between xx(1) and xx(n)
!  ! 
!  !  input: xx,yy,n,x
!  !  output: y,ier
!
!  implicit none
!
!  integer n,ier,i,jl,ju,jm,j
!  real xx(n),yy(n),x,y,d
!
!  ier = 0
!
!  ! Make sure xx is increasing or decreasing monotonically
!  do i=2,n
!     if (xx(n).gt.xx(1).and.xx(i).lt.xx(i-1)) then
!        write(*,*) ' lintp: xx is not increasing monotonically '
!        write(*,*) n,xx
!        stop
!     endif
!     if (xx(n).lt.xx(1).and.xx(i).gt.xx(i-1)) then
!        write(*,*) ' lintp: xx is not decreasing monotonically '
!        write(*,*) n,xx
!        stop
!     endif
!  enddo
!
!  ! Set ier=1 if out of range
!  if (xx(n).gt.xx(1)) then
!     if (x.lt.xx(1).or.x.gt.xx(n)) ier=1
!  else
!     if (x.gt.xx(1).or.x.lt.xx(n)) ier=1
!  endif
!  if (ier.eq.1) then
!     write(*,*) ' Error: ier.eq.1'
!     print *,'n,x ',n,x
!     print *,'xx(1:n) ',xx(1:n)
!     stop
!  endif
!
!  ! initialize lower and upper values
!  jl=1
!  ju=n
!
!  ! if not done compute a midpoint
!10 if (ju-jl.gt.1) then
!     jm=(ju+jl)/2
!     ! now replace lower or upper limit
!     if ((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm))) then
!        jl=jm
!     else
!        ju=jm
!     endif
!     ! try again
!     go to 10
!  endif
!
!  ! this is the j
!  j=jl      ! if x.le.xx(1) then j=1
!  ! if x.gt.x(j).and.x.le.x(j+1) then j=j
!  ! if x.gt.x(n) then j=n-1
!  d=xx(j+1)-xx(j)
!  y=(yy(j)*(xx(j+1)-x)+yy(j+1)*(x-xx(j)))/d
!
!end subroutine lintp


!-------------------------------------------------------------------------------
subroutine lintp2aIM(x,y,v,nx,ny,x1,y1,v1)
  !-----------------------------------------------------------------------------
  !  This sub program takes 2-d interplation. x is 2-D and y is 1-D.
  !
  !  Input: x,y,v,nx,ny,x1,y1
  !  Output: v1

  implicit none               

  integer nx,ny,j,j1,i,i1,i2,i3
  real x(nx,ny),y(ny),v(nx,ny),x1,y1,v1,a,a1,b,x1d(1000)   ! max(nx)=1000
  real q00,q01,q10,q11

  call locate1IM(y,ny,y1,j)
  j1=j+1
  if (j.eq.0.or.j1.gt.ny) then
     b=1.
     if (j.eq.0) j=j1
     if (j1.gt.ny) j1=j
  else
     b=(y1-y(j))/(y(j+1)-y(j))
  endif

  x1d(1:nx)=x(1:nx,j)

  call locate1IM(x1d,nx,x1,i)
  i1=i+1
  if (i.eq.0.or.i1.gt.nx) then
     a=1.
     if (i.eq.0) i=i1
     if (i1.gt.nx) i1=i
  else
     a=(x1-x1d(i))/(x1d(i+1)-x1d(i))
  endif

  x1d(1:nx)=x(1:nx,j1)
  
  call locate1IM(x1d,nx,x1,i2)
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

end subroutine lintp2aIM

!-------------------------------------------------------------------------------
subroutine lintp2IM(x,y,v,nx,ny,x1,y1,v1)
!-------------------------------------------------------------------------------
!  Routine does 2-D interpolation.  x and y must be increasing or decreasing
!  monotonically
!
  real x(nx),y(ny),v(nx,ny)
  
  call locate1IM(x,nx,x1,i)
  if (i.gt.(nx-1)) i=nx-1      ! extrapolation if out of range
  if (i.lt.1) i=1              ! extrapolation if out of range
  i1=i+1
  a=(x1-x(i))/(x(i1)-x(i))
  
  call locate1IM(y,ny,y1,j)
  if (j.gt.(ny-1)) j=ny-1      ! extrapolation if out of range
  if (j.lt.1) j=1              ! extrapolation if out of range
  j1=j+1
  b=(y1-y(j))/(y(j1)-y(j))
  
  q00=(1.-a)*(1.-b)
  q01=(1.-a)*b
  q10=a*(1.-b)
  q11=a*b
  v1=q00*v(i,j)+q01*v(i,j1)+q10*v(i1,j)+q11*v(i1,j1)
  
  return
end subroutine lintp2IM



!--------------------------------------------------------------------------
subroutine locate1IM(xx,n,x,j)
  !--------------------------------------------------------------------------
  !  Routine return a value of j such that x is between xx(j) and xx(j+1).
  !  xx must be increasing or decreasing monotonically.
  !  If xx is increasing:
  !     If x=xx(m), j=m-1 so if x=xx(1), j=0  and if x=xx(n), j=n-1
  !     If x < xx(1), j=0  and if x > xx(n), j=n
  !  If xx is decreasing:
  !     If x=xx(m), j=m so if x=xx(1), j=1  and if x=xx(n), j=n
  !     If x > xx(1), j=0  and if x < xx(n), j=n
  !
  !  Input: xx,n,x
  !  Output: j

  implicit none

  integer n,j,i,jl,ju,jm
  real xx(n),x

  ! Make sure xx is increasing or decreasing monotonically
  do i=2,n
     if (xx(n).gt.xx(1).and.xx(i).lt.xx(i-1)) then
        write(*,*) ' locate1IM: xx is not increasing monotonically '
        write(*,*) n, (xx(j),j=1,n)
        call CON_STOP('CIMI stopped in locate1IM')
     endif
     if (xx(n).lt.xx(1).and.xx(i).gt.xx(i-1)) then
        write(*,*) ' locate1IM: xx is not decreasing monotonically '
        write(*,*) ' n, xx  ',n,xx
        call CON_STOP('CIMI stopped in locate1IM')
     endif
  enddo

  jl=0
  ju=n+1
  test: do
     if (ju-jl.le.1) exit test
     jm=(ju+jl)/2
     if ((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm))) then
        jl=jm
     else
        ju=jm
     endif
  end do test
  j=jl

end subroutine locate1IM


!-------------------------------------------------------------------------------
        subroutine lintp3IM(x,y,z,v,nx,ny,nz,x1,y1,z1,v1)
!-------------------------------------------------------------------------------
!  This sub program takes 3-d interplation.
!
      real x(nx),y(ny),z(nz),v(nx,ny,nz)

      call locate1IM(x,nx,x1,i)
      if (i.gt.(nx-1)) i=nx-1      ! extrapolation if out of range
      if (i.lt.1) i=1              ! extrapolation if out of range

      call locate1IM(y,ny,y1,j)
      if (j.gt.(ny-1)) j=ny-1      ! extrapolation if out of range
      if (j.lt.1) j=1              ! extrapolation if out of range

      call locate1IM(z,nz,z1,k)
      if (k.gt.(nz-1)) k=nz-1      ! extrapolation if out of range
      if (k.lt.1) k=1              ! extrapolation if out of range


      i1=i+1
      j1=j+1
      k1=k+1
      a=(x1-x(i))/(x(i1)-x(i))
      b=(y1-y(j))/(y(j1)-y(j))
      c=(z1-z(k))/(z(k1)-z(k))

      q000=(1.-a)*(1.-b)*(1.-c)*v(i,j,k)
      q001=(1.-a)*(1.-b)*c*v(i,j,k1)
      q010=(1.-a)*b*(1.-c)*v(i,j1,k)
      q011=(1.-a)*b*c*v(i,j1,k1)
      q100=a*(1.-b)*(1.-c)*v(i1,j,k)
      q101=a*(1.-b)*c*v(i1,j,k1)
      q110=a*b*(1.-c)*v(i1,j1,k)
      q111=a*b*c*v(i1,j1,k1)
      v1=q000+q001+q010+q011+q100+q101+q110+q111

      end subroutine lintp3IM


!-----------------------------------------------------------------------------

      subroutine tridagIM(a,b,c,r,u,n,ier)
!-----------------------------------------------------------------------------

      parameter (nmax=100)
      real gam(nmax),a(n),b(n),c(n),r(n),u(n),bet
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
      do 11 j=2,n
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
11    continue
!
!    back substitution
!
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue

      end subroutine tridagIM

!-------------------------------------------------------------------------------












!Old CLOSED SUBROUTINE
!!--------------------------------------------------------------------------
!subroutine closed(n1,n2,yy,dx,ss)
!  !--------------------------------------------------------------------------
!  ! Routine does numerical integration using closed form.
!  ! 
!  ! Input: n1,n2,yy,dx
!  ! Output: ss
!
!  implicit none
!
!  integer n1,n2,i
!  real yy(n2),dx(n2),ss
!
!  ss=0.
!  do i=n1,n2
!     ss=ss+yy(i)*dx(i)
!  enddo
!
!end subroutine closed

