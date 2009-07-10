^CFG COPYRIGHT UM
!==============================================================================
module ModUser
  use ModReadParam, ONLY: lStringLine
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_init_session,               &
       IMPLEMENTED3 => user_set_ics,                    &
       IMPLEMENTED4 => user_initial_perturbation,       &
       IMPLEMENTED5 => user_face_bcs,                   &
       IMPLEMENTED6 => user_get_log_var,                &
       IMPLEMENTED7 => user_get_b0,                     &
       IMPLEMENTED8 => user_update_states,              &
       IMPLEMENTED9 => user_specify_refinement,         &
       IMPLEMENTED10=> user_set_boundary_cells,         &
       IMPLEMENTED11=> user_set_plot_var

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'Empirical Solar Wind and MHD Turbulence'

  character(len=lStringLine) :: NameModel
  character(len=lStringLine) :: NameModelSpectrum

contains
  !============================================================================
  subroutine user_read_inputs
    use EEE_ModMain,    ONLY: EEE_set_parameters
    use ModMain
    use ModProcMH,      ONLY: iProc
    use ModReadParam,   ONLY: read_line, read_command, read_var
    use ModIO,          ONLY: write_prefix, write_myname, iUnitOut
    use ModMagnetogram, ONLY: set_parameters_magnetogram
    implicit none

    character (len=100) :: NameCommand
    !--------------------------------------------------------------------------
    UseUserInitSession = .true.

    if(iProc == 0 .and. lVerbose > 0)then
       call write_prefix;
       write(iUnitOut,*)'User read_input TURBULENCE CORONA starts'
    endif

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)
       case("#PFSSM")
          call read_var('UseUserB0'  ,UseUserB0)
          if(UseUserB0)then
             call set_parameters_magnetogram
             call read_var('dt_UpdateB0',dt_UpdateB0)
             DoUpdateB0 = dt_updateb0 > 0.0
          end if

       case("#ARCH","#TD99FLUXROPE","#GL98FLUXROPE")
          call EEE_set_parameters(NameCommand)

       case("#EMPIRICALSW")
          call read_var('NameModel',NameModel)

       case('#USERINPUTEND')
          if(iProc == 0 .and. lVerbose > 0)then
             call write_prefix;
             write(iUnitOut,*)'User read_input TURBULENCE CORONA ends'
          endif
          EXIT

       case default
          if(iProc == 0) then
             call write_myname; write(*,*) &
                  'ERROR: Invalid user defined #COMMAND in user_read_inputs. '
             write(*,*) '--Check user_read_inputs for errors'
             write(*,*) '--Check to make sure a #USERINPUTEND command was used'
             write(*,*) '  *Unrecognized command was: '//NameCommand
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if
       end select
    end do

  end subroutine user_read_inputs

  !============================================================================

  subroutine user_init_session
    use EEE_ModMain,    ONLY: EEE_initialize
    use ModIO,          ONLY: write_prefix, iUnitOut
    use ModMagnetogram, ONLY: read_magnetogram_file
    use ModMain,        ONLY: UseUserB0
    use ModPhysics,     ONLY: BodyNDim_I,BodyTDim_I,g
    use ModProcMH,      ONLY: iProc
    use ModReadParam,   ONLY: i_line_command
    implicit none
    !--------------------------------------------------------------------------
    if(iProc == 0)then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'user_init_session:'
       call write_prefix; write(iUnitOut,*) ''
    end if

    if(i_line_command("#PFSSM", iSessionIn = 1) < 0)then
       write(*,*) 'In session 1, a magnetogram file has to be read via #PFSSM'
       call stop_mpi('ERROR: Correct PARAM.in!')
    end if
    if(i_line_command("#PFSSM") > 0)then
       call read_magnetogram_file
    end if

    if(i_line_command("#EMPIRICALSW", iSessionIn = 1) < 0)then
       write(*,*) 'An empirical model has to be set via #EMPIRICALSW'
       call stop_mpi('ERROR: Correct PARAM.in!')
    end if
    if(i_line_command("#EMPIRICALSW") > 0)then
       call set_empirical_model(trim(NameModel),BodyTDim_I(1))
    end if

 
    call EEE_initialize(BodyNDim_I(1),BodyTDim_I(1),g)

    if(iProc == 0)then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'user_init_session finished'
       call write_prefix; write(iUnitOut,*) ''
    end if

  end subroutine user_init_session

  !============================================================================
  subroutine user_face_bcs(VarsGhostFace_V)
    use EEE_ModMain,   ONLY: EEE_get_state_BC
    use ModSize,       ONLY: East_,West_,South_,North_,Bot_,Top_,nDim
    use ModMain,       ONLY: time_accurate,x_,y_,z_, UseRotatingFrame, n_step, Iteration_Number
    use ModVarIndexes !ONLY: nVar,Ew_,rho_,Ux_,Uy_,Uz_,Bx_,By_,Bz_,P_
    use ModAdvance,    ONLY: State_VGB
    use ModPhysics,    ONLY: inv_gm1,OmegaBody,Si2No_V, &
         UnitB_,UnitU_,UnitRho_,UnitP_
    use ModNumConst,   ONLY: cTolerance,cTiny
    use ModFaceBc, ONLY: FaceCoords_D, VarsTrueFace_V, TimeBc, &
         iFace, jFace, kFace, iSide, iBlockBc
    implicit none

    real, intent(out):: VarsGhostFace_V(nVar)

    integer:: iCell,jCell,kCell

    real:: DensCell,PresCell,GammaCell,TBase,B1dotR  
    real, dimension(3):: RFace_D,B1_D,U_D,B1t_D,B1n_D

    real :: RhoCME,UCME_D(nDim),BCME_D(nDim),pCME
    real :: BCMEn,BCMEn_D(nDim),UCMEn,UCMEn_D(nDim),UCMEt_D(nDim)
    !--------------------------------------------------------------------------

    RFace_D  = FaceCoords_D/sqrt(sum(FaceCoords_D**2))

    U_D (x_:z_)  = VarsTrueFace_V(Ux_:Uz_)
    B1_D(x_:z_)  = VarsTrueFace_V(Bx_:Bz_)
    B1dotR       = dot_product(RFace_D,B1_D)
    B1n_D(x_:z_) = B1dotR*RFace_D(x_:z_)
    B1t_D        = B1_D-B1n_D

    !\
    ! Update BCs for velocity and induction field::
    !/
    VarsGhostFace_V(Ux_:Uz_) = -U_D(x_:z_)
    VarsGhostFace_V(Bx_:Bz_) = B1t_D(x_:z_)!-B1n_D(x_:z_)

    !\
    ! Compute the perturbed state of the eruptive event at RFace_D::
    !/
    call EEE_get_state_BC(RFace_D,RhoCME,UCME_D,BCME_D,pCME,TimeBc, &
         n_step,iteration_number)

    RhoCME = RhoCME*Si2No_V(UnitRho_)
    UCME_D = UCME_D*Si2No_V(UnitU_)
    BCME_D = BCME_D*Si2No_V(UnitB_)
    pCME = pCME*Si2No_V(UnitP_)

    !\
    ! Fix the normal component of the CME field to BCMEn_D at the Sun::
    !/
    BCMEn   = dot_product(RFace_D,BCME_D)
    BCMEn_D = BCMEn*RFace_D
    VarsGhostFace_V(Bx_:Bz_) = VarsGhostFace_V(Bx_:Bz_) + BCMEn_D

    !\
    ! Fix the tangential components of the CME velocity at the Sun
    !/
    UCMEn   = dot_product(RFace_D,UCME_D)
    UCMEn_D = UCMEn*RFace_D
    UCMEt_D = UCME_D-UCMEn_D   
    VarsGhostFace_V(Ux_:Uz_) = VarsGhostFace_V(Ux_:Uz_) + 2.0*UCMEt_D

    !\
    ! Update BCs for the mass density, EnergyRL,
    ! and pressure::
    !/
    iCell = iFace; jCell = jFace; kCell = kFace
    select case(iSide)
    case(East_)
       iCell  = iFace
    case(West_)
       iCell  = iFace-1
    case(South_)
       jCell  = jFace
    case(North_)
       jCell  = jFace-1
    case(Bot_)
       kCell  = kFace
    case(Top_)
       kCell  = kFace-1
    case default
       write(*,*)'ERROR: iSide = ',iSide
       call stop_mpi('incorrect iSide value in user_face_bcs')
    end select

    call get_plasma_parameters_cell(iCell,jCell,kCell,iBlockBc,&
         DensCell,PresCell,GammaCell)
    VarsGhostFace_V(Rho_) = &
         max(-VarsTrueFace_V(Rho_) + 2.0*(DensCell+RhoCME), &
         VarsTrueFace_V(Rho_))
    TBase = (PresCell+pCME)/(DensCell+RhoCME)
    VarsGhostFace_V(P_) = &
         max(VarsGhostFace_V(Rho_)*TBase, &
         VarsTrueFace_V(P_))
    VarsGhostFace_V(Ew_) = &!max(-VarsTrueFace_V(Ew_)+ &
         VarsGhostFace_V(Rho_)*TBase &
         *(1.0/(GammaCell-1.0)-inv_gm1)

    !\
    ! Apply corotation
    !/
    if (.not.UseRotatingFrame) then
       VarsGhostFace_V(Ux_) = VarsGhostFace_V(Ux_) &
            - 2*OmegaBody*FaceCoords_D(y_)
       VarsGhostFace_V(Uy_) = VarsGhostFace_V(Uy_) &
            + 2*OmegaBody*FaceCoords_D(x_)
    end if
  end subroutine user_face_bcs

  !============================================================================

  subroutine get_plasma_parameters_cell(iCell,jCell,kCell,iBlock,&
       DensCell,PresCell,GammaCell)

    ! This subroutine computes the cell values for density and pressure 
    ! assuming an isothermal atmosphere
    
    use ModGeometry,   ONLY: x_BLK,y_BLK,z_BLK,R_BLK
    use ModNumConst
    use ModPhysics,    ONLY: GBody,BodyRho_I,Si2No_V,UnitTemperature_
    use ModExpansionFactors,  ONLY: UMin,T0
    implicit none

    integer, intent(in)  :: iCell,jCell,kCell,iBlock
    real, intent(out)    :: DensCell,PresCell,GammaCell
    real :: UFinal       !The solar wind speed at the far end of the Parker spiral,
                         !which originates from the given cell
    real :: URatio       !The coronal based values for temperature density 
                         !are scaled as functions of UFinal/UMin ratio
    real :: Temperature
    !--------------------------------------------------------------------------

    call get_gamma_emp(x_BLK(iCell,jCell,kCell,iBlock),&
         y_BLK(iCell,jCell,kCell,iBlock),&
         z_BLK(iCell,jCell,kCell,iBlock),&
         GammaCell)
    call get_bernoulli_integral(x_BLK(iCell,jCell,kCell,iBlock)/&
         R_BLK(iCell,jCell,kCell,iBlock),&
         y_BLK(iCell,jCell,kCell,iBlock)/R_BLK(iCell,jCell,kCell,iBlock),&
         z_BLK(iCell,jCell,kCell,iBlock)/R_BLK(iCell,jCell,kCell,iBlock),UFinal)
    URatio=UFinal/UMin

    !This is the temperature variation
    Temperature = T0/(min(URatio,2.0))*Si2No_V(UnitTemperature_)

    DensCell  = ((1.0/URatio)**2) &          !This is the density variation
         *BodyRho_I(1)*exp(-GBody/Temperature &
         *(1.0/max(R_BLK(iCell,jCell,kCell,iBlock),0.90)-1.0))

    PresCell = DensCell*Temperature

  end subroutine get_plasma_parameters_cell

  !============================================================================

  subroutine user_initial_perturbation
    use ModMain, ONLY: nBLK,unusedBLK,x_,y_,z_,n_step,iteration_number
    use ModVarIndexes
    use ModAdvance,   ONLY: State_VGB 
    use ModPhysics,   ONLY: Si2No_V,UnitB_
    use ModGeometry
    use ModMagnetogram, ONLY: get_magnetogram_field
   
    implicit none

    logical :: oktest,oktest_me

    !--------------------------------------------------------------------------
 
        call set_oktest('user_initial_perturbation',oktest,oktest_me)

        call set_wave_spectrum
           
        write(*,*) 'SC: Finished initializing wave spectrum'
 
     
   end subroutine user_initial_perturbation
  !============================================================================

  subroutine user_set_ics
    use ModMain,      ONLY: globalBLK,nI,nJ,nK
    use ModVarIndexes
    use ModAdvance,   ONLY: State_VGB 
    use ModPhysics,   ONLY: inv_gm1,BodyTDim_I
    use ModGeometry
    implicit none

    integer :: i,j,k,iBLK
    logical :: oktest,oktest_me
    real :: Dens_BLK,Pres_BLK,Gamma_BLK
    real :: x,y,z,R,ROne,Rmax,U0
    !--------------------------------------------------------------------------
    call set_oktest('user_set_ics',oktest,oktest_me)

    iBLK = globalBLK

    select case(TypeGeometry)
    case('cartesian')
       Rmax = max(2.1E+01,sqrt(x2**2+y2**2+z2**2))
    case('spherical_lnr')
       Rmax = max(2.1E+01,exp(XyzMax_D(1)))
    end select

    ! The sqrt is for backward compatibility with older versions of the Sc
    U0 = 4.0*sqrt(2.0E+6/BodyTDim_I(1))

    do k=1,nK; do j=1,nJ; do i=1,nI
       x = x_BLK(i,j,k,iBLK)
       y = y_BLK(i,j,k,iBLK)
       z = z_BLK(i,j,k,iBLK)
       R = max(R_BLK(i,j,k,iBLK),cTolerance)
       ROne = max(1.0,R)
       State_VGB(Bx_:Bz_,i,j,k,iBLK) = 0.0
       call get_plasma_parameters_cell(i,j,k,iBLK,&
            Dens_BLK,Pres_BLK,Gamma_BLK)
       State_VGB(rho_,i,j,k,iBLK) = Dens_BLK
       State_VGB(P_,i,j,k,iBLK)   = Pres_BLK
       State_VGB(RhoUx_,i,j,k,iBLK) = Dens_BLK &
            *U0*((ROne-1.0)/(Rmax-1.0))*x/R
       State_VGB(RhoUy_,i,j,k,iBLK) = Dens_BLK &
            *U0*((ROne-1.0)/(Rmax-1.0))*y/R
       State_VGB(RhoUz_,i,j,k,iBLK) = Dens_BLK &
            *U0*((ROne-1.0)/(Rmax-1.0))*z/R
       State_VGB(Ew_,i,j,k,iBLK) = Pres_BLK &
            *(1.0/(Gamma_BLK-1.0)-inv_gm1) 
    end do; end do; end do

  end subroutine user_set_ics

  !============================================================================
  subroutine user_get_b0(xInput,yInput,zInput,B0_D)
    use EEE_ModMain,    ONLY: EEE_get_B0
    use ModPhysics,     ONLY: Io2No_V,Si2No_V,UnitB_
    use ModMagnetogram, ONLY: get_magnetogram_field

    implicit none

    real, intent(in):: xInput,yInput,zInput
    real, intent(out), dimension(3):: B0_D
   
    real :: x_D(3),B_D(3)
    !--------------------------------------------------------------------------

    call get_magnetogram_field(xInput,yInput,zInput,B0_D)
    B0_D = B0_D*Si2No_V(UnitB_)

    x_D = (/ xInput, yInput, zInput /)
    call EEE_get_B0(x_D,B_D)
    B0_D = B0_D + B_D*Si2No_V(UnitB_)

  end subroutine user_get_b0

  !============================================================================

  subroutine user_update_states(iStage,iBlock)
    use ModMain,    ONLY: iteration_number
    use ModVarIndexes
    use ModSize
    use ModAdvance, ONLY: State_VGB, B0_DGB
    use ModMain,    ONLY: nStage
    use ModPhysics, ONLY: inv_gm1
    use ModGeometry,ONLY: R_BLK
    use ModEnergy,  ONLY: calc_energy_cell
    use ModExpansionFactors, ONLY: gammaSS

    implicit none

    integer,intent(in)                  :: iStage,iBlock
    integer                             :: i,j,k
    real                                :: DensCell,PresCell,GammaCell,Beta
    !--------------------------------------------
    call update_states_MHD(iStage,iBlock)
    !\
    ! Begin update of pressure and relaxation energy::
    !/
    !  if (iStage/=nStage) return
    do k=1,nK; do j=1,nJ; do i=1,nI
       call get_plasma_parameters_cell(i,j,k,iBlock,&
            DensCell,PresCell,GammaCell)
       if(R_BLK(i,j,k,iBlock)>2.5)&
            GammaCell=GammaCell-(GammaCell-gammaSS)*max(0.0, &
            -1.0 + 2*State_VGB(P_,i,j,k,iBlock)/&
            (State_VGB(P_   ,i,j,k,iBlock)+sum(&
            (State_VGB(Bx_:Bz_ ,i,j,k,iBlock)+B0_DGB(:,i,j,k,iBlock))**2)&
            *0.25*(R_BLK(i,j,k,iBlock)/2.5)**1.50))
       State_VGB(P_   ,i,j,k,iBlock)=(GammaCell-1.0)*      &
            (inv_gm1*State_VGB(P_,i,j,k,iBlock) + State_VGB(Ew_,i,j,k,iBlock))
       State_VGB(Ew_,i,j,k,iBlock)= State_VGB(P_,i,j,k,iBlock) &
            *(1.0/(GammaCell-1.0)-inv_gm1)
    end do; end do; end do
    call calc_energy_cell(iBlock)
    !\
    ! End update of pressure and relaxation energy::
    !/
    
  end subroutine user_update_states

!=======================================================================
   subroutine calc_cutoff_freq(FreqCutOff_C,iBLK)

     ! This subroutine calculates the cut-off frequency for Alfven waves, which is the ion cyclotron frequency, 
     ! \omega_{c.o.}=zeB/m. For ions z=1. 
     
     use ModMain,              ONLY: nI,nJ,nK
     use ModVarIndexes
     use ModAdvance,           ONLY: State_VGB
     use ModConst,             ONLY: cElectronCharge, cProtonMass
     use ModSize

     implicit none
     
     integer                                       :: i,j,k
     real, dimension(nI,nJ,nK,nBLK)                :: Btot_CB
     real, intent(out), dimension(nI,nj,nK)        :: FreqCutOff_C
     integer, intent(in)                           :: iBLK
     ! -----------------------------------------------------------------
     
     do i=1,nI; do j=1,nJ ; do k=1,nK
        Btot_CB(i,j,k,iBLK)=(State_VGB(Bx_,i,j,k,iBLK)**2+  &
             State_VGB(By_,i,j,k,iBLK)**2+State_VGB(Bz_,i,j,k,iBLK)**2)**0.5
        FreqCutOff_C(i,j,k)=cElectronCharge*Btot_CB(i,j,k,iBLK)/cProtonMass
     end do; end do ; end do
    

   end subroutine calc_cutoff_freq
!=======================================================================
 subroutine set_wave_spectrum

! This subroutine sets the Alfven wave spectrum at the surface of the sun, up to 1.01 Rs.
! The frequancy range of the spectrum is determined by Alfven waves frequencie observed at the chromosphere 
! and around 1AU. The maximal energy density is determined by solar observations from about 1000km above the limb
! by Hinode. Thus the wave spectrum initial conditions are set only up tp 1.01Rs (~1000km).
! The spectrum can then by further modofoed empirically in order to match the observed spectrum at 1AU.

! Further developement of this model will enable setting wave spectrum properties via PARAM.in.
! As explained in ModEquationMhdTurb.f90, the current version only allows advection of wave energy with the background fluid.
! Further development will include evolution of the wave spectrum (via dissipation/growth, plasma expansion/contraction etc)
  
   use ModMain,      ONLY: nI,nJ,nK,x_,y_,z_,nBLK
   use ModVarIndexes
   use ModAdvance,   ONLY: State_VGB, B0_DGB
   use ModGeometry
   use ModSize
   use ModConst,     ONLY: cMu,cLightSpeed,cElectronCharge, cGEV,cAU
   use ModNumConst,  ONLY: cPi,cTiny
   use ModPhysics,     ONLY: Io2No_V,Si2No_V,UnitB_
   use ModMagnetogram, ONLY: get_magnetogram_field
   
   implicit none

   ! Spatial grid variables
   integer                      :: i,j,k,iBLK
   real                         :: x,y,z,R
   
   ! 1D Frequency grid variables
   real                         :: FreqMin, FreqMax, dLogFreq
   real, dimension(I50_-I01_+1) :: WaveEnergy_I, WaveLogFreq_I
   integer                      :: iFreq
   real,dimension(nI,nJ,nK)     :: FreqCutOff_C

   ! Initial spectrum model parameters
   !the intensity of outward travelling waves (initial condition)
   real,parameter               :: Alpha=1.0/10.0
   real,parameter               :: Lambda0=4.0/10.0  ![AU]
   real,parameter               :: KolmogPower=-5.0/3.0
   real,parameter               :: ConstCoeff=(1-Alpha)*54/(4*cPi*Lambda0)
   !ConstCoeff is a constant coeffiecient appearing in the initial spectrum
   real                         :: bCoeff ,rCoeff ! B and r dependence of initial spectrum
   real,dimension(nI,nJ,nK,nBLK):: EnergyCoeff_CB ! product of 3 previous coefficients
   real                         :: EnergyCoeffMax, VAlfvenCell, B0Cell,RhoCell
   real, dimension(3)           :: B0_D, B_D
   
   ! ------------------------------------------------------------------

   ! \
   ! Set frequency grid
   ! /
   
   FreqMin = 6e-6  ! minimal frequency of Alfven waves spectrum
   !This is the frequancy of the variable I01_
   FreqMax = 0.001 ! maximal frequency of Alfven waves spectrum
   !This is the frequency of the variable I50_
   
   ! calculate frequency interval on a natural logarithmic scale 
   dLogFreq = (log(FreqMax)-log(FreqMin))/(I50_-I01_) 
  
   ! Divide the spectrum into frequncy groups on a log scale
   do iFreq = 1,I50_-I01_+1

      WaveLogFreq_I(iFreq)=log(FreqMin)+(iFreq-1)*dLogFreq

   end do
    
   ! \
   ! Calculate wave energy  coefficients in all cells 
   ! /
   do iBLK=1,nBLK

      do k=1,nK; do j=1,nJ; do i=1,nI
      
         if (R_BLK(i,j,k,iBLK) .ge. 1.0) then

            x=x_BLK(i,j,k,iBLK)
            y=y_BLK(i,j,k,iBLK)
            z=z_BLK(i,j,k,iBLK)
            
            !Get B0, the initial potential field
            call get_magnetogram_field(x,y,z,B0_D)
            B0_D(1) = B0_D(1)*Si2No_V(UnitB_)
            B0_D(2)= B0_D(2)*Si2No_V(UnitB_)
            B0_D(3)= B0_D(3)*Si2No_V(UnitB_)

            !B_D(1) = B0_DGB(x_,i,j,k,iBLK)*Si2No_V(UnitB_)
            !B_D(2)= B0_DGB(y_,i,j,k,iBLK)*Si2No_V(UnitB_)
            !B_D(3)= B0_DGB(z_,i,j,k,iBLK)*Si2No_V(UnitB_)
            
            B0Cell=((B0_D(1))**2+(B0_D(2))**2 + &
                 (B0_D(3))**2)**0.5
                  
  
            ! The wave energy spectrum is described by the power law:
            ! I_{+}= Coeff*bCoeff*rCoeff*(w/Va)^{-5/3}
            ! where Va is the Alfven speed Va=B^2/sqrt{\mu_0\rho}
            ! ConstCoeff=54/4\pi\lambda0
            ! bCoeff=B^2(eBc/1GeV)^{-1/3}
            ! rCoeff=(r/1AU)^{-1}   
            
            bCoeff=B0Cell**2*(cElectronCharge*B0Cell*cLightSpeed/cGEV)&
                 **(-(1.0/3.0))
         
            rCoeff=cAU/R_BLK(i,j,k,iBLK)
         
            RhoCell=State_VGB(rho_,i,j,k,iBLK)
            VAlfvenCell=B0Cell/(RhoCell*cMu)**0.5
         
            EnergyCoeff_CB(i,j,k,iBLK)=ConstCoeff*bCoeff*rCoeff*(VAlfvenCell)**(-KolmogPower)
          
         end if

      end do; end do; end do
   end do
   
   !\
   ! Initialize spectrum in all frequency groups
   !/

   EnergyCoeffMax=maxval(EnergyCoeff_CB)
   do iBLK=1,nBLK
      do i=1,nI; do j=1,nJ; do k=1,nK
         if (R_BLK(i,j,k,iBLK)>1.0) then
            
            ! Normalize EnergyCoeff_CB according to maximum value in the spatial domain
            EnergyCoeff_CB(i,j,k,iBLK)=EnergyCoeff_CB(i,j,k,iBLK)/EnergyCoeffMax
            
            ! Start filling frequency groups
            do iFreq=1,I50_-I01_+1
   
               WaveEnergy_I(iFreq)=EnergyCoeff_CB(i,j,k,iBLK)*exp(WaveLogFreq_I(iFreq))**(KolmogPower)
               State_VGB(iFreq-1+I01_,i,j,k,iBLK)=WaveEnergy_I(iFreq)
               
            end do
         else
            State_VGB(I01_:I50_,i,j,k,iBLK)=cTiny
         end if
   
      end do; end do ; end do
   end do

   !\
   ! Calculate cut-off frequancy for all cells
   !/
   do iBLK=1,nBLK
      call calc_cutoff_freq(FreqCutOff_C,iBLK)
      write(*,*) 'Maximum cut-off = ',maxval(FreqCutOff_C)
   end do
  end subroutine set_wave_spectrum
!========================================================================

  subroutine user_get_log_var(VarValue,TypeVar,Radius)

    use ModIO,         ONLY: write_myname
    use ModMain,       ONLY: unusedBLK,nBLK,x_,y_,z_
    use ModVarIndexes !ONLY: Ew_,Bx_,By_,Bz_,rho_,rhoUx_,rhoUy_,rhoUz_,P_ 
    use ModGeometry,   ONLY: R_BLK
    use ModAdvance,    ONLY: State_VGB,tmp1_BLK,B0_DGB
    use ModPhysics,    ONLY: inv_gm1,&
         No2Si_V,UnitEnergydens_,UnitX_,UnitU_,UnitRho_

    real, intent(out):: VarValue
    character (LEN=10), intent(in):: TypeVar 
    real, intent(in), optional :: Radius
    !
    integer:: iBLK
    real:: unit_energy,unit_mass
    real, external:: integrate_BLK
    !--------------------------------------------------------------------------
    unit_energy = 1.0e7*No2Si_V(UnitEnergydens_)*No2Si_V(UnitX_)**3
    unit_mass   = 1.0e3*No2Si_V(UnitRho_)*No2Si_V(UnitX_)**3
    !\
    ! Define log variable to be saved::
    !/
    select case(TypeVar)
    case('em_t','Em_t','em_r','Em_r')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(:,:,:,iBLK) = & 
               (B0_DGB(x_,:,:,:,iBLK)+State_VGB(Bx_,:,:,:,iBLK))**2+&
               (B0_DGB(y_,:,:,:,iBLK)+State_VGB(By_,:,:,:,iBLK))**2+&
               (B0_DGB(z_,:,:,:,iBLK)+State_VGB(Bz_,:,:,:,iBLK))**2
       end do
       VarValue = unit_energy*0.5*integrate_BLK(1,tmp1_BLK)
    case('ek_t','Ek_t','ek_r','Ek_r')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(:,:,:,iBLK) = &
               (State_VGB(rhoUx_,:,:,:,iBLK)**2 +&
               State_VGB(rhoUy_,:,:,:,iBLK)**2 +&
               State_VGB(rhoUz_,:,:,:,iBLK)**2)/&
               State_VGB(rho_  ,:,:,:,iBLK)             
       end do
       VarValue = unit_energy*0.5*integrate_BLK(1,tmp1_BLK)
    case('et_t','Et_t','et_r','Et_r')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(:,:,:,iBLK) = State_VGB(P_,:,:,:,iBLK)
       end do
       VarValue = unit_energy*inv_gm1*integrate_BLK(1,tmp1_BLK)
    case('ew_t','Ew_t','ew_r','Ew_r')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(:,:,:,iBLK) = State_VGB(Ew_,:,:,:,iBLK)
       end do
       VarValue = unit_energy*integrate_BLK(1,tmp1_BLK)
    case('ms_t','Ms_t')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(:,:,:,iBLK) = &
               State_VGB(rho_,:,:,:,iBLK)/R_BLK(:,:,:,iBLK)
       end do
       VarValue = unit_mass*integrate_BLK(1,tmp1_BLK)
    case('vol','Vol')
       tmp1_BLK(:,:,:,iBLK) = 1.0
       VarValue = integrate_BLK(1,tmp1_BLK)
    case default
       VarValue = -7777.
       call write_myname;
       write(*,*) 'Warning in set_user_logvar: unknown logvarname = ',TypeVar
    end select
  end subroutine user_get_log_var
  
 !===========================================================================
  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModSize, ONLY: nI, nJ, nK
    use ModAdvance, ONLY: State_VGB
    use ModVarIndexes
    use ModPhysics, ONLY: No2Si_V,UnitEnergydens_,UnitX_

    implicit none

    integer,          intent(in)   :: iBlock
    character(len=*), intent(in)   :: NameVar
    logical,          intent(in)   :: IsDimensional
    real,             intent(out)  :: PlotVar_G(-1:nI+2, -1:nJ+2, -1:nK+2)
    real,             intent(out)  :: PlotVarBody
    logical,          intent(out)  :: UsePlotVarBody
    character(len=*), intent(inout):: NameTecVar
    character(len=*), intent(inout):: NameTecUnit
    character(len=*), intent(inout):: NameIdlUnit
    logical,          intent(out)  :: IsFound

    character (len=*), parameter :: NameSub = 'user_set_plot_var'
    real :: unit_energy, IntIwCell,dw
    integer :: i,j,k, I_Index
    logical :: IsError
    !-------------------------------------------------------------------    
    !UsePlotVarBody = .true.
    !PlotVarBody = 0.0
    IsFound=.true.

    !\                                                                              
    ! Define plot variable to be saved::
    !/ 
    ! IntIwCell is the integral of the wave energy over frequency
    select case(NameVar)
    case('IntIw')
       PlotVar_G(i,j,k)=IntIwCell
        case default
       IsFound= .false.
    end select
  end subroutine user_set_plot_var

!=========================================================================== 
  subroutine user_specify_refinement(iBlock, iArea, DoRefine)

    use ModSize,     ONLY: nI, nJ, nK
    use ModAdvance,  ONLY: State_VGB, Bx_, By_, Bz_, B0_DGB
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, far_field_BCs_BLK
    use ModNumConst, ONLY: cTiny
    use ModMain,     ONLY: x_, y_, z_

    integer, intent(in) :: iBlock, iArea
    logical,intent(out) :: DoRefine

    real :: rDotB_G(1:nI,1:nJ,0:nK+1)
    integer :: i,j,k
    character (len=*), parameter :: NameSub = 'user_specify_refinement'
    !-------------------------------------------------------------------

    if(far_field_BCs_BLK(iBlock))then
       DoRefine = .false.
       RETURN
    end if

    ! Calculate r.B in all physical cells and ghost cells 
    ! in the Z/Theta direction to find current sheet 
    ! passing between blocks
    do k=0, nK+1; do j=1, nJ; do i=1, nI
       rDotB_G(i,j,k) = x_BLK(i,j,k,iBlock)   &
            * (B0_DGB(x_,i,j,k,iBlock) + State_VGB(Bx_,i,j,k,iBlock)) &
            +              y_BLK(i,j,k,iBlock)   &
            * (B0_DGB(y_,i,j,k,iBlock) + State_VGB(By_,i,j,k,iBlock)) &
            +              z_BLK(i,j,k,iBlock)   &
            * (B0_DGB(z_,i,j,k,iBlock) + State_VGB(Bz_,i,j,k,iBlock))
    end do; end do; end do;

    DoRefine = maxval(rDotB_G) > cTiny .and. minval(rDotB_G) < -cTiny

  end subroutine user_specify_refinement

  !============================================================================
  subroutine user_set_boundary_cells(iBLK)
    use ModGeometry,      ONLY: ExtraBc_, IsBoundaryCell_GI, r_Blk
    use ModBoundaryCells, ONLY: SaveBoundaryCells
    use ModPhysics,       ONLY: rBody
    implicit none

    integer, intent(in):: iBLK

    character (len=*), parameter :: Name='user_set_boundary_cells'
    !--------------------------------------------------------------------------
    IsBoundaryCell_GI(:,:,:,ExtraBc_) = r_Blk(:,:,:,iBLK) < rBody

    if(SaveBoundaryCells) return
    call stop_mpi('Set SaveBoundaryCells=.true. in PARAM.in file')

  end subroutine user_set_boundary_cells

end module ModUser

