!^CFG COPYRIGHT UM
!==============================================================================
module ModUser
  use EEE_ModMain
  use ModMagnetogram
  use ModExpansionFactors
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPMENENTED2 => user_set_ics,                    &
       IMPLEMENTED3 => user_initial_perturbation,       &
       IMPLEMENTED4 => user_face_bcs,                   &
       IMPLEMENTED5 => user_get_log_var,                &
       IMPLEMENTED6 => user_get_b0,                     &
       IMPLEMENTED7 => user_update_states,              &
       IMPLEMENTED8 => user_specify_initial_refinement

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'EMPIRICAL SC - Cohen, Sokolov'

contains
  !============================================================================
  subroutine user_read_inputs
    use ModMain
    use ModProcMH,    ONLY: iProc
    use ModReadParam
    use ModIO,        ONLY: write_prefix, write_myname, iUnitOut
    use ModPhysics,   ONLY: BodyNDim_I,BodyTDim_I,g

    integer:: i
    character (len=100) :: NameCommand
    character (len=lStringLine)   :: NameModel
    !-------------------------------------------------------------------------

    call EEE_initialize(BodyNDim_I(1),BodyTDim_I(1),g)

    if(iProc==0.and.lVerbose > 0)then
       call write_prefix; write(iUnitOut,*)'User read_input HELIOSPHERE starts'
    endif
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#PFSSM")
          call read_var('UseUserB0'  ,UseUserB0)
          if (UseUserB0)then
             call read_magnetogram_file
             call read_var('dt_UpdateB0',dt_UpdateB0)
             DoUpdateB0 = dt_updateb0 > 0.0
          endif
       case("#TD99FLUXROPE")
          call read_var('UseTD99Perturbation' ,UseTD99Perturbation)
          call read_var('UseVariedCurrent'    ,UseVariedCurrent)
          call read_var('CurrentStartTime'    ,CurrentStartTime)
          call read_var('CurrentRiseTime '    ,CurrentRiseTime)
          call read_var('DoTD99FluxRope'      ,DoTD99FluxRope)
          call read_var('DoEquilItube'        ,DoEquilItube)
          call read_var('DoRevCurrent'        ,DoRevCurrent)
          call read_var('aratio_TD99'         ,aratio_TD99)
          call read_var('Itube_TD99'          ,Itube_TD99)
          call read_var('Rtube_TD99'          ,Rtube_TD99)
          call read_var('atube_TD99'          ,atube_TD99)
          call read_var('d_TD99'              ,d_TD99)
          call read_var('Mass_TD99'           ,Mass_TD99)
          call read_var('LongitudeTD99'       ,LongitudeTD99)
          call read_var('LatitudeTD99'        ,LatitudeTD99)
          call read_var('OrientationTD99'     ,OrientationTD99)
          call read_var('DoBqField'           ,DoBqField)
          call read_var('q_TD99'              ,q_TD99)
          call read_var('L_TD99'              ,L_TD99)
       case("#GL98FLUXROPE")
          call read_var('UseFluxRope',     UseFluxRope)
          call read_var('cme_a',           cme_a)
          call read_var('cme_r1',          cme_r1)
          call read_var('cme_r0',          cme_r0)
          call read_var('cme_a1',          cme_a1)
          call read_var('cme_alpha',       cme_alpha)
          call read_var('cme_rho1',        cme_rho1)
          call read_var('cme_rho2',        cme_rho2)
          call read_var('ModulationRho',   ModulationRho)
          call read_var('ModulationP',     ModulationP)
          call read_var('LongitudeGL98',   LongitudeGL98)
          call read_var('LatitudeGL98',    LatitudeGL98)
          call read_var('OrientationGL98', OrientationGL98)
       case("#EMPIRICALSW")
          call read_var('NameModel',NameModel)
          call set_empirical_model(trim(NameModel))
       case('#USERINPUTEND')
          if(iProc==0.and.lVerbose > 0)then
             call write_prefix;
             write(iUnitOut,*)'User read_input HELIOSPHERE ends'
          endif
          EXIT
       case default
          if(iProc==0) then
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
  subroutine user_face_bcs(VarsGhostFace_V)

    use ModSize,       ONLY: East_,West_,South_,North_,Bot_,Top_
    use ModMain,       ONLY: time_accurate,x_,y_,z_, UseRotatingFrame, &
         n_step, Iteration_Number
    use ModVarIndexes, ONLY: nVar,Ew_,rho_,Ux_,Uy_,Uz_,Bx_,By_,Bz_,P_

    use ModGeometry,   ONLY: R_BLK
    use ModAdvance,    ONLY: State_VGB
    use ModPhysics,    ONLY: inv_gm1,OmegaBody,Si2No_V,UnitB_,UnitU_,UnitRho_
    use ModNumConst,   ONLY: cTolerance,cTiny

    use ModBlockData, ONLY: use_block_data, put_block_data, get_block_data

    use ModFaceBc, ONLY: FaceCoords_D, VarsTrueFace_V, TimeBc, &
         iFace, jFace, kFace, iSide, iBlockBc

    implicit none

    real, intent(out):: VarsGhostFace_V(nVar)

    integer:: iCell,jCell,kCell

    real:: DensCell,PresCell,GammaCell, B1dotR  
    real, dimension(3):: RFace_D,B1_D,U_D,B1t_D,B1n_D

    real :: RhoCME,UCME_D(nDim),BCME_D(nDim),pCME
    real :: BCMEn, BCMEn_D(nDim)
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

    BCMEn   = dot_product(RFace_D,BCME_D)
    BCMEn_D = BCMEn*RFace_D
    !\
    ! Fix the normal component of the CME field to BCMEn_D at the Sun::
    !/
    VarsGhostFace_V(Bx_:Bz_) = VarsGhostFace_V(Bx_:Bz_) + BCMEn_D

    ! CME velocity update should be added here

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
    VarsGhostFace_V(P_) = &
         max(VarsGhostFace_V(Rho_)*(PresCell+pCME)/(DensCell+RhoCME), &
         VarsTrueFace_V(P_))
    VarsGhostFace_V(Ew_) = &!max(-VarsTrueFace_V(Ew_)+ &
         VarsGhostFace_V(Rho_)*(PresCell+pCME)/(DensCell+RhoCME) &
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
    use ModPhysics,    ONLY: g,inv_g,GBody,BodyTdim_I
    use ModExpansionFactors,  ONLY: UMin,T0
    implicit none

    integer, intent(in)  :: iCell,jCell,kCell,iBlock
    real, intent(out)    :: DensCell,PresCell,GammaCell
    real :: UFinal       !The solar wind speed at the far end of the Parker spiral,
                         !which originates from the given cell
    real :: URatio       !The coronal based values for temperature density 
                         !are scaled as functions of UFinal/UMin ratio
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
    DensCell  = ((1.0/URatio)**2)*&          !This is the density variation
         exp(-GBody*g*&
         (min(URatio,2.0)*BodyTdim_I(1)/T0)*&!This is the temperature variation
         (1.0/max(R_BLK(iCell,jCell,kCell,iBlock),0.90)&
         -1.0))

    PresCell  = inv_g*DensCell*&
         T0/(min(URatio,2.0)*BodyTdim_I(1))  !This is the temperature variation
  end subroutine get_plasma_parameters_cell

  !============================================================================

  subroutine user_initial_perturbation
    use ModMain, ONLY: nI,nJ,nK,nBLK,unusedBLK,x_,y_,z_,n_step,iteration_number
    use ModVarIndexes
    use ModAdvance,   ONLY: State_VGB 
    use ModPhysics,   ONLY: Si2No_V,UnitU_,UnitRho_,UnitP_,UnitB_
    use ModGeometry
    use ModEnergy,    ONLY: calc_energy_cell
    implicit none

    integer :: i,j,k,iBLK
    logical :: oktest,oktest_me
    real :: x_D(nDim),Rho,U_D(nDim),B_D(nDim),p

    real :: Mass=0.0
    !--------------------------------------------------------------------------
    call set_oktest('user_initial_perturbation',oktest,oktest_me)

    do iBLK=1,nBLK
       if(unusedBLK(iBLK))CYCLE
       do k=1,nK; do j=1,nJ; do i=1,nI

          x_D(x_) = x_BLK(i,j,k,iBLK)
          x_D(y_) = y_BLK(i,j,k,iBLK)
          x_D(z_) = z_BLK(i,j,k,iBLK)

          call EEE_get_state_init(x_D,Rho,U_D,B_D,p, &
               n_step,iteration_number)

          Rho = Rho*Si2No_V(UnitRho_)
          U_D = U_D*Si2No_V(UnitU_)
          B_D = B_D*Si2No_V(UnitB_)
          p = p*Si2No_V(UnitP_)

          !\
          ! Add the eruptive event state to the solar wind
          !/
          State_VGB(Rho_,i,j,k,iBLK) = State_VGB(Rho_,i,j,k,iBLK) + Rho
          State_VGB(RhoUx_,i,j,k,iBLK) = &
               State_VGB(RhoUx_,i,j,k,iBLK) + Rho*U_D(x_)
          State_VGB(RhoUy_,i,j,k,iBLK) = &
               State_VGB(RhoUy_,i,j,k,iBLK) + Rho*U_D(y_)
          State_VGB(RhoUz_,i,j,k,iBLK) = &
               State_VGB(RhoUz_,i,j,k,iBLK) + Rho*U_D(z_)
          State_VGB(Bx_,i,j,k,iBLK) = State_VGB(Bx_,i,j,k,iBLK) + B_D(x_)
          State_VGB(By_,i,j,k,iBLK) = State_VGB(By_,i,j,k,iBLK) + B_D(y_)
          State_VGB(Bz_,i,j,k,iBLK) = State_VGB(Bz_,i,j,k,iBLK) + B_D(z_)
          State_VGB(P_,i,j,k,iBLK) = State_VGB(P_,i,j,k,iBLK) + p

          !\
          ! Calculate the mass added to the eruptive event
          !/
          Mass = Mass + Rho/vInv_CB(i,j,k,iBLK)
       end do; end do; end do

       !\
       ! Update the total energy::
       !/
       call calc_energy_cell(iBLK)

    end do

  end subroutine user_initial_perturbation

  !============================================================================

  subroutine user_set_ics
    use ModMain,      ONLY: globalBLK,nI,nJ,nK
    use ModVarIndexes
    use ModAdvance,   ONLY: State_VGB 
    use ModPhysics,   ONLY: inv_gm1
    use ModGeometry
    use ModEnergy,    ONLY: calc_energy_cell
    implicit none

    integer :: i,j,k,iBLK
    logical :: oktest,oktest_me
    real :: Dens_BLK,Pres_BLK,Gamma_BLK
    real :: x,y,z,R,ROne,Rmax
    !--------------------------------------------------------------------------
    call set_oktest('user_set_ics',oktest,oktest_me)

    iBLK = globalBLK

    do k=1,nK; do j=1,nJ; do i=1,nI
       x = x_BLK(i,j,k,iBLK)
       y = y_BLK(i,j,k,iBLK)
       z = z_BLK(i,j,k,iBLK)
       R = sqrt(x**2+y**2+z**2+cTolerance**2)
       ROne  = max(1.0,R)
       Rmax  = max(2.1E+01,sqrt(x2**2+y2**2+z2**2))
       State_VGB(Bx_,i,j,k,iBLK) = 0.0
       State_VGB(By_,i,j,k,iBLK) = 0.0
       State_VGB(Bz_,i,j,k,iBLK) = 0.0
       call get_plasma_parameters_cell(i,j,k,iBLK,&
            Dens_BLK,Pres_BLK,Gamma_BLK)
       State_VGB(rho_,i,j,k,iBLK) = Dens_BLK
       State_VGB(P_,i,j,k,iBLK)   = Pres_BLK
       State_VGB(RhoUx_,i,j,k,iBLK) = Dens_BLK &
            *4.0*((ROne-1.0)/(Rmax-1.0))*x/R
       State_VGB(RhoUy_,i,j,k,iBLK) = Dens_BLK &
            *4.0*((ROne-1.0)/(Rmax-1.0))*y/R
       State_VGB(RhoUz_,i,j,k,iBLK) = Dens_BLK &
            *4.0*((ROne-1.0)/(Rmax-1.0))*z/R
       State_VGB(Ew_,i,j,k,iBLK) = Pres_BLK &
            *(1.0/(Gamma_BLK-1.0)-inv_gm1) 
    end do; end do; end do

  end subroutine user_set_ics

  !============================================================================
  subroutine user_get_b0(xInput,yInput,zInput,B0_D)
    use ModPhysics,  ONLY: Io2No_V,UnitB_
    implicit none
    real, intent(in):: xInput,yInput,zInput
    real, intent(out), dimension(3):: B0_D
    call get_magnetogram_field(xInput,yInput,zInput,B0_D)
    B0_D = B0_D*Io2No_V(UnitB_)
  end subroutine user_get_b0

  !============================================================================

  subroutine user_update_states(iStage,iBlock)
    use ModVarIndexes
    use ModSize
    use ModAdvance, ONLY: State_VGB, B0xCell_BLK, B0yCell_BLK, B0zCell_BLK
    use ModMain,    ONLY: nStage
    use ModPhysics, ONLY: inv_gm1
    use ModGeometry,ONLY: R_BLK
    use ModEnergy,  ONLY: calc_energy_cell
    use ModExpansionFactors, ONLY: gammaSS
    implicit none
    integer,intent(in):: iStage,iBlock
    integer:: i,j,k
    real:: DensCell,PresCell,GammaCell,Beta
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
            (State_VGB(P_   ,i,j,k,iBlock)+(&
            (State_VGB(Bx_   ,i,j,k,iBlock)+B0xCell_BLK(i,j,k,iBlock))**2+&
            (State_VGB(By_   ,i,j,k,iBlock)+B0yCell_BLK(i,j,k,iBlock))**2+&
            (State_VGB(Bz_   ,i,j,k,iBlock)+B0zCell_BLK(i,j,k,iBlock))**2)&
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

  !========================================================================
  subroutine user_get_log_var(VarValue,TypeVar,Radius)

    use ModProcMH,     ONLY: nProc
    use ModIO,         ONLY: dn_output,logfile_,write_myname
    use ModMain,       ONLY: unusedBLK,nBLK,iteration_number,   &
         x_,y_,z_
    use ModVarIndexes, ONLY: Ew_,Bx_,By_,Bz_,rho_,rhoUx_,rhoUy_,rhoUz_,P_ 
    use ModGeometry,   ONLY: R_BLK
    use ModAdvance,    ONLY: State_VGB,tmp1_BLK,B0xCell_BLK,    &
         B0yCell_BLK,B0zCell_BLK
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
               (B0xcell_BLK(:,:,:,iBLK)+State_VGB(Bx_,:,:,:,iBLK))**2+&
               (B0ycell_BLK(:,:,:,iBLK)+State_VGB(By_,:,:,:,iBLK))**2+&
               (B0zcell_BLK(:,:,:,iBLK)+State_VGB(Bz_,:,:,:,iBLK))**2
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
  !----------------------------------------------------------
  subroutine user_specify_initial_refinement(iBLK,refineBlock,lev,DxBlock, &
       xCenter,yCenter,zCenter,rCenter,                        &
       minx,miny,minz,minR,maxx,maxy,maxz,maxR,found)
    use ModMain,ONLY:time_loop,nI,nJ,nK
    use ModAMR,ONLY:InitialRefineType
    use ModNumConst
    use ModAdvance,ONLY:&
         State_VGB,Bx_,By_,Bz_,B0xCell_BLK,B0yCell_BLK,B0zCell_BLK
    use ModGeometry
    use ModPhysics,ONLY:rBody
    implicit none
    logical,intent(out) :: refineBlock, found
    integer, intent(in) :: lev
    real, intent(in)    :: DxBlock
    real, intent(in)    :: xCenter,yCenter,zCenter,rCenter
    real, intent(in)    :: minx,miny,minz,minR
    real, intent(in)    :: maxx,maxy,maxz,maxR
    integer, intent(in) :: iBLK

    character (len=*), parameter :: Name='user_specify_initial_refinement'
    real::BDotRMin,BDotRMax,critx
    integer::i,j,k
    !-------------------------------------------------------------------
    select case (InitialRefineType)
    case ('helio_init')
       if(.not.time_loop)then
          !refine to have resolution not worse 4.0 and
          !refine the body intersecting blocks
          ! Refine all blocks time through (starting with 1 block)
          if (lev <= 4) then
             refineBlock = .true.
          else
             critx=(XyzMax_D(1)-XyzMin_D(1))/(2.0**real(lev-2))
             if ( rCenter < 1.10*rBody + critx ) then
                refineBlock = .true.
             else
                refineBlock = .false.
             end if
          endif
       elseif(dx_BLK(iBLK)<0.20.or.far_field_BCs_BLK(iBLK))then
          refineBlock=.false. !Do not refine body or outer boundary
       else
          !refine heliosheath
          BDotRMin=0.0
          do k=0,nK+1;do j=1,nJ
             BDotRMin=min( BDotRMin,minval(&
                  (B0xCell_BLK(1:nI,j,k,iBLK)+&
                  State_VGB(Bx_,1:nI,j,k,iBLK))*&
                  x_BLK(1:nI,j,k,iBLK)+&
                  (B0yCell_BLK(1:nI,j,k,iBLK)+&
                  State_VGB(By_,1:nI,j,k,iBLK))*&
                  y_BLK(1:nI,j,k,iBLK)+&
                  (B0zCell_BLK(1:nI,j,k,iBLK)+&
                  State_VGB(Bz_,1:nI,j,k,iBLK))*&
                  z_BLK(1:nI,j,k,iBLK)))
          end do;end do
          BDotRMax=0.0
          do k=0,nK+1;do j=1,nJ
             BDotRMax=max( BDotRMax,maxval(&
                  (B0xCell_BLK(1:nI,j,k,iBLK)+&
                  State_VGB(Bx_,1:nI,j,k,iBLK))*&
                  x_BLK(1:nI,j,k,iBLK)+&
                  (B0yCell_BLK(1:nI,j,k,iBLK)+&
                  State_VGB(By_,1:nI,j,k,iBLK))*&
                  y_BLK(1:nI,j,k,iBLK)+&
                  (B0zCell_BLK(1:nI,j,k,iBLK)+&
                  State_VGB(Bz_,1:nI,j,k,iBLK))*&
                  z_BLK(1:nI,j,k,iBLK)))
          end do;end do
          refineBlock =BDotRMin<-cTiny.and.&
               BDotRMax>cTiny
       end if
       found=.true.
    end select
  end subroutine user_specify_initial_refinement

end module ModUser

