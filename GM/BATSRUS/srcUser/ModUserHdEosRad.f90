!^CFG COPYRIGHT UM
!============================================================================
module ModUser
  ! This is the default user module which contains empty methods defined
  ! in ModUserEmpty.f90

  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_update_states,              &
       IMPLEMENTED2 => user_calc_sources

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'HYDRO + Grey Diffusion'
contains
!============================================================================
  subroutine user_update_states(iStage,iBlock)
    use ModVarIndexes
    use ModSize
    use ModAdvance, ONLY: State_VGB
    use ModMain,    ONLY: nStage
    use ModPhysics
    use ModEnergy,  ONLY: calc_energy_cell
    use ModEos
    implicit none

    integer,intent(in):: iStage,iBlock
    integer:: i,j,k
    real:: PressureSI,EInternal,EInternalSI,RhoSI

    !------------------------------------------------------------------------
    call update_states_MHD(iStage,iBlock)
    !\
    ! Begin update of pressure and relaxation energy::
    !/
    !  if (iStage/=nStage) return
    do k=1,nK; do j=1,nJ; do i=1,nI
       !Total plasma internal energy, P_gas/(\gamma -1) + ExtraEInt
       ! transformed to SI
       EInternalSI = No2Si_V(UnitEnergyDens_) *( inv_gm1 &
            *(State_VGB(p_,i,j,k,iBlock)-State_VGB(ERad_,i,j,k,iBlock)/3.0) &
            + State_VGB(ExtraEInt_,i,j,k,iBlock) )

       !Density, transformed to SI
       RhoSI = No2Si_V(UnitRho_)*State_VGB(Rho_,i,j,k,iBlock)
       
       
       !Apply the EOS, get gas kinetic pressure in SI
       call eos(&
            UDensityTotal=EInternalSI,& !Input total energy density SI,[J/m^3]
            Rho=RhoSI,                & !Input mass density, SI [kg/m^3] 
            iMaterial=0,              & !Input: sort of material
            PTotalOut=PressureSI      ) !Output, OPTIONAL, pressure, SI [Pa]
  
       !Put total pressure (gas kinetic + radiation) and
       ! ExtraEInt = Total plasma internal energy - P_gas/(\gamma -1)
       State_VGB(p_,i,j,k,iBlock) = PressureSI*Si2No_V(UnitP_) &
            + State_VGB(ERad_,i,j,k,iBlock)/3.0
       State_VGB(ExtraEInt_,i,j,k,iBlock) = Si2No_V(UnitEnergyDens_)*&
            (EInternalSI - PressureSI*inv_gm1)
    end do; end do; end do
    call calc_energy_cell(iBlock)
    !\
    ! End update of pressure and relaxation energy::
    !/
  end subroutine user_update_states

  !==========================================================================

  subroutine user_calc_sources
    use ModAdvance,    ONLY: State_VGB, Source_VC, &
         uDotArea_XI, uDotArea_YI, uDotArea_ZI
    use ModGeometry,   ONLY: vInv_CB
    use ModMain,       ONLY: GlobalBLK
    use ModSize,       ONLY: nI, nJ, nK
    use ModVarIndexes, ONLY: ERad_, p_, Energy_
    implicit none

    integer :: iBlock, i, j, k, iFluid
    real :: PRadDivU

    character (len=*), parameter :: NameSub = 'user_calc_sources'
    !------------------------------------------------------------------------

    iBlock = GlobalBlk

    iFluid = 1

    ! add -Prad*Div(U) to radiation energy density
    do k=1,nK; do j=1,nJ; do i=1,nI

       PRadDivU = State_VGB(ERad_,i,j,k,iBlock)/3.0 &
            *vInv_CB(i,j,k,iBlock)&
            *(uDotArea_XI(i+1,j,k,iFluid) - uDotArea_XI(i,j,k,iFluid) &
             +uDotArea_YI(i,j+1,k,iFluid) - uDotArea_YI(i,j,k,iFluid) &
             +uDotArea_ZI(i,j,k+1,iFluid) - uDotArea_ZI(i,j,k,iFluid))

       ! add -Prad*Div(U) to radiation energy density
       Source_VC(ERad_,i,j,k) = Source_VC(ERad_,i,j,k) - PRadDivU

       ! add + 1/2 *Prad*Div(U) to the total energy density
       Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + 0.5*PRadDivU

       ! add + 1/3 *Prad*Div(U) to the total pressure
       Source_VC(p_,i,j,k) = Source_VC(p_,i,j,k) + PRadDivU/3.0
    end do; end do; end do

  end subroutine user_calc_sources

end module ModUser
