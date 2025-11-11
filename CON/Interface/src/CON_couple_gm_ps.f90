!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE GM
!^CMP FILE PS

!
! Couple GM and PS components both ways.
!
module CON_couple_gm_ps

  use CON_coupler
  use CON_transfer_data, ONLY: transfer_real_array

  use GM_wrapper, ONLY: GM_put_from_ps
  use PS_wrapper, ONLY: PS_get_for_gm

  implicit none

  private ! except

  public :: couple_gm_ps_init ! initialize both couplings
  public :: couple_gm_ps      ! couple GM to PS
  public :: couple_ps_gm      ! couple PS to GM

  ! revision history:
  ! 11/17/2018 D.Welling <dwelling@uta.edu> - initial version

  logical :: IsInitialized = .false.

  ! Size of the 2D cylindrical structured PS grid
  integer, save :: iSize, jSize

  ! Is GM in multifluid mode?
  logical, save :: DoMultiFluidPSCoupling

  ! Number of fluids in GM (useful as multifluid can have more than 2 fluids)
  integer, save :: nDensityGM

contains
  !============================================================================

  subroutine couple_gm_ps_init
    !DESCRIPTOIN:
    ! Store PS grid size, determine number of GM fluids.
    use ModProcessVarName, ONLY: process_var_name

    integer :: nSpeedGm, nPGm, nPparGm, nWaveGm, nMaterialGm, nChargeStateAllGm

    ! Set initialization status:
    !--------------------------------------------------------------------------
    if(IsInitialized) RETURN
    IsInitialized = .true.

    if(.not. (is_proc(PS_) .or. is_proc(GM_))) RETURN

    ! This works for a regular PS grid only
    iSize = Grid_C(PS_) % nCoord_D(1)
    jSize = Grid_C(PS_) % nCoord_D(2)

    ! Get information on number of fluids used in GM:
    call process_var_name(Grid_C(GM_)%NameVar, nDensityGm, nSpeedGm, &
         nPGm, nPparGm, nWaveGm, nMaterialGm, nChargeStateAllGm)

    DoMultiFluidPSCoupling = nDensityGm > 1
    ! write(*,*)'!!!! DEBUG: DoMultiFluidPSCoupling = ', DoMultiFluidPSCoupling

  end subroutine couple_gm_ps_init
  !============================================================================

  subroutine couple_gm_ps(tSimulation)

    use CON_world, ONLY: get_comp_info
    use CON_comp_param, ONLY: lNameVersion

    real, intent(in) :: tSimulation     ! simulation time at coupling

    ! Couple between two components:
    !    Global Magnetosphere (GM) source
    !      (PS) target
    !
    ! Send field line volumes & field magnitudes to PS from GM.

    character(len=*), parameter:: NameSub = 'couple_gm_ps'
    !--------------------------------------------------------------------------

    call CON_stop(NameSub//': GM to PS coupling not implemented yet.')

  end subroutine couple_gm_ps
  !============================================================================

  subroutine couple_ps_gm(tSimulation)

    real, intent(in) :: tSimulation     ! simulation time at coupling

    ! Couple between two components:
    !    Plasmasphere  (PS) source
    !    Global Magnetosphere (GM) target
    !
    ! Send pressure from PS to GM.

    ! Number of variables to pass
    integer:: nVarPsGm

    ! Names of variables to pass
    character(len=100) :: NameVar

    ! Buffer for the variables on the 2D PS grid
    real, allocatable:: Buffer_IIV(:,:,:)

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'couple_ps_gm'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    if(DoTest)write(*,*)NameSub,' starting, iProc=',i_proc()

    if(DoMultiFluidPSCoupling)then
       NameVar='p:rho:Hpp:Opp:Hprho:Oprho'
       nVarPsGm=6
    else
       NameVar='p:rho'
       nVarPsGm=2
    end if

    allocate(Buffer_IIV(iSize,jSize,nVarPsGm))
    if(is_proc(PS_)) &
         call PS_get_for_gm(Buffer_IIV, iSize, jSize, nVarPsGm, NameVar)
    call transfer_real_array(PS_, GM_, size(Buffer_IIV), Buffer_IIV)
    if(is_proc(GM_)) &
         call GM_put_from_ps(Buffer_IIV, iSize, jSize, nVarPsGm, NameVar)
    deallocate(Buffer_IIV)

    if(DoTest)write(*,*)NameSub,': finished iProc=', i_proc()

  end subroutine couple_ps_gm
  !============================================================================

end module CON_couple_gm_ps
!==============================================================================
