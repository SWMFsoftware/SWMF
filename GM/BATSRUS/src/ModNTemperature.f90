!^CFG COPYRIGHT UM
!============================================================================
module ModNTemperature

  implicit none
  save

  private !except

  ! Back-end module for implicitly updating the N temperature
  ! diffusion-relaxation model.

  integer, public :: nTemperature
  integer, allocatable :: iTemperature_I(:), iRelaxationT_I(:)
  real, allocatable :: &
       Temperature_VGB(:,:,:,:,:), &
       HeatConductionCoef_VG(:,:,:,:), &
       SpecificHeat_VC(:,:,:,:), &
       RelaxationCoef_VC(:,:,:,:)

contains

  !==========================================================================
  subroutine init_NTemperature

    use ModSize, ONLY: nI, nJ, nK

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = 'init_NTemperature'
    !------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    if(nTemperature<1) &
       call stop_mpi('ERROR in init_NTemperature: nTemperature must be > 0')

    if(.not.allocated(Temperature_VGB)) then
       allocate( &
            iTemperature_I(nTemperature), &
            Temperature_VGB(nTemperature,0:nI+1,0:nJ+1,0:nK+1,nBlk), &
            HeatConductionCoef_VG(nTemperature,0:nI+1,0:nJ+1,0:nK+1), &
            SpecificHeat_VC(nTemperature,1:nI,1:nJ,1:nK) )

       if(nTemperature>1)then
          allocate( &
               RelaxationCoef_VC(2:nTemperature,1:nI,1:nJ,1:nK), &
               iRelaxationT_I(2:nTemperature) )
       end if
    end if

  end subroutine init_NTemperature

end module ModNTemperature
