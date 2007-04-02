module ModUser

  use ModUserEmpty,               &
       IMPLEMENTED1 => user_set_ics

  include 'user_module.h' !list of public methods
 
  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'Multi-ion fluid, G. Toth and Y. Ma'

contains

  !=====================================================================
  subroutine user_set_ics
    use ModMain,     ONLY: globalBLK
    use ModGeometry, ONLY: r_BLK
    use ModVarIndexes
    use ModAdvance,  ONLY: State_VGB
    use ModPhysics,  ONLY: rBody, BodyRho_I, BodyP_I

    implicit none

    integer :: iBlock
    !--------------------------------------------------------------------------
    iBlock = globalBLK

    where(r_BLK(:,:,:,iBlock) < 3*rBody)
       State_VGB(Rho_  ,:,:,:,iBlock) = BodyRho_I(1)
       State_VGB(OpRho_,:,:,:,iBlock) = BodyRho_I(2)
       State_VGB(RhoUx_,:,:,:,iBlock) = 0.0
       State_VGB(RhoUy_,:,:,:,iBlock) = 0.0
       State_VGB(RhoUz_,:,:,:,iBlock) = 0.0
       State_VGB(OpRhoUx_,:,:,:,iBlock) = 0.0
       State_VGB(OpRhoUy_,:,:,:,iBlock) = 0.0
       State_VGB(OpRhoUz_,:,:,:,iBlock) = 0.0
       State_VGB(P_,      :,:,:,iBlock) = BodyP_I(1)
       State_VGB(OpP_,    :,:,:,iBlock) = BodyP_I(2)
    end where

  end subroutine user_set_ics

end module ModUser

