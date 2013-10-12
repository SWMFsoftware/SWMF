!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
! Note this code is adapted from GITM's vertical solver.
! Some variables are different from the rest of the polar
! wind code. I will try to list how they line up.
! dAlt = drbnd, nAlts = nDim, Altitude=ALTD
! These variables were read from other parts of the program before
! now they must be passed in. 




!\
! ------------------------------------------------------------
! calc_rusanov
! ------------------------------------------------------------
!/

subroutine calc_rusanov(nDim,Var_C, VarBottom, VarTop, cMax_C, DiffVar_C)
 
  use ModParameters
  implicit none
  




  integer, intent(in)    :: nDim
  real, intent(in)    :: Var_C(MaxGrid), VarBottom, VarTop, cMax_C(MaxGrid)
  real, intent(out)   :: DiffVar_C(MaxGrid)

  real, dimension(-1:MaxGrid+2):: Var_G, cMax_G
  real, dimension(1:MaxGrid+1) :: VarLeft_F, VarRight_F, DiffFlux_F
  !------------------------------------------------------------
  
  ! Add true ghost cells
  Var_G(1:nDim)         = Var_C(1:nDim)
  Var_G(-1:0)           = VarBottom
  Var_G(nDim+1: nDim+2) = VarTop

  cMax_G(1:nDim) = cMax_C(1:nDim)
  cMax_G(0)      = cMax_C(1)
  cMax_G(nDim+1) = cMax_C(nDim)

  call calc_facevalues(nDim,Var_G, VarLeft_F, VarRight_F)

  ! Rusanov/Lax-Friedrichs diffusive term
  
  !DiffFlux_F(1:nDim+1) = 0.5 * max(cMax_G(0:nDim),cMax_G(1:nDim+1)) &
  !     * (VarRight_F(1:nDim+1) - VarLeft_F(1:nDim+1))
 
  !Global timestep version
  DiffFlux_F(1:nDim+1) = 0.5 * (VarRight_F(1:nDim+1) - VarLeft_F(1:nDim+1))

  DiffVar_C(1:nDim) = DiffFlux_F(2:nDim+1) - DiffFlux_F(1:nDim)

end subroutine calc_rusanov


!\
! ------------------------------------------------------------
! calc_facevalues
! ------------------------------------------------------------
!/

subroutine calc_facevalues(nDim, Var_G, VarLeft_F, VarRight_F)

  use ModParameters
  use ModPWOM, Only: Beta
  implicit none

  integer, intent(in)    :: nDim
  real,    intent(in) :: Var_G(-1:nDim+2)
  real,    intent(out):: VarLeft_F(1:nDim+1), VarRight_F(1:nDim+1)

  real :: dVarUp, dVarDown, dVarLimited_G(0:MaxGrid+1)

  integer :: i
  !----------------------------------------------------------------------------

  do i=0,nDim+1

     dVarUp            = (Var_G(i+1) - Var_G(i))  
     dVarDown          = (Var_G(i)   - Var_G(i-1))
     dVarLimited_G(i) = Limiter_mc(dVarUp, dVarDown)

  end do

  do i=1,nDim+1

     VarLeft_F(i)  = Var_G(i-1) + 0.5*dVarLimited_G(i-1) 
     VarRight_F(i) = Var_G(i)   - 0.5*dVarLimited_G(i) 

  end do

contains
  !---------------------------------------------------
  real function limiter_mc(dUp, dDown)

    real :: dUp, dDown

!    limiter_mc = 0.0
!    RETURN

    if (dUp > 0.0) then
       if (dDown > 0.0) then
          limiter_mc = min(Beta*dUp,Beta*dDown,(dUp+dDown)*0.5)
       else
          limiter_mc = 0.0
       endif
    else
       if (dDown < 0.0) then
          limiter_mc = max(Beta*dUp,Beta*dDown,(dUp+dDown)*0.5)
       else
          limiter_mc = 0.0
       endif
    endif

  end function limiter_mc

end subroutine calc_facevalues
