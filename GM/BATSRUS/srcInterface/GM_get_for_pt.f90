!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine GM_get_for_pt(IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, &
     Data_VI)

  ! Get magnetic field data from GM to PT

  use ModProcMH,  ONLY: iProc
  use ModPhysics, ONLY: Si2No_V, UnitX_, No2Si_V, UnitB_, UnitU_, UnitP_, UnitRho_, iUnitCons_V
  use ModAdvance, ONLY: State_VGB, Rho_, RhoUx_, RhoUz_, Bx_, Bz_, nVar, p_, NamePrimitiveVar
  use ModVarIndexes, ONLY: nVar, NamePrimitiveVar
  use ModB0,      ONLY: UseB0, get_b0
  use BATL_lib,   ONLY: nDim, MaxDim, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, find_grid_block
  use ModInterpolate, ONLY: trilinear

  implicit none

  logical,          intent(in):: IsNew   ! true for new point array
  character(len=*), intent(in):: NameVar ! List of variables
  integer,          intent(in):: nVarIn  ! Number of variables in Data_VI
  integer,          intent(in):: nDimIn  ! Dimensionality of positions
  integer,          intent(in):: nPoint  ! Number of points in Xyz_DI

  real, intent(in) :: Xyz_DI(nDimIn,nPoint)  ! Position vectors
  real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array

  real:: Xyz_D(MaxDim), b_D(MaxDim), u_D(MaxDim), rho_D, p_D
  real:: Dist_D(MaxDim), Dx1, Dx2, Dy1, Dy2, Dz1, Dz2, State_V(nVar)
  integer:: iCell_D(MaxDim), i1, i2, j1, j2, k1, k2

  integer:: iPoint, iBlock, iProcFound

  character(len=*), parameter :: NameSub='GM_get_for_pt'
  !--------------------------------------------------------------------------
  if(nVarIn /= nVar .or. NameVar /= NamePrimitiveVar )then
     write(*,*)NameSub,' ERROR: nVarIn, NameVar=', nVarIn, NameVar
     call stop_mpi(NameSub//': Requested variables are inconsistent with the current equation set.')
  end if

  ! We should have second order accurate magnetic field in the ghost cells

  do iPoint = 1, nPoint

     Xyz_D = Xyz_DI(:,iPoint)*Si2No_V(UnitX_)
     call find_grid_block(Xyz_D, iProcFound, iBlock, iCell_D, Dist_D, &
          UseGhostCell = .true.)

     if(iProcFound /= iProc)then
        write(*,*)NameSub,' ERROR: Xyz_D, iProcFound=', Xyz_D, iProcFound
        call stop_mpi(NameSub//' could not find position on this proc')
     end if

     State_V = trilinear(State_VGB(:,:,:,:,iBlock), &
         nVar, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, Xyz_D, iCell_D = iCell_D, Dist_D = Dist_D)

     if(UseB0)then
        call get_b0(Xyz_D, b_D)
        State_V(Bx_:Bz_) = State_V(Bx_:Bz_) + b_D
     else
        b_D = 0.0
     end if

     Data_VI(:,iPoint) = State_V*No2Si_V(iUnitCons_V)

  end do

end subroutine GM_get_for_pt
