!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module GM_couple_pt

  implicit none

  private ! except
  public:: GM_get_for_pt

contains

  subroutine GM_get_for_pt(IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, &
       Data_VI)

    ! Get magnetic field data from GM to PT

    use ModProcMH,  ONLY: iProc
    use ModPhysics, ONLY: Si2No_V, No2Si_V, iUnitCons_V, UnitX_
    use ModAdvance, ONLY: State_VGB, nVar, Bx_, Bz_
    use ModVarIndexes, ONLY: nVar
    use ModB0,      ONLY: UseB0, get_b0
    use BATL_lib,   ONLY: MaxDim, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, find_grid_block
    use ModInterpolate, ONLY: trilinear

    logical,          intent(in):: IsNew   ! true for new point array
    character(len=*), intent(in):: NameVar ! List of variables
    integer,          intent(in):: nVarIn  ! Number of variables in Data_VI
    integer,          intent(in):: nDimIn  ! Dimensionality of positions
    integer,          intent(in):: nPoint  ! Number of points in Xyz_DI

    real, intent(in) :: Xyz_DI(nDimIn,nPoint)  ! Position vectors
    real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array

    real:: Xyz_D(MaxDim), b_D(MaxDim)
    real:: Dist_D(MaxDim), State_V(nVar)
    integer:: iCell_D(MaxDim)

    integer:: iPoint, iBlock, iProcFound

    character(len=*), parameter :: NameSub='GM_get_for_pt'
    !--------------------------------------------------------------------------

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
            nVar, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, Xyz_D, &
            iCell_D = iCell_D, Dist_D = Dist_D)

       if(UseB0)then
          call get_b0(Xyz_D, b_D)
          State_V(Bx_:Bz_) = State_V(Bx_:Bz_) + b_D
       else
          b_D = 0.0
       end if

       Data_VI(:,iPoint) = State_V*No2Si_V(iUnitCons_V)

    end do

  end subroutine GM_get_for_pt

end module GM_couple_pt
