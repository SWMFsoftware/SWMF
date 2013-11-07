!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine GM_get_for_pt(IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, &
     Data_VI)

  ! Get magnetic field data from GM to PT

  use ModProcMH,  ONLY: iProc
  use ModPhysics, ONLY: Si2No_V, UnitX_, No2Si_V, UnitB_
  use ModAdvance, ONLY: State_VGB, Bx_, Bz_
  use ModB0,      ONLY: UseB0, get_b0
  use BATL_lib,   ONLY: nDim, MaxDim, MaxJ, MaxK, find_grid_block

  implicit none

  logical,          intent(in):: IsNew   ! true for new point array
  character(len=*), intent(in):: NameVar ! List of variables
  integer,          intent(in):: nVarIn  ! Number of variables in Data_VI
  integer,          intent(in):: nDimIn  ! Dimensionality of positions
  integer,          intent(in):: nPoint  ! Number of points in Xyz_DI

  real, intent(in) :: Xyz_DI(nDimIn,nPoint)  ! Position vectors
  real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array

  real:: Xyz_D(MaxDim), b_D(MaxDim)
  real:: Dist_D(MaxDim), Dx1, Dx2, Dy1, Dy2, Dz1, Dz2
  integer:: iCell_D(MaxDim), i1, i2, j1, j2, k1, k2

  integer:: iPoint, iBlock, iProcFound

  character(len=*), parameter :: NameSub='GM_get_for_pt'
  !--------------------------------------------------------------------------
  if(nVarIn /= 3 .or. NameVar /= 'Bx By Bz')then
     write(*,*)NameSub,' ERROR: nVarIn, NameVar=', nVarIn, NameVar
     call stop_mpi(NameSub//' is only implemented for transferring B field')
  end if

  ! We should have second order accurate magnetic field in the ghost cells

  do iPoint = 1, nPoint

     Xyz_D = Xyz_DI(:,iPoint)*Si2No_V(UnitX_)
     call find_grid_block(Xyz_D, iProcFound, iBlock, iCell_D, Dist_D)

     if(iProcFound /= iProc)then
        write(*,*)NameSub,' ERROR: Xyz_D, iProcFound=', Xyz_D, iProcFound
        call stop_mpi(NameSub//' could not find position on this proc')
     end if

     if(UseB0)then
        call get_b0(Xyz_D, b_D)
     else
        b_D = 0.0
     end if

     if(Dist_D(1)>0) then
        i1=iCell_D(1)
        i2=i1+1
        Dx1=Dist_D(1)
        Dx2=1.0-Dx1
     else
        i1=iCell_D(1)-1
        i2=i1+1
        Dx1=1.0+Dist_D(1)
        Dx2=1.0-Dx1
     end if

     if(Dist_D(2)>0) then
        j1=iCell_D(2)
        j2=min(j1+1, MaxJ)
        Dy1=Dist_D(2)
        Dy2=1.0-Dy1
     else
        j1=iCell_D(2)-1
        j2=min(j1+1, MaxJ)
        Dy1=1.0+Dist_D(2)
        Dy2=1.0-Dy1
     end if

     if(Dist_D(3)>0) then
        k1=iCell_D(3)
        k2=min(k1+1, MaxK)
        Dz1=Dist_D(3)
        Dz2=1.0-Dz1
     else
        k1=iCell_D(3)-1
        k2=min(k1+1, MaxK)
        Dz1=1.0+Dist_D(3)
        Dz2=1.0-Dz1
     end if

     b_D = b_D + Dz2*(Dy2*(Dx2*State_VGB(Bx_:Bz_,i1,j1,k1,iBlock)   &
          +                Dx1*State_VGB(Bx_:Bz_,i2,j1,k1,iBlock))  &
          +           Dy1*(Dx2*State_VGB(Bx_:Bz_,i1,j2,k1,iBlock)   &
          +                Dx1*State_VGB(Bx_:Bz_,i2,j2,k1,iBlock))) &
          +      Dz1*(Dy2*(Dx2*State_VGB(Bx_:Bz_,i1,j1,k2,iBlock)   &
          +                Dx1*State_VGB(Bx_:Bz_,i2,j1,k2,iBlock))  &
          +           Dy1*(Dx2*State_VGB(Bx_:Bz_,i1,j2,k2,iBlock)   &
          +                Dx1*State_VGB(Bx_:Bz_,i2,j2,k2,iBlock)))

     Data_VI(:,iPoint) = b_D*No2Si_V(UnitB_)

  end do

end subroutine GM_get_for_pt
