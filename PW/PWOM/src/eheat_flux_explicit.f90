!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine PW_eheat_flux_explicit(nCell,&
     Rgas, DtIn,&
     OldState_GV,&
     RhoSource_C, RhoUSource_C, eSource_C,&
     HeatCon_G,NewT_G)

  use ModCommonVariables,ONLY: Ar12,Ar23,CellVolume_C,Gmin1,DrBnd,Gamma
  
  implicit none

  integer, intent(in)      :: nCell
  real, intent(in)         :: Rgas,DtIn
  real, intent(in)         :: OldState_GV(-1:nCell+2,3)
  real, dimension(nCell), intent(in)  :: RhoSource_C, RhoUSource_C, eSource_C
  real, intent(in)  :: HeatCon_G(0:nCell+1)
  real, intent(out) :: NewT_G(-1:nCell+2)
  
  
  integer,parameter    :: Rho_=1,U_=2,P_=3
  real, allocatable    :: T_G(:),Conduction_G(:),Diffusion_C(:),Conduction_F(:)
  real, allocatable    :: DivU_C(:),LeftU_F(:),RightU_F(:),&
                          LeftT_F(:),RightT_F(:),u_F(:),tFlux_F(:),&
                          DivUT_C(:),GradT_F(:)
  integer :: i,iCell
  real :: Coeff
  !---------------------------------------------------------------------------
  if (.not.allocated(T_G)) then
     allocate(T_G(-1:nCell+2),Conduction_G(-1:nCell+1),Conduction_F(0:nCell), &
              GradT_F(0:nCell+1),Diffusion_C(nCell),      &
              DivU_C(1:nCell),DivUT_C(1:nCell), &
              LeftU_F(1:nCell+1),RightU_F(1:nCell+1),&
              LeftT_F(1:nCell+1), RightT_F(1:nCell+1),&
              u_F(1:nCell+1),tFlux_F(1:nCell+1))
  endif

  ! get temperature from pressure and density, and calculate heat flow
  ! source term.
  T_G(-1:nCell+2)=OldState_GV(-1:nCell+2,P_)/Rgas/OldState_GV(-1:nCell+2,Rho_)
  
  do iCell = 0,nCell+1
     GradT_F(iCell)=(T_G(iCell+1)-T_G(iCell)) / (DrBnd)
  enddo

  !get facevalues
  do iCell = 0,nCell
     Conduction_F(iCell) = (HeatCon_G(iCell+1)+HeatCon_G(iCell)) / 2.0 &
          * GradT_F(iCell)
  enddo
  
  do iCell = 1,nCell
     Diffusion_C(iCell) = &
          (Ar23(iCell)*Conduction_F(iCell)-Ar12(iCell)*Conduction_F(iCell-1))&
          / (CellVolume_C(iCell))
!     if (iCell == nCell)  then
!        write(*,*) 'Gmin1/Rgas/OldState_GV(nCell,Rho_)*Diffusion_C(nCell):',&
!             Gmin1/Rgas/OldState_GV(nCell,Rho_)*Diffusion_C(nCell)
!        write(*,*) 'Ar23(nCell)        :',Ar23(nCell)
!        write(*,*)'Conduction_F(nCell)      :',Conduction_F(nCell)
!        write(*,*)'Ar12(nCell)         :',Ar12(nCell)
!        write(*,*)'Conduction_F(nCell-1)    :',Conduction_F(nCell-1)
!        write(*,*)'CellVolume_C(nCell) :',CellVolume_C(nCell)
!        write(*,*)'GradT_F(nCell+1)    :',GradT_F(nCell+1)
!        write(*,*)'T_G(nCell+1)        :',T_G(nCell+1)
!        write(*,*)'OldState_GV(nCell+1,P_)',OldState_GV(nCell+1,P_)
!        write(*,*)'------------------------'
!     end if
  enddo
  
  ! Get terms that need to be flux limited
  call calc_facevalues(nCell, OldState_GV(-1:nCell+2,U_), LeftU_F, RightU_F)
  call calc_facevalues(nCell, T_G(-1:nCell+2), LeftT_F, RightT_F)
  Coeff = 0.47*DrBnd/DtIn
  
  u_F(1:nCell+1)  = 0.5*( LeftU_F(1:nCell+1) + RightU_F(1:nCell+1) )
  
  tFlux_F(1:nCell+1)  = &
       0.5*( LeftU_F(1:nCell+1)*LeftT_F(1:nCell+1) &
            + RightU_F(1:nCell+1)*RightT_F(1:nCell+1) )!&
!       - Coeff * (RightT_F(1:nCell+1) - LeftT_F(1:nCell+1))
  
  do iCell = 1,nCell
     DivU_C (iCell) = &
          (Ar23(iCell)*U_F(iCell+1) - Ar12(iCell)*U_F(iCell)) &
          / (CellVolume_C(iCell))
     DivUT_C(iCell) = &
          (Ar23(iCell)*tFlux_F(iCell+1) - Ar12(iCell)*tFlux_F(iCell)) &
          / (CellVolume_C(iCell))
  enddo

  !Update Temperature state
  NewT_G(1:nCell) = T_G(1:nCell) + DtIn *(                                &
       Gmin1/Rgas/OldState_GV(1:nCell,Rho_)*eSource_C(1:nCell)            &
       +Diffusion_C(1:nCell)/OldState_GV(1:nCell,Rho_)                     &
       -DivUT_C(1:nCell)                                                  &
       -T_G(1:nCell)*RhoSource_C(1:nCell)/OldState_GV(1:nCell,Rho_)       &
       -T_G(1:nCell)*(Gamma-2.0)*DivU_C(1:nCell) )

!  write(*,*) 'OLD T:',T_G(nCell-1:nCell)
!  write(*,*) 'RHS  :', DtIn *(                                    &
!       Gmin1/Rgas/OldState_GV(nCell,Rho_)*                                  &
!           (Diffusion_C(nCell) + eSource_C(nCell))                       &
!       -DivUT_C(nCell)                        &
!       -T_G(nCell)*RhoSource_C(nCell)/OldState_GV(nCell,Rho_)           &
!       -T_G(nCell)*(Gamma-2.0)*DivU_C(nCell) )
!!
!  write(*,*) 'Gmin1/Rgas/OldState_GV(nCell,Rho_)*Diffusion_C(nCell) :',&
!       Gmin1/Rgas/OldState_GV(nCell,Rho_)*Diffusion_C(nCell)
!  write(*,*) 'Gmin1/Rgas/OldState_GV(nCell,Rho_)*eSource_C(nCell)    :',&
!       Gmin1/Rgas/OldState_GV(nCell,Rho_)*eSource_C(nCell)
!  write(*,*) '-DivUT_C(nCell)                                       :',-DivUT_C(nCell) 
!  write(*,*) '-T_G(nCell)*RhoSource_C(nCell)/OldState_GV(nCell,Rho_) :',&
!       -T_G(nCell)*RhoSource_C(nCell)/OldState_GV(nCell,Rho_)
!  write(*,*) '-T_G(nCell)*(Gamma-2.0)*DivU_C(nCell)                 :',-T_G(nCell)*(Gamma-2.0)*DivU_C(nCell) 
!  write(*,*) '-------------------'
!  write(*,*) NewT_G(nCell),T_G(nCell),DivU_C(nCell),DivUT_C(nCell),&
!       OldState_GV(nCell,Rho_),RhoSource_C(nCell),eSource_C(nCell),Diffusion_C(nCell)


  deallocate(T_G,Conduction_G,Conduction_F,     &
              GradT_F,Diffusion_C,      &
              DivU_C,DivUT_C, &
              LeftU_F,RightU_F,&
              LeftT_F, RightT_F,&
              u_F,tFlux_F)

end subroutine PW_eheat_flux_explicit
