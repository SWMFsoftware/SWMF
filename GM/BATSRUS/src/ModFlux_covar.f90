!^CFG COPYRIGHT UM
!^CFG FILE NOT CARTESIAN

Module ModFlux
  use ModVarIndexes,ONLY:Energy_
  implicit none
  integer, parameter :: nFlux=Energy_+2, MaxStrip=10

  ! symmetric face variables

  real,dimension(3,MaxStrip)::FaceArea_DI
  real::FaceArea2Min
  real, dimension(MaxStrip) :: v_B0x,v_B0y,v_B0z,B0n
  real, dimension(MaxStrip) :: v_max_hf          ! Maximum perturbation velocity

  ! left face 
  real, dimension(MaxStrip) :: v_rho_lf,  v_p_lf 
  real, dimension(MaxStrip) :: v_Ux_lf,v_B1x_lf,v_Uy_lf,v_B1y_lf,v_Uz_lf,v_B1z_lf
  real, dimension(MaxStrip) :: v_Un_lf,B1n_lf

  ! right face	
  real, dimension(MaxStrip) :: v_rho_rf, v_p_rf  
  real, dimension(MaxStrip) :: v_Ux_rf,v_B1x_rf,v_Uy_rf,v_B1y_rf,v_Uz_rf,v_B1z_rf
  real, dimension(MaxStrip) :: v_Un_rf,B1n_rf

  ! flux
  real, dimension(nFlux,MaxStrip) :: Full_Flux     ! Flux

  ! These common blocks improve speed slightly
  common/igors/ v_rho_lf, v_Ux_lf,v_Uy_lf,v_Uz_lf,&
       v_B1x_lf, v_B1y_lf, v_B1z_lf, v_p_lf, &
       v_B0x, v_B0y, v_B0z, &
       v_rho_rf, v_Ux_rf, v_Uy_rf, v_Uz_rf,&
       v_B1x_rf, v_B1y_rf, v_B1z_rf, v_p_rf,&
       v_Un_lf,v_Un_rf,B1n_lf,B1n_rf,B0n, &
       Full_Flux,v_max_hf,&
       FaceArea_DI,FaceArea2Min

end module ModFlux











