!^CFG COPYRIGHT UM
!=========================================================================!
!BOP
!MODULE: ionosphere_bc - converts IE  potential to GM  boundary velocity 
!DESCRIPTION:               
!             Ionosphere/magnetosphere coupling                           
!  The data for potential at the GMIE_ grid are transformed  to the       
!  ionosphere velocities by ionosphere_magBCs. Then for any point near the
!  sphere of the radius of IONO_Radius_Mag_Boundary, the subroutine       
!  calc_inner_BC_velocities gives the three cartesian components of the   
!  ionosphere velocities using bilinear interpolation                     
!EOP
!=============================================================================

!BOP
!INTERFACE:
subroutine ionosphere_magBCs(&
     PHI_BC, ETh_BC,EPs_BC,UR_BC, UTh_BC, UPs_BC, Radius_BC,              &
     Theta_BC, Psi_BC, nTheta, nPsi)
  !USES:
  use ModPhysics,ONLY:UnitSI_B,UnitSI_U
  use ModInnerBC,ONLY:&
       IONO_NORTH_THETA,IONO_NORTH_PSI,& 
       IONO_SOUTH_THETA,IONO_SOUTH_PSI,& 
       IONO_PI,IONO_Radius, IONO_Toler
!EOP
  implicit none

  integer, intent(in)  :: nTheta, nPsi
  real, intent(in)     :: Radius_BC

  real, dimension(nTheta,nPsi),intent(in)     :: &
        PHI_BC,Theta_BC, Psi_BC
  real, dimension(nTheta,nPsi), intent(out)   :: &
       ETh_BC, EPs_BC ,UR_BC, UTh_BC, UPs_BC                      
          
  real, dimension(nTheta) :: dTheta
  real, dimension(nPsi)   :: dPsi

  logical :: north
  integer :: i, j
  real :: dTheta_l, dPsi_l, dTheta_l2, dPsi_l2
  real :: dd, dd2 
  real :: cosTheta, sinTheta, cosPhi, sinPhi
  real :: xx, yy, zz
  real :: ER, Ex, Ey, Ez
  real :: Ux, Uy, Uz
  real :: BB, Bx, By,Bz
 

  real, dimension(3) :: Mag_Loc, Mag_B


  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------
  call set_oktest('ionosphere_magBCs',DoTest,DotestMe)

  if(DoTestMe)write(*,*)'ionosphere_magBCs starting with iono_pi=',IONO_PI

  if (Theta_BC(1,1) < IONO_PI/4.0) then
     north = .true.
  else
     north = .false.
  end if

  dTheta_l=(Theta_BC(nTheta,1)-Theta_BC(1,1))/real(nTheta-1)
  dPsi_l=(Psi_BC(1,nPsi)-Psi_BC(1,1))/real(nPsi-1)
  dTheta_l2=dTheta_l*dTheta_l
  dPsi_l2=dPsi_l*dPsi_l
  dd=dPsi_l
  dd2=dPsi_l2

  if(DoTestMe)write(*,*)'dTheta_l,dPsi_l=',dTheta_l,dPsi_l

  ! Determine the potential at the magnetospheric inner boundary.



  ! Compute the electric field at the magnetospheric inner boundary.

  if(DoTestMe)write(*,*)'ionosphere_magBCs compute electric field'

  do j = 1, nPsi
     if (j > 1 .and. j < nPsi ) then 
        do i = 2, nTheta-1
           sinTheta = sin(Theta_BC(i,j))
           ETh_BC(i,j) = -(PHI_BC(i+1,j)-PHI_BC(i-1,j))/ &
                         (2.00*dd*Radius_BC)
           EPs_BC(i,j) = -(PHI_BC(i,j+1)-PHI_BC(i,j-1))/ &
                         (2.00*dd*Radius_BC*sinTheta)
        end do
        ETh_BC(1,j) = -(PHI_BC(2,j)-PHI_BC(1,j))/ &
                      (dd*Radius_BC)
        EPs_BC(1,j) = EPs_BC(2,j)
        ETh_BC(nTheta,j) = -(PHI_BC(nTheta,j)-PHI_BC(nTheta-1,j))/ &
                           (dd*Radius_BC)
        EPs_BC(nTheta,j) = EPs_BC(nTheta-1,j)
     else if (j == 1) then
        do i = 2, nTheta-1
           sinTheta = sin(Theta_BC(i,j))
           ETh_BC(i,j) = -(PHI_BC(i+1,j)-PHI_BC(i-1,j))/ &
                         (2.00*dd*Radius_BC)
           EPs_BC(i,j) = -(PHI_BC(i,j+1)-PHI_BC(i,nPsi-1))/ &
                         (2.00*dd*Radius_BC*sinTheta)
        end do
        ETh_BC(1,j) = -(PHI_BC(2,j)-PHI_BC(1,j))/ &
                      (dd*Radius_BC)
        EPs_BC(1,j) = EPs_BC(2,j)
        ETh_BC(nTheta,j) = -(PHI_BC(nTheta,j)-PHI_BC(nTheta-1,j))/ &
                           (dd*Radius_BC)
        EPs_BC(nTheta,j) = EPs_BC(nTheta-1,j)
     else
        do i = 2, nTheta-1
           sinTheta = sin(Theta_BC(i,j))
           ETh_BC(i,j) = -(PHI_BC(i+1,j)-PHI_BC(i-1,j))/ &
                         (2.00*dd*Radius_BC)
           EPs_BC(i,j) = -(PHI_BC(i,2)-PHI_BC(i,j-1))/ &
                         (2.00*dd*Radius_BC*sinTheta)
        end do
        ETh_BC(1,j) = -(PHI_BC(2,j)-PHI_BC(1,j))/ &
                      (dd*Radius_BC)
        EPs_BC(1,j) = EPs_BC(2,j)
        ETh_BC(nTheta,j) = -(PHI_BC(nTheta,j)-PHI_BC(nTheta-1,j))/ &
                           (dd*Radius_BC)
        EPs_BC(nTheta,j) = EPs_BC(nTheta-1,j)
     end if
  end do

  ! Compute the convection velocities at the magnetospheric inner boundary.

  if(DoTestMe)write(*,*)'ionosphere_magBCs compute velocities'

  do j = 1, nPsi
     do i = 1, nTheta
        cosTheta = cos(Theta_BC(i,j))
        sinTheta = sin(Theta_BC(i,j))
        cosPhi = cos(Psi_BC(i,j))
        sinPhi = sin(Psi_BC(i,j))

        if (north .and. i == nTheta) then
           ER = 0.00
        else if (.not.north .and. i == 1) then
           ER = 0.00
        else
           ER = -0.50*(sinTheta/(cosTheta+IONO_Toler**2))*ETh_BC(i,j)
        end if
      
        Ex = ER*sinTheta*cosPhi + ETh_BC(i,j)*cosTheta*cosPhi - &
             EPs_BC(i,j)*sinPhi
        Ey = ER*sinTheta*sinPhi + ETh_BC(i,j)*cosTheta*sinPhi + &
             EPs_BC(i,j)*cosPhi
        Ez = ER*cosTheta - ETh_BC(i,j)*sinTheta
        
        Mag_Loc(1) = Radius_BC * sin(Theta_BC(i,j))* &
                                 cos(Psi_BC(i,j)) / IONO_Radius
        Mag_Loc(2) = Radius_BC * sin(Theta_BC(i,j))* &
                                 sin(Psi_BC(i,j)) / IONO_Radius
        Mag_Loc(3) = Radius_BC * cos(Theta_BC(i,j)) / IONO_Radius

        call get_b0(Mag_Loc(1),Mag_Loc(2),Mag_Loc(3),Mag_B)

        Mag_B = Mag_B * UnitSI_B

        BB = sum(Mag_B**2)

        bx = Mag_B(1)
        by = Mag_B(2)
        bz = Mag_B(3)
        
        Ux = (Ey*bz - Ez*by)/BB
        Uy = (Ez*bx - Ex*bz)/BB
        Uz = (Ex*by - Ey*bx)/BB
        
        xx = sinTheta*cosPhi
        yy = sinTheta*sinPhi
        zz = cosTheta

        UR_BC(i,j) = Ux*xx + Uy*yy + Uz*zz
        UTh_BC(i,j) = ((Ux*xx + Uy*yy)*zz - &
                       Uz*(xx**2+yy**2)) / &
                      sqrt(xx**2+yy**2+IONO_TOLER**2)
        UPs_BC(i,j) = (Uy*xx - Ux*yy)*sinTheta / &
                      (xx**2+yy**2+IONO_TOLER**2)

        UR_BC(i,j) = UR_BC(i,j) / UnitSI_U
        UTh_BC(i,j) = UTh_BC(i,j) / UnitSI_U
        UPs_BC(i,j) = UPs_BC(i,j) / UnitSI_U

     end do  
  end do

end subroutine ionosphere_magBCs

!==========================================================================
subroutine calc_inner_bc_velocity(nIter,tSimulation,Xyz_D,B1_D,B0_D,u_D)

  use ModIonoPotential
  use ModMain,           ONLY: TypeCoordSystem, nDim
  use CON_axes,          ONLY: transform_matrix
  use ModCoordTransform, ONLY: xyz_to_dir, cross_product
  use CON_planet_field,  ONLY: map_planet_field

  implicit none
  integer, intent(in) :: nIter
  real, intent(in)    :: tSimulation
  real, intent(in)    :: Xyz_D(nDim)    ! Position vector
  real, intent(in)    :: B1_D(nDim)     ! Magnetic field perturbation
  real, intent(in)    :: B0_D(nDim)     ! Magnetic field of planet
  real, intent(out)   :: u_D(nDim)      ! Velocity vector (output)

  real, parameter :: Epsilon = 0.01 ! Perturbation of X, Y or Z

  real :: XyzEpsilon_D(nDim)       ! Points shifted by Epsilon
  real :: XyzIono_D(nDim)          ! Mapped point on the ionosphere
  real :: Theta, Phi               ! Mapped point colatitude, longitude
  real :: ThetaNorm, PhiNorm       ! Normalized colatitude, longitude
  real :: Dist1, Dist2             ! Distance from ionosphere grid point

  real :: Potential_DI(nDim, 2)    ! Potential at the shifted positions
  real :: eField_D(nDim)           ! Electric field
  real :: b_D(nDim)                ! Magnetic field
  real :: B2                       ! Magnetic field squared

  integer :: iDim, iSide, iTheta, iPhi, iHemisphere

  character(len=*), parameter :: NameSub = 'calc_inner_bc_velocity'
  logical :: DoTest, DoTestMe
  real :: tSimulationLast=-1.0
  real, save :: SmgGm_DD(nDim,nDim)
  !-------------------------------------------------------------------------

  call set_oktest(NameSub, DoTest, DoTestMe)

  if(DoTestMe)write(*,*)NameSub,' Xyz_D=',Xyz_D

  ! Calculate conversion matrix between GM and SMG coordinates
  if( tSimulationLast /= tSimulation ) then
     tSimulationLast = tSimulation
     SmgGm_DD = transform_matrix(tSimulation, TypeCoordSystem, 'SMG')
  end if

  ! Map points to obtain potential
  do iSide = 1, 2
     do iDim = 1, nDim

        ! Perturb the iDim coordinate
        XyzEpsilon_D = Xyz_D
        if(iSide == 1)then
           XyzEpsilon_D(idim) = XyzEpsilon_D(idim) - Epsilon
        else
           XyzEpsilon_D(idim) = XyzEpsilon_D(idim) + Epsilon
        end if

        ! Transform into SMG coordinates
        XyzEpsilon_D = matmul(SmgGm_DD, XyzEpsilon_D)

        ! Map down to the ionosphere at radius rIonosphere
        call map_planet_field(tSimulation, XyzEpsilon_D, 'SMG NORM', &
             rIonosphere, XyzIono_D, iHemisphere)

        ! Calculate angular coordinates
        call xyz_to_dir(XyzIono_D, Theta, Phi)

        ! Interpolate potential
        ThetaNorm = Theta / dThetaIono
        PhiNorm   = Phi   / dPhiIono

        iTheta    = floor(ThetaNorm) + 1
        iPhi      = floor(PhiNorm)   + 1

        if(iTheta<1 .or. iTheta > nThetaIono .or. &
             iPhi < 1 .or. iPhi > nPhiIono)then
           write(*,*)NameSub,' PhiNorm, ThetaNorm=',PhiNorm,ThetaNorm
           write(*,*)NameSub,' Phi, Theta=',Phi,Theta
           write(*,*)NameSub,' nPhi, nTheta=',nPhiIono,nThetaIono
           write(*,*)NameSub,' iPhi, iTheta=',iPhi,iTheta
           call stop_mpi(NameSub//' index out of bounds')
        end if

        Dist1     = ThetaNorm - (iTheta - 1)
        Dist2     = PhiNorm   - (iPhi   - 1)

        Potential_DI(iDim, iSide) = &
             (1 - Dist1)*( (1-Dist2) * IonoPotential_II(iTheta  , iPhi  )  &
             +             Dist2     * IonoPotential_II(iTheta  , iPhi+1)) &
             + Dist1    *( (1-Dist2) * IonoPotential_II(iTheta+1, iPhi  )  &
             +             Dist2     * IonoPotential_II(iTheta+1, iPhi+1))

        if(DoTestMe)then
           write(*,*)NameSub,' iDim, iSide  =',iDim, iSide
           write(*,*)NameSub,' XyzEpsilon_D =',XyzEpsilon_D
           write(*,*)NameSub,' XyzIono_D    =',XyzIono_D
           write(*,*)NameSub,' Theta, Phi   =',Theta,Phi
           write(*,*)NameSub,' iTheta, iPhi =',iTheta,iPhi
           write(*,*)NameSub,' Dist1, Dist2 =',Dist1,Dist2
           write(*,*)NameSub,' Potential_DI =',Potential_DI(iDim, iSide)
        end if

     end do
  end do

  b_D = B1_D + B0_D
  B2  = sum(b_D**2)
  
  ! E = -grad(Potential)
  eField_D = - (Potential_DI(:,2) - Potential_DI(:,1))/(2*Epsilon)

  ! U = (E x B) / B^2
  u_D = cross_product(eField_D, b_D) / B2

  if(DoTestMe)then
     write(*,*)NameSub,' b_D=',b_D
     write(*,*)NameSub,' E_D=',eField_D
     write(*,*)NameSub,' u_D=',u_D
  endif

  ! Subtract the radial component of the velocity
  u_D = u_D - Xyz_D * sum(Xyz_D * u_D) / sum(Xyz_D**2)

  if(DoTestMe)then
     write(*,*)NameSub,' Final u_D=',u_D
     !!! if(maxval(abs(u_D))>0.00001)call stop_mpi(NameSub)
  end if

end subroutine calc_inner_bc_velocity

!==========================================================================
subroutine calc_inner_bc_velocity2(nIter,tSimulation,Xyz_D,B1_D,B0_D,u_D)

  use ModIonoPotential
  use ModMain,           ONLY: TypeCoordSystem, nDim
  use CON_axes,          ONLY: transform_matrix
  use ModCoordTransform, ONLY: xyz_to_dir, cross_product
  use CON_planet_field,  ONLY: map_planet_field

  implicit none
  integer, intent(in) :: nIter
  real, intent(in)    :: tSimulation
  real, intent(in)    :: Xyz_D(nDim)    ! Position vector
  real, intent(in)    :: B1_D(nDim)     ! Magnetic field perturbation
  real, intent(in)    :: B0_D(nDim)     ! Magnetic field of planet
  real, intent(out)   :: u_D(nDim)      ! Velocity vector (output)

  real :: XyzIono_D(nDim)      ! Mapped point on the ionosphere
  real :: Theta, Phi           ! Mapped point colatitude, longitude
  real :: ThetaNorm, PhiNorm   ! Normalized colatitude, longitude
  real :: Dist1, Dist2         ! Distance from ionosphere grid point

  real :: dPotential_D(2)      ! Gradient of potential at the mapped position
  real :: DdirDxyz_DD(2,3)     ! Jacobian matrix between Theta, Phi and Xyz_D
  real :: eField_D(nDim)       ! Electric field
  real :: b_D(nDim)            ! Magnetic field
  real :: B2                   ! Magnetic field squared

  integer :: iTheta, iPhi, iHemisphere

  character(len=*), parameter :: NameSub = 'calc_inner_bc_velocity'
  logical :: DoTest, DoTestMe
  real :: tSimulationLast=-1.0
  !-------------------------------------------------------------------------

  call set_oktest(NameSub, DoTest, DoTestMe)

  if(DoTestMe)write(*,*)NameSub,' Xyz_D=',Xyz_D

  ! Map down to the ionosphere at radius rIonosphere. Result is in SMG.
  ! Also obtain the Jacobian matrix between Theta,Phi and Xyz_D
  call map_planet_field(tSimulation, Xyz_D, TypeCoordSystem//' NORM', &
       rIonosphere, XyzIono_D, iHemisphere, .true., DdirDxyz_DD)

  ! Calculate angular coordinates
  call xyz_to_dir(XyzIono_D, Theta, Phi)

  ! Interpolate the spherical gradients of the electric potential
  ThetaNorm = Theta / dThetaIono
  PhiNorm   = Phi   / dPhiIono

  iTheta    = floor(ThetaNorm) + 1
  iPhi      = floor(PhiNorm)   + 1

  if(iTheta<1 .or. iTheta > nThetaIono .or. &
       iPhi < 1 .or. iPhi > nPhiIono)then
     write(*,*)NameSub,' PhiNorm, ThetaNorm=',PhiNorm,ThetaNorm
     write(*,*)NameSub,' Phi, Theta=',Phi,Theta
     write(*,*)NameSub,' nPhi, nTheta=',nPhiIono,nThetaIono
     write(*,*)NameSub,' iPhi, iTheta=',iPhi,iTheta
     call stop_mpi(NameSub//' index out of bounds')
  end if

  Dist1     = ThetaNorm - (iTheta - 1)
  Dist2     = PhiNorm   - (iPhi   - 1)

  dPotential_D = &
       (1 - Dist1)*( (1-Dist2) * dIonoPotential_DII(:, iTheta  , iPhi  )  &
       +             Dist2     * dIonoPotential_DII(:, iTheta  , iPhi+1)) &
       + Dist1    *( (1-Dist2) * dIonoPotential_DII(:, iTheta+1, iPhi  )  &
       +             Dist2     * dIonoPotential_DII(:, iTheta+1, iPhi+1))

  if(DoTestMe)then
     write(*,*)NameSub,' Xyz_D        =',Xyz_D
     write(*,*)NameSub,' XyzIono_D    =',XyzIono_D
     write(*,*)NameSub,' Theta, Phi   =',Theta,Phi
     write(*,*)NameSub,' iTheta, iPhi =',iTheta,iPhi
     write(*,*)NameSub,' Dist1, Dist2 =',Dist1,Dist2
     write(*,*)NameSub,' dPotential_D =',dPotential_D
  end if

  ! E = -grad(Potential) = - dPotential/d(Theta,Phi) * d(Theta,Phi)/d(x,y,z)
  eField_D = - matmul( dPotential_D, DdirDxyz_DD)

  ! Magnetic field
  b_D = B1_D + B0_D
  B2  = sum(b_D**2)

  ! U = (E x B) / B^2
  u_D = cross_product(eField_D, b_D) / B2

  if(DoTestMe)then
     write(*,*)NameSub,' b_D=',b_D
     write(*,*)NameSub,' E_D=',eField_D
     write(*,*)NameSub,' u_D=',u_D
  endif

  ! Subtract the radial component of the velocity
  u_D = u_D - Xyz_D * sum(Xyz_D * u_D) / sum(Xyz_D**2)

  if(DoTestMe)then
     write(*,*)NameSub,' Final u_D=',u_D
     !!! if(maxval(abs(u_D))>0.00001)call stop_mpi(NameSub)
  end if

end subroutine calc_inner_bc_velocity2

!============================================================================  

!BOP
!INTERFACE:
subroutine calc_inner_BC_velocities(iter,time_now,XFace,YFace,ZFace,&
     BxOutside,ByOutside,BzOutside,B0faceX,B0FaceY,B0FaceZ,&
     VxFace,VyFace,VzFace)
!USES:
  use ModNumConst
  use ModInnerBC
!DESCRIPTION:
  !\
  ! This subroutine should calculate -x -y and -z components of the "ionosphere" velocities
  ! as a function of the cartesian ponit coordinates x,y,z, as well as the 
  ! induction and static magnetic field
  ! Iteration number and time can be also used.
  !/
!EOP
  implicit none

  integer,intent(in) :: iter
  real, intent(in) :: time_now
  real,intent(in) :: XFace,YFace,ZFace,&
       BxOutside,ByOutside,BzOutside,B0faceX,B0FaceY,B0FaceZ
  real,intent(out) :: VxFace,VyFace,VzFace
  real:: VrFace,VThetaFace,VphiFace
  integer::i0,j0
  real:: RFace, cosTheta, sinTheta, cosPhi, sinPhi
  real :: ThetaFace, PhiFace, dTheta, dPsi, dd

  real :: xC1, yC1, xC2, yC2, xC3, yC3, xC4, yC4, &
       fC1, fC2, fC3, fC4, &
       a_bilin, b_bilin, c_bilin, d_bilin

  RFace=sqrt(XFace*XFace+YFace*YFace+ZFace*ZFace)
  cosTheta = ZFace/RFace
  sinTheta = sqrt(XFace**2+YFace**2)/RFace
  cosPhi = XFace/sqrt(XFace**2+YFace**2+cTolerance**2)
  sinPhi = YFace/sqrt(XFace**2+YFace**2+cTolerance**2)

  if (ZFace > cZero) then
     !\
     ! Northern Hemisphere Ionosphere BC.
     !/
     dTheta = (IONO_NORTH_Theta(IONO_nTheta,1)- &
          IONO_NORTH_Theta(1,1))/real(IONO_nTheta-1)
     dPsi = (IONO_NORTH_Psi(1,IONO_nPsi)- &
          IONO_NORTH_Psi(1,1))/real(IONO_nPsi-1)
     dd = dPsi

     if (XFace == cZero .and. YFace == cZero) then
        PhiFace = cZero
     else
        PhiFace = atan2(sinPhi, cosPhi)
     end if
     if (PhiFace < cZero) PhiFace = PhiFace + cTwoPi
     ThetaFace = atan2(sinTheta, cosTheta)

     i0 = ThetaFace/dd + 1
     j0 = PhiFace/dd + 1

     xC1 = IONO_NORTH_Theta(i0,j0)
     yC1 = IONO_NORTH_Psi(i0,j0)
     xC2 = IONO_NORTH_Theta(i0,j0+1)
     yC2 = IONO_NORTH_Psi(i0,j0+1)
     xC3 = IONO_NORTH_Theta(i0+1,j0+1)
     yC3 = IONO_NORTH_Psi(i0+1,j0+1)
     xC4 = IONO_NORTH_Theta(i0+1,j0)
     yC4 = IONO_NORTH_Psi(i0+1,j0)

     fC1 = IONO_NORTH_UR_BC(i0,j0)
     fC2 = IONO_NORTH_UR_BC(i0,j0+1)
     fC3 = IONO_NORTH_UR_BC(i0+1,j0+1)
     fC4 = IONO_NORTH_UR_BC(i0+1,j0)


     a_bilin = fC1
     b_bilin = (fC4-fC1)/dd
     c_bilin = (fC2-fC1)/dd
     d_bilin = (fC3+fC1-fC4-fC2)/(dd*dd)
     VrFace = a_bilin + &
          b_bilin*(ThetaFace-xC1) + &
          c_bilin*(PhiFace-yC1) + &
          d_bilin*(ThetaFace-xC1)*(PhiFace-yC1)

     fC1 = IONO_NORTH_UTh_BC(i0,j0)
     fC2 = IONO_NORTH_UTh_BC(i0,j0+1)
     fC3 = IONO_NORTH_UTh_BC(i0+1,j0+1)
     fC4 = IONO_NORTH_UTh_BC(i0+1,j0)
     a_bilin = fC1
     b_bilin = (fC4-fC1)/dd
     c_bilin = (fC2-fC1)/dd
     d_bilin = (fC3+fC1-fC4-fC2)/(dd*dd)
     VthetaFace = a_bilin + &
          b_bilin*(ThetaFace-xC1) + &
          c_bilin*(PhiFace-yC1) + &
          d_bilin*(ThetaFace-xC1)*(PhiFace-yC1)

     fC1 = IONO_NORTH_UPs_BC(i0,j0)
     fC2 = IONO_NORTH_UPs_BC(i0,j0+1)
     fC3 = IONO_NORTH_UPs_BC(i0+1,j0+1)
     fC4 = IONO_NORTH_UPs_BC(i0+1,j0)
     a_bilin = fC1
     b_bilin = (fC4-fC1)/dd
     c_bilin = (fC2-fC1)/dd
     d_bilin = (fC3+fC1-fC4-fC2)/(dd*dd)
     VphiFace = a_bilin + &
          b_bilin*(ThetaFace-xC1) + &
          c_bilin*(PhiFace-yC1) + &
          d_bilin*(ThetaFace-xC1)*(PhiFace-yC1)

  else
     !\
     ! Southern Hemisphere Ionosphere BC.
     !/
     dTheta = (IONO_SOUTH_Theta(IONO_nTheta,1)- &
          IONO_SOUTH_Theta(1,1))/real(IONO_nTheta-1)
     dPsi = (IONO_SOUTH_Psi(1,IONO_nPsi)- &
          IONO_SOUTH_Psi(1,1))/real(IONO_nPsi-1)
     dd = dPsi

     if (XFace == cZero .and. YFace == cZero) then
        PhiFace = cZero
     else
        PhiFace = atan2(sinPhi, cosPhi)
     end if
     if (PhiFace < cZero) PhiFace = PhiFace + cTwoPi
     ThetaFace = cPi - &
          atan2(sinTheta, -cosTheta)

     i0 = (ThetaFace-cHalfPi)/dd + 1
     j0 = PhiFace/dd + 1

     xC1 = IONO_SOUTH_Theta(i0,j0)
     yC1 = IONO_SOUTH_Psi(i0,j0)
     xC2 = IONO_SOUTH_Theta(i0,j0+1)
     yC2 = IONO_SOUTH_Psi(i0,j0+1)
     xC3 = IONO_SOUTH_Theta(i0+1,j0+1)
     yC3 = IONO_SOUTH_Psi(i0+1,j0+1)
     xC4 = IONO_SOUTH_Theta(i0+1,j0)
     yC4 = IONO_SOUTH_Psi(i0+1,j0)

     fC1 = IONO_SOUTH_UR_BC(i0,j0)
     fC2 = IONO_SOUTH_UR_BC(i0,j0+1)
     fC3 = IONO_SOUTH_UR_BC(i0+1,j0+1)
     fC4 = IONO_SOUTH_UR_BC(i0+1,j0)
     a_bilin = fC1
     b_bilin = (fC4-fC1)/dd
     c_bilin = (fC2-fC1)/dd
     d_bilin = (fC3+fC1-fC4-fC2)/(dd*dd)
     VrFace = a_bilin + &
          b_bilin*(ThetaFace-xC1) + &
          c_bilin*(PhiFace-yC1) + &
          d_bilin*(ThetaFace-xC1)*(PhiFace-yC1)

     fC1 = IONO_SOUTH_UTh_BC(i0,j0)
     fC2 = IONO_SOUTH_UTh_BC(i0,j0+1)
     fC3 = IONO_SOUTH_UTh_BC(i0+1,j0+1)
     fC4 = IONO_SOUTH_UTh_BC(i0+1,j0)
     a_bilin = fC1
     b_bilin = (fC4-fC1)/dd
     c_bilin = (fC2-fC1)/dd
     d_bilin = (fC3+fC1-fC4-fC2)/(dd*dd)
     VthetaFace = a_bilin + &
          b_bilin*(ThetaFace-xC1) + &
          c_bilin*(PhiFace-yC1) + &
          d_bilin*(ThetaFace-xC1)*(PhiFace-yC1)

     fC1 = IONO_SOUTH_UPs_BC(i0,j0)
     fC2 = IONO_SOUTH_UPs_BC(i0,j0+1)
     fC3 = IONO_SOUTH_UPs_BC(i0+1,j0+1)
     fC4 = IONO_SOUTH_UPs_BC(i0+1,j0)
     a_bilin = fC1
     b_bilin = (fC4-fC1)/dd
     c_bilin = (fC2-fC1)/dd
     d_bilin = (fC3+fC1-fC4-fC2)/(dd*dd)
     VphiFace = a_bilin + &
          b_bilin*(ThetaFace-xC1) + &
          c_bilin*(PhiFace-yC1) + &
          d_bilin*(ThetaFace-xC1)*(PhiFace-yC1)
  end if ! South or North

 
  VrFace = cZero
  VxFace  = VrFace*cosPhi*sinTheta + VthetaFace*cosTheta*cosPhi - VphiFace*sinPhi 
  VyFace  = VrFace*sinPhi*sinTheta + VthetaFace*cosTheta*sinPhi + VphiFace*cosPhi 
  VzFace  = VrFace*cosTheta - VthetaFace*sinTheta

end subroutine calc_inner_BC_velocities

!=============================================================================

real function logvar_ionosphere(NameLogvar)
  use ModProcMH,  ONLY: iProc
  use ModIO,      ONLY: write_myname
  use ModInnerBc, ONLY: IONO_NORTH_PHI_BC,IONO_SOUTH_PHI_BC
  implicit none
  character (len=*), intent(in) :: NameLogvar
  integer :: nWarn = 0 ! warn multiple times to catch multiple log variables
  !---------------------------------------------------------------------------
  select case(NameLogvar)
  case('cpcpn','cpcp_n','cpcp_north','cpcpnorth')
     logvar_ionosphere = maxval(IONO_NORTH_PHI_BC)-minval(IONO_NORTH_PHI_BC)
  case('cpcps','cpcp_s','cpcp_south','cpcpsouth')
     logvar_ionosphere = maxval(IONO_SOUTH_PHI_BC)-minval(IONO_SOUTH_PHI_BC)
  case default
     if(nWarn < 2 .and. iProc==0)then
        call write_myname;
        write(*,'(a)')'WARNING in logvar_ionosphere: unknown NameLogvar='//&
             trim(NameLogvar)//', returning -777.77'
        nWarn = nWarn + 1
     end if
     logvar_ionosphere = -777.77
  end select

end function logvar_ionosphere
