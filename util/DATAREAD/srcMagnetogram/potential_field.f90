module ModPotentialField

  implicit none

  logical :: DoReadMagnetogram = .true.
  logical :: UseCosTheta = .true. 

  integer         :: nR = 90, nTheta = 90, nPhi = 90
  real, parameter :: rMin = 1.0, rMax = 2.5
  real :: PhiShift

  logical :: UseBr = .true.

  real, dimension(:), allocatable :: &
       Radius_I, Theta_I, Phi_I, SinTheta_I, &
       dRadius_I, dTheta_I, dPhi_I, &
       RadiusNode_I, ThetaNode_I, PhiNode_I, SinThetaNode_I, &
       dRadiusNode_I, dThetaNode_I, dPhiNode_I

  real, allocatable:: Br_II(:,:), Potential_C(:,:,:), Rhs_C(:,:,:), &
       B0_DG(:,:,:,:), DivB_C(:,:,:), PlotVar_VG(:,:,:,:)

contains

  !=================================================================
  subroutine read_magnetogram

    ! Read the raw magnetogram file into a 2d array

    ! Name of input file
    character (len=*), parameter:: NameFileIn = 'fitsfile.dat'

    real, parameter:: BrMax = 3500.0 ! Saturation level of MDI
    integer, parameter:: iUnit = 9   

    integer:: iError
    integer:: nCarringtonRotation
    integer:: nTheta0, nPhi0, nThetaRatio, nPhiRatio
    integer:: iTheta, iPhi, iTheta0, iTheta1, iPhi0, iPhi1
    real :: BrAverage
    character (len=100) :: String

    real, allocatable:: Br0_II(:,:)
    !----------------------------------------------------------
    open(iUnit, file=NameFileIn, status='old', iostat=iError)
    if(iError /= 0)then
       write(*,*) 'Error: could not open input file ',NameFileIn
       stop
    end if
    do 
       read(iUnit,'(a)', iostat = iError ) String
       if(index(String,'#CR')>0)then
          read(iUnit,*) nCarringtonRotation
       endif
       if(index(String,'#ARRAYSIZE')>0)then
          read(iUnit,*) nPhi0
          read(iUnit,*) nTheta0
       endif
       if(index(String,'#START')>0) EXIT
    end do
    
    write(*,*)'nCarringtonRotation, nTheta0, nPhi0: ',&
         nCarringtonRotation, nTheta0, nPhi0

    allocate(Br0_II(nTheta0,nPhi0))

    ! fitsfile.dat is in longitude, latitude
    do iTheta = nTheta0, 1, -1
       do iPhi = 1, nPhi0
          read(iUnit,*) Br0_II(iTheta,iPhi)
          if (abs(Br0_II(iTheta,iPhi)) > BrMax) &
               Br0_II(iTheta,iPhi) = sign(BrMax, Br0_II(iTheta,iPhi))
       end do
    end do

    if(nTheta0 > nTheta)then
       ! Set integer coarsening ratio
       nThetaRatio = nTheta0 / nTheta
       nTheta      = nTheta0 / nThetaRatio
    else
       nThetaRatio = 1
       nTheta      = nTheta0
    end if

    if(nPhi0 > nPhi)then
       nPhiRatio = nPhi0 / nPhi
       nPhi      = nPhi0 / nPhiRatio
    else
       nPhiRatio = 1
       nPhi      = nPhi0
    end if

    allocate(Br_II(nTheta,nPhi))

    do iPhi = 1, nPhi
       iPhi0 = nPhiRatio*(iPhi-1) + 1
       iPhi1 = iPhi0 + nPhiRatio - 1
       do iTheta = 1, nTheta
          iTheta0 = nThetaRatio*(iTheta-1) + 1
          iTheta1 = iTheta0 + nThetaRatio - 1
          
          Br_II(iTheta,iPhi) = sum( Br0_II(iTheta0:iTheta1, iPhi0:iPhi1)) &
               / (nThetaRatio*nPhiRatio)
       end do
    end do

    ! remove monopole
    BrAverage = sum(Br_II)/(nTheta*nPhi)
    Br_II = Br_II - BrAverage

    deallocate(Br0_II)

    ! longitude = 0 is shifted by -PhiShift (BATSRUS will shift by -PhiShift)
    PhiShift = 360.0/nPhi0 *0.5*(nPhiRatio-1)

    close(iUnit)

  end subroutine read_magnetogram

  !============================================================================

  subroutine init_potential_field

    use ModConst, ONLY: cPi, cTwoPi

    integer :: iR, iTheta, iPhi
    real:: dR, dTheta, dPhi, dZ, z
    !--------------------------------------------------------------------------

    allocate( &
         Radius_I(0:nR+1), Theta_I(0:nTheta+1), Phi_I(0:nPhi+1), &
         dRadius_I(nR), dTheta_I(nTheta), dPhi_I(nPhi), &
         SinTheta_I(0:nTheta+1), SinThetaNode_I(nTheta+1), &
         RadiusNode_I(nR+1), ThetaNode_I(nTheta+1), PhiNode_I(nPhi+1), &
         dRadiusNode_I(nR+1), dThetaNode_I(nTheta+1), dPhiNode_I(nPhi+1))

    ! nR is the number of mesh cells in radial direction
    ! cell centered radial coordinate
    dR = (rMax - rMin)/nR
    do iR = 0, nR+1
       Radius_I(iR) = rMin + (iR - 0.5)*dR
    end do
    ! node based radial coordinate
    do iR = 1, nR+1
       RadiusNode_I(iR) = rMin + (iR - 1)*dR
    end do
    dRadius_I = RadiusNode_I(2:nR+1) - RadiusNode_I(1:nR)
    dRadiusNode_I = Radius_I(1:nR+1) - Radius_I(0:nR)

    if(UseCosTheta)then
       dZ = 2.0/nTheta
       do iTheta = 1, nTheta
          z = 1 - (iTheta - 0.5)*dZ
          Theta_I(iTheta) = acos(z)
       end do
       Theta_I(0)        = -Theta_I(1)
       Theta_I(nTheta+1) = cTwoPi - Theta_I(nTheta)
       do iTheta = 1, nTheta + 1
          z = max(-1.0, min(1.0, 1 - (iTheta-1)*dZ))
          ThetaNode_I(iTheta) = acos(z)
       end do
    else
       dTheta = cPi/nTheta
       do iTheta = 0, nTheta+1
          Theta_I(iTheta) = (iTheta - 0.5)*dTheta
       end do

       do iTheta = 1, nTheta+1
          ThetaNode_I(iTheta) = (iTheta - 1)*dTheta
       end do
    end if
    SinTheta_I = sin(Theta_I)
    SinThetaNode_I = sin(ThetaNode_I)

    dTheta_I     = ThetaNode_I(2:nTheta+1) - ThetaNode_I(1:nTheta)
    dThetaNode_I = Theta_I(1:nTheta+1)     - Theta_I(0:nTheta)

    dPhi = cTwoPi/nPhi
    do iPhi = 0, nPhi+1
       Phi_I(iPhi) = (iPhi - 0.5)*dPhi
    end do
    do iPhi = 1, nPhi+1
       PhiNode_I(iPhi) = (iPhi - 1)*dPhi
    end do
    dPhi_I = PhiNode_I(2:nPhi+1) - PhiNode_I(1:nPhi)
    dPhiNode_I = Phi_I(1:nPhi+1) - Phi_I(0:nPhi)

    allocate( &
         Potential_C(nR,nTheta,nPhi), &
         Rhs_C(nR,nTheta,nPhi), &
         B0_DG(3,nR+1,nTheta+1,nPhi+1), &
         DivB_C(nR,nTheta,nPhi), &
         PlotVar_VG(6,nR+1,nTheta+1,nPhi+1))

    Potential_C       =   0.0
    Rhs_C             =   0.0

  end subroutine init_potential_field

  !============================================================================

  subroutine save_potential_field

    use ModInterpolate, ONLY: trilinear
    use ModNumConst,    ONLY: cHalfPi
    use ModPlotFile,    ONLY: save_plot_file

    real, allocatable :: Potential_G(:,:,:)
    real, allocatable :: Br_G(:,:,:), rBtheta_G(:,:,:), rBphi_G(:,:,:)
    real, allocatable :: B_DN(:,:,:,:)

    real :: r, Theta, Phi
    integer :: iR, iTheta, iPhi
    !--------------------------------------------------------------------------

    allocate(Potential_G(0:nR+1,0:nTheta+1,0:nPhi+1))
    Potential_G = 0.0
    call set_boundary(Potential_C, Potential_G)

    ! Staggered components of the magnetic field
    allocate( &
         Br_G(1:nR+1,0:nTheta+1,0:nPhi+1), &
         rBtheta_G(0:nR+1,1:nTheta+1,0:nPhi+1), &
         rBphi_G(0:nR+1,0:nTheta+1,1:nPhi+1) )

    do iPhi = 0, nPhi+1
       do iTheta = 0, nTheta+1
          do iR = 1, nR+1
             Br_G(iR,iTheta,iPhi) = &
                  (Potential_G(iR,iTheta,iPhi)-Potential_G(iR-1,iTheta,iPhi)) &
                  / dRadiusNode_I(iR)
          end do
       end do
    end do

    do iPhi = 0, nPhi+1
       do iTheta = 1, nTheta+1
          do iR = 0, nR+1
             rBtheta_G(iR,iTheta,iPhi) = &
                  (Potential_G(iR,iTheta,iPhi)-Potential_G(iR,iTheta-1,iPhi)) &
                  / dThetaNode_I(iTheta)
          end do
       end do
    end do

    do iPhi = 1, nPhi+1
       do iTheta = 0, nTheta+1
          do iR = 0, nR+1
             rBphi_G(iR,iTheta,iPhi) = &
                  (Potential_G(iR,iTheta,iPhi)-Potential_G(iR,iTheta,iPhi-1)) &
                  / (SinTheta_I(iTheta)*dPhiNode_I(iPhi))
          end do
       end do
    end do

    ! The magnetic field on the final grid used in BATSRUS
    ! Note, the coordinates are in longitutude and latitude, but the
    ! magnetic field is in Bphi, Btheta = -Blatitude
    allocate(B_DN(3,1:nR+1,1:nPhi+1,1:nTheta))

    do iPhi = 1, nPhi+1
       do iTheta = 1, nTheta
          do iR = 1, nR+1
             r = RadiusNode_I(iR)
             Theta = Theta_I(iTheta)
             Phi = Phi_I(iPhi)

             B_DN(1,iR,iPhi,nTheta+1-iTheta) = Br_G(iR,iTheta,iPhi)
             B_DN(3,iR,iPhi,nTheta+1-iTheta) = &
                  trilinear(rBtheta_G, 0, nR+1, 1, nTheta+1, 0, nPhi+1, &
                  (/r,Theta,Phi/), x_I=Radius_I, y_I=ThetaNode_I, &
                  z_I=Phi_I, DoExtrapolate=.false.) / r
             B_DN(2,iR,iPhi,nTheta+1-iTheta) = &
                  trilinear(rBphi_G, 0, nR+1, 0, nTheta+1, 1, nPhi+1, &
                  (/r,Theta,Phi/), x_I=Radius_I, y_I=Theta_I, &
                  z_I=PhiNode_I, DoExtrapolate=.false.) / r

          end do
       end do
    end do

    call save_plot_file('potentialfield.out', TypeFileIn='real8', &
         StringHeaderIn = 'Radius [Rs] Longitude [Rad] Latitude [Rad] B [G]', &
         nameVarIn = 'Radius Longitude Latitude Br Bphi Btheta' &
         //' Ro_PFSSM Rs_PFSSM PhiShift nRExt', &
         ParamIn_I = (/ rMin, rMax, PhiShift, 0.0 /), &
         nDimIn=3, VarIn_VIII=B_DN, &
         Coord1In_I=RadiusNode_I, &
         Coord2In_I=PhiNode_I, &
         Coord3In_I=cHalfPi-Theta_I(nTheta:1:-1))

  end subroutine save_potential_field

  !============================================================================

  subroutine set_boundary(x_C, x_G)

    real, intent(in):: x_C(nR,nTheta,nPhi)
    real, intent(inout):: x_G(0:nR+1,0:nTheta+1,0:nPhi+1)

    integer :: iPhi, jPhi
    !--------------------------------------------------------------------------
    ! Current solution inside
    x_G(1:nR,1:nTheta,1:nPhi) = x_C

    ! The slope is forced to be Br at the inner boundary
    if(UseBr)then
       x_G(0,1:nTheta,1:nPhi) = x_C(1,:,:) - dRadiusNode_I(1)*Br_II
    else
       x_G(0,1:nTheta,1:nPhi) = x_C(1,:,:)
    end if

    ! The potential is zero at the outer boundary
    x_G(nR+1,:,:) = -x_G(nR,:,:)

    ! Symmetric in Theta but shifted by nPhi/2
    do iPhi = 1, nPhi
       jPhi = modulo(iPhi-1 + nPhi/2, nPhi) + 1
       x_G(:,0,iPhi)        = x_G(:,1,jPhi)
       x_G(:,nTheta+1,iPhi) = x_G(:,nTheta,jPhi)
    end do

    ! Periodic around phi
    x_G(:,:,0)      = x_G(:,:,nPhi)
    x_G(:,:,nPhi+1) = x_G(:,:,1)

  end subroutine set_boundary

end module ModPotentialField

!==============================================================================

module ModB0Matvec

  implicit none

contains

  !============================================================================

  subroutine matvec(x_C, y_C, n)

    use ModPotentialField, ONLY: B0_DG

    integer, intent(in) :: n
    real, intent(in)    :: x_C(n)
    real, intent(out)   :: y_C(n)
    !--------------------------------------------------------------------------

    ! Calculate y = laplace x in two steps
    call get_gradient(x_C, B0_DG)
    call get_divergence(B0_DG, y_C)

  end subroutine matvec

  !============================================================================

  subroutine get_gradient(x_C, Grad_DG)

    use ModPotentialField, ONLY: nR, nTheta, nPhi, Radius_I, SinTheta_I, &
         dRadiusNode_I, dThetaNode_I, dPhiNode_I, Br_II, set_boundary

    real, intent(in):: x_C(nR,nTheta,nPhi)
    real, intent(out):: Grad_DG(3,nR+1,nTheta+1,nPhi+1)

    real, allocatable, save :: x_G(:,:,:)

    integer:: iR, iTheta, iPhi
    !--------------------------------------------------------------------------
    if(.not.allocated(x_G))then
       allocate(x_G(0:nR+1,0:nTheta+1,0:nPhi+1))
       ! Initialize so that corners are all set
       x_G = 0.0
    end if

    call set_boundary(x_C, x_G)

    ! This initialization is only for the corners
    Grad_DG = 0.0

    do iPhi = 1, nPhi
       do iTheta = 1, nTheta
          do iR = 1, nR+1
             Grad_DG(1,iR,iTheta,iPhi) = &
                  (x_G(iR,iTheta,iPhi) - x_G(iR-1,iTheta,iPhi)) &
                  / dRadiusNode_I(iR)
          end do
       end do
    end do

    do iPhi = 1, nPhi
       do iTheta = 1, nTheta+1
          do iR = 1, nR
             Grad_DG(2,iR,iTheta,iPhi) = &
                  (x_G(iR,iTheta,iPhi) - x_G(iR,iTheta-1,iPhi)) &
                  / (Radius_I(iR)*dThetaNode_I(iTheta))
          end do
       end do
    end do

    do iPhi = 1, nPhi+1
       do iTheta = 1, nTheta
          do iR = 1, nR
             Grad_DG(3,iR,iTheta,iPhi) = &
                  (x_G(iR,iTheta,iPhi) - x_G(iR,iTheta,iPhi-1)) &
                  / (Radius_I(iR)*SinTheta_I(iTheta)*dPhiNode_I(iPhi))
          end do
       end do
    end do

  end subroutine get_gradient

  !============================================================================

  subroutine get_divergence(b_DG, DivB_C)

    use ModPotentialField, ONLY: nR, nTheta, nPhi, Radius_I, dRadius_I, &
         dTheta_I, dPhi_I, SinTheta_I, RadiusNode_I, SinThetaNode_I

    real, intent(in) :: b_DG(3,nR+1,nTheta+1,nPhi+1)
    real, intent(out):: DivB_C(nR,nTheta,nPhi)

    integer:: iR, iTheta, iPhi
    !--------------------------------------------------------------------------
    do iPhi = 1, nPhi
       do iTheta = 1, nTheta
          do iR = 1, nR
             DivB_C(iR,iTheta,iPhi) = &
                  ( RadiusNode_I(iR+1)**2*b_DG(1,iR+1,iTheta,iPhi)   &
                  - RadiusNode_I(iR)**2  *b_DG(1,iR  ,iTheta,iPhi) ) &
                  / (Radius_I(iR)**2 *dRadius_I(iR)) &
                  + &
                  ( SinThetaNode_I(iTheta+1)*b_DG(2,iR,iTheta+1,iPhi)   &
                  - SinThetaNode_I(iTheta)  *b_DG(2,iR,iTheta  ,iPhi) ) &
                  / (Radius_I(iR)*SinTheta_I(iTheta)*dTheta_I(iTheta)) &
                  + &
                  ( b_DG(3,iR,iTheta,iPhi+1) &
                  - b_DG(3,iR,iTheta,iPhi) ) &
                  / (Radius_I(iR)*SinTheta_I(iTheta)*dPhi_I(iPhi))
          end do
       end do
    end do

  end subroutine get_divergence

end module ModB0Matvec

!==============================================================================

program potential_field

  ! Solve 3D potential field with given Br at inner boundary,
  ! radial field at outer boundary.

  use ModPotentialField
  use ModB0Matvec, ONLY: get_gradient, get_divergence, matvec
  use ModLinearSolver, ONLY: gmres, bicgstab
  use ModPlotFile, ONLY: save_plot_file

  implicit none

  integer :: nKrylov=400, nIter=10000
  real    :: Tolerance = 0.0001
  integer :: n, iError, iTheta
  !--------------------------------------------------------------------------

  if(DoReadMagnetogram) call read_magnetogram

  call init_potential_field

  if(.not.DoReadMagnetogram)then
     allocate(Br_II(nTheta,nPhi))
     do iTheta = 1, nTheta
        Br_II(iTheta,:) = cos(Theta_I(iTheta))
     end do
  end if

  n = nR*nTheta*nPhi
  UseBr = .true.
  call matvec(Potential_C, Rhs_C, n)
  Rhs_C = -Rhs_C

  UseBr = .false.
  !    call gmres(matvec, Rhs_C, Potential_C, .false., n, &
  !         nKrylov, Tolerance, 'rel', nIter, iError, DoTest=.true.)

  call bicgstab(matvec, Rhs_C, Potential_C, .false., n, &
       Tolerance, 'rel', nIter, iError, DoTest=.true.)

  UseBr = .true.
  write(*,*)'nIter, Tolerance, iError=', nIter, Tolerance, iError

  PlotVar_VG = 0.0

  call get_gradient(Potential_C, B0_DG)
  PlotVar_VG(1:3,:,:,:) = B0_DG

  ! report maximum divb
  call get_divergence(B0_DG, DivB_C)
  write(*,*) 'max(abs(divb)) = ', maxval(abs(DivB_C))
  PlotVar_VG(4,1:nR,1:nTheta,1:nPhi) = DivB_C

  PlotVar_VG(5,1:nR,1:nTheta,1:nPhi) = Rhs_C

  PlotVar_VG(6,1:nR,1:nTheta,1:nPhi) = Potential_C

  call save_plot_file('potentialtest.out', &
       StringHeaderIn='potential field', &
       NameVarIn='r theta phi br btheta bphi divb rhs pot', &
       Coord1In_I=RadiusNode_I, &
       Coord2In_I=ThetaNode_I, &
       Coord3In_I=PhiNode_I, &
       VarIn_VIII=PlotVar_VG)

  deallocate(PlotVar_VG)

  call save_potential_field

end program potential_field
!==============================================================================
subroutine CON_stop(String)

  character(len=*), intent(in):: String

  write(*,*) 'ERROR:', String
  stop

end subroutine CON_stop
!==============================================================================
