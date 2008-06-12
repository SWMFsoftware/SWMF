!^CFG COPYRIGHT UM
!==============================================================================
module EEE_ModShearFlow
  use EEE_ModCommonVariables
  implicit none
  save

  private

  public :: set_parameters_shearflow,get_shearflow

  logical, public :: UseShearFlow=.false.
  real :: FlowAmplitude, Longitude, Latitude

  real :: xFlow_D(3)

contains

  !============================================================================

  subroutine set_parameters_shearflow
    use ModReadParam, ONLY: read_var
    implicit none
    !--------------------------------------------------------------------------
    call read_var('UseShearFlow',  UseShearFlow)
    call read_var('FlowAmplitude', FlowAmplitude)
    call read_var('Longitude',     Longitude)
    call read_var('Latitude',      Latitude)

  end subroutine set_parameters_shearflow

  !============================================================================

  subroutine get_shearflow(x_D,Time,U_D)

    !DESCRIPTION:
    ! Boundary shear flow that conserves the radial magnetic flux

    use EEE_ModCommonVariables
    use EEE_ModGetB0,   ONLY: EEE_get_B0
    use ModMagnetogram, ONLY: get_magnetogram_field
    use ModNumConst,    ONLY: cTolerance, cDegToRad
    implicit none

    real, intent(in)  :: x_D(3), Time
    real, intent(out) :: U_D(3)

    real :: R, Theta, Phi, Xy, SinTheta, CosTheta, SinPhi, CosPhi, UnitR_D(3)
    real :: dTheta, dPhi
    real :: B_D(3), B0_D(3), FullBn, FullBnL, FullBnR
    real :: ShearProfileL, ShearProfileR
    real :: UTheta, UPhi
    logical, save :: DoFirst=.true.
    !--------------------------------------------------------------------------
    call CON_stop('The shear flow implementation is not yet finished')

    if(DoFirst)then
       DoFirst = .false.
       xFlow_D(1) = cos(Latitude*cDegToRad)*cos(Longitude*cDegToRad)
       xFlow_D(2) = cos(Latitude*cDegToRad)*sin(Longitude*cDegToRad)
       xFlow_D(3) = sin(Latitude*cDegToRad)
    end if

    R        = sqrt(sum(x_D**2))
    Theta    = acos(x_D(3)/R)
    Phi      = atan2(x_D(2),x_D(1))
    Xy       = max(sqrt(x_D(1)**2+x_D(2)**2),cTolerance)
    SinTheta = Xy/R
    CosTheta = x_D(3)/R
    SinPhi   = x_D(2)/Xy
    CosPhi   = x_D(1)/Xy

    call get_magnetogram_field(x_D(1),x_D(2),x_D(3),B0_D)
    call EEE_get_B0(x_D,B_D)
    B0_D = (B0_D + B_D)*Si2No_V(UnitB_)

    UnitR_D = x_D/R
    FullBn = dot_product(UnitR_D,B0_D)

    dTheta = cTiny
    dPhi   = cTiny

    if(abs(FullBn)<cTolerance)then
       ! No flow at polarity inversion lines
       U_D = 0.0
       return
    else
       ShearProfileL = shear_profile(R,Theta,Phi-0.5*dPhi,Time,FullBnL)
       ShearProfileR = shear_profile(R,Theta,Phi+0.5*dPhi,Time,FullBnR)
       if(FullBnL*FullBnR<0.0 .and. abs(FullBnL)>cTolerance &
            .and. abs(FullBnR)>cTolerance)then
          UTheta = 0.0
       else
          UTheta = 1.0/FullBn &
               *(ShearProfileR-ShearProfileL)/(R*SinTheta*dPhi)
       end if

       ShearProfileL = shear_profile(R,Theta-0.5*dTheta,Phi,Time,FullBnL)
       ShearProfileR = shear_profile(R,Theta+0.5*dTheta,Phi,Time,FullBnR)
       if(FullBnL*FullBnR<0.0 .and. abs(FullBnL)>cTolerance &
            .and. abs(FullBnR)>cTolerance)then
          UPhi = 0.0
       else
          UPhi = -1.0/FullBn &
               *(ShearProfileR-ShearProfileL)/(R*dTheta)
       end if
    end if

    U_D(1) = UTheta*CosTheta*CosPhi - UPhi*SinPhi
    U_D(2) = UTheta*CosTheta*SinPhi + UPhi*CosPhi
    U_D(3) =-UTheta*SinTheta

    U_D = U_D*No2Si_V(UnitU_)

  end subroutine get_shearflow

  !============================================================================

  real function shear_profile(R,Theta,Phi,Time,FullBn)
    use EEE_ModCommonVariables
    use EEE_ModGetB0,   ONLY: EEE_get_B0
    use ModMagnetogram, ONLY: get_magnetogram_field
    implicit none

    real, intent(in)  :: R, Theta, Phi, Time
    real, intent(out) :: FullBn

    real :: Xy, x_D(3), UnitR_D(3)
    real :: B0_D(3), B_D(3)
    real :: Del_D(3), Mask, Ramp
    real :: MaxBnAR, MaxBnARDim, BnRatio, ShearAngle, Angle
    !--------------------------------------------------------------------------
    Xy = R*sin(Theta)
    x_D(1) = Xy*cos(Phi)
    x_D(2) = Xy*sin(Phi)
    x_D(3) = R*cos(Theta)

    call get_magnetogram_field(x_D(1),x_D(2),x_D(3),B0_D)
    call EEE_get_B0(x_D,B_D)
    B0_D = (B0_D + B_D)*Si2No_V(UnitB_)

    UnitR_D = x_D/R
    FullBn = dot_product(UnitR_D,B0_D)

    Del_D = xFlow_D - UnitR_D
    Angle = 2.0*asin(0.5*sqrt(dot_product(Del_D,Del_D)))
    Mask = exp(-(Angle/ShearAngle)**2)

    Ramp = 1.0

    MaxBnAR = MaxBnARDim*Io2No_V(UnitB_)
    BnRatio = FullBn/MaxBnAR

    shear_profile = FlowAmplitude*FullBn**3*exp(-BnRatio**2)*Mask*Ramp

  end function shear_profile

end module EEE_ModShearFlow
