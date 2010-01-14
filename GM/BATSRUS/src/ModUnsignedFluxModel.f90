!^CFG COPYRIGHT UM
!==============================================================================
module ModUnsignedFluxModel

  ! based on Abett (2007)

  implicit none
  save

  private ! except

  ! Public methods
  public :: get_coronal_heat_factor
  public :: get_coronal_heating

  !Bill Abbet model, if .true.
  logical, public :: UseUnsignedFluxModel = .false.

  ! Normalized value of Heating constant
  real, public :: HeatFactor = 0.0

  ! Cgs value of total power input from coronal heating
  real, public :: TotalCoronalHeatingCgs = 1.0e+28

  ! Exponential Scale height to truncate heating function
  real, public :: DecayLength = 1.0


contains

  !============================================================================

  subroutine get_coronal_heat_factor

    use ModAdvance,     ONLY: State_VGB, B0_DGB, Bx_, Bz_
    use ModGeometry,    ONLY: vInv_CB, true_BLK, true_cell
    use ModMagnetogram, ONLY: nTheta, nPhi, dSinTheta, dPhi, &
         get_magnetogram_field
    use ModMain,        ONLY: nI, nJ, nK, nBlock, UnusedBLK
    use ModMpi,         ONLY: MPI_REAL, MPI_SUM
    use ModNumConst,    ONLY: cHalfPi
    use ModPhysics,     ONLY: Si2No_V, No2Si_V, UnitX_, UnitB_, UnitT_, &
         UnitEnergyDens_, rBody
    use ModProcMH,      ONLY: nProc, iComm

    integer :: i, j, k, iBlock
    integer :: iTheta, iPhi, iError
    real :: UnsignedFlux, UnsignedFluxCgs, dAreaCgs
    real :: HeatFunction, HeatFunctionVolume, HeatFunctionVolumePe
    real :: x, y, z, Theta, Phi, SinTheta, CosTheta, SinPhi, CosPhi
    real :: FullB_D(3), B0_D(3), BrSi, BrCgs, SumUnsignedBrCgs

    real, save :: TotalCoronalHeating
    logical, save :: DoFirst = .true.

    real, parameter :: HeatExponent = 1.1488, HeatCoef = 89.4
    !--------------------------------------------------------------------------

    if(DoFirst)then

       ! uniform cell area on sphere
       dAreaCgs = rBody**2*dSinTheta*dPhi*No2Si_V(UnitX_)**2*1e4
       SumUnsignedBrCgs = 0.0

       do iTheta = 0, nTheta
          Theta = cHalfPi - asin((real(iTheta) + 0.5)*dSinTheta - 1.0)
          SinTheta = sin(Theta)
          CosTheta = cos(Theta)
          do iPhi = 1, nPhi
             Phi=(real(iPhi)-0.5)*dPhi
             SinPhi = sin(Phi)
             CosPhi = cos(Phi)

             x = rBody*SinTheta*CosPhi
             y = rBody*SinTheta*SinPhi
             z = rBody*CosTheta

             call get_magnetogram_field(x, y, z, B0_D)
             BrSi = (x*B0_D(1) + y*B0_D(2) + z*B0_D(3))/rBody
             BrCgs = BrSi*1e4
             SumUnsignedBrCgs = SumUnsignedBrCgs + abs(BrCgs)
          end do
       end do

       UnsignedFluxCgs = SumUnsignedBrCgs*dAreaCgs

       TotalCoronalHeatingCgs = HeatCoef*UnsignedFluxCgs**HeatExponent

       TotalCoronalHeating = TotalCoronalHeatingCgs*1e-7 &
            *Si2No_V(UnitEnergyDens_)*Si2No_V(UnitX_)**3/Si2No_V(UnitT_)

       DoFirst = .false.
    end if

    HeatFunctionVolume = 0
    do iBlock = 1, nBlock
       if(unusedBLK(iBlock)) CYCLE

       if(true_BLK(iBlock)) then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             call get_heat_function(i, j, k, iBlock, HeatFunction)
             HeatFunctionVolume = HeatFunctionVolume &
                  + HeatFunction/vInv_CB(i,j,k,iBlock)
          end do; end do; end do
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(true_cell(i,j,k,iBlock))then
                call get_heat_function(i, j, k, iBlock, HeatFunction)
                HeatFunctionVolume = HeatFunctionVolume &
                     + HeatFunction/vInv_CB(i,j,k,iBlock)
             end if
          end do; end do; end do
       end if
    end do

    if(nProc>1)then
       HeatFunctionVolumePe = HeatFunctionVolume
       call MPI_allreduce(HeatFunctionVolumePe, HeatFunctionVolume, 1, &
            MPI_REAL, MPI_SUM, iComm, iError)
    end if

    HeatFactor = TotalCoronalHeating/HeatFunctionVolume

  end subroutine get_coronal_heat_factor

  !============================================================================

  subroutine get_coronal_heating(i, j, k, iBlock, CoronalHeating)

    integer, intent(in) :: i, j, k, iBlock
    real, intent(out) :: CoronalHeating

    real :: HeatFunction
    !--------------------------------------------------------------------------

    call get_heat_function(i, j, k, iBlock, HeatFunction)

    CoronalHeating = HeatFactor*HeatFunction

  end subroutine get_coronal_heating

  !============================================================================

  subroutine get_heat_function(i, j, k, iBlock, HeatFunction)

    use ModAdvance, ONLY: State_VGB, B0_DGB, Bx_, Bz_
    use ModGeometry, ONLY: r_BLK

    integer, intent(in) :: i, j, k, iBlock
    real, intent(out) :: HeatFunction

    real :: Bmagnitude, B_D(3)
    !--------------------------------------------------------------------------

    B_D = B0_DGB(:,i,j,k,iBlock) + State_VGB(Bx_:Bz_,i,j,k,iBlock)
    Bmagnitude = sqrt(sum(B_D**2))

    HeatFunction = Bmagnitude*exp(-(r_BLK(i,j,k,iBlock)-1.0)/DecayLength)

  end subroutine get_heat_function

end module ModUnsignedFluxModel
