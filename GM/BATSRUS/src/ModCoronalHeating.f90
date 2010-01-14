module ModCoronalHeating
  use ModUnsignedFluxModel
  use ModNonWKBHeating
  use ModMain,      ONLY: nBLK, nI, nJ, nK
  use ModReadParam, ONLY: lStringLine
  implicit none
  SAVE

  logical,public:: UseCoronalHeating = .false.
  character(len=lStringLine) :: NameModel, TypeCoronalHeating

  ! Variables and parameters for various heating models
  !\
  ! Exponential:
  !/  

  ! quantitative parameters for exponential heating model

  real :: HeatingAmplitude, HeatingAmplitudeCgs = 6.07e-7
  real :: DecayLengthExp = 0.7  ! in Solar Radii units
  logical :: UseExponentialHeating = .false.

  ! parameters for high-B transition (in Gauss)
  ! idea is to grossly approx 'Active Region' heating values
  logical :: UseArComponent = .false.
  real :: ArHeatB0 = 30.0
  real :: DeltaArHeatB0 = 5.0
  real :: ArHeatFactorCgs = 4.03E-05  ! cgs energy density = [ergs cm-3 s-1]
  
  ! open closed heat is if you want to have different heating
  ! between open and closed field lines. The routines for the WSA
  ! model in the ExpansionFactors module essentially do this, so will
  ! need to call them
  logical :: DoOpenClosedHeat = .false.
  real :: WsaT0 = 3.50
  
  !\
  ! Abbet's model
  !/
  
  
  ! Normalization constant for Abbet Model
  real :: HeatNormalization = 1.0


  
  ! long scale height heating (Ch = Coronal Hole)
  logical :: DoChHeat = .false.
  real :: HeatChCgs = 5.0e-7
  real :: DecayLengthCh = 0.7
  
contains
  subroutine read_corona_heating
    use ModReadParam,   ONLY: read_var
    !----------------------------------
    call read_var('TypeCoronalHeating', TypeCoronalHeating)

    !Initialize logicals
    UseCoronalHeating = .true.

    UseUnsignedFluxModel = .false.
    UseCranmerHeating      = .false.
    UseExponentialHeating= .false.


    select case(TypeCoronalHeating)
    case('F','none')
       UseCoronalHeating = .false.
       
    case('exponential')
       UseExponentialHeating = .true.
       call read_var('DecayLengthExp', DecayLengthExp)
       call read_var('HeatingAmplitudeCgs', HeatingAmplitudeCgs)

    case('unsignedflux','Abbet')
       UseUnsignedFluxModel = .true.
       call read_var('DecayLength', DecayLength)
       call read_var('HeatNormalization', HeatNormalization)
    case('Cranmer','NonWKB')
       UseCranmerHeating     = .true.
    case default
       call stop_mpi('Read_corona_heating: unknown TypeCoronalHeating = ' &
            //TypeCoronalHeating)
    end select
  end subroutine read_corona_heating
  !===========================
  subroutine read_active_region_heating
    use ModReadParam,   ONLY: read_var
    !---------------------------------
    call read_var('UseArComponent', UseArComponent)
    if(UseArComponent) then
       call read_var('ArHeatFactorCgs', ArHeatFactorCgs)
       call read_var('ArHeatB0', ArHeatB0)
       call read_var('DeltaArHeatB0', DeltaArHeatB0)
    endif
  end subroutine read_active_region_heating
  !===========================
  subroutine read_longscale_heating
    use ModReadParam,   ONLY: read_var
    use ModReadParam,   ONLY: read_var
    !---------------------------------
    call read_var('DoChHeat', DoChHeat)
    if(.not.DoChHeat)return
    call read_var('HeatChCgs', HeatChCgs)
    call read_var('DecayLengthCh', DecayLengthCh)
  end subroutine read_longscale_heating
  !===========================
  
  subroutine get_cell_heating(i, j, k, iBlock, CoronalHeating)

    use ModGeometry,       ONLY: r_BLK, x_BLK, y_BLK, z_BLK
    use ModPhysics,        ONLY: Si2No_V, No2Si_V, UnitEnergyDens_, UnitT_, &
         No2Io_V, UnitB_,UnitRho_,UnitX_,UnitU_
    use ModExpansionFactors, ONLY: UMin
    use ModMain,       ONLY: x_, y_, z_
    use ModVarIndexes, ONLY: Bx_, By_, Bz_,Rho_,RhoUx_,RhoUz_
    use ModAdvance,    ONLY: State_VGB, B0_DGB

    integer, intent(in) :: i, j, k, iBlock
    real, intent(out) :: CoronalHeating
    
    real :: HeatCh

    ! parameters for open/closed. This uses WSA model to determine if 
    ! cell is in an 'open' or 'closed' field region.
    !
    ! ** NOTE ** WSA does field line tracing on an auxiliary grid.
    ! should really be using the computational domain, but global 
    ! feild line tracing for this purpose is not easily implemented
    real :: Ufinal
    real :: UminIfOpen = 290.0
    real :: UmaxIfOpen = 300.0
    real :: Weight
    
    ! local variables for ArHeating (Active Region Heating)
    real :: FractionB, Bcell

    real :: RhoSI, RSI, UMagSI, BMagSI, QHeatSI
    !--------------------------------------------------------------------------
    
    if(UseUnsignedFluxModel)then
       
       call get_coronal_heating(i, j, k, iBlock, CoronalHeating)
       CoronalHeating = CoronalHeating * HeatNormalization
       
    elseif(UseExponentialHeating)then
    
       CoronalHeating = HeatingAmplitude &
            *exp(- max(r_BLK(i,j,k,iBlock) - 1.0, 0.0) / DecayLengthExp)
            
    elseif(UseCranmerHeating)then

       RSI   = r_BLK(i,j,k,iBlock) * No2Si_V(UnitX_)
       RhoSI = State_VGB(Rho_,i,j,k,iBlock) * No2Si_V(UnitRho_)
       UMagSI= sqrt( sum( State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock)**2 ) )/&
                          State_VGB(Rho_,i,j,k,iBlock) * No2Si_V(UnitU_)
       BmagSI= sqrt( sum( State_VGB(Bx_   :   Bz_,i,j,k,iBlock)**2 ) ) &
                                                       * No2Si_V(UnitB_)
       call Cranmer_heating_function(&
            RSI=RSI,      &
            RhoSI = RhoSI,&
            UMagSI=UMagSI,&
            BMagSI=BMagSI,&
            QHeatSI=QHeatSI)
       CoronalHeating = QHeatSI *  Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)
    else
       CoronalHeating = 0.0
    end if
  
    if(DoOpenClosedHeat)then
       ! If field is less than 1.05 times the minimum speed, mark as closed
       ! Interpolate between 1.05 and 1.10 for smoothness
       UminIfOpen = UMin*1.05
       UmaxIfOpen = UMin*1.1
       call get_bernoulli_integral(x_BLK( i, j, k, iBlock)/&
            R_BLK( i, j, k, iBlock),&
            y_BLK( i, j, k, iBlock)/R_BLK( i, j, k, iBlock),&
            z_BLK( i, j, k, iBlock)/R_BLK( i, j, k, iBlock), UFinal)
       if(UFinal <= UminIfOpen) then
          Weight = 0.0
       else if (UFinal >= UmaxIfOpen) then
          Weight = 1.0
       else 
          Weight = (UFinal - UminIfOpen)/(UmaxIfOpen - UminIfOpen)
       end if
       
       CoronalHeating = (1.0 - Weight) * CoronalHeating
    end if

    if(DoChHeat) then
       HeatCh = HeatChCgs * 0.1 * Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)
       CoronalHeating = CoronalHeating + HeatCh &
            *exp(- max(r_BLK(i,j,k,iBlock) - 1.0, 0.0) / DecayLengthCh)
    end if

    if(UseExponentialHeating.and.UseArComponent) then

       Bcell = No2Io_V(UnitB_) * sqrt(&
            ( B0_DGB(x_,i,j,k,iBlock) + State_VGB(Bx_,i,j,k,iBlock))**2 &
            +(B0_DGB(y_,i,j,k,iBlock) + State_VGB(By_,i,j,k,iBlock))**2 &
            +(B0_DGB(z_,i,j,k,iBlock) + State_VGB(Bz_,i,j,k,iBlock))**2)
          
       FractionB = 0.5*(1.0+tanh((Bcell - ArHeatB0)/DeltaArHeatB0))
       CoronalHeating = max(CoronalHeating, & 
                 FractionB * ArHeatFactorCgs * Bcell &
                 * 0.1 * Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_))

    endif
                       
  end subroutine get_cell_heating

  !============================================================================
  
end module ModCoronalHeating
