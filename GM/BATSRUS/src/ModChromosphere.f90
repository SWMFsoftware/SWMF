!This code is a copyright protected software (c) 2002- University of Michigan

!==============================================================================
module ModChromosphere
  !Here all parameters relating to chromosphere and tansition region are
  !collected

  use ModSize, ONLY: nI, nJ, nK
  implicit none
  save
  
  !The use of short-scale exponential heat function with 
  !the decay length = (30 m/K)*TeCromosphere SI
  logical:: UseChromosphereHeating    = .false. 
 

  real   :: NumberDensChromosphereCgs = 1.0e+12 ![cm^{-3}  
  real   :: TeChromosphereSi = 1.0e4            ![K]

  !\
  ! TRANSITION REGION
  !/
  ! 
  
  logical :: DoExtendTransitionRegion = .false.

  !The following variables are meaningful if
  !DoExtendTransitionRegion=.true.

  real :: TeModSi = 3.0E+5                !K
  real :: DeltaTeModSi = 1E+4             !K

  !The following variable is meaningful if 
  !DoExtendTransitionRegion = .false. . Al long as
  !the unextended transition region cannot be resolved
  !we set the 'corona base temperature' equal to the
  !temperature at the top of the transition region and 
  !use the integral ralationship across the transition
  !region to find the number density at the top of the
  !transition region

  real :: TeTransitionRegionTopSi = 4.0e+5 ![K]

  ! Electron temperature in K:
  real :: TeSi_C(nI,nJ,nK)

contains
  !================================

  subroutine read_chromosphere
    use ModReadParam, ONLY: read_var
    !-------------------------------
    call read_var('UseChromosphereHeating'   , UseChromosphereHeating)
    call read_var('NumberDensChromosphereCgs', NumberDensChromosphereCgs)
    call read_var('TeChromosphereSi',          TeChromosphereSi      )
  end subroutine read_chromosphere

  !============================================================================

  real function extension_factor(TeSi)
    real, intent(in) :: TeSi    !Dimensionless
    
    real :: FractionSpitzer
    !--------------------------------
    FractionSpitzer = 0.5*(1.0+tanh((TeSi-TeModSi)/DeltaTeModSi))
    
    extension_factor = FractionSpitzer + &
         (1.0 - FractionSpitzer)*(TeModSi/TeSi)**2.5
  end function extension_factor

  !============================================================================

  subroutine get_tesi_c(iBlock, TeSi_C)

    use ModAdvance,    ONLY: UseElectronPressure, UseIdealEos
    use ModAdvance,    ONLY: State_VGB, p_, Pe_, Rho_
    use ModSize,       ONLY: nI, nJ, nK
    use ModPhysics,    ONLY: No2Si_V, UnitTemperature_, &
         AverageIonCharge, PePerPtotal
    use ModMultifluid, ONLY: MassIon_I
    use ModUser,       ONLY: user_material_properties

    integer, intent(in)  :: iBlock
    real,    intent(out) :: TeSi_C(1:nI, 1:nJ, 1:nK)

    integer:: i, j, k
    !--------------------------------------------------------------------------
    if(UseIdealEos)then
       if(UseElectronPressure)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             TeSi_C(i,j,k) = State_VGB(Pe_,i,j,k,iBlock) &
                  /State_VGB(Rho_,i,j,k,iBlock)
          end do; end do; end do
          TeSi_C = TeSi_C * No2Si_V(UnitTemperature_ ) * &
               MassIon_I(1)/AverageIonCharge
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             TeSi_C(i,j,k) = State_VGB(p_,i,j,k,iBlock) &
                  /State_VGB(Rho_,i,j,k,iBlock)
          end do; end do; end do
          TeSi_C = TeSi_C * No2Si_V(UnitTemperature_ ) * &
               MassIon_I(1)/AverageIonCharge * PePerPtotal
       end if
    else
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          call user_material_properties( &
               State_VGB(:,i,j,k,iBlock), TeOut=TeSi_C(i,j,k))
       end do; end do; end do
    end if
  end subroutine get_tesi_c

end module ModChromosphere
