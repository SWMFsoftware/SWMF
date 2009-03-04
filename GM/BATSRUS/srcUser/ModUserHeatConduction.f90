!^CFG COPYRIGHT UM
!==============================================================================
module ModUser
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_init_session,               &
       IMPLEMENTED3 => user_set_ics,                    &
       IMPLEMENTED4 => user_update_states,              &
       IMPLEMENTED5 => user_set_plot_var,               &
       IMPLEMENTED6 => user_material_properties

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'heat conduction'

  integer :: iHeatConductionTest
  real :: HeatConductionCoef, AmplitudeTemperature, T0

contains

  !============================================================================

  subroutine user_init_session

    use ModGeometry, ONLY: x1, x2, XyzMin_D, XyzMax_D
    use ModMain,     ONLY: nI
    use ModParallel, ONLY: proc_dims

    character(len=*), parameter :: NameSub = 'user_init_session'
    !--------------------------------------------------------------------------

    if(x2 /= -x1) &
       call stop_mpi(NameSub//' : xMin should be -xMax in the #GRID command')

    ! Shift grid to the right by half a grid cell size,
    ! so that x=0 is at a cell center
    x1 = x1 - x1/(2*nI*proc_dims(1))
    x2 = x2 + x2/(2*nI*proc_dims(2))

    XyzMin_D(1) = x1
    XyzMax_D(1) = x2

    select case(iHeatConductionTest)
    case(1)
       HeatConductionCoef = 0.1
       AmplitudeTemperature = 10.0
       T0 = 10.0
    case(2)
       HeatConductionCoef = 0.1
       AmplitudeTemperature = 10.0
       T0 = 3.0
    end select

  end subroutine user_init_session

  !============================================================================

  subroutine user_update_states(iStage, iBlock)

    integer, intent(in) :: iStage, iBlock

    character(len=*), parameter :: NameSub = 'user_update_states'
    !--------------------------------------------------------------------------

    ! No call to update_states_MHD to nullify the effect of the hydro solver
    ! call update_states_MHD(iStage,iBlock)

  end subroutine user_update_states

  !============================================================================

  subroutine user_read_inputs

    use ModIO,          ONLY: write_prefix, write_myname, iUnitOut
    use ModMain,        ONLY: lVerbose
    use ModProcMH,      ONLY: iProc
    use ModReadParam,   ONLY: read_line, read_command, read_var

    character (len=100) :: NameCommand
    !--------------------------------------------------------------------------
    if(iProc == 0 .and. lVerbose > 0)then
       call write_prefix;
       write(iUnitOut,*)'User read_input starts'
    endif

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)
       case("#HEATCONDUCTIONTEST")
          call read_var('iHeatConductionTest',iHeatConductionTest)

       case('#USERINPUTEND')
          if(iProc == 0 .and. lVerbose > 0)then
             call write_prefix;
             write(iUnitOut,*)'User read_input ends'
          endif
          EXIT

       case default
          if(iProc == 0) then
             call write_myname; write(*,*) &
                  'ERROR: Invalid user defined #COMMAND in user_read_inputs. '
             write(*,*) '--Check user_read_inputs for errors'
             write(*,*) '--Check to make sure a #USERINPUTEND command was used'
             write(*,*) '  *Unrecognized command was: '//NameCommand
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if
       end select
    end do

  end subroutine user_read_inputs

  !============================================================================

  subroutine user_set_ics

    use ModAdvance,    ONLY: State_VGB
    use ModGeometry,   ONLY: x_Blk, y_Blk, Dx_Blk, y2
    use ModMain,       ONLY: GlobalBlk, nI, nJ, nK, x_, Time_Simulation
    use ModNumConst,   ONLY: cPi
    use ModPhysics,    ONLY: gm1
    use ModVarIndexes, ONLY: Rho_, RhoUx_, p_

    integer :: iBlock, i, j, k
    real :: x, Dx, Temperature, Spread
    real :: r, Lambda

    character(len=*), parameter :: NameSub = "user_set_ics"
    !--------------------------------------------------------------------------
    iBlock = GlobalBlk

    select case(iHeatConductionTest)
    case(1)
       ! temperature spike at x=0, otherwise set the temperature to zero
       do i = -1, nI+2
          if(Time_Simulation == 0.0)then
             if(x_Blk(i,1,1,iBlock) > -0.5*Dx_Blk(iBlock) .and. &
                  x_Blk(i,1,1,iBlock) < 0.5*Dx_Blk(iBlock)) then
                Temperature = T0 + AmplitudeTemperature/Dx_Blk(iBlock)
             else
                Temperature = T0
             end if
          else
             Spread = 4.0*HeatConductionCoef*Time_Simulation
             Temperature = T0 + AmplitudeTemperature/(sqrt(cPi*Spread)) &
                  *exp(-x_Blk(i,1,1,iBlock)**2/Spread)
          end if
          do k = -1, nK+2; do j = -1, nJ+2
             State_VGB(Rho_,i,j,k,iBlock) = 1.0
             State_VGB(RhoUx_:p_-1,i,j,k,iBlock) = 0.0
             State_VGB(p_,i,j,k,iBlock) = Temperature
          end do; end do
       end do
    case(2)
       Lambda = -(3.831705970/y2)**2
       Spread = 4.0*HeatConductionCoef*Time_Simulation
       do j = -1, nJ+2; do i = -1, nI+2
          r = abs(y_BLK(i,j,1,iBlock))
          Temperature = T0 + AmplitudeTemperature &
               *exp(Lambda*HeatConductionCoef*Time_Simulation) &
               *bessj0(sqrt(-Lambda)*r) &
               /(sqrt(cPi*Spread))*exp(-x_Blk(i,j,1,iBlock)**2/Spread)
          do k = -1, nK+2
             State_VGB(Rho_,i,j,k,iBlock) = 1.0
             State_VGB(RhoUx_:p_-1,i,j,k,iBlock) = 0.0
             State_VGB(p_,i,j,k,iBlock) = Temperature
          end do
       end do; end do
    end select

  end subroutine user_set_ics

  !============================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModAdvance,    ONLY: State_VGB
    use ModGeometry,   ONLY: x_Blk, y_Blk, Dx_Blk, y2
    use ModMain,       ONLY: nI, nJ, nK, Time_Simulation
    use ModNumConst,   ONLY: cPi
    use ModVarIndexes, ONLY: p_, Rho_

    integer,          intent(in)   :: iBlock
    character(len=*), intent(in)   :: NameVar
    logical,          intent(in)   :: IsDimensional
    real,             intent(out)  :: PlotVar_G(-1:nI+2, -1:nJ+2, -1:nK+2)
    real,             intent(out)  :: PlotVarBody
    logical,          intent(out)  :: UsePlotVarBody
    character(len=*), intent(inout):: NameTecVar
    character(len=*), intent(inout):: NameTecUnit
    character(len=*), intent(inout):: NameIdlUnit
    logical,          intent(out)  :: IsFound

    real :: Spread, Lambda, r, Temperature
    integer :: i, j, k

    character (len=*), parameter :: NameSub = 'user_set_plot_var'
    !--------------------------------------------------------------------------

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0

    IsFound = .true.
    select case(NameVar)
    case('te')
       PlotVar_G = State_VGB(p_,:,:,:,iBlock)/State_VGB(Rho_,:,:,:,iBlock)
    case('te0')
       select case(iHeatConductionTest)
       case(1)
          if(Time_Simulation == 0.0)then
             PlotVar_G = State_VGB(p_,:,:,:,iBlock) &
                  /State_VGB(Rho_,:,:,:,iBlock)
          else
             Spread = 4.0*HeatConductionCoef*Time_Simulation
             PlotVar_G = T0 + AmplitudeTemperature/(sqrt(cPi*Spread)) &
                  *exp(-x_Blk(:,:,:,iBlock)**2/Spread)
          end if
       case(2)
          Lambda = -(3.831705970/y2)**2
          Spread = 4.0*HeatConductionCoef*Time_Simulation
          do j = -1, nJ+2; do i = -1, nI+2
             r = abs(y_BLK(i,j,1,iBlock))
             Temperature = T0 + AmplitudeTemperature &
                  *exp(Lambda*HeatConductionCoef*Time_Simulation) &
                  *bessj0(sqrt(-Lambda)*r) &
                  /(sqrt(cPi*Spread))*exp(-x_Blk(i,j,1,iBlock)**2/Spread)
             do k = -1, nK+2
                PlotVar_G(i,j,k) = Temperature
             end do
          end do; end do
       end select
    case default
       IsFound = .false.
    end select

  end subroutine user_set_plot_var

  !============================================================================

  subroutine user_material_properties(State_V, EinternalSiIn, &
       TeSiIn, EinternalSiOut, TeSiOut, PressureSiOut, CvSiOut, &
       AbsorptionOpacitySiOut, RosselandMeanOpacitySiOut, &
       HeatConductionCoefSiOut)

    ! The State_V vector is in normalized units

    use ModPhysics,    ONLY: gm1, inv_gm1, No2Si_V, Si2No_V, &
         UnitRho_, UnitP_, UnitEnergyDens_, UnitTemperature_, &
         UnitX_, UnitT_, UnitU_
    use ModVarIndexes, ONLY: nVar, Rho_, p_

    real, intent(in) :: State_V(nVar)
    real, optional, intent(in)  :: EinternalSiIn             ! [J/m^3]
    real, optional, intent(in)  :: TeSiIn                    ! [K]
    real, optional, intent(out) :: EinternalSiOut            ! [J/m^3]
    real, optional, intent(out) :: TeSiOut                   ! [K]
    real, optional, intent(out) :: AbsorptionOpacitySiOut    ! [1/m]
    real, optional, intent(out) :: RosselandMeanOpacitySiOut ! [1/m]
    real, optional, intent(out) :: CvSiOut                   ! [J/(K*m^3)]
    real, optional, intent(out) :: PressureSiOut             ! [Pa]
    real, optional, intent(out) :: HeatConductionCoefSiOut   ! [Jm^2/(Ks)]

    real :: Rho, Pressure, Tmat
    real :: RhoSi, pSi, TmatSi

    character(len=*), parameter :: NameSub = 'user_material_properties'
    !--------------------------------------------------------------------------

    Rho = State_V(Rho_)
    RhoSi = Rho*No2Si_V(Rho_)

    if(present(EinternalSiIn))then
       pSi = EinternalSiIn*gm1
       Pressure = pSi*Si2No_V(UnitP_)
       Tmat = Pressure/Rho
       TmatSi = Tmat*No2Si_V(UnitTemperature_)
    elseif(present(TeSiIn))then
       TmatSi = TeSiIn
       Tmat = TmatSi*Si2No_V(UnitTemperature_)
       Pressure = Rho*Tmat
       pSi = Pressure*No2Si_V(UnitP_)
    else
       Pressure = State_V(p_)
       pSi = Pressure*No2Si_V(UnitP_)
       Tmat = Pressure/Rho
       TmatSi = Tmat*No2Si_V(UnitTemperature_)
    end if

    if(present(EinternalSiOut)) EinternalSiOut = pSi*inv_gm1
    if(present(TeSiOut)) TeSiOut = TmatSi
    if(present(PressureSiOut)) PressureSiOut = pSi

    if(present(CvSiOut)) CvSiOut = inv_gm1*Rho &
         *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_)

    if(present(HeatConductionCoefSiOut)) &
       HeatConductionCoefSiOut = HeatConductionCoef &
            *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_) &
            *No2Si_V(UnitU_)*No2Si_V(UnitX_)

    if(present(AbsorptionOpacitySiOut)) AbsorptionOpacitySiOut = 0.0
    if(present(RosselandMeanOpacitySiOut)) RosselandMeanOpacitySiOut = 0.0

  end subroutine user_material_properties

  !============================================================================

  real function bessj0(x)

    ! bessj0 from numerical recipes (changed to free format)

    real :: x, z, ax, xx

    real y,p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6
    data p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4, &
         -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1, &
         .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
    data r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,651619640.7d0, &
         -11214424.18d0,77392.33017d0,-184.9052456d0/, &
         s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0, &
         9494680.718d0,59272.64853d0,267.8532712d0,1.d0/

    !--------------------------------------------------------------------------

    if(abs(x).lt.8.)then
       y=x**2
       bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))))) &
            /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
    else
       ax=abs(x)
       z=8./ax
       y=z**2
       xx=ax-.785398164
       bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y &
            *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
    endif

  end function bessj0

end module ModUser
