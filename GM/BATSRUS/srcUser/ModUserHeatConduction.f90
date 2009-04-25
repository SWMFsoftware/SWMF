!^CFG COPYRIGHT UM
!==============================================================================
module ModUser
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_init_session,               &
       IMPLEMENTED3 => user_set_ics,                    &
       IMPLEMENTED4 => user_set_outerbcs,               &
       IMPLEMENTED5 => user_update_states,              &
       IMPLEMENTED6 => user_set_plot_var,               &
       IMPLEMENTED7 => user_material_properties

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'heat conduction'

  character(len=20) :: TypeProblem
  real :: HeatConductionCoef, AmplitudeTemperature, T0

  character(len=100) :: NameRefFile
  integer :: nCellRef
  integer :: &
       iRhoRef  = 1, &
       iTmatRef = 2, &
       iUrRef   = 3, &
       nVarRef  = 3
  real, allocatable :: rRef_C(:), StateRef_VC(:,:)

contains

  !============================================================================

  subroutine user_init_session

    use ModIoUnit,  ONLY: UnitTmp_
    use ModMain,    ONLY: Time_Simulation
    use ModProcMH,  ONLY: iProc

    integer :: iCell, iError
    integer :: nStepRef, nDimRef, nParamRef, nVarRef
    real :: TimeRef, GammaRef
    character(len=500) :: StringHeaderRef, NameVarRef

    character(len=*), parameter :: NameSub = 'user_init_session'
    !--------------------------------------------------------------------------

    if(TypeProblem=='gaussian' .or. TypeProblem=='rz')then
       if(Time_Simulation <= 0.0)then
          if(iProc == 0) write(*,*) NameSub// &
               ' : starting simulation time should be larger than 0'
          call stop_mpi('reset time with #TIMESIMULATION')
       end if
    end if

    select case(TypeProblem)
    case('gaussian')
       HeatConductionCoef = 0.1
       AmplitudeTemperature = 10.0
       T0 = 10.0

    case('rz')
       HeatConductionCoef = 0.1
       AmplitudeTemperature = 10.0
       T0 = 3.0

    case('rmtv')
       NameRefFile = 'rmtv_initial.out'

       open(UnitTmp_, FILE=NameRefFile, STATUS="old", IOSTAT=iError)

       if(iError /= 0)call stop_mpi(NameSub // &
            " could not open reference file="//NameRefFile)

       read(UnitTmp_,'(a)') StringHeaderRef
       read(UnitTmp_,*) nStepRef, TimeRef, nDimRef, nParamRef, nVarRef
       read(UnitTmp_,*) nCellRef
       read(UnitTmp_,*) GammaRef
       read(UnitTmp_,'(a)') NameVarRef

       allocate( rRef_C(nCellRef), StateRef_VC(nVarRef,nCellRef) )

       do iCell = 1, nCellRef
          read(UnitTmp_,*) rRef_C(iCell), StateRef_VC(:,iCell)
       end do

       close(UnitTmp_)

       rRef_C = rRef_C
       StateRef_VC(iRhoRef,:)  = StateRef_VC(iRhoRef,:)           ! g/cm^3
       StateRef_VC(iTmatRef,:) = StateRef_VC(iTmatRef,:)*1.0e-3   ! keV
       StateRef_VC(iUrRef,:)   = StateRef_VC(iUrRef,:)*1.0e-8     ! cm/sh

       do iCell = 1, nCellRef
          StateRef_VC(iTmatRef,iCell) = &
               max(StateRef_VC(iTmatRef,iCell),1.0e-12)
       end do

    case default
       call stop_mpi(NameSub//' : undefined problem type='//TypeProblem)
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
          call read_var('TypeProblem',TypeProblem)

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
    use ModPhysics,    ONLY: gm1, No2Io_V, UnitTemperature_
    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUy_, RhoUz_, ExtraEint_, p_

    integer :: iBlock, i, j, k, iCell
    real :: x, y, Dx, Temperature, Spread, Pressure
    real :: r, Lambda, Weight1, Weight2, Rho, Tmat, Ur, p, RhoU_D(2)

    character(len=*), parameter :: NameSub = "user_set_ics"
    !--------------------------------------------------------------------------
    iBlock = GlobalBlk

    select case(TypeProblem)
    case('gaussian')
       do i = -1, nI+2
          Spread = 4.0*HeatConductionCoef*Time_Simulation
          Temperature = T0 + AmplitudeTemperature/(sqrt(cPi*Spread)) &
               *exp(-x_Blk(i,1,1,iBlock)**2/Spread)
          do k = -1, nK+2; do j = -1, nJ+2
             State_VGB(Rho_,i,j,k,iBlock) = 1.0
             State_VGB(RhoUx_:p_-1,i,j,k,iBlock) = 0.0
             State_VGB(p_,i,j,k,iBlock) = Temperature
          end do; end do
       end do

    case('rz')
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

    case('rmtv')
       do j=1,nJ; do i=1,nI
          x = x_Blk(i,j,0,iBlock)
          y = y_Blk(i,j,0,iBlock)
          r = sqrt(x**2+y**2)

          do iCell = 1, nCellRef
             if(rRef_C(iCell) >= r) EXIT
          end do
          if(iCell == 1) call stop_mpi(NameSub // &
               " Reference solution does not cover the left boundary")

          if(iCell > nCellRef)then
             ! Cell is beyond the last point of Reference input: use last cell
             iCell   = nCellRef
             Weight1 = 0.0
             Weight2 = 1.0
          else
             ! Assign weights for linear interpolation between
             ! iCell-1 and iCell
             Weight1 = (rRef_C(iCell) - r) &
                  /    (rRef_C(iCell) - rRef_C(iCell-1))
             Weight2 = 1.0 - Weight1
          end if

          Rho  = ( Weight1*StateRef_VC(iRhoRef, iCell-1) &
               +   Weight2*StateRef_VC(iRhoRef, iCell) )
          Tmat = ( Weight1*StateRef_VC(iTmatRef, iCell-1) &
               +   Weight2*StateRef_VC(iTmatRef, iCell) )
          Ur   = ( Weight1*StateRef_VC(iUrRef, iCell-1) &
               +   Weight2*StateRef_VC(iUrRef, iCell) )

          p = Rho*Tmat

          RhoU_D(1) = Rho*Ur*x/r
          RhoU_D(2) = Rho*Ur*y/r

          do k=1,nk
             State_VGB(Rho_,i,j,k,iBlock) = Rho
             State_VGB(RhoUx_:RhoUy_,i,j,k,iBlock) = RhoU_D
             State_VGB(RhoUz_,i,j,k,iBlock) = 0.0
             State_VGB(ExtraEint_,i,j,k,iBlock) = 0.0
             State_VGB(p_,i,j,k,iBlock) = p
          end do
       end do; end do

    case default
       call stop_mpi(NameSub//' : undefined problem type='//TypeProblem)
    end select

  end subroutine user_set_ics

  !============================================================================

  subroutine user_set_outerbcs(iBlock,iSide, TypeBc, IsFound)

    use ModAdvance,    ONLY: State_VGB
    use ModGeometry,   ONLY: x_Blk, y_Blk
    use ModImplicit,   ONLY: StateSemi_VGB
    use ModMain,       ONLY: nI, nJ, nK
    use ModVarIndexes, ONLY: Rho_, RhoUx_, p_

    integer,          intent(in)  :: iBlock, iSide
    character(len=20),intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    integer :: i, j, k
    real :: r, Temperature

    character (len=*), parameter :: NameSub = 'user_set_outerbcs'
    !--------------------------------------------------------------------------

    if(.not. (iSide==2 .or. iSide==4) )then
       write(*,*) NameSub//' : user boundary not defined at iSide = ', iSide
       call stop_mpi(NameSub)
    end if

    select case(TypeBc)
    case('user')
       select case(iSide)
       case(2) ! z-direction in rz-geometry
          do k = -1, nK+2; do j = -1, nJ+2
             Temperature = State_VGB(p_,nI,j,k,iBlock) &
                  /State_VGB(Rho_,nI,j,k,iBlock)
             do i = nI+1, nI+2
                r = sqrt(x_Blk(i,j,1,iBlock)**2+y_Blk(i,j,1,iBlock)**2)
                State_VGB(Rho_,i,j,k,iBlock) = 1.0/r**(19.0/9.0)
                State_VGB(RhoUx_:p_-1,i,j,k,iBlock) = &
                     State_VGB(RhoUx_:p_-1,nI,j,k,iBlock)
                State_VGB(p_,i,j,k,iBlock) = State_VGB(Rho_,i,j,k,iBlock) &
                     *Temperature
             end do
          end do; end do
       case(4) ! r-direction in rz-geometry
          do k = -1, nK+2; do i = -1, nI+2
             Temperature = State_VGB(p_,i,nJ,k,iBlock) &
                  /State_VGB(Rho_,i,nJ,k,iBlock)
             do j = nJ+1, nJ+2
                r = sqrt(x_Blk(i,j,1,iBlock)**2+y_Blk(i,j,1,iBlock)**2)
                State_VGB(Rho_,i,j,k,iBlock) = 1.0/r**(19.0/9.0)
                State_VGB(RhoUx_:p_-1,i,j,k,iBlock) = &
                     State_VGB(RhoUx_:p_-1,i,nJ,k,iBlock)
                State_VGB(p_,i,j,k,iBlock) = State_VGB(Rho_,i,j,k,iBlock) &
                     *Temperature
             end do
          end do; end do
       end select
    case('usersemi')
       select case(iSide)
       case(2)
          StateSemi_VGB(1,nI+1,:,:,iBlock) = StateSemi_VGB(1,nI,:,:,iBlock)
       case(4)
          StateSemi_VGB(1,:,nJ+1,:,iBlock) = StateSemi_VGB(1,:,nJ,:,iBlock)
       end select
    end select

    IsFound = .true.

  end subroutine user_set_outerbcs

  !============================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModAdvance,    ONLY: State_VGB
    use ModConst,      ONLY: cKToEv
    use ModGeometry,   ONLY: x_Blk, y_Blk, Dx_Blk, y2
    use ModMain,       ONLY: nI, nJ, nK, Time_Simulation
    use ModNumConst,   ONLY: cPi
    use ModPhysics,    ONLY: No2Si_V, UnitTemperature_
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
    case('t0','temp0')
       select case(TypeProblem)
       case('gaussian')
          Spread = 4.0*HeatConductionCoef*Time_Simulation
          PlotVar_G = T0 + AmplitudeTemperature/(sqrt(cPi*Spread)) &
               *exp(-x_Blk(:,:,:,iBlock)**2/Spread)
       case('rz')
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

    real :: Rho, Pressure, Temperature
    real :: RhoSi, pSi, TemperatureSi

    character(len=*), parameter :: NameSub = 'user_material_properties'
    !--------------------------------------------------------------------------

    Rho = State_V(Rho_)
    RhoSi = Rho*No2Si_V(Rho_)

    if(present(EinternalSiIn))then
       pSi = EinternalSiIn*gm1
       Pressure = pSi*Si2No_V(UnitP_)
       Temperature = Pressure/Rho
       TemperatureSi = Temperature*No2Si_V(UnitTemperature_)
    elseif(present(TeSiIn))then
       TemperatureSi = TeSiIn
       Temperature = TemperatureSi*Si2No_V(UnitTemperature_)
       Pressure = Rho*Temperature
       pSi = Pressure*No2Si_V(UnitP_)
    else
       Pressure = State_V(p_)
       pSi = Pressure*No2Si_V(UnitP_)
       Temperature = Pressure/Rho
       TemperatureSi = Temperature*No2Si_V(UnitTemperature_)
    end if

    if(present(EinternalSiOut)) EinternalSiOut = pSi*inv_gm1
    if(present(TeSiOut)) TeSiOut = TemperatureSi
    if(present(PressureSiOut)) PressureSiOut = pSi

    if(present(CvSiOut)) CvSiOut = inv_gm1*Rho &
         *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_)

    if(present(HeatConductionCoefSiOut))then
       select case(TypeProblem)
       case('rmtv')
          HeatConductionCoefSiOut = Temperature**6.5/Rho**2 &
               *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_) &
               *No2Si_V(UnitU_)*No2Si_V(UnitX_)
       case default
          HeatConductionCoefSiOut = HeatConductionCoef &
               *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_) &
               *No2Si_V(UnitU_)*No2Si_V(UnitX_)
       end select
    end if

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
