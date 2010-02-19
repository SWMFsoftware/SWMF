!^CFG COPYRIGHT UM
!========================================================================
module ModUserTD99
  use ModConst
  use ModMain,       ONLY: x_,y_,z_,nByteReal
  use ModVarIndexes, ONLY: Ux_,Uy_,Uz_,Bx_,By_,Bz_
  use ModPhysics,    ONLY: No2Si_V,No2Io_V,UnitX_,UnitRho_,UnitU_,UnitB_,UnitJ_
  implicit none
  save
  !----------------------------------------------------------------------
  ! Additional normalization units::
  real:: No2SiUnitI, No2IoUnitI
  !----------------------------------------------------------------------
  ! Variables related to the position of the flux rope::
  real, parameter:: Li_TD99=cHalf
  real:: LongitudeTD99,LatitudeTD99,OrientationTD99
  !----------------------------------------------------------------------
  ! Variables related to the flux rope properties::
  real:: Itube_TD99,Rtube_TD99,atube_TD99,d_TD99,aratio_TD99
  real:: Mass_TD99,InvH0_TD99,Rho0_TD99
  real:: ItubeSaved
  !----------------------------------------------------------------------
  ! Variables related to the properties of the strapping field, Bq::
  integer:: nStepSaved=-1
  real, parameter:: AlphaRamp=9.52381E-04         !in [-]
  real, parameter:: VTransX=1.500E+03             !in [m/s]
  real, parameter:: VTransY=-2.900E+04            !in [m/s]
  real, parameter:: UVorCMax0=2.5                 !in units of 100 km/s
  real, parameter:: BqZMax0=3.768210E+01          !in [Gauss]
  real:: BqZMax,BqZMaxSaved,UVorCMax
  real:: q_TD99,L_TD99
  !----------------------------------------------------------------------
  ! Logical variables related to the magnetic field computation::
  logical:: UseTD99Perturbation=.false.
  logical:: DoTD99FluxRope=.false.
  logical:: DoEquilItube=.false.
  logical:: DoRevCurrent=.false.
  logical:: DoBqField=.false.
  logical:: UseVariedCurrent=.false.
  logical:: DoMaintainEpot=.false.
  real   :: CurrentRiseTime,CurrentStartTime
  !----------------------------------------------------------------------
  ! Declare the rotational matrix of coordinate transformation::
  real, dimension(3,3):: RotateTD99_DD

contains

  !=====================================================================!

  real function varied_current(Time)
    real, intent(in):: Time
    varied_current = min(&
         max(Time-CurrentStartTime,cZero)/CurrentRiseTime,cOne)*&
         ItubeSaved
  end function varied_current

  !=====================================================================!

  real function time_increment(Time)
    !--------------------------------------------------------------------
    real, intent(in):: Time
    time_increment = min(&
         max(Time-CurrentStartTime,cZero)/CurrentRiseTime,cOne)
  end function time_increment

  !=====================================================================!

  subroutine get_transformed_TD99fluxrope(RFace_D,BFRope_D,UVorT_D,&
       RhoFRope,Time)
    use ModCoordTransform, ONLY: rot_matrix_x,rot_matrix_y,&
         rot_matrix_z
    !--------------------------------------------------------------------
    real, dimension(3), intent(in):: RFace_D
    real, dimension(3), intent(out):: BFRope_D
    real, dimension(3), intent(out):: UVorT_D
    !--------------------------------------------------------------------
    real, intent(in), optional:: Time
    real, intent(out), optional:: RhoFRope
    !--------------------------------------------------------------------
    logical:: DoFirstCallTD99=.true.
    real:: atemp,Itemp
    real:: UVorR
    real, dimension(3):: B1FRopeTemp_D
    real, dimension(3):: R1Face_D,B1FRope_D,B1qField_D
    real, dimension(3):: UVorC_D,U1VorC_D
    !--------------------------------------------------------------------
    !\
    ! Initialize the TD99 model parameters once::
    !/
    !--------------------------------------------------------------------
    if (DoFirstCallTD99) then
       call init_TD99_parameters
       DoFirstCallTD99=.false.
    endif
    if (present(Time).and.UseVariedCurrent) &
         Itube_TD99 = varied_current(Time)
    !--------------------------------------------------------------------
    !\
    ! Check if the potential electric field needs to be applied.
    ! Maintain Epot for (CurrentStartTime < Time < CurrentRiseTime)::
    !/
    !--------------------------------------------------------------------
    if (present(Time).and.DoBqField) then
       DoMaintainEpot = (&
            (Time.gt.CurrentStartTime).and.&
            (Time.lt.CurrentRiseTime))
    else
       DoMaintainEpot = .false.
    endif
    !--------------------------------------------------------------------
    !\
    ! Compute the flux rope, and transform coordinates and vectors to
    ! position the flux rope in the desired way::
    !/
    !--------------------------------------------------------------------
    R1Face_D = matmul(RotateTD99_DD,RFace_D)
    if (DoTD99FluxRope) then
       call compute_TD99_FluxRope(R1Face_D,B1FRope_D,RhoFRope)
       if (DoRevCurrent) then
          Itemp = Itube_TD99; Itube_TD99 = -Itemp
          atemp = atube_TD99; atube_TD99 = aratio_TD99*atemp
          call compute_TD99_FluxRope(R1Face_D,B1FRopeTemp_D)
          B1FRope_D = B1FRope_D+B1FRopeTemp_D
          Itube_TD99 = Itemp
          atube_TD99 = atemp
       endif
    else
       B1FRope_D = cZero
       RhoFRope  = cZero
    endif
    U1VorC_D = cZero
    if (DoBqField) then
       if (present(Time)) then
          call compute_TD99_BqField(R1Face_D,B1qField_D,&
               U1VorC_D,DoMaintainEpot,Time)
       else
          call compute_TD99_BqField(R1Face_D,B1qField_D,&
               U1VorC_D,DoMaintainEpot)
       endif
       B1FRope_D = B1FRope_D+B1qField_D
    endif
    BFRope_D = matmul(B1FRope_D,RotateTD99_DD)
    UVorC_D  = matmul(U1VorC_D,RotateTD99_DD)
    !--------------------------------------------------------------------
    !\
    ! Compute the tangential component of the velocity field, UVorT_D,
    ! associated with the potential electric field, Epot::
    !/
    !--------------------------------------------------------------------
    UVorR          = dot_product(RFace_D,UVorC_D)
    UVorT_D(x_:z_) = UVorC_D-RFace_D(x_:z_)*UVorR
    !--------------------------------------------------------------------
  end subroutine get_transformed_TD99fluxrope

  !=====================================================================!

  subroutine init_TD99_parameters
    use ModProcMH,         ONLY: iProc
    use ModCoordTransform, ONLY: rot_matrix_x,rot_matrix_y,&
         rot_matrix_z
    use ModIO,             ONLY: iUnitOut, write_prefix
    !--------------------------------------------------------------------
    real:: AlphaRope,LInduct,WFRope,FootSepar
    !--------------------------------------------------------------------
    !\
    ! Define the SI normalization units here::
    !/
    !--------------------------------------------------------------------
    No2SiUnitI   = No2Si_V(UnitJ_)*No2Si_V(UnitX_)**2 ! in [A]
    !--------------------------------------------------------------------
    !\
    ! Define the I/O normalization units here::
    !/
    !--------------------------------------------------------------------
    No2IoUnitI = 1.0E+6*No2SiUnitI      ! in [microA]
    !--------------------------------------------------------------------
    !\
    ! Compute the magnetic energy, WFRope, associated with the portion
    ! of the flux rope current that is above the solar surface::
    !/
    !--------------------------------------------------------------------
    InvH0_TD99 = cGravitation*Msun/Rsun/No2Si_V(UnitU_)**2       ! in [-]
    AlphaRope  = cTwo*acos(d_TD99/Rtube_TD99)             ! in [rad]
    FootSepar  = Rtube_TD99*sin(AlphaRope/cTwo)/cE6       ! in [Mm]
    LInduct    = cMu*(AlphaRope/cTwo/cPi)*Rtube_TD99*log(cTwo**3*&
         (Rtube_TD99-d_TD99)/atube_TD99-cTwo+cQuarter)    ! in [H]
    WFRope     = LInduct*Itube_TD99**2/cTwo*cE6*cE1       ! in [ergs]
    !--------------------------------------------------------------------
    ! Compute the average density inside the flux rope assuming that the
    ! total amount of prominence mass is Mass_TD99 (=10^16g=10^13kg)::
    !--------------------------------------------------------------------
    Rho0_TD99  = Mass_TD99/(AlphaRope*Rtube_TD99*cPi*atube_TD99**2)
                                                          ! in [kg/m^3]
    !--------------------------------------------------------------------
    !\
    ! Define the normalized model parameters here::
    !/
    !--------------------------------------------------------------------
    ! Flux rope::
    !--------------------------------------------------------------------
    Rtube_TD99 = Rtube_TD99/No2Si_V(UnitX_)
    atube_TD99 = atube_TD99/No2Si_V(UnitX_)
    Itube_TD99 = Itube_TD99/No2SiUnitI
    Rho0_TD99  = Rho0_TD99/No2Si_V(UnitRho_)
    !--------------------------------------------------------------------
    ! Save the maximum value of the current for possible use in
    ! varied_current case::
    !--------------------------------------------------------------------
    ItubeSaved = Itube_TD99 
    !--------------------------------------------------------------------
    ! Strapping field::
    !--------------------------------------------------------------------
    d_TD99     = d_TD99/No2Si_V(UnitX_)
    L_TD99     = L_TD99/No2Si_V(UnitX_)
    q_TD99     = q_TD99/(No2Si_V(UnitB_)*No2Si_V(UnitX_)**2)
    !--------------------------------------------------------------------
    !\
    ! Construct the rotational matrix, RotateTD99_DD, to position the
    ! flux rope in the desired way on the solar surface::
    !/
    !--------------------------------------------------------------------
    RotateTD99_DD = matmul(rot_matrix_y(-cPi/cTwo),&
         rot_matrix_x(-OrientationTD99*cDegToRad))
    RotateTD99_DD = matmul(RotateTD99_DD,          &
         rot_matrix_y(LatitudeTD99*cDegToRad))
    RotateTD99_DD = matmul(RotateTD99_DD,          &
         rot_matrix_z(-LongitudeTD99*cDegToRad))
    !--------------------------------------------------------------------
    !\
    ! Print some stuff on the screen::
    !/
    !--------------------------------------------------------------------
    if (iProc==0) then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) '>>>>>>>>>>>>>>>>>>>                   <<<<<<<<<<<<<<<<<<<<<'
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) '    Twisted Flux Rope Model by Titov & Demoulin, 1999.     '
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) '>>>>>>>>>>>>>>>>>>>                   <<<<<<<<<<<<<<<<<<<<<'
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) '>>>>     Normalization Units (in MKS) in the model.    <<<<'
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'X0 = ',No2Si_V(UnitX_),'[m] = ',No2Si_V(UnitX_)/cE6,'[Mm] = Rsun'
       call write_prefix; write(iUnitOut,*) 'B0 = ',No2Si_V(UnitB_),'[T] = ',No2Si_V(UnitB_)*cE2*cE2,'[Gauss]'
       call write_prefix; write(iUnitOut,*) 'I0 = ',No2Si_V(UnitJ_),'[A] = ',No2Si_V(UnitJ_)*cE6,'[microA]'
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) '>>>>       Normalized values of model parameters.      <<<<'
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'd_TD99      = ',d_TD99,'[X0] ',              &
            d_TD99*No2Si_V(UnitX_)/cE6,'[Mm]'
       call write_prefix; write(iUnitOut,*) 'Rtube_TD99  = ',Rtube_TD99,'[X0] ',          &
            Rtube_TD99*No2Si_V(UnitX_)/cE6,'[Mm]'
       call write_prefix; write(iUnitOut,*) 'atube_TD99  = ',atube_TD99,'[X0] ',          &
            atube_TD99*No2Si_V(UnitX_)/cE6,'[Mm]'
       call write_prefix; write(iUnitOut,*) 'atube/Rtube = ',atube_TD99/Rtube_TD99,'[-]'
       call write_prefix; write(iUnitOut,*) 'Itube_TD99  = ',Itube_TD99,'[I0] ',          &
            Itube_TD99*No2SiUnitI,'[A]'
       call write_prefix; write(iUnitOut,*) 'aratio_TD99 = ',aratio_TD99,'[-]'
       call write_prefix; write(iUnitOut,*) 'Mass_TD99   = ',Mass_TD99*cE3,'[g] ',        &
            'InvH0_TD99  = ',InvH0_TD99,'[1/X0]' 
       call write_prefix; write(iUnitOut,*) 'Rho0_TD99   = ',Rho0_TD99,'[Rho0] = ',       &
            Rho0_TD99*No2Io_V(UnitRho_),'[g/cm^3]'
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'q_TD99      = ',q_TD99,'[B0*X0^2] ',         &
            q_TD99*No2Si_V(UnitB_)*No2Si_V(UnitX_)**2,'[T m^2]'
       call write_prefix; write(iUnitOut,*) 'L_TD99      = ',L_TD99,'[X0] ',              &
            L_TD99*No2Si_V(UnitX_)/cE6,'[Mm]'
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'Free energy of flux rope is ',WFRope,'Ergs.'
       call write_prefix; write(iUnitOut,*) 'Separation of flux rope ends is ',FootSepar, &
            'Mm, or',cPi*FootSepar*cE6/(cTwo*Rsun)*cRadToDeg,'deg.'
       call write_prefix; write(iUnitOut,*) ''
       if (UseVariedCurrent) then
          call write_prefix; write(iUnitOut,*) '>>>>>       UseVariedCurrent is set to .true.!!!      <<<<<'
          call write_prefix; write(iUnitOut,*) 'CurrentStartTime = ',CurrentStartTime,'[s]'
          call write_prefix; write(iUnitOut,*) 'CurrentRiseTime  = ',CurrentRiseTime,'[s]'
          call write_prefix; write(iUnitOut,*) ''
       endif
       call write_prefix; write(iUnitOut,*) '>>>>>>>>>>>>>>>>>>>>>>>> Action!!! <<<<<<<<<<<<<<<<<<<<<<<<'
       call write_prefix; write(iUnitOut,*) ''
    endif
    if (DoEquilItube) then
       !----------------------------------------------------------------
       !\
       ! Compute the equilibrium toroidal current, Itube_TD99, based
       ! on the force balance in direction normal to the surface of
       ! the flux tube.
       !/
       !----------------------------------------------------------------
       Itube_TD99 = cTwo*cFour*cPi*q_TD99*L_TD99*Rtube_TD99* &
            (L_TD99**2+Rtube_TD99**2)**(-(cOne+cHalf))     / &
            (alog(cTwo*cFour*Rtube_TD99/atube_TD99)        - &
            (cOne+cHalf)+Li_TD99/cTwo)                           ! in [-]
       WFRope    = LInduct*(Itube_TD99*No2SiUnitI)**2/cTwo*cE6*cE1 ! in [ergs]
    endif
    if (DoEquilItube.and.iProc==0) then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) '>>>>>>>>>>>>>>>>>>>      Wait...      <<<<<<<<<<<<<<<<<<<<<'
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'The strapping field, Bq, is added and the EQUILIBRIUM value'
       call write_prefix; write(iUnitOut,*) 'of Itube_TD99 is computed!!!'
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'The value of Itube_TD99 is reset to :: ',Itube_TD99
       call write_prefix; write(iUnitOut,*) 'The free energy of the flux rope is :: ',WFRope,'Ergs.'
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) '>>>>>>>>>>>>>>>>>>>   Now Action!!!   <<<<<<<<<<<<<<<<<<<<<'
       call write_prefix; write(iUnitOut,*) ''
    endif
    !--------------------------------------------------------------------
  end subroutine init_TD99_parameters

  !=====================================================================!

  subroutine compute_TD99_BqField(RFace_D,BqField_D,UVorC_D,&
       DoMaintainEpot,TimeNow)
    use ModMain,    ONLY: n_step,iteration_number
    !--------------------------------------------------------------------
    logical, intent(in):: DoMaintainEpot
    real, intent(in), optional:: TimeNow
    real, intent(in), dimension(3):: RFace_D
    real, intent(out), dimension(3):: BqField_D
    real, intent(out), dimension(3):: UVorC_D
    !--------------------------------------------------------------------
    ! Variables related to coordinates::
    !--------------------------------------------------------------------
    real:: R2Plus,R2Mins
    real, dimension(3):: RPlus_D,RMins_D
    !--------------------------------------------------------------------
    ! Variables related to computations of potential electric field::
    !--------------------------------------------------------------------
    real:: BqZOverBqZ0,BqZFunction
    real, dimension(3):: EpotC_D,UTranC_D
    real, dimension(3):: GradBqZ_D,GradPsiC_D
    !--------------------------------------------------------------------
    !\
    ! Compute the locations, RMins_D and RPlus_D, of the two magnetic
    ! charges, -/+q::
    !/
    !--------------------------------------------------------------------
    if (present(TimeNow)) then
       RPlus_D(x_) = RFace_D(x_)-L_TD99 - &
            VTransX*(TimeNow-CurrentStartTime)/No2Si_V(UnitX_)
       RMins_D(x_) = RFace_D(x_)+L_TD99 + &
            VTransX*(TimeNow-CurrentStartTime)/No2Si_V(UnitX_)
       RPlus_D(y_) = RFace_D(y_)        - &
            VTransY*(TimeNow-CurrentStartTime)/No2Si_V(UnitX_)
       RMins_D(y_) = RFace_D(y_)        + &
            VTransY*(TimeNow-CurrentStartTime)/No2Si_V(UnitX_)
    else
       RPlus_D(x_) = RFace_D(x_)-L_TD99
       RMins_D(x_) = RFace_D(x_)+L_TD99
       RPlus_D(y_) = RFace_D(y_)
       RMins_D(y_) = RPlus_D(y_)
    endif
    RPlus_D(z_) = RFace_D(z_)+d_TD99-cOne
    RMins_D(z_) = RPlus_D(z_)
    R2Plus = sqrt(dot_product(RPlus_D,RPlus_D))
    R2Mins = sqrt(dot_product(RMins_D,RMins_D))
    !--------------------------------------------------------------------
    !\
    ! Compute the field of the strapping magnetic field, BqField_D::
    !/
    !--------------------------------------------------------------------
    BqField_D(x_:z_) = q_TD99*     &
         (RPlus_D(x_:z_)/R2Plus**3-&
          RMins_D(x_:z_)/R2Mins**3)
    !--------------------------------------------------------------------
    !\
    ! Update the values of BqZMax and BqZMaxSaved once at each iteration::
    !/
    !--------------------------------------------------------------------
    if (n_step/=nStepSaved) then
       nStepSaved  = n_step
       BqZMaxSaved = BqZMax
       BqZMax  = cZero
    end if
    if (iteration_number==1) &
         BqZMaxSaved = BqZMax0/No2Io_V(UnitB_)
    BqZMax = max(abs(BqField_D(z_)),BqZMax)
    !--------------------------------------------------------------------
    !--------------------------------------------------------------------
    ! Apply Epot only if DoMaintainEpot=.true.!!!
    !--------------------------------------------------------------------
    !--------------------------------------------------------------------
    if (DoMaintainEpot) then
       !-----------------------------------------------------------------
       !\
       ! Compute the gradient of the z-component of the Bq field on the
       ! solar surface in Cartesian geometry -- GradBqZ_D(x_:z_)::
       !/
       !-----------------------------------------------------------------
       GradBqZ_D(x_:y_) = cOne       *&
            cThree*q_TD99*RMins_D(z_)*&
            (RMins_D(x_:y_)/R2Mins**5-&
             RPlus_D(x_:y_)/R2Plus**5)
       GradBqZ_D(z_)    = cZero
       !    GradBqZ_D(x_:z_) = GradBqZ_D(x_:z_)  - &
       !         RFace_D(x_:z_)*dot_product(RFace_D,GradBqZ_D)
       !-----------------------------------------------------------------
       !\
       ! Compute the gradient of the scalar potential in Cartesian
       ! geometry -- GradPsiC_D::
       !/
       !-----------------------------------------------------------------
       BqZOverBqZ0 = min(cOne,abs(BqField_D(z_))/BqZMaxSaved)
       BqZFunction = max(cZero,cOne-BqZOverBqZ0**2)
       if (BqZOverBqZ0.gt.cZero) then
          GradPsiC_D(x_:z_) = GradBqZ_D(x_:z_)*cTwo * &
               AlphaRamp*BqZFunction*BqZOverBqZ0    * &
               exp(-BqZFunction)
       else
          GradPsiC_D(x_:z_) = cZero
       endif
       !-----------------------------------------------------------------
       !\
       ! Compute the potential electric field on the solar surface to
       ! be applied at one of the spots -- EpotC_D(x_:z_).
       ! This is given by:
       ! EpotC_D(x_:z_) = sign(BqField_D(z_))*BqZOverBqZ0*GradPsiC_D::
       !/
       !-----------------------------------------------------------------
       EpotC_D(x_:z_) = sign(cOne,BqField_D(z_)) * &
            GradPsiC_D(x_:z_)*BqZOverBqZ0
       !-----------------------------------------------------------------
       !\
       ! Compute the plasma velocity, UVorC_D, associated with the static
       ! eletric field, EpotC_D.
       ! This is given by:
       ! UVorC_D(x_:z_) = sign(BqField_D(z_))*GradPsiC_D(x_:z_) X Ez,
       ! where Ez = (0,0,1)
       !/
       !-----------------------------------------------------------------
       UVorC_D(x_)    =  GradPsiC_D(y_)*sign(cOne,BqField_D(z_))
       UVorC_D(y_)    = -GradPsiC_D(x_)*sign(cOne,BqField_D(z_))
       UVorC_D(z_)    =  cZero
       UVorC_D(x_:z_) = UVorC_D(x_:z_)*UVorCMax0
       if (iteration_number.gt.2) &
            UVorC_D(x_:z_) = UVorC_D(x_:z_)*UVorCMax0/UVorCMax
       !-----------------------------------------------------------------
       !\
       ! Compute the translational velocity, UTranC_D, at which the two
       ! q-sources are moved with respect to each other.
       !/
       !-----------------------------------------------------------------
       if (present(TimeNow)) then
          UTranC_D(x_) = (VTransX/No2Si_V(UnitU_))*BqZOverBqZ0 * &
               sign(cOne,BqField_D(z_))*exp(-BqZFunction)
          UTranC_D(y_) = (VTransY/No2Si_V(UnitU_))*BqZOverBqZ0 * &
               sign(cOne,BqField_D(z_))*exp(-BqZFunction)
          UTranC_D(z_) = cZero
       else
          UTranC_D(x_:z_) = cZero
       endif
       !-----------------------------------------------------------------
       ! Add the translational velocity, UTranC_D, to UVorC_D::
       !-----------------------------------------------------------------
       UVorC_D = UVorC_D+UTranC_D
       !-----------------------------------------------------------------
    else
       !-----------------------------------------------------------------
       EpotC_D(x_:z_) = cZero
       UVorC_D(x_:z_) = cZero
       !-----------------------------------------------------------------
    endif
    !--------------------------------------------------------------------
  end subroutine compute_TD99_BqField

  !=====================================================================!

  subroutine compute_TD99_FluxRope(RFace_D,BFRope_D,RhoFRope)
    !-------------------------------------------------------------------!
    !\__                                                             __/!
    !    Twisted Magnetic Field Configuration by Titov & Demoulin '99   !
    !                                                                   !
    ! An instability that causes a CME eruption is expected to occur at !
    ! R > L*sqrt(2). For a detailed description of the initial state    !
    ! refer to A&A, 1999, v.351, pp.707-720                             !
    !                                                                   !
    ! ___  This module was written by Ilia Roussev on June 10, 2002 ___ !
    !/                                                                 \!
    !-------------------------------------------------------------------!
    real, intent(in), dimension(3):: RFace_D
    real, intent(out), dimension(3):: BFRope_D
    !--------------------------------------------------------------------
    real, intent(out), optional:: RhoFRope
    !--------------------------------------------------------------------
    real:: xxx,yyy,zzz,R2Face
    real:: RhoTB_TD99,Rperp_TD99
    real:: CHIin_TD99,CHIex_TD99
    real:: xUVx_TD99,xUVy_TD99,xUVz_TD99
    real:: ThetaUVy_TD99,ThetaUVz_TD99
    real:: RperpUVx_TD99,RperpUVy_TD99,RperpUVz_TD99
    real:: Kappa_TD99,dKappadx_TD99,dKappadr_TD99
    real:: KappaA_TD99,dKappaAdx_TD99,dKappaAdr_TD99
    !--------------------------------------------------------------------
    ! Complete elliptic integrals related variables::
    real:: K_elliptic, E_elliptic
    real:: K_ellipticA, E_ellipticA
    !--------------------------------------------------------------------
    ! Vector potential related variables::
    real:: Ak_TD99,dAkdk_TD99
    real:: AkA_TD99,dAkdkA_TD99,d2Akdk2A_TD99
    real:: AI_TD99,dAIdx_TD99,dAIdr_TD99
    real:: AIin_TD99,dAIindx_TD99,dAIindr_TD99
    real:: AIex_TD99,dAIexdx_TD99,dAIexdr_TD99
    ! Flux-rope related variables::
    real:: BIphix_TD99,BIphiy_TD99,BIphiz_TD99
    !--------------------------------------------------------------------
    !\
    ! Assign X,Y,Z coordinates at which to compute the magnetic field::
    !/
    !--------------------------------------------------------------------
    xxx = RFace_D(x_)
    yyy = RFace_D(y_)
    zzz = RFace_D(z_)
    R2Face = sqrt(dot_product(RFace_D,RFace_D))
    !--------------------------------------------------------------------
    !\
    ! Compute Rperp_TD99 and RhoTB_TD99::
    !/
    !--------------------------------------------------------------------
    Rperp_TD99 = sqrt(yyy**2+(zzz+d_TD99-cOne)**2)
    RhoTB_TD99 = sqrt(xxx**2+(Rperp_TD99-Rtube_TD99)**2)
    !--------------------------------------------------------------------
    !\
    ! Define the Heaviside step function in the internal region 
    ! (RhoTB_TD99<atube_TD99), CHIin_TD99, and the external one
    ! (RhoTB_TD99>atube_TD99), CHIex_TD99::
    !/
    !--------------------------------------------------------------------
    if (RhoTB_TD99.lt.atube_TD99) then
       CHIin_TD99 = cOne
       CHIex_TD99 = cZero
    else
       CHIin_TD99 = cZero
       CHIex_TD99 = cOne
    endif
    !--------------------------------------------------------------------
    !\
    ! Add the prominence material inside the flux rope, assuming that the
    ! total amount mass is 10^13kg, and that the desnity scale-height is
    ! the same as the pressure scale-height, 1/InvH0 (i.e., iso-thermal
    ! atmoshpere)::
    !/
    !--------------------------------------------------------------------
    if (present(RhoFRope)) &
         RhoFRope = Rho0_TD99*                 &
         exp(-cE1*(RhoTB_TD99/atube_TD99)**6)* &
         exp(-InvH0_TD99*abs(R2Face-cOne))    
    !--------------------------------------------------------------------
    !\
    ! Compute the field produced by the ring current, Itube_TD99, both
    ! inside and outside the torus, BI_TD99 = BFRope_D(x_:z_)::
    !/
    !--------------------------------------------------------------------
    ThetaUVy_TD99 = -(zzz+d_TD99-cOne)/Rperp_TD99
    ThetaUVz_TD99 = yyy/Rperp_TD99
    !--------------------------------------------------------------------
    !\
    ! Compute the toroidal field (BIphix_TD99, BIphiy_TD99, BIphiz_TD99)
    ! produced by the azimuthal current Iphi. This is needed to ensure
    ! that the flux rope configuration is force free. 
    !/
    !--------------------------------------------------------------------
    BIphix_TD99 = cZero
    BIphiy_TD99 = abs(Itube_TD99)/(cTwo*cPi*atube_TD99**2)*   &
         sqrt(CHIin_TD99*cTwo*(atube_TD99**2-RhoTB_TD99**2))* &
         ThetaUVy_TD99
    BIphiz_TD99 = abs(Itube_TD99)/(cTwo*cPi*atube_TD99**2)*   &
         sqrt(CHIin_TD99*cTwo*(atube_TD99**2-RhoTB_TD99**2))* &
         ThetaUVz_TD99
    !--------------------------------------------------------------------
    !\
    ! Compute the components of the unit vector in the plane of symmetry
    ! x=0::
    !/
    !--------------------------------------------------------------------
    RperpUVx_TD99 = cZero
    RperpUVy_TD99 = yyy/Rperp_TD99
    RperpUVz_TD99 = (zzz+d_TD99-cOne)/Rperp_TD99
    !--------------------------------------------------------------------
    !\
    ! Compute the components of the unit vector pointing in the positive 
    ! x-direction::
    !/
    !--------------------------------------------------------------------
    xUVx_TD99 = cOne
    xUVy_TD99 = cZero
    xUVz_TD99 = cZero
    !--------------------------------------------------------------------
    !\
    ! Define two model parameters, Kappa_TD99 and KappaA_TD99::
    !/
    !--------------------------------------------------------------------
    Kappa_TD99 = cTwo*sqrt(Rperp_TD99*Rtube_TD99 / &
         ((Rperp_TD99+Rtube_TD99)**2+xxx**2))
    KappaA_TD99 = cTwo*sqrt(Rperp_TD99*Rtube_TD99/ &
         (cFour*Rperp_TD99*Rtube_TD99+atube_TD99**2))
    !--------------------------------------------------------------------
    !\
    ! Truncate the value of Kappa_TD99::
    !/
    !--------------------------------------------------------------------
    if (abs(cOne-Kappa_TD99).lt.cTiny/cE1) &
         Kappa_TD99 = cOne-cTiny/cE1
    !--------------------------------------------------------------------
    !\
    ! Compute the vector potential in the internal, AIin_TD99, and
    ! external (outside the current torus), AIex_TD99, regions::   
    !/
    !--------------------------------------------------------------------
    call calc_elliptic_int_1kind(Kappa_TD99,K_elliptic)
    call calc_elliptic_int_2kind(Kappa_TD99,E_elliptic)
    Ak_TD99       = ((cTwo-Kappa_TD99**2)*K_elliptic      - &
         cTwo*E_elliptic)/Kappa_TD99
    dAkdk_TD99    = (cTwo-Kappa_TD99**2)*E_elliptic       / &
         (Kappa_TD99**2*(cOne-Kappa_TD99**2))             - &
         cTwo*K_elliptic/Kappa_TD99**2
    call calc_elliptic_int_1kind(KappaA_TD99,K_ellipticA)
    call calc_elliptic_int_2kind(KappaA_TD99,E_ellipticA)
    AkA_TD99      = ((cTwo-KappaA_TD99**2)*K_ellipticA    - &
         cTwo*E_ellipticA)/KappaA_TD99
    dAkdkA_TD99   = (cTwo-KappaA_TD99**2)*E_ellipticA     / &
         (KappaA_TD99**2*(cOne-KappaA_TD99**2))           - &
         cTwo*K_ellipticA/KappaA_TD99**2
    d2Akdk2A_TD99 = (((cFour+cThree)*KappaA_TD99**2-cFour - &
         KappaA_TD99**4)*E_ellipticA/(cOne-KappaA_TD99**2)+ &
         (cFour-(cOne+cFour)*KappaA_TD99**2)*K_ellipticA) / &
         (KappaA_TD99**3*(cOne-KappaA_TD99**2))
    !--------------------------------------------------------------------
    !\
    ! Define AIin_TD99 and AIex_TD99::
    !/
    !--------------------------------------------------------------------
    AIex_TD99     = Itube_TD99/(cTwo*cPi)*sqrt(Rtube_TD99 / &
         Rperp_TD99)*Ak_TD99
    AIin_TD99     = Itube_TD99/(cTwo*cPi)*sqrt(Rtube_TD99 / &
         Rperp_TD99)*(AkA_TD99+dAkdkA_TD99*(Kappa_TD99    - &
         KappaA_TD99))
    !--------------------------------------------------------------------
    !\
    ! Compute the vector potential, AI_TD99, of the magnetic field 
    ! produced by the ring current Itube_TD99 in the whole space::
    !/
    !--------------------------------------------------------------------
    AI_TD99 = CHIin_TD99*AIin_TD99+CHIex_TD99*AIex_TD99
    !--------------------------------------------------------------------
    !\
    ! Derive the BI_TD99 field from the corresponding vector potential,
    ! AI_TD99 (this involves the comp. of some nasty derivatives)::
    !/
    !--------------------------------------------------------------------
    dKappadx_TD99  = -xxx*Kappa_TD99/(xxx**2+(Rperp_TD99+Rtube_TD99)**2)
    dKappadr_TD99  = Kappa_TD99*(Rtube_TD99**2-Rperp_TD99**2+xxx**2)  / &
         (cTwo*Rperp_TD99*((Rtube_TD99+Rperp_TD99)**2+xxx**2))
    dKappaAdx_TD99 = cZero 
    dKappaAdr_TD99 = KappaA_TD99*atube_TD99**2/(cTwo*Rperp_TD99*(cFour* &
         Rperp_TD99*Rtube_TD99+atube_TD99**2))
    !--------------------------------------------------------------------
    !\
    ! Derivative of AIin_TD99 with respect to `x` and `rperp`:: 
    !/
    !--------------------------------------------------------------------
    dAIindx_TD99   = Itube_TD99/(cTwo*cPi)*sqrt(Rtube_TD99/Rperp_TD99)* &
         (dAkdkA_TD99*dKappadx_TD99)
    dAIindr_TD99   = Itube_TD99/(cTwo*cPi)*sqrt(Rtube_TD99/Rperp_TD99)* &
         (dAkdkA_TD99*dKappadr_TD99+d2Akdk2A_TD99*dKappaAdr_TD99      * &
         (Kappa_TD99-KappaA_TD99))-AIin_TD99/(cTwo*Rperp_TD99)
    !--------------------------------------------------------------------
    !\
    ! Derivative of AIex_TD99 with respect to `x` and `rperp`::
    !/
    !--------------------------------------------------------------------
    dAIexdx_TD99   = Itube_TD99/(cTwo*cPi)*sqrt(Rtube_TD99/Rperp_TD99)* &
         (dAkdk_TD99*dKappadx_TD99)
    dAIexdr_TD99   = Itube_TD99/(cTwo*cPi)*sqrt(Rtube_TD99/Rperp_TD99)* &
         (dAkdk_TD99*dKappadr_TD99)-AIex_TD99/(cTwo*Rperp_TD99)
    !--------------------------------------------------------------------
    !\
    ! Derivatives of AI with respect to `x` and `rperp`::
    !/
    !--------------------------------------------------------------------
    dAIdx_TD99 = CHIin_TD99*dAIindx_TD99+CHIex_TD99*dAIexdx_TD99
    dAIdr_TD99 = CHIin_TD99*dAIindr_TD99+CHIex_TD99*dAIexdr_TD99
    !--------------------------------------------------------------------
    !\              
    ! Obtain the BI_TD99 field in the whole space from the corresponding
    ! vector potential, AI_TD99 -->
    ! BI_TD99 = curl(AI_TD99*ThetaUV_TD99) = BFRope_D(x_:z_)::
    !/
    !--------------------------------------------------------------------
    BFRope_D(x_) = -dAIdx_TD99*RperpUVx_TD99+ &
         (dAIdr_TD99+AI_TD99/Rperp_TD99)*xUVx_TD99
    BFRope_D(y_) = -dAIdx_TD99*RperpUVy_TD99+ &
         (dAIdr_TD99+AI_TD99/Rperp_TD99)*xUVy_TD99
    BFRope_D(z_) = -dAIdx_TD99*RperpUVz_TD99+ &
         (dAIdr_TD99+AI_TD99/Rperp_TD99)*xUVz_TD99
    !--------------------------------------------------------------------
    ! Add the field of the azimuthal current, Iphi::
    !--------------------------------------------------------------------
    BFRope_D(x_) = BFRope_D(x_)+BIphix_TD99
    BFRope_D(y_) = BFRope_D(y_)+BIphiy_TD99
    BFRope_D(z_) = BFRope_D(z_)+BIphiz_TD99
    !--------------------------------------------------------------------
  contains
    !====================================================================
    subroutine calc_elliptic_int_1kind(Kappa,K_elliptic)
      !------------------------------------------------------------------
      real, intent(in):: Kappa
      real, intent(out):: K_elliptic
      !------------------------------------------------------------------
      integer:: iN
      real,parameter:: pK_LIMIT1 = cSqrtTwo/cTwo
      real,parameter:: pK_LIMIT2 = 0.9930000000000000000000000000000 
      real:: DESIRED_CEI_ACCURACY
      real:: pK,pK1,K_ell_sum,K_ell_sum_old
      real:: TwoN_1FactOverNFact
      !------------------------------------------------------------------
      !\
      ! Compute the complete elliptic integral of 1st kind from the series
      ! representations given by ...
      ! see formulae 8.113.1 (for 0<k<0.701) and 8.113.2 (for 0.701=<k<1)
      ! therein::
      !/
      ! The stability is ensured up to pK = 0.9935 (sin**-1=83.5deg)!!!
      !\
      ! Set the desired accuracy for the integral computation::KappaA_TD99
      !/
      !------------------------------------------------------------------
      if (nByteReal==8) then
         DESIRED_CEI_ACCURACY = cOne/cE15
      else
         DESIRED_CEI_ACCURACY = cE1**2/cE9
      endif
      !------------------------------------------------------------------
      !\
      ! Initialize some variables::
      !/
      !------------------------------------------------------------------
      iN                  = 1
      TwoN_1FactOverNFact = cOne 
      pK                  = Kappa
      pK1                 = sqrt(cOne-pK**2)
      !------------------------------------------------------------------
      !\
      ! Compute the CEI of 1st kind::
      !/
      !------------------------------------------------------------------
      if (abs(pK).lt.pK_LIMIT1) then
         K_ell_sum_old = cZero
         K_ell_sum     = (cOne+pK**2/cFour)*cPi/cTwo
         do while (abs(K_ell_sum-K_ell_sum_old).gt.DESIRED_CEI_ACCURACY)
            iN                  = iN+1
            TwoN_1FactOverNFact = TwoN_1FactOverNFact*(cTwo*iN-cOne)/iN
            K_ell_sum_old       = K_ell_sum
            K_ell_sum           = K_ell_sum+(TwoN_1FactOverNFact/cTwo**iN)**2* &
                 pK**(cTwo*iN)*cPi/cTwo
         enddo
      else
         if (abs(pK).lt.pK_LIMIT2) then
            K_ell_sum_old = cZero
            K_ell_sum     = (cOne+((cOne-pK1)/(cOne+pK1))**2/cFour)*&
                 cPi/(cOne+pK1)
            do while (abs(K_ell_sum-K_ell_sum_old).gt.DESIRED_CEI_ACCURACY)
               iN                  = iN+1
               TwoN_1FactOverNFact = TwoN_1FactOverNFact*(cTwo*iN-cOne)/iN
               K_ell_sum_old       = K_ell_sum
               K_ell_sum           = K_ell_sum+(TwoN_1FactOverNFact/cTwo**iN)**2* &
                    ((cOne-pK1)/(cOne+pK1))**(cTwo*iN)*cPi/(cOne+pK1)
            enddo
         else
            K_ell_sum     = alog(cFour/pK1)+(alog(cFour/pK1)-cOne)*pK1**2/cFour        + &
                 (alog(cFour/pK1)-cOne-cOne/cThree/cTwo)*pK1**4*(cThree/cTwo/cFour)**2 + &
                 (alog(cFour/pK1)-cOne-cOne/cThree/cTwo-cOne/(cTwo+cThree)/cThree)     * &
                 pK1**6*(cThree*(cTwo+cThree)/cTwo/cFour/(cTwo+cFour))**2              + &
                 (alog(cFour/pK1)-cOne-cOne/cThree/cTwo-cOne/(cTwo+cThree)/cThree      - &
                 cOne/(cThree+cFour)/cFour)*pK1**8*(cThree*(cTwo+cThree)*(cThree+cFour)/ &
                 cTwo/cFour/(cTwo+cFour)/(cFour+cFour))**2
         endif
      endif
      K_elliptic = K_ell_sum
      !------------------------------------------------------------------
    end subroutine calc_elliptic_int_1kind

    !====================================================================

    subroutine calc_elliptic_int_2kind(Kappa,E_elliptic)
      !------------------------------------------------------------------
      real, intent(in):: Kappa
      real, intent(out):: E_elliptic
      !------------------------------------------------------------------
      integer:: iN
      real,parameter:: pK_LIMIT1 = cSqrtTwo/cTwo
      real,parameter:: pK_LIMIT2 = 0.9990000000000000000000000000000 
      real:: DESIRED_CEI_ACCURACY
      real:: pK,pK1,E_ell_sum,E_ell_sum_old
      real:: TwoN_1FactOverNFact,TwoN_3FactOverNFact
      !------------------------------------------------------------------
      !\
      ! Compute the complete elliptic integral of 2nd kind from the series
      ! representations given by ...
      ! see formulae 8.114.1 (for 0<k<0.701) and 8.114.2 (for 0.701=<k<1)
      ! therein::
      !/
      ! The stability is ensured up to pK = 0.9993 (sin**-1=88.0deg)!!!
      !\
      ! Set the desired accuracy for the integral computation::
      !/
      !------------------------------------------------------------------    
      if (nByteReal==8) then
         DESIRED_CEI_ACCURACY = cOne/cE15
      else
         DESIRED_CEI_ACCURACY = cE1**2/cE9
      endif
      !------------------------------------------------------------------
      !\
      ! Initialize some variables::
      !/
      !------------------------------------------------------------------
      iN                  = 1
      TwoN_1FactOverNFact = cOne
      TwoN_3FactOverNFact = cOne
      pK                  = Kappa
      pK1                 = sqrt(cOne-pK**2)
      !------------------------------------------------------------------
      !\
      ! Compute the CEI of 2nd kind::
      !/
      !------------------------------------------------------------------
      if (abs(pK).lt.pK_LIMIT1) then
         E_ell_sum_old = cZero
         E_ell_sum     = (cOne-pK**2/cFour)*cPi/cTwo
         do while (abs(E_ell_sum-E_ell_sum_old).gt.DESIRED_CEI_ACCURACY)
            iN                  = iN+1
            TwoN_1FactOverNFact = TwoN_1FactOverNFact*(cTwo*iN-cOne)/iN
            E_ell_sum_old       = E_ell_sum
            E_ell_sum           = E_ell_sum-(TwoN_1FactOverNFact/cTwo**iN)**2/ &
                 (cTwo*iN-cOne)*pK**(cTwo*iN)*cPi/cTwo
         enddo
      else
         if (abs(pK).lt.pK_LIMIT2) then
            E_ell_sum_old = cZero
            E_ell_sum     = (cOne+((cOne-pK1)/(cOne+pK1))**2/cFour)* &
                 cPi*(cOne+pK1)/cFour
            do while (abs(E_ell_sum-E_ell_sum_old).gt.DESIRED_CEI_ACCURACY)
               iN                  = iN+1
               TwoN_3FactOverNFact = TwoN_3FactOverNFact*(cTwo*iN-cThree)/iN
               E_ell_sum_old       = E_ell_sum
               E_ell_sum           = E_ell_sum+(TwoN_3FactOverNFact/cTwo**iN)**2* &
                    ((cOne-pK1)/(cOne+pK1))**(cTwo*iN)*cPi*(cOne+pK1)/cFour
            enddo
         else
            E_ell_sum     = cOne+(alog(cFour/pK1)-cOne/cTwo)*pK1**2/cTwo                + &
                 (alog(cFour/pK1)-cOne-cOne/cThree/cFour)*pK1**4*(cThree/cFour/cFour)   + &
                 (alog(cFour/pK1)-cOne-cOne/cThree/cTwo-cOne/(cTwo+cThree)/(cTwo+cFour))* &
                 pK1**6*(cThree/cTwo/cFour)**2*(cTwo+cThree)/(cTwo+cFour)               + &
                 (alog(cFour/pK1)-cOne-cOne/cThree/cTwo-cOne/(cTwo+cThree)/cThree       - &
                 cOne/(cThree+cFour)/(cFour+cFour))*pK1**8*(cThree*(cTwo+cThree)        / &
                 cTwo/cFour/(cTwo+cFour))**2*(cThree+cFour)/(cFour+cFour)
         endif
      endif
      E_elliptic = E_ell_sum
      !------------------------------------------------------------------
    end subroutine calc_elliptic_int_2kind
    !--------------------------------------------------------------------  
  end subroutine compute_TD99_FluxRope
  !----------------------------------------------------------------------
end module ModUserTD99

!========================================================================

module ModUser
  use ModNumConst, ONLY: cHalf,cTwo,cThree,&
       cFour,cE1,cHundred,cHundredth,cZero,&
       cOne
  use ModMain,     ONLY: UseUserB0,UseUserHeating
  use ModSize,     ONLY: nI,nJ,nK,gcn,nBLK
  use ModUserTD99  ! To include TD99 flux rope.

  use ModUserEmpty, ONLY:               &
!!!       user_read_inputs,                &
!!!       user_init_session,               &
       user_set_ics,                    &
!!!       user_initial_perturbation,       &
!!!       user_set_boundary_cells,        &
!!!       user_face_bcs,                   &
       user_set_outerbcs,               &
!!!       user_specify_initial_refinement, &
!!!       user_amr_criteria,               &
       user_write_progress,             &
!!!       user_get_log_var,                &
       user_set_plot_var,               &
!!!       user_calc_sources,               &
!!!       user_get_b0,                     &
!!!       user_update_states,              &
       user_io_units

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'HELIOSPHERE, Manchester, Roussev'

  ! Cooling related parameters
  logical :: UseCooling   = .true.
  real    :: Temp0Cooling = 5.0e+7   ! [K]      target temperature
  real    :: TauCooling   = 100.0    ! [second] relaxation time

  !\
  ! PFSSM related variables::
  !/
  ! Maximum order of spherical harmonics
  integer, parameter:: N_PFSSM=90

  ! Weights of the spherical harmonics
  real, dimension(N_PFSSM+1,N_PFSSM+1):: g_nm,h_nm

  ! Temporary variable
  real, dimension(N_PFSSM+1):: FactRatio1

  ! Number of header lines in the file
  integer:: iHead_PFSSM

  ! Logical for Initialization 
  logical:: DoFirst=.true.

  ! Name of the input file
  character (LEN=32):: File_PFSSM

  ! Header string in the input file
  character (LEN=80):: Head_PFSSM

  ! Rotation angle around z-axis, in degrees,
  ! from the coordinate system of the component
  ! towards the coordinate system of the magnetic 
  ! map in the positive direction - hence, is positive. 
  ! Is set automatically to be equal to the 
  ! H(eliographic) L(ongitude) of C(entral) M(eridian)
  ! of the M(ap) minus 180 Deg, if in the input file 
  ! Phi_Shift is negative
  !
  real:: Phi_Shift

  ! Units of the magnetic field in the file including corrections
  ! relative to the magnetic units in the SC (Gauss)
  real :: UnitB


  ! Rs - radius of outer source surface where field is taken to be 0.
  ! Ro - radius of inner boundary for the potential
  ! H  - height of ??
  !  
  ! Typical values are Rs_PFSSM=31.0; Ro_PFSSM=1.0; H_PFSSM=0.1::
  real:: R_PFSSM,Rs_PFSSM,Ro_PFSSM,H_PFSSM

  !\
  ! Parameters related to the empirical heating::
  !/
  real:: InvH0,DegF_Ratio,Dens_Ratio
  real:: DegFrm1,Tnot=cOne
  real:: MaxB0_1,MaxB0_2,Bnot
  !\
  ! Gibson & Low 1998 related variables::
  !/
  logical:: DoFirst_GL=.true.
  real:: Mrope_GL98

contains

  !========================================================================
  !  SUBROUTINE user_read_inputs
  !========================================================================

  !\
  ! This subroutine allows the user to add input commands to the PARAM.in
  ! file that are specific to an application.  Although the user does not
  ! have to read input in the same manner as BATSRUS (with a #COMMAND
  ! followed by the data), the following example is set up to do this.
  ! This method is encouraged because it allows flexibility and it also
  ! will echo PARAM.in to standard out the same as the rest of the input
  ! file.
  !
  ! In the PARAM.in file user commands need to be enclosed in the
  ! #USERINPUTBEGIN and #USERINPUTEND commands as follows
  !
  ! #USERINPUTBEGIN
  !
  ! #USERSPECIFIEDCOMMANDS
  ! ...
  !
  ! #USERINPUTEND
  !
  ! The leading #USERINPUTBEGIN tells BATSRUS to call user_read_inputs.
  ! The user_read_inputs routine should return to the standard read_inputs
  ! on reading the command #USERINPUTEND.
  !
  ! The reading routine read_var can be used by the user to read all data
  ! types. BATSRUS uses a fortran 90 interface to overload the function so
  ! that it can read real, integer, logical or character strings. The syntax
  ! is read_var('characterstring',variable).  Where characterstring is
  ! typically just the variable name or a description of the variable's
  ! meaning.  For example: read_var('A random character',RandomCharacter)
  ! The advantage of using read_var is that it advances the line (iline)
  ! automatically and it also echos the command correctly back to standard
  ! out the way all other commands are echoed in read_inputs.
  !
  ! The user should declare their specific variables in ModUser above.
  !
  ! As with other user subroutines DO NOT MODIFY ANY GLOBAL VARIABLE DEFINED
  ! IN THE MODULES INCLUDED IN THIS SUBROUTINE UNLESS SPECIFIED!!
  subroutine user_read_inputs
    use ModMain
    use ModProcMH,    ONLY: iProc
    use ModReadParam
    use ModIO,        ONLY: write_prefix, write_myname, iUnitOut

    integer:: i
    character (len=100) :: NameCommand
    !-------------------------------------------------------------------------

    if(iProc==0.and.lVerbose > 0)then
       call write_prefix; write(iUnitOut,*)'User read_input HELIOSPHERE starts'
    endif
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#USER_FLAGS")
          call read_var('UseUserInnerBCs'         ,UseUserInnerBCs)
          call read_var('UseUserSource'           ,UseUserSource)
          call read_var('UseUserPerturbation'     ,UseUserPerturbation)
          call read_var('UseUserOuterBcs'         ,UseUserOuterBcs)
          call read_var('UseUserICs'              ,UseUserICs)
          call read_var('UseUserSpecifyRefinement',UseUserSpecifyRefinement)
          call read_var('UseUserLogFiles'         ,UseUserLogFiles)
          call read_var('UseUserWritePlot'        ,UseUserWritePlot)
          call read_var('UseUserAMR'              ,UseUserAMR)
          call read_var('UseUserEchoInput'        ,UseUserEchoInput)
          call read_var('UseUserB0'               ,UseUserB0)
          call read_var('UseUserInitSession'     ,UseUserInitSession)
          call read_var('UseUserUpdateStates'     ,UseUserUpdateStates)
       case("#COOLING")
          call read_var('UseCooling', UseCooling)
          if(UseCooling)then
             call read_var('Temp0Cooling',Temp0Cooling)
             call read_var('TauCooling',  TauCooling)
          end if
       case("#USEUSERHEATING")
          call read_var('UseUserHeating'          ,UseUserHeating)
       case("#PFSSM")
          call read_var('UseUserB0'  ,UseUserB0)
          if (UseUserB0) then
             call read_var('Ro_PFSSM'   ,Ro_PFSSM)
             call read_var('Rs_PFSSM'   ,Rs_PFSSM)
             call read_var('H_PFSSM'    ,H_PFSSM)
             call read_var('File_PFSSM' ,File_PFSSM)
             call read_var('iHead_PFSSM',iHead_PFSSM)
             call read_var('Phi_Shift'  ,Phi_Shift)
             call read_var('UnitB'      ,UnitB)
             call read_var('dt_UpdateB0',dt_UpdateB0)
             DoUpdateB0 = dt_updateb0 > 0.0
          endif
       case("#AWHEAT")
          call read_var('Bnot        ',Bnot)
          call read_var('Tnot        ',Tnot)
          call read_var('DegFrm1     ',DegFrm1)
          call read_var('DegF_Ratio  ',DegF_Ratio)
          call read_var('Dens_Ratio  ',Dens_Ratio)
       case("#TD99FLUXROPE")
          call read_var('UseTD99Perturbation' ,UseTD99Perturbation)
          call read_var('UseVariedCurrent'    ,UseVariedCurrent)
          call read_var('CurrentStartTime'    ,CurrentStartTime)
          call read_var('CurrentRiseTime '    ,CurrentRiseTime)
          call read_var('DoTD99FluxRope'      ,DoTD99FluxRope)
          call read_var('DoEquilItube'        ,DoEquilItube)
          call read_var('DoRevCurrent'        ,DoRevCurrent)
          call read_var('aratio_TD99'         ,aratio_TD99)
          call read_var('Itube_TD99'          ,Itube_TD99)
          call read_var('Rtube_TD99'          ,Rtube_TD99)
          call read_var('atube_TD99'          ,atube_TD99)
          call read_var('d_TD99'              ,d_TD99)
          call read_var('Mass_TD99'           ,Mass_TD99)
          call read_var('LongitudeTD99'       ,LongitudeTD99)
          call read_var('LatitudeTD99'        ,LatitudeTD99)
          call read_var('OrientationTD99'     ,OrientationTD99)
          call read_var('DoBqField'           ,DoBqField)
          call read_var('q_TD99'              ,q_TD99)
          call read_var('L_TD99'              ,L_TD99)
       case('#USERINPUTEND')
          if(iProc==0.and.lVerbose > 0)then
             call write_prefix;
             write(iUnitOut,*)'User read_input HELIOSPHERE ends'
          endif
          EXIT
       case default
          if(iProc==0) then
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

  !=====================================================================
  subroutine user_init_session

    ! Modify gravity to get better temperature !!! Do this only once

    use ModPhysics, ONLY: gBody
    character (len=*), parameter :: Name='user_init_session'
    !-------------------------------------------------------------------

    gBody  = gBody/Tnot

  end subroutine user_init_session

  !========================================================================
  subroutine user_set_boundary_cells(iBLK)

    ! Set the boundary cell information IsBoundaryCell_GI(:,:,:,ExtraBc_) 
    ! for a sphere of radius 1 around the origin.
    ! Allow resolution change.

    use ModGeometry
    use ModNumConst
    use ModParallel, ONLY: neiLEV,neiLtop, neiLbot,&
         neiLeast,neiLwest,neiLnorth,neiLsouth
    use ModMPCells,ONLY:DoOneCoarserLayer
    use ModBoundaryCells,ONLY:SaveBoundaryCells
    implicit none
    integer, intent(in):: iBLK

    integer:: i2,j2,k2,i,j,k
    real, parameter:: R=1.0, R2=R*R

    real :: DxHalf, DyHalf, DzHalf, DxQuarter, DyQuarter, DzQuarter
    real :: DxDown, DyDown, DzDown, DxUp, DyUp, DzUp

    !--------------------------------------------------------------------------
    IsBoundaryCell_GI(:,:,:,ExtraBc_) = R_BLK(:,:,:,iBLK)<R
    if(SaveBoundaryCells)return
    !\
    ! Neighbor solution block refinement levels::
    ! ( 0=neighbors at same level, 
    !  -1=neighbors at lower level,
    !  +1=neighbors at higher level,
    !  NOBLK=no neighbors).
    !/
    if(all(neiLEV(:,iBLK)==0))return

    ! Distance between fine and coarse cell centers
    DxHalf = 0.5*dx_BLK(iBLK)
    DyHalf = 0.5*dy_BLK(iBLK)
    DzHalf = 0.5*dz_BLK(iBLK)

    DxQuarter = 0.5*DxHalf
    DyQuarter = 0.5*DyHalf
    DzQuarter = 0.5*DzHalf

    ! The jump from the +++fine cell to the first coarse cell
    DxUp   = -DxHalf
    DyUp   = -DyHalf
    DzUp   = -DzHalf
    DxDown = -DxHalf
    DyDown = -DyHalf
    DzDown = -DzHalf

    if(.not.DoOneCoarserLayer)then
       ! The jump from the +++fine cell to the second coarse cell
       ! The coarce cell centers are twice as far as the local cell size
       DxUp   = DxUp   + 4*DxHalf
       DyUp   = DyUp   + 4*DyHalf
       DzUp   = DzUp   + 4*DzHalf
       DxDown = DxDown - 4*DxHalf
       DyDown = DyDown - 4*DyHalf
       DzDown = DzDown - 4*DzHalf
    endif

    ! Make any ghostcell covered by finer cell become a body cell, if
    ! any of the finer cell centers is inside the body
    ! Note: this includes all 8 fine cells !
    !       In the original implementation only the 4 closeby cells 
    !       were checked, which is good for typical face restricted
    !       message passing, but it is wrong if all 8 cells are restricted.

    if(neiLeast(iBLK)  == -1) call set_coarse_ghost(   0,   0,  1, nJ, 1, nK)
    if(neiLwest(iBLK)  == -1) call set_coarse_ghost(nI+1, nI+1, 1, nJ, 1, nK)
    if(neiLsouth(iBLK) == -1) call set_coarse_ghost(1, nI,    0,    0, 1, nK)
    if(neiLnorth(iBLK) == -1) call set_coarse_ghost(1, nI, nJ+1, nJ+1, 1, nK)
    if(neiLbot(iBLK)   == -1) call set_coarse_ghost(1, nI, 1, nJ,    0,    0)
    if(neiLtop(iBLK)   == -1) call set_coarse_ghost(1, nI, 1, nJ, nK+1, nK+1)

    ! Make all the ghostcells covered by a coarser cell become body cells,
    ! if the coarser cell center is inside the body

    if(neiLeast(iBLK)==+1)then
       i2=0
       do k2=2,nK,2; do j2=2,nJ,2
          IsBoundaryCell_GI(i2,j2-1:j2,k2-1:k2,ExtraBc_)=   &
               ( (x_BLK(i2,j2,k2,iBLK) - DxHalf)**2         &
               + (y_BLK(i2,j2,k2,iBLK) - DyHalf)**2         &
               + (z_BLK(i2,j2,k2,iBLK) - DzHalf)**2 ) < R2
          IsBoundaryCell_GI(i2-1,j2-1:j2,k2-1:k2,ExtraBc_)= &
               ( (x_BLK(i2,j2,k2,iBLK) + DxDown)**2         &
               + (y_BLK(i2,j2,k2,iBLK) - DyHalf)**2         &
               + (z_BLK(i2,j2,k2,iBLK) - DzHalf)**2 ) < R2
       end do; end do
    end if

    if(neiLwest(iBLK)==+1)then
       i2=nI+2
       do k2=2,nK,2; do j2=2,nJ,2
          IsBoundaryCell_GI(i2-1,j2-1:j2,k2-1:k2,ExtraBc_)= &
               ( (x_BLK(i2,j2,k2,iBLK) - DxHalf)**2         &
               + (y_BLK(i2,j2,k2,iBLK) - DyHalf)**2         &
               + (z_BLK(i2,j2,k2,iBLK) - DzHalf)**2 ) < R2
          IsBoundaryCell_GI(i2,j2-1:j2,k2-1:k2,ExtraBc_)=   &
               ( (x_BLK(i2,j2,k2,iBLK) + DxUp  )**2         &
               + (y_BLK(i2,j2,k2,iBLK) - DyHalf)**2         &
               + (z_BLK(i2,j2,k2,iBLK) - DzHalf)**2 ) < R2
       end do; end do
    end if

    if(neiLsouth(iBLK)==+1)then
       j2=0
       do k2=2,nK,2; do i2=2,nI,2
          IsBoundaryCell_GI(i2-1:i2,j2,k2-1:k2,ExtraBc_)=   &
               ( (x_BLK(i2,j2,k2,iBLK) - DxHalf)**2         &
               + (y_BLK(i2,j2,k2,iBLK) - DyHalf)**2         &
               + (z_BLK(i2,j2,k2,iBLK) - DzHalf)**2) < R2
          IsBoundaryCell_GI(i2-1:i2,j2-1,k2-1:k2,ExtraBc_)= &
               ( (x_BLK(i2,j2,k2,iBLK) - DxHalf)**2         &
               + (y_BLK(i2,j2,k2,iBLK) + DyDown)**2         &
               + (z_BLK(i2,j2,k2,iBLK) - DzHalf)**2) < R2
       end do; end do
    end if

    if(neiLnorth(iBLK)==+1)then
       j2=nJ+2
       do k2=2,nK,2; do i2=2,nI,2
          IsBoundaryCell_GI(i2-1:i2,j2-1,k2-1:k2,ExtraBc_)=         &
               ( (x_BLK(i2,j2,k2,iBLK) - DxHalf)**2       &
               + (y_BLK(i2,j2,k2,iBLK) - DyHalf)**2       &
               + (z_BLK(i2,j2,k2,iBLK) - DzHalf)**2) < R2
          IsBoundaryCell_GI(i2-1:i2,j2,k2-1:k2,ExtraBc_)=           &
               ( (x_BLK(i2,j2,k2,iBLK) - DxHalf)**2       &
               + (y_BLK(i2,j2,k2,iBLK) + DyUp  )**2       &
               + (z_BLK(i2,j2,k2,iBLK) - DzHalf)**2 ) < R2
       end do; end do
    end if

    if(neiLbot(iBLK)==+1)then
       k2=0
       do j2=2,nJ,2; do i2=2,nI,2
          IsBoundaryCell_GI(i2-1:i2,j2-1:j2,k2,ExtraBc_)=           &
               ( (x_BLK(i2,j2,k2,iBLK) - DxHalf)**2       &
               + (y_BLK(i2,j2,k2,iBLK) - DyHalf)**2       &
               + (z_BLK(i2,j2,k2,iBLK) - DzHalf)**2 ) < R2
          IsBoundaryCell_GI(i2-1:i2,j2-1:j2,k2-1,ExtraBc_)=         &
               ( (x_BLK(i2,j2,k2,iBLK) - DxHalf)**2       &
               + (y_BLK(i2,j2,k2,iBLK) - DyHalf)**2       &
               + (z_BLK(i2,j2,k2,iBLK) + DzDown)**2 ) < R2
       end do; end do
    end if

    if(neiLtop(iBLK)==+1)then
       k2=nK+2
       do j2=2,nJ,2; do i2=2,nI,2
          IsBoundaryCell_GI(i2-1:i2,j2-1:j2,k2-1,ExtraBc_)=         &
               ( (x_BLK(i2,j2,k2,iBLK) - DxHalf)**2      &
               + (y_BLK(i2,j2,k2,iBLK) - DyHalf)**2      &
               + (z_BLK(i2,j2,k2,iBLK) - DzHalf)**2) < R2
          IsBoundaryCell_GI(i2-1:i2,j2-1:j2,k2,ExtraBc_)=           &
               ( (x_BLK(i2,j2,k2,iBLK) - DxHalf)**2       &
               + (y_BLK(i2,j2,k2,iBLK) - DyHalf)**2       &
               + (z_BLK(i2,j2,k2,iBLK) + DzUp  )**2 ) < R2
       end do; end do
    end if

  contains

    subroutine set_coarse_ghost(iMin, iMax, jMin, jMax, kMin, kMax)

      ! Set boundary cell logical for the fine ghost cells 
      ! in the region determined by the arguments

      integer, intent(in) :: iMin, iMax, jMin, jMax, kMin, kMax
      !-------------------------------------------------------------------

      ! Get the distance of the fine cell closest to the origin
      do k=kMin,kMax; do j=jMin,jMax; do i=iMin, iMax
         IsBoundaryCell_GI( i, j, k, ExtraBc_) = &
              ( (abs(x_BLK(i, j, k, iBLK)) - DxQuarter)**2 &
              + (abs(y_BLK(i, j, k, iBLK)) - DyQuarter)**2 &
              + (abs(z_BLK(i, j, k, iBLK)) - DzQuarter)**2 ) < R2
      end do; enddo; enddo

    end subroutine set_coarse_ghost

  end subroutine user_set_boundary_cells

  !========================================================================
  subroutine user_face_bcs(VarsGhostFace_V)

    use ModSize,       ONLY: nDim,East_,West_,South_,North_,Bot_,    &
         Top_
    use ModMain,       ONLY: time_accurate,UseUserHeating,x_,y_,z_,  &
         UseRotatingFrame
    use ModVarIndexes, ONLY: rho_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, P_, &
         Ew_ => ScalarFirst_

    use ModGeometry,   ONLY: R_BLK
    use ModAdvance,    ONLY: nVar, State_VGB
    use ModPhysics,    ONLY: CosThetaTilt,SinThetaTilt,g,inv_g,      &
         inv_gm1,OmegaBody,No2Io_V(UnitB_)
    use ModNumConst,   ONLY: cZero,cHalf,cOne,cTwo,cTolerance,cTiny

    use ModBlockData, ONLY: use_block_data, put_block_data, get_block_data
    use ModFaceBc

    implicit none

    real, intent(out):: VarsGhostFace_V(nVar)
    !\
    ! User declared local variables go here::
    !/
    integer:: iCell,jCell,kCell
    real:: XFace,YFace,ZFace,RFace
    real:: VxFaceOutside,VyFaceOutside,VzFaceOutside
    real:: BxFaceOutside,ByFaceOutside,BzFaceOutside
    real:: VrFaceOutside,VthetaFaceOutside,VphiFaceOutside, &
         VrFaceInside,VthetaFaceInside,VphiFaceInside,    &
         BrFaceOutside,BthetaFaceOutside,BphiFaceOutside, &
         BrFaceInside,BthetaFaceInside,BphiFaceInside
    real:: CosTheta,SinTheta,CosPhi,SinPhi
    real, dimension(1:3):: location,v_phi
    real:: XFaceT,YFaceT,ZFaceT,sin2Theta_coronal_hole
    real:: CosThetaT,SinThetaT,CosPhiT,SinPhiT
    real:: DensCell,PresCell,GammaCell
    real:: B1dotR,BdotR,UdotR,SignB0n
    real, dimension(3):: RFace_D,B1_D,U_D,&
         B1t_D,B1n_D,Bt_D,Bn_D,Ut_D,Un_D
    !---------------------------------------------------------------------
    ! Variables related to the TD99 flux rope::
    !---------------------------------------------------------------------
    real:: BFRdotR,RhoFRope=cZero
    real, dimension(3):: BFRope_D,BFRn_D,BFRt_D,UVorT_D
    !--------------------------------------------------------------------------
    !\
    ! Calculation of boundary conditions should start here::
    !/
    !--------------------------------------------------------------------------
    !
    XFace = FaceCoords_D(1)
    YFace = FaceCoords_D(2)
    ZFace = FaceCoords_D(3)
    RFace = sqrt(XFace**2+YFace**2+ZFace**2)
    if (.not.UseUserHeating) then
       !\
       ! Apply some tricks to incorporate velocity
       ! shear at the boundary::
       !/
       RFace_D(x_)  = XFace/RFace
       RFace_D(y_)  = YFace/RFace
       RFace_D(z_)  = ZFace/RFace
       U_D (x_:z_)  = VarsTrueFace_V(Ux_:Uz_)
       UdotR        = dot_product(RFace_D,U_D)
       Un_D(x_:z_)  = UdotR*RFace_D(x_:z_)
       Ut_D(x_:z_)  = U_D(x_:z_)-Un_D(x_:z_)
       B1_D(x_:z_)  = VarsTrueFace_V(Bx_:Bz_)
       B1dotR       = dot_product(RFace_D,B1_D)
       B1n_D(x_:z_) = B1dotR*RFace_D(x_:z_)
       B1t_D        = B1_D-B1n_D
       BdotR        = dot_product(RFace_D,B0Face_D)
       Bn_D(x_:z_)  = BdotR*RFace_D(x_:z_)
       Bt_D(x_:z_)  = B0Face_D(x_:z_)-Bn_D(x_:z_)
       SignB0n      = sign(cOne,BdotR)
       !\
       ! Update BCs for velocity and induction field::
       !/
       VarsGhostFace_V(Ux_:Uz_) = -U_D(x_:z_)
       VarsGhostFace_V(Bx_:Bz_) = B1t_D(x_:z_)!-B1n_D(x_:z_)
       !\
       ! Compute the magnetic field of TD99 flux rope at RFace_D::
       !/
       if (DoTD99FluxRope.or.DoBqField) then
          if(UseVariedCurrent .or. .not.use_block_data(iBlock))then
             call get_transformed_TD99fluxrope(RFace_D,BFRope_D,&
                  UVorT_D,RhoFRope,time_now)
             if(.not. UseVariedCurrent)then
                ! Store calculated values
                call put_block_data(iBlock, RhoFRope)
                call put_block_data(iBlock, 3, BFRope_D)
                if(DoBqField)call put_block_data(iBlock, 3, UVorT_D)
             endif
          else
             ! Reuse calculated values
             call get_block_data(iBlock, RhoFRope)
             call get_block_data(iBlock, 3, BFRope_D)
             if(DoBqField)then
                call get_block_data(iBlock, 3, UVorT_D)
             else
                UVorT_D=0.0
             endif
          endif

          !\
          ! Compute the normal, BFRn_D, and tangential, BFRt_D,
          ! field components of the flux rope::
          !/
          BFRdotR       = dot_product(RFace_D,BFRope_D)
          BFRn_D(x_:z_) = BFRdotR*RFace_D(x_:z_)
          BFRt_D        = BFRope_D-BFRn_D
          !\
          ! Fix the normal component of the flux rope's field
          ! to BFRn_D at the Sun::
          !/
          VarsGhostFace_V(Bx_:Bz_) = VarsGhostFace_V(Bx_:Bz_)+&
               BFRn_D(x_:z_)!*cTwo
       endif
       !\
       ! Update BCs for the mass density, pressure and wave energy
       !/
       iCell = iFace; jCell = jFace; kCell = kFace
       select case(iSide)
       case(East_)
          iCell  = iFace
       case(West_)
          iCell  = iFace-1
       case(South_)
          jCell  = jFace
       case(North_)
          jCell  = jFace-1
       case(Bot_)
          kCell  = kFace
       case(Top_)
          kCell  = kFace-1
       case default
          call stop_mpi('user_face_bcs')
       end select
       call get_plasma_parameters_cell(iCell,jCell,kCell,iBlock,&
            DensCell,PresCell,GammaCell)
       VarsGhostFace_V(Rho_) = max(-VarsTrueFace_V(rho_     )+ &
            cTwo*(DensCell+RhoFRope),VarsTrueFace_V(rho_))
       VarsGhostFace_V(P_  ) = max(-VarsTrueFace_V(P_       )+ &
            cTwo*PresCell,VarsTrueFace_V(P_  ))
       VarsGhostFace_V(Ew_) = max(-VarsTrueFace_V(Ew_)+             & 
            cTwo*PresCell*(cOne/(GammaCell-cOne)-inv_gm1),          &
            VarsTrueFace_V(Ew_))
    else
       !\
       ! Rotate to spherical coordinates
       !/
       CosTheta = ZFace/RFace
       SinTheta = sqrt(XFace**2+YFace**2)/RFace
       CosPhi   = XFace/sqrt(XFace**2+YFace**2+cTolerance**2)
       SinPhi   = YFace/sqrt(XFace**2+YFace**2+cTolerance**2)
       VxFaceOutside = VarsTrueFace_V(Ux_)
       VyFaceOutside = VarsTrueFace_V(Uy_)
       VzFaceOutside = VarsTrueFace_V(Uz_)
       BxFaceOutside = VarsTrueFace_V(Bx_)
       ByFaceOutside = VarsTrueFace_V(By_)
       BzFaceOutside = VarsTrueFace_V(Bz_)
       VrFaceOutside = (VxFaceOutside*XFace      +&
            VyFaceOutside*YFace                  +&
            VzFaceOutside*ZFace)/RFace
       VthetaFaceOutside = ((VxFaceOutside*XFace +&
            VyFaceOutside*YFace)*ZFace           -&
            VzFaceOutside*(XFace**2+YFace**2))   /&
            (sqrt(XFace**2+YFace**2+cTolerance**2)*RFace)
       VphiFaceOutside = (VyFaceOutside*XFace    -&
            VxFaceOutside*YFace)*SinTheta        /&
            ((XFace**2+YFace**2+cTolerance**2)/RFace)
       BrFaceOutside = (BxFaceOutside*XFace      +&
            ByFaceOutside*YFace                  +&
            BzFaceOutside*ZFace)/RFace
       BthetaFaceOutside = ((BxFaceOutside*XFace +&
            ByFaceOutside*YFace)*ZFace           -&
            BzFaceOutside*(XFace**2+YFace**2))   /&
            (sqrt(XFace**2+YFace**2+cTolerance**2)*RFace)
       BphiFaceOutside = (ByFaceOutside*XFace    -&
            BxFaceOutside*YFace)*SinTheta        /&
            ((XFace**2+YFace**2+cTolerance**2)/RFace)
       RFace = cOne
       call coronal_hole_boundary(RFace,sin2Theta_coronal_hole)
       XFaceT    =  CosThetaTilt*XFace+SinThetaTilt*ZFace
       YFaceT    = YFace
       ZFaceT    = -SinThetaTilt*XFace+CosThetaTilt*ZFace
       CosThetaT = ZFaceT/RFace
       SinThetaT = sqrt(XFaceT**2+YFaceT**2)/RFace
       CosPhiT   = XFaceT/sqrt(XFaceT**2+YFaceT**2+cTolerance**2)
       SinPhiT   = YFaceT/sqrt(XFaceT**2+YFaceT**2+cTolerance**2)
       if (SinThetaT*SinThetaT > sin2Theta_coronal_hole) then
          !\
          ! At the base of closed field regions::
          !/
          VrFaceInside     = -VrFaceOutside
          VthetaFaceInside = -VthetaFaceOutside
          VphiFaceInside   = -VphiFaceOutside
          BrFaceInside     =  BrFaceOutside
          BthetaFaceInside =  BthetaFaceOutside
          BphiFaceInside   =  BphiFaceOutside
          VarsGhostFace_V(rho_     ) = VarsTrueFace_V(rho_     )
          VarsGhostFace_V(P_       ) = VarsTrueFace_V(P_       )
          VarsGhostFace_V(Ew_) = VarsTrueFace_V(Ew_)
       else
          !\
          ! At the base of open field regions::
          !/
          VrFaceInside     = cZero
          VthetaFaceInside = cZero
          VphiFaceInside   = cZero
          BrFaceInside     = cZero
          BthetaFaceInside = cZero 
          BphiFaceInside   = cZero
          VarsGhostFace_V(Rho_) = cOne
          VarsGhostFace_V(P_  ) = inv_g
          VarsGhostFace_V(Ew_ ) = VarsTrueFace_V(Ew_)
       endif
       !\
       ! Rotate back to cartesian coordinates::
       !/
       VarsGhostFace_V(Ux_) = VrFaceInside*XFace/RFace+&
            VthetaFaceInside*CosTheta*CosPhi          -&
            VphiFaceInside*SinPhi 
       VarsGhostFace_V(Uy_) = VrFaceInside*YFace/RFace+&
            VthetaFaceInside*CosTheta*SinPhi          +&
            VphiFaceInside*CosPhi
       VarsGhostFace_V(Uz_) = VrFaceInside*ZFace/RFace-&
            VthetaFaceInside*SinTheta
       VarsGhostFace_V(Bx_) = BrFaceInside*XFace/RFace+&
            BthetaFaceInside*CosTheta*CosPhi          -&
            BphiFaceInside*SinPhi
       VarsGhostFace_V(By_) = BrFaceInside*YFace/RFace+&
            BthetaFaceInside*CosTheta*SinPhi          +&
            BphiFaceInside*CosPhi
       VarsGhostFace_V(Bz_) = BrFaceInside*ZFace/RFace-&
            BthetaFaceInside*SinTheta
    endif
    !\
    ! Apply corotation:: Currently works only for the first body.
    !/
    if (.not.UseRotatingFrame) then
       !\
       ! The program is called which calculates the cartesian 
       ! corotation velocity::
       !/
       VarsGhostFace_V(Ux_) = VarsGhostFace_V(Ux_) -&
            cTwo*OmegaBody*YFace
       VarsGhostFace_V(Uy_) = VarsGhostFace_V(Uy_) +&
            cTwo*OmegaBody*XFace
    end if
  end subroutine user_face_bcs

  !========================================================================

  subroutine user_calc_sources

    !\
    ! This subroutine is used to calculate sources for the MHD equations.  
    ! The routine is called for each block separately, indexed by GlobalBlk
    ! Note that Source_VC(Energy_) and Source_VC(p_) must both be loaded. 
    !
    ! Calculates a location based heating if UseUserHeating is true. 
    ! Also fills in ModAdvance::tHeat0 for point implicit source term,
    ! but this is not used any longer.
    !/

    use ModMain,       ONLY: nI, nJ, nK, GlobalBlk, PROCtest, BLKtest, &
         UseUserHeating
    use ModAdvance,    ONLY: State_VGB, Source_VC, Rho_, p_, Energy_, tHeat0, &
         qHeat_BLK
    use ModGeometry,   ONLY: x_BLK, y_BLK, z_BLK, R_BLK, true_cell
    use ModPhysics,    ONLY: g, CosThetaTilt, SinThetaTilt, tHeat, &
         No2Si_V(UnitTemperature_), No2Si_V(UnitT_)
    use ModProcMH,     ONLY: iProc 
    !\
    ! Local variables
    !/
    integer:: i, j, k, iBlock
    real   :: CosTheta, SinTheta, CosPhi, SinPhi, Sin2Theta_coronal_hole
    real   :: xT, yT, zT
    real   :: Heating, Cooling, Temp, Temp0, Tau
    logical:: DoTest, DoTestMe
    !--------------------------------------------------------------------------
    if (.not.(UseUserHeating .or. UseCooling)) RETURN

    if (iProc==PROCtest.and.globalBLK==BLKtest) then
       call set_oktest('user_calc_sources', DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    iBlock = GlobalBlk

    if(UseUserHeating)then
       do k=1,nK; do j=1,nJ; do i=1,nI
          !\
          ! Colatitude dependent heating source term
          !/
          XT =  CosThetaTilt*x_BLK(i,j,k,iBlock)+&
               SinThetaTilt*z_BLK(i,j,k,iBlock)
          YT =  y_BLK(i,j,k,iBlock)
          ZT = -SinThetaTilt*x_BLK(i,j,k,iBlock)+&
               CosThetaTilt*z_BLK(i,j,k,iBlock)
          CosTheta = ZT/(R_BLK(i,j,k,iBlock)+cTolerance)
          SinTheta = sqrt(XT**2+YT**2)             /&
               (R_BLK(i,j,k,iBlock)+cTolerance)
          CosPhi = XT/sqrt(XT**2+YT**2+cTolerance**2)
          SinPhi = YT/sqrt(XT**2+YT**2+cTolerance**2)
          call coronal_hole_boundary(R_BLK(i,j,k,iBlock),&
               sin2Theta_coronal_hole)
          if (SinTheta**2 < sin2Theta_coronal_hole) then
             Theat0(i,j,k) = Theat
          else
             Theat0(i,j,k) = cOne
          end if
          Heating = State_VGB(rho_,i,j,k,iBlock)     *&
               qheat_BLK(i,j,k,iBlock)*(Theat0(i,j,k) &
               -g*State_VGB(P_,i,j,k,iBlock)/State_VGB(rho_,i,j,k,iBlock))

          Source_VC(P_     ,i,j,k) = Source_VC(P_,i,j,k) + (g-1)*Heating
          Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + Heating
       end do; end do; end do
    end if

    if(UseCooling)then
       ! Apply exponential cooling for tempretures T > T0 
       ! with time constant tau:
       ! dE/dt = -rho*max(0,T-T0)/tau 

       ! Convert to normalized units
       Temp0 = Temp0Cooling / No2Si_V(UnitTemperature_)
       Tau   = TauCooling   / No2Si_V(UnitT_)

       do k=1,nK; do j=1,nJ; do i=1,nI

          if(.not.true_cell(i,j,k,iBlock)) CYCLE

          Temp = State_VGB(P_,i,j,k,iBlock)/ State_VGB(Rho_,i,j,k,iBlock)

          if (Temp < Temp0) CYCLE

          Cooling = State_VGB(Rho_,i,j,k,iBlock)*(Temp - Temp0) / Tau

          Source_VC(P_     ,i,j,k) = Source_VC(P_,i,j,k)      - Cooling * (g-1)
          Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) - Cooling

       end do; end do; end do
    end if

  end subroutine user_calc_sources

  !========================================================================

  subroutine get_plasma_parameters_cell(iCell,jCell,kCell,iBlock,&
       DensCell,PresCell,GammaCell)
    !--------------------------------------------------------------------------
    !
    ! This module computes the cell values for density and pressure assuming
    ! a politropic equation of state with variable gamma = [2+n(T)]/n(T),
    ! where n(T)=n0+n1*T^2.
    ! This subroutine is written by ILR on May 29, 2003.
    ! Last updated is made by IVS and ILR on Nov 2, 2004.
    !--------------------------------------------------------------------------
    !
    use ModVarIndexes, ONLY: Bx_,By_,Bz_,P_
    use ModGeometry,   ONLY: x_BLK,y_BLK,z_BLK
    use ModNumConst,   ONLY: cZero,cHalf,cOne,cTwo,cThree,cFour,cTiny
    use ModAdvance,    ONLY: B0_DGB,  &
         State_VGB
    use ModMain,       ONLY:x_,y_,z_
    use ModPhysics,    ONLY: g,inv_g,No2Io_V(UnitB_)
    implicit none
    !--------------------------------------------------------------------------
    integer, intent(in):: iCell,jCell,kCell,iBlock
    real, intent(out):: DensCell,PresCell,GammaCell
    !--------------------------------------------------------------------------
    real, parameter:: n0=cFour
    real:: AAc,BBc,Fn1,Fg1
    real:: BrCell,B2Cell
    real:: BetaCell,BetaFactor
    real:: XCell,YCell,ZCell,RCell
    real:: Temp_Ratio,TempCell,DegFrmCell
    real:: DegF_Modulation,Dens_Modulation,Temp_Modulation
    !--------------------------------------------------------------------------
    ! Set MaxB0 stuff
    MaxB0_1 = Bnot/No2Io_V(UnitB_)
    MaxB0_2 = 2.00E+01/No2Io_V(UnitB_)
    !\
    ! Get cell coordinates and radial distance from the Sun::
    !/
    XCell = x_BLK(iCell,jCell,kCell,iBlock)
    YCell = y_BLK(iCell,jCell,kCell,iBlock)
    ZCell = z_BLK(iCell,jCell,kCell,iBlock)
    RCell = sqrt(XCell**2+YCell**2+ZCell**2)
    if (RCell>0.5) then
       !\
       ! Get the absolute value of the radial component of the magnetic filed::
       !/
       BrCell = abs(&
            (XCell*B0_DGB(x_,iCell,jCell,kCell,iBlock) +&
            YCell*B0_DGB(y_,iCell,jCell,kCell,iBlock) +&
            ZCell*B0_DGB(z_,iCell,jCell,kCell,iBlock))/&
            RCell)
       B2Cell = sqrt(sum(B0_DGB(:,iCell,jCell,kCell,iBlock)**2))
       !\
       ! Modulate the degrees of freedom so that in the CS (slow wind) the
       ! number of degrees of freedom is LARGER (~27) than in the open field
       ! regions (~13; fast wind).  This will result in a more isothermal
       ! solution in the helmet streamer belt than in the open field regions.
       ! Therefore, the kinetic gas pressure at the Sun in open field regions
       ! is somewhat GREATER than that in the closed field regions.  However,
       ! by imposing the density and temperature modulation described below,
       ! we increase the kinetic gas pressure in closed field regions.  All
       ! this is done to achieve a good agreement with observations at 1AU....
       !/
       DegF_Modulation = DegF_Ratio/&
            (cOne+min(cOne,BrCell*RCell**3/MaxB0_1)*&
            (DegF_Ratio-cOne))
       !\
       ! Modulate the temperature and density in order to have higher
       ! temperature (~2.5) and lower density (~1/5) in open field regions
       ! [abs(Br)<6G] compared to closed ones [abs(Br)~0].  This is done
       ! to reproduce the observed bi-modal structure of the solar wind.
       !/
       Temp_Ratio = Dens_Ratio/cTwo
       Temp_Modulation = cOne+min(cOne,BrCell*RCell**3/MaxB0_1)*&
            (Temp_Ratio-cOne)
       Dens_Modulation = cOne/&
            (cOne+min(cOne,BrCell*RCell**3/MaxB0_1)*&
            (Dens_Ratio-cOne))
       !\
       ! Increase artificially the mass density in regions where
       ! B2Cell>MaxB0_2 (~25G), so that to maintain a reasonable value
       ! of the Alfven speed for B2Cell>MaxB0_2.
       !/
       if (B2Cell>MaxB0_2) then
          Dens_Modulation = min(cTwo+cHalf,&
               (B2Cell/MaxB0_2)**2/Dens_Ratio)
          Temp_Modulation = cOne/Dens_Modulation
       end if
    else
       DegF_Modulation = cOne
       Temp_Modulation = cOne
       Dens_Modulation = cOne
    end if
    !\
    ! Obtain the temperature at the given cell::
    !/
    Fg1 = (cOne+cHalf*n0+cHalf*cHalf*DegFrm1*DegF_Modulation)/&
         (cOne+cHalf*n0+cHalf*cHalf*DegFrm1)
    Fn1 = DegFrm1*DegF_Modulation
    AAc = cFour*(cOne+cHalf*n0)/Fn1
    BBc = cFour*(cOne+cHalf*n0+cHalf*cHalf*Fn1)/Fn1/RCell
    TempCell = cTwo*BBc/(sqrt(AAc**2+cFour*BBc)+AAc) !=1 as long as RCell=1
    !\
    ! Compute the plasma beta, BetaCell, in order to diminish the
    ! turbulent heating in the current sheet, where BetaCell is
    ! large::
    !/
    BetaCell = cTwo*&
         State_VGB(P_ ,iCell,jCell,kCell,iBlock)/    &
         max(cTiny, sum((State_VGB(Bx_:Bz_,iCell,jCell,kCell,iBlock)+    &
                            B0_DGB(:       ,iCell,jCell,kCell,iBlock))**2))
         
    !\
    ! Compute a Beta-multiplier, which will be used to
    ! diminish the turbulent heating in the current sheet::
    !/
    BetaFactor = max(cZero,cOne-BetaCell*(RCell/2.50E+01))
    !\
    ! Use BetaFactor to set degrees of freedom close to n0(=4)
    ! in the heliospheric current sheet (BetaFactor approx 0),
    ! so that GammaCell=1.5 there::
    !
    ! Also, below we assume a quadratic dependence of the
    ! degrees of freedom on the plasma temperature.
    !/
    DegFrmCell = n0+(Fn1*TempCell**2)*BetaFactor
    GammaCell  = (DegFrmCell+cTwo)/DegFrmCell 
    !\
    ! Obtain the mass density and pressure in the given cell::
    !/
    PresCell  = inv_g*Dens_Modulation*Temp_Modulation/Fg1*&
         (TempCell**(cOne+cHalf*n0))*exp(-cHalf*Fn1*(cOne-TempCell))
    DensCell  = g*PresCell*Fg1/TempCell/Temp_Modulation
  end subroutine get_plasma_parameters_cell

  !========================================================================
  subroutine user_initial_perturbation

    use ModMain,      ONLY: nI,nJ,nK,nBLK,                           &
         unusedBLK,UseUserHeating,UseUserB0,gcn,x_,y_,z_
    use ModIO,        ONLY: restart
    use ModVarIndexes,ONLY: Rho_, RhoUx_, RhoUy_, RhoUz_, Bx_, By_, Bz_, P_, &
         Ew_ => ScalarFirst_ 
    use ModAdvance,   ONLY: State_VGB,B0_DGB,tmp1_BLK,tmp2_BLK
    use ModProcMH,    ONLY: iProc,nProc,iComm
    use ModNumConst,  ONLY: cZero,cQuarter,cHalf,cOne,cTwo,cE1,cE9,  &
         cTolerance,cThree
    use ModConst,     ONLY: Rsun,Msun,cGravitation
    use ModGeometry,  ONLY: x_BLK,y_BLK,z_BLK,R_BLK,cV_BLK,x2,y2,z2
    use ModPhysics,   ONLY: Gbody,g,inv_g,gm1,inv_gm1,ModulationP,   &
         ModulationRho,UseFluxRope,rot_period_dim,OmegaBody,Rbody,   &
         No2Si_V(UnitU_),No2Si_V(UnitRho_),No2Si_V(UnitX_),&
         No2Io_V(UnitEnergydens_),No2Io_V(UnitT_),No2Io_V(UnitB_)

    use ModEnergy, ONLY: calc_energy_cell
    !\
    ! Variables required by this user subroutine::
    !/
    integer:: i,j,k,iBLK,iError
    logical:: oktest,oktest_me
    real:: volume
    real:: xx,yy,zz,RR,ROne,Rmax
    real:: rho_GL98,p_GL98
    real:: Bx_GL98,By_GL98,Bz_GL98
    real:: Dens_BLK,Pres_BLK,Gamma_BLK
    real, dimension(3):: R_GL98_D,B_GL98_D
    real, dimension(3):: R_TD99_D,B_TD99_D,U_TD99  ! To include TD99 flux rope.
    real:: Rho_TD99=cZero                          ! To include TD99 flux rope.
    !
    !--------------------------------------------------------------------------
    !\
    ! Variable meanings:
    !/
    !--------------------------------------------------------------------------
    !
    call set_oktest('user_initial_perturbation',oktest,oktest_me)
    !\
    ! Initialize some auxilary variables::
    !/
    Rbody  = cOne
    Mrope_GL98 = cZero

    InvH0 = cGravitation*Msun/Rsun/No2Si_V(UnitU_)**2
    do iBLK=1,nBLK
       if (unusedBLK(iBLK)) CYCLE   
       do k=1,nK;do j=1,nJ; do i=1,nI
          if ((.not.UseUserHeating).and.(.not.restart)) then
             xx = x_BLK(i,j,k,iBLK)
             yy = y_BLK(i,j,k,iBLK)
             zz = z_BLK(i,j,k,iBLK)
             RR = sqrt(xx**2+yy**2+zz**2+cTolerance**2)
             ROne  = max(cOne,RR)
             Rmax  = max(2.1E+01,sqrt(x2**2+y2**2+z2**2))
             State_VGB(Bx_      ,i,j,k,iBLK) = cZero
             State_VGB(By_      ,i,j,k,iBLK) = cZero
             State_VGB(Bz_      ,i,j,k,iBLK) = cZero
             call get_plasma_parameters_cell(i,j,k,iBLK,&
                  Dens_BLK,Pres_BLK,Gamma_BLK)
             State_VGB(rho_     ,i,j,k,iBLK) = Dens_BLK
             State_VGB(P_       ,i,j,k,iBLK) = Pres_BLK
             State_VGB(rhoUx_   ,i,j,k,iBLK) = Dens_BLK*&
                  4.0E+01*((ROne-cOne)/(Rmax-cOne))*xx/RR
             State_VGB(rhoUy_   ,i,j,k,iBLK) = Dens_BLK*&
                  4.0E+01*((ROne-cOne)/(Rmax-cOne))*yy/RR
             State_VGB(rhoUz_   ,i,j,k,iBLK) = Dens_BLK*&
                  4.0E+01*((ROne-cOne)/(Rmax-cOne))*zz/RR
             State_VGB(Ew_,i,j,k,iBLK) = &
                  Pres_BLK*(cOne/(Gamma_BLK-cOne)-inv_gm1)
          endif
          !----------------------------------------------------------------
          !\
          ! Add Titov & Demoulin (TD99) flux rope here:: 
          !/
          !----------------------------------------------------------------
          if (UseTD99Perturbation.and.restart) then
             !-------------------------------------------------------------
             !\
             ! Assign the coordinates at which to compute the field::
             !/
             !-------------------------------------------------------------
             R_TD99_D(x_) = x_BLK(i,j,k,iBLK)
             R_TD99_D(y_) = y_BLK(i,j,k,iBLK)
             R_TD99_D(z_) = z_BLK(i,j,k,iBLK)
             !-------------------------------------------------------------
             !\
             ! Computed the magnetic field::
             !/
             !-------------------------------------------------------------
             if (.not.UseVariedCurrent) &
                  call get_transformed_TD99fluxrope(R_TD99_D,B_TD99_D,&
                  U_TD99,Rho_TD99)
             !-------------------------------------------------------------
             !\
             ! Add the flux rope field to the induction field, B1::
             !/
             !-------------------------------------------------------------
             State_VGB(rho_,i,j,k,iBLK)          = &
                  State_VGB(rho_,i,j,k,iBLK)+Rho_TD99
             State_VGB(Bx_ ,i,j,k,iBLK)          = &
                  State_VGB(Bx_ ,i,j,k,iBLK)+B_TD99_D(x_)
             State_VGB(By_ ,i,j,k,iBLK)          = &
                  State_VGB(By_ ,i,j,k,iBLK)+B_TD99_D(y_)
             State_VGB(Bz_ ,i,j,k,iBLK)          = &
                  State_VGB(Bz_ ,i,j,k,iBLK)+B_TD99_D(z_)
             !-------------------------------------------------------------
          endif
          !----------------------------------------------------------------

          !----------------------------------------------------------------
          if (UseFluxRope) then
             R_GL98_D(x_) = x_BLK(i,j,k,iBLK)
             R_GL98_D(y_) = y_BLK(i,j,k,iBLK)
             R_GL98_D(z_) = z_BLK(i,j,k,iBLK)
             call add_GL98_fluxrope(R_GL98_D,rho_GL98,p_GL98,&
                  B_GL98_D)
             State_VGB(Bx_,i,j,k,iBLK)          = &
                  State_VGB(Bx_,i,j,k,iBLK)+B_GL98_D(x_)
             State_VGB(By_,i,j,k,iBLK)          = &
                  State_VGB(By_,i,j,k,iBLK)+B_GL98_D(y_)
             State_VGB(Bz_,i,j,k,iBLK)          = &
                  State_VGB(Bz_,i,j,k,iBLK)+B_GL98_D(z_)
             !\
             ! Add just `ModulationRho' times of the CME mass
             ! to the mass density::
             !/
             if ((State_VGB(rho_,i,j,k,iBLK)+     &
                  ModulationRho*rho_GL98) < cOne*&
                  State_VGB(rho_,i,j,k,iBLK)) then
                State_VGB(rho_,i,j,k,iBLK)       = &
                     cOne*State_VGB(rho_,i,j,k,iBLK)
             else
                State_VGB(rho_,i,j,k,iBLK)       = &
                     State_VGB(rho_,i,j,k,iBLK)  + &
                     ModulationRho*rho_GL98
             endif
             !\
             ! Add just `ModulationP' times of the CME pressure
             ! to the kinetic pressure::
             !/
             if ((State_VGB(P_,i,j,k,iBLK)+&
                  ModulationP*p_GL98) < cOne* &
                  State_VGB(P_,i,j,k,iBLK)) then
                State_VGB(P_,i,j,k,iBLK)     = &
                     cOne*State_VGB(P_,i,j,k,iBLK)
             else 
                State_VGB(P_,i,j,k,iBLK)     = &
                     State_VGB(P_,i,j,k,iBLK)+ &
                     ModulationP*p_GL98
             endif
             !\
             ! Calculate the mass added to the flux rope::
             !/
             Mrope_GL98 = Mrope_GL98+ModulationRho*&
                  rho_GL98*cV_BLK(iBLK)
          endif
       end do; end do; end do
       !\
       ! Update the total energy::
       !/
       call calc_energy_cell(iBLK)
    end do
    !\
    ! Write out some statistics::
    !/
    call post_init_stat

  end subroutine user_initial_perturbation

  !========================================================================

  subroutine add_GL98_fluxrope(R_GL98_D,rho_GL98,p_GL98,B_GL98_D)
    !
    !---------------------------------------------------------------------------
    ! PARAMETER LIST: cme_a, cme_r1, cme_r0, cme_a1, cme_rho1, cme_rho2, B1_dim,
    !                  RHOsun,Vscl
    ! Definition of Parameters used for the initial state
    !   cme_a    = contraction distance as in   r --> r -a
    !   cme_r1   = distance of flux rope from sun center = 1.2
    !   cme_r0   = radius of flux rope
    !   cme_a1   = constant for setting pressure in flux rope
    !   Rscl     = 1.0  scaled radius of the sun
    !   RHOscl   = 1.0  scaled density of RHOsun
    !   SSPscl   = 1.0  scaled soundspeed of the sun
    !   rho1scl  = uniform backround density of the solution before contraction
    !   rho2scl  = background powerlaw density added to after contraction
    !   B1scl    = magnetic field strength parameter of the flux rope
    !   Gsun     = gravitational acceleration at sun's surface = 2.734e4 cm/s**2
    !   Vscl     = V/SSPsun     
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////
    !=====================================================================
    !\
    ! Calculates magnetic field, pressure and density for a coronal flux rope 
    ! capable of self-similar expansion and erupting in a CME.
    ! The analytical solution is taken from Gibson and Low 
    ! Astrophysical Journal, Vol 493, p. 460.
    !
    ! Written by Chip Manchester Jan 18 2001
    ! Rewritten by Chip Nov 29 2001 for flux rope injection
    !/
    !   Bug fixes
    !   March 18       dpres_1dr1 = cme_a1*(dA2dr*dr2dr1 + dA2dth*dth2dr1)
    !   change to..... dpres_1dr1 = a1scl*(dA2dr*dr2dr1 + dA2dth*dth2dr1)
    !   without above fix, same as used for runs 12, 13, 14
    !   Feb   07, 2002 Br1 is changed to Br1**2 in density calc thanks to Ilia
    !   Feb   12, 2002 expression for ga0r0 is fixed     
    !   Feb   22, 2002 derivative of B2*2/dr1 is fixed
    !   Feb   22, 2002 angles in 2nd coordinate system are redefined
    !
    !---------------------------------------------------------------------------
    !
    use ModMain,           ONLY: x_,y_,z_
    use ModNumConst,       ONLY: cZero,cOne,cTwo,cFour,cPi,&
         cDegToRad
    use ModProcMH,         ONLY: iProc
    use ModCoordTransform, ONLY: rot_matrix_x,rot_matrix_y,&
         rot_matrix_z
    use ModPhysics,        ONLY: Gbody,g,inv_g,gm1,inv_gm1,&
         cme_a,cme_r1,cme_r0,cme_a1,a1scl,rho1scl,rho2scl, &
         SSPscl,cme_B1_dim,cme_alpha,ModulationRho,        &
         ModulationP,OrientationGL98,LatitudeGL98,LongitudeGL98
    use ModIO,             ONLY: iUnitOut, write_prefix
    implicit none
    !  
    real, dimension(3), intent(in):: R_GL98_D
    real, dimension(3), intent(out):: B_GL98_D
    real, intent(out):: rho_GL98,p_GL98
    !\
    ! User declared local variables go here::
    !/
    real:: x,y,z,  x_1,y_1,z_1,  x_2,y_2,z_2
    real:: r,cos_theta,sin_theta,cos_phi,sin_phi
    real:: r_1,cos_theta1,sin_theta1,cos_phi1,sin_phi1,lambda
    real:: r_2,cos_theta2,sin_theta2,cos_phi2,sin_phi2
    real:: dr2dr1,dth2dr1,cos_thmax,sin_thmax,dsin_thmaxdr
    real:: Br,Btheta,Bphi
    real:: Br1,Btheta1,Bphi1 
    real:: Br2,Btheta2,Bphi2
    real:: Bx_1,By_1,Bz_1
    real:: Bx_2,By_2,Bz_2
    real:: Br_r0,Btheta_r0
    real:: dBr1dr,dBtheta1dr,dBphi1dr
    real:: dBr_r0dr,dBtheta_r0dr 
    real:: dBr2dr1,dBtheta2dr1,dBphi2dr1
    real:: dBr2dr2,dBth2dr2,dBphi2dr2
    real:: dBr2dth2,dBth2dth2,dBphi2dth2
    real:: A2,dA2dr,dA2dth, d2A2dr2,d2A2drdth,d2A2dth2
    real:: pres_1,dpres_1dr1,F_grav,alpha0,ga0r0,delta 
    real, dimension(3):: R1_GL98_D,B1_GL98_D
    real, dimension(3,3):: RotateGL98_DD
    !------------------------------------------------------------------------
    if (DoFirst_GL) then

       DoFirst_GL=.false.
       !\
       ! Construct the rotational matrix RotateGL98_DD,
       !/
       RotateGL98_DD  = matmul( &
            rot_matrix_z(OrientationGL98*cDegToRad),&
            rot_matrix_y((LatitudeGL98-90)*cDegToRad))
       RotateGL98_DD = matmul(RotateGL98_DD, &
            rot_matrix_z(-LongitudeGL98*cDegToRad))

       if(iProc==0)then
          call write_prefix; write(iUnitOut,*) ''
          call write_prefix; write(iUnitOut,*) &
               '>>>>>>>>>>>>>>>>>>>                  <<<<<<<<<<<<<<<<<<<<<'
          call write_prefix; write(iUnitOut,*) &
               '            Initial Perturbation Is Initiated!!!'
          call write_prefix; write(iUnitOut,*) &
               '>>>>>>>>>>>>>>>>>>>                  <<<<<<<<<<<<<<<<<<<<<'
          call write_prefix; write(iUnitOut,*) ''
          call write_prefix; write(iUnitOut,*) 'B1_dim = ',cme_B1_dim
          call write_prefix; write(iUnitOut,*) 'cme_a  = ',cme_a
          call write_prefix; write(iUnitOut,*) 'cme_r1 = ',cme_r1
          call write_prefix; write(iUnitOut,*) 'cme_r0 = ',cme_r0
          call write_prefix; write(iUnitOut,*) 'cme_a1 = ',cme_a1
          call write_prefix; write(iUnitOut,*) 'ModulationRho = ',ModulationRho
          call write_prefix; write(iUnitOut,*) 'ModulationP   = ',ModulationP
       end if
    end if
    !
    delta = cOne/(cTwo*(cOne+cFour))
    !\
    ! Compute R1_GL98_D::
    !/
    R1_GL98_D = matmul(RotateGL98_DD, R_GL98_D)
    !\
    ! CALCULATE CELL CENTER FOR GLOBAL CARTESIAN COORDINATES 
    !/
    x = R1_GL98_D(x_)
    y = R1_GL98_D(y_)
    z = R1_GL98_D(z_)
    !\
    ! CALCULATE CELL CENTER FOR GLOBAL SPHERICAL COORDINATES 
    !/
    r = sqrt(x**2 + y**2 + z**2)
    cos_theta = z/r
    sin_theta = sqrt(x**2 + y**2)/r
    cos_phi   = x/sqrt(x**2 + y**2)
    sin_phi   = y/sqrt(x**2 + y**2)
    if (r <= delta) then 
       r = delta
       x = delta*sin_theta*cos_phi
       y = delta*sin_theta*sin_phi
       z = delta*cos_theta
    end if
    !\
    ! CALCULATE CELL CENTER FOR TRANSFORMED SPHERICAL COORDINATES 
    ! stretching transformation of variables r --> r - a
    !/
    lambda = r + cme_a
    r_1    = lambda
    cos_theta1 = cos_theta
    sin_theta1 = sin_theta 
    cos_phi1   = cos_phi
    sin_phi1   = sin_phi
    !\
    ! CALCULATE CELL CENTER FOR TRANSFORMED CARTESIAN COORDINATES 
    !/
    x_1 = lambda*sin_theta1*cos_phi1
    y_1 = lambda*sin_theta1*sin_phi1
    z_1 = lambda*cos_theta1
    !---------------------------FLUX ROPE REGION---------------------
    ! CALCULATE CELL CENTER CARTESIAN COORDINATES for CME FLUX ROPE
    ! stretching transformation r = r --> r - a
    !----------------------------------------------------------------
    x_2 = x_1
    y_2 = y_1 
    z_2 = z_1 - cme_r1 
    !\
    ! CALCULATE CELL CENTER SPHERICAL COORDINATES for CME FLUX ROPE
    !/
    r_2 = sqrt(x_2**2 + y_2**2 + z_2**2)
    cos_theta2 = x_2/r_2
    sin_theta2 = sqrt(z_2**2 + y_2**2)/r_2
    cos_phi2   = y_2/sqrt(z_2**2 + y_2**2)
    sin_phi2   = z_2/sqrt(z_2**2 + y_2**2)
    if (r_2 <= delta) then
       r_2 = delta
       y_2 = delta*sin_theta2*cos_phi2
       z_2 = delta*sin_theta2*sin_phi2
       x_2 = delta*cos_theta2
    end if
    alpha0 = 5.763854/cme_r0
    ga0r0 = sin(alpha0*cme_r0)/(alpha0*cme_r0) - cos(alpha0*cme_r0)
    A2 = (cFour*cPi*a1scl/alpha0**2)*((cme_r0**2/ga0r0)             * &
         (sin(alpha0*r_2)/(alpha0*r_2) - cos(alpha0*r_2)) - r_2**2) * &
         sin_theta2**2
    dA2dr = ((cFour*cPi*a1scl/alpha0**2)*((cme_r0**2/ga0r0)         * & 
         (cos(alpha0*r_2)/r_2 - sin(alpha0*r_2)/(alpha0*r_2**2)     + & 
         alpha0*sin(alpha0*r_2)) - cTwo*r_2))*sin_theta2**2
    dA2dth = (8.0*cPi*a1scl/alpha0**2)*((cme_r0**2/ga0r0)           * & 
         (sin(alpha0*r_2)/(alpha0*r_2) - cos(alpha0*r_2)) - r_2**2) * &
         sin_theta2*cos_theta2 
    d2A2dr2 = (cFour*cPi*a1scl/alpha0**2)*sin_theta2**2             * &
         ( (cme_r0**2/ga0r0)*(cTwo*sin(alpha0*r_2)/(alpha0*r_2**3)  - &
         cTwo*cos(alpha0*r_2)/(r_2**2) - alpha0*sin(alpha0*r_2)/r_2 + &
         (alpha0**2)*cos(alpha0*r_2)) - cTwo)  
    d2A2drdth = (8.0*cPi*a1scl/alpha0**2)*sin_theta2*cos_theta2     * &
         ((cme_r0**2/ga0r0)*(cos(alpha0*r_2)/r_2                    - & 
         sin(alpha0*r_2)/(alpha0*r_2**2)                            + &
         alpha0*sin(alpha0*r_2)) - cTwo*r_2) 
    d2A2dth2 = (8.0*cPi*a1scl/alpha0**2)*((cme_r0**2/ga0r0)         * &
         (sin(alpha0*r_2)/(alpha0*r_2) - cos(alpha0*r_2)) - r_2**2) * &
         (cos_theta2**2 - sin_theta2**2)
    dr2dr1  =  (x_2/r_2)*sin_theta1*cos_phi1                        + &
         (y_2/r_2)*sin_theta1*sin_phi1                              + &
         (z_2/r_2)*cos_theta1
    dth2dr1 =  (cOne/r_2)*(-sin_theta2*sin_theta1*cos_phi1          + &
         cos_theta2*cos_phi2*sin_theta1*sin_phi1                    + &
         cos_theta2*sin_phi2*cos_theta1)
    !\
    ! Derivatives of field components in flux rope spherical coordinates
    !/
    dBr2dr2 = -cTwo*dA2dth/(sin_theta2*r_2**3)           + &
         d2A2drdth/(sin_theta2*r_2**2)
    dBr2dth2 = -cos_theta2*dA2dth/(r_2**2*sin_theta2**2) + & 
         d2A2dth2/(sin_theta2*r_2**2)
    dBth2dr2 = dA2dr/(sin_theta2*r_2**2)                 - &
         d2A2dr2/(r_2 *sin_theta2)
    dBth2dth2 = cos_theta2*dA2dr/(r_2*sin_theta2**2)     - &
         d2A2drdth/(r_2*sin_theta2)
    dBphi2dr2 = alpha0*dA2dr/(r_2*sin_theta2)            - &
         alpha0*A2 /(sin_theta2*r_2**2)
    dBphi2dth2 = alpha0*dA2dth/(r_2*sin_theta2)          - & 
         alpha0*cos_theta2*A2/(r_2*sin_theta2**2)
    !\
    ! Total derivative of the flux rope field components in terms of `r1'
    !/
    dBr2dr1     = dBr2dr2  *dr2dr1 + dBr2dth2  *dth2dr1
    dBtheta2dr1 = dBth2dr2 *dr2dr1 + dBth2dth2 *dth2dr1
    dBphi2dr1   = dBphi2dr2*dr2dr1 + dBphi2dth2*dth2dr1
    !\
    ! Magnetic field components in the flux rope spherical coordinates
    !/
    Br2     = dA2dth/(sin_theta2*r_2**2)
    Btheta2 = -dA2dr/(sin_theta2*r_2)
    Bphi2   = alpha0*A2/(sin_theta2*r_2)
    !\
    ! Magnetic field components in the second cartesian coordinates
    ! X-COMPONENT OF MAGNETIC FIELD
    !/
    Bx_2 = Br2*cos_theta2            - &
         Btheta2*sin_theta2
    ! Y-COMPONENT OF MAGNETIC FIELD
    By_2 = Br2*sin_theta2*cos_phi2   + &
         Btheta2*cos_theta2*cos_phi2 - &
         Bphi2*sin_phi2
    ! Z-COMPONENT OF MAGNETIC FIELD
    Bz_2 = Br2*sin_theta2*sin_phi2   + &
         Btheta2*cos_theta2*sin_phi2 + &
         Bphi2*cos_phi2
    !\
    ! Define the magnetic field in the global cartesian coordinates
    ! INSIDE THE MAGNETIC FLUX ROPE REGION
    !/
    if (sqrt(x_2**2 + y_2**2 + z_2**2) <= cme_r0) then
       Bx_1 = Bx_2 
       By_1 = By_2 
       Bz_1 = Bz_2 
       !\
       ! Magnetic field components in global sperical coordinates
       !/
       Br1     = Bx_1*sin_theta1*cos_phi1 + &
            By_1*sin_theta1*sin_phi1      + &
            Bz_1*cos_theta1
       Btheta1 = Bx_1*cos_theta1*cos_phi1 + &
            By_1*cos_theta1*sin_phi1      - &
            Bz_1*sin_theta1
       Bphi1   = Bx_1*(-cOne)*sin_phi1    + &
            By_1*cos_phi1
       !\
       ! Compute kinetic gas pressure
       !/
       pres_1     = inv_g*rho1scl*SSPscl**2 + a1scl*A2
       dpres_1dr1 = a1scl*(dA2dr*dr2dr1 + dA2dth*dth2dr1)
       !\
       ! MAGNETIC FIELD transformed with stretching transformation
       !/
       Br     = Br1    *(lambda/r)**2
       Btheta = Btheta1*(lambda/r)
       Bphi   = Bphi1  *(lambda/r)
       !\
       ! Magnetic field components in global cartesian coordinates
       !/
       B1_GL98_D(x_) = Br*sin_theta1*cos_phi1   + &
            Btheta*cos_theta1*cos_phi1    - &
            Bphi*sin_phi1
       !\
       ! Y-COMPONENT OF MAGNETIC FIELD:: (x,y,z)_BATSRUS -> (x,y,z)
       !/
       B1_GL98_D(y_) = Br*sin_theta1*sin_phi1   + &
            Btheta*cos_theta1*sin_phi1    + &
            Bphi*cos_phi1
       !\
       ! Z-COMPONENT OF MAGNETIC FIELD:: (x,y,z)_BATSRUS -> (x,y,z)
       !/
       B1_GL98_D(z_) = Br*cos_theta1-Btheta*sin_theta1
       !\
       ! Transform back to the original coordinates
       ! given by R_GL98_D:: 
       !/
       B_GL98_D  = matmul(B1_GL98_D, RotateGL98_DD)
       !\
       ! PLASMA DENSITY with stretching transformation
       !/
       F_grav   = (abs(Gbody)/(r**2) + cme_alpha*r)
       rho_GL98 = (cOne/ F_grav)*(((lambda/r)**2)                   * &
            ((lambda/r)**2 - cOne)*(dpres_1dr1 + (cOne/(cFour*cPi)) * &
            (Br2*dBr2dr1 + Btheta2*dBtheta2dr1 + Bphi2*dBphi2dr1))  + &
            cTwo*lambda*cme_a*pres_1/r**3 + cme_a*lambda            / &
            (cFour*cPi*r**3)*(cOne - cTwo*(lambda/r)**2)*Br1**2     + &
            ((lambda/r)**2)*((cme_a/r)**2 +cTwo*cme_a/r)*(Btheta1**2+ &
            Bphi1**2)/(cFour*cPi*lambda))
       !\
       ! Add background density
       !/
       rho_GL98 = rho_GL98 + rho2scl/r**3
       !\
       ! PLASMA PRESSURE with contraction transformation
       !/
       p_GL98 = pres_1*(lambda/r)**2                               - &
            (cOne/(cTwo*cFour*cPi))*((lambda/r)**2)*((lambda/r)**2 - &
            cOne)*Br1**2 
       !\
       ! Add background pressure
       !/
       p_GL98 = p_GL98 + abs(Gbody)*rho2scl/(cFour*r**4)
    else
       B_GL98_D = cZero; rho_GL98 = cZero; p_GL98 = cZero
    endif

  end subroutine add_GL98_fluxrope

  !========================================================================

  subroutine post_init_stat
    use ModAdvance,    ONLY: State_VGB,rho_,P_
    use ModProcMH,     ONLY: iProc
    use ModPhysics,    ONLY: Gbody
    use ModIO,         ONLY: iUnitOut, write_prefix
    implicit none
    real, external :: maxval_blk, minval_blk
    real :: pMin, pMax, RhoMin, RhoMax
    !---------------------------------------------------------------------
    !\
    ! Post-initialization statistics::
    !/
    !  call write_prefix; write(iUnitOut,*) 'Mass in the flux rope on processor::',&
    !       iProc,Mrope_GL98
    if (iProc==0) then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) &
            '>>>>>>>>>>>>>>>>>>> Pressure and Density Log <<<<<<<<<<<<<<<<<<<<<'
       call write_prefix; write(iUnitOut,*) 'At PE=0'
       call write_prefix; write(iUnitOut,*) 'The value of MaxB0_1 is :: ',MaxB0_1
       call write_prefix; write(iUnitOut,*) 'The value of MaxB0_2 is :: ',MaxB0_2
       call write_prefix; write(iUnitOut,*) 'The value of Bnot    is :: ',Bnot
       pMin = minval_BLK(1,State_VGB(P_,:,:,:,:))
       pMax = maxval_BLK(1,State_VGB(P_,:,:,:,:))
       call write_prefix; write(iUnitOut,*) 'The min,max P is        :: ',&
            pMin, pMax
       RhoMin = minval_BLK(1,State_VGB(Rho_,:,:,:,:))
       RhoMax = maxval_BLK(1,State_VGB(Rho_,:,:,:,:))
       call write_prefix; write(iUnitOut,*) 'The min,max Rho is      :: ',&
            RhoMin, RhoMax
       call write_prefix; write(iUnitOut,*) 'The value of Tnot  is   :: ',Tnot
       call write_prefix; write(iUnitOut,*) 'The value of Gbody is   :: ',Gbody
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) &
            '>>>>>>>>>>>>>>>>>>>                          <<<<<<<<<<<<<<<<<<<<<'
       call write_prefix; write(iUnitOut,*) ''
    end if
  end subroutine post_init_stat

  !========================================================================
  subroutine get_atmosphere_orig(i,j,k,iBLK,Dens_BLK,Pres_BLK)
    !
    !---------------------------------------------------------------------------
    !
    ! This module computes the background atmosphere in the 
    ! presence of gravity for an isothermal plasma::
    ! The subroutine is written by ILR on Feb 3, 2003.
    !
    !---------------------------------------------------------------------------
    !
    use ModGeometry,   ONLY: x_BLK,y_BLK,z_BLK
    use ModNumConst,   ONLY: cHalf,cOne
    use ModPhysics,    ONLY: inv_g
    implicit none

    integer, intent(in):: i,j,k,iBLK
    real, intent(out):: Dens_BLK,Pres_BLK
    real, parameter:: Rho0=cOne
    real:: xx,yy,zz,RR,BBr
    !\
    ! Get the coordinates and radial distance from the Sun::
    !/
    xx = x_BLK(i,j,k,iBLK)
    yy = y_BLK(i,j,k,iBLK)
    zz = z_BLK(i,j,k,iBLK)
    RR = sqrt(xx**2+yy**2+zz**2)
    if (RR > cHalf) then
       Dens_BLK = Rho0/RR**2
       Pres_BLK = Rho0*inv_g/RR**2
    else
       Dens_BLK = Rho0/cHalf**2
       Pres_BLK = Rho0*inv_g/cHalf**2
    endif

  end subroutine get_atmosphere_orig
  !=================================================================
  ! SUBROUTINE get_hlcmm
  ! Read H(eliographic) L(ongitude) of the C(entral) M(eridian) of 
  ! the M(ap). Assign Phi_Shift=HLCMM-180
  subroutine get_hlcmm
    use ModProcMH,ONLY:iProc
    use ModIO,       ONLY: iUnitOut, write_prefix
    use CON_axes,    ONLY: dLongitudeHgrDeg
    use ModNumConst
    implicit none
    real::HLCMM     !Heliographic Longitudee of the central meridian of map
    integer::iHLCMM !The same, but HLCMM is integer at WSO magnetograms
    integer::iErrorRead,iPosition
    iPosition=index(Head_PFSSM,'Centered')	
    if (iPosition>0.and.Phi_Shift<cZero)then	
       Head_PFSSM(1:len(Head_PFSSM)-iPosition)=&
            Head_PFSSM(iPosition+1:len(Head_PFSSM))
       iPosition=index(Head_PFSSM,':')
       Head_PFSSM(1:len(Head_PFSSM)-iPosition)=&
            Head_PFSSM(iPosition+1:len(Head_PFSSM))
       iPosition=index(Head_PFSSM,':')
       Head_PFSSM(1:len(Head_PFSSM)-iPosition)=&
            Head_PFSSM(iPosition+1:len(Head_PFSSM))
       read(Head_PFSSM,'(i3)',iostat=iErrorRead)iHLCMM
       if(iErrorRead>0)call stop_mpi(&
            'Can nod find HLCMM, '//File_PFSSM//&
            ' is not a true WSO magnetogram')
       Phi_Shift=modulo(iHLCMM-180-dLongitudeHgrDeg, 360.0) 
       if(iProc==0)then
          call write_prefix;write(iUnitOut,*)'Phi_Shift=',Phi_Shift
       end if
       return
    end if
    iPosition=index(Head_PFSSM,'Central')
    if(iPosition>0.and.Phi_Shift<cZero)then
       Head_PFSSM(1:len(Head_PFSSM)-iPosition)=&
            Head_PFSSM(iPosition+1:len(Head_PFSSM))
       iPosition=index(Head_PFSSM,':')
       Head_PFSSM(1:len(Head_PFSSM)-iPosition)=&
            Head_PFSSM(iPosition+1:len(Head_PFSSM))
       read(Head_PFSSM,*,iostat=iErrorRead)HLCMM
       if(iErrorRead>0)call stop_mpi(&
            'Can nod find HLCMM, '//File_PFSSM//&
            ' is not a true MDI magnetogram')
       Phi_Shift=modulo(HLCMM-180-dLongitudeHgrDeg, 360.0) 
       if(iProc==0)then
          call write_prefix;write(iUnitOut,*)'Phi_Shift=',Phi_Shift
       end if
    end if
  end subroutine get_hlcmm

  !============================================================================
  subroutine user_get_b0(xx,yy,zz,B0_PFSSM)
    !
    !---------------------------------------------------------------------------
    ! This subroutine computes PFSS (Potential Field Source Surface)
    ! field model components in spherical coordinates at user-specified
    ! r (solar radii units, r>1) and theta, phi (both radians).
    ! The subroutine requires a file containing the spherical
    ! harmonic coefficients g(n,m) & h(n.m) obtained from a separate analysis 
    ! of a photospheric synoptic map, in a standard ("radial")
    ! format and normalization used by Stanford.
    !
    ! The PFSS field model assumes no currents in the corona and
    ! a pure radial field at the source surface, here R=2.5 Rsun
    !
    ! Get solar coefficients from Todd Hoeksema's files:
    !    1. Go to http://solar.stanford.edu/~wso/forms/prgs.html
    !    2. Fill in name and email as required
    !    3. Chose Carrington rotation (leave default 180 center longitude)
    ! For most requests of integer CRs with order < 20, result will come back
    ! immediately on the web.
    !    4. Count header lines before 1st (0,0) coefficient -this will be asked!
    !---------------------------------------------------------------------------
    ! Notes:
    !
    ! In the calling routine you must initialize one variable: istart=0 (it is a 
    ! flag used to tell the subroutine to read the coefficient file the first 
    ! time only). The first time around (DoFirst=0), the subroutine will ask for
    ! the coefficient file name, the order of the expansion to use (N_PFSSM=40 or 
    ! less*, but the coeff file can contain more orders than you use), and the 
    ! number of lines in the coefficient file header. (*note computation time 
    ! increases greatly with order used).
    !
    ! The source surface surface radius has been set at Rs=2.5*Ro in the 
    ! subroutine. PFSS fields at R>Rs are radial.(br,bthet,bphi) are the resulting
    ! components. Note the units of the B fields will differ with observatory used
    ! for the coefficients. Here we assume use of the wso coefficients so units are
    ! microT. The computation of the B fields is taken mainly from Altschuler, 
    ! Levine, Stix, and Harvey, "High Resolutin Mapping of the Magnetic Field of
    ! the Solar Corona," Solar Physics 51 (1977) pp. 345-375. It uses Schmidt
    ! normalized Legendre polynomials and the normalization is explained in the 
    ! paper. The field expansion in terms of the Schmidt normalized Pnm and dPnm's
    ! is best taken from Todd Hoeksema's notes which can be downloaded from the Web
    ! http://quake.stanford.edu/~wso/Description.ps
    ! The expansions  used to get include radial factors to make the field become
    ! purely radial at the source surface. The expans. in Altschuler et al assumed
    ! that the the coefficient g(n,m) and h(n,m) were the LOS coefficients -- but 
    ! the g(n,m) and h(n,m) now available are radial (according to Janet Luhman). 
    ! Therefore, she performs an initial correction to the g(n,m) and h(n,m) to 
    ! make them consistent with the the expansion. There is no reference for this
    ! correction.
    !---------------------------------------------------------------------------

    use ModMain,     ONLY: Time_Simulation, TypeCoordSystem
    use ModNumConst, ONLY: cZero,cHalf,cOne,cTwo,cThree, &
         cE1,cE9,cTolerance,cTiny,cDegToRad,cPi
    use ModProcMH,   ONLY: iProc
    use ModPhysics,  ONLY: No2Io_V, UnitB_, UnitT_, OmegaBody
    use ModIO,       ONLY: iUnitOut, write_prefix
    use ModIoUnit,   ONLY: io_unit_new
    use CON_axes,    ONLY: transform_matrix, XyzPlanetHgi_D

    implicit none

    real, intent(in):: xx,yy,zz
    real, intent(out), dimension(3):: B0_PFSSM
    integer:: iError
    integer:: i,n,m
    real:: gtemp,htemp
    real:: c_n
    real:: SinPhi,CosPhi
    real:: sinmPhi,cosmPhi
    real:: CosTheta,SinTheta
    real:: stuff1,stuff2,stuff3
    real:: SumR,SumT,SumP,SumPsi
    real:: Rin_PFSSM,Theta_PFSSM,Phi_PFSSM
    real:: Br_PFSSM,Btheta_PFSSM,Bphi_PFSSM,Psi_PFSSM
    real, dimension(N_PFSSM+1,N_PFSSM+1):: p_nm,dp_nm

    integer :: iUnit

    !\
    ! Optimization by G. Toth::
    !/
    integer, parameter:: MaxInt=10000
    real, save:: Sqrt_I(MaxInt)
    real:: Xy
    real:: SinThetaM, SinThetaM1
    integer:: delta_m0
    real, dimension(-1:N_PFSSM+2):: RoRsPower_I, RoRPower_I, rRsPower_I
    !--------------------------------------------------------------------------
    !\
    ! Calculate cell-centered spherical coordinates::
    !/
    Rin_PFSSM   = sqrt(xx**2+yy**2+zz**2)
    !\
    ! Avoid calculating B0 inside a critical radius = 0.5*Rsun
    !/
    if (Rin_PFSSM < 9.00E-01) then
       B0_PFSSM = cZero
       RETURN
    end if
    Theta_PFSSM = acos(zz/Rin_PFSSM)
    Xy          = sqrt(xx**2+yy**2+cTiny**2)
    Phi_PFSSM   = atan2(yy,xx)
    SinTheta    = Xy/Rin_PFSSM
    CosTheta    = zz/Rin_PFSSM
    SinPhi      = yy/Xy
    CosPhi      = xx/Xy
    !\
    ! Set the source surface radius::
    ! The inner boundary in the simulations starts at a height
    ! H_PFSSM above that of the magnetic field measurements!
    !/
    if ((Rin_PFSSM+H_PFSSM) > Rs_PFSSM) then 
       R_PFSSM = Rs_PFSSM
    else
       R_PFSSM = Rin_PFSSM+H_PFSSM 
    endif

    if (DoFirst) then
       !\
       ! Initialize once g(n+1,m+1) & h(n+1,m+1) by reading a file
       ! created from Web data::
       !/ 
       DoFirst=.false.
       if (iProc==0) then
          call write_prefix; write(iUnitOut,*) 'Norder = ',N_PFSSM
          call write_prefix; write(iUnitOut,*) 'Entered coefficient file name :: ',File_PFSSM
          call write_prefix; write(iUnitOut,*) 'Entered number of header lines:: ',iHead_PFSSM
       endif
       !\
       ! Formats adjusted for wso CR rad coeffs::
       !/
       iUnit = io_unit_new()
       open(iUnit,file=File_PFSSM,status='old',iostat=iError)
       if (iHead_PFSSM /= 0) then
          do i=1,iHead_PFSSM
             read(iUnit,'(a)') Head_PFSSM
             if(Phi_Shift<-cTiny)call get_hlcmm	
          enddo
       endif
       !\
       ! Initialize all coefficient arrays::
       !/
       g_nm(:,:) = cZero; h_nm(:,:)  = cZero
       p_nm(:,:) = cZero; dp_nm(:,:) = cZero
       !\
       ! Read file with coefficients, g_nm and h_nm::
       !/
       do
          read(iUnit,*,iostat=iError) n,m,gtemp,htemp
          if (iError /= 0) EXIT
          if (n > N_PFSSM .or. m > N_PFSSM) CYCLE
          g_nm(n+1,m+1) = gtemp
          h_nm(n+1,m+1) = htemp
       enddo
       close(iUnit)
       !\
       ! Add correction factor for radial, not LOS, coefficients::
       ! Note old "coefficients" file are LOS, all new coeffs and 
       ! files are radial)
       !/
       do n=0,N_PFSSM
          stuff1 = cOne/real(n+1+(n/(Rs_PFSSM**(2*n+1))))
          do m=0,n
             g_nm(n+1,m+1) = g_nm(n+1,m+1)*stuff1
             h_nm(n+1,m+1) = h_nm(n+1,m+1)*stuff1
          enddo
       enddo
       !\
       ! Calculate sqrt(integer) from 1 to 10000::
       !/
       do m=1,MaxInt
          Sqrt_I(m) = sqrt(real(m))
       end do
       !\
       ! Calculate the ratio sqrt(2m!)/(2^m*m!)::
       !/
       factRatio1(:) = cZero; factRatio1(1) = cOne
       do m=1,N_PFSSM
          factRatio1(m+1) = factRatio1(m)*Sqrt_I(2*m-1)/Sqrt_I(2*m)
       enddo

    end if
    !\
    ! Transform Phi_PFSSM from the component's frame to the magnetogram's frame.
    !/
    Phi_PFSSM = Phi_PFSSM - Phi_Shift*cDegToRad

    !\
    ! Calculate powers of the ratios of radii
    !/
    rRsPower_I(-1) = Rs_PFSSM/R_PFSSM ! This one can have negative power.
    rRsPower_I(0)  = cOne
    RoRsPower_I(0) = cOne
    RoRPower_I(0)  = cOne
    do m=1,N_PFSSM+2
       RoRsPower_I(m) = RoRsPower_I(m-1) * (Ro_PFSSM/Rs_PFSSM)
       RoRPower_I(m)  = RoRPower_I(m-1)  * (Ro_PFSSM/R_PFSSM)
       rRsPower_I(m)  = rRsPower_I(m-1)  * (R_PFSSM /Rs_PFSSM)
    end do
    !\
    ! Calculate polynomials with appropriate normalization
    ! for Theta_PFSSMa::
    !/
    SinThetaM  = cOne
    SinThetaM1 = cOne

    do m=0,N_PFSSM
       if (m == 0) then
          delta_m0 = 1
       else
          delta_m0 = 0
       endif
       !\
       ! Eq.(27) from Altschuler et al. 1976::
       !/
       p_nm(m+1,m+1) = factRatio1(m+1)*Sqrt_I((2-delta_m0)*(2*m+1))* &
            SinThetaM
       !\
       ! Eq.(28) from Altschuler et al. 1976::
       !/
       if (m < N_PFSSM) p_nm(m+2,m+1) = p_nm(m+1,m+1)*Sqrt_I(2*m+3)* &
            CosTheta
       !\
       ! Eq.(30) from Altschuler et al. 1976::
       !/
       dp_nm(m+1,m+1) = factRatio1(m+1)*Sqrt_I((2-delta_m0)*(2*m+1))*&
            m*CosTheta*SinThetaM1
       !\
       ! Eq.(31) from Altschuler et al. 1976::
       !/
       if (m < N_PFSSM) &
            dp_nm(m+2,m+1) = Sqrt_I(2*m+3)*(CosTheta*&
            dp_nm(m+1,m+1)-SinTheta*p_nm(m+1,m+1))

       SinThetaM1 = SinThetaM
       SinThetaM  = SinThetaM*SinTheta

    enddo
    do m=0,N_PFSSM-2; do n=m+2,N_PFSSM
       !\
       ! Eq.(29) from Altschuler et al. 1976::
       !/
       stuff1         = Sqrt_I(2*n+1)/Sqrt_I(n**2-m**2)
       stuff2         = Sqrt_I(2*n-1)
       stuff3         = Sqrt_I((n-1)**2-m**2)/Sqrt_I(2*n-3)
       p_nm(n+1,m+1)  = stuff1*(stuff2*CosTheta*p_nm(n,m+1)-  &
            stuff3*p_nm(n-1,m+1))
       !\
       ! Eq.(32) from Altschuler et al. 1976::
       !/
       dp_nm(n+1,m+1) = stuff1*(stuff2*(CosTheta*dp_nm(n,m+1)-&
            SinTheta*p_nm(n,m+1))-stuff3*dp_nm(n-1,m+1))
    enddo; enddo
    !\
    ! Apply Schmidt normalization::
    !/
    do m=0,N_PFSSM; do n=m,N_PFSSM
       !\
       ! Eq.(33) from Altschuler et al. 1976::
       !/
       stuff1 = cOne/Sqrt_I(2*n+1)
       !\
       ! Eq.(34) from Altschuler et al. 1976::
       !/
       p_nm(n+1,m+1)  = p_nm(n+1,m+1)*stuff1
       dp_nm(n+1,m+1) = dp_nm(n+1,m+1)*stuff1
    enddo; enddo
    !\
    ! Truncate the value of SinTheta::
    !/
    if (SinTheta == cZero) SinTheta = cOne/(cE9*cE1)
    !\
    ! Initialize the values of SumR,SumT,SumP, and SumPsi::
    !/
    SumR = cZero; SumT   = cZero
    SumP = cZero; SumPsi = cZero
    !\
    ! Leave out monopole (n=0) term::
    !/
    g_nm(1,1) = cZero
    !\
    ! Calculate B for (R_PFSSM,Phi_PFSSM)::
    ! Also calculate magnetic potential Psi_PFSSM
    !/
    do m=0,N_PFSSM
       cosmPhi  = cos(m*Phi_PFSSM)
       sinmPhi  = sin(m*Phi_PFSSM)
       do n=m,N_PFSSM
          !\
          ! c_n corresponds to Todd's c_l::
          !/
          c_n    = -RoRsPower_I(n+2)
          !\
          ! Br_PFSSM = -d(Psi_PFSSM)/dR_PFSSM::
          !/
          stuff1 = (n+1)*RoRPower_I(n+2)-c_n*n*rRsPower_I(n-1)
          stuff2 = g_nm(n+1,m+1)*cosmPhi+h_nm(n+1,m+1)*sinmPhi
          SumR   = SumR + p_nm(n+1,m+1)*stuff1*stuff2
          !\
          ! Bt_PFSSM = -(1/R_PFSSM)*d(Psi_PFSSM)/dTheta_PFSSM::
          !/
          stuff1 = RoRPower_I(n+2)+c_n*rRsPower_I(n-1)
          SumT   = SumT-dp_nm(n+1,m+1)*stuff1*stuff2
          !\
          ! Psi_PFSSM::
          !/
          SumPsi = SumPsi+R_PFSSM*p_nm(n+1,m+1)*stuff1*stuff2
          !\
          ! Bp_PFSSM = -(1/R_PFSSM)*d(Psi_PFSSM)/dPhi_PFSSM::
          !/
          stuff2 = g_nm(n+1,m+1)*sinmPhi-h_nm(n+1,m+1)*cosmPhi
          SumP   = SumP + p_nm(n+1,m+1)*m/SinTheta*stuff1*stuff2
       enddo
    enddo
    !\
    ! Compute (Br_PFSSM,Btheta_PFSSM,Bphi_PFSSM) and Psi_PFSSM::
    !/
    Psi_PFSSM    = SumPsi
    Br_PFSSM     = SumR
    if (Rin_PFSSM > Rs_PFSSM) Br_PFSSM = Br_PFSSM*(Rs_PFSSM/Rin_PFSSM)**2
    Btheta_PFSSM = SumT
    Bphi_PFSSM   = SumP
    !\
    ! Magnetic field components in global Cartesian coordinates::
    ! Set B0xCell_BLK::
    !/
    B0_PFSSM(1) = Br_PFSSM*SinTheta*CosPhi+    &
         Btheta_PFSSM*CosTheta*CosPhi-&
         Bphi_PFSSM*SinPhi
    !\
    ! Set B0yCell_BLK::
    !/
    B0_PFSSM(2) = Br_PFSSM*SinTheta*SinPhi+    &
         Btheta_PFSSM*CosTheta*SinPhi+&
         Bphi_PFSSM*CosPhi
    !\
    ! Set B0zCell_BLK::
    !/
    B0_PFSSM(3) = Br_PFSSM*CosTheta-           &
         Btheta_PFSSM*SinTheta
    !\
    ! Apply field strength normalization::
    ! UnitB contains the units of the CR file relative to 1 Gauss
    ! and a possible correction factor (e.g. 1.8 or 1.7).
    !/
    B0_PFSSM(1:3) = B0_PFSSM(1:3)*(UnitB/No2Io_V(UnitB_))

  end subroutine user_get_b0

  !========================================================================
  subroutine user_specify_initial_refinement(iBLK,RefineBlock,lev, &
       dxBlock,xCenter,yCenter,zCenter,rCenter,minx,miny,minz,minR,&
       maxx,maxy,maxz,maxR,IsFound)
    use ModSize,       ONLY: nI,nJ,nK,gcn
    use ModVarIndexes, ONLY: Bx_,By_,Bz_,P_
    use ModAdvance,    ONLY: B0_DGB,State_VGB
    use ModAMR,        ONLY: InitialRefineType
    use ModMain,       ONLY: x_,y_,z_,iteration_number,     &
         time_loop
    use ModGeometry,   ONLY: XyzMin_D,XyzMax_D,x_BLK,y_BLK, &
         z_BLK,R_BLK,dx_BLK,dy_BLK,dz_BLK
    use ModNumConst,   ONLY: cTiny,cHundredth,cEighth,cHalf,&
         cQuarter,cOne,cTwo,cFour,cE1,cE2,cZero
    use ModPhysics,    ONLY: No2Io_V, UnitB_

    logical, intent(out):: RefineBlock,IsFound
    integer, intent(in):: lev
    real, intent(in):: dxBlock
    real, intent(in):: xCenter,yCenter,zCenter,rCenter
    real, intent(in):: minx,miny,minz,minR
    real, intent(in):: maxx,maxy,maxz,maxR
    integer, intent(in):: iBLK

    logical:: DoRefineInitCS=.false.
    logical:: ResolveNullPount=.false.
    logical:: DoCallUserB0
    logical:: IsInRangeAR,IsInRangeCS
    real:: dsMin,dsNewMin
    real:: critx,critvRdotR0,critxCenter
    real:: RminRv,RdotRv,RminRn,RdotR0,R2Cell
    real, dimension(3):: RCell_D,RvCell_D,RnCell_D,R0Cell_D

    integer:: i,j,k
    real, parameter:: &
         XLoc_V=-9.6722621E-01,YLoc_V=-4.2230144E-02,ZLoc_V=-2.5038001E-01
    real, parameter:: &
         XLoc_0=-9.6800017E-01,YLoc_0=-1.6896725E-02,ZLoc_0=-2.5038001E-01
    real:: XCell,YCell,ZCell,RCell,RCentre
    real:: B0xCell,B0yCell,B0zCell
    real:: BIxCell,BIyCell,BIzCell
    real:: DensCell,PresCell,GammaCell
    logical, dimension(3):: IsGhostCell_D
    real, dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn):: Br_D,Beta_D
    !------------------------------------------------------------------------
    dsMin = min(&
         dx_BLK(iBLK),dy_BLK(iBLK),dz_BLK(iBLK))
    RCell_D (x_) = xCenter
    RCell_D (y_) = yCenter
    RCell_D (z_) = zCenter
    RvCell_D(x_) = XLoc_V
    RvCell_D(y_) = YLoc_V
    RvCell_D(z_) = ZLoc_V
    !At 01:35UT on Oct 27, 2003:
    !  RnCell_D(x_) = -1.0156200
    !  RnCell_D(y_) = -0.034179699
    !  RnCell_D(z_) = -0.18133120
    !At 09:35UT on Oct 28, 2003:
    RnCell_D(x_) = -1.083020
    RnCell_D(y_) = -0.106262
    RnCell_D(z_) = -0.220000
    R0Cell_D(x_) = XLoc_0
    R0Cell_D(y_) = YLoc_0
    R0Cell_D(z_) = ZLoc_0
    R2Cell = sqrt(dot_product(RCell_D,RCell_D))
    RdotRv = dot_product(RCell_D,RvCell_D)/&
         max(cTiny,R2Cell)
    RdotR0 = dot_product(RCell_D,R0Cell_D)/&
         max(cTiny,R2Cell)
    RminRv = sqrt(dot_product(&
         RCell_D-RvCell_D,RCell_D-RvCell_D))
    RminRn = sqrt(dot_product(&
         RCell_D-RnCell_D,RCell_D-RnCell_D))
    RefineBlock = .false.
    IsInRangeAR = .false.; IsInRangeCS = .false.
    !\
    ! Construct a logical switch DoRefineInitCS whether or
    ! not to refine the current sheet initially AND the active
    ! region of interest::
    !/
    if (DoRefineInitCS) then
       DoCallUserB0 = (lev>9.and.RminRv<0.19.and.R2Cell>0.98).or.&
            (lev>3.and.lev<6)
    else
       DoCallUserB0 = (lev>9.and.RminRv<0.19.and.R2Cell>0.98).or.&
            !For 01:35UT on Oct 27, 2003: (lev>12.and.RminRn<2.00E-01.and.R2Cell>1.035)
            !For 09:35UT on Oct 28, 2003: (lev>11.and.RminRn<4.00E-01.and.R2Cell>1.07)
       (lev>11.and.RminRn<0.4.and.R2Cell>1.07)
    endif
    !\
    ! For initial refinement, get the radial magnetic field
    ! to refine the desired blocks in the AR and in the CS,
    ! as long as DoCallUserB0 = .true.
    !/
    if (.not.time_loop) then
       if (DoCallUserB0) then
          call set_b0(iBLK)
          do k=1-gcn,nK+gcn
             IsGhostCell_D(3)=k<1.or.k>nK
             do j=1-gcn,nJ+gcn
                IsGhostCell_D(2)=j<1.or.j>nJ           
                do i=1-gcn,nI+gcn
                   IsGhostCell_D(1)=i<1.or.i>nI
                   if (count(IsGhostCell_D)>1) then
                      Br_D(i,j,k)   = huge(cOne)
                      Beta_D(i,j,k) = cZero
                      CYCLE
                   end if
                   call get_plasma_parameters_cell(i,j,k,iBLK,&
                        DensCell,PresCell,GammaCell)
                   XCell = x_BLK(i,j,k,iBLK)
                   YCell = y_BLK(i,j,k,iBLK)
                   ZCell = z_BLK(i,j,k,iBLK)
                   RCell = sqrt(XCell**2+YCell**2+ZCell**2)
                   B0xCell = B0_DGB(x_,i,j,k,iBLK)
                   B0yCell = B0_DGB(y_,i,j,k,iBLK)
                   B0zCell = B0_DGB(z_,i,j,k,iBLK)
                   Br_D(i,j,k) = abs(   &
                        (XCell*B0xCell+ &
                        YCell*B0yCell+ &
                        ZCell*B0zCell)/&
                        RCell)
                   Beta_D(i,j,k) = cTwo*PresCell/max(cTiny,&
                        (B0xCell**2+B0yCell**2+B0zCell**2))
                end do
             end do
          end do
          !\
          ! Construct refinement criteria to refine blocks in the
          ! AR (IsInRangeAR) and in the CS (IsInRangeCS)::
          !/
          if (lev>3.and.lev<6.and.DoRefineInitCS) &
               IsInRangeCS = (minval(Br_D)<1.50E-05)
          if (lev>9) &
               IsInRangeAR = (minval(Br_D)<2.50E+00/No2Io_V(UnitB_))
          !For 01:35UT on Oct 27, 2003: if (lev>12) &
          !For 09:35UT on Oct 28, 2003: if (lev>11) &
          if (lev>11) &
               IsInRangeAR = (minval(Br_D)<2.50E+00/No2Io_V(UnitB_)).and.&
               !For 01:35UT on Oct 27, 2003: (R2Cell<1.020E+00).or.(maxval(Beta_D)>0.01)
               !For 09:35UT on Oct 28, 2003: (R2Cell<1.025E+00).or.(maxval(Beta_D)>0.055)
          (R2Cell<1.025E+00).or.(maxval(Beta_D)>0.055)
       endif
    else
       if (.not.ResolveNullPount) then
          dsNewMin = (cOne+cQuarter)/cFour+cHundredth
          if (dsMin>dsNewMin) then
             do k=1-gcn,nK+gcn
                IsGhostCell_D(3)=k<1.or.k>nK
                do j=1-gcn,nJ+gcn
                   IsGhostCell_D(2)=j<1.or.j>nJ           
                   do i=1-gcn,nI+gcn
                      IsGhostCell_D(1)=i<1.or.i>nI
                      if (count(IsGhostCell_D)>1) then
                         Br_D(i,j,k) = huge(cOne)
                         CYCLE
                      end if
                      XCell = x_BLK(i,j,k,iBLK)
                      YCell = y_BLK(i,j,k,iBLK)
                      ZCell = z_BLK(i,j,k,iBLK)
                      RCell = sqrt(XCell**2+YCell**2+ZCell**2)
                      B0xCell = B0_DGB(x_,i,j,k,iBLK)
                      B0yCell = B0_DGB(y_,i,j,k,iBLK)
                      B0zCell = B0_DGB(z_,i,j,k,iBLK)
                      BIxCell = State_VGB(Bx_,i,j,k,iBLK)
                      BIyCell = State_VGB(By_,i,j,k,iBLK)
                      BIzCell = State_VGB(Bz_,i,j,k,iBLK)
                      Br_D(i,j,k) = abs(             &
                           (XCell*(B0xCell+BIxCell)+ &
                           YCell*(B0yCell+BIyCell)+ &
                           ZCell*(B0zCell+BIzCell))/&
                           RCell)
                   end do
                end do
             end do
             !\
             ! Construct refinement criteria to refine blocks in the
             ! current sheet only (IsInRangeCS = .true.)::
             !/
             RCentre = cEighth*&
                  (R_BLK( 1, 1, 1,iBLK)+R_BLK( 1, 1,nK,iBLK)+&
                  R_BLK( 1,nJ, 1,iBLK)+R_BLK( 1,nJ,nK,iBLK)+&
                  R_BLK(nI, 1, 1,iBLK)+R_BLK(nI, 1,nK,iBLK)+&
                  R_BLK(nI,nJ, 1,iBLK)+R_BLK(nI,nJ,nK,iBLK))
             IsInRangeCS = (RCentre**3*minval(Br_D)<3.50E-01)
          end if
       else
          !For 01:35UT on Oct 27, 2003: if (RminRn<1.50E-01.and.R2Cell>1.035) then
          !For 09:35UT on Oct 28, 2003: if (RminRn<1.50E-01.and.R2Cell>1.060) then
          if (RminRn<1.50E-01.and.R2Cell>1.060) then
             do k=1-gcn,nK+gcn
                IsGhostCell_D(3)=k<1.or.k>nK
                do j=1-gcn,nJ+gcn
                   IsGhostCell_D(2)=j<1.or.j>nJ           
                   do i=1-gcn,nI+gcn
                      IsGhostCell_D(1)=i<1.or.i>nI
                      if (count(IsGhostCell_D)>1) then
                         Beta_D(i,j,k) = cZero
                         CYCLE
                      end if
                      B0xCell  = B0_DGB(x_,i,j,k,iBLK)
                      B0yCell  = B0_DGB(y_,i,j,k,iBLK)
                      B0zCell  = B0_DGB(z_,i,j,k,iBLK)
                      BIxCell  = State_VGB(Bx_   ,i,j,k,iBLK)
                      BIyCell  = State_VGB(By_   ,i,j,k,iBLK)
                      BIzCell  = State_VGB(Bz_   ,i,j,k,iBLK)
                      PresCell = State_VGB(P_,i,j,k,iBLK)
                      Beta_D(i,j,k) = cTwo*PresCell/max(cTiny,&
                           ((B0xCell+BIxCell)**2+ &
                           (B0yCell+BIyCell)**2+ &
                           (B0zCell+BIzCell)**2))
                   end do
                end do
             end do
             !\
             ! Construct refinement criteria to refine blocks in the
             ! current sheet only (IsInRangeCS = .true.)::
             !/
             !For 01:35UT on Oct 27, 2003: IsInRangeCS = (maxval(Beta_D)>0.015)
             !For 09:35UT on Oct 28, 2003: IsInRangeCS = (maxval(Beta_D)>0.045)
             IsInRangeCS = (maxval(Beta_D)>0.045)
          end if
       endif
    end if

    select case (InitialRefineType)
    case('UserHELIO','UserHelio','userhelio')
       if (lev<3) then
          RefineBlock = .true.
       else if (lev<10) then
          critx=(XyzMax_D(1)-XyzMin_D(1))/(2.0**real(lev-2))
          if (rCenter<1.10+critx) then
             RefineBlock = .true.
          end if
       else if (lev<11) then
          RefineBlock = (RminRv<0.29).and.(R2Cell>0.96)
       else if (lev<12) then
          RefineBlock = (RminRv<0.23).and.(R2Cell>0.96)
       else
          RefineBlock = (RminRv<0.17).and.(R2Cell>0.96)
       endif
       IsFound = .true.
    case('UserAR8210')
       RefineBlock =  (RminRv<0.36).and.(R2Cell>0.95)
       IsFound = .true.
    case('UserCME')
       if (lev<4) then
          critxCenter = cHalf*(XyzMax_D(1)-XyzMin_D(1))/cTwo**real(3-lev)
          critvRdotR0 = 0.96
          if ((RdotR0>critvRdotR0).and.(abs(xCenter)>critxCenter)) &
               RefineBlock = .true.
       else
          RefineBlock = .false.
       endif
       IsFound = .true.
    case('UserAR486','userar486','USERAR486')
       if (.not.time_loop) then
          if (lev<4) then
             RefineBlock = .true.
          else if (lev<10) then
             !Block is in the CS or intersects the body::
             RefineBlock = IsInRangeCS.or.(minR<=1.00)
          else 
             !Block is in the AR::
             RefineBlock = IsInRangeAR
          endif
          !Do not refine inside the body::
          RefineBlock = RefineBlock.and.(maxR>1.00)
       else
          !Block is in the CS::
          RefineBlock = IsInRangeCS
       end if
       IsFound = .true.
    endselect
  end subroutine user_specify_initial_refinement
  !========================================================================
  !========================================================================
  !  SUBROUTINE user_amr_criteria
  !========================================================================
  !
  !\
  ! This subroutine allows the user to add a refinement type
  ! based on a geometric criteria or physical criteria.  The `case'
  ! specified in the #AMRCRITERIA file will be read here.
  !/
  subroutine user_amr_criteria(iBLK, userCriteria, TypeCriteria, IsFound)
    use ModMain
    use ModAdvance
    use ModGeometry, ONLY:x_BLK,y_BLK,z_BLK,R_BLK,&
         dx_BLK,dy_BLK,dz_BLK,true_cell
    use ModPhysics
    use ModConst
    !\ 
    ! Variables required by this user subroutine::
    !/
    integer, intent(in):: iBLK
    logical, intent(out):: IsFound
    real, intent(out):: userCriteria
    character (len=20),intent(in):: TypeCriteria
    !\
    ! Local variables::
    !/
    logical:: IsInRange
    integer:: i,j,k
    real:: dsMin,dsMax,dsTwo
    real:: XCell,YCell,ZCell,RCell,RCenter
    real:: B0xCell,B0yCell,B0zCell,MinBr,MaxBr
    real:: BIxCell,BIyCell,BIzCell
    real, dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn):: Br_D
    logical,dimension(3)::IsGhostCell_D
    !\
    ! Find the radial location of the center of the block and
    ! the min/max cell size::
    !/
    RCenter = cEighth*&
         (R_BLK( 1, 1, 1,iBLK)+R_BLK( 1, 1,nK,iBLK)+&
         R_BLK( 1,nJ, 1,iBLK)+R_BLK( 1,nJ,nK,iBLK)+&
         R_BLK(nI, 1, 1,iBLK)+R_BLK(nI, 1,nK,iBLK)+&
         R_BLK(nI,nJ, 1,iBLK)+R_BLK(nI,nJ,nK,iBLK))
    dsMin = min(dx_BLK(iBLK),dy_BLK(iBLK),dz_BLK(iBLK))
    dsMax = max(dx_BLK(iBLK),dy_BLK(iBLK),dz_BLK(iBLK))
    dsTwo = dsMin*dsMax

    select case (TypeCriteria)
    case('UserCS','USERCS','usercs')
       !\
       ! Get the radial magnetic field in iBLK::
       !/
       do k=1-gcn,nK+gcn
          IsGhostCell_D(3)=k<1.or.k>nK
          do j=1-gcn,nJ+gcn
             IsGhostCell_D(2)=j<1.or.j>nJ           
             do i=1-gcn,nI+gcn
                IsGhostCell_D(1)=i<1.or.i>nI
                if (count(IsGhostCell_D)>1) then
                   Br_D(i,j,k) = huge(cOne)
                   CYCLE
                end if
                XCell   = x_BLK(i,j,k,iBLK)
                YCell   = y_BLK(i,j,k,iBLK)
                ZCell   = z_BLK(i,j,k,iBLK)
                RCell   = R_BLK(i,j,k,iBLK)
                B0xCell = B0_DGB(x_,i,j,k,iBLK)
                B0yCell = B0_DGB(y_,i,j,k,iBLK)
                B0zCell = B0_DGB(z_,i,j,k,iBLK)
                BIxCell = State_VGB(Bx_,i,j,k,iBLK)
                BIyCell = State_VGB(By_,i,j,k,iBLK)
                BIzCell = State_VGB(Bz_,i,j,k,iBLK)
                Br_D(i,j,k) = abs(&
                     (XCell*(B0xCell+BIxCell)+ &
                     YCell*(B0yCell+BIyCell)+ &
                     ZCell*(B0zCell+BIzCell))/&
                     RCell)
             end do
          end do
       end do
       !\
       ! Find the minimum of abs(Br) in iBLK::
       !/
       MinBr = minval(Br_D)
       !\
       ! Construct refine criteria based on ds2, RCenter, and MinBr::
       !/
       IsInRange = (RCenter<1.50*Rs_PFSSM).and.&
            (MinBr<3.0E-05) 
       if (IsInRange) then
          userCriteria = dsTwo*RCenter*exp(-MinBr)
       else
          userCriteria = dsTwo*RCenter*exp(-MinBr)/cE6
       end if
       IsFound = .true.
    endselect
  end subroutine user_amr_criteria

  !========================================================================
  subroutine user_update_states(iStage,iBlock)
    use ModVarIndexes, Ew_ => ScalarFirst_
    use ModSize
    use ModAdvance, ONLY: State_VGB
    use ModMain,    ONLY: nStage
    use ModPhysics, ONLY: inv_gm1
    use ModEnergy,  ONLY: calc_energy_cell

    implicit none
    integer,intent(in):: iStage,iBlock
    integer:: i,j,k
    real:: DensCell,PresCell,GammaCell

    call update_states_MHD(iStage,iBlock)
    !\
    ! Begin update of pressure and relaxation energy::
    !/
    !  if (iStage/=nStage) return
    do k=1,nK; do j=1,nJ; do i=1,nI
       call get_plasma_parameters_cell(i,j,k,iBlock,&
            DensCell,PresCell,GammaCell)
       State_VGB(P_   ,i,j,k,iBlock)=           &
            (GammaCell-cOne)*                   &
            (inv_gm1*State_VGB(P_,i,j,k,iBlock)&
            + State_VGB(Ew_,i,j,k,iBlock)& 
       )
       State_VGB(Ew_,i,j,k,iBlock)=       & 
            State_VGB(P_,i,j,k,iBlock)*(cOne/(GammaCell - cOne) - inv_gm1)
    end do; end do; end do
    call calc_energy_cell(iBlock)
    !\
    ! End update of pressure and relaxation energy::
    !/
  end subroutine user_update_states


  !========================================================================
  !  SUBROUTINE user_get_log_var
  !========================================================================
  !
  ! This subroutine allows the user to write to the log files variables
  ! (normalized and not) which are problem specific.
  !
  ! The variables specific to the problem are loaded from ModUser.
  !
  !========================================================================
  subroutine user_get_log_var(VarValue,TypeVar,Radius)
    use ModProcMH,     ONLY: nProc
    use ModIO,         ONLY: dn_output,logfile_,write_myname
    use ModMain,       ONLY: unusedBLK,nBLK,iteration_number,   &
         x_,y_,z_
    use ModVarIndexes, ONLY: Bx_, By_, Bz_, rho_, rhoUx_, rhoUy_, rhoUz_, P_, &
         Ew_ => ScalarFirst_
    use ModGeometry,   ONLY: R_BLK
    use ModAdvance,    ONLY: State_VGB,tmp1_BLK,B0_DGB
    use ModPhysics,    ONLY: inv_gm1,No2Si_V(UnitEnergydens_),No2Si_V(UnitX_),&
         No2Si_V(UnitU_),No2Si_V(UnitRho_)
    use ModNumConst,   ONLY: cOne,cHalf,cE1,cE3,cE6
    real, intent(out):: VarValue
    character (LEN=10), intent(in):: TypeVar 
    real, intent(in), optional :: Radius

    integer:: iBLK
    real:: unit_energy,unit_mass
    real, external:: integrate_BLK
    !--------------------------------------------------------------------------
    unit_energy = cE1*cE6*No2Si_V(UnitEnergydens_)*No2Si_V(UnitX_)**3
    unit_mass   = cE3*No2Si_V(UnitRho_)*No2Si_V(UnitX_)**3
    !\
    ! Define log variable to be saved::
    !/
    select case(TypeVar)
    case('em_t','Em_t','em_r','Em_r')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(:,:,:,iBLK) = &
               (B0_DGB(x_,:,:,:,iBLK)+State_VGB(Bx_,:,:,:,iBLK))**2+&
               (B0_DGB(y_,:,:,:,iBLK)+State_VGB(By_,:,:,:,iBLK))**2+&
               (B0_DGB(z_,:,:,:,iBLK)+State_VGB(Bz_,:,:,:,iBLK))**2
       end do
       VarValue = unit_energy*cHalf*integrate_BLK(1,tmp1_BLK)
    case('ek_t','Ek_t','ek_r','Ek_r')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(:,:,:,iBLK) = &
               (State_VGB(rhoUx_,:,:,:,iBLK)**2 +&
               State_VGB(rhoUy_,:,:,:,iBLK)**2 +&
               State_VGB(rhoUz_,:,:,:,iBLK)**2)/&
               State_VGB(rho_  ,:,:,:,iBLK)             
       end do
       VarValue = unit_energy*cHalf*integrate_BLK(1,tmp1_BLK)
    case('et_t','Et_t','et_r','Et_r')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(:,:,:,iBLK) = State_VGB(P_,:,:,:,iBLK)
       end do
       VarValue = unit_energy*inv_gm1*integrate_BLK(1,tmp1_BLK)
    case('ew_t','Ew_t','ew_r','Ew_r')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(:,:,:,iBLK) = State_VGB(Ew_,:,:,:,iBLK)
       end do
       VarValue = unit_energy*integrate_BLK(1,tmp1_BLK)
    case('ms_t','Ms_t')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(:,:,:,iBLK) = &
               State_VGB(rho_,:,:,:,iBLK)/R_BLK(:,:,:,iBLK)
       end do
       VarValue = unit_mass*integrate_BLK(1,tmp1_BLK)
    case('vol','Vol')
       tmp1_BLK(:,:,:,iBLK) = cOne
       VarValue = integrate_BLK(1,tmp1_BLK)
    case default
       VarValue = -7777.
       call write_myname;
       write(*,*) 'Warning in user_get_log_var: unknown logvarname = ',TypeVar
    end select
  end subroutine user_get_log_var

end module ModUser

