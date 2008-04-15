!^CFG COPYRIGHT UM
!========================================================================
module ModUserTD99
  use ModConst
  use ModMain,       ONLY: x_,y_,z_,nByteReal
  use ModVarIndexes, ONLY: Ux_,Uy_,Uz_,Bx_,By_,Bz_
  use ModPhysics,    ONLY: Si2No_V, No2Si_V,No2Io_V,UnitX_,UnitRho_,UnitU_,UnitB_,UnitJ_

  implicit none
  save
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
    real:: AlphaRope,LInduct,WFRope,FootSepar,ItubeDim
    !--------------------------------------------------------------------
    !\
    ! Compute the magnetic energy, WFRope, associated with the portion
    ! of the flux rope current that is above the solar surface::
    !/
    !--------------------------------------------------------------------
    InvH0_TD99 = cGravitation*Msun/Rsun*Si2No_V(UnitU_)**2       ! in [-]
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
    Rtube_TD99 = Rtube_TD99*Si2No_V(UnitX_)
    atube_TD99 = atube_TD99*Si2No_V(UnitX_)
    ItubeDim   = Itube_TD99
    Itube_TD99 = ItubeDim*Si2No_V(UnitJ_)*Si2No_V(UnitX_)**2 ! in [A]
    Rho0_TD99  = Rho0_TD99*Si2No_V(UnitRho_)
    !--------------------------------------------------------------------
    ! Save the maximum value of the current for possible use in
    ! varied_current case::
    !--------------------------------------------------------------------
    ItubeSaved = Itube_TD99 
    !--------------------------------------------------------------------
    ! Strapping field::
    !--------------------------------------------------------------------
    d_TD99     = d_TD99*Si2No_V(UnitX_)
    L_TD99     = L_TD99*Si2No_V(UnitX_)
    q_TD99     = q_TD99*Si2No_V(UnitB_)*Si2No_V(UnitX_)**2
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
            ItubeDim,'[A]'
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
       WFRope    = LInduct*(ItubeDim)**2/cTwo*cE6*cE1 ! in [ergs]
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
            VTransX*(TimeNow-CurrentStartTime)*Si2No_V(UnitX_)
       RMins_D(x_) = RFace_D(x_)+L_TD99 + &
            VTransX*(TimeNow-CurrentStartTime)*Si2No_V(UnitX_)
       RPlus_D(y_) = RFace_D(y_)        - &
            VTransY*(TimeNow-CurrentStartTime)*Si2No_V(UnitX_)
       RMins_D(y_) = RFace_D(y_)        + &
            VTransY*(TimeNow-CurrentStartTime)*Si2No_V(UnitX_)
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
          UTranC_D(x_) = (VTransX*Si2No_V(UnitU_))*BqZOverBqZ0 * &
               sign(cOne,BqField_D(z_))*exp(-BqZFunction)
          UTranC_D(y_) = (VTransY*Si2No_V(UnitU_))*BqZOverBqZ0 * &
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

    end subroutine calc_elliptic_int_2kind

  end subroutine compute_TD99_FluxRope

end module ModUserTD99

!==============================================================================
!==============================================================================
!==============================================================================

module ModUser
  use ModNumConst, ONLY: cHalf,cTwo,cThree,&
       cFour,cE1,cHundred,cHundredth,cZero,&
       cOne
  use ModMain,     ONLY: UseUserB0
  use ModSize,     ONLY: nI,nJ,nK,gcn,nBLK
  use ModUserTD99  ! To include TD99 flux rope.
  use ModMagnetogram
  use ModExpansionFactors
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_initial_perturbation,       &
       IMPLEMENTED3 => user_face_bcs,                   &
       IMPLEMENTED4 => user_get_log_var,                &
       IMPLEMENTED5 => user_get_b0,                     &
       IMPLEMENTED6 => user_update_states,              &
       IMPLEMENTED7 => user_specify_initial_refinement

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'EMPIRICAL SC - Cohen, Sokolov'
  !\
  ! Parameters related to the empirical heating::
  !/
contains
  !============================================================================
  subroutine user_read_inputs
    use ModMain
    use ModProcMH,    ONLY: iProc
    use ModReadParam
    use ModIO,        ONLY: write_prefix, write_myname, iUnitOut

    integer:: i
    character (len=100) :: NameCommand
    character (len=lStringLine)   :: NameModel
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
       case("#PFSSM")
          call read_var('UseUserB0'  ,UseUserB0)
          if (UseUserB0)then
             call read_magnetogram_file
             call read_var('dt_UpdateB0',dt_UpdateB0)
             DoUpdateB0 = dt_updateb0 > 0.0
          endif
          !       case("#AWHEAT")
          !          call read_var('Bnot        ',Bnot)
          !          call read_var('Tnot        ',Tnot)
          !          call read_var('DegFrm1     ',DegFrm1)
          !          call read_var('DegF_Ratio  ',DegF_Ratio)
          !          call read_var('Dens_Ratio  ',Dens_Ratio)
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
       case("#EMPIRICALSW")
          call read_var('NameModel',NameModel)
          call set_empirical_model(trim(NameModel))
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
  
  !========================================================================
  subroutine user_face_bcs(VarsGhostFace_V)

    use ModSize,       ONLY: East_,West_,South_,North_,Bot_,Top_
    use ModMain,       ONLY: time_accurate,x_,y_,z_, UseRotatingFrame
    use ModVarIndexes, ONLY: nVar,Ew_,rho_,Ux_,Uy_,Uz_,Bx_,By_,Bz_,P_

    use ModGeometry,   ONLY: R_BLK
    use ModAdvance,    ONLY: State_VGB
    use ModPhysics,    ONLY: inv_gm1, OmegaBody
    use ModNumConst,   ONLY: cTolerance,cTiny

    use ModBlockData, ONLY: use_block_data, put_block_data, get_block_data

    use ModFaceBc, ONLY: FaceCoords_D, VarsTrueFace_V, TimeBc, &
         iFace, jFace, kFace, iSide, iBlockBc

    implicit none

    real, intent(out):: VarsGhostFace_V(nVar)

    integer:: iCell,jCell,kCell

    real:: DensCell,PresCell,GammaCell, B1dotR  
    real, dimension(3):: RFace_D,B1_D,U_D,B1t_D,B1n_D

    ! Variables related to the TD99 flux rope::
    real:: BFRdotR,RhoFRope=cZero
    real, dimension(3):: BFRope_D,BFRn_D,BFRt_D,UVorT_D

    !--------------------------------------------------------------------------

    RFace_D  = FaceCoords_D/sqrt(sum(FaceCoords_D**2))

    U_D (x_:z_)  = VarsTrueFace_V(Ux_:Uz_)
    B1_D(x_:z_)  = VarsTrueFace_V(Bx_:Bz_)
    B1dotR       = dot_product(RFace_D,B1_D)
    B1n_D(x_:z_) = B1dotR*RFace_D(x_:z_)
    B1t_D        = B1_D-B1n_D

    !\
    ! Update BCs for velocity and induction field::
    !/
    VarsGhostFace_V(Ux_:Uz_) = -U_D(x_:z_)
    VarsGhostFace_V(Bx_:Bz_) = B1t_D(x_:z_)!-B1n_D(x_:z_)
    !\
    ! Compute the magnetic field of TD99 flux rope at RFace_D::
    !/
    if (DoTD99FluxRope.or.DoBqField) then
       call get_transformed_TD99fluxrope(RFace_D,BFRope_D,&
            UVorT_D,RhoFRope,TimeBc)
       if(.not.DoBqField)UVorT_D=0.0
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
            BFRn_D(x_:z_)
    else
       RhoFRope=cZero
    end if
    !\
    ! Update BCs for the mass density, EnergyRL,
    ! and pressure::
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
       write(*,*)'ERROR: iSide = ',iSide
       call stop_mpi('incorrect iSide value in user_face_bcs')
    end select

    call get_plasma_parameters_cell(iCell,jCell,kCell,iBlockBc,&
         DensCell,PresCell,GammaCell)
    VarsGhostFace_V(rho_     ) = max(-VarsTrueFace_V(rho_     )+ &
         2*(DensCell+RhoFRope),&!+RhoFRope)
         VarsTrueFace_V(rho_))
    VarsGhostFace_V(P_       ) =max(VarsGhostFace_V(rho_     )*&
         PresCell/(DensCell+RhoFRope),&
                                !max(-VarsTrueFace_V(P_       )+ &
                                !2*PresCell,&
         VarsTrueFace_V(P_  ))
    VarsGhostFace_V(Ew_) = &!max(-VarsTrueFace_V(Ew_)+ &  
         VarsGhostFace_V(rho_     )/(DensCell+RhoFRope)*&  !2* ???
         PresCell*(1.0/(GammaCell-cOne)-inv_gm1)

    !\
    ! Apply corotation
    !/
    if (.not.UseRotatingFrame) then

       VarsGhostFace_V(Ux_) = VarsGhostFace_V(Ux_) -&
            2*OmegaBody*FaceCoords_D(y_)
       VarsGhostFace_V(Uy_) = VarsGhostFace_V(Uy_) +&
            2*OmegaBody*FaceCoords_D(x_)
    end if
  end subroutine user_face_bcs

  !============================================================================

  subroutine get_plasma_parameters_cell(iCell,jCell,kCell,iBlock,&
       DensCell,PresCell,GammaCell)

    ! This subroutine computes the cell values for density and pressure 
    ! assuming an isothermal atmosphere
    
    use ModVarIndexes
    use ModGeometry,   ONLY: x_BLK,y_BLK,z_BLK,R_BLK
    use ModNumConst
    use ModPhysics,    ONLY: g,inv_g,GBody,BodyTdim_I
    use ModExpansionFactors,  ONLY: UMin,T0
    implicit none

    integer, intent(in)  :: iCell,jCell,kCell,iBlock
    real, intent(out)    :: DensCell,PresCell,GammaCell
    real :: UFinal       !The solar wind speed at the far end of the Parker spiral,
                         !which originates from the given cell
    real :: URatio       !The coronal based values for temperature density 
                         !are scaled as functions of UFinal/UMin ratio
    !--------------------------------------------------------------------------

    call get_gamma_emp(x_BLK(iCell,jCell,kCell,iBlock),&
         y_BLK(iCell,jCell,kCell,iBlock),&
         z_BLK(iCell,jCell,kCell,iBlock),&
         GammaCell)
    call get_bernoulli_integral(x_BLK(iCell,jCell,kCell,iBlock)/&
         R_BLK(iCell,jCell,kCell,iBlock),&
         y_BLK(iCell,jCell,kCell,iBlock)/R_BLK(iCell,jCell,kCell,iBlock),&
         z_BLK(iCell,jCell,kCell,iBlock)/R_BLK(iCell,jCell,kCell,iBlock),UFinal)
    URatio=UFinal/UMin
    DensCell  = ((cOne/URatio)**2)*&               !This is the density variation
         exp(-GBody*g*&
         (min(URatio,2.0)*BodyTdim_I(1)/T0)*&!This is the temperature variation
         (cOne/max(R_BLK(iCell,jCell,kCell,iBlock),0.90)&
         -cOne))

    PresCell  = inv_g*DensCell*&
         T0/(min(URatio,2.0)*BodyTdim_I(1))  !This is the temperature variation
  end subroutine get_plasma_parameters_cell

  !========================================================================

  subroutine user_initial_perturbation
    use ModMain,      ONLY: nI, nJ, nK, nBLK, unusedBLK, gcn, x_, y_, z_
    use ModIO,        ONLY: restart
    use ModVarIndexes
    use ModAdvance,   ONLY: State_VGB 
    use ModNumConst
    use ModPhysics,   ONLY:inv_gm1
    use ModGeometry
    use ModEnergy,    ONLY: calc_energy_cell
    implicit none

    !\
    ! Variables required by this user subroutine::
    !/
    integer:: i,j,k,iBLK,iError
    logical:: oktest,oktest_me
    real:: Dens_BLK,Pres_BLK,Gamma_BLK
    real:: xx,yy,zz,RR,ROne,Rmax,Speed
    real, dimension(3):: R_TD99_D,B_TD99_D,U_TD99  ! To include TD99 flux rope.
    real:: Rho_TD99=cZero                          ! To include TD99 flux rope.
    !
    !---------------------------------------------------------------------------
    !\
    ! Variable meanings:
    !
    !
    !/
    !---------------------------------------------------------------------------
    !
    call set_oktest('user_initial_perturbation',oktest,oktest_me)
    do iBLK=1,nBLK
       if(unusedBLK(iBLK))CYCLE
       if ((.not.restart)) then   
          do k=1,nK;do j=1,nJ; do i=1,nI
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
                  4.0E+00*((ROne-cOne)/(Rmax-cOne))*xx/RR
             State_VGB(rhoUy_   ,i,j,k,iBLK) = Dens_BLK*&
                  4.0E+00*((ROne-cOne)/(Rmax-cOne))*yy/RR
             State_VGB(rhoUz_   ,i,j,k,iBLK) = Dens_BLK*&
                  4.0E+00*((ROne-cOne)/(Rmax-cOne))*zz/RR
             State_VGB(Ew_,i,j,k,iBLK) = Pres_BLK   *&
                  (cOne/(Gamma_BLK-cOne)-inv_gm1) 
          end do;end do; end do
       elseif(UseTD99Perturbation)then
          !----------------------------------------------------------------
          !\
          ! Add Titov & Demoulin (TD99) flux rope here:: 
          !/
          !----------------------------------------------------------------
          do k=1,nK; do j=1,nJ; do i=1,nI
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
             if (.not.UseVariedCurrent) then
                call get_transformed_TD99fluxrope(R_TD99_D,B_TD99_D,&
                     U_TD99,Rho_TD99)
             else
                B_TD99_D=cZero
             end if
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
          end do; end do; end do
       endif
       !----------------------------------------------------------------


       !\
       ! Update the total energy::
       !/
       call calc_energy_cell(iBLK)
    end do
  end subroutine user_initial_perturbation

  !============================================================================
  subroutine user_get_b0(xInput,yInput,zInput,B0_D)
    use ModPhysics,  ONLY: Io2No_V,UnitB_
    implicit none
    real, intent(in):: xInput,yInput,zInput
    real, intent(out), dimension(3):: B0_D
    call get_magnetogram_field(xInput,yInput,zInput,B0_D)
    B0_D = B0_D*Io2No_V(UnitB_)
  end subroutine user_get_b0

  !===========================================================================
  subroutine user_update_states(iStage,iBlock)
    use ModVarIndexes
    use ModSize
    use ModAdvance, ONLY: State_VGB, B0xCell_BLK, B0yCell_BLK, B0zCell_BLK
    use ModMain,    ONLY: nStage
    use ModPhysics, ONLY: inv_gm1
    use ModGeometry,ONLY: R_BLK
    use ModEnergy,  ONLY: calc_energy_cell
    use ModExpansionFactors, ONLY: gammaSS
    implicit none
    integer,intent(in):: iStage,iBlock
    integer:: i,j,k
    real:: DensCell,PresCell,GammaCell,Beta
    call update_states_MHD(iStage,iBlock)
    !\
    ! Begin update of pressure and relaxation energy::
    !/
    !  if (iStage/=nStage) return
    do k=1,nK; do j=1,nJ; do i=1,nI
       call get_plasma_parameters_cell(i,j,k,iBlock,&
            DensCell,PresCell,GammaCell)
       if(R_BLK(i,j,k,iBlock)>2.5)&
            GammaCell=GammaCell-(GammaCell-gammaSS)*max(0.0, &
            -1.0 + 2*State_VGB(P_,i,j,k,iBlock)/&
            (State_VGB(P_   ,i,j,k,iBlock)+(&
            (State_VGB(Bx_   ,i,j,k,iBlock)+B0xCell_BLK(i,j,k,iBlock))**2+&
            (State_VGB(By_   ,i,j,k,iBlock)+B0yCell_BLK(i,j,k,iBlock))**2+&
            (State_VGB(Bz_   ,i,j,k,iBlock)+B0zCell_BLK(i,j,k,iBlock))**2)&
            *cQuarter*(R_BLK(i,j,k,iBlock)/2.5)**1.50))
       State_VGB(P_   ,i,j,k,iBlock)=(GammaCell-cOne)*      &
            (inv_gm1*State_VGB(P_,i,j,k,iBlock) + State_VGB(Ew_,i,j,k,iBlock))
       State_VGB(Ew_,i,j,k,iBlock)= State_VGB(P_,i,j,k,iBlock) &
            *(cOne/(GammaCell-cOne)-inv_gm1)
    end do; end do; end do
    call calc_energy_cell(iBlock)
    !\
    ! End update of pressure and relaxation energy::
    !/
  end subroutine user_update_states

  !========================================================================
  subroutine user_get_log_var(VarValue,TypeVar,Radius)

    use ModProcMH,     ONLY: nProc
    use ModIO,         ONLY: dn_output,logfile_,write_myname
    use ModMain,       ONLY: unusedBLK,nBLK,iteration_number,   &
         x_,y_,z_
    use ModVarIndexes, ONLY: Ew_,Bx_,By_,Bz_,rho_,rhoUx_,rhoUy_,rhoUz_,P_ 
    use ModGeometry,   ONLY: R_BLK
    use ModAdvance,    ONLY: State_VGB,tmp1_BLK,B0xCell_BLK,    &
         B0yCell_BLK,B0zCell_BLK
    use ModPhysics,    ONLY: inv_gm1,&
         No2Si_V,UnitEnergydens_,UnitX_,UnitU_,UnitRho_
    use ModNumConst,   ONLY: cOne,cHalf,cE1,cE3,cE6

    real, intent(out):: VarValue
    character (LEN=10), intent(in):: TypeVar 
    real, intent(in), optional :: Radius
    !
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
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(:,:,:,iBLK) = &
               (B0xcell_BLK(:,:,:,iBLK)+State_VGB(Bx_,:,:,:,iBLK))**2+&
               (B0ycell_BLK(:,:,:,iBLK)+State_VGB(By_,:,:,:,iBLK))**2+&
               (B0zcell_BLK(:,:,:,iBLK)+State_VGB(Bz_,:,:,:,iBLK))**2
       end do
       VarValue = unit_energy*cHalf*integrate_BLK(1,tmp1_BLK)
    case('ek_t','Ek_t','ek_r','Ek_r')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(:,:,:,iBLK) = &
               (State_VGB(rhoUx_,:,:,:,iBLK)**2 +&
               State_VGB(rhoUy_,:,:,:,iBLK)**2 +&
               State_VGB(rhoUz_,:,:,:,iBLK)**2)/&
               State_VGB(rho_  ,:,:,:,iBLK)             
       end do
       VarValue = unit_energy*cHalf*integrate_BLK(1,tmp1_BLK)
    case('et_t','Et_t','et_r','Et_r')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
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
          if (unusedBLK(iBLK)) CYCLE
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
       write(*,*) 'Warning in set_user_logvar: unknown logvarname = ',TypeVar
    end select
  end subroutine user_get_log_var
  !----------------------------------------------------------
  subroutine user_specify_initial_refinement(iBLK,refineBlock,lev,DxBlock, &
       xCenter,yCenter,zCenter,rCenter,                        &
       minx,miny,minz,minR,maxx,maxy,maxz,maxR,found)
    use ModMain,ONLY:time_loop,nI,nJ,nK
    use ModAMR,ONLY:InitialRefineType
    use ModNumConst
    use ModAdvance,ONLY:&
         State_VGB,Bx_,By_,Bz_,B0xCell_BLK,B0yCell_BLK,B0zCell_BLK
    use ModGeometry
    use ModPhysics,ONLY:rBody
    implicit none
    logical,intent(out) :: refineBlock, found
    integer, intent(in) :: lev
    real, intent(in)    :: DxBlock
    real, intent(in)    :: xCenter,yCenter,zCenter,rCenter
    real, intent(in)    :: minx,miny,minz,minR
    real, intent(in)    :: maxx,maxy,maxz,maxR
    integer, intent(in) :: iBLK

    character (len=*), parameter :: Name='user_specify_initial_refinement'
    real::BDotRMin,BDotRMax,critx
    integer::i,j,k
    !-------------------------------------------------------------------
    select case (InitialRefineType)
    case ('helio_init')
       if(.not.time_loop)then
          !refine to have resolution not worse 4.0 and
          !refine the body intersecting blocks
          ! Refine all blocks time through (starting with 1 block)
          if (lev <= 4) then
             refineBlock = .true.
          else
             critx=(XyzMax_D(1)-XyzMin_D(1))/(2.0**real(lev-2))
             if ( rCenter < 1.10*rBody + critx ) then
                refineBlock = .true.
             else
                refineBlock = .false.
             end if
          endif
       elseif(dx_BLK(iBLK)<0.20.or.far_field_BCs_BLK(iBLK))then
          refineBlock=.false. !Do not refine body or outer boundary
       else
          !refine heliosheath
          BDotRMin=cZero
          do k=0,nK+1;do j=1,nJ
             BDotRMin=min( BDotRMin,minval(&
                  (B0xCell_BLK(1:nI,j,k,iBLK)+&
                  State_VGB(Bx_,1:nI,j,k,iBLK))*&
                  x_BLK(1:nI,j,k,iBLK)+&
                  (B0yCell_BLK(1:nI,j,k,iBLK)+&
                  State_VGB(By_,1:nI,j,k,iBLK))*&
                  y_BLK(1:nI,j,k,iBLK)+&
                  (B0zCell_BLK(1:nI,j,k,iBLK)+&
                  State_VGB(Bz_,1:nI,j,k,iBLK))*&
                  z_BLK(1:nI,j,k,iBLK)))
          end do;end do
          BDotRMax=cZero
          do k=0,nK+1;do j=1,nJ
             BDotRMax=max( BDotRMax,maxval(&
                  (B0xCell_BLK(1:nI,j,k,iBLK)+&
                  State_VGB(Bx_,1:nI,j,k,iBLK))*&
                  x_BLK(1:nI,j,k,iBLK)+&
                  (B0yCell_BLK(1:nI,j,k,iBLK)+&
                  State_VGB(By_,1:nI,j,k,iBLK))*&
                  y_BLK(1:nI,j,k,iBLK)+&
                  (B0zCell_BLK(1:nI,j,k,iBLK)+&
                  State_VGB(Bz_,1:nI,j,k,iBLK))*&
                  z_BLK(1:nI,j,k,iBLK)))
          end do;end do
          refineBlock =BDotRMin<-cTiny.and.&
               BDotRMax>cTiny
       end if
       found=.true.
    end select
  end subroutine user_specify_initial_refinement

end module ModUser

