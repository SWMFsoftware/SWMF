
!========================================================================
Module ModUser
  use ModNumConst, ONLY: cHalf,cTwo,cThree,&
       cFour,cE1,cHundred,cHundredth,cZero,&
       cOne
  use ModMain,     ONLY: UseUserB0,UseUserHeating
  use ModSize,     ONLY: nI,nJ,nK,gcn,nBLK
  !  use ModUserTD99  ! To include TD99 flux rope.
  use ModUserEmpty, ONLY:               &
!!!       user_read_inputs,                &
!!!       user_init_session,               &
       user_set_ics,                    &
!!!       user_initial_perturbation,       &
       user_set_boundary_cells,         &
!!!       user_face_bcs,                   &
       user_set_outerbcs,               &
!!!       user_specify_initial_refinement, &
!!!       user_amr_criteria,               &
       user_write_progress,             &
!!!       user_get_log_var,                &
       user_set_plot_var,               &
       user_calc_sources,               &
!!!       user_get_b0,                     &
!!!       user_update_states,              &
       user_io_units

  include 'user_module.h' !list of public methods
 
  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'HELIOSPHERE, Manchester, Roussev'

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
       case("#USEUSERHEATING")
          call read_var('UseUserHeating'          ,UseUserHeating)
       case("#PFSSM")
          call read_var('UseUserB0'  ,UseUserB0)
          if (UseUserB0)then
             call SC_read_magnetogram_file
             call SC_set_expansion_factors
             call read_var('dt_UpdateB0',dt_UpdateB0)
             DoUpdateB0 = dt_updateb0 > 0.0
          endif
       case("#AWHEAT")
          call read_var('Bnot        ',Bnot)
          call read_var('Tnot        ',Tnot)
          call read_var('DegFrm1     ',DegFrm1)
          call read_var('DegF_Ratio  ',DegF_Ratio)
          call read_var('Dens_Ratio  ',Dens_Ratio)
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
  subroutine user_face_bcs(iFace,jFace,kFace,iBlock,iSide,iBoundary, &
       iter,time_now,FaceCoords_D,VarsTrueFace_V,VarsGhostFace_V,    &
       B0Face_D,UseIonosphereHere,UseRotatingBcHere)
    use ModSize,       ONLY: nDim,East_,West_,South_,North_,Bot_,    &
         Top_
    use ModMain,       ONLY: time_accurate,UseUserHeating,x_,y_,z_,  &
         UseRotatingFrame
    use ModVarIndexes, ONLY: rho_,Ux_,Uy_,Uz_,Bx_,By_,Bz_,P_,Ew_ 

    use ModGeometry,   ONLY: R_BLK
    use ModAdvance,    ONLY: nFaceValueVars,State_VGB
    use ModPhysics,    ONLY: CosThetaTilt,SinThetaTilt,g,inv_g,      &
         inv_gm1,OmegaBody,unitUSER_B
    use ModNumConst,   ONLY: cZero,cHalf,cOne,cTwo,cTolerance,cTiny

    use ModBlockData, ONLY: use_block_data, put_block_data, get_block_data
    implicit none
    !\
    ! Variables required by this user subroutine
    !/
    integer, intent(in):: iFace,jFace,kFace,iBlock,iSide,&
         iBoundary,iter
    real, intent(in):: time_now
    real, dimension(nDim), intent(in):: FaceCoords_D,    &
         B0Face_D
    real, dimension(nFaceValueVars), intent(in)::        &
         VarsTrueFace_V
    logical, intent(in):: UseIonosphereHere,             &
         UseRotatingBcHere
    real, dimension(nFaceValueVars), intent(out)::       &
         VarsGhostFace_V
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

    !\
    ! Calculation of boundary conditions should start here::
    !/
    !--------------------------------------------------------------------------
    !
    XFace = FaceCoords_D(1)
    YFace = FaceCoords_D(2)
    ZFace = FaceCoords_D(3)
    RFace = sqrt(XFace**2+YFace**2+ZFace**2)

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
       call stop_mpi('user_face_bcs')
    end select
    call get_plasma_parameters_cell(iCell,jCell,kCell,iBlock,&
         DensCell,PresCell,GammaCell)
    VarsGhostFace_V(rho_     ) = max(-VarsTrueFace_V(rho_     )+ &
         cTwo*(DensCell),&!+RhoFRope)
         VarsTrueFace_V(rho_))
    VarsGhostFace_V(P_       ) = max(-VarsTrueFace_V(P_       )+ &
         cTwo*PresCell,VarsTrueFace_V(P_  ))
    VarsGhostFace_V(Ew_) = max(-VarsTrueFace_V(Ew_)+ &  
         cTwo*PresCell*(cOne/(GammaCell-cOne)-inv_gm1),&
         VarsTrueFace_V(Ew_))

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
    use ModAdvance,    ONLY: B0xCell_BLK,B0yCell_BLK,B0zCell_BLK,  &
         State_VGB
    use ModPhysics,    ONLY: g,inv_g,unitUSER_B
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
    MaxB0_1 = Bnot/unitUSER_B
    MaxB0_2 = 2.00E+01/unitUSER_B
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
            (XCell*B0xCell_BLK(iCell,jCell,kCell,iBlock) +&
            YCell*B0yCell_BLK(iCell,jCell,kCell,iBlock) +&
            ZCell*B0zCell_BLK(iCell,jCell,kCell,iBlock))/&
            RCell)
       B2Cell = sqrt(&
            B0xCell_BLK(iCell,jCell,kCell,iBlock)**2+&
            B0yCell_BLK(iCell,jCell,kCell,iBlock)**2+&
            B0zCell_BLK(iCell,jCell,kCell,iBlock)**2)
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
         max(cTiny, &
         ((B0xCell_BLK(iCell,jCell,kCell,iBlock)+    &
         State_VGB(Bx_,iCell,jCell,kCell,iBlock))**2+&
         (B0yCell_BLK(iCell,jCell,kCell,iBlock)+    &
         State_VGB(By_,iCell,jCell,kCell,iBlock))**2+&
         (B0zCell_BLK(iCell,jCell,kCell,iBlock)+    &
         State_VGB(Bz_,iCell,jCell,kCell,iBlock))**2))
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
  !  SUBROUTINE user_initial_perturbation
  !========================================================================
  !\
  ! This subroutine allows the user to add a perturbation to a solutions
  ! read in from a restart file.  The idea is to allow the user to "start"
  ! some process on top of the already converged solution. The routine loops
  ! over all blocks, but only needs to load the perturbation where appropriate.
  ! As with all user subroutines, the variables declared in ModUser are 
  ! available here.  Again, as with other user subroutines DO NOT MODIFY ANY 
  ! GLOBAL VARIABLE DEFINED IN THE MODULES INCLUDED IN THIS SUBROUTINE UNLESS 
  ! SPECIFIED!!
  !
  ! The user should load the global variables:
  !   rho_BLK,rhoUx_BLk,rhoUy_BLK,rhoUz_BLK,p_BLK,E_BLK
  !   Bx_BLK, By_BLK, Bz_BLK
  !
  ! Note that in most cases E (energy) and P (pressure) should both be loaded.
  !/
  subroutine user_initial_perturbation
    use ModMain,      ONLY: nI,nJ,nK,nBLK,                           &
         unusedBLK,UseUserHeating,UseUserB0,gcn,x_,y_,z_
    use ModIO,        ONLY: restart
    use ModVarIndexes,ONLY: rho_,rhoUx_,rhoUy_,rhoUz_,Bx_,By_,Bz_,P_,Ew_
    use ModAdvance,   ONLY: State_VGB,B0xCell_BLK,B0yCell_BLK,       &
         B0zCell_BLK,tmp1_BLK,tmp2_BLK
    use ModProcMH,    ONLY: iProc,nProc,iComm
    use ModNumConst,  ONLY: cZero,cQuarter,cHalf,cOne,cTwo,cE1,cE9,  &
         cTolerance,cThree
    use ModConst,     ONLY: Rsun,Msun,cGravitation
    use ModGeometry,  ONLY: x_BLK,y_BLK,z_BLK,R_BLK,cV_BLK,x2,y2,z2
    use ModPhysics,   ONLY: Gbody,g,inv_g,gm1,inv_gm1,ModulationP,   &
         ModulationRho,UseFluxRope,rot_period_dim,OmegaBody,Rbody,   &
         unitSI_U,unitSI_rho,unitSI_x,unitUSER_energydens,           &
         unitUSER_t,unitUSER_B,Body_rho_dim
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
    !---------------------------------------------------------------------------
    !\
    ! Variable meanings:
    !
    !
    !/
    !---------------------------------------------------------------------------
    !
    call set_oktest('user_initial_perturbation',oktest,oktest_me)
    !\
    ! Initialize some auxilary variables::
    !/
    Rbody  = cOne
    Mrope_GL98 = cZero

    InvH0 = cGravitation*Msun/Rsun/unitSI_U**2
    do iBLK=1,nBLK
       if ((.not.UseUserHeating).and.(.not.restart)) then
          if (unusedBLK(iBLK)) CYCLE   
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
                  4.0E+01*((ROne-cOne)/(Rmax-cOne))*xx/RR
             State_VGB(rhoUy_   ,i,j,k,iBLK) = Dens_BLK*&
                  4.0E+01*((ROne-cOne)/(Rmax-cOne))*yy/RR
             State_VGB(rhoUz_   ,i,j,k,iBLK) = Dens_BLK*&
                  4.0E+01*((ROne-cOne)/(Rmax-cOne))*zz/RR
             State_VGB(Ew_,i,j,k,iBLK) = Pres_BLK   *& 
                  (cOne/(Gamma_BLK-cOne)-inv_gm1)
          end do;end do; end do

       endif

       !\
       ! Update the total energy::
       !/
       call calc_energy(iBLK)
    end do
    !\
    ! Write out some statistics::
    !/
    call post_init_stat

  end subroutine user_initial_perturbation

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

  !============================================================================
  subroutine user_get_b0(xInput,yInput,zInput,B0_D)
    use ModPhysics,  ONLY: unitUSER_B
    implicit none
    real, intent(in):: xInput,yInput,zInput
    real, intent(out), dimension(3):: B0_D
    call SC_get_magnetogram_field(xInput,yInput,zInput,B0_D)
    B0_D = B0_D/unitUSER_B
  end subroutine user_get_b0

  !========================================================================
  subroutine user_specify_initial_refinement(iBLK,RefineBlock,lev, &
       dxBlock,xCenter,yCenter,zCenter,rCenter,minx,miny,minz,minR,&
       maxx,maxy,maxz,maxR,IsFound)
    use ModSize,       ONLY: nI,nJ,nK,gcn
    use ModVarIndexes, ONLY: Bx_,By_,Bz_,P_
    use ModAdvance,    ONLY: B0xCell_BLK,B0yCell_BLK,       &
         B0zCell_BLK,State_VGB
    use ModAMR,        ONLY: InitialRefineType
    use ModMain,       ONLY: x_,y_,z_,iteration_number,     &
         time_loop
    use ModGeometry,   ONLY: XyzMin_D,XyzMax_D,x_BLK,y_BLK, &
         z_BLK,R_BLK,dx_BLK,dy_BLK,dz_BLK
    use ModNumConst,   ONLY: cTiny,cHundredth,cEighth,cHalf,&
         cQuarter,cOne,cTwo,cFour,cE1,cE2,cZero
    use ModPhysics,    ONLY: unitUSER_B

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
                   B0xCell = B0xCell_BLK(i,j,k,iBLK)
                   B0yCell = B0yCell_BLK(i,j,k,iBLK)
                   B0zCell = B0zCell_BLK(i,j,k,iBLK)
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
               IsInRangeAR = (minval(Br_D)<2.50E+00/unitUSER_B)
          !For 01:35UT on Oct 27, 2003: if (lev>12) &
          !For 09:35UT on Oct 28, 2003: if (lev>11) &
          if (lev>11) &
               IsInRangeAR = (minval(Br_D)<2.50E+00/unitUSER_B).and.&
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
                      B0xCell = B0xCell_BLK(i,j,k,iBLK)
                      B0yCell = B0yCell_BLK(i,j,k,iBLK)
                      B0zCell = B0zCell_BLK(i,j,k,iBLK)
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
                      B0xCell  = B0xCell_BLK(i,j,k,iBLK)
                      B0yCell  = B0yCell_BLK(i,j,k,iBLK)
                      B0zCell  = B0zCell_BLK(i,j,k,iBLK)
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
    real::Rs_PFSSM
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
                B0xCell = B0xCell_BLK(i,j,k,iBLK)
                B0yCell = B0yCell_BLK(i,j,k,iBLK)
                B0zCell = B0zCell_BLK(i,j,k,iBLK)
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
    use ModVarIndexes
    use ModSize
    use ModAdvance, ONLY: State_VGB
    use ModMain,    ONLY: nStage
    use ModPhysics, ONLY: inv_gm1
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
       State_VGB(P_ ,i,j,k,iBlock) = (GammaCell-cOne)* &
            (inv_gm1*State_VGB(P_,i,j,k,iBlock) + State_VGB(Ew_,i,j,k,iBlock))
       State_VGB(Ew_,i,j,k,iBlock) = &
            State_VGB(P_,i,j,k,iBlock)*(cOne/(GammaCell - cOne) - inv_gm1)
    end do; end do; end do
    call calc_energy(iBlock)
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
  subroutine user_get_log_var(VarValue,TypeVar)
    use ModProcMH,     ONLY: nProc
    use ModIO,         ONLY: dn_output,logfile_,write_myname
    use ModMain,       ONLY: unusedBLK,nBLK,iteration_number,   &
         x_,y_,z_
    use ModVarIndexes, ONLY: Bx_,By_,Bz_,rho_,rhoUx_,rhoUy_,rhoUz_,P_,Ew_
    use ModGeometry,   ONLY: R_BLK
    use ModAdvance,    ONLY: State_VGB,tmp1_BLK,B0xCell_BLK,    &
         B0yCell_BLK,B0zCell_BLK
    use ModPhysics,    ONLY: inv_gm1,unitSI_energydens,unitSI_x,&
         unitSI_U,unitSI_rho
    use ModNumConst,   ONLY: cOne,cHalf,cE1,cE3,cE6
    real, intent(out):: VarValue
    character (LEN=10), intent(in):: TypeVar 
    !
    integer:: iBLK
    real:: unit_energy,unit_mass
    real, external:: integrate_BLK
    !--------------------------------------------------------------------------
    unit_energy = cE1*cE6*unitSI_energydens*unitSI_x**3
    unit_mass   = cE3*unitSI_rho*unitSI_x**3
    !\
    ! Define log variable to be saved::
    !/
    select case(TypeVar)
    case('em_t','Em_t','em_r','Em_r')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(:,:,:,iBLK) = &
               (B0xcell_BLK(:,:,:,iBLK)+State_VGB(Bx_,:,:,:,iBLK))**2+&
               (B0ycell_BLK(:,:,:,iBLK)+State_VGB(By_,:,:,:,iBLK))**2+&
               (B0zcell_BLK(:,:,:,iBLK)+State_VGB(Bz_,:,:,:,iBLK))**2
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
       write(*,*) 'Warning in set_user_logvar: unknown logvarname = ',TypeVar
    end select
  end subroutine user_get_log_var


end module ModUser

