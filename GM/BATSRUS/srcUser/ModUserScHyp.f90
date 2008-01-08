!^CFG COPYRIGHT UM
!========================================================================
module ModUser
  use ModNumConst, ONLY: cHalf,cTwo,cThree,&
       cFour,cE1,cHundred,cHundredth,cZero,&
       cOne
  use ModMain,     ONLY: UseUserB0
  use ModSize,     ONLY: nI,nJ,nK,gcn,nBLK
  use ModMagnetogram
  use ModExpansionFactors
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_set_ics,                    &
       IMPLEMENTED3 => user_face_bcs,                   &
       IMPLEMENTED5 => user_get_b0,                     &
       IMPLEMENTED6 => user_update_states,              &
       IMPLEMENTED7 => user_specify_initial_refinement

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.1
  character (len=*), parameter :: &
       NameUserModule = 'EMPIRICAL SC - Cohen, Sokolov. van der Holst'

contains

  !==========================================================================
  subroutine user_read_inputs
    use ModMain
    use ModProcMH,    ONLY: iProc
    use ModReadParam
    use ModIO,        ONLY: write_prefix, write_myname, iUnitOut

    integer:: i
    character (len=100) :: NameCommand
    character (len=lStringLine)   :: NameModel
    !------------------------------------------------------------------------

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

  !==========================================================================
  subroutine user_set_ics
    use ModMain,      ONLY: nI, nJ, nK, nBLK, unusedBLK
    use ModVarIndexes
    use ModAdvance,   ONLY: State_VGB 
    use ModNumConst
    use ModPhysics,   ONLY: inv_gm1
    use ModGeometry
    use ModEnergy,    ONLY: calc_energy_cell
    implicit none

    integer :: i,j,k,iBLK
    logical :: oktest,oktest_me
    real    :: Denscell,Prescell,Gammacell
    real    :: xx,yy,zz,RR,ROne,Rmax
    !------------------------------------------------------------------------

    call set_oktest('user_set_ics',oktest,oktest_me)
    do iBLK=1,nBLK
       if(unusedBLK(iBLK))CYCLE
       do k=1,nK;do j=1,nJ; do i=1,nI
          xx = x_BLK(i,j,k,iBLK)
          yy = y_BLK(i,j,k,iBLK)
          zz = z_BLK(i,j,k,iBLK)
          RR = sqrt(xx**2+yy**2+zz**2+cTolerance**2)
          ROne  = max(cOne,RR)
          Rmax  = max(2.1E+01,sqrt(xx**2+yy**2+zz**2))
          State_VGB(Bx_      ,i,j,k,iBLK) = cZero
          State_VGB(By_      ,i,j,k,iBLK) = cZero
          State_VGB(Bz_      ,i,j,k,iBLK) = cZero
          State_VGB(Hyp_     ,i,j,k,iBLK) = cZero
          call get_plasma_parameters_cell(xx,yy,zz,RR,&
               Denscell,Prescell,Gammacell)
          State_VGB(rho_     ,i,j,k,iBLK) = Denscell
          State_VGB(P_       ,i,j,k,iBLK) = Prescell
          State_VGB(rhoUx_   ,i,j,k,iBLK) = Denscell*&
               4.0E+00*((ROne-cOne)/(Rmax-cOne))*xx/RR
          State_VGB(rhoUy_   ,i,j,k,iBLK) = Denscell*&
               4.0E+00*((ROne-cOne)/(Rmax-cOne))*yy/RR
          State_VGB(rhoUz_   ,i,j,k,iBLK) = Denscell*&
               4.0E+00*((ROne-cOne)/(Rmax-cOne))*zz/RR
          State_VGB(Ew_      ,i,j,k,iBLK) = Prescell*&
               (cOne/(Gammacell-cOne)-inv_gm1) 
       end do;end do; end do

       !\
       ! Update the total energy::
       !/
       call calc_energy_cell(iBLK)
    end do
  end subroutine user_set_ics
  
  !==========================================================================
  subroutine user_face_bcs(VarsGhostFace_V)

    use ModSize,       ONLY: East_,West_,South_,North_,Bot_,Top_
    use ModMain,       ONLY: time_accurate,x_,y_,z_, UseRotatingFrame, &
                             UseHyperbolicDivb
    use ModVarIndexes, ONLY: nVar,Ew_,rho_,Ux_,Uy_,Uz_,Bx_,By_,Bz_,P_

    use ModGeometry,   ONLY: R_BLK
    use ModAdvance,    ONLY: State_VGB
    use ModPhysics,    ONLY: inv_gm1, OmegaBody
    use ModNumConst,   ONLY: cTolerance,cTiny

    use ModBlockData, ONLY: use_block_data, put_block_data, get_block_data

    use ModFaceBc, ONLY: FaceCoords_D, VarsTrueFace_V, TimeBc, &
         iBlockBc, B0Face_D

    implicit none

    real, intent(out):: VarsGhostFace_V(nVar)

    integer :: iHyp
    real    :: Denscell, Prescell, Gammacell
    real    :: RR, FullBnorm
    real, dimension(3) :: RFace_D, B1_D, U_D, B1t_D, B1n_D, &
                          FullB_D, FullBunit_D
    !------------------------------------------------------------------------

    RR = sqrt(sum(FaceCoords_D**2))
    RFace_D = FaceCoords_D/RR

    U_D (x_:z_)  = VarsTrueFace_V(Ux_:Uz_)
    B1_D(x_:z_)  = VarsTrueFace_V(Bx_:Bz_)
    B1n_D(x_:z_) = dot_product(RFace_D,B1_D)*RFace_D(x_:z_)
    B1t_D        = B1_D-B1n_D

    !\
    ! Update BCs for induction field::
    !/
    VarsGhostFace_V(Bx_:Bz_) = B1t_D(x_:z_)

    !\
    ! Update BCs for velocity
    !/
    FullB_D = VarsGhostFace_V(Bx_:Bz_) + B0Face_D
    FullBnorm = sqrt(sum(FullB_D**2))
    if (FullBnorm<cTolerance) then
       ! HD: extrapolate in radial direction
       VarsGhostFace_V(Ux_:Uz_) = dot_product(U_D,RFace_D)*RFace_D
    else
       ! MHD: extrapolate in field line direction
       FullBunit_D = FullB_D/FullBnorm
       VarsGhostFace_V(Ux_:Uz_) = dot_product(U_D,FullBunit_D)*FullBunit_D
    end if

    !\
    ! Update BCs for the mass density, EnergyRL, and pressure::
    !/
    call get_plasma_parameters_cell(FaceCoords_D(x_),FaceCoords_D(y_), &
                                    FaceCoords_D(z_),RR, &
                                    Denscell,Prescell,Gammacell)
    VarsGhostFace_V(rho_) = max(-VarsTrueFace_V(rho_) + 2*Denscell, &
                                VarsTrueFace_V(rho_))
    VarsGhostFace_V(P_)   = max(VarsGhostFace_V(rho_)*Prescell/Denscell, &
                                VarsTrueFace_V(P_))
    VarsGhostFace_V(Ew_)  = &!max(-VarsTrueFace_V(Ew_)+ &  
         VarsGhostFace_V(rho_)*Prescell/Denscell &
         *(1.0/(Gammacell-cOne)-inv_gm1)

    !\
    ! Apply corotation
    !/
    if (.not.UseRotatingFrame) then

       VarsGhostFace_V(Ux_) = VarsGhostFace_V(Ux_) -&
            2*OmegaBody*FaceCoords_D(y_)
       VarsGhostFace_V(Uy_) = VarsGhostFace_V(Uy_) +&
            2*OmegaBody*FaceCoords_D(x_)
    end if

    if(UseHyperbolicDivb)then
       iHyp = Bz_ + 1
       VarsGhostFace_V(iHyp) = cZero
    end if

  end subroutine user_face_bcs

  !==========================================================================
  subroutine get_plasma_parameters_cell(xx,yy,zz,RR, &
                                        Denscell,Prescell,Gammacell)

    ! This subroutine computes the values for density and pressure 
    ! assuming an isothermal atmosphere
    
    use ModVarIndexes
    use ModNumConst
    use ModPhysics,           ONLY: g,inv_g,GBody,BodyTdim_I
    use ModExpansionFactors,  ONLY: UMin,T0
    implicit none

    real, intent(in)     :: xx, yy, zz, RR
    real, intent(out)    :: Denscell,Prescell,Gammacell

    real :: UFinal       !The solar wind speed at the far end of the Parker spiral,
                         !which originates from the given cell
    real :: URatio       !The coronal based values for temperature density 
                         !are scaled as functions of UFinal/UMin ratio
    !------------------------------------------------------------------------

    call get_gamma_emp(xx,yy,zz,Gammacell)
    call get_bernoulli_integral(xx/RR,yy/RR,zz/RR,UFinal)
    URatio=UFinal/UMin
    Denscell  = ((cOne/URatio)**2)*&         !This is the density variation
         exp(-GBody*g*&
         (min(URatio,2.0)*BodyTdim_I(1)/T0)*&!This is the temperature variation
         (cOne/max(RR,0.90)-cOne))

    Prescell  = inv_g*Denscell*&
         T0/(min(URatio,2.0)*BodyTdim_I(1))  !This is the temperature variation
  end subroutine get_plasma_parameters_cell

  !==========================================================================
  subroutine user_get_b0(xInput,yInput,zInput,B0_D)
    use ModPhysics,  ONLY: Io2No_V,UnitB_
    implicit none

    real, intent(in):: xInput,yInput,zInput
    real, intent(out), dimension(3):: B0_D
    !------------------------------------------------------------------------

    call get_magnetogram_field(xInput,yInput,zInput,B0_D)
    B0_D = B0_D*Io2No_V(UnitB_)
  end subroutine user_get_b0

  !==========================================================================
  subroutine user_update_states(iStage,iBlock)
    use ModVarIndexes
    use ModSize
    use ModAdvance, ONLY: State_VGB, B0xCell_BLK, B0yCell_BLK, B0zCell_BLK
    use ModMain,    ONLY: nStage
    use ModPhysics, ONLY: inv_gm1
    use ModGeometry,ONLY: x_BLK, y_BLK, z_BLK, R_BLK
    use ModEnergy,  ONLY: calc_energy_cell
    implicit none

    integer, intent(in) :: iStage, iBlock

    integer :: i, j, k
    real    :: Denscell,Prescell,Gammacell,Beta
    !------------------------------------------------------------------------

    call update_states_MHD(iStage,iBlock)
    !\
    ! Begin update of pressure and relaxation energy::
    !/
    !  if (iStage/=nStage) return
    do k=1,nK; do j=1,nJ; do i=1,nI
       call get_plasma_parameters_cell( &
               x_BLK(i,j,k,iBlock),y_BLK(i,j,k,iBlock), &
               z_BLK(i,j,k,iBlock),R_BLK(i,j,k,iBlock), &
               Denscell,Prescell,Gammacell)
       if(R_BLK(i,j,k,iBlock)>2.5)&
            Gammacell=Gammacell-(Gammacell-1.1)*max(0.0, &
            -1.0 + 2*State_VGB(P_,i,j,k,iBlock)/&
            (State_VGB(P_    ,i,j,k,iBlock)+(&
            (State_VGB(Bx_   ,i,j,k,iBlock)+B0xCell_BLK(i,j,k,iBlock))**2+&
            (State_VGB(By_   ,i,j,k,iBlock)+B0yCell_BLK(i,j,k,iBlock))**2+&
            (State_VGB(Bz_   ,i,j,k,iBlock)+B0zCell_BLK(i,j,k,iBlock))**2)&
            *cQuarter*(R_BLK(i,j,k,iBlock)/2.5)**1.50))
       State_VGB(P_   ,i,j,k,iBlock)=(Gammacell-cOne)*      &
            (inv_gm1*State_VGB(P_,i,j,k,iBlock) + State_VGB(Ew_,i,j,k,iBlock))
       State_VGB(Ew_,i,j,k,iBlock)= State_VGB(P_,i,j,k,iBlock) &
            *(cOne/(Gammacell-cOne)-inv_gm1)
    end do; end do; end do
    call calc_energy_cell(iBlock)
    !\
    ! End update of pressure and relaxation energy::
    !/
  end subroutine user_update_states

  !==========================================================================
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
       elseif(lev<=6)then
          !refine heliosheath up to 6 levels
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
       else
          refineBlock=.false.
       end if
       found=.true.
    end select
  end subroutine user_specify_initial_refinement

end module ModUser

