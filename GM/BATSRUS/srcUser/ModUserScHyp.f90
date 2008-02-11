!^CFG COPYRIGHT UM
!============================================================================

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
       IMPLEMENTED2 => user_initial_perturbation,       &
       IMPLEMENTED3 => user_face_bcs,                   &
       IMPLEMENTED4 => user_get_b0,                     &
       IMPLEMENTED5 => user_update_states,              &
       IMPLEMENTED6 => user_specify_initial_refinement, &
       IMPLEMENTED7 => user_set_outerbcs

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.1
  character (len=*), parameter :: &
       NameUserModule = 'EMPIRICAL SC - Cohen, Sokolov, van der Holst'

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
       case("#PFSSM")
          call read_var('UseUserB0'  ,UseUserB0)
          if (UseUserB0)then
             call read_magnetogram_file
             call read_var('dt_UpdateB0',dt_UpdateB0)
             DoUpdateB0 = dt_updateb0 > 0.0
          endif
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
  subroutine user_set_outerbcs(iBlock, iSide, TypeBc, IsFound)
    use ModAdvance,  ONLY: State_VGB,B0xCell_Blk,B0yCell_Blk,B0zCell_Blk,Hyp_
    use ModMain,     ONLY: UseHyperbolicDivb
    use ModNumConst, ONLY: cTolerance
    use ModGeometry, ONLY: x_Blk, y_Blk, z_Blk, r_Blk
    use ModPhysics,  ONLY: inv_gm1
    use ModSize
    use ModVarIndexes
    implicit none

    integer,          intent(in)  :: iBlock, iSide
    character(len=20),intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    character (len=*), parameter :: Name='user_set_outerbcs'

    integer :: i, j, k, iMin, jMin, kMin, iMax, jMax, kMax
    real, dimension(nDim) :: B1_D,B1n_D,B1t_D, FullB_D
    real, dimension(nDim) :: RhoU_D,RhoUn_D,RhoUt_D,r2RhoUn_D,GRhoU_D
    real, dimension(nDim) :: FaceCoords_D, Normal_D, Unit_D
    real :: RFace, Dens, Pres, Gamma, B1dotR, FullBNorm, RhoUdotR
    !------------------------------------------------------------------------
    IsFound = .true.

    if(iSide/=East_) call stop_mpi('Wrong iSide in user_set_outerBCs')

    iMin=1-gcn;iMax=0
    jMin=1-gcn;jMax=nJ+gcn
    kMin=1-gcn;kMax=nK+gcn

    do k=kMin,kMax; do j=jMin,jMax
       FaceCoords_D(x_)=0.5*sum(x_Blk(0:1,j,k,iBlock))
       FaceCoords_D(y_)=0.5*sum(y_Blk(0:1,j,k,iBlock))
       FaceCoords_D(z_)=0.5*sum(z_Blk(0:1,j,k,iBlock))
       RFace = sqrt(sum(FaceCoords_D**2))

       call get_plasma_parameters_cell(FaceCoords_D(x_),FaceCoords_D(y_), &
            FaceCoords_D(z_),RFace,Dens,Pres,Gamma)

       State_VGB(Rho_,iMax,j,k,iBlock) = max( 2.0*Dens &
            - State_VGB(Rho_,iMax+1,j,k,iBlock), &
            State_VGB(Rho_,iMax+1,j,k,iBlock) )
       do i=iMax-1,iMin,-1
          State_VGB(Rho_,i,j,k,iBlock) = 2.0*State_VGB(Rho_,i+1,j,k,iBlock) &
            - State_VGB(Rho_,i+2,j,k,iBlock)
       end do

       do i=iMin,iMax
          call get_plasma_parameters_cell(x_Blk(i,j,k,iBlock), &
               y_Blk(i,j,k,iBlock),z_Blk(i,j,k,iBlock),r_Blk(i,j,k,iBlock), &
               Dens,Pres,Gamma)
          State_VGB(P_,i,j,k,iBlock) = State_VGB(Rho_,i,j,k,iBlock) &
               *Pres/Dens
          State_VGB(Ew_,i,j,k,iBlock) = State_VGB(P_,i,j,k,iBlock) &
               *(1.0/(Gamma-1.0)-inv_gm1)
       end do

       Normal_D(x_)=x_Blk(1,j,k,iBlock)/r_Blk(1,j,k,iBlock)
       Normal_D(y_)=y_Blk(1,j,k,iBlock)/r_Blk(1,j,k,iBlock)
       Normal_D(z_)=z_Blk(1,j,k,iBlock)/r_Blk(1,j,k,iBlock)

       B1_D   = State_VGB(Bx_:Bz_,1,j,k,iBlock)
       B1dotR = dot_product(Normal_D,B1_D)
       B1n_D  = B1dotR*Normal_D
       B1t_D  = B1_D-B1n_D

       do i=iMin,iMax
          State_VGB(Bx_:Bz_,i,j,k,iBlock) = B1t_D
       end do

       RhoU_D   = State_VGB(RhoUx_:RhoUz_,1,j,k,iBlock)
       RhoUdotR = dot_product(Normal_D,RhoU_D)
       RhoUn_D  = RhoUdotR*Normal_D
       RhoUt_D  = RhoU_D - RhoUn_D
       r2RhoUn_D = r_Blk(1,j,k,iBlock)**2*RhoUn_D

       do i=iMin,iMax
          GRhoU_D = RhoUt_D + r2RhoUn_D/r_Blk(i,j,k,iBlock)**2
          FullB_D(x_) = State_VGB(Bx_,i,j,k,iBlock) &
               + B0xCell_BLK(i,j,k,iBlock)
          FullB_D(y_) = State_VGB(By_,i,j,k,iBlock) &
               + B0yCell_BLK(i,j,k,iBlock)
          FullB_D(z_) = State_VGB(Bz_,i,j,k,iBlock) &
               + B0zCell_BLK(i,j,k,iBlock)
          FullBNorm = sqrt(sum(FullB_D**2))
          if(FullBNorm<cTolerance)then
             ! HD: extrapolate in radial direction
             Unit_D(x_)=x_Blk(i,j,k,iBlock)/r_Blk(i,j,k,iBlock)
             Unit_D(y_)=y_Blk(i,j,k,iBlock)/r_Blk(i,j,k,iBlock)
             Unit_D(z_)=z_Blk(i,j,k,iBlock)/r_Blk(i,j,k,iBlock)
          else
             ! MHD: extrapolate in field line direction
             Unit_D = FullB_D/FullBNorm
          end if
          State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = &
               dot_product(GRhoU_D,Unit_D)*Unit_D
       end do

    end do; end do

    if(UseHyperbolicDivb)then
       State_VGB(Hyp_,iMin:iMax,jMin:jMax,kMin:kMax,iBlock) = 0.0
    end if

  end subroutine user_set_outerbcs

  !==========================================================================
  subroutine user_face_bcs(VarsGhostFace_V)

    use ModSize,       ONLY: East_,West_,South_,North_,Bot_,Top_
    use ModMain,       ONLY: time_accurate,x_,y_,z_, UseRotatingFrame, &
                             UseHyperbolicDivb
    use ModVarIndexes, ONLY: nVar,Ew_,rho_,Ux_,Uy_,Uz_,Bx_,By_,Bz_,P_

    use ModGeometry,   ONLY: R_BLK
    use ModAdvance,    ONLY: State_VGB, Hyp_
    use ModPhysics,    ONLY: inv_gm1, OmegaBody
    use ModNumConst,   ONLY: cTolerance,cTiny

    use ModBlockData, ONLY: use_block_data, put_block_data, get_block_data

    use ModFaceBc, ONLY: FaceCoords_D, VarsTrueFace_V, TimeBc, &
         iBlockBc, iSide !!!, B0Face_D

    implicit none

    real, intent(out):: VarsGhostFace_V(nVar)

    real :: DensCell, PresCell, GammaCell
    real :: RR, B1dotR  
    real, dimension(3):: RFace_D,B1_D,U_D,B1t_D,B1n_D

!!!    real    :: RR, B1dotR, FullBnorm, UdotRFace
!!!    real, dimension(3) :: RFace_D, B1_D, U_D, B1t_D, B1n_D, &
!!!                         FullB_D, FullBunit_D, UFace_D
    !------------------------------------------------------------------------

    RR = sqrt(sum(FaceCoords_D**2))
    RFace_D = FaceCoords_D/RR

    U_D (x_:z_)  = VarsTrueFace_V(Ux_:Uz_)
    B1_D(x_:z_)  = VarsTrueFace_V(Bx_:Bz_)
    B1dotR       = dot_product(RFace_D,B1_D)
    B1n_D(x_:z_) = B1dotR*RFace_D(x_:z_)
    B1t_D        = B1_D-B1n_D

! A field aligned flow extrapolation is more conform the MHD equations.
! However, in the current heating scenario the nonzero inflow will create
! a nonzero Ew*vr flux and therefor change the coronal heating properties.
! We leave this boundary condition here for future heating explorations.
!    !\
!    ! Update BCs for induction field::
!    !/
!    VarsGhostFace_V(Bx_:Bz_) = B1t_D(x_:z_)
!
!    !\
!    ! Update BCs for velocity
!    !/
!    FullB_D(x_:z_) = VarsGhostFace_V(Bx_:Bz_) + B0Face_D(x_:z_)
!    FullBnorm = sqrt(sum(FullB_D**2))
!    if(FullBnorm<cTolerance)then
!       ! HD: extrapolate in radial direction
!       UFace_D = dot_product(U_D,RFace_D)*RFace_D
!    else
!       ! MHD: extrapolate in field line direction
!       FullBunit_D = FullB_D/FullBnorm
!       UFace_D = dot_product(U_D,FullBunit_D)*FullBunit_D
!    endif
!    UdotRFace = dot_product(UFace_D,RFace_D)
!    if(UdotRFace < 0.0)then
!       VarsGhostFace_V(Ux_:Uz_) = 0.0
!    else
!       VarsGhostFace_V(Ux_:Uz_) = UFace_D(x_:z_)
!    endif

    !\
    ! Update BCs for velocity and induction field::
    !/
    VarsGhostFace_V(Ux_:Uz_) = -U_D(x_:z_)
    VarsGhostFace_V(Bx_:Bz_) = B1t_D(x_:z_)!-B1n_D(x_:z_)
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
    VarsGhostFace_V(Ew_)  = VarsGhostFace_V(rho_)*Prescell/Denscell &
         *(1.0/(Gammacell-cOne)-inv_gm1)
    !\
    ! Apply corotation
    !/
    if (.not.UseRotatingFrame) then
       VarsGhostFace_V(Ux_) = VarsGhostFace_V(Ux_) &
            - 2*OmegaBody*FaceCoords_D(y_)
       VarsGhostFace_V(Uy_) = VarsGhostFace_V(Uy_) &
            + 2*OmegaBody*FaceCoords_D(x_)
    end if

    if(UseHyperbolicDivb)then
       VarsGhostFace_V(Hyp_) = 0.0
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
  subroutine user_initial_perturbation
    use ModMain,      ONLY: nI, nJ, nK, nBLK, unusedBLK, UseHyperbolicDivb
    use ModIO,        ONLY: restart
    use ModVarIndexes
    use ModAdvance,   ONLY: State_VGB, Hyp_
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
    !------------------------------------------------------------------------
    call set_oktest('user_initial_perturbation',oktest,oktest_me)

    select case(TypeGeometry)
    case('cartesian')
       Rmax=sqrt(x2**2+y2**2+z2**2)
    case('spherical_lnr')
       Rmax=exp(XyzMax_D(1))
    case('spherical')
       Rmax=XyzMax_D(1)
    case default
       call stop_mpi('user_initial_perturbation is not implemented ' &
                     //'for geometry= '//TypeGeometry)
    end select
    Rmax=max(2.1E+01,Rmax)

    do iBLK=1,nBLK
       if(unusedBLK(iBLK))CYCLE
       if ((.not.restart)) then   
          do k=1,nK;do j=1,nJ; do i=1,nI
             xx = x_BLK(i,j,k,iBLK)
             yy = y_BLK(i,j,k,iBLK)
             zz = z_BLK(i,j,k,iBLK)
             RR = sqrt(xx**2+yy**2+zz**2+cTolerance**2)
             ROne  = max(cOne,RR)
             State_VGB(Bx_      ,i,j,k,iBLK) = cZero
             State_VGB(By_      ,i,j,k,iBLK) = cZero
             State_VGB(Bz_      ,i,j,k,iBLK) = cZero
             if(UseHyperbolicDivb)then
                State_VGB(Hyp_  ,i,j,k,iBLK) = 0.0
             end if
             call get_plasma_parameters_cell(xx,yy,zz,RR,&
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
       endif
       !\
       ! Update the total energy::
       !/
       call calc_energy_cell(iBLK)
    end do
  end subroutine user_initial_perturbation

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
    real    :: Gammacell
    !------------------------------------------------------------------------
    call update_states_MHD(iStage,iBlock)
    !\
    ! Begin update of pressure and relaxation energy::
    !/
    !  if (iStage/=nStage) return
    do k=1,nK; do j=1,nJ; do i=1,nI
       call get_gamma_emp(x_BLK(i,j,k,iBlock),y_BLK(i,j,k,iBlock), &
            z_BLK(i,j,k,iBlock),Gammacell)
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
    use ModAMR,ONLY:InitialRefineType, initial_refine_levels
    use ModNumConst
    use ModAdvance,ONLY:&
         State_VGB,Bx_,By_,Bz_,B0xCell_BLK,B0yCell_BLK,B0zCell_BLK
    use ModGeometry
    use ModPhysics,ONLY:rBody
    use ModOctree
    implicit none

    logical,intent(out) :: refineBlock, found
    integer, intent(in) :: lev
    real, intent(in)    :: DxBlock
    real, intent(in)    :: xCenter,yCenter,zCenter,rCenter
    real, intent(in)    :: minx,miny,minz,minR
    real, intent(in)    :: maxx,maxy,maxz,maxR
    integer, intent(in) :: iBLK

    character (len=*), parameter :: Name='user_specify_initial_refinement'
    real :: BDotRMin,BDotRMax,CritR, MinRBlk
    integer :: i,j,k, levmin, levelHS, levBLK
    !------------------------------------------------------------------------

    select case (InitialRefineType)
    case ('helio_init')

       select case(TypeGeometry)
       case('cartesian')
          if(.not.time_loop)then
             !refine to have resolution not worse 4.0 and
             !refine the body intersecting blocks
             ! Refine all blocks time through (starting with 1 block)
             if (lev <= 4) then
                refineBlock = .true.
             else
                critr=(XyzMax_D(1)-XyzMin_D(1))/(2.0**real(lev-2))
                if ( rCenter < 1.10*rBody + critr ) then
                   refineBlock = .true.
                else
                   refineBlock = .false.
                end if
             endif
          elseif(dx_BLK(iBLK)<0.20.or.far_field_BCs_BLK(iBLK))then
             refineBlock = .false. !Do not refine body or outer boundary
          else
             call refine_heliosheath
          end if
       case('spherical_lnr')
          levBLK=global_block_ptrs(iBLK,iProc+1)%ptr%LEV
          levmin=1
          levelHS=2
          if(.not.time_loop)then
             !refine to have resolution not worse than levmin
             if(lev<=levmin)then
                refineBlock = .true.
             else
                refineBlock = .false.
             end if
          elseif( levBLK<=levmin ) then
             refineBlock = .true.
          elseif( levBLK<=levelHS-1 )then
             MinRBlk = exp(XyzStart_Blk(1,iBlk)-0.5*Dx_Blk(iBlk))
             CritR=rBody + (exp(XyzMax_D(1))-exp(XyzMin_D(1))) &
                  /(2.0**real(levBLK+2-levmin))
             if( MinRBlk < CritR )then
                refineBlock = .true.
             else
                call refine_heliosheath
             endif
          else
             refineBlock=.false.
          end if
       case default
          call stop_mpi('user_specify_initial_refinement is ' &
               //'not implemented for geometry= '//TypeGeometry)
       end select
       found=.true.
    end select

    contains
      subroutine refine_heliosheath

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
      end subroutine refine_heliosheath
  end subroutine user_specify_initial_refinement

end module ModUser
