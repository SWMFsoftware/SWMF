!^CFG COPYRIGHT UM
!============================================================================

module ModUser
  use ModExpansionFactors
  use ModMagnetogram
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_set_ics,                    &
       IMPLEMENTED3 => user_get_b0,                     &
       IMPLEMENTED4 => user_update_states,              &
       IMPLEMENTED5 => user_specify_initial_refinement, &
       IMPLEMENTED6 => user_set_outerbcs,               &
       IMPLEMENTED7 => user_set_plot_var

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.1
  character (len=*), parameter :: &
       NameUserModule = 'EMPIRICAL SC on spherical grids'

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
    use ModMain,     ONLY: UseHyperbolicDivb, UseRotatingFrame
    use ModNumConst, ONLY: cTolerance
    use ModGeometry, ONLY: x_Blk, y_Blk, z_Blk, r_Blk
    use ModPhysics,  ONLY: inv_gm1, OmegaBody
    use ModSize
    use ModVarIndexes
    use ModFaceValue, ONLY: UseLogRhoLimiter
    implicit none

    integer,          intent(in)  :: iBlock, iSide
    character(len=20),intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    character (len=*), parameter :: Name='user_set_outerbcs'

    integer :: i, j, k, iMin, jMin, kMin, iMax, jMax, kMax
    real, dimension(nDim) :: B1_D, B1t_D, FullB_D, FullBt_D
    real, dimension(nDim) :: RhoU_D, RhoUt_D
    real, dimension(nDim) :: Normal_D
    real :: xx, yy, zz, RR, Dens, Pres, Gamma
    real :: B1n, FullBn, RhoUn, r2RhoUn
    !------------------------------------------------------------------------
    IsFound = .true.

    if(iSide/=East_) call stop_mpi('Wrong iSide in user_set_outerBCs')

    iMin=1-gcn;iMax=0
    jMin=1-gcn;jMax=nJ+gcn
    kMin=1-gcn;kMax=nK+gcn

    do k=kMin,kMax; do j=jMin,jMax
       xx = 0.5*sum(x_Blk(0:1,j,k,iBlock))
       yy = 0.5*sum(y_Blk(0:1,j,k,iBlock))
       zz = 0.5*sum(z_Blk(0:1,j,k,iBlock))
       RR = 0.5*sum(r_Blk(0:1,j,k,iBlock))

       call plasma_at_base(xx,yy,zz,RR,Dens,Pres)
       call my_gamma_emp(xx,yy,zz,Gamma)

       if(UseLogRhoLimiter)then
          State_VGB(Rho_,1,j,k,iBlock)=log(State_VGB(Rho_,1,j,k,iBlock))
          Dens = log(Dens)
       end if
       State_VGB(Rho_,iMax,j,k,iBlock) = 2.0*Dens &
            - State_VGB(Rho_,iMax+1,j,k,iBlock)
       do i=iMax-1,iMin,-1
          State_VGB(Rho_,i,j,k,iBlock) = 2.0*State_VGB(Rho_,i+1,j,k,iBlock) &
               - State_VGB(Rho_,i+2,j,k,iBlock)
       end do
       if(UseLogRhoLimiter)then
          State_VGB(Rho_,-1:1,j,k,iBlock) = &
               exp(State_VGB(Rho_,-1:1,j,k,iBlock))
          Dens = exp(Dens)
       end if

       do i=iMin,iMax
          State_VGB(P_,i,j,k,iBlock) = State_VGB(Rho_,i,j,k,iBlock) &
               *Pres/Dens
          State_VGB(Ew_,i,j,k,iBlock) = State_VGB(Rho_,i,j,k,iBlock) &
               *Pres/Dens*(1.0/(Gamma-1.0)-inv_gm1)
       end do

       Normal_D(x_)=x_Blk(1,j,k,iBlock)/r_Blk(1,j,k,iBlock)
       Normal_D(y_)=y_Blk(1,j,k,iBlock)/r_Blk(1,j,k,iBlock)
       Normal_D(z_)=z_Blk(1,j,k,iBlock)/r_Blk(1,j,k,iBlock)

       B1_D  = State_VGB(Bx_:Bz_,1,j,k,iBlock)
       B1n   = dot_product(Normal_D,B1_D)
       B1t_D = B1_D - B1n*Normal_D

       do i=iMin,iMax
          State_VGB(Bx_:Bz_,i,j,k,iBlock) = B1t_D
       end do

       RhoU_D = State_VGB(RhoUx_:RhoUz_,1,j,k,iBlock)
       r2RhoUn = r_Blk(1,j,k,iBlock)**2*dot_product(Normal_D,RhoU_D)

       do i=iMin,iMax
          FullB_D(x_) = State_VGB(Bx_,i,j,k,iBlock) &
               + B0xCell_BLK(i,j,k,iBlock)
          FullB_D(y_) = State_VGB(By_,i,j,k,iBlock) &
               + B0yCell_BLK(i,j,k,iBlock)
          FullB_D(z_) = State_VGB(Bz_,i,j,k,iBlock) &
               + B0zCell_BLK(i,j,k,iBlock)
          FullBn = dot_product(Normal_D,FullB_D)

          if(abs(FullBn) < cTolerance)then
             ! No flow at polarity inversion lines
             State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
          else
             FullBt_D = FullB_D - FullBn*Normal_D
             ! No backflow
             RhoUn = max(0.0,r2RhoUn/r_Blk(i,j,k,iBlock)**2)
             RhoUt_D = (RhoUn/FullBn)*FullBt_D
             State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = &
                  RhoUt_D + RhoUn*Normal_D
          end if
       end do
    end do; end do

    !\
    ! Apply corotation
    !/
    if(.not.UseRotatingFrame)then
       State_VGB(RhoUx_,iMin:iMax,jMin:jMax,kMin:kMax,iBlock) = &
            State_VGB(RhoUx_,iMin:iMax,jMin:jMax,kMin:kMax,iBlock) &
            -2.0*OmegaBody*y_Blk(iMin:iMax,jMin:jMax,kMin:kMax,iBlock) &
            *State_VGB(Rho_,iMin:iMax,jMin:jMax,kMin:kMax,iBlock)
       State_VGB(RhoUy_,iMin:iMax,jMin:jMax,kMin:kMax,iBlock) = &
            State_VGB(RhoUy_,iMin:iMax,jMin:jMax,kMin:kMax,iBlock) &
            -2.0*OmegaBody*x_Blk(iMin:iMax,jMin:jMax,kMin:kMax,iBlock) &
            *State_VGB(Rho_,iMin:iMax,jMin:jMax,kMin:kMax,iBlock)
    end if

    if(UseHyperbolicDivb)then
       State_VGB(Hyp_,iMin:iMax,jMin:jMax,kMin:kMax,iBlock) = 0.0
    end if

  end subroutine user_set_outerbcs

  !==========================================================================
  subroutine user_set_ics
    use ModMain,      ONLY: globalBLK, nI, nJ, nK, UseHyperbolicDivb
    use ModVarIndexes
    use ModAdvance,   ONLY: State_VGB, Hyp_
    use ModNumConst,  ONLY: cTolerance
    use ModPhysics,   ONLY: inv_gm1, GBody
    use ModGeometry, ONLY: x_Blk, y_Blk, z_Blk, r_Blk
    implicit none

    integer :: iBLK, i, j, k
    integer :: iCrita, iCritb, IterCount
    logical :: oktest,oktest_me
    real :: xx, yy, zz, rr, Dens, Pres, Gamma
    real :: Ur(1:nI), Ur0, Ur1, del, UBase, rTransonic, UEscape

    real, parameter :: Epsilon = 1.0e-6
    !------------------------------------------------------------------------
    call set_oktest('user_set_ics',oktest,oktest_me)

    iBLK = globalBLK

    UEscape = sqrt(-GBody*2.0)

    !\
    ! Initialize MHD wind with Parker's solution
    !/
    rTransonic = 0.25*UEscape**2
    if(.not.(rTransonic>exp(1.0))) call stop_mpi('sonic point inside Sun')

    if(r_BLK(1,1,1,iBLK) > rTransonic)then
       ! sonic point below subdomain
       iCrita = 0
       iCritb = 1
    elseif(r_BLK(nI,1,1,iBLK) <= rTransonic)then
       ! sonic point above subdomain
       iCrita = nI
       iCritb = nI + 1
    else
       ! sonic point inside subdomain
       i = 1
       do
          i = i + 1
          if(r_BLK(i,1,1,iBLK) > rTransonic .or. i > nI) EXIT
       end do
       iCrita = i - 1
       iCritb = i
    end if

    !\
    ! start construction of Parker solution from
    ! sonic point outwards
    !
    ! Part I: from inside sonic point to base
    !/
    Ur0 = 1.0
    do i=iCrita,1,-1
       rr = r_BLK(i,1,1,iBLK)
       IterCount = 0
       do
          IterCount = IterCount + 1
          Ur1 = (UEscape**2/(4.0*rr))**2 &
               *exp(0.5*(Ur0**2 + 3.0 - UEscape**2/rr))
          del = abs(Ur1 - Ur0)
          if(del < Epsilon)then
             Ur(i) = Ur1
             Ur0 = Ur1
             EXIT
          else
             if(IterCount < 1000)then
                Ur0 = Ur1
                CYCLE
             else
                call stop_mpi('PARKER > 1000 it.')
             end if
          end if
       end do
    end do

    !\
    ! Part II: from outside sonic point to infinity
    !/
    Ur0 = 1.0
    do i = iCritb,nI,+1
       rr = r_BLK(i,1,1,iBLK)
       IterCount = 0
       do
          IterCount = IterCount + 1
          Ur1 = sqrt(UEscape**2/rr - 3.0 &
               + 2.0*log(16.0*Ur0*rr**2/UEscape**4))
          del = abs(Ur1 - Ur0)
          if(del < Epsilon)then
             Ur(i) = Ur1
             Ur0 = Ur1
             EXIT
          else
             if(IterCount < 1000)then
                Ur0 = Ur1
                CYCLE
             else
                call stop_mpi('PARKER > 1000 it.')
             end if
          end if
       end do
    end do
    !\
    ! construct solution which obeys
    !   rho x u_r x r^2 = constant
    !/
    UBase=rTransonic**2*exp(1.5-2.0*rTransonic)

    do k=1,nK;do j=1,nJ; do i=1,nI
       xx = x_BLK(i,j,k,iBLK)
       yy = y_BLK(i,j,k,iBLK)
       zz = z_BLK(i,j,k,iBLK)
       rr = r_BLK(i,j,k,iBLK)
       call plasma_at_base(xx,yy,zz,rr,Dens,Pres)
       State_VGB(Rho_,i,j,k,iBLK) = Dens*UBase/(rr**2*Ur(i))
       State_VGB(RhoUx_,i,j,k,iBLK) = State_VGB(Rho_,i,j,k,iBLK) &
            *Ur(i)*xx/rr
       State_VGB(RhoUy_,i,j,k,iBLK) = State_VGB(Rho_,i,j,k,iBLK) &
            *Ur(i)*yy/rr
       State_VGB(RhoUz_,i,j,k,iBLK) = State_VGB(Rho_,i,j,k,iBLK) &
            *Ur(i)*zz/rr
       State_VGB(Bx_:Bz_,i,j,k,iBLK) = 0.0
       if(UseHyperbolicDivb)then
          State_VGB(Hyp_,i,j,k,iBLK) = 0.0
       end if
       State_VGB(P_,i,j,k,iBLK) = State_VGB(Rho_,i,j,k,iBLK)*Pres/Dens
       call my_gamma_emp(xx,yy,zz,Gamma)
       State_VGB(Ew_,i,j,k,iBLK) = State_VGB(P_,i,j,k,iBLK) &
            *(1.0/(Gamma-1.0)-inv_gm1) 
    end do;end do; end do

  end subroutine user_set_ics

  !==========================================================================
  subroutine plasma_at_base(xx,yy,zz,rr,Dens,Pres)
    use ModPhysics,  ONLY: inv_g, BodyTdim_I
    implicit none

    real, intent(in) :: xx, yy, zz, rr
    real, intent(out) :: Dens, Pres

    real :: UFinal, URatio
    !------------------------------------------------------------------------
!    call get_bernoulli_integral(xx/rr,yy/rr,zz/rr,UFinal)
!    URatio = UFinal/UMin
!    Dens  = (1.0/URatio)**2
!    Pres  = inv_g*Dens*T0/(min(URatio,2.0)*BodyTdim_I(1))
    Dens  = 1.0
    Pres  = inv_g*Dens*T0/BodyTdim_I(1)

  end subroutine plasma_at_base

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
    use ModAdvance, ONLY: State_VGB
    use ModMain,    ONLY: nStage
    use ModPhysics, ONLY: inv_gm1
    use ModGeometry,ONLY: x_BLK, y_BLK, z_BLK, R_BLK
    use ModEnergy,  ONLY: calc_energy_cell
    use ModExpansionFactors, ONLY:  GammaSS
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
       call my_gamma_emp(x_BLK(i,j,k,iBlock),y_BLK(i,j,k,iBlock), &
            z_BLK(i,j,k,iBlock),Gammacell)
       call correct_gamma(i,j,k,iBlock,Gammacell)

       State_VGB(P_,i,j,k,iBlock) = (Gammacell-1.0) &
            *(inv_gm1*State_VGB(P_,i,j,k,iBlock) + State_VGB(Ew_,i,j,k,iBlock))
       State_VGB(Ew_,i,j,k,iBlock) = State_VGB(P_,i,j,k,iBlock) &
            *(1.0/(Gammacell-1.0)-inv_gm1)
    end do; end do; end do
    call calc_energy_cell(iBlock)
    !\
    ! End update of pressure and relaxation energy::
    !/
  end subroutine user_update_states

  !==========================================================================
  subroutine correct_gamma(i,j,k,iBlock,Gamma)
    use ModAdvance, ONLY: State_VGB, B0xCell_BLK, B0yCell_BLK, B0zCell_BLK
    use ModGeometry,ONLY: R_BLK
    use ModVarIndexes
    implicit none

    integer, intent(in) :: i, j, k, iBlock
    real, intent(inout) :: Gamma
    !------------------------------------------------------------------------
    !\
    ! Change the empirical gamma near the current sheet based on the
    ! plasma beta
    !/
    if(R_BLK(i,j,k,iBlock) > Rs_PFSSM) &
         Gamma = Gamma - (Gamma - GammaSS)*max(0.0, &
         -1.0 + 2.0*State_VGB(P_,i,j,k,iBlock)/(State_VGB(P_,i,j,k,iBlock) &
         +((State_VGB(Bx_,i,j,k,iBlock) + B0xCell_BLK(i,j,k,iBlock))**2 &
         + (State_VGB(By_,i,j,k,iBlock) + B0yCell_BLK(i,j,k,iBlock))**2 &
         + (State_VGB(Bz_,i,j,k,iBlock) + B0zCell_BLK(i,j,k,iBlock))**2) &
         *0.25*(R_BLK(i,j,k,iBlock)/Rs_PFSSM)**1.5))

  end subroutine correct_gamma

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
       case('spherical_lnr')
          ! assumes 3x4x2 blocks on root level
          levBLK=global_block_ptrs(iBLK,iProc+1)%ptr%LEV
          levmin=3
          levelHS=4
          if(.not.time_loop)then
             !refine to have resolution not worse than levmin
             if(lev<=levmin)then
                refineBlock = .true.
             else
                MinRBlk = exp(XyzStart_Blk(1,iBlk)-0.5*Dx_Blk(iBlk))
                CritR=rBody + (exp(XyzMax_D(1))-exp(XyzMin_D(1))) &
                     /(2.0**(lev+2-levmin))
                if( MinRBlk < CritR )then
                   refineBlock = .true.
                else
                   refineBlock = .false.
                endif
             end if
          elseif( levBLK<=levelHS-1 )then
             call refine_heliosheath
          else
             refineBlock = .false.
          end if
       case default
          call stop_mpi('user_specify_initial_refinement is ' &
               //'not implemented for geometry= '//TypeGeometry)
       end select
       found=.true.
    end select

    contains
      subroutine refine_heliosheath

        BDotRMin=0.0
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
        BDotRMax=0.0
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

  !=========================================================================
  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModAdvance, ONLY: State_VGB, B0xCell_BLK, B0yCell_BLK, B0zCell_BLK
    use ModSize, ONLY: nI, nJ, nK
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, r_BLK
    use ModVarIndexes
    use ModExpansionFactors, ONLY:  GammaSS
    implicit none

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

    character (len=*), parameter :: Name='user_set_plot_var'

    integer :: i,j,k
    real :: Gamma
    !------------------------------------------------------------------------

    IsFound=.true.

    select case(NameVar)
    case('g0')
       do k=-1,nK+2; do j=-1,nJ+2; do i=-1,nI+2
          call my_gamma_emp(x_BLK(i,j,k,iBlock),y_BLK(i,j,k,iBlock), &
               z_BLK(i,j,k,iBlock),Gamma)

          PlotVar_G(i,j,k) = Gamma
       end do; end do; end do
    case('g1')
       do k=-1,nK+2; do j=-1,nJ+2; do i=-1,nI+2
          call my_gamma_emp(x_BLK(i,j,k,iBlock),y_BLK(i,j,k,iBlock), &
               z_BLK(i,j,k,iBlock),Gamma)
          call correct_gamma(i,j,k,iBlock,Gamma)

          PlotVar_G(i,j,k) = Gamma
       end do; end do; end do
    end select

    UsePlotVarBody=.true.
    PlotVarBody=1.0

  end subroutine user_set_plot_var

  !==========================================================================
  subroutine my_gamma_emp(xx,yy,zz,gammaOut)

    ! Subroutine my_gamma_emp
    ! Provides the distribution of the polytropic index, complying with
    ! the WSA or Fisk semi-empirical models

    use ModExpansionFactors
    use ModNumConst
    implicit none

    real, intent(in) :: xx,yy,zz
    real, intent(out)   :: gammaOut 
    real :: RR,Uf,BernoulliFactor
    real, parameter :: gammaIH=1.5
    real, parameter :: R1=2.50,R2=12.50
    integer,parameter::nPowerIndex=2
    !------------------------------------------------------------------------
    !\
    ! Calculate cell-centered spherical coordinates::
    RR   = sqrt(xx**2+yy**2+zz**2)
    !\
    ! Avoid calculating inside a critical radius = 0.5*Rsun
    !/
    if (RR <max(Ro_PFSSM-dR*nRExt,0.90*Ro_PFSSM)) then 
       gammaOut= gammaSS
       RETURN
    end if

    ! Calculate gamma
    if(RR >= R2)then
       gammaOut=gammaIH
    else if(RR >= R1)then
       gammaOut=gammaSS+(RR-R1)*(gammaIH-gammaSS)/(R2-R1)
    else
       call get_bernoulli_integral(xx,yy,zz,Uf)
       BernoulliFactor=(cHalf*Uf**2+cSunGravitySI)/&
            (T0*cBoltzmann/cProtonMass)&
            *(R1-RR)*&
            & (Ro_PFSSM/RR)**nPowerIndex/ (R1-Ro_PFSSM)+ GammaSS&
            &/(GammaSS-cOne)*(cOne- (R1-RR)*(Ro_PFSSM/RR)&
            &**nPowerIndex/ (R1-Ro_PFSSM))
       gammaOut = BernoulliFactor/(BernoulliFactor-cOne)
    end if

  end subroutine my_gamma_emp

end module ModUser
