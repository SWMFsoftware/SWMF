!^CFG COPYRIGHT UM
!^CFG FILE USERFILES
!========================================================================
!========================================================================
!  This file contains routines that let the user make customizations to
!  BATSRUS without having to make changes in the "kernel" of the code.
!  Each of the routines is accessed in a run but setting the appropriate
!  input flag in the PARAM.in file.  Also included here is a ModUser
!  module.  Variables that are needed in more than one user routine (in
!  otherwords, if they must be global) can be defined here and accessed
!  in the other user subroutines.
!
!  In each of the routines there are sections for the user to place 
!  variable definitions and to put code to compute the desired quantity.
!  Note that in many routines (user_set_inner_BCs for example) the 
!  solution to a problem will be greatly affected by changes the user
!  introduces.  Often we suggest using one of the BATSRUS defaults rather
!  then coding a unique option.  However, if the user chooses to create
!  a user routine, we suggest that you follow the format of the examples
!  give below very carefully.  At the very least, the user should take
!  care to set the global variables that the routine requires to be set
!  and no others.  This will happen if the user adhears to the examples 
!  given and does not modify the looping structure or the sections of the 
!  subroutine which are documented as parts not to change.
!========================================================================
!========================================================================

!========================================================================
!  MODULE MODUSER
!========================================================================
!\
! This module is where the user should define any variables that are
! needed in more than a single routine.  In otherwords, variables which
! need to be global to all the subroutines in this file.  
!/
Module ModUser
  use ModNumConst, ONLY: cHalf,cTwo,cThree,&
       cFour,cE1,cHundred,cHundredth,cZero,&
       cOne
  use ModMain,     ONLY: UseUserB0,UseUserHeating
  use ModSize,     ONLY: nI,nJ,nK,gcn,nBLK
  implicit none
  !\
  ! PFSSM related variables::
  !/
  integer, parameter:: N_PFSSM=30
  integer, parameter:: UNIT_PFSSM=30
  real, parameter:: unitPFSSM_B=cHundredth
  integer:: iHead_PFSSM
  logical:: DoFirst=.true.
  logical:: DoSecond=.true.
  character (LEN=32):: File_PFSSM
  character (LEN=80):: Head_PFSSM
  real, dimension(N_PFSSM+1,N_PFSSM+1):: g_nm,h_nm
  real, dimension(N_PFSSM+1):: factRatio1
  !\
  ! Rs_PFSSM=2.5; Ro_PFSSM=1.0::
  !/
  real:: R_PFSSM,Rs_PFSSM,Ro_PFSSM
  !\
  ! Parameters related to the empirical heating::
  !/
  logical:: DoStaticICs=.false.
  real:: InvH0,TemRatio
  real:: DegFrm1,Tnot=1.0
  real:: MaxB0,Bnot
  !\
  ! Variables related to energies::
  !/
  real:: Emag_0, Ekin_0, Ethe_0
  !\
  ! Gibson & Low 1998 related variables::
  !/
  logical:: DoFirst_GL=.true.
  real:: Mrope_GL98
  !\
  ! The following are needed in user_sources::
  !/
  real, dimension(1:nI,1:nJ,1:nK):: &
       Srho,SrhoUx,SrhoUy,SrhoUz,SBx,SBy,SBz,Sp,SEr,SE
end Module ModUser
!========================================================================
!  SUBROUTINE USER_READ_INPUTS
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
  use ModIO, ONLY: write_prefix, write_myname, iUnitOut
  use ModUser,      ONLY: DoStaticICs,UseUserB0,Ro_PFSSM,&
       Rs_PFSSM,File_PFSSM,iHead_PFSSM,TemRatio,DegFrm1, &
       Tnot,Bnot
  implicit none

  integer:: i
  character (len=100) :: NameCommand
  !---------------------------------------------------------------------------

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
        call read_var('UseUserSetPhysConst'     ,UseUserSetPhysConst)
        call read_var('UseUserUpdateStates'     ,UseUserUpdateStates)
     case("#USEUSERHEATING")
	call read_var('UseUserHeating'          ,UseUserHeating)
     case("#PFSSM")
        call read_var('UseUserB0'  ,UseUserB0)
        call read_var('Ro_PFSSM'   ,Ro_PFSSM)
        call read_var('Rs_PFSSM'   ,Rs_PFSSM)
        call read_var('File_PFSSM' ,File_PFSSM)
        call read_var('iHead_PFSSM',iHead_PFSSM)
        call read_var('dt_UpdateB0',dt_UpdateB0)
     case("#AWHEAT")
        call read_var('Bnot     '  ,Bnot)
        call read_var('Tnot     '  ,Tnot)
        call read_var('DegFrm1  '  ,DegFrm1)
        call read_var('TemRatio '  ,TemRatio)
        call read_var('DoStaticICs',DoStaticICs)
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

subroutine user_heat_source
  use ModMain,     ONLY: globalBLK
  use ModAdvance,  ONLY: qheat_BLK
  use ModNumConst, ONLY: cZero
  implicit none
  qheat_BLK(:,:,:,globalBLK) = cZero
  call heat_source_averages
end subroutine user_heat_source

!========================================================================
!  SET_EXTRA_BOUNDARY_CELLS
!  Allows to define boundary conditions at the user defined boundary.
!========================================================================
subroutine set_extra_boundary_cells(iBLK)
  use ModGeometry	
  implicit none
  integer, intent(in):: iBLK
!  SHOULD define IsBoundaryCell_GI(:,:,:,ExtraBc_) using
!  a boundary condition for iBLK block
!  EXAMPLE: OUTER SPHERICAL BOUNDARY of radius of 100.
!  IsBoundaryCell_GI(:,:,:,ExtraBc_) = R_BLK(:,:,:,iBLK)<100.
end subroutine set_extra_boundary_cells
!========================================================================

!========================================================================
!  SUBROUTINE USER_SET_INNER_BCS
!========================================================================
!\
! This subroutine allows the user to apply boundary conditions to the inner
! body which are problem specific and cannot be created using the predefined
! options in BATSRUS.  The available options in BATSRUS have been designed
! to be self consistent and reasonably robust.  We generally recommend that
! you use on of those or a variant that is very close to one of them.  They
! can be considered reasonably safe.
!
! An example of a reasonable variant would be to use a modification of the
! "ionosphere" boundary where the density is fixed at the boundary to a 
! value that is a function of latitude.
!
! This routine is called for a single inner boundary face.  Since BATSRUS is
! is block cartesian, the values inside the boundary face must be passed back
! in cartesian coordinates.  Values that must be set are:
!
!  RhoFaceInside, pFaceInside, VxFaceInside, VyFaceInside, VzFaceInside
!  BxFaceInside, ByFaceInside, BzFaceInside
!
! Typically the boundary conditions are applied for the spherical coordinates
! and then transformed to the cartesian ones.
!
! As with all user subroutines, the variables declared in ModUser are 
! available here.  Again, as with other user subroutines DO NOT MODIFY 
! ANY GLOBAL VARIABLE DEFINED IN THE MODULES INCLUDED IN THIS SUBROUTINE 
! UNLESS SPECIFIED!!
!/
subroutine user_face_bcs(iFace,jFace,kFace,iBlock,iSide,iBoundary,&
     iter,time_now,FaceCoords_D,VarsTrueFace_V,VarsGhostFace_V,   &
     B0Face_D,UseIonosphereHere,UseCorotationHere)
  use ModSize,       ONLY: nDim,East_,West_,South_,North_,Bot_,   &
       Top_
  use ModMain,       ONLY: UseUserHeating,UseInertial
  use ModVarIndexes, ONLY: &
 !       EnergyRL_,&     !^CFG UNCOMMENT IF ALWAVES
	rho_,Ux_,Uy_,Uz_,Bx_,By_,Bz_,P_

  use ModGeometry,   ONLY: R_BLK
  use ModAdvance,    ONLY: nFaceValueVars
  use ModPhysics,    ONLY: cosTHETAtilt,sinTHETAtilt,g,inv_g,     &
       inv_gm1,OMEGAbody
  use ModNumConst,   ONLY: cZero,cHalf,cOne,cTwo,cTolerance
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
       UseCorotationHere
  real, dimension(nFaceValueVars), intent(out)::       &
       VarsGhostFace_V
  !\
  ! User declared local variables go here::
  !/
  integer:: iCell,jCell,kCell
  real:: XFace,YFace,ZFace,RFace,dRFace
  real:: VxFaceOutside,VyFaceOutside,VzFaceOutside
  real:: BxFaceOutside,ByFaceOutside,BzFaceOutside
  real:: VrFaceOutside,VthetaFaceOutside,VphiFaceOutside,&
         VrFaceInside,VthetaFaceInside,VphiFaceInside,   &
         BrFaceOutside,BthetaFaceOutside,BphiFaceOutside,&
         BrFaceInside,BthetaFaceInside,BphiFaceInside
  real:: cosTheta,sinTheta,cosPhi,sinPhi
  real, dimension(1:3):: location,v_phi
  real:: XFaceT,YFaceT,ZFaceT,sin2Theta_coronal_hole
  real:: cosThetaT,sinThetaT,cosPhiT,sinPhiT
  real:: BodyDens,BodyPres,BodyGamma
  !
  !---------------------------------------------------------------------------
  !\
  ! Calculation of boundary conditions should start here::
  !/
  !---------------------------------------------------------------------------
  !
  XFace = FaceCoords_D(1)
  YFace = FaceCoords_D(2)
  ZFace = FaceCoords_D(3)
  VxFaceOutside = VarsTrueFace_V(Ux_)
  VyFaceOutside = VarsTrueFace_V(Uy_)
  VzFaceOutside = VarsTrueFace_V(Uz_)
  BxFaceOutside = VarsTrueFace_V(Bx_)
  ByFaceOutside = VarsTrueFace_V(By_)
  BzFaceOutside = VarsTrueFace_V(Bz_)
  !\
  ! Rotate to spherical coordinates
  !/
  RFace    = sqrt(XFace**2+YFace**2+ZFace**2)
  cosTheta = ZFace/RFace
  sinTheta = sqrt(XFace**2+YFace**2)/RFace
  cosPhi   = XFace/sqrt(XFace**2+YFace**2+cTolerance**2)
  sinPhi   = YFace/sqrt(XFace**2+YFace**2+cTolerance**2)
  VrFaceOutside = (VxFaceOutside*XFace      +&
       VyFaceOutside*YFace                  +&
       VzFaceOutside*ZFace)/RFace
  VthetaFaceOutside = ((VxFaceOutside*XFace +&
       VyFaceOutside*YFace)*ZFace           -&
       VzFaceOutside*(XFace**2+YFace**2))   /&
       (sqrt(XFace**2+YFace**2+cTolerance**2)*RFace)
  VphiFaceOutside = (VyFaceOutside*XFace    -&
       VxFaceOutside*YFace)*sinTheta        /&
       ((XFace**2+YFace**2+cTolerance**2)/RFace)
  BrFaceOutside = (BxFaceOutside*XFace      +&
       ByFaceOutside*YFace                  +&
       BzFaceOutside*ZFace)/RFace
  BthetaFaceOutside = ((BxFaceOutside*XFace +&
       ByFaceOutside*YFace)*ZFace           -&
       BzFaceOutside*(XFace**2+YFace**2))   /&
       (sqrt(XFace**2+YFace**2+cTolerance**2)*RFace)
  BphiFaceOutside = (ByFaceOutside*XFace    -&
       BxFaceOutside*YFace)*sinTheta        /&
       ((XFace**2+YFace**2+cTolerance**2)/RFace)
  if (.not.UseUserHeating) then
     VrFaceInside     = -VrFaceOutside
     VthetaFaceInside = -VthetaFaceOutside
     VphiFaceInside   = -VphiFaceOutside
     BrFaceInside     =  BrFaceOutside
     BthetaFaceInside =  BthetaFaceOutside
     BphiFaceInside   =  BphiFaceOutside
     !\
     ! Begin update of BCs for the mass density, EnergyRL,
     ! and pressure in the case of gravity::
     !/
     iCell = iFace; jCell = jFace; kCell = kFace
     select case(iSide)
     case(East_)
        iCell  = iFace
        dRFace = R_BLK(iCell  ,jCell,kCell,iBlock)-&
                 R_BLK(iCell-1,jCell,kCell,iBlock)
     case(West_)
        iCell  = iFace-1
        dRFace = R_BLK(iCell  ,jCell,kCell,iBlock)-&
                 R_BLK(iCell+1,jCell,kCell,iBlock)
     case(South_)
        jCell  = jFace
        dRFace = R_BLK(iCell,jCell  ,kCell,iBlock)-&
                 R_BLK(iCell,jCell-1,kCell,iBlock)
     case(North_)
        jCell  = jFace-1
        dRFace = R_BLK(iCell,jCell  ,kCell,iBlock)-&
                 R_BLK(iCell,jCell+1,kCell,iBlock)
     case(Bot_)
        kCell  = kFace
        dRFace = R_BLK(iCell,jCell,kCell  ,iBlock)-&
                 R_BLK(iCell,jCell,kCell-1,iBlock)
     case(Top_)
        kCell  = kFace-1
        dRFace = R_BLK(iCell,jCell,kCell  ,iBlock)-&
                 R_BLK(iCell,jCell,kCell+1,iBlock)
     case default
        call stop_mpi('User_face_BCs')
     end select
     call get_atmosphere_body
     VarsGhostFace_V(rho_     ) = -VarsTrueFace_V(rho_     )+&
          cTwo*BodyDens
     VarsGhostFace_V(P_   ) = -VarsTrueFace_V(P_   )+&
          cTwo*BodyPres
   !  VarsGhostFace_V(EnergyRL_) = -VarsTrueFace_V(EnergyRL_)+&  !^CFG UNCOMMENT IF ALWAVES
   !       cTwo*BodyPres*(cOne/(BodyGamma-cOne)-inv_gm1)         !^CFG UNCOMMENT IF ALWAVES
     !\
     ! End update of BCs for the mass density, EnergyRL,
     ! and pressure in the case of gravity::
     !/
  else
     RFace = cOne
     call coronal_hole_boundary(RFace,sin2Theta_coronal_hole)
     XFaceT    =  cosTHETAtilt*XFace+sinTHETAtilt*ZFace
     YFaceT    = YFace
     ZFaceT    = -sinTHETAtilt*XFace+cosTHETAtilt*ZFace
     cosThetaT = ZFaceT/RFace
     sinThetaT = sqrt(XFaceT**2+YFaceT**2)/RFace
     cosPhiT   = XFaceT/sqrt(XFaceT**2+YFaceT**2+cTolerance**2)
     sinPhiT   = YFaceT/sqrt(XFaceT**2+YFaceT**2+cTolerance**2)
     if (sinThetaT*sinThetaT.gt.sin2Theta_coronal_hole) then
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
        VarsGhostFace_V(P_   ) = VarsTrueFace_V(P_   )
 !       VarsGhostFace_V(EnergyRL_) = VarsTrueFace_V(EnergyRL_)!^CFG UNCOMMENT IF ALWAVES
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
        VarsGhostFace_V(rho_     ) = cOne
        VarsGhostFace_V(P_   ) = inv_g
       ! VarsGhostFace_V(EnergyRL_) = VarsTrueFace_V(EnergyRL_)!^CFG UNCOMMENT IF ALWAVES
     endif
  endif
  !\
  ! Rotate back to cartesian coordinates::
  !/
  VarsGhostFace_V(Ux_) = VrFaceInside*XFace/RFace+&
       VthetaFaceInside*cosTheta*cosPhi          -&
       VphiFaceInside*sinPhi 
  VarsGhostFace_V(Uy_) = VrFaceInside*YFace/RFace+&
       VthetaFaceInside*cosTheta*sinPhi          +&
       VphiFaceInside*cosPhi
  VarsGhostFace_V(Uz_) = VrFaceInside*ZFace/RFace-&
       VthetaFaceInside*sinTheta
  VarsGhostFace_V(Bx_) = BrFaceInside*XFace/RFace+&
       BthetaFaceInside*cosTheta*cosPhi          -&
       BphiFaceInside*sinPhi
  VarsGhostFace_V(By_) = BrFaceInside*YFace/RFace+&
       BthetaFaceInside*cosTheta*sinPhi          +&
       BphiFaceInside*cosPhi
  VarsGhostFace_V(Bz_) = BrFaceInside*ZFace/RFace-&
       BthetaFaceInside*sinTheta
  !\
  ! Apply corotation:: Currently works only for the first body.
  !/
  if (UseInertial) then
     !\
     ! The program is called which calculates the cartesian 
     ! corotation velocity::
     !/
     VarsGhostFace_V(Ux_) = VarsGhostFace_V(Ux_) -&
          cTwo*OMEGAbody*YFace
     VarsGhostFace_V(Uy_) = VarsGhostFace_V(Uy_) +&
          cTwo*OMEGAbody*XFace
  end if
  
Contains

  subroutine get_atmosphere_body
    !---------------------------------------------------------------------------
    !
    ! This module computes the boundary values for pressure and density
    ! for a politropic equation of state with variable gamma = 1+n/2,
    ! where n=n0+n1*T^pow::
    ! The subroutine is written by ILR on May 29, 2003.
    !
    !---------------------------------------------------------------------------
    !
    use ModVarIndexes, ONLY: Bx_,By_,Bz_
    use ModAdvance,    ONLY: B0xCell_BLK,B0yCell_BLK,&
         B0zCell_BLK,State_VGB
    use ModGeometry,   ONLY: x_BLK,y_BLK,z_BLK
    use ModNumConst,   ONLY: cHalf,cOne,cTwo,cFour
    use ModUser,       ONLY: TemRatio,MaxB0,DegFrm1
    use ModPhysics,    ONLY: g,inv_g
    implicit none
    real, parameter:: n0=cFour
    real:: XCell,YCell,ZCell,RCell
    real:: AAc,BBc,Fn1,Fg1
    real:: BodyTemp,BodyDegF
    real:: BrBody,TemMod
    !\
    ! Get the radial field component at the body::
    !/
    XCell = x_BLK(iCell,jCell,kCell,iBlock)
    YCell = y_BLK(iCell,jCell,kCell,iBlock)
    ZCell = z_BLK(iCell,jCell,kCell,iBlock)
    RCell = sqrt(XCell**2+YCell**2+ZCell**2)
    BrBody = (XCell*(B0xCell_BLK(iCell,jCell,kCell,iBlock)  +&
                   State_VGB(Bx_,iCell,jCell,kCell,iBlock)) +&
              YCell*(B0yCell_BLK(iCell,jCell,kCell,iBlock)  +&
                   State_VGB(By_,iCell,jCell,kCell,iBlock)) +&
              ZCell*(B0zCell_BLK(iCell,jCell,kCell,iBlock)  +&
                   State_VGB(Bz_,iCell,jCell,kCell,iBlock)))/&
              RCell
    TemMod = cOne+min(cOne,abs(BrBody*RCell**3/MaxB0))*&
          (TemRatio-cOne)
    Fg1 = (cOne+cHalf*n0+cHalf*cHalf*DegFrm1*TemMod)  /&
          (cOne+cHalf*n0+cHalf*cHalf*DegFrm1)
    Fn1 = DegFrm1*TemMod
    !\
    ! Find the solution for T at RCell::
    !/
    AAc = cFour*(cOne+cHalf*n0)/Fn1
    BBc = cFour*(cOne+cHalf*n0+cHalf*cHalf*Fn1)/RCell/Fn1
    BodyTemp = cTwo*BBc/(sqrt(AAc**2+cFour*BBc)+AAc)
    BodyDegF = n0+Fn1*BodyTemp
    !\
    ! Get the mass density & pressure at the Body::
    !/
    BodyPres  = inv_g*(BodyTemp**(cOne+cHalf*n0))*exp(-Fn1*&
         cHalf*(cOne-BodyTemp))/Fg1
    BodyDens  = g*BodyPres/BodyTemp/TemMod
    BodyGamma = (BodyDegF+cTwo)/BodyDegF
  end subroutine get_atmosphere_body
  
end subroutine user_face_bcs

subroutine user_calc_sources
  use ModVarIndexes, ONLY:&
        !EnergyRL_,& !^CFG UNCOMMENT IF ALWAVES 
	 rho_,rhoUx_,rhoUy_,rhoUz_,Bx_,By_,Bz_,P_
  use ModAdvance,    ONLY: Source_VC,Energy_,     &
       Flux_VX,Flux_VY,Flux_VZ 
  use ModUser,       ONLY: Srho,SrhoUx,SrhoUy,SrhoUz,&
       SBx,SBy,SBz,SP,SEr,SE
  use ModNumConst,   ONLY: cZero
  implicit none
  
  Srho   = cZero
  SrhoUx = cZero
  SrhoUy = cZero
  SrhoUz = cZero
  SBx    = cZero
  SBy    = cZero
  SBz    = cZero
  SP     = cZero
  SEr    = cZero
  SE     = cZero
  call user_sources
  Source_VC(rho_       ,:,:,:) = Srho  +Source_VC(rho_,:,:,:)
  Source_VC(rhoUx_     ,:,:,:) = SrhoUx+Source_VC(rhoUx_,:,:,:)
  Source_VC(rhoUy_     ,:,:,:) = SrhoUy+Source_VC(rhoUy_,:,:,:)
  Source_VC(rhoUz_     ,:,:,:) = SrhoUz+Source_VC(rhoUz_,:,:,:)
  Source_VC(Bx_        ,:,:,:) = SBx   +Source_VC(Bx_,:,:,:)
  Source_VC(By_        ,:,:,:) = SBy   +Source_VC(By_,:,:,:)
  Source_VC(Bz_        ,:,:,:) = SBz   +Source_VC(Bz_,:,:,:)
  Source_VC(P_     ,:,:,:) = SP    +Source_VC(P_,:,:,:)
  Source_VC(Energy_,:,:,:) = SE    +Source_VC(Energy_,:,:,:)
 ! Source_VC(EnergyRL_  ,:,:,:) = SEr   +Source_VC(EnergyRL_  ,:,:,:) !^CFG UNCOMMENT IF ALWAVES
end subroutine user_calc_sources
!========================================================================
!  SUBROUTINE USER_SOURCES
!========================================================================
!\
! This subroutine is used to calculate sources for the MHD equations.  The
! routine is called for each block separately so that the user would typically
! need only to code the source term calculation for a single block (in other
! words inside the the k,j,i loop below).  As with all user subroutines, the
! variables declared in ModUser are available here.  Again, as with other
! user subroutines DO NOT MODIFY ANY GLOBAL VARIABLE DEFINED IN THE MODULES
! INCLUDED IN THIS SUBROUTINE UNLESS SPECIFIED!!
!
! The user should load the global variables:
!      Srho,SrhoUx,SrhoUy,SrhoUz,SBx,SBy,SBz,SE,SP,SEr
!
! Note that SE (energy) and SP (pressure) must both be loaded if the code is 
! going to use both the primitive and the conservative MHD equation advance  
! (see the USER MANUAL and the DESIGN document).  If using only primitive SP 
! must be loaded.  If using only conservative SE must be loaded.  The safe
! approach is to load both.
!/
subroutine user_sources
  use ModMain,       ONLY: nI,nJ,nK,globalBLK,PROCtest,BLKtest,  &
       UseInertial,UseUserHeating,gcn
  use ModVarIndexes, ONLY: &
       !EnergyRL_, &    !^CFG UNCOMMENT IF ALWAVES
       rho_,rhoUx_,rhoUy_,rhoUz_,Bx_,By_,Bz_,P_
  use ModAdvance,    ONLY: State_VGB,StateOld_VCB,Source_VC,      &
       B0xCell_BLK,B0yCell_BLK,B0zCell_BLK,UDotFA_X,UDotFA_Y,    &
       UDotFA_Z, Flux_VX,Flux_VY,Flux_VZ,   &
       Theat0,qheat_BLK
  use ModGeometry,   ONLY: x_BLK,y_BLK,z_BLK,R_BLK,VolumeInverse_I
  use ModConst,      ONLY: cZero,cHalf,cOne,cTwo,cTolerance
  use ModPhysics,    ONLY: g,OMEGAbody,cosTHETAtilt,sinTHETAtilt,&
       Theat
  use ModProcMH,     ONLY: iProc 
  use ModUser,       ONLY: SrhoUx,SrhoUy,SrhoUz,SBx,SBy,SBz,SP,SE
  implicit none
  !\
  ! Variables required by this user subroutine::
  !/
  integer:: i,j,k
  logical:: oktest,oktest_me
  !
  !---------------------------------------------------------------------------
  !\
  ! Variable meanings:
  !   Srho: Source terms for the continuity equation
  !   SE,SP: Source terms for the energy (conservative) and presure
  !          (primative) equations
  !   SrhoUx,SrhoUy,SrhoUz:  Source terms for the momentum equation
  !   SBx,SBy,SBz:  Souce terms for the magnetic field equations
  !   SEr: Source terms for the relaxation energy equation
  !/
  !---------------------------------------------------------------------------
  !
  if (iProc==PROCtest.and.globalBLK==BLKtest) then
     call set_oktest('user_sources',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if
  do k=1,nK; do j=1,nJ; do i=1,nI
     !\
     ! User coding of the source terms should be inside this loop
     ! Note that user source terms should be added to the the
     ! values already loaded in the arrays. For example::
     ! Srho(i,j,k) = Srho(i,j,k) + user source term for rho
     !/
     if (UseUserHeating) call calc_OLD_heating
     if (.not.UseInertial) then
        !\
        ! Note that the SrhoUx, SrhoUy, SBx, SBy corresponds to Omega x B::
        ! This is used only in the rot. frame and is yet under investigation.
        !/
        SrhoUx(i,j,k) = SrhoUx(i,j,k) +&
             cTwo*OMEGAbody*State_VGB(rhoUy_,i,j,k,globalBLK)
        SrhoUy(i,j,k) = SrhoUy(i,j,k) -&
             cTwo*OMEGAbody*State_VGB(rhoUx_,i,j,k,globalBLK)
        SBx(i,j,k)    = SBx(i,j,k)    +&
             OMEGAbody*B0yCell_BLK(i,j,k,globalBLK)
        SBy(i,j,k)    = SBy(i,j,k)    -&
             OMEGAbody*B0xCell_BLK(i,j,k,globalBLK)
        SE(i,j,k)     = SE(i,j,k)     +&
             State_VGB(Bx_,i,j,k,globalBLK)        *&
             OMEGAbody*B0yCell_BLK(i,j,k,globalBLK)-&
             State_VGB(By_,i,j,k,globalBLK)        *&
             OMEGAbody*B0xCell_BLK(i,j,k,globalBLK)
     endif
  end do; end do; end do
  
Contains

    subroutine calc_OLD_heating
    implicit none
    real:: cosTheta,sinTheta,cosPhi,sinPhi,&
         sin2Theta_coronal_hole,XT,YT,ZT
    !\
    ! Compute Heliosphere source terms::
    !/
    XT =  cosTHETAtilt*x_BLK(i,j,k,globalBLK)+&
         sinTHETAtilt*z_BLK(i,j,k,globalBLK)
    YT =  y_BLK(i,j,k,globalBLK)
    ZT = -sinTHETAtilt*x_BLK(i,j,k,globalBLK)+&
          cosTHETAtilt*z_BLK(i,j,k,globalBLK)
    cosTheta = ZT/(R_BLK(i,j,k,globalBLK)+cTolerance)
    sinTheta = sqrt(XT**2+YT**2)             /&
         (R_BLK(i,j,k,globalBLK)+cTolerance)
    cosPhi = XT/sqrt(XT**2+YT**2+cTolerance**2)
    sinPhi = YT/sqrt(XT**2+YT**2+cTolerance**2)
    call coronal_hole_boundary(R_BLK(i,j,k,globalBLK),&
         sin2Theta_coronal_hole)
    if (sinTheta*sinTheta.lt.sin2Theta_coronal_hole) then
       Theat0(i,j,k) = Theat
    else
       Theat0(i,j,k) = cOne
    end if
    SP(i,j,k)  = SP(i,j,k)+&
        (g-cOne)*State_VGB(rho_,i,j,k,globalBLK)     *&
         qheat_BLK(i,j,k,globalBLK)*(Theat0(i,j,k)-g *&
         State_VGB(P_,i,j,k,globalBLK)           /&
         State_VGB(rho_,i,j,k,globalBLK))
    SE(i,j,k)  = SE(i,j,k)+&
         State_VGB(rho_,i,j,k,globalBLK)             *&
         qheat_BLK(i,j,k,globalBLK)*(Theat0(i,j,k)-g *&
         State_VGB(P_,i,j,k,globalBLK)           /&
         State_VGB(rho_,i,j,k,globalBLK))
  end subroutine calc_OLD_heating
  
end subroutine user_sources

!========================================================================
! SUBROUTINE GAMMA_IBLOCK
!========================================================================
!
! This subroutine allows the user to work with a vaiable polytropic index.
!
! The subroutine was written by ILR on Sep 25, 2003.
!
!========================================================================
subroutine Gamma_iBlock(i,j,k,iBLK,Gamma_BLK)
  use ModVarIndexes, ONLY: Bx_,By_,Bz_
  use ModAdvance,    ONLY: State_VGB,B0xCell_BLK,B0yCell_BLK,&
       B0zCell_BLK
  use ModGeometry,   ONLY: x_BLK,y_BLK,z_BLK
  use ModNumConst,   ONLY: cTolerance,cHalf,cOne,cTwo,cFour
  use ModUser,       ONLY: TemRatio,MaxB0,DegFrm1
  implicit none
  integer, intent(in):: i,j,k,iBLK
  real, intent(out):: Gamma_BLK
  real, parameter:: n0=cFour
  real:: BBr,TemMod
  real:: xx,yy,zz,RR
  real:: AAc,BBc,Fn1
  real:: Temp_BLK,DegF_BLK
  !\
  ! Get the coordinates and radial distance from the Sun::
  !/
  xx = x_BLK(i,j,k,iBLK)
  yy = y_BLK(i,j,k,iBLK)
  zz = z_BLK(i,j,k,iBLK)
  RR = sqrt(xx**2+yy**2+zz**2+cTolerance**2)
  if (RR.gt.cHalf) then
     !\
     ! Get the radial field component::
     !/
     BBr = (xx*(B0xCell_BLK(i,j,k,iBLK)+State_VGB(Bx_,i,j,k,iBLK)) +&
            yy*(B0yCell_BLK(i,j,k,iBLK)+State_VGB(By_,i,j,k,iBLK)) +&
            zz*(B0zCell_BLK(i,j,k,iBLK)+State_VGB(Bz_,i,j,k,iBLK)))/&
            RR
     TemMod = cOne+min(cOne,abs(BBr*RR**3/MaxB0))*(TemRatio-cOne)
  else
     TemMod = cOne
  endif
  Fn1 = DegFrm1*TemMod
  !\
  ! Find the solution for T as a function of R::
  !/
  AAc = cFour*(cOne+cHalf*n0)/Fn1
  BBc = cFour*(cOne+cHalf*n0+cHalf*cHalf*Fn1)/RR/Fn1
  Temp_BLK = cTwo*BBc/(sqrt(AAc**2+cFour*BBc)+AAc)
  DegF_BLK = n0+Fn1*Temp_BLK
  !\
  ! Get Gamma_BLK at distance R::
  !/
  Gamma_BLK = (DegF_BLK+cTwo)/DegF_BLK
end subroutine Gamma_iBlock

!========================================================================
!  SUBROUTINE USER_INITIAL_PERTURBATION
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
  use ModMain,      ONLY: nI,nJ,nK,nBLK,globalBLK,PROCtest,       &
       BLKtest,unusedBLK,UseInertial,UseUserHeating,UseUserB0,gcn,&
       x_,y_,z_
  use ModIO,        ONLY: restart
  use ModVarIndexes,ONLY:&
	!EnergyRL_ ,&  !^CFG UNCOMMENT IF ALWAVES
	rho_,rhoUx_,rhoUy_,rhoUz_,Bx_,By_,Bz_,P_
  use ModAdvance,   ONLY: State_VGB,B0xCell_BLK,B0yCell_BLK,      &
       B0zCell_BLK,tmp1_BLK,tmp2_BLK
  use ModProcMH,    ONLY: iProc,iComm
  use ModNumConst,  ONLY: cZero,cQuarter,cHalf,cOne,cTwo,cE1,cE9, &
       cTolerance,cThree
  use ModConst,     ONLY: Rsun,Msun,cGravitation
  use ModGeometry,  ONLY: x_BLK,y_BLK,z_BLK,R_BLK,cV_BLK,x2,y2,z2
  use ModPhysics,   ONLY: Gbody,g,inv_g,gm1,inv_gm1,ModulationP,  &
       ModulationRho,UseFluxRope,rot_period_dim,OMEGAbody,        &
       unitSI_U,unitSI_rho,unitSI_x,unitUSER_energydens,          &
       unitUSER_t,unitUSER_B,Body_T_dim,Body_rho_dim
  use ModUser,      ONLY: Mrope_GL98,InvH0,DoStaticICs,MaxB0,Tnot,&
       Bnot,DegFrm1,Emag_0,Ekin_0,Ethe_0
  use ModMpi
  implicit none
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
  real, external:: maxval_BLK
  real, dimension(3):: R_GL98_D,B_GL98_D
  !
  !---------------------------------------------------------------------------
  !\
  ! Variable meanings:
  !
  !
  !/
  !---------------------------------------------------------------------------
  !
  if (iProc==PROCtest.and.globalBLK==BLKtest) then
     call set_oktest('user_initial_perturbation',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if
  !\
  ! Find the maximum value of B0 at time zero::
  !/
  do iBLK=1,nBLK
     if (unusedBLK(iBLK)) CYCLE
     do k=1,nK; do j=1,nJ; do i=1,nI
        if (R_BLK(i,j,k,iBLK).ge.cOne) then
           tmp1_BLK(i,j,k,iBLK) = sqrt(    &
                B0xCell_BLK(i,j,k,iBLK)**2+&
                B0yCell_BLK(i,j,k,iBLK)**2+&
                B0zCell_BLK(i,j,k,iBLK)**2)
        else
           tmp1_BLK(i,j,k,iBLK) = cZero
        endif
     enddo; enddo; enddo
  end do
  !\
  ! Initialize some auxilary variables::
  !/
  Mrope_GL98 = cZero; volume = cZero
  Emag_0 = cZero; Ekin_0 = cZero; Ethe_0 = cZero
  !\
  ! Set the value of MaxB0::
  !/
  if (UseUserB0) then
     !     MaxB0 = 7.1873379E+00/1.526
     !     MaxB0 = 7.1873379E+00/2.0
     MaxB0 = Bnot/unitUSER_B
  else
     MaxB0 = maxval_BLK(2,tmp1_BLK)
  endif
  InvH0 = cGravitation*Msun/Rsun/unitSI_U**2
  call MPI_BCAST(MaxB0     ,1,MPI_REAL,0,iComm,iError)
  call MPI_BCAST(InvH0     ,1,MPI_REAL,0,iComm,iError)
  call MPI_BCAST(Mrope_GL98,1,MPI_REAL,0,iComm,iError)
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
           call get_atmosphere_BLK(i,j,k,iBLK,Dens_BLK,  &
                Pres_BLK,Gamma_BLK)
           State_VGB(rho_     ,i,j,k,iBLK) = Dens_BLK
           State_VGB(P_   ,i,j,k,iBLK) = Pres_BLK
           if (DoStaticICs) then
              State_VGB(rhoUx_   ,i,j,k,iBLK) = cZero
              State_VGB(rhoUy_   ,i,j,k,iBLK) = cZero
              State_VGB(rhoUz_   ,i,j,k,iBLK) = cZero
           else
              State_VGB(rhoUx_   ,i,j,k,iBLK) = Dens_BLK*&
                   4.0E+01*((ROne-cOne)/(Rmax-cOne))*xx/RR
              State_VGB(rhoUy_   ,i,j,k,iBLK) = Dens_BLK*&
                   4.0E+01*((ROne-cOne)/(Rmax-cOne))*yy/RR
              State_VGB(rhoUz_   ,i,j,k,iBLK) = Dens_BLK*&
                   4.0E+01*((ROne-cOne)/(Rmax-cOne))*zz/RR
           endif
           !State_VGB(EnergyRL_,i,j,k,iBLK) = Pres_BLK   *& !^CFG UNCOMMENT IF ALWAVES
           !     (cOne/(Gamma_BLK-cOne)-inv_gm1)            !^CFG UNCOMMENT IF ALWAVES
        endif
        if (UseFluxRope.and.restart) then
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
                ModulationRho*rho_GL98).lt.cOne*&
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
                ModulationP*p_GL98).lt.cOne* &
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
        !\
        ! Compute initial energies integrated over volume::
        !/
        volume = volume+cV_BLK(iBLK)
        Emag_0 = Emag_0+cV_BLK(iBLK)*cHalf*(&
             (B0xcell_BLK(i,j,k,iBLK)+State_VGB(Bx_,i,j,k,iBLK))**2+&
             (B0ycell_BLK(i,j,k,iBLK)+State_VGB(By_,i,j,k,iBLK))**2+&
             (B0zcell_BLK(i,j,k,iBLK)+State_VGB(Bz_,i,j,k,iBLK))**2)
        Ekin_0 = Ekin_0+cV_BLK(iBLK)*cHalf*(&
             (State_VGB(rhoUx_,i,j,k,iBLK)**2 +&
              State_VGB(rhoUy_,i,j,k,iBLK)**2 +&
              State_VGB(rhoUz_,i,j,k,iBLK)**2)/&
              State_VGB(rho_  ,i,j,k,iBLK))
        Ethe_0 = Ethe_0+cV_BLK(iBLK)*(&
             inv_gm1*State_VGB(P_,i,j,k,iBLK))
     end do; end do; end do
     !\
     ! Update the total energy::
     !/
     globalBLK=iBLK
     call correctE
  end do
  !\
  ! Compute initial energies averaged over volume::
  !/
  Emag_0 = Emag_0/volume
  Ekin_0 = Ekin_0/volume
  Ethe_0 = Ethe_0/volume
  !\
  ! Write out some statistics::
  !/
  call post_init_stat
  
end subroutine user_initial_perturbation

subroutine user_set_physics_constants
  use ModPhysics, ONLY: FaceState_VI,CellState_VI
end subroutine user_set_physics_constants

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
       ModulationP,cRot_x_GL98,cRot_y_GL98,cRot_z_GL98
  use ModUser,           ONLY: DoFirst_GL
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
  real, dimension(3,3):: Txz_GL98_DD,Txzy_GL98_DD
  if (iProc==0.and.DoFirst_GL) then
     write(*,*) ''
     write(*,*) '>>>>>>>>>>>>>>>>>>>                  <<<<<<<<<<<<<<<<<<<<<'
     write(*,*) '            Initial Perturbation Is Initiated!!!'
     write(*,*) '>>>>>>>>>>>>>>>>>>>                  <<<<<<<<<<<<<<<<<<<<<'
     write(*,*) ''
     write(*,*) 'cme_init called by processor',iProc
     write(*,*) 'B1_dim = ',cme_B1_dim
     write(*,*) 'cme_a  = ',cme_a
     write(*,*) 'cme_r1 = ',cme_r1
     write(*,*) 'cme_r0 = ',cme_r0
     write(*,*) 'cme_a1 = ',cme_a1
     write(*,*) 'ModulationRho = ',ModulationRho
     write(*,*) 'ModulationP   = ',ModulationP
     DoFirst_GL=.false.
  end if
  !
  delta = cOne/(cTwo*(cOne+cFour))
  !\
  ! Construct the rotational matrix Txzy_GL98_DD,
  ! and compute R1_GL98_D::
  !/
  Txz_GL98_DD  = matmul(rot_matrix_x(cRot_x_GL98*cDegToRad),&
                        rot_matrix_z(cRot_z_GL98*cDegToRad))
  Txzy_GL98_DD = matmul(rot_matrix_y(cRot_y_GL98*cDegToRad),&
                        Txz_GL98_DD)
  R1_GL98_D = matmul(Txzy_GL98_DD,R_GL98_D)
  !\
  ! CALCULATE CELL CENTER FOR GLOBAL CARTESIAN COORDINATES 
  ! (z,x,y)_BATSRUS -> (x,y,z)
  !/
  x = R1_GL98_D(z_)
  y = R1_GL98_D(x_)
  z = R1_GL98_D(y_)
  !\
  ! CALCULATE CELL CENTER FOR GLOBAL SPHERICAL COORDINATES 
  !/
  r = sqrt(x**2 + y**2 + z**2)
  cos_theta = z/r
  sin_theta = sqrt(x**2 + y**2)/r
  cos_phi   = x/sqrt(x**2 + y**2)
  sin_phi   = y/sqrt(x**2 + y**2)
  if (r.le.delta) then 
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
  if (r_2.le.delta) then
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
  if (sqrt(x_2**2 + y_2**2 + z_2**2).le.cme_r0) then
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
     ! Mind that (z,x,y)_BATSRUS -> (x,y,z)!!!
     ! X-COMPONENT OF MAGNETIC FIELD:: (x,y,z)_BATSRUS -> (x,y,z)
     !/
     B1_GL98_D(z_) = Br*sin_theta1*cos_phi1   + &
          Btheta*cos_theta1*cos_phi1    - &
          Bphi*sin_phi1
     !\
     ! Y-COMPONENT OF MAGNETIC FIELD:: (x,y,z)_BATSRUS -> (x,y,z)
     !/
     B1_GL98_D(x_) = Br*sin_theta1*sin_phi1   + &
          Btheta*cos_theta1*sin_phi1    + &
          Bphi*cos_phi1
     !\
     ! Z-COMPONENT OF MAGNETIC FIELD:: (x,y,z)_BATSRUS -> (x,y,z)
     !/
     B1_GL98_D(y_) = Br*cos_theta1-Btheta*sin_theta1
     !\
     ! Transform back to the original coordinates
     ! given by R_GL98_D:: 
     !/
     B_GL98_D  = matmul(transpose(Txzy_GL98_DD),B1_GL98_D)
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
  endif
  
end subroutine add_GL98_fluxrope

subroutine post_init_stat
  use ModVarIndexes, ONLY: rho_,P_
  use ModAdvance,    ONLY: State_VGB
  use ModProcMH,     ONLY: iProc
  use ModUser,       ONLY: InvH0,Mrope_GL98,MaxB0,Tnot,Bnot
  use ModPhysics,    ONLY: Gbody
  implicit none
  !\
  ! Post-initialization statistics::
  !/
  !  write(*,*) 'Mass in the flux rope on processor::',&
  !       iProc,Mrope_GL98
  if (iProc==0) then
     write(*,*) ''
     write(*,*) '>>>>>>>>>>>>>>>>>>> Pressure and Density Log <<<<<<<<<<<<<<<<<<<<<'
     write(*,*) ''
     write(*,*) 'The value of MaxB0 is :: ',MaxB0
     write(*,*) 'The value of Bnot  is :: ',Bnot
     write(*,*) 'The min,max P is      :: ',&
          minval(State_VGB(P_,:,:,:,:)),&
          maxval(State_VGB(P_,:,:,:,:))
     write(*,*) 'The min,max Rho is    :: ',&
          minval(State_VGB(rho_,:,:,:,:))  ,&
          maxval(State_VGB(rho_,:,:,:,:))
     write(*,*) 'The value of Tnot  is :: ',Tnot
     write(*,*) 'The value of InvH0 is :: ',InvH0
     write(*,*) 'The value of Gbody is :: ',Gbody
     write(*,*) ''
     write(*,*) '>>>>>>>>>>>>>>>>>>>                          <<<<<<<<<<<<<<<<<<<<<'
     write(*,*) ''
  end if
end subroutine post_init_stat

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
  if (RR.gt.cHalf) then
     Dens_BLK = Rho0/RR**2
     Pres_BLK = Rho0*inv_g/RR**2
  else
     Dens_BLK = Rho0/cHalf**2
     Pres_BLK = Rho0*inv_g/cHalf**2
  endif
  
end subroutine get_atmosphere_orig

subroutine get_atmosphere_poli(i,j,k,iBLK,Dens_BLK,Pres_BLK)
  !
  !---------------------------------------------------------------------------
  !
  ! This module computes the background atmosphere in the
  ! presence of gravity for a politropic equation of state::
  ! Use gamma=[1.05,1.10] only!!!
  ! The subroutine is written by ILR on May 15, 2003.
  !
  !---------------------------------------------------------------------------
  !
  use ModVarIndexes, ONLY: Bx_,By_,Bz_
  use ModAdvance,    ONLY: State_VGB,B0xCell_BLK,B0yCell_BLK,&
       B0zCell_BLK
  use ModGeometry,   ONLY: x_BLK,y_BLK,z_BLK
  use ModNumConst,   ONLY: cTolerance,cHalf,cOne
  use ModUser,       ONLY: InvH0,TemRatio,MaxB0
  use ModPhysics,    ONLY: inv_g
  implicit none
  
  integer, intent(in):: i,j,k,iBLK
  real, intent(out):: Dens_BLK,Pres_BLK
  real, parameter:: Rho0=cOne
  real, parameter:: Gamma_TRG=1.05E+00
  real:: xx,yy,zz,RR
  real:: RadDep,TemMod,BBr
  real:: g_mn1_tr,inv_g_tr
  !\
  ! Get the coordinates and radial distance from the Sun::
  !/
  g_mn1_tr = Gamma_TRG-cOne
  inv_g_tr = cOne/Gamma_TRG
  xx = x_BLK(i,j,k,iBLK)
  yy = y_BLK(i,j,k,iBLK)
  zz = z_BLK(i,j,k,iBLK)
  RR = sqrt(xx**2+yy**2+zz**2+cTolerance**2)
  if (RR.gt.cHalf) then
     !\
     ! Get the radial field component::
     !/
     BBr = (xx*(B0xCell_BLK(i,j,k,iBLK)+State_VGB(Bx_,i,j,k,iBLK)) +&
            yy*(B0yCell_BLK(i,j,k,iBLK)+State_VGB(By_,i,j,k,iBLK)) +&
            zz*(B0zCell_BLK(i,j,k,iBLK)+State_VGB(Bz_,i,j,k,iBLK)))/&
            RR
     TemMod = cOne+min(cOne,abs(BBr*RR**3/MaxB0)) *&
          (TemRatio-cOne)
     RadDep = (cOne+g_mn1_tr*inv_g_tr*InvH0*(cOne /&
          RR-cOne))**(cOne/g_mn1_tr)
     !\
     ! Get the mass density at the body::
     ! The Lat-Long dependence is determined by abs(Br)*R**3. For a
     ! simple dipole field this gives a cos(Latitude) dependence.
     !/
     !\
     ! simple cos(Latitude) dependence::
     !/
     Dens_BLK = Rho0*RadDep/TemMod
     Pres_BLK = inv_g*RadDep**Gamma_TRG
  else
     Dens_BLK = Rho0
     Pres_BLK = inv_g
  endif
end subroutine get_atmosphere_poli

subroutine get_atmosphere_BLK(i,j,k,iBLK,Dens_BLK,Pres_BLK,Gamma_BLK)
  !
  !---------------------------------------------------------------------------
  !
  ! This module computes the background atmosphere in the
  ! presence of gravity for a politropic equation of state
  ! with variable gamma = 1+n/2, where n=n0+n1*T^pow::
  ! The subroutine is written by ILR on May 29, 2003.
  !
  !---------------------------------------------------------------------------
  !
  use ModAdvance,    ONLY: State_VGB,B0xCell_BLK,B0yCell_BLK,&
       B0zCell_BLK
  use ModGeometry,   ONLY: x_BLK,y_BLK,z_BLK
  use ModNumConst,   ONLY: cTolerance,cHalf,cOne,cTwo,cFour
  use ModUser,       ONLY: TemRatio,MaxB0,DegFrm1
  use ModPhysics,    ONLY: g,inv_g
  implicit none
  integer, intent(in):: i,j,k,iBLK
  real, intent(out):: Dens_BLK,Pres_BLK,Gamma_BLK
  real, parameter:: n0=cFour
  real:: BBr,TemMod
  real:: xx,yy,zz,RR
  real:: AAc,BBc,Fn1,Fg1
  real:: Temp_BLK,DegF_BLK
  !\
  ! Get the coordinates and radial distance from the Sun::
  !/
  xx = x_BLK(i,j,k,iBLK)
  yy = y_BLK(i,j,k,iBLK)
  zz = z_BLK(i,j,k,iBLK)
  RR = sqrt(xx**2+yy**2+zz**2+cTolerance**2)
  if (RR.gt.cHalf) then
     !\
     ! Get the radial field component::
     !/
     BBr = (xx*B0xCell_BLK(i,j,k,iBLK) +&
            yy*B0yCell_BLK(i,j,k,iBLK) +&
            zz*B0zCell_BLK(i,j,k,iBLK))/&
            RR
     TemMod = cOne+min(cOne,abs(BBr*RR**3/MaxB0)) *&
          (TemRatio-cOne)
  else
     TemMod = cOne
  endif
  Fg1 = (cOne+cHalf*n0+cHalf*cHalf*DegFrm1*TemMod)/&
        (cOne+cHalf*n0+cHalf*cHalf*DegFrm1)
  Fn1 = DegFrm1*TemMod
  !\
  ! Find the solution for T as a function of R::
  !/
  AAc = cFour*(cOne+cHalf*n0)/Fn1
  BBc = cFour*(cOne+cHalf*n0+cHalf*cHalf*Fn1)/RR/Fn1
  Temp_BLK = cTwo*BBc/(sqrt(AAc**2+cFour*BBc)+AAc)
  DegF_BLK = n0+Fn1*Temp_BLK
  !\
  ! Get the mass density, pressure and gamma at distance R::
  !/
  Pres_BLK  = inv_g*(Temp_BLK**(cOne+cHalf*n0))   *&
       exp(-cHalf*Fn1*(cOne-Temp_BLK))/Fg1
  Dens_BLK  = g*Pres_BLK/Temp_BLK/TemMod
  Gamma_BLK = (DegF_BLK+cTwo)/DegF_BLK
end subroutine get_atmosphere_BLK

subroutine get_Br_BLK(xx,yy,zz,Br_BLK)
  use ModMain,     ONLY: Time_Simulation
  use ModNumConst, ONLY: cZero,cHalf,cOne,cTwo,cThree,  &
       cE1,cE9,cTolerance,cTiny
  use ModProcMH,   ONLY: iProc
  use ModUser,     ONLY: DoSecond,File_PFSSM,Head_PFSSM,&
       N_PFSSM,UNIT_PFSSM,R_PFSSM,Rs_PFSSM,Ro_PFSSM,    &
       unitPFSSM_B,iHead_PFSSM,g_nm,h_nm,factRatio1,    &
       MaxB0
  use ModPhysics,  ONLY: unitUSER_B,OMEGAbody,unitUSER_t
  implicit none
  
  real, intent(in):: xx,yy,zz
  real, intent(out):: Br_BLK
  integer:: iError
  integer:: i,n,m
  real:: gtemp,htemp
  real:: delta_m0,c_n
  real:: sinPhi,cosPhi
  real:: sinmPhi,cosmPhi
  real:: cosTheta,sinTheta
  real:: stuff1,stuff2,stuff3
  real:: sumr,sumt,sump,sumpsi
  real:: Rin_PFSSM,Theta_PFSSM,Phi_PFSSM
  real, dimension(N_PFSSM+1,N_PFSSM+1):: p_nm,dp_nm
  !\
  ! Calculate cell-centered spherical coordinates::
  !/
  Rin_PFSSM   = sqrt(xx**2+yy**2+zz**2)
  !\
  ! Avoid calculating B0 inside a critical radius = 0.5*Rsun
  !/
  if (Rin_PFSSM.lt.9.00E-01) then
     Br_BLK = cZero
     RETURN
  end if
  Theta_PFSSM = acos(zz/Rin_PFSSM)
  Phi_PFSSM   = atan2(yy/sqrt(xx**2+yy**2+cTiny**2),&
                      xx/sqrt(xx**2+yy**2+cTiny**2))
  sinTheta    = sqrt(xx**2+yy**2+cTiny**2)/Rin_PFSSM
  cosTheta    = zz/Rin_PFSSM
  sinPhi      = yy/sqrt(xx**2+yy**2+cTiny**2)
  cosPhi      = xx/sqrt(xx**2+yy**2+cTiny**2)
  Rin_PFSSM   = Ro_PFSSM
  R_PFSSM     = Rin_PFSSM     
  !\
  ! Update Phi_PFSSM for solar rotation::
  !/
  !  Phi_PFSSM = Phi_PFSSM-OMEGAbody*(Time_Simulation/unitUSER_t)
  !
  if (DoSecond) then
     !\
     ! Initialize once g(n+1,m+1) & h(n+1,m+1) by reading a file
     ! created from Web data::
     !/ 
     DoSecond=.false.
     if (iProc==0) then
        write(*,*) '!!!This is the call from ICs!!!'
        write(*,*) 'Norder = ',N_PFSSM
        write(*,*) 'Entered coefficient file name :: ',File_PFSSM
        write(*,*) 'Entered number of header lines:: ',iHead_PFSSM
     endif
     !\
     ! Formats adjusted for wso CR rad coeffs::
     !/
     open(UNIT_PFSSM,file=File_PFSSM,status='old',iostat=iError)
     if (iHead_PFSSM.ne.0) then
        do i=1,iHead_PFSSM
           read(UNIT_PFSSM,'(a)') Head_PFSSM
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
        read(UNIT_PFSSM,*,iostat=iError) n,m,gtemp,htemp
        if (iError.ne.0) EXIT
        if (n.gt.N_PFSSM.or.m.gt.N_PFSSM) CYCLE
        g_nm(n+1,m+1) = gtemp
        h_nm(n+1,m+1) = htemp
     enddo
     close(UNIT_PFSSM)
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
     ! Calculate the ratio sqrt(2m!)/(2^m*m!)::
     !/
     factRatio1(:) = cZero; factRatio1(1) = cOne
     do m=1,N_PFSSM
        factRatio1(m+1) = factRatio1(m)*&
             sqrt(cOne-cHalf/real(m))
     enddo
  endif
  !\
  ! Calculate polynomials with appropriate normalization
  ! for Theta_PFSSMa::
  !/
  do m=0,N_PFSSM
     if (m.eq.0) then
        delta_m0 = cOne
     else
        delta_m0 = cZero
     endif
     !\
     ! Eq.(27) from Altschuler et al. 1976::
     !/
     p_nm(m+1,m+1) = factRatio1(m+1)*sqrt((cTwo-delta_m0)    *&
          real(2*m+1))*sinTheta**m
     !\
     ! Eq.(28) from Altschuler et al. 1976::
     !/
     if (m.lt.N_PFSSM) &
          p_nm(m+2,m+1)  = p_nm(m+1,m+1)*sqrt(real(2*m+3))   *&
          cosTheta
     !\
     ! Eq.(30) from Altschuler et al. 1976::
     !/
     dp_nm(m+1,m+1)      = factRatio1(m+1)*sqrt((cTwo        -&
          delta_m0)*real(2*m+1))*m*cosTheta*sinTheta**(m-1)
     !\
     ! Eq.(31) from Altschuler et al. 1976::
     !/
     if (m.lt.N_PFSSM) &
          dp_nm(m+2,m+1) = sqrt(real(2*m+3))*(cosTheta       *&
          dp_nm(m+1,m+1)-sinTheta*p_nm(m+1,m+1))
  enddo
  do m=0,N_PFSSM-2; do n=m+2,N_PFSSM
     !\
     ! Eq.(29) from Altschuler et al. 1976::
     !/
     stuff1         = sqrt(real(2*n+1)/real(n**2-m**2))
     stuff2         = sqrt(real(2*n-1))
     stuff3         = sqrt(real((n-1)**2-m**2)/real(2*n-3))
     p_nm(n+1,m+1)  = stuff1*(stuff2*cosTheta*p_nm(n,m+1)    -&
          stuff3*p_nm(n-1,m+1))
     !\
     ! Eq.(32) from Altschuler et al. 1976::
     !/
     dp_nm(n+1,m+1) = stuff1*(stuff2*(cosTheta*dp_nm(n,m+1)  -&
          sinTheta*p_nm(n,m+1))-stuff3*dp_nm(n-1,m+1))
  enddo; enddo
  !\
  ! Apply Schmidt normailization::
  !/
  do m=0,N_PFSSM; do n=m,N_PFSSM
     !\
     ! Eq.(33) from Altschuler et al. 1976::
     !/
     stuff1 = cOne/sqrt(real(2*n+1))
     !\
     ! Eq.(34) from Altschuler et al. 1976::
     !/
     p_nm(n+1,m+1)  = p_nm(n+1,m+1)*stuff1
     dp_nm(n+1,m+1) = dp_nm(n+1,m+1)*stuff1
  enddo; enddo
  !\
  ! Truncate the value of sinTheta::
  !/
  if (sinTheta.eq.cZero) sinTheta = cOne/(cE9*cE1)
  !\
  ! Initialize the values of sumr,sumt,sump, and sumpsi::
  !/
  sumr = cZero; sumt   = cZero
  sump = cZero; sumpsi = cZero
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
        c_n    = -(Ro_PFSSM/Rs_PFSSM)**(n+2)
        !\
        ! Br_PFSSM = -d(Psi_PFSSM)/dR_PFSSM::
        !/
        stuff1 = (real(n)+cOne)*(Ro_PFSSM/R_PFSSM)**(n+2)    -&
             c_n*real(n)*(R_PFSSM/Rs_PFSSM)**(n-1)
        stuff2 = g_nm(n+1,m+1)*cosmPhi+h_nm(n+1,m+1)*sinmPhi
        sumr   = sumr+p_nm(n+1,m+1)*stuff1*stuff2
        !\
        ! Bt_PFSSM = -(1/R_PFSSM)*d(Psi_PFSSM)/dTheta_PFSSM::
        !/
        stuff1 = (Ro_PFSSM/R_PFSSM)**(n+2)+c_n*(R_PFSSM      /&
             Rs_PFSSM)**(n-1)
        sumt   = sumt-dp_nm(n+1,m+1)*stuff1*stuff2
        !\
        ! Psi_PFSSM::
        !/
        sumpsi = sumpsi+R_PFSSM*p_nm(n+1,m+1)*stuff1*stuff2
        !\
        ! Bp_PFSSM = -(1/R_PFSSM)*d(Psi_PFSSM)/dPhi_PFSSM::
        !/
        stuff2 = g_nm(n+1,m+1)*sinmPhi-h_nm(n+1,m+1)*cosmPhi
        sump   = sump+p_nm(n+1,m+1)*real(m)/sinTheta*stuff1  *&
             stuff2
     enddo
  enddo
  !\
  ! Compute Br_BLK at R_PFSSM = Ro_PFSSM = cOne::
  ! Apply field strength normalization::
  ! Note that unitUSER_B is in Gauss,
  ! while unit_PFSSM_B is in microT=0.01*Gauss::
  !/
  !  Br_BLK = sqrt(sumr**2+sumt**2+sump**2)*unitPFSSM_B/unitUSER_B
  Br_BLK = abs(sumr)*unitPFSSM_B/unitUSER_B
end subroutine get_Br_BLK
!============================================================================
subroutine get_user_b0(xx,yy,zz,B0_PFSSM)
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
  !
  use ModMain,     ONLY: Time_Simulation
  use ModNumConst, ONLY: cZero,cHalf,cOne,cTwo,cThree, &
       cE1,cE9,cTolerance,cTiny,cPi
  use ModProcMH,   ONLY: iProc
  use ModUser,     ONLY: DoFirst,File_PFSSM,Head_PFSSM,&
       N_PFSSM,UNIT_PFSSM,R_PFSSM,Rs_PFSSM,Ro_PFSSM,   &
       unitPFSSM_B,iHead_PFSSM,g_nm,h_nm,factRatio1
  use ModPhysics,  ONLY: unitUSER_B,OMEGAbody,unitUSER_t
  implicit none
  
  real, intent(in):: xx,yy,zz
  real, intent(out), dimension(3):: B0_PFSSM
  integer:: iError
  integer:: i,n,m
!  real, parameter:: Phi_Shift=-cPi*1.39E+02/1.80E+02
  real, parameter:: Phi_Shift=cPi
  real:: gtemp,htemp
  real:: delta_m0,c_n
  real:: sinPhi,cosPhi
  real:: sinmPhi,cosmPhi
  real:: cosTheta,sinTheta
  real:: stuff1,stuff2,stuff3
  real:: sumr,sumt,sump,sumpsi
  real:: Rin_PFSSM,Theta_PFSSM,Phi_PFSSM
  real:: Br_PFSSM,Btheta_PFSSM,Bphi_PFSSM,Psi_PFSSM
  real, dimension(N_PFSSM+1,N_PFSSM+1):: p_nm,dp_nm
  !\
  ! Calculate cell-centered spherical coordinates::
  !/
  Rin_PFSSM   = sqrt(xx**2+yy**2+zz**2)
  !\
  ! Avoid calculating B0 inside a critical radius = 0.5*Rsun
  !/
  if (Rin_PFSSM.lt.9.00E-01) then
     B0_PFSSM = cZero
     RETURN
  end if
  Theta_PFSSM = acos(zz/Rin_PFSSM)
  Phi_PFSSM   = atan2(yy/sqrt(xx**2+yy**2+cTiny**2),&
                      xx/sqrt(xx**2+yy**2+cTiny**2))
  sinTheta    = sqrt(xx**2+yy**2+cTiny**2)/Rin_PFSSM
  cosTheta    = zz/Rin_PFSSM
  sinPhi      = yy/sqrt(xx**2+yy**2+cTiny**2)
  cosPhi      = xx/sqrt(xx**2+yy**2+cTiny**2)
  !\
  ! Update Phi_PFSSM for central meridian and solar rotation::
  !/
  Phi_PFSSM = Phi_PFSSM-OMEGAbody*(Time_Simulation/unitUSER_t)
  Phi_PFSSM = Phi_PFSSM-Phi_Shift
  !\
  ! Set the source surface radius::
  !/
  if (Rin_PFSSM.gt.Rs_PFSSM) then 
     R_PFSSM = Rs_PFSSM
  else
     R_PFSSM = Rin_PFSSM     
  endif
  if (DoFirst) then
     !\
     ! Initialize once g(n+1,m+1) & h(n+1,m+1) by reading a file
     ! created from Web data::
     !/ 
     DoFirst=.false.
     if (iProc==0) then
        write(*,*) 'Norder = ',N_PFSSM
        write(*,*) 'Entered coefficient file name :: ',File_PFSSM
        write(*,*) 'Entered number of header lines:: ',iHead_PFSSM
     endif
     !\
     ! Formats adjusted for wso CR rad coeffs::
     !/
     open(UNIT_PFSSM,file=File_PFSSM,status='old',iostat=iError)
     if (iHead_PFSSM.ne.0) then
        do i=1,iHead_PFSSM
           read(UNIT_PFSSM,'(a)') Head_PFSSM
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
        read(UNIT_PFSSM,*,iostat=iError) n,m,gtemp,htemp
        if (iError.ne.0) EXIT
        if (n.gt.N_PFSSM.or.m.gt.N_PFSSM) CYCLE
        g_nm(n+1,m+1) = gtemp
        h_nm(n+1,m+1) = htemp
     enddo
     close(UNIT_PFSSM)
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
     ! Calculate the ratio sqrt(2m!)/(2^m*m!)::
     !/
     factRatio1(:) = cZero; factRatio1(1) = cOne
     do m=1,N_PFSSM
        factRatio1(m+1) = factRatio1(m)*&
             sqrt(cOne-cHalf/real(m))
     enddo
  endif
  !\
  ! Calculate polynomials with appropriate normalization
  ! for Theta_PFSSMa::
  !/
  do m=0,N_PFSSM
     if (m.eq.0) then
        delta_m0 = cOne
     else
        delta_m0 = cZero
     endif
     !\
     ! Eq.(27) from Altschuler et al. 1976::
     !/
     p_nm(m+1,m+1) = factRatio1(m+1)*sqrt((cTwo-delta_m0)    *&
          real(2*m+1))*sinTheta**m
     !\
     ! Eq.(28) from Altschuler et al. 1976::
     !/
     if (m.lt.N_PFSSM) &
          p_nm(m+2,m+1)  = p_nm(m+1,m+1)*sqrt(real(2*m+3))   *&
          cosTheta
     !\
     ! Eq.(30) from Altschuler et al. 1976::
     !/
     dp_nm(m+1,m+1)      = factRatio1(m+1)*sqrt((cTwo        -&
          delta_m0)*real(2*m+1))*m*cosTheta*sinTheta**(m-1)
     !\
     ! Eq.(31) from Altschuler et al. 1976::
     !/
     if (m.lt.N_PFSSM) &
          dp_nm(m+2,m+1) = sqrt(real(2*m+3))*(cosTheta       *&
          dp_nm(m+1,m+1)-sinTheta*p_nm(m+1,m+1))
  enddo
  do m=0,N_PFSSM-2; do n=m+2,N_PFSSM
     !\
     ! Eq.(29) from Altschuler et al. 1976::
     !/
     stuff1         = sqrt(real(2*n+1)/real(n**2-m**2))
     stuff2         = sqrt(real(2*n-1))
     stuff3         = sqrt(real((n-1)**2-m**2)/real(2*n-3))
     p_nm(n+1,m+1)  = stuff1*(stuff2*cosTheta*p_nm(n,m+1)    -&
          stuff3*p_nm(n-1,m+1))
     !\
     ! Eq.(32) from Altschuler et al. 1976::
     !/
     dp_nm(n+1,m+1) = stuff1*(stuff2*(cosTheta*dp_nm(n,m+1)  -&
          sinTheta*p_nm(n,m+1))-stuff3*dp_nm(n-1,m+1))
  enddo; enddo
  !\
  ! Apply Schmidt normailization::
  !/
  do m=0,N_PFSSM; do n=m,N_PFSSM
     !\
     ! Eq.(33) from Altschuler et al. 1976::
     !/
     stuff1 = cOne/sqrt(real(2*n+1))
     !\
     ! Eq.(34) from Altschuler et al. 1976::
     !/
     p_nm(n+1,m+1)  = p_nm(n+1,m+1)*stuff1
     dp_nm(n+1,m+1) = dp_nm(n+1,m+1)*stuff1
  enddo; enddo
  !\
  ! Truncate the value of sinTheta::
  !/
  if (sinTheta.eq.cZero) sinTheta = cOne/(cE9*cE1)
  !\
  ! Initialize the values of sumr,sumt,sump, and sumpsi::
  !/
  sumr = cZero; sumt   = cZero
  sump = cZero; sumpsi = cZero
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
        c_n    = -(Ro_PFSSM/Rs_PFSSM)**(n+2)
        !\
        ! Br_PFSSM = -d(Psi_PFSSM)/dR_PFSSM::
        !/
        stuff1 = (real(n)+cOne)*(Ro_PFSSM/R_PFSSM)**(n+2)    -&
             c_n*real(n)*(R_PFSSM/Rs_PFSSM)**(n-1)
        stuff2 = g_nm(n+1,m+1)*cosmPhi+h_nm(n+1,m+1)*sinmPhi
        sumr   = sumr+p_nm(n+1,m+1)*stuff1*stuff2
        !\
        ! Bt_PFSSM = -(1/R_PFSSM)*d(Psi_PFSSM)/dTheta_PFSSM::
        !/
        stuff1 = (Ro_PFSSM/R_PFSSM)**(n+2)+c_n*(R_PFSSM      /&
             Rs_PFSSM)**(n-1)
        sumt   = sumt-dp_nm(n+1,m+1)*stuff1*stuff2
        !\
        ! Psi_PFSSM::
        !/
        sumpsi = sumpsi+R_PFSSM*p_nm(n+1,m+1)*stuff1*stuff2
        !\
        ! Bp_PFSSM = -(1/R_PFSSM)*d(Psi_PFSSM)/dPhi_PFSSM::
        !/
        stuff2 = g_nm(n+1,m+1)*sinmPhi-h_nm(n+1,m+1)*cosmPhi
        sump   = sump+p_nm(n+1,m+1)*real(m)/sinTheta*stuff1  *&
             stuff2
     enddo
  enddo
  !\
  ! Compute (Br_PFSSM,Btheta_PFSSM,Bphi_PFSSM) and Psi_PFSSM::
  !/
  Psi_PFSSM    = sumpsi
  Br_PFSSM     = sumr
  if (Rin_PFSSM.gt.Rs_PFSSM) Br_PFSSM = Br_PFSSM              *&
       (Rs_PFSSM/Rin_PFSSM)**2
  Btheta_PFSSM = sumt
  Bphi_PFSSM   = sump
  !\
  ! Magnetic field components in global Cartesian coordinates::
  ! Set B0xCell_BLK::
  !/
  B0_PFSSM(1) = Br_PFSSM*sinTheta*cosPhi+&
       Btheta_PFSSM*cosTheta*cosPhi     -&
       Bphi_PFSSM*sinPhi
  !\
  ! Set B0yCell_BLK::
  !/
  B0_PFSSM(2) = Br_PFSSM*sinTheta*sinPhi+&
       Btheta_PFSSM*cosTheta*sinPhi     +&
       Bphi_PFSSM*cosPhi
  !\
  ! Set B0zCell_BLK::
  !/
  B0_PFSSM(3) = Br_PFSSM*cosTheta       -&
       Btheta_PFSSM*sinTheta
  !\
  ! Apply field strength normalization::
  ! Note that unitUSER_B is in Gauss,
  ! while unit_PFSSM_B is in microT=0.01*Gauss::
  !/
  B0_PFSSM(1) = 3.0*B0_PFSSM(1)*unitPFSSM_B/unitUSER_B
  B0_PFSSM(2) = 3.0*B0_PFSSM(2)*unitPFSSM_B/unitUSER_B
  B0_PFSSM(3) = 3.0*B0_PFSSM(3)*unitPFSSM_B/unitUSER_B
end subroutine get_user_b0
!========================================================================
!========================================================================
!  SUBROUTINE USER_SET_OUTER_BCS 
!========================================================================
! This subroutine allows the user to apply boundary conditions to the outer
! body which are problem specific and cannot be created using the predefined
! options in BATSRUS.
! The variables specific to the problem are loaded from ModUser
! Any of the outer boundary specified in BATSRUS can be used here.
  subroutine user_set_outerBCs(iBlock,TypeBc,IsFound)
  use ModMain
  use ModPhysics
  use ModProcMH
  use ModAdvance 
  use ModGeometry
  use ModSetOuterBC
  use ModUser
  implicit none
  integer,intent(in)::iBlock
  ! Variables required by this user subroutine
  logical,intent(out) :: IsFound
  character (len=20),intent(in) :: TypeBc
end subroutine user_set_outerBCs
!========================================================================
!========================================================================
!  SUBROUTINE USER_SPECIFY_INITIAL_REFINEMENT
!========================================================================
!
!\
! This subroutine allows the user to add an initial refinement type
! based on a geometric criteria.  The `case' specified in the PARAM.in file
! will be read here.
! As with all user subroutines, the variables declared in ModUser are
! available here.  Again, as with other user subroutines DO NOT MODIFY ANY
! GLOBAL VARIABLE DEFINED IN THE MODULES INCLUDED IN THIS SUBROUTINE UNLESS
! SPECIFIED!!
!
!/
! 
subroutine user_specify_initial_refinement(iBLK,refineBlock,lev,dxBlock,   &
              xCenter,yCenter,zCenter,rCenter,                        &
              minx,miny,minz,minR,maxx,maxy,maxz,maxR,IsFound)
  use ModPhysics
  use ModNumConst
  use ModGeometry
  use ModAMR
  use ModUser 
  implicit none
  logical,intent(out) :: refineBlock, IsFound
  integer, intent(in) :: lev
  real, intent(in) :: dxBlock
  real, intent(in) :: xCenter,yCenter,zCenter,rCenter
  real, intent(in) :: minx,miny,minz,minR
  real, intent(in) :: maxx,maxy,maxz,maxR
  integer,intent(in) :: iBLK
end subroutine user_specify_initial_refinement
!========================================================================
!========================================================================
!  SUBROUTINE USER_AMR_CRITERIA
!========================================================================
!
!\
! This subroutine allows the user to add a refinement type
! based on a geometric criteria or physical criteria.  The `case'
! specified in the #AMRCRITERIA file will be read here.
!
!/
subroutine user_amr_criteria(iBLK, userCriteria, TypeCriteria, IsFound)
  use ModMain
  use ModAdvance
  use ModGeometry, ONLY:x_BLK,y_BLK,z_BLK,R_BLK,dx_BLK,dy_BLK,dz_BLK,true_cell
  use ModPhysics
  use ModConst
  use ModUser
 implicit none
 ! Variables required by this user subroutine
  logical ,intent(out):: IsFound
  integer, intent(in) :: iBLK
  integer, intent(out) :: userCriteria
  character (len=20),intent(in) :: TypeCriteria
end subroutine user_amr_criteria
!========================================================================
!========================================================================
!  SUBROUTINE USER_SET_ICs
! (It will include set_ICs_global.f90
!!\
! Calculates the initial conditions of the grid for the Global Heliosphere
!
! Written by Merav Opher Feb 14  2002
!/
! OMEGAbody is the rotation frequency of the Sun
!========================================================================

! This subroutine allows the user to apply initial conditions to the domain
! which are problem specific and cannot be created using the predefined
! options in BATSRUS.
! The variables specific to the problem are loaded from ModUser
subroutine user_set_ICs
  use ModMain, ONLY : nI,nJ,nK,gcn,globalBLK
  use ModAdvance
  use ModGeometry, ONLY : x2,y2,z2,x_BLK,y_BLK,z_BLK,R_BLK,true_cell
  use ModIO, ONLY : restart
  use ModPhysics
  use ModUser
  implicit none
end subroutine user_set_ICs
!========================================================================
!========================================================================
subroutine user_update_states(iStage,iBlock)
  use ModVarIndexes
  use ModSize
  use ModAdvance, ONLY: State_VGB
  use ModMain,    ONLY: nStage
  use ModPhysics, ONLY: inv_gm1
  implicit none
  integer,intent(in):: iStage,iBlock
  integer::i,j,k
  real:: Gamma_BLK
  
  call update_states_MHD(iStage,iBlock)
  !\
  ! Begin update of pressure and relaxation energy::
  !/
  !  if (iStage/=nStage) return
  !  do k=1-gcn,nK+gcn; do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
  do k=1,nK; do j=1,nJ; do i=1,nI
     call Gamma_iBlock(i,j,k,iBlock,Gamma_BLK)
     State_VGB(P_   ,i,j,k,iBlock)=           &
          (Gamma_BLK-cOne)                       *&
          (inv_gm1*State_VGB(P_,i,j,k,iBlock)&
          !+State_VGB(EnergyRL_,i,j,k,iBlock)&  !^CFG UNCOMMENT IF ALWAVES
	  )
    ! State_VGB(EnergyRL_,i,j,k,iBlock)=           &!^CFG UNCOMMENT IF ALWAVES
    !      State_VGB(P_,i,j,k,iBlock)*(cOne   /&    !^CFG UNCOMMENT IF ALWAVES
    !      (Gamma_BLK-cOne)-inv_gm1)                !^CFG UNCOMMENT IF ALWAVES 
  end do; end do; end do
  call correctE
  !\
  ! End update of pressure and relaxation energy::
  !/
end subroutine user_update_states
!========================================================================
!========================================================================
!  SUBROUTINE USER_WRITE_PROGRESS
!========================================================================

!
! This subroutine allows the user to write to the log files variables
! (normalized and not) which are problem specific.
!
! The variables specific to the problem are loaded from ModUser
subroutine user_write_progress
    use ModProcMH
   use ModMain
   use ModPhysics
   use ModUser
   implicit none
end subroutine user_write_progress
!========================================================================
!========================================================================
!  SUBROUTINE USER_GET_LOG_VAR
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
  use ModIO,         ONLY: dn_output,logfile_
  use ModMain,       ONLY: unusedBLK,nBLK,iteration_number
  use ModVarIndexes, ONLY: Bx_,By_,Bz_,rho_,rhoUx_,rhoUy_,   &
       rhoUz_,P_
  use ModAdvance,    ONLY: tmp1_BLK,B0xCell_BLK,B0yCell_BLK, &
       B0zCell_BLK,State_VGB
  use ModPhysics,    ONLY: inv_gm1,unitSI_x,unitSI_energydens 
  use ModNumConst,   ONLY: cOne,cHalf,cE1,cE6 
  use ModUser,       ONLY: Emag_0,Ekin_0,Ethe_0
  implicit none
  real, intent(out):: VarValue
  character (LEN=10), intent(in):: TypeVar 
  !
  integer:: iBLK
  real:: volume
  real, external:: integrate_BLK
  !
  tmp1_BLK = cOne; volume = integrate_BLK(nProc,tmp1_BLK)
  !\
  ! Compute Emag_0,Ekin_0,Ethe_0 one time at
  ! iteration number==dn_output(logfile_)::
  !/
  if (iteration_number==dn_output(logfile_)) then
     do iBLK=1,nBLK
        if (unusedBLK(iBLK)) cycle
        tmp1_BLK(:,:,:,iBLK) = &
             (B0xcell_BLK(:,:,:,iBLK)+State_VGB(Bx_,:,:,:,iBLK))**2+&
             (B0ycell_BLK(:,:,:,iBLK)+State_VGB(By_,:,:,:,iBLK))**2+&
             (B0zcell_BLK(:,:,:,iBLK)+State_VGB(Bz_,:,:,:,iBLK))**2
     end do
     Emag_0 = cHalf*integrate_BLK(nProc,tmp1_BLK)/volume
     do iBLK=1,nBLK
        if (unusedBLK(iBLK)) cycle
        tmp1_BLK(:,:,:,iBLK) = &
             (State_VGB(rhoUx_,:,:,:,iBLK)**2 +&
              State_VGB(rhoUy_,:,:,:,iBLK)**2 +&
              State_VGB(rhoUz_,:,:,:,iBLK)**2)/&
              State_VGB(rho_  ,:,:,:,iBLK)             
     end do
     Ekin_0 = cHalf*integrate_BLK(nProc,tmp1_BLK)/volume
     do iBLK=1,nBLK
        if (unusedBLK(iBLK)) cycle
        tmp1_BLK(:,:,:,iBLK) = &
             inv_gm1*State_VGB(P_,:,:,:,iBLK)
     end do
     Ethe_0 = integrate_BLK(nProc,tmp1_BLK)/volume
  endif
  !\
  ! Define log variable to be saved::
  !/
  select case(TypeVar)
  case('em_r','Em_r')
     do iBLK=1,nBLK
        if (unusedBLK(iBLK)) cycle
        tmp1_BLK(:,:,:,iBLK) = &
             (B0xcell_BLK(:,:,:,iBLK)+State_VGB(Bx_,:,:,:,iBLK))**2+&
             (B0ycell_BLK(:,:,:,iBLK)+State_VGB(By_,:,:,:,iBLK))**2+&
             (B0zcell_BLK(:,:,:,iBLK)+State_VGB(Bz_,:,:,:,iBLK))**2
     end do
     VarValue = cHalf*integrate_BLK(nProc,tmp1_BLK) 
     VarValue = VarValue*cE1*cE6*unitSI_energydens*unitSI_x**3
  case('ek_r','Ek_r')
     do iBLK=1,nBLK
        if (unusedBLK(iBLK)) cycle
        tmp1_BLK(:,:,:,iBLK) = &
             (State_VGB(rhoUx_,:,:,:,iBLK)**2 +&
              State_VGB(rhoUy_,:,:,:,iBLK)**2 +&
              State_VGB(rhoUz_,:,:,:,iBLK)**2)/&
              State_VGB(rho_  ,:,:,:,iBLK)             
     end do
     VarValue = cHalf*integrate_BLK(nProc,tmp1_BLK)
     VarValue = VarValue*cE1*cE6*unitSI_energydens*unitSI_x**3
  case('et_r','Et_r')
     do iBLK=1,nBLK
        if (unusedBLK(iBLK)) cycle
        tmp1_BLK(:,:,:,iBLK) = &
             inv_gm1*State_VGB(P_,:,:,:,iBLK)
     end do
     VarValue = integrate_BLK(nProc,tmp1_BLK)
     VarValue = VarValue*cE1*cE6*unitSI_energydens*unitSI_x**3
  case default
     VarValue = -7777.
     write(*,*) 'Warning in set_user_logvar: unknown logvarname = ',TypeVar
  end select
end subroutine user_get_log_var
!========================================================================
