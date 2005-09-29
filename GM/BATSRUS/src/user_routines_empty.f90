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

module ModUser
  use ModSize, ONLY: nI,nJ,nK
! Define variables global to this file only here
!
! Example:
!
!  real :: test_real
!  integer :: test_integer
!  character (len=100) :: test_character
  real:: MaxB0,Rs_PFSSM
  real,dimension(1:nI,1:nJ,1:nK)::&
       Srho,SrhoUx,SrhoUy,SrhoUz,SBx,SBy,SBz,SE,Sp
end module
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
!
!/

subroutine user_read_inputs
  use ModInterface
  use ModProcMH
  use ModReadParam
  use ModUser
  use ModMain
  use ModIO, ONLY: write_prefix, write_myname, iUnitOut

  implicit none
  integer::i
  character (len=100) :: NameCommand
  !---------------------------------------------------------------------------

  if(iProc==0.and.lVerbose > 0)then
     call write_prefix; write(iUnitOut,*)'User read_input starts'
  endif
  do
     if(.not.read_line() ) EXIT
     if(.not.read_command(NameCommand)) CYCLE

     ! It is sufficient to leave only the read commands for those
     ! parameters you want to switch on or modify. The default value for all 
     !UseUserLogical is .false.
     select case(NameCommand)
!     case("#USER_FLAGS")
!        call read_var('UseUserInnerBCs'         ,UseUserInnerBCs)
!        call read_var('UseUserSource'           ,UseUserSource)
!        call read_var('UseUserPerturbation'     ,UseUserPerturbation)
!        call read_var('UseUserOuterBcs'         ,UseUserOuterBcs)
!        call read_var('UseUserICs'              ,UseUserICs)
!        call read_var('UseUserSpecifyRefinement',UseUserSpecifyRefinement)
!        call read_var('UseUserLogFiles'         ,UseUserLogFiles)
!        call read_var('UseUserWritePlot'        ,UseUserWritePlot)
!        call read_var('UseUserAMR'              ,UseUserAMR)
!        call read_var('UseUserEchoInput'        ,UseUserEchoInput)
!        call read_var('UseUserB0'               ,UseUserB0)
     case("#PFSSM")
        call read_var('UseUserB0' ,UseUserB0)
!        call read_var('Ro_PFSSM'  ,Ro_PFSSM)
!        call read_var('Rs_PFSSM'  ,Rs_PFSSM)
!        call read_var('File_PFSSM',File_PFSSM)
!        call read_var('Head_PFSSM',Head_PFSSM)
     case('#USERINPUTEND')
        if(iProc==0.and.lVerbose > 0)then
           call write_prefix; write(iUnitOut,*)'User read_input ends'
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
  implicit none
  call heat_source_averages
! Define your own heating function here
end subroutine user_heat_source

subroutine user_set_boundary_cells(iBLK)
  use ModGeometry	
  implicit none
  integer,intent(in)::iBLK
!  SHOULD define IsBoundaryCell_GI(:,:,:,ExtraBc_) using
!  a boundary condition for iBLK block
!  EXAMPLE: OUTER SPHERICAL BOUNDARY of radius of 100.
!  IsBoundaryCell_GI(:,:,:,ExtraBc_) = R_BLK(:,:,:,iBLK)<100.
end subroutine user_set_boundary_cells

!========================================================================
!========================================================================
!! The following routines are used to promote backward compatability of
!! user routines.  They should probably be located somewhere else, but
!! for now will be left here.
!========================================================================
!========================================================================
!========================================================================
!  SUBROUTINE USER_FACE_BCS
!========================================================================
subroutine user_face_bcs(iFace,jFace,kFace,iBlock,iSide,iBoundary,&
     iter,time_now, FaceCoords_D,&
     VarsTrueFace_V,VarsGhostFace_V,&
     B0Face_D,  UseIonosphereHere,UseRotatingBcHere)
  use ModUser
  use ModMain
  use ModAdvance
  use ModPhysics	
  implicit none
  
  integer,intent(in)::iFace,jFace,kFace,iBlock,iSide,iBoundary,iter
  real,intent(in)::time_now
  real,dimension(nDim),intent(in)::FaceCoords_D,B0Face_D
  real,dimension(nFaceValueVars),intent(in)::VarsTrueFace_V
  real,dimension(nFaceValueVars),intent(out)::VarsGhostFace_V
  logical,intent(in)::UseIonosphereHere,UseRotatingBcHere
  
  
  if( UseUserInnerBcs) then
!For backward compatibility only! Remove this, if you need not an old version
     call user_set_innerBCs(iter,time_now, &
          FaceCoords_D(x_), FaceCoords_D(y_), FaceCoords_D(z_),&
          !
          VarsTrueFace_V(rho_),&
          VarsTrueFace_V(Ux_ ),&
          VarsTrueFace_V(Uy_ ),&
          VarsTrueFace_V(Uz_ ),&
          VarsTrueFace_V(Bx_ ),&
          VarsTrueFace_V(By_ ),&
          VarsTrueFace_V(Bz_ ),&
          VarsTrueFace_V(P_),&
          !
          VarsGhostFace_V(rho_),& 
          VarsGhostFace_V(Ux_ ),&
          VarsGhostFace_V(Uy_ ),&
          VarsGhostFace_V(Uz_ ),&  
          VarsGhostFace_V(Bx_ ),&
          VarsGhostFace_V(By_ ),&
          VarsGhostFace_V(Bz_ ),&  
          VarsGhostFace_V(P_),&
          !
          B0Face_D(x_),B0Face_D(y_),B0Face_D(z_),&
          UseIonosphereHere,UseRotatingBcHere,&
          FaceState_VI(rho_,iBoundary), FaceState_VI(P_,iBoundary),&
          iBoundary)
  end if
end subroutine user_face_bcs




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

subroutine user_set_innerBCs(iter,time_now, XFace, YFace, ZFace,&
            RhoFaceOutside,VxFaceOutside,VyFaceOutside,VzFaceOutside,&
            BxFaceOutside,ByFaceOutside,BzFaceOutside,pFaceOutside ,&
            RhoFaceInside,VxFaceInside, VyFaceInside, VzFaceInside,&
            BxFaceInside,ByFaceInside,BzFaceInside, pFaceInside,&
            B0xFace,B0yFace,B0zFace,&
            UseIonosphereHere,UseRotatingBcHere,&
            RhoBodyHere, pBodyHere,BodyNumberHere)

  use ModMain
  use ModPhysics
  use ModUser
  use ModProcMH
  use ModGeometry

  implicit none

  ! Variables required by this user subroutine

  ! Arguments
  integer, intent(in) :: iter,BodyNumberHere
  real, intent(in) :: time_now, XFace, YFace, ZFace,&
            RhoFaceOutside,VxFaceOutside,VyFaceOutside,VzFaceOutside,&
            BxFaceOutside,ByFaceOutside,BzFaceOutside,pFaceOutside
  real, intent(out) ::  RhoFaceInside,VxFaceInside, VyFaceInside, VzFaceInside,&
            BxFaceInside,ByFaceInside,BzFaceInside, pFaceInside 
  real, intent (in) ::  B0xFace,B0yFace,B0zFace
  real, intent (in) ::  RhoBodyHere,pBodyHere
  logical,intent (in) :: UseIonosphereHere
  logical,intent (in) :: UseRotatingBcHere
  logical :: oktest,oktest_me
  real:: xBodyHere,yBodyHere,zBodyHere

  if(BodyNumberHere==body2_)then         !^CFG IF SECONDBODY BEGIN
	xBodyHere=xBody2
	yBodyHere=yBody2
	zBodyHere=zBody2
  else                                   !^CFG END SECONDBODY
	xBodyHere=cZero
	yBodyHere=cZero
	zBodyHere=cZero
  end if                                 !^CFG IF SECONDBODY

  !\
  ! Variable meanings:
  !
  ! iter, time_now:  current iteration and current time.  
  ! RhoFaceOutside,pFaceOutside: density and pressure just outside the body
  ! VxFaceOutside,VyFaceOutside,VzFaceOutside: velocity components outside body
  ! BxFaceOutside,ByFaceOutside,BzFaceOutside: B field components outside body
  !
  ! RhoFaceInside,pFaceInside: density and pressure inside the body used to
  !                            set boundary conditions
  ! VxFaceInside,VyFaceInside,VzFaceInside: velocity components inside body
  ! BxFaceInside,ByFaceInside,BzFaceOutside: B field components inside body
  !
  ! BodyNumberHere: the number of the body for which we are applying boundary
  !                 conditions. 
  ! B0xFace,B0yFace,B0zFace - B0 field at the center of the face
  ! XFace, YFace, ZFace: x,y,z location of the current face relative to the
  !                      body center
  ! RhoBodyHere: value set in read_input as the density to assign the body
  ! pBodyHere: value set in read_inputs as the pressure to assign the body
  !
  ! UseIonosphereHere, UseRotatingBcHere:  whether to use the ionosphere and
  !                                        to apply corotation.  
  !/

  ! User declared local variables go here


  !---------------------------------------------------------------------------

  if(iProc==PROCtest .and. globalBLK==BLKtest)then
     call set_oktest('user_set_innerBCs',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if

  ! Calculation of boundary conditions should start here

  !\
  ! The following is an example of the calculation of the inner boundary
  ! condition for the option "ionosphere".  The option UseRotatingBc is
  ! here used, but the UseIonosphere and the calls associated with calling
  ! the ionosphere potential solve as a boundary condition are not shown.
  !
  ! For this case the following definitions would need to be made above
  ! in the user declared local variable section.
  !
  !  real :: VrFaceOutside,VthetaFaceOutside,VphiFaceOutside,&
  !       VrFaceInside,VthetaFaceInside,VphiFaceInside,&
  !       BrFaceOutside,BthetaFaceOutside,BphiFaceOutside, &
  !       BrFaceInside,BthetaFaceInside,BphiFaceInside
  !  real :: cosTheta, sinTheta, cosPhi, sinPhi, RFace
  !  real:: PressureJumpLimit=0.0,DensityJumpLimit=0.1
  !  real, dimension(1:3) :: location, v_phi
  !
  ! Code would be :
  !
  !  !Rotate to spherical coordinates
  !  RFace = sqrt(XFace**2 + YFace**2 + ZFace**2)
  !  cosTheta = ZFace/RFace
  !  sinTheta = sqrt(XFace**2+YFace**2)/RFace
  !  cosPhi = XFace/sqrt(XFace**2+YFace**2+cTolerance**2)
  !  sinPhi = YFace/sqrt(XFace**2+YFace**2+cTolerance**2)
  !  VrFaceOutside = (VxFaceOutside*XFace + &
  !  	    VyFaceOutside*YFace + &
  !  	    VzFaceOutside*ZFace) / RFace
  !  VthetaFaceOutside = ((VxFaceOutside*XFace + &
  !  	    VyFaceOutside*YFace) * ZFace - &
  !  	    VzFaceOutside*(XFace**2+YFace**2)) / &
  !  	    (sqrt(XFace**2+YFace**2+cTolerance**2)*RFace)
  !  VphiFaceOutside = (VyFaceOutside*XFace - &
  !  	    VxFaceOutside*YFace)*sinTheta / &
  !  	    ((XFace**2+YFace**2+cTolerance**2)/RFace)
  !  BrFaceOutside = (BxFaceOutside*XFace + &
  !  	    ByFaceOutside*YFace + &
  !  	    BzFaceOutside*ZFace) / RFace
  !  BthetaFaceOutside = ((BxFaceOutside*XFace + &
  !  	    ByFaceOutside*YFace) * ZFace - &
  !  	    BzFaceOutside*(XFace**2+YFace**2)) / &
  !  	    (sqrt(XFace**2+YFace**2+cTolerance**2)*RFace)
  !  BphiFaceOutside = (ByFaceOutside*XFace - &
  !         BxFaceOutside*YFace)*sinTheta / &
  !  	    ((XFace**2+YFace**2+cTolerance**2)/RFace)
  !
  !  !Apply boundary conditions
  !  VrFaceInside     = - VrFaceOutside
  !  VthetaFaceInside = - VthetaFaceOutside
  !  VphiFaceInside   = - VphiFaceOutside
  !  BrFaceInside     = cZero
  !  BthetaFaceInside = BthetaFaceOutside
  !  BphiFaceInside   = BphiFaceOutside
  !  RhoFaceInside    = RhoFaceOutside+sign(cOne,RhoBodyHere-RhoFaceOutside)*&
  !       min(abs(RhoBodyHere-RhoFaceOutside),DensityJumpLimit*RhoFaceOutside)
  !  pFaceInside      = pFaceOutside+sign(cOne,pBodyHere-pFaceOutside)*&
  !       min(abs(pBodyHere-pFaceOutside),PressureJumpLimit*pFaceOutside)
  !  
  !  ! Rotate back to cartesian coordinates
  !  VxFaceInside  = VrFaceInside*XFace/RFace + VthetaFaceInside*cosTheta*cosPhi - &
  !       VphiFaceInside*sinPhi 
  !  VyFaceInside  = VrFaceInside*YFace/RFace + VthetaFaceInside*cosTheta*sinPhi + &
  !       VphiFaceInside*cosPhi 
  !  VzFaceInside  = VrFaceInside*ZFace/RFace - VthetaFaceInside*sinTheta
  !  
  !  BxFaceInside  = BrFaceInside*XFace/RFace + BthetaFaceInside*cosTheta*cosPhi - &
  !       BphiFaceInside*sinPhi
  !  ByFaceInside  = BrFaceInside*YFace/RFace + BthetaFaceInside*cosTheta*sinPhi + &
  !       BphiFaceInside*cosPhi
  !  BzFaceInside  = BrFaceInside*ZFace/RFace - BthetaFaceInside*sinTheta
  !  
  !  ! apply corotation.  currently works only for the first body 
  !  if (UseRotatingBcHere) then
  !  
  !     location(1) = XFace 
  !     location(2) = YFace 
  !     location(3) = ZFace
  !     !\
  !     ! The program is called which calculates the cartesian corotation velocity vector
  !     ! v_phi as a function of the radius-vector "location"
  !     !/
  !     call calc_corotation_velocities(iter,time_now,location,v_phi)
  !  
  !     VxFaceInside = VxFaceInside + cTwo*v_phi(1)
  !     VyFaceInside = VyFaceInside + cTwo*v_phi(2)   
  !     VzFaceInside = VzFaceInside + cTwo*v_phi(3)
  !  end if
  !/

end subroutine user_set_innerBCs


!========================================================================
!  SUBROUTINE USER_SET_OUTER_BCS 
!========================================================================
! This subroutine allows the user to apply boundary conditions to the outer
! body which are problem specific and cannot be created using the predefined
! options in BATSRUS.
! The variables specific to the problem are loaded from ModUser
! Any of the outer boundary specified in BATSRUS can be used here.

  subroutine user_set_outerBCs(iBlock,TypeBc,found)

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

  logical,intent(out) :: found
  character (len=20),intent(in) :: TypeBc
end subroutine user_set_outerBCs


subroutine user_calc_sources
  use ModAdvance
  use ModUser
  use ModNumConst
  implicit none

  Srho=cZero
  SrhoUx=cZero
  SrhoUy=cZero
  SrhoUz=cZero
  SBx=cZero
  SBy=cZero
  SBz=cZero
  SP=cZero
  SE=cZero
  call user_sources
  Source_VC(rho_,:,:,:)=Srho+ Source_VC(rho_,:,:,:)
  Source_VC(rhoUx_,:,:,:)=SrhoUx+Source_VC(rhoUx_,:,:,:)
  Source_VC(rhoUy_,:,:,:)=SrhoUy+Source_VC(rhoUy_,:,:,:)
  Source_VC(rhoUz_,:,:,:)=SrhoUz+ Source_VC(rhoUz_,:,:,:)
  Source_VC(Bx_,:,:,:)=SBx+ Source_VC(Bx_,:,:,:)
  Source_VC(By_,:,:,:)=SBy+ Source_VC(By_,:,:,:)
  Source_VC(Bz_,:,:,:)=SBz+ Source_VC(Bz_,:,:,:)
  Source_VC(P_,:,:,:)=SP+ Source_VC(P_,:,:,:)
  Source_VC(Energy_,:,:,:)=SE+ Source_VC(Energy_,:,:,:)
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
!      Srho,SrhoUx,SrhoUy,SrhoUz,SBx,SBy,SBz,SE,SP
!
! Note that SE (energy) and SP (pressure) must both be loaded if the code is 
! going to use both the primitive and the conservative MHD equation advance  
! (see the USER MANUAL and the DESIGN document).  If using only primitive SP 
! must be loaded.  If using only conservative SE must be loaded.  The safe
! approach is to load both.
!/

subroutine user_sources
  use ModMain
  use ModAdvance
  use ModGeometry
  use ModPhysics
  use ModProcMH
  use ModUser

  implicit none

  ! Variables required by this user subroutine
  integer :: i,j,k
  logical :: oktest,oktest_me

  !\
  ! Variable meanings:
  !   Srho: Source terms for the continuity equation
  !   SE,SP: Source terms for the energy (conservative) and presure (primative)
  !          equations
  !   SrhoUx,SrhoUy,SrhoUz:  Source terms for the momentum equation
  !   SBx,SBy,SBz:  Souce terms for the magnetic field equations
  !/

  ! User declared local variables go here

  !---------------------------------------------------------------------------
  if(iProc==PROCtest .and. globalBLK==BLKtest)then
     call set_oktest('user_sources',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if

  do k = 1, nK
     do j = 1, nJ
        do i = 1, nI

           ! User coding of the source terms should be inside this loop
           ! Note that user source terms should be added to the the
           ! values already loaded in the arrays.  For example
           !     Srho(i,j,k) = Srho(i,j,k) + user source term for rho

           Srho(i,j,k)   = Srho(i,j,k)   + 0.0
           SrhoUx(i,j,k) = SrhoUx(i,j,k) + 0.0
           SrhoUy(i,j,k) = SrhoUy(i,j,k) + 0.0
           SrhoUz(i,j,k) = SrhoUz(i,j,k) + 0.0
           SBx(i,j,k)    = SBx(i,j,k)    + 0.0
           SBy(i,j,k)    = SBy(i,j,k)    + 0.0
           SBz(i,j,k)    = SBz(i,j,k)    + 0.0
           SE(i,j,k)     = SE(i,j,k)     + 0.0
           SP(i,j,k)     = SP(i,j,k)     + 0.0

        end do     ! end the i loop
     end do     ! end the j loop
  end do     ! end the k loop

end subroutine user_sources



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
!   CellCenteredVars_VGB
!
! Note that in most cases E (energy) and P (pressure) should both be loaded.
!/

subroutine user_initial_perturbation
  use ModMain
  use ModAdvance
  use ModGeometry
  use ModPhysics
  use ModProcMH
  use ModUser

  implicit none

  ! Variables required by this user subroutine
  integer :: i,j,k,iBLK
  logical :: oktest,oktest_me

  !\
  ! Variable meanings:
  !
  ! Current MHD state :
  !   CellCenteredVars_VGB
  !
  !/

  ! User declared local variables go here

  !---------------------------------------------------------------------------
  call set_oktest('user_initial_perturbation',oktest,oktest_me)

  do iBLK = 1,nBLK
     if (unusedBLK(iBLK)) CYCLE   
     do k=1-gcn,nK+gcn
        do j=1-gcn,nJ+gcn
           do i=1-gcn,nI+gcn

           ! User coding of the perturbation to the initial state should
           ! be inside this loop.  The perturbation can be added to the
           ! current state of the solution from the restart file or can
           ! be over written.  The example here adds a perturbation.
           !     State_VGB(rho_,i,j,k) = State_VGB(rho_,i,j,k) + perturbation

   	      State_VGB(rho_,i,j,k,iBLK)   = State_VGB(rho_,i,j,k,iBLK)   + 0.0
   	      State_VGB(rhoUx_,i,j,k,iBLK) = State_VGB(rhoUx_,i,j,k,iBLK) + 0.0
   	      State_VGB(rhoUy_,i,j,k,iBLK) = State_VGB(rhoUy_,i,j,k,iBLK) + 0.0
   	      State_VGB(rhoUz_,i,j,k,iBLK) = State_VGB(rhoUz_,i,j,k,iBLK) + 0.0
   	      State_VGB(Bx_,i,j,k,iBLK)    = State_VGB(Bx_,i,j,k,iBLK)    + 0.0
   	      State_VGB(By_,i,j,k,iBLK)    = State_VGB(By_,i,j,k,iBLK)    + 0.0
   	      State_VGB(Bz_,i,j,k,iBLK)    = State_VGB(Bz_,i,j,k,iBLK)    + 0.
   	      State_VGB(P_,i,j,k,iBLK)     = State_VGB(P_,i,j,k,iBLK)     + 0.0

           end do     ! end the i loop
        end do     ! end the j loop
     end do     ! end the k loop
  end do     ! end the iBLK loop

end subroutine user_initial_perturbation
!==================================================
subroutine user_get_b0(xx,yy,zz,B0_PFSSM)
  implicit none
  real, intent(in):: xx,yy,zz
  real, intent(out), dimension(3):: B0_PFSSM
end subroutine user_get_b0
!==================================================
subroutine user_set_physics_constants
  use ModPhysics,ONLY:FaceState_VI,CellState_VI
end subroutine user_set_physics_constants
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
              minx,miny,minz,minR,maxx,maxy,maxz,maxR,found)
  
  use ModPhysics
  use ModNumConst
  use ModGeometry
  use ModAMR
  use ModUser 

  implicit none

  logical,intent(out) :: refineBlock, found
  integer, intent(in) :: lev
  real, intent(in) :: dxBlock
  real, intent(in) :: xCenter,yCenter,zCenter,rCenter
  real, intent(in) :: minx,miny,minz,minR
  real, intent(in) :: maxx,maxy,maxz,maxR
  integer,intent(in) :: iBLK
end subroutine user_specify_initial_refinement
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

subroutine user_amr_criteria(iBLK, userCriteria, TypeCriteria, found)

!
  use ModMain
  use ModAdvance
  use ModGeometry, ONLY:x_BLK,y_BLK,z_BLK,R_BLK,dx_BLK,dy_BLK,dz_BLK,true_cell
  use ModPhysics
  use ModConst
  use ModUser

 implicit none

 ! Variables required by this user subroutine
  logical ,intent(inout):: found
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
  implicit none
  integer,intent(in)::iStage,iBlock
  call update_states_MHD(iStage,iBlock)	
end subroutine user_update_states

!========================================================================

!========================================================================
!  SUBROUTINE USER_WRITE_PROGRESS
!========================================================================

!
! This subroutine allows the user to write to the log files variables (normalized and not)
! which are problem specific.
!
! The variables specific to the problem are loaded from ModUser


subroutine user_write_progress
 
   use ModProcMH
   use ModMain
   use ModPhysics
   use ModUser
   implicit none

end subroutine user_write_progress

subroutine user_get_log_var(VarValue,TypeVar)

  use ModProcMH
  use ModMain
  use ModVarIndexes
  use ModAdvance
  use ModPhysics
  use ModNumConst
  use ModUser
  implicit none
  real, intent(out):: VarValue
  character (LEN=10), intent(in):: TypeVar

end subroutine user_get_log_var

!real function rho_BLK(i,j,k,iBlock)
!  use ModAdvance
!  implicit none
!  integer::i,j,k,iBlock
!  rho_BLK=State_VGB(rho_,i,j,k,iBlock)
!end function rho_BLK 
!real function rhoUx_BLK(i,j,k,iBlock)
!  use ModAdvance
!  implicit none
!  integer::i,j,k,iBlock
!  rhoUx_BLK=State_VGB(rhoUx_,i,j,k,iBlock)
!end function rhoUx_BLK 
!real function rhoUy_BLK(i,j,k,iBlock)
!  use ModAdvance
!  implicit none
!  integer::i,j,k,iBlock
!  rhoUy_BLK=State_VGB(rhoUy_,i,j,k,iBlock)
!end function rhoUy_BLK 
!real function rhoUz_BLK(i,j,k,iBlock)
!  use ModAdvance
!  implicit none
!  integer::i,j,k,iBlock
!  rhoUz_BLK=State_VGB(rhoUz_,i,j,k,iBlock)
!end function rhoUz_BLK
!real function Bx_BLK(i,j,k,iBlock)
!  use ModAdvance
!  implicit none
!  integer::i,j,k,iBlock
!  Bx_BLK=State_VGB(Bx_,i,j,k,iBlock)
!end function Bx_BLK 
!real function By_BLK(i,j,k,iBlock)
!  use ModAdvance
!  implicit none
!  integer::i,j,k,iBlock
!  By_BLK=State_VGB(By_,i,j,k,iBlock)
!end function By_BLK 
!real function Bz_BLK(i,j,k,iBlock)
!  use ModAdvance
!  implicit none
!  integer::i,j,k,iBlock
!  Bz_BLK=State_VGB(Bz_,i,j,k,iBlock)
!end function Bz_BLK
!real function P_BLK(i,j,k,iBlock)
!  use ModAdvance
!  implicit none
!  integer::i,j,k,iBlock
!  P_BLK=State_VGB(P_,i,j,k,iBlock)
!end function P_BLK
