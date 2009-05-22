!^CFG COPYRIGHT UM
!==============================================================================
module ModUser

  ! This file contains the CCMC version of "specify_refinement.f90"
  ! Other than the usual #GRID command, the PARAM.in file must contain 
  ! the following commands in the GM section:
  !
  ! #USERINPUTBEGIN ----------------
  ! 
  ! ! Define the name of the grid here
  ! CCMCGRID
  ! ror_shock
  ! 
  ! #USERINPUTEND  -----------------
  ! 
  ! ! Define the initial resolution
  ! #GRIDRESOLUTION
  ! 0.5                    Resolution
  ! initial                 NameArea
  ! 
  ! ! A "user" defined area with the highest resolution during the run
  ! #GRIDRESOLUTION
  ! 0.125                    Resolution
  ! user                     NameArea

  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_specify_initial_refinement

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 0.9
  character (len=*), parameter :: &
       NameUserModule = 'CCMC grids, Kuznetsova and Toth'

  ! Name of the grid 
  character(len=30):: NameGrid = "ror6_short"

contains

  !===========================================================================
  subroutine user_read_inputs

    use ModReadParam,   ONLY: read_line, read_command, read_var
    character (len=100) :: NameCommand
    character(len=*), parameter :: NameSub = "ModUser::user_read_inputs"
    !--------------------------------------------------------------------------
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)
       case("#CCMCGRID")
          call read_var("NameGrid", NameGrid)
       case('#USERINPUTEND')
          EXIT
       case default
          call stop_mpi(NameSub//': unknown command name='//trim(NameCommand))
       end select
    end do
  end subroutine user_read_inputs

  !===========================================================================

!!! Above this everything is the same as in ModUserCCMC.f90 of 10/26/2008
!!! except name of subroutine user_specify_refinement(iBlk, iArea, DoRefine)

  subroutine user_specify_initial_refinement(iBLK, DoRefine, lev, &
       Dx0,xxx0,yyy0,zzz0,RR0,minx0,miny0,minz0,minRblk0, &
       maxx0,maxy0,maxz0,maxRblk0, &
       IsFound)

    use ModSize,     ONLY: nI, nJ, nK
    use ModMain,     ONLY: body1, UseRotatingBc
    use ModGeometry, ONLY: XyzMin_D, XyzMax_D, XyzStart_BLK, &
         dy_BLK, dz_BLK, TypeGeometry, x1, x2, x_BLK, y_BLK, z_BLK, dx_BLK
    use ModPhysics,  ONLY: rBody, rCurrents
    use ModAMR,      ONLY: nRefineLevel, InitialRefineType !!!
    use ModNumConst, ONLY: cPi, cTiny

!!!    integer, intent(in) :: iBlk, iArea

    integer, intent(in) :: iBLK, lev                                   !!!
    logical, intent(out):: DoRefine, IsFound                           !!!

    ! These are ignored all
    real,    intent(in) :: Dx0,xxx0,yyy0,zzz0,RR0, &                   !!!
         minx0, miny0, minz0, minRblk0, maxx0, maxy0, maxz0, maxRblk0  !!!

!!!    integer :: lev
    real :: xxx,yyy,zzz,RR, xxPoint,yyPoint,zzPoint
    real :: RRR, RRR1, RRR2, RRR3, RRR4,RRR5
    real :: minRblk,maxRblk,xx1,xx2,yy1,yy2,zz1,zz2, tmpminRblk,tmpmaxRblk
    real :: critx
    real :: x_off,x_off1,x_off2,curve_rad,curve_rad1,curve_rad2
    real :: mach,frac,bound1,bound2
    real :: minx, miny, minz
    real :: SizeMax

    real,parameter::cRefinedTailCutoff=0.75*0.875

    character (len=*), parameter :: NameSub = 'user_specify_refinement'

    logical :: oktest, oktest_me

    !--------------------------------------------------------------------------
    call set_oktest('initial_refinement',oktest,oktest_me)

!!!    lev = nRefineLevel
    IsFound = .true.                     !!!
    NameGrid = InitialRefineType         !!!

!!! Below this everything is the same as in ModUserCCMC.f90 of 10/26/2008

    DoRefine = .false.


    ! Block min, max radius values
    xx1 = 0.5*(x_BLK( 0, 0, 0,iBLK)+x_BLK(   1,   1  , 1,iBLK))
    xx2 = 0.5*(x_BLK(nI,nJ,nK,iBLK)+x_BLK(nI+1,nJ+1,nK+1,iBLK))
    yy1 = 0.5*(y_BLK( 0, 0, 0,iBLK)+y_BLK(   1,   1,   1,iBLK))
    yy2 = 0.5*(y_BLK(nI,nJ,nK,iBLK)+y_BLK(nI+1,nJ+1,nK+1,iBLK))
    zz1 = 0.5*(z_BLK( 0, 0, 0,iBLK)+z_BLK(   1,   1,   1,iBLK))
    zz2 = 0.5*(z_BLK(nI,nJ,nK,iBLK)+z_BLK(nI+1,nJ+1,nK+1,iBLK))

    select case(TypeGeometry)
    case('cartesian')
       SizeMax=dx_BLK(iBLK)
       ! Block center coordinates
       xxx = 0.5*(x_BLK(nI,nJ,nK,iBLK)+x_BLK(1,1,1,iBLK)) 
       yyy = 0.5*(y_BLK(nI,nJ,nK,iBLK)+y_BLK(1,1,1,iBLK))
       zzz = 0.5*(z_BLK(nI,nJ,nK,iBLK)+z_BLK(1,1,1,iBLK))         
       RR = sqrt( xxx*xxx + yyy*yyy + zzz*zzz )
       minRblk = sqrt(&
            minmod(xx1,xx2)**2 + minmod(yy1,yy2)**2 + minmod(zz1,zz2)**2)
       maxRblk = sqrt((max(abs(xx1),abs(xx2)))**2 + &
            (max(abs(yy1),abs(yy2)))**2 + &
            (max(abs(zz1),abs(zz2)))**2)
       if(body1.and.maxRblk<rBody) RETURN
    case('spherical')
       SizeMax=max(dx_BLK(iBLK),dy_BLK(iBLK)*maxRblk)
       minRblk = XyzStart_BLK(1,iBLK)-0.5*dx_BLK(iBLK)            
       maxRblk = minRblk+nI*dx_BLK(iBLK)                            
       RR=0.5*(minRblk+maxRblk)                                   
       xxx=RR*sin(XyzStart_BLK(3,iBLK)+0.5*(nK-1)*dz_BLK(iBLK))*& 
            cos(XyzStart_BLK(2,iBLK)+0.5*(nJ-1)*dy_BLK(iBLK))       
       yyy=RR*sin(XyzStart_BLK(3,iBLK)+0.5*(nK-1)*dz_BLK(iBLK))*& 
            sin(XyzStart_BLK(2,iBLK)+0.5*(nJ-1)*dy_BLK(iBLK))       
       zzz=RR*cos(XyzStart_BLK(3,iBLK)+0.5*(nK-1)*dz_BLK(iBLK))
    case('spherical_lnr')                                          
       SizeMax=max(dx_BLK(iBLK)*maxRblk,dy_BLK(iBLK)*maxRblk)
       minRblk =exp( XyzStart_BLK(1,iBLK)-0.5*dx_BLK(iBLK))            
       maxRblk =exp(nI*dx_BLK(iBLK))*minRblk                            
       RR=0.5*(minRblk+maxRblk)                                   
       xxx=RR*sin(XyzStart_BLK(3,iBLK)+0.5*(nK-1)*dz_BLK(iBLK))*& 
            cos(XyzStart_BLK(2,iBLK)+0.5*(nJ-1)*dy_BLK(iBLK))       
       yyy=RR*sin(XyzStart_BLK(3,iBLK)+0.5*(nK-1)*dz_BLK(iBLK))*& 
            sin(XyzStart_BLK(2,iBLK)+0.5*(nJ-1)*dy_BLK(iBLK))       
       zzz=RR*cos(XyzStart_BLK(3,iBLK)+0.5*(nK-1)*dz_BLK(iBLK))
    case default
       call stop_mpi('Unknown TypeGeometry = '//TypeGeometry)
    end select

    select case (NameGrid)
    case ('none')
       ! Refine no blocks
    case ('all')
       ! Refine all used blocks
       DoRefine = .true.

    case ('2Dbodyfocus')
       ! 

    case ('3Dbodyfocus')
       ! Refine, focusing on body
       if (maxRblk > Rbody) then
          if (lev <= 2) then
             ! Refine all blocks first two times through
             DoRefine = .true.
          else
             ! Refine blocks intersecting body
             if (minRblk < 1.5*Rbody) DoRefine = .true.
          end if
       end if

    case('xplanefocus')
       ! concentrate refinement around x=0 plane
       if (lev <= 2) then
          ! Refine all blocks first two times through
          DoRefine = .true.
       else
          ! Refine blocks intersecting x=0 plane
          if(xx1<=0.and.0<=xx2)DoRefine = .true.
       end if
    case ('spherefocus')
       ! Refine, focusing spherically - like
       ! for a body, but can be used when there is
       ! no body
       if (lev <= 2) then
          ! Refine all blocks first two times through
          DoRefine = .true.
       else
          ! Refine blocks intersecting R=4.5 sphere
          if (minRblk < 4.5) DoRefine = .true.
       end if

    case ('helio_init')
       ! Refine all blocks time through (starting with 1 block)
       if (lev <= 4) then
          DoRefine = .true.
       else
          critx=(XyzMax_D(1)-XyzMin_D(1))/(2.0**real(lev-2))
          if ( RR < 1.10*rBody + critx ) then
             DoRefine = .true.
          else
             DoRefine = .false.
          end if
       endif

    case('helio_z=4')
       ! Refine all blocks time through (starting with 1 block)
       if (lev <= 2) then
          DoRefine = .true.
       else
          !         if ((xxx < 60.0).and.(xxx >-60.0)) then
          !           if ((yyy < 60.0).and.(yyy >-60.0)) then
          if ((SizeMax <=2.0).and.(SizeMax >1.0) &
               .and.(zzz<=10.0).and.(zzz>=-10.0)) then
             DoRefine =.true.
          endif
          !         endif
          !        endif
       endif

    case ('all_then_focus')
       !\
       ! Initial refinement criteria for heliospheric flow problem.
       !/
       ! ***** CME MILESTONE 3 RUN (March 1999) *****
       ! Starts from one block and refines until 512 block
       ! focusing grid is made.
       if (lev <= 2) then
          ! Refine all blocks 2 times through
          DoRefine = .true.
       else
          if (lev <= 14) then
             critx=(XyzMax_D(1)-XyzMin_D(1))/(2.0**real(lev-1))
             if ( RR < 1.10*rBody + critx ) then
                DoRefine = .true.
             else
                DoRefine = .false.
             end if
          end if
       endif

    case ('points')
       ! Refine around given points
       if(lev<=3)then
          ! Refine all blocks first 3 times through
          DoRefine = .true.
       else
          ! List as many points as desired
          xxPoint =  0.
          yyPoint =  0.
          zzPoint =  0.
          RR = sqrt( (xxx-xxPoint)**2 + (yyy-yyPoint)**2 + (zzz-zzPoint)**2 )
          if (RR < 2.*real(nI)*dx_BLK(iBLK)) DoRefine = .true.
       end if

                  case ('ds1')
           ! 6 x 6 x 6 block for High resolution Runs On Request
           ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
           if (maxRblk > Rbody) then

              if (dx_BLK(iBLK) > 8.) then
                 ! Refine all blocks with dx greater than 8
                DoRefine = .true.
              else

                 if (dx_BLK(iBLK) > 4. .and. xxx > -159.0) DoRefine = .true.

                 if (dx_BLK(iBLK) > 2. .and. xxx > -39.0) DoRefine = .true.

                 if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                   abs(yyy)<48..and. abs(zzz)<24.) DoRefine = .true.
!    1 Re
                 if (dx_BLK(iBLK)>1. .and. xxx > -27. .and. xxx < 33.0 .and. &
                    abs(yyy)<48..and. abs(zzz)<36.)DoRefine = .true.

               if (dx_BLK(iBLK)>1. .and. xxx > -51. .and. xxx < -27 .and. &
                    abs(yyy)<48..and. abs(zzz)<24.)DoRefine = .true.

                 if (dx_BLK(iBLK)>1. .and. xxx < -51. .and. xxx > -75. .and. &
                    abs(yyy)<24..and. abs(zzz)<12.)DoRefine = .true.
!    0.5 Re
                  if (dx_BLK(iBLK)>0.5 .and. xxx < 21.0 .and. xxx > -27. .and. &
                    abs(yyy)<30. .and. abs(zzz)<21.)DoRefine = .true.

                  if (dx_BLK(iBLK)>0.5 .and. xxx < -27.0 .and. xxx > - 45. .and. &
                    abs(yyy)<30. .and. abs(zzz)<15.)DoRefine = .true.

                  if (dx_BLK(iBLK)>0.5 .and. xxx < - 45. .and. xxx> - 63.0 .and. &
                    abs(yyy)<12..and. abs(zzz)<6.)DoRefine = .true.
!    0.25 Re
                 if (dx_BLK(iBLK)>0.25 .and. xxx < 18. .and. xxx> -6.0 .and. &
                    abs(yyy)<9. .and. abs(zzz)<9.)DoRefine = .true.

                  if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> 9.0 .and. &
                    abs(yyy)<15. .and. abs(zzz)<12.)DoRefine = .true.

                   if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                    abs(yyy)<21. .and. abs(zzz)<15.)DoRefine = .true.

                   if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -36.0 .and. &
                    abs(yyy)<21. .and. abs(zzz)<12.)DoRefine = .true.

!!!

                  if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                    abs(yyy)<6..and. abs(zzz)<6.)DoRefine = .true.
!  0.125 Re
!  0.125 Re

!!!                   if (dx_BLK(iBLK)>0.125 .and.  xxx > -24. .and.  xxx < -7.5 .and. &
!!!                  abs(yyy) < 19.5 .and. abs(yyy) > 13.5  .and. &
!!!                   abs(zzz)< 6.)DoRefine = .true.
!!!
!!!                   if (dx_BLK(iBLK)>0.125 .and.  xxx > -7.5 .and.  xxx < -2.5 .and. &
!!!                  abs(yyy) < 19.5 .and. abs(yyy) > 12.  .and. &
!!!                   abs(zzz)< 6)DoRefine = .true.


                   if (dx_BLK(iBLK)>0.125 .and.  xxx > 4.5 .and.  xxx < 12 .and. &
                  abs(yyy) < 6.  .and.  zzz < 4. .and. &
                   zzz > -10.)DoRefine = .true.


                  RRR=sqrt((xxx+3.)**2+1.3*zzz**2+(5.+0.77*abs(yyy))**2-25.)

                  if (dx_BLK(iBLK)>0.125 .and.  xxx > -3. .and.  xxx < 14. .and. &
                  RRR >10. .and.  RRR < 17. .and.  zzz < 6. .and. zzz > -12. &
                   )DoRefine = .true.
!
!
!!!! corrections on 04/14/07

                  if (dx_BLK(iBLK)>0.125 .and.  xxx > 1. .and.  xxx < 14. .and. &
                  RRR >10. .and.  RRR < 20. .and. zzz < 6. .and. zzz > -12. &
                   )DoRefine = .true.
!!!!!         end corrections 04/14/07
! 0.0625

!!old                  if (dx_BLK(iBLK)>0.0625 .and.  xxx > -18. .and.  xxx < -7.5 .and. &
!!old                  abs(yyy) < 18. .and. abs(yyy) > 15.  .and. &
!!old                   abs(zzz)< 3.)DoRefine = .true.

!!                   if (dx_BLK(iBLK)>0.0625 .and.  xxx > -7.5 .and.  xxx < -0.75 .and. &
!!                  abs(yyy) < 18. .and. abs(yyy) > 13.5  .and. &
!!                   abs(zzz)< 3.)DoRefine = .true.

               if (dx_BLK(iBLK)>0.0625 .and.  xxx > 6. .and.  xxx < 12. .and. &
                  abs(yyy) < 6.  .and.  zzz < 3.5 .and. &
                   zzz > -8.5)DoRefine = .true.


               if (dx_BLK(iBLK)>0.0625 .and.  xxx > 0 .and.  xxx < 12. .and. &
                 RRR > 11. .and. RRR < 15.5 .and.  zzz < 3.5 .and. &
                  zzz > -8.5)DoRefine = .true.
!!! corrections 04/14/07
               if (dx_BLK(iBLK)>0.0625 .and.  xxx > 1. .and.  xxx < 14. .and. &
                 RRR > 11. .and. RRR < 18.5 .and.  zzz < 3.5 .and. &
                  zzz > -8.5)DoRefine = .true.

!!! end 04/14/07
!!!
!!!
!!! corrections 04/15/07
               if (dx_BLK(iBLK)>0.0625 .and.  xxx > 6. .and.  xxx < 12. .and. &
                 abs(yyy) < 7.5 .and.  zzz < 3.5 .and. &
                  zzz > -8.5)DoRefine = .true.
!!! end 04/15/07
              end if
           end if





         case ('tn1')
           ! 6 x 6 x 6 block for High resolution Runs On Request
           ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
           if (maxRblk > Rbody) then

              if (dx_BLK(iBLK) > 8.) then
                 ! Refine all blocks with dx greater than 8
                DoRefine = .true.
              else

                 if (dx_BLK(iBLK) > 4. .and. xxx > -159.0) DoRefine = .true.

                 if (dx_BLK(iBLK) > 2. .and. xxx > -39.0) DoRefine = .true.

                 if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                   abs(yyy)<48..and. abs(zzz)<24.) DoRefine = .true.
!    1 Re
                 if (dx_BLK(iBLK)>1. .and. xxx > -27. .and. xxx < 33.0 .and. &
                    abs(yyy)<48..and. abs(zzz)<36.)DoRefine = .true.

               if (dx_BLK(iBLK)>1. .and. xxx > -51. .and. xxx < -27 .and. &
                    abs(yyy)<48..and. abs(zzz)<24.)DoRefine = .true.

                 if (dx_BLK(iBLK)>1. .and. xxx < -51. .and. xxx > -75. .and. &
                    abs(yyy)<24..and. abs(zzz)<12.)DoRefine = .true.
!    0.5 Re
                  if (dx_BLK(iBLK)>0.5 .and. xxx < 21.0 .and. xxx > -27. .and. &
                    abs(yyy)<30. .and. abs(zzz)<21.)DoRefine = .true.

                  if (dx_BLK(iBLK)>0.5 .and. xxx < -27.0 .and. xxx > - 45. .and. &
                    abs(yyy)<30. .and. abs(zzz)<15.)DoRefine = .true.

                  if (dx_BLK(iBLK)>0.5 .and. xxx < - 45. .and. xxx> - 63.0 .and. &
                    abs(yyy)<12..and. abs(zzz)<6.)DoRefine = .true.
!    0.25 Re
                 if (dx_BLK(iBLK)>0.25 .and. xxx < 18. .and. xxx> -6.0 .and. &
                    abs(yyy)<9. .and. abs(zzz)<9.)DoRefine = .true.

                  if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> 9.0 .and. &
                    abs(yyy)<15. .and. abs(zzz)<12.)DoRefine = .true.

                   if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                    abs(yyy)<21. .and. abs(zzz)<15.)DoRefine = .true.

                   if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -36.0 .and. &
                    abs(yyy)<21. .and. abs(zzz)<12.)DoRefine = .true.

!!!

                  if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                    abs(yyy)<6..and. abs(zzz)<6.)DoRefine = .true.
!  0.125 Re
                  if (dx_BLK(iBLK)>0.125 .and.  xxx > 4.5 .and.  xxx < 12 .and. &
                  abs(yyy) < 6.  .and.  abs(zzz) < 9.)DoRefine = .true.


                  RRR=sqrt((xxx+3.)**2+1.3*zzz**2+(5.+0.77*abs(yyy))**2-25.)

                  if (dx_BLK(iBLK)>0.125 .and.  xxx > -3. .and.  xxx < 14. .and. &
                  RRR >10. .and.  RRR < 17. .and.  abs(zzz) < 9. )DoRefine = .true.
!
!
!!!! corrections on 04/14/07

                  if (dx_BLK(iBLK)>0.125 .and.  xxx > 1. .and.  xxx < 14. .and. &
                  RRR >10. .and.  RRR < 20. .and. abs(zzz) < 9. )DoRefine = .true.
!!!!!         end corrections 04/14/07
! 0.0625

!!old                  if (dx_BLK(iBLK)>0.0625 .and.  xxx > -18. .and.  xxx < -7.5 .and. &
!!old                  abs(yyy) < 18. .and. abs(yyy) > 15.  .and. &
!!old                   abs(zzz)< 3.)DoRefine = .true.

!!                   if (dx_BLK(iBLK)>0.0625 .and.  xxx > -7.5 .and.  xxx < -0.75 .and. &
!!                  abs(yyy) < 18. .and. abs(yyy) > 13.5  .and. &
!!                   abs(zzz)< 3.)DoRefine = .true.

               if (dx_BLK(iBLK)>0.0625 .and.  xxx > 6. .and.  xxx < 12. .and. &
                  abs(yyy) < 6.  .and.  abs(zzz) < 4.5)DoRefine = .true.


               if (dx_BLK(iBLK)>0.0625 .and.  xxx > 0 .and.  xxx < 12. .and. &
                 RRR > 11. .and. RRR < 15.5 .and. abs(zzz) < 4.5)DoRefine = .true.
!!! corrections 04/14/07
               if (dx_BLK(iBLK)>0.0625 .and.  xxx > 1. .and.  xxx < 14. .and. &
                 RRR > 11. .and. RRR < 18.5 .and.  abs(zzz) < 4.5)DoRefine = .true.

!!! end 04/14/07
!!!
!!!
!!! corrections 04/15/07
               if (dx_BLK(iBLK)>0.0625 .and.  xxx > 6. .and.  xxx < 12. .and. &
                 abs(yyy) < 7.5 .and.  abs(zzz) < 4.5 )DoRefine = .true.
!!! end 04/15/07
              end if
           end if





    case ('ror6_short')
       ! 6 x 6 x 6 block for pressure pulse high resolution stady
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -87.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 21.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -51. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 9.0 .and. xxx > -3. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.                     

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true. 

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -27.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.      


             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. abs(xxx) < 3.5 .and. &
                  abs(yyy)<3.5 .and. abs(zzz)<3.5) DoRefine = .true.


             ! 2 Re:  X>-84
             ! 1 Re:  |Z|<24, |Y|<24, -48<X<24
             ! 0.5Re:  -24<X<12
             ! 0.25Re
             ! 0.25Re


             ! Refine tail sheet

             ! This will make it:
             !    


          end if
       end if
       !
    case ('msh_grid1')
       ! 6 x 6 x 6 block for pressure pulse high resolution stady
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.


             if (dx_BLK(iBLK) > 2. .and. xxx > -87.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<48.)  DoRefine = .true.


             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<36..and. abs(zzz)<36.) DoRefine = .true.


             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -51. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.



             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
                  abs(yyy)<30..and. abs(zzz)<30.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > 15. .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -27.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.




             ! 2 Re:  X>-84
             ! 1 Re:  |Z|<24, |Y|<24, -48<X<24
             ! 0.5Re:  -24<X<12
             ! 0.25Re
             ! 0.25Re


             ! Refine tail sheet

             ! This will make it:
             !


          end if
       end if
       !
    case ('msh_grid2')
       ! 6 x 6 x 6 block for pressure pulse high resolution stady
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.


             if (dx_BLK(iBLK) > 2. .and. xxx > -87.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<48.)  DoRefine = .true.


             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<36..and. abs(zzz)<36.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -51. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
                  abs(yyy)<30..and. abs(zzz)<30.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > 15. .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -27.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx > -6. .and. &
                  abs(yyy)<27..and. abs(zzz)<27.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 18. .and. xxx > 12. .and. &
                  abs(yyy)<21..and. abs(zzz)<21.) DoRefine = .true.




             ! 2 Re:  X>-84
             ! 1 Re:  |Z|<24, |Y|<24, -48<X<24
             ! 0.5Re:  -24<X<12
             ! 0.25Re
             ! 0.25Re


             ! Refine tail sheet

             ! This will make it:
             !


          end if
       end if
       !

    case ('jm6')
       ! 6 x 6 x 6 block for pressure pulse high resolution stady
       ! Box sizes: -228<X<60, -96<Z<96, -96<Y<96
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -132.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<48.)  DoRefine = .true.

             if (dx_BLK(iBLK) > 4. .and. xxx > -84.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -12.0 .and. &
                  abs(yyy)<72..and. abs(zzz)<72.)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -36.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<48.)  DoRefine = .true.                   

             if (dx_BLK(iBLK) > 2. .and. xxx > -84.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -24. .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -48. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -12. .and. &
                  abs(yyy)<36..and. abs(zzz)<36.) DoRefine = .true.
             if (dx_BLK(iBLK)>1. .and. xxx > 0. .and. &
                  abs(yyy)<48..and. abs(zzz)<36.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.5 .and. xxx > -12. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

          end if
       end if

    case ('ror_shock')
       ! 6 x 6 x 6 block for pressure pulse high resolution stady
       ! Box sizes: -264<X<24, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -168.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -48.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -96.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -24. .and. xxx < 24.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<48.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -48. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 24.0 .and. xxx > -18. .and. &
                  abs(yyy)<36..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 18. .and. xxx > -18. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 24. .and. xxx> -12.0 .and. &
                  abs(yyy)<30..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 12.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. xxx < 15.0 .and. xxx> 0.0  .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. xxx < 12.0 .and. xxx> 0.0  .and. &
                  abs(yyy)<18..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. xxx < 9.0 .and. xxx> 0.0  .and. &
                  abs(yyy)<24..and. abs(zzz)<6.) DoRefine = .true.




             ! 2 Re:  X>-84
             ! 1 Re:  |Z|<24, |Y|<24, -48<X<24
             ! 0.5Re:  -24<X<12
             ! 0.25Re
             ! 0.25Re


             ! Refine tail sheet

             ! Refine tail sheet

             ! This will make it:
             !


          end if
       end if


    case ('ror6_long')
       ! 6 x 6 x 6 block for overview Runs On Request (long tail added)
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 21.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 9.0 .and. xxx > -3. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.                     

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true. 

             if (dx_BLK(iBLK)>0.5 .and. xxx < 0. .and. xxx> -99. .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.      

             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. abs(xxx) < 5. .and. &
                  abs(yyy)<5..and. abs(zzz)<5.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.06125 .and. abs(xxx) < 4.5 .and. &
                  abs(yyy)<4.5.and. abs(zzz)<4.5) DoRefine = .true.

             ! This will make it:
             !                                            

             ! Tail
             ! 4 Re:  X>-159
             ! 2 Re:  X>-135
             ! 1 Re:  |Z|<12, |Y|<24, -111<X<0
             ! 0.5Re: |Z|<6, |Y|<12, -99<X<0

             ! Magnetopause/Shock
             ! 1 Re:   |Z|<24, |Y|<24, -15 <X<21 
             ! 0.5Re:  |Z|<12, |Y|<12, 9<X<15
             ! 0.5Re:  |Z|<18, |Y|<18, -3<X<9

             ! Near Earth

             ! 0.5Re:  |Z|<12, |Y|<15, -9<X<15  
             ! 0.25Re: |Z|<6, |Y|<6, |X|<6,

          end if
       end if

    case ('mn1')
       ! 6 x 6 x 6 block for overview Runs On Request (long tail added)
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -168.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -72.0)  DoRefine = .true.


             if (dx_BLK(iBLK)>1. .and. xxx > -60. .and. xxx < 36.0 .and. &
                  abs(yyy)<60. .and. abs(zzz)<48.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.5 .and. xxx < 30. .and. xxx> -54. .and. &
                  abs(yyy)<48..and. abs(zzz)<36.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

!             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 24. .and. yyy > -24. .and. &
!                  yyy<30. .and. abs(zzz)<5.) DoRefine = .true.
!
!             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 24. .and. yyy > -24. .and. &
!                  yyy<30. .and. abs(zzz)<5.) DoRefine = .true.



          end if
       end if

    case ('mn2')
       ! 6 x 6 x 6 block for overview Runs On Request (long tail added)
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -168.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -72.0)  DoRefine = .true.


             if (dx_BLK(iBLK)>1. .and. xxx > 0. .and. xxx < 36.0 .and. &
                  yyy> -60. .and. yyy < 36.   .and. abs(zzz)<36.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -12. .and. xxx < 0 .and. &
                 yyy> -72. .and. yyy < 36.  .and. abs(zzz)<36.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -60. .and. xxx < -12 .and. &
                 yyy> -168. .and. yyy < 36.   .and. abs(zzz)<36.) DoRefine = .true.
  

             if (dx_BLK(iBLK)>0.5 .and. xxx < 30. .and. xxx> 0. .and. &
                  yyy> -48. .and. yyy < 24. .and. abs(zzz)<30.) DoRefine = .true.

            if (dx_BLK(iBLK)>0.5 .and. xxx < 0. .and. xxx> -12. .and. &
                 yyy> -60. .and. yyy < 30.  .and. abs(zzz)<30.) DoRefine = .true.

            if (dx_BLK(iBLK)>0.5 .and. xxx < -10. .and. xxx> -54. .and. &
                 yyy> -156. .and. yyy < 30.  .and. abs(zzz)<30.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

!             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 24. .and. yyy > -24. .and. &
!                  yyy<30. .and. abs(zzz)<5.) DoRefine = .true.
!
!             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 24. .and. yyy > -24. .and. &
!                  yyy<30. .and. abs(zzz)<5.) DoRefine = .true.



          end if
       end if


    case ('mn2_old')
       ! 6 x 6 x 6 block for overview Runs On Request (long tail added)
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -168.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -72.0)  DoRefine = .true.


             if (dx_BLK(iBLK)>1. .and. xxx > 0. .and. xxx < 36.0 .and. &
                  yyy> -60. .and. yyy < 36.   .and. abs(zzz)<36.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -10. .and. xxx < 0 .and. &
                 yyy> -72. .and. yyy < 36.  .and. abs(zzz)<36.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -60. .and. xxx < -10 .and. &
                 yyy> -96. .and. yyy < 36.   .and. abs(zzz)<36.) DoRefine = .true.



             if (dx_BLK(iBLK)>0.5 .and. xxx < 30. .and. xxx> 0. .and. &
                  yyy> -48. .and. yyy < 24. .and. abs(zzz)<30.) DoRefine = .true.

            if (dx_BLK(iBLK)>0.5 .and. xxx < 0. .and. xxx> -10. .and. &
                 yyy> -60. .and. yyy < 30.  .and. abs(zzz)<30.) DoRefine = .true.

            if (dx_BLK(iBLK)>0.5 .and. xxx < -10. .and. xxx> -54. .and. &
                 yyy> -96. .and. yyy < 30.  .and. abs(zzz)<30.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

!             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 24. .and. yyy > -24. .and. &
!                  yyy<30. .and. abs(zzz)<5.) DoRefine = .true.
!
!             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 24. .and. yyy > -24. .and. &
!                  yyy<30. .and. abs(zzz)<5.) DoRefine = .true.



          end if
       end if
 case ('ror6_RB_lr')
       ! 6 x 6 x 6 block for overview Runs On Request (long tail added)
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -27. .and. xxx < 21.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -87. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 9.0 .and. xxx > -3. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx > -15. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 0. .and. xxx> -51. .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 12. .and. &
                  abs(yyy)<12. .and. abs(zzz)<12. ) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. abs(xxx) < 5. .and. &
                  abs(yyy)<5..and. abs(zzz)<5.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.06125 .and. abs(xxx) < 4.5 .and. &
                  abs(yyy)<4.5.and. abs(zzz)<4.5) DoRefine = .true.

          end if
       end if



    case ('ror6_long2')
       ! 6 x 6 x 6 block for overview Runs On Request (long tail added)
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 21.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 9.0 .and. xxx > -3. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.                     

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true. 

             if (dx_BLK(iBLK)>0.5 .and. xxx < 0. .and. xxx> -99. .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.      

             ! increase near-Earth region to 8 R_E                    
             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 8.0 .and. &
                  abs(yyy)<8..and. abs(zzz)<8.) DoRefine = .true.

             ! increase resolution in near-Earth tail                    
             !                  if (dx_BLK(iBLK)>0.25 .and. xxx > 40 .and. xxx < -6.0 .and. &
             !                    abs(yyy)<5..and. abs(zzz)<2.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. abs(xxx) < 5. .and. &
                  abs(yyy)<5..and. abs(zzz)<5.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.06125 .and. abs(xxx) < 4.5 .and. &
                  abs(yyy)<4.5.and. abs(zzz)<4.5) DoRefine = .true.

             ! This will make it:
             !                                            

             ! Tail
             ! 4 Re:  X>-159
             ! 2 Re:  X>-135
             ! 1 Re:  |Z|<12, |Y|<24, -111<X<0
             ! 0.5Re: |Z|<6, |Y|<12, -99<X<0

             ! Magnetopause/Shock
             ! 1 Re:   |Z|<24, |Y|<24, -15 <X<21 
             ! 0.5Re:  |Z|<12, |Y|<12, 9<X<15
             ! 0.5Re:  |Z|<18, |Y|<18, -3<X<9

             ! Near Earth

             ! 0.5Re:  |Z|<12, |Y|<15, -9<X<15  
             ! 0.25Re: |Z|<6, |Y|<6, |X|<6,

          end if
       end if


    case ('ror6_long3')
       ! 6 x 6 x 6 block for overview Runs On Request (long tail added)
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 21.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 9.0 .and. xxx > -3. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.                     

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true. 

             if (dx_BLK(iBLK)>0.5 .and. xxx < 0. .and. xxx> -99. .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.      

             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             ! increase resolution in near-Earth tail                    
             if (dx_BLK(iBLK)>0.25 .and. xxx > -30 .and. xxx < -6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. abs(xxx) < 5. .and. &
                  abs(yyy)<5..and. abs(zzz)<5.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.06125 .and. abs(xxx) < 4.5 .and. &
                  abs(yyy)<4.5.and. abs(zzz)<4.5) DoRefine = .true.

             ! This will make it:
             !                                            

             ! Tail
             ! 4 Re:  X>-159
             ! 2 Re:  X>-135
             ! 1 Re:  |Z|<12, |Y|<24, -111<X<0
             ! 0.5Re: |Z|<6, |Y|<12, -99<X<0

             ! Magnetopause/Shock
             ! 1 Re:   |Z|<24, |Y|<24, -15 <X<21 
             ! 0.5Re:  |Z|<12, |Y|<12, 9<X<15
             ! 0.5Re:  |Z|<18, |Y|<18, -3<X<9

             ! Near Earth

             ! 0.5Re:  |Z|<12, |Y|<15, -9<X<15  
             ! 0.25Re: |Z|<6, |Y|<6, |X|<6,
          endif
       endif

    case ('ror6_hr_3')
       ! 6 x 6 x 6 block for overview Runs On Request (long tail added)
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else
 
RRR=sqrt(xxx*xxx+yyy*yyy+zzz*zzz)

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.
                   
             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 21.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 9.0 .and. xxx > -3. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.
                   
             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 0. .and. xxx> -99. .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. RRR < 18. ) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. RRR < 12. ) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR< 7.5) DoRefine = .true.

             if (dx_BLK(iBLK)>0.06125 .and. RRR< 5.25) DoRefine = .true.


          end if
       end if

    case ('ror6_hr_4')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -99.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 0. .and. xxx> -33.0 .and. &
                  abs(yyy)<18..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 51.0 .and. &
                  abs(yyy)<9..and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 24.0 .and. &
                  abs(yyy)<15..and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx> 0.0 .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.

              RRR=sqrt(xxx*xxx+yyy*yyy+zzz*zzz)

             if (dx_BLK(iBLK)>0.5 .and. RRR < 18. ) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. RRR < 12. ) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR< 7.5) DoRefine = .true.

             if (dx_BLK(iBLK)>0.06125 .and. RRR< 5.25) DoRefine = .true.

          end if
       end if



    case ('mag6_rt')
       ! 6 x 6 x 6 block for real time runs
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -87.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 21.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -51. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx > -3. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -27.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.



          end if
       end if


    case ('ror6_rb')
       ! 6 x 6 x 6 block for overview Runs On Request (long tail added)
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 21.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 9.0 .and. xxx > -3. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 0. .and. xxx> -99. .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 9. .and. xxx> -24. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             !                  if (dx_BLK(iBLK)>0.25 .and. xxx < 3. .and. xxx > -18. .and. &
             !                    yyy < 15. .and. yyy > -18. .and. abs(zzz)<9.) DoRefine = .true.
             !
             !                    if (dx_BLK(iBLK)>0.125 .and. xxx < 0. .and. xxx > -15. .and. &
             !                  yyy < -7. .and. yyy > -14. .and. zzz < 4. .and. zzz > -7.) &
             !                     DoRefine = .true.
             !


             ! This will make it:
             !

             ! Tail
             ! 4 Re:  X>-159
             ! 2 Re:  X>-135
             ! 1 Re:  |Z|<12, |Y|<24, -111<X<0
             ! 0.5Re: |Z|<6, |Y|<12, -99<X<0

             ! Magnetopause/Shock
             ! 1 Re:   |Z|<24, |Y|<24, -15 <X<21
             ! 0.5Re:  |Z|<12, |Y|<12, 9<X<15
             ! 0.5Re:  |Z|<18, |Y|<18, -3<X<9

             ! Near Earth

             ! 0.5Re:  |Z|<12, |Y|<15, -9<X<15
             ! 0.25Re: |Z|<6, |Y|<6, |X|<6,

          end if
       end if


    case ('ror6_tail')
       ! 6 x 6 x 6 block for overview Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 21.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 9.0 .and. xxx > -3. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.                     

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true. 

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -99.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 93.0 .and. &
                  abs(yyy)<9..and. abs(zzz)<3.) DoRefine = .true.      


             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.


             ! This will make it:
             !                                            

             ! Tail
             ! 4 Re:   X>-159
             ! 2 Re:   X>-135
             ! 1 Re:   |Z|<12, |Y|<24, -111<X<0
             ! 0.5Re:  |Z|<6, |Y|<12, -99<X<0
             ! 0.25Re: |Z|<3, |Y|<9, -93<X<0

             ! Magnetopause/Shock
             ! 1 Re:   |Z|<24, |Y|<24, -15 <X<21 
             ! 0.5Re:  |Z|<12, |Y|<12, 9<X<15
             ! 0.5Re:  |Z|<18, |Y|<18, -3<X<9

             ! Near Earth

             ! 0.5Re:  |Z|<12, |Y|<15, -9<X<15  
             ! 0.25Re: |Z|<6, |Y|<6, |X|<6,   


          end if
       end if
    case ('ror6_kn_2')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -96<Z<96, -96<Y<96
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<48.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -51. .and. xxx < 33.0 .and. &
                  abs(yyy)<36..and. abs(zzz)<36.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -63. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.5 .and. xxx < 6.0 .and. xxx > -39. .and. &
                  abs(yyy)<36..and. abs(zzz)<24.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -39. .and. &
                  abs(yyy)<36..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.25 .and. xxx < 6. .and. &
                  xxx > - 21.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 21.0 .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx> 0.0 .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.




             ! This will make it:
             !
             ! Tail
             ! 4 Re:  X>-159
             ! 2 Re:  X>-135
             ! 1 Re:  |Z|<12, |Y|<24, -63<X<0
             ! 1 Re:  |Z|<24, |Y|<36, -51<X<0
             ! 0.5Re: |Z|<12, |Y|<24, -39<X<0
             ! 0.25   |Z|<6, |Y|<21, -21<X<3

             ! Nightside magnetopause


             ! Magnetopause/Shock
             ! 1 Re:   |Z|<24, |Y|<24, -15 <X< 33
             ! 0.5Re:  |Z|<18, |Y|<18, -3<X<15
             ! 0.5Re:  |Z|<12, |Y|<12, -9<X<21
             ! 0.25Re: |Z|<9, |Y|<9, -6<X<15
             ! 0.25Re: |Z|<12, |Y|<12, 0<X<12
             ! 0.25Re: |Z|<15, |Y|<15, -3<X<9

             ! Near Earth
             ! 0.25Re: |Z|<6, |Y|<6, |X|<6


          end if
       end if

    case ('ror6_kn_1')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -96<Z<96, -96<Y<96
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<48.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -51. .and. xxx < 33.0 .and. &
                  abs(yyy)<36..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -63. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -39. .and. &
                  abs(yyy)<36..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.25 .and. xxx < 6. .and. &
                  xxx > - 21.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<6.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx> 0.0 .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and.  xxx < 3. .and. xxx > -18. .and. &
                  yyy < 19.5 .and. yyy > 10. .and. & 
                  abs(zzz)< 3) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and.  xxx < 0. .and. xxx > -15. .and. &
                  yyy < 17. .and. yyy > 13. .and. & 
                  abs(zzz)< 1.5) DoRefine = .true.



             ! This will make it:
             !
             ! Tail
             ! 4 Re:  X>-159
             ! 2 Re:  X>-135
             ! 1 Re:  |Z|<12, |Y|<24, -63<X<0
             ! 1 Re:  |Z|<24, |Y|<36, -51<X<0
             ! 0.5Re: |Z|<12, |Y|<24, -39<X<0
             ! 0.25   |Z|<6, |Y|<21, -21<X<3

             ! Nightside magnetopause


             ! Magnetopause/Shock
             ! 1 Re:   |Z|<24, |Y|<24, -15 <X< 33
             ! 0.5Re:  |Z|<18, |Y|<18, -3<X<15
             ! 0.5Re:  |Z|<12, |Y|<12, -9<X<21
             ! 0.25Re: |Z|<9, |Y|<9, -6<X<15
             ! 0.25Re: |Z|<12, |Y|<12, 0<X<12
             ! 0.25Re: |Z|<15, |Y|<15, -3<X<9

             ! Near Earth
             ! 0.25Re: |Z|<6, |Y|<6, |X|<6


          end if
       end if

    case ('ror6_mp')
       ! 6 x 6 x 6 block for overview Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 21.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -3. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.                     

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true. 

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -99.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx> 0.0 .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> 0.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.



             ! This will make it:
             !                                            

             ! Tail
             ! 4 Re:  X>-159
             ! 2 Re:  X>-135
             ! 1 Re:  |Z|<12, |Y|<12, -111<X<0
             ! 0.5Re: |Z|<6, |Y|<6, -99<X<0

             ! Magnetopause/Shock
             ! 1 Re:   |Z|<24, |Y|<24, -15 <X<21 
             ! 0.5Re:  |Z|<18, |Y|<18, -3<X<15
             ! 0.25Re: |Z|<12, |Y|<12, 0<X<12
             ! 0.25Re: |Z|<15, |Y|<15, 0<X<9

             ! Near Earth

             ! 0.5Re:  |Z|<12, |Y|<15, -9<X<15  
             ! 0.25Re: |Z|<6, |Y|<6, |X|<6

             !    


          end if
       end if

    case ('ror6_mp2')
       ! 6 x 6 x 6 block for overview Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -3. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.                     

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true. 

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -99.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx> 0.0 .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> 0.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> 0.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             ! This will make it:
             !                                            
             ! Tail
             ! 4 Re:  X>-159
             ! 2 Re:  X>-135
             ! 1 Re:  |Z|<12, |Y|<12, -111<X<0
             ! 0.5Re: |Z|<6, |Y|<6, -99<X<0

             ! Magnetopause/Shock
             ! 1 Re:   |Z|<24, |Y|<24, -15 <X< 33
             ! 0.5Re:  |Z|<18, |Y|<18, -3<X<15
             ! 0.5Re:  |Z|<12, |Y|<12, -9<X<21
             ! 0.25Re: |Z|<9, |Y|<9, 0<X<15
             ! 0.25Re: |Z|<12, |Y|<12, 0<X<12
             ! 0.25Re: |Z|<15, |Y|<15, 0<X<9

             ! Near Earth
             ! 0.25Re: |Z|<6, |Y|<6, |X|<6

          end if
       end if
    case ('ror6_mp3')
       ! 6 x 6 x 6 block for High Resolution Runs On Request (magnetopause/shock/cusp)
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -99.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx> 0.0 .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -6.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             ! This will make it:
             !
             ! Tail
             ! 4 Re:  X>-159
             ! 2 Re:  X>-135
             ! 1 Re:  |Z|<12, |Y|<12, -111<X<0
             ! 0.5Re: |Z|<6, |Y|<6, -99<X<0

             ! Magnetopause/Shock
             ! 1 Re:   |Z|<24, |Y|<24, -15 <X< 33
             ! 0.5Re:  |Z|<18, |Y|<18, -3<X<15
             ! 0.5Re:  |Z|<12, |Y|<12, -9<X<21
             ! 0.25Re: |Z|<9, |Y|<9, -6<X<15
             ! 0.25Re: |Z|<12, |Y|<12, 0<X<12
             ! 0.25Re: |Z|<15, |Y|<15, -6<X<9

             ! Near Earth
             ! 0.25Re: |Z|<6, |Y|<6, |X|<6

          end if
       end if

    case ('ror6_mp4')
       ! 6 x 6 x 6 block for High Resolution Runs On Request (magnetopause/shock/cusp)
       ! Box sizes: -384<X<192, -192<Z<192, -192<Y<192
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.

             ! add cylinder in night side and half sphere in dayside with radius=50 
             !                 if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -30. .and. &
             !                    (yyy*yyy+zzz*zzz)<2500.) DoRefine = .true.

             ! add elliptic area in front of -30 R_E up to 45 R_E in X and with
             ! 170 R_E radous in in Y-Z plane
             if (dx_BLK(iBLK)>1. .and. xxx > -40. .and. &
                  ((xxx+40.)*(xxx+40.)+(yyy*yyy+zzz*zzz)/4.)<7225.) DoRefine = .true.



             !                if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
             !                   abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             !                if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
             !                    abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             !                if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -99.0 .and. &
             !                    abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx> 0.0 .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -6.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             ! This will make it:
             !
             ! Tail
             ! 4 Re:  X>-159
             ! 2 Re:  X>-135
             ! 1 Re:  |Z|<12, |Y|<12, -111<X<0
             ! 0.5Re: |Z|<6, |Y|<6, -99<X<0

             ! Magnetopause/Shock
             ! 1 Re:   |Z|<24, |Y|<24, -15 <X< 33
             ! 0.5Re:  |Z|<18, |Y|<18, -3<X<15
             ! 0.5Re:  |Z|<12, |Y|<12, -9<X<21
             ! 0.25Re: |Z|<9, |Y|<9, -6<X<15
             ! 0.25Re: |Z|<12, |Y|<12, 0<X<12
             ! 0.25Re: |Z|<15, |Y|<15, -6<X<9

             ! Near Earth
             ! 0.25Re: |Z|<6, |Y|<6, |X|<6

          end if
       end if

    case ('mp4_hr')
       ! 6 x 6 x 6 block for High Resolution Runs On Request (magnetopause/shock/cusp)
       ! Box sizes: -384<X<192, -192<Z<192, -192<Y<192
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -24.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -120.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -24. .and. xxx < 48 .and. &
                  abs(yyy)<48..and. abs(zzz)<48.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -96. .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx > 0. .and. xxx < 36 .and. &
                  abs(yyy)<30. .and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 1. .and. xxx > -9. .and. &
                  abs(yyy)<24. .and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -30.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.


             ! add elliptic area in front of -30 R_E up to 45 R_E in X and with
             ! 170 R_E radous in in Y-Z plane
             if (dx_BLK(iBLK)>1. .and. xxx > -24. .and. &
                  ((xxx+40.)*(xxx+40.)+(yyy*yyy+zzz*zzz)/4.)<7225.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx > -3. .and. xxx < 18 .and. &
                  abs(yyy)<21..and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -6.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -12.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.


          end if
       end if

    case ('mp5_hr')
       ! 6 x 6 x 6 block for High Resolution Runs On Request (magnetopause/shock/cusp)
       ! Box sizes: -384<X<192, -192<Z<192, -192<Y<192
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -24.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -120.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -24. .and. xxx < 48 .and. &
                  abs(yyy)<48..and. abs(zzz)<48.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -96. .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx > 0. .and. xxx < 36 .and. &
                  abs(yyy)<30. .and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 1. .and. xxx > -9. .and. &
                  abs(yyy)<24. .and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -30.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.


             ! add elliptic area in front of -30 R_E up to 45 R_E in X and with
             ! 170 R_E radous in in Y-Z plane
             if (dx_BLK(iBLK)>1. .and. xxx > -24. .and. &
                  ((xxx+40.)*(xxx+40.)+(yyy*yyy+zzz*zzz)/4.)<7225.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx > -3. .and. xxx < 27. .and. &
                  abs(yyy)<21..and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -6.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -12.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.


          end if
       end if


    case ('ror6_hr')
       ! 6 x 6 x 6 block for High resolution Runs On Request 
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -99.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 51.0 .and. &
                  abs(yyy)<9..and. abs(zzz)<3.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx> 0.0 .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and.  abs(xxx) < 4.5 .and. &
                  abs(yyy) < 4.5 .and. abs(zzz)< 4.5) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and.  abs(xxx) < 3. .and. &
                  abs(yyy) < 3. .and. abs(zzz)< 3.) DoRefine = .true.

             ! This will make it:
             !
             ! Tail
             ! 4 Re:  X>-159
             ! 2 Re:  X>-135
             ! 1 Re:  |Z|<12, |Y|<24, -111<X<0
             ! 0.5Re: |Z|<6, |Y|<12, -99<X<0
             ! 0.25   |Z|<3, |Y|<9, -51<X<0

             ! Magnetopause/Shock
             ! 1 Re:   |Z|<24, |Y|<24, -15 <X< 33
             ! 0.5Re:  |Z|<18, |Y|<18, -3<X<15
             ! 0.5Re:  |Z|<12, |Y|<12, -9<X<21
             ! 0.25Re: |Z|<9, |Y|<9, -6<X<15
             ! 0.25Re: |Z|<12, |Y|<12, 0<X<12
             ! 0.25Re: |Z|<15, |Y|<15, -3<X<9

             ! Near Earth
             ! 0.25Re: |Z|<6, |Y|<6, |X|<6
             ! 0.125Re: |Z|<4.5, |Y|< 4.56, |X|< 4.5

          end if
       end if

    case ('AS1')
       ! 6 x 6 x 6 block 
       ! Box sizes: -267<X<21, -48<Z<48, -48<Y<48
       !
       ! Customized grid requested by  Andrey Samsonov
       ! Research: Interplanetary shock interaction with magnetosheath. Inflow boundary moved close
       ! to allow better resolution of the interplanetary shock.
       ! Resolution in magnetosheath up to 1/8 Re

       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -171.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -51.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -147.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -39. .and. xxx < 21.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -75. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21.0 .and. xxx > -21. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -51.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 36.0 .and. &
                  abs(yyy)<9..and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 24.0 .and. &
                  abs(yyy)<15..and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and.  xxx < 21.  .and.  xxx > -12. .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.125 .and.  abs(xxx) < 6. .and. &
                  abs(yyy) < 6. .and. abs(zzz)< 6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and.  xxx < 21. .and. xxx > -6. .and. &
                  abs(yyy) < 8. .and. abs(zzz)< 8.) DoRefine = .true.

          end if
       end if

    case ('bl1')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -351<X<33, -192<Z<192, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else
             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             ! add elliptic area in front of -30 R_E up to 45 R_E in X and with
             ! 170 R_E radous in in Y-Z plane
             if (dx_BLK(iBLK)>1. .and. xxx > -24. .and. &
                  ((xxx+40.)*(xxx+40.)+(yyy*yyy+zzz*zzz)/4.)<7225.) DoRefine = .true.



             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -99.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 0. .and. xxx> -33.0 .and. &
                  abs(yyy)<18..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 51.0 .and. &
                  abs(yyy)<9..and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 24.0 .and. &
                  abs(yyy)<15..and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx> 0.0 .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             !                  if (dx_BLK(iBLK)>0.125 .and.  abs(xxx) < 4.5 .and. &
             !                    abs(yyy) < 4.5 .and. abs(zzz)< 4.5) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and.  abs(xxx) < 6. .and. &
                  abs(yyy) < 6. .and. abs(zzz)< 6.) DoRefine = .true.


             !                if (dx_BLK(iBLK)>0.0625 .and.  abs(xxx) < 6. .and. &
             !                    abs(yyy) < 6. .and. abs(zzz)< 6.) DoRefine = .true.

             ! This will make it:
             !
             ! Tail
             ! 4 Re:  X>-159
             ! 2 Re:  X>-135
             ! 1 Re:  |Z|<12, |Y|<24, -111<X<0
             ! 0.5Re: |Z|<6, |Y|<12, -99<X<0
             ! 0.25   |Z|<3, |Y|<9, -51<X<0

             ! Magnetopause/Shock
             ! 1 Re:   |Z|<24, |Y|<24, -15 <X< 33
             ! 0.5Re:  |Z|<18, |Y|<18, -3<X<15
             ! 0.5Re:  |Z|<12, |Y|<12, -9<X<21
             ! 0.25Re: |Z|<9, |Y|<9, -6<X<15
             ! 0.25Re: |Z|<12, |Y|<12, 0<X<12
             ! 0.25Re: |Z|<15, |Y|<15, -3<X<9

             ! Near Earth
             ! 0.25Re: |Z|<6, |Y|<6, |X|<6
             ! 0.125Re: |Z|<4.5, |Y|< 4.56, |X|< 4.5

          end if
       end if

    case ('ror6_nykyri')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -99.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 0. .and. xxx> -33.0 .and. &
                  abs(yyy)<18..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 51.0 .and. &
                  abs(yyy)<9..and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 24.0 .and. &
                  abs(yyy)<15..and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx> 0.0 .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

! Katarina Nykyri

           if (dx_BLK(iBLK)>1. .and. xxx > -27. .and. xxx < 33.0 .and. &
                  abs(yyy)<30. .and. abs(zzz)<30.) DoRefine = .true.

          
           if (dx_BLK(iBLK)>0.5  .and. xxx > -21. .and. xxx < 12.0 .and. &
                  abs(yyy)<24. .and. abs(zzz)<24.) DoRefine = .true.

           if (dx_BLK(iBLK)>0.25 .and. xxx > -18. .and. xxx < 9.0 .and. &
                  yyy > -6. .and. yyy<15. .and.  zzz>3. .and. zzz < 18.) DoRefine = .true.

          if (dx_BLK(iBLK)>0.125 .and. xxx > -14. .and. xxx < 5.0 .and. &
                  yyy > -1. .and. yyy<10. .and.  zzz>7. .and. zzz < 16.) DoRefine = .true.


             !                  if (dx_BLK(iBLK)>0.125 .and.  abs(xxx) < 4.5 .and. &
             !                    abs(yyy) < 4.5 .and. abs(zzz)< 4.5) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and.  abs(xxx) < 6. .and. &
                  abs(yyy) < 6. .and. abs(zzz)< 6.) DoRefine = .true.


             !                if (dx_BLK(iBLK)>0.0625 .and.  abs(xxx) < 6. .and. &
             !                    abs(yyy) < 6. .and. abs(zzz)< 6.) DoRefine = .true.

             ! This will make it:
             !
             ! Tail
             ! 4 Re:  X>-159
             ! 2 Re:  X>-135
             ! 1 Re:  |Z|<12, |Y|<24, -111<X<0
             ! 0.5Re: |Z|<6, |Y|<12, -99<X<0
             ! 0.25   |Z|<3, |Y|<9, -51<X<0

             ! Magnetopause/Shock
             ! 1 Re:   |Z|<24, |Y|<24, -15 <X< 33
             ! 0.5Re:  |Z|<18, |Y|<18, -3<X<15
             ! 0.5Re:  |Z|<12, |Y|<12, -9<X<21
             ! 0.25Re: |Z|<9, |Y|<9, -6<X<15
             ! 0.25Re: |Z|<12, |Y|<12, 0<X<12
             ! 0.25Re: |Z|<15, |Y|<15, -3<X<9

             ! Near Earth
             ! 0.25Re: |Z|<6, |Y|<6, |X|<6
             ! 0.125Re: |Z|<4.5, |Y|< 4.56, |X|< 4.5

          end if
       end if



    case ('ror6_hr_1')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -99.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 0. .and. xxx> -33.0 .and. &
                  abs(yyy)<18..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 51.0 .and. &
                  abs(yyy)<9..and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 24.0 .and. &
                  abs(yyy)<15..and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx> 0.0 .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             !                  if (dx_BLK(iBLK)>0.125 .and.  abs(xxx) < 4.5 .and. &
             !                    abs(yyy) < 4.5 .and. abs(zzz)< 4.5) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and.  abs(xxx) < 6. .and. &
                  abs(yyy) < 6. .and. abs(zzz)< 6.) DoRefine = .true.


             !                if (dx_BLK(iBLK)>0.0625 .and.  abs(xxx) < 6. .and. &
             !                    abs(yyy) < 6. .and. abs(zzz)< 6.) DoRefine = .true.

             ! This will make it:
             !
             ! Tail
             ! 4 Re:  X>-159
             ! 2 Re:  X>-135
             ! 1 Re:  |Z|<12, |Y|<24, -111<X<0
             ! 0.5Re: |Z|<6, |Y|<12, -99<X<0
             ! 0.25   |Z|<3, |Y|<9, -51<X<0

             ! Magnetopause/Shock
             ! 1 Re:   |Z|<24, |Y|<24, -15 <X< 33
             ! 0.5Re:  |Z|<18, |Y|<18, -3<X<15
             ! 0.5Re:  |Z|<12, |Y|<12, -9<X<21
             ! 0.25Re: |Z|<9, |Y|<9, -6<X<15
             ! 0.25Re: |Z|<12, |Y|<12, 0<X<12
             ! 0.25Re: |Z|<15, |Y|<15, -3<X<9

             ! Near Earth
             ! 0.25Re: |Z|<6, |Y|<6, |X|<6
             ! 0.125Re: |Z|<4.5, |Y|< 4.56, |X|< 4.5

          end if
       end if

    case ('ror6_hr_2')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -735<X<33, -192<Z<192, -192<Y<192
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -351.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -20. .and. &
                  abs(yyy)<96. .and. abs(zzz)<96. )  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -207.0 .and. xxx< -15. .and. &
                  abs(yyy)<144. .and. abs(zzz)<144.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -99.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 0. .and. xxx> -33.0 .and. &
                  abs(yyy)<18..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 51.0 .and. &
                  abs(yyy)<9..and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 24.0 .and. &
                  abs(yyy)<15..and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx> 0.0 .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.


             ! This will make it:
             !
             ! Tail
             ! 4 Re:  X>-351
             ! 2 Re:  X>-207
             ! 1 Re:  |Z|<12, |Y|<24, -111<X<0
             ! 0.5Re: |Z|<6, |Y|<12, -99<X<0
             ! 0.25   |Z|<3, |Y|<9, -51<X<0

             ! Magnetopause/Shock
             ! 1 Re:   |Z|<24, |Y|<24, -15 <X< 33
             ! 0.5Re:  |Z|<18, |Y|<18, -3<X<15
             ! 0.5Re:  |Z|<12, |Y|<12, -9<X<21
             ! 0.25Re: |Z|<9, |Y|<9, -6<X<15
             ! 0.25Re: |Z|<12, |Y|<12, 0<X<12
             ! 0.25Re: |Z|<15, |Y|<15, -3<X<9

             ! Near Earth
             ! 0.25Re: |Z|<6, |Y|<6, |X|<6
             ! 0.125Re: |Z|<4.5, |Y|< 4.56, |X|< 4.5

          end if
       end if

    case ('ror6_Lis_Ros')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<36..and. abs(zzz)<36.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21.0 .and. xxx > -9. .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -99.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 0. .and. xxx> -33.0 .and. &
                  abs(yyy)<18..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 51.0 .and. &
                  abs(yyy)<9..and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 24.0 .and. &
                  abs(yyy)<15..and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx> 0.0 .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -6.0 .and. &
                  abs(yyy)<21. .and. abs(zzz)<21.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             !                  if (dx_BLK(iBLK)>0.125 .and.  abs(xxx) < 4.5 .and. &
             !                    abs(yyy) < 4.5 .and. abs(zzz)< 4.5) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and.  abs(xxx) < 6. .and. &
                  abs(yyy) < 6. .and. abs(zzz)< 6.) DoRefine = .true.


             !                if (dx_BLK(iBLK)>0.0625 .and.  abs(xxx) < 6. .and. &
             !                    abs(yyy) < 6. .and. abs(zzz)< 6.) DoRefine = .true.


          end if
       end if



    case ('ESS-261lr')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -21. .and. xxx < 33.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<36.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21.0 .and. xxx > -15. .and. &
                  abs(yyy)<36. .and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -99.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 0. .and. xxx> -33.0 .and. &
                  abs(yyy)<18..and. abs(zzz)<6.) DoRefine = .true.

             !
             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.


          end if
       end if

    case ('ESS-261hr')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -21. .and. xxx < 33.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<36.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21.0 .and. xxx > -15. .and. &
                  abs(yyy)<36. .and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -99.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 0. .and. xxx> -33.0 .and. &
                  abs(yyy)<18..and. abs(zzz)<6.) DoRefine = .true.

             !
!!! Uncomment for high resolution run

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 51.0 .and. &
                  abs(yyy)<9..and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 24.0 .and. &
                  abs(yyy)<15..and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. abs(xxx)  <  12. .and.  &
                  abs(yyy)< 21. .and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 18. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -3.0 .and. &
                  abs(yyy)<21. .and. abs(zzz)<15.) DoRefine = .true.

!!! End of high resolution run

             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and.  xxx > -9. .and.  xxx < 3. .and. &
                  abs(yyy) < 18 .and. abs(zzz)< 6) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and.  xxx > 2. .and.  xxx < 8. .and. &
                  abs(yyy) < 15 .and. abs(zzz)< 6) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and.  xxx > 7.5 .and.  xxx < 12. .and. &
                  abs(yyy) < 12 .and. abs(zzz)< 6) DoRefine = .true.


          end if
       end if

    case ('as1')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.
             !    1 Re
             if (dx_BLK(iBLK)>1. .and. xxx > -27. .and. xxx < 33.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<36.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -51. .and. xxx < -27 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < -51. .and. xxx > -75. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.
             !    0.5 Re
             if (dx_BLK(iBLK)>0.5 .and. xxx < 21.0 .and. xxx > -27. .and. &
                  abs(yyy)<30. .and. abs(zzz)<21.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < -27.0 .and. xxx > - 45. .and. &
                  abs(yyy)<30. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < - 45. .and. xxx> - 63.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             !    0.25 Re

             if (dx_BLK(iBLK)>0.25 .and. xxx < 18. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> 9.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                  abs(yyy)<21. .and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -36.0 .and. &
                  abs(yyy)<21. .and. abs(zzz)<18.) DoRefine = .true.

!!!
!!!

             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.


             !  0.125 Re

!!!                   if (dx_BLK(iBLK)>0.125 .and.  xxx > -24. .and.  xxx < -7.5 .and. &
!!!                  abs(yyy) < 19.5 .and. abs(yyy) > 13.5  .and. &
!!!                   abs(zzz)< 6.) DoRefine = .true.
!!!
!!!                   if (dx_BLK(iBLK)>0.125 .and.  xxx > -7.5 .and.  xxx < -2.5 .and. &
!!!                  abs(yyy) < 19.5 .and. abs(yyy) > 12.  .and. &
!!!                   abs(zzz)< 6) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and.  xxx > 6. .and.  xxx < 11 .and. &
                  abs(yyy) < 9.  .and. &
                  abs(zzz)< 10.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.125 .and.  xxx > 4.5 .and.  xxx < 12 .and. &
                  abs(yyy) < 6.  .and. &
                  abs(zzz)< 6.) DoRefine = .true.


             RRR=sqrt((xxx+3.)**2+1.3*zzz**2+(5.+0.77*abs(yyy))**2-25.)

             if (dx_BLK(iBLK)>0.125 .and.  xxx > -10. .and.  xxx < 14. .and. &
                  RRR >10. .and.  RRR < 17. .and. &
                  abs(zzz)< 15. .and. abs(yyy) < 18.) DoRefine = .true.
             !
             !
!!!! corrections on 04/14/07

             if (dx_BLK(iBLK)>0.125 .and.  xxx > -10. .and.  xxx < 14. .and. &
                  RRR >10. .and.  RRR < 20. .and. &
                  abs(zzz)< 6. .and. abs(yyy) < 18.) DoRefine = .true.
!!!!!         end corrections 04/14/07
             if (dx_BLK(iBLK)>0.125 .and. xxx > -10.5 .and.  xxx < 7.5 .and.  abs(yyy) < 12. .and. &
                  abs(zzz) < 10.5) DoRefine = .true.
!!!!!  correction on northern hemisphere
             if (dx_BLK(iBLK)>0.125 .and. xxx > -10.5 .and.  xxx < 10.5  .and.  yyy < 18. .and. &
                  yyy > 0. .and. zzz > 0. .and. zzz  < 15.) DoRefine = .true.


          end if
       end if


    case ('jb3')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.
             !    1 Re
             if (dx_BLK(iBLK)>1. .and. xxx > -27. .and. xxx < 33.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<36.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -51. .and. xxx < -27 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < -51. .and. xxx > -75. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.
             !    0.5 Re
             if (dx_BLK(iBLK)>0.5 .and. xxx < 21.0 .and. xxx > -27. .and. &
                  abs(yyy)<30. .and. abs(zzz)<21.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < -27.0 .and. xxx > - 45. .and. &
                  abs(yyy)<30. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < - 45. .and. xxx> - 63.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             !    0.25 Re
             if (dx_BLK(iBLK)>0.25 .and. xxx < 18. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> 9.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                  abs(yyy)<21. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -36.0 .and. &
                  abs(yyy)<21. .and. abs(zzz)<12.) DoRefine = .true.

!!!

             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.


             !  0.125 Re

!!!                   if (dx_BLK(iBLK)>0.125 .and.  xxx > -24. .and.  xxx < -7.5 .and. &
!!!                  abs(yyy) < 19.5 .and. abs(yyy) > 13.5  .and. &
!!!                   abs(zzz)< 6.) DoRefine = .true.
!!!
!!!                   if (dx_BLK(iBLK)>0.125 .and.  xxx > -7.5 .and.  xxx < -2.5 .and. &
!!!                  abs(yyy) < 19.5 .and. abs(yyy) > 12.  .and. &
!!!                   abs(zzz)< 6) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and.  xxx > 6. .and.  xxx < 11 .and. &
                  abs(yyy) < 9.  .and. &
                  abs(zzz)< 10.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.125 .and.  xxx > 4.5 .and.  xxx < 12 .and. &
                  abs(yyy) < 6.  .and. &
                  abs(zzz)< 6.) DoRefine = .true.


             RRR=sqrt((xxx+3.)**2+1.3*zzz**2+(5.+0.77*abs(yyy))**2-25.)

             if (dx_BLK(iBLK)>0.125 .and.  xxx > -3. .and.  xxx < 14. .and. &
                  RRR >10. .and.  RRR < 17. .and. &
                  abs(zzz)< 6.) DoRefine = .true.
             !
             !
!!!! corrections on 04/14/07

             if (dx_BLK(iBLK)>0.125 .and.  xxx > 1. .and.  xxx < 14. .and. &
                  RRR >10. .and.  RRR < 20. .and. &
                  abs(zzz)< 6.) DoRefine = .true.
!!!!!         end corrections 04/14/07
             ! 0.0625

             !!old                  if (dx_BLK(iBLK)>0.0625 .and.  xxx > -18. .and.  xxx < -7.5 .and. &
             !!old                  abs(yyy) < 18. .and. abs(yyy) > 15.  .and. &
             !!old                   abs(zzz)< 3.) DoRefine = .true.

             !!                   if (dx_BLK(iBLK)>0.0625 .and.  xxx > -7.5 .and.  xxx < -0.75 .and. &
             !!                  abs(yyy) < 18. .and. abs(yyy) > 13.5  .and. &
             !!                   abs(zzz)< 3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and.  xxx > 6. .and.  xxx < 12. .and. &
                  abs(yyy) < 6.  .and. &
                  abs(zzz)< 4.5) DoRefine = .true.


             if (dx_BLK(iBLK)>0.0625 .and.  xxx > 0 .and.  xxx < 12. .and. &
                  RRR > 11. .and. RRR < 15.5 .and. &
                  abs(zzz)< 4.5) DoRefine = .true.
!!! corrections 04/14/07
             if (dx_BLK(iBLK)>0.0625 .and.  xxx > 1. .and.  xxx < 14. .and. &
                  RRR > 11. .and. RRR < 18.5 .and. &
                  abs(zzz)< 4.5) DoRefine = .true.
!!! end 04/14/07
!!!
!!! corrections 04/15/07
             if (dx_BLK(iBLK)>0.0625 .and.  xxx > 6. .and.  xxx < 12. .and. &
                  abs(yyy) < 7.5 .and. &
                  abs(zzz)< 4.5) DoRefine = .true.
!!! end 04/15/07

          end if
       end if




    case ('jb1')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.
             !    1 Re
             if (dx_BLK(iBLK)>1. .and. xxx > -27. .and. xxx < 33.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<36.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -51. .and. xxx < -27 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < -51. .and. xxx > -75. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.
             !    0.5 Re
             if (dx_BLK(iBLK)>0.5 .and. xxx < 21.0 .and. xxx > -27. .and. &
                  abs(yyy)<30. .and. abs(zzz)<21.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < -27.0 .and. xxx > - 45. .and. &
                  abs(yyy)<30. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < - 45. .and. xxx> - 63.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             !    0.25 Re

             if (dx_BLK(iBLK)>0.25 .and. xxx < 18. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> 9.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                  abs(yyy)<21. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -36.0 .and. &
                  abs(yyy)<21. .and. abs(zzz)<12.) DoRefine = .true.

!!!

             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.


             !  0.125 Re

!!!                   if (dx_BLK(iBLK)>0.125 .and.  xxx > -24. .and.  xxx < -7.5 .and. &
!!!                  abs(yyy) < 19.5 .and. abs(yyy) > 13.5  .and. &
!!!                   abs(zzz)< 6.) DoRefine = .true.
!!!
!!!                   if (dx_BLK(iBLK)>0.125 .and.  xxx > -7.5 .and.  xxx < -2.5 .and. &
!!!                  abs(yyy) < 19.5 .and. abs(yyy) > 12.  .and. &
!!!                   abs(zzz)< 6) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and.  xxx > 6. .and.  xxx < 11 .and. &
                  abs(yyy) < 9.  .and. &
                  abs(zzz)< 10.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.125 .and.  xxx > 7.5 .and.  xxx < 12 .and. &
                  abs(yyy) < 6.  .and. &
                  abs(zzz)< 6.) DoRefine = .true.


             RRR=sqrt((xxx+3.)**2+1.3*zzz**2+(5.+0.77*abs(yyy))**2-25.)

             if (dx_BLK(iBLK)>0.125 .and.  xxx > -3. .and.  xxx < 14. .and. &
                  RRR >10. .and.  RRR < 17. .and. &
                  abs(zzz)< 6.) DoRefine = .true.

             ! 0.0625

             !!old                  if (dx_BLK(iBLK)>0.0625 .and.  xxx > -18. .and.  xxx < -7.5 .and. &
             !!old                  abs(yyy) < 18. .and. abs(yyy) > 15.  .and. &
             !!old                   abs(zzz)< 3.) DoRefine = .true.

             !!                   if (dx_BLK(iBLK)>0.0625 .and.  xxx > -7.5 .and.  xxx < -0.75 .and. &
             !!                  abs(yyy) < 18. .and. abs(yyy) > 13.5  .and. &
             !!                   abs(zzz)< 3.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.0625 .and.  xxx > 0 .and.  xxx < 13. .and. &
                  RRR > 11. .and. RRR < 15.5 .and. &
                  abs(zzz)< 4.5) DoRefine = .true.

          end if
       end if

    case ('ESS261hr2')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.
             !    1 Re
             if (dx_BLK(iBLK)>1. .and. xxx > -27. .and. xxx < 33.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<36.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -51. .and. xxx < -27 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < -51. .and. xxx > -75. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.
             !    0.5 Re
             if (dx_BLK(iBLK)>0.5 .and. xxx < 21.0 .and. xxx > -27. .and. &
                  abs(yyy)<30. .and. abs(zzz)<21.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < -27.0 .and. xxx > - 45. .and. &
                  abs(yyy)<30. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < - 45. .and. xxx> - 63.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.


             !    0.25 Re

             if (dx_BLK(iBLK)>0.25 .and. xxx < 18. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> 9.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                  abs(yyy)<21. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -36.0 .and. &
                  abs(yyy)<21. .and. abs(zzz)<12.) DoRefine = .true.

!!! 

             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             !  0.125 Re

             if (dx_BLK(iBLK)>0.125 .and.  xxx > -24. .and.  xxx < -7.5 .and. &
                  abs(yyy) < 19.5 .and. abs(yyy) > 13.5  .and. &
                  abs(zzz)< 6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and.  xxx > -7.5 .and.  xxx < -2.5 .and. &
                  abs(yyy) < 19.5 .and. abs(yyy) > 12.  .and. &
                  abs(zzz)< 6) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and.  xxx > 6. .and.  xxx < 11 .and. &
                  abs(yyy) < 9.  .and. &
                  abs(zzz) > 5. .and.  abs(zzz)< 10.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.125 .and.  xxx > 7.5 .and.  xxx < 12 .and. &
                  abs(yyy) < 6.  .and. &
                  abs(zzz)< 6.) DoRefine = .true.


             RRR=sqrt((xxx+3.)**2+1.3*zzz**2+(5.+0.77*abs(yyy))**2-25.)

             if (dx_BLK(iBLK)>0.125 .and.  xxx > -3. .and.  xxx < 16. .and. &
                  RRR >12. .and.  RRR < 18. .and. &
                  abs(zzz)< 6.) DoRefine = .true.



             ! 0.0625

             !!old                  if (dx_BLK(iBLK)>0.0625 .and.  xxx > -18. .and.  xxx < -7.5 .and. &
             !!old                  abs(yyy) < 18. .and. abs(yyy) > 15.  .and. &
             !!old                   abs(zzz)< 3.) DoRefine = .true.

             !!                   if (dx_BLK(iBLK)>0.0625 .and.  xxx > -7.5 .and.  xxx < -0.75 .and. &
             !!                  abs(yyy) < 18. .and. abs(yyy) > 13.5  .and. &
             !!                   abs(zzz)< 3.) DoRefine = .true.


             !!                  if (dx_BLK(iBLK)>0.0625 .and.  xxx > -3. .and.  xxx < 16. .and. &
             !!                  RRR > 12.5 .and. RRR < 16.5 .and. &
             !!                   abs(zzz)< 3.) DoRefine = .true.

          end if
       end if

             case ('ror6_small_new6')
           ! 6 x 6 x 6 block for High resolution Runs On Request
           ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
           if (maxRblk > Rbody) then

              if (dx_BLK(iBLK) > 8.) then
                 ! Refine all blocks with dx greater than 8
                 DoRefine  = .true.
              else

                 if (dx_BLK(iBLK) > 4. .and. xxx > -159.0 .and. &
                 abs(yyy)<96..and. abs(zzz)<96.)  DoRefine  = .true.

                 if (dx_BLK(iBLK) > 2. .and. xxx > -39.0 .and. &
              abs(yyy)<96..and. abs(zzz)<96.)  DoRefine  = .true.

                 if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                   abs(yyy)<48..and. abs(zzz)<48.)  DoRefine  = .true.

                 if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                    abs(yyy)<24..and. abs(zzz)<24.) DoRefine  = .true.

                 if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -123. .and. &
                    abs(yyy)<36..and. abs(zzz)<24.) DoRefine  = .true.

                  if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
                    abs(yyy)<18..and. abs(zzz)<18.) DoRefine  = .true.

                  if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                    abs(yyy)<12..and. abs(zzz)<12.) DoRefine  = .true.

                  if (dx_BLK(iBLK)>0.5 .and. xxx < 18. .and. xxx> -111.0 .and. &
                    abs(yyy)<30..and. abs(zzz)<18.) DoRefine  = .true.

                  if (dx_BLK(iBLK)>0.25 .and. xxx < 6. .and. xxx> -6. .and. &
                    abs(yyy)<8. .and. abs(zzz)<8.) DoRefine  = .true.

              end if
           end if

  case ('ror6_tail_new6')
           ! 6 x 6 x 6 block for High resolution Runs On Request
           ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
           if (maxRblk > Rbody) then

              if (dx_BLK(iBLK) > 8.) then
                 ! Refine all blocks with dx greater than 8
                 DoRefine  = .true.
              else

                 if (dx_BLK(iBLK) > 4. .and. xxx > -159.0 .and. &
                 abs(yyy)<96..and. abs(zzz)<96.)  DoRefine  = .true.

                 if (dx_BLK(iBLK) > 2. .and. xxx > -39.0 .and. &
              abs(yyy)<96..and. abs(zzz)<96.)  DoRefine  = .true.

                 if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                   abs(yyy)<48..and. abs(zzz)<48.)  DoRefine  = .true.

                 if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                    abs(yyy)<24..and. abs(zzz)<24.) DoRefine  = .true.

                 if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -123. .and. &
                    abs(yyy)<36..and. abs(zzz)<24.) DoRefine  = .true.

                  if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
                    abs(yyy)<18..and. abs(zzz)<18.) DoRefine  = .true.

                  if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                    abs(yyy)<12..and. abs(zzz)<12.) DoRefine  = .true.

                  if (dx_BLK(iBLK)>0.5 .and. xxx < 18. .and. xxx> -111.0 .and. &
                    abs(yyy)<30..and. abs(zzz)<18.) DoRefine  = .true.

                  if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. xxx> -87. .and. &
                    abs(yyy)<27. .and. abs(zzz)<6.) DoRefine  = .true.

                  if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx> -3.0 .and. &
                    abs(yyy)<15. .and. abs(zzz)<15.) DoRefine  = .true.

                 if (dx_BLK(iBLK)>0.25 .and. xxx < 18. .and. xxx> -6.0 .and. &
                    abs(yyy)<12. .and. abs(zzz)<9.) DoRefine  = .true.

                    if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. xxx> -26.0 .and. &
                    abs(yyy)<21. .and. abs(zzz)<9.) DoRefine  = .true.

!!                  if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 10.0 .and. &
!!                    abs(yyy)<12..and. abs(zzz)<12.) DoRefine  = .true.

                  if (dx_BLK(iBLK)>0.125 .and.  abs(xxx) < 12. .and. &
                    abs(yyy) < 8. .and. abs(zzz)< 6. ) DoRefine  = .true.


                    if (dx_BLK(iBLK)>0.125 .and. xxx < 0. .and. &
                     xxx > - 24.0 .and. &
                    abs(yyy)< 18. .and. abs(zzz)<6.) DoRefine  = .true.

                    if (dx_BLK(iBLK)>0.125 .and. xxx < 8. .and. &
                     xxx > - 6.0 .and. &
                    abs(yyy)< 12. .and. abs(zzz)<6.) DoRefine  = .true.

!!
                    if (dx_BLK(iBLK)>0.125 .and. xxx < -6. .and. &
                     xxx > - 48.0 .and. &
                    abs(yyy)< 21. .and. abs(zzz)<3.) DoRefine  = .true.

                    if (dx_BLK(iBLK)>0.125 .and. xxx < -40. .and. &
                     xxx > - 60.0 .and. &
                    abs(yyy)< 18. .and. abs(zzz)<3.) DoRefine  = .true.

!!                    if (dx_BLK(iBLK)>0.125 .and. xxx < -70. .and. &
!!                     xxx > - 87.0 .and. &
!!                    abs(yyy)< 15. .and. abs(zzz)<3.) DoRefine  = .true.


!!                  if (dx_BLK(iBLK)>0.0625 .and. xxx < -8. .and. &
!!                     xxx > - 20.0 .and. &
!!                    abs(yyy)<12. .and. abs(zzz)<3.) DoRefine  = .true.
!!
!!                      if (dx_BLK(iBLK)>0.0625 .and. xxx < -8. .and. &
!!                     xxx > - 50.0 .and. &
!!                    abs(yyy)<9. .and. abs(zzz)<1.5) DoRefine  = .true.
!!
!!
!!                  if (dx_BLK(iBLK)>0.03125 .and. xxx < -10. .and. &
!!                     xxx > - 16.0 .and. &
!!                    abs(yyy)<9. .and. abs(zzz)<0.75) DoRefine  = .true.
              end if
           end if





    case ('ror6_sawt_jb')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -63.  .and. &
                  abs(yyy)<36..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 9. .and. xxx> -51.0 .and. &
                  abs(yyy)<30..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -99.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 0. .and. xxx> -33.0 .and. &
                  abs(yyy)<18..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. xxx> -39. .and. &
                  abs(yyy)<24. .and. abs(zzz)<6.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 51.0 .and. &
                  abs(yyy)<9..and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 24.0 .and. &
                  abs(yyy)<15..and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx> 0.0 .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and.  abs(xxx) < 4.5 .and. &
                  abs(yyy) < 4.5 .and. abs(zzz)< 4.5) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and.  abs(xxx) < 3. .and. &
                  abs(yyy) < 3. .and. abs(zzz)< 3.) DoRefine = .true.

             ! This will make it:
             !
             ! Tail
             ! 4 Re:  X>-159
             ! 2 Re:  X>-135
             ! 1 Re:  |Z|<12, |Y|<24, -111<X<0
             ! 0.5Re: |Z|<6, |Y|<12, -99<X<0
             ! 0.25   |Z|<3, |Y|<9, -51<X<0

             ! Magnetopause/Shock
             ! 1 Re:   |Z|<24, |Y|<24, -15 <X< 33
             ! 0.5Re:  |Z|<18, |Y|<18, -3<X<15
             ! 0.5Re:  |Z|<12, |Y|<12, -9<X<21
             ! 0.25Re: |Z|<9, |Y|<9, -6<X<15
             ! 0.25Re: |Z|<12, |Y|<12, 0<X<12
             ! 0.25Re: |Z|<15, |Y|<15, -3<X<9

             ! Near Earth
             ! 0.25Re: |Z|<6, |Y|<6, |X|<6
             ! 0.125Re: |Z|<4.5, |Y|< 4.56, |X|< 4.5

          end if
       end if


    case ('ror6_cusp_pd')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21.0 .and. xxx > -9. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -99.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 0. .and. xxx> -33.0 .and. &
                  abs(yyy)<18..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 51.0 .and. &
                  abs(yyy)<9..and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 24.0 .and. &
                  abs(yyy)<15..and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 18. .and. xxx> -6.  .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.           

             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and.  abs(xxx) < 4.5 .and. &
                  abs(yyy) < 4.5 .and. abs(zzz)< 4.5) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and.  abs(xxx) < 3. .and. &
                  abs(yyy) < 3. .and. abs(zzz)< 3.) DoRefine = .true.! This will make it:
             !
             ! Tail
             ! 4 Re:  X>-159
             ! 2 Re:  X>-135
             ! 1 Re:  |Z|<12, |Y|<24, -111<X<0
             ! 0.5Re: |Z|<6, |Y|<12, -99<X<0
             ! 0.5Re: |Z|<6, |Y|<18, -33<X<0      
             ! 0.25   |Z|<3, |Y|<9, -51<X<0
             ! 0.25   |Z|<3, |Y|<15, -24<X<0

             ! Magnetopause/Shock
             ! 1 Re:   |Z|<24, |Y|<24, -15 <X< 33
             ! 0.5Re:  |Z|<18, |Y|<18, -9<X<121  
             ! 0.25Re: |Z|<12, |Y|<12, -6<X<18
             ! 0.25Re: |Z|<15, |Y|<15, -3<X<9

             ! Near Earth
             ! 0.25Re: |Z|<6, |Y|<6, |X|<6
             ! 0.125Re: |Z|<4.5, |Y|< 4.56, |X|< 4.5

          end if
       end if

    case ('magcon1')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<24.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -99.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 0 .and. xxx> -51.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0 .and. xxx> -33.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 51.0 .and. &
                  abs(yyy)<9..and. abs(zzz)<3.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx> 0.0 .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.



             ! This will make it:
             !
             ! Tail
             ! 4 Re:  X>-159
             ! 2 Re:  X>-135
             ! 1 Re:  |Z|<12, |Y|<24, -111<X<0
             ! 0.5Re: |Z|<6, |Y|<12, -99<X<0
             ! 0.5Re:  |Z|<12, |Y|<24, -51<X<0
             ! 0.25   |Z|<3, |Y|<9, -51<X<0
             ! 0.25   |Z|<6, |Y|<12, -33<X<0

             ! Magnetopause/Shock
             ! 1 Re:   |Z|<24, |Y|<24, -15 <X< 33
             ! 0.5Re:  |Z|<18, |Y|<18, -3<X<15
             ! 0.5Re:  |Z|<12, |Y|<12, -9<X<21
             ! 0.25Re: |Z|<9, |Y|<9, -6<X<15
             ! 0.25Re: |Z|<12, |Y|<12, 0<X<12
             ! 0.25Re: |Z|<15, |Y|<15, -3<X<9

             ! Near Earth
             ! 0.25Re: |Z|<6, |Y|<6, |X|<6

          end if
       end if


    case ('ror6_tail_new2')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<48.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<48..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15. .and. xxx> -87.0 .and. &
                  abs(yyy)<36..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. &
                  xxx > - 63.0 .and. &
                  abs(yyy)<24.and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. xxx < 0. .and. &
                  xxx > - 51.0 .and. &
                  abs(yyy)< 18. .and. abs(zzz)<1.5) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. xxx < -6. .and. &
                  xxx > - 36.0 .and. &
                  abs(yyy)<12. .and. abs(zzz)<0.75) DoRefine = .true.

             if (dx_BLK(iBLK)>0.03125 .and. xxx < -6. .and. &
                  xxx > - 26.0 .and. &
                  abs(yyy)<3.75 .and. abs(zzz)<0.375) DoRefine = .true.


             if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx> 0.0 .and. &
                  abs(yyy)<12. .and. abs(zzz)<12.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -3.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.
             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 6.0 .and. &
                  abs(yyy)<6..and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and.  abs(xxx) < 4.5 .and. &
                  abs(yyy) < 4.5 .and. abs(zzz)< 4.5) DoRefine = .true.


             ! This will make it:
             ! This will make it:
             !
             ! Tail
             ! 4 Re:  X>-159
             ! 2 Re:  X>-135
             ! 1 Re:  |Z|<24, |Y|<24, -111<X<0
             ! 0.5Re: |Z|<12, |Y|<18, -87<X<0
             ! 0.25   |Z|<6, |Y|<15, -69<X<0
             ! 0.125   |Z|<4.5, |Y|<12, -63<X<0
             ! 0.0625  |Z|<3.,  |Y|<10.5 -60<X< -4.5
             ! 0.03125  |Z|<1.5  |Y|< 9  -57<X< - 5

             ! Magnetopause/Shock
             ! 1 Re:   |Z|<24, |Y|<24, -15 <X< 33
             ! 0.5Re:  |Z|<18, |Y|<18, -3<X<15
             ! 0.5Re:  |Z|<12, |Y|<12, -9<X<21
             ! 0.25Re: |Z|<9, |Y|<9, -6<X<15
             ! 0.25Re: |Z|<12, |Y|<12, 0<X<12
             ! 0.25Re: |Z|<15, |Y|<15, -3<X<9

             ! Near Earth
             ! 0.25Re: |Z|<6, |Y|<6, |X|<6
             ! 0.125Re: |Z|<4.5, |Y|< 4.56, |X|< 4.5

          end if
       end if

    case ('ror6_tail_new5')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -351<X<33, -96<Z<96, -96<Y<96
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0 .and. &
                  abs(yyy)<96..and. abs(zzz)<96.)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0 .and. &
                  abs(yyy)<96..and. abs(zzz)<96.)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<48.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -123. .and. &
                  abs(yyy)<36..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 18. .and. xxx> -111.0 .and. &
                  abs(yyy)<30..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. xxx> -99. .and. &
                  abs(yyy)<27. .and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx> -3.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 18. .and. xxx> -6.0 .and. &
                  abs(yyy)<12. .and. abs(zzz)<9.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. xxx> -26.0 .and. &
                  abs(yyy)<21. .and. abs(zzz)<9.) DoRefine = .true.


             !!                  if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 12.0 .and. &
             !!                    abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and.  abs(xxx) < 10. .and. &
                  abs(yyy) < 8. .and. abs(zzz)< 6. ) DoRefine = .true.


             if (dx_BLK(iBLK)>0.125 .and. xxx < 0. .and. &
                  xxx > - 24.0 .and. &
                  abs(yyy)< 18. .and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. xxx < 8. .and. &
                  xxx > - 6.0 .and. &
                  abs(yyy)< 12. .and. abs(zzz)<6.) DoRefine = .true.

             !!
             if (dx_BLK(iBLK)>0.125 .and. xxx < -6. .and. &
                  xxx > - 48.0 .and. &
                  abs(yyy)< 21. .and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. xxx < -40. .and. &
                  xxx > - 72.0 .and. &
                  abs(yyy)< 18. .and. abs(zzz)<3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. xxx < -70. .and. &
                  xxx > - 87.0 .and. &
                  abs(yyy)< 15. .and. abs(zzz)<3.) DoRefine = .true.


             !!                  if (dx_BLK(iBLK)>0.0625 .and. xxx < -8. .and. &
             !!                     xxx > - 20.0 .and. &
             !!                    abs(yyy)<12. .and. abs(zzz)<3.) DoRefine = .true.
             !!
             !!                      if (dx_BLK(iBLK)>0.0625 .and. xxx < -8. .and. &
             !!                     xxx > - 50.0 .and. &
             !!                    abs(yyy)<9. .and. abs(zzz)<1.5) DoRefine = .true.
             !!
             !!
             !!                  if (dx_BLK(iBLK)>0.03125 .and. xxx < -10. .and. &
             !!                     xxx > - 16.0 .and. &
             !!                    abs(yyy)<9. .and. abs(zzz)<0.75) DoRefine = .true.
          end if
       end if

    case ('tail_winter5')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0 .and. &
                  abs(yyy)<96..and. abs(zzz)<96.)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0 .and. &
                  abs(yyy)<96..and. abs(zzz)<96.)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<48.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -123. .and. &
                  abs(yyy)<36..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 18. .and. xxx> -87.0 .and. &
                  abs(yyy)<30..and. abs(zzz)<18.) DoRefine = .true.

             !!   0.25

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. xxx> -75. .and. &
                  abs(yyy)<27. .and. zzz < 6. .and.  zzz > -9.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx> -3.0 .and. &
                  abs(yyy)<15. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 18. .and. xxx> -6.0 .and. &
                  abs(yyy)<12. .and. abs(zzz)<9.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. xxx> -51.0 .and. &
                  abs(yyy)<21. .and. zzz < 6. .and. zzz > -12.) DoRefine = .true.


             !!                  if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 12.0 .and. &
             !!  
             !!                    abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.

             !!   0.125

!!!!             if (dx_BLK(iBLK)>0.125 .and.  abs(xxx) < 10. .and. &
!!!!                  abs(yyy) < 8. .and. abs(zzz)< 6. ) DoRefine = .true.


!!!!             if (dx_BLK(iBLK)>0.125 .and. xxx < 8. .and. &
!!!!                  xxx > - 6.0 .and. &
!!!!                  abs(yyy)< 12. .and. abs(zzz)<6.) DoRefine = .true.

             !!
!!!!
!!!!             if (dx_BLK(iBLK)>0.125 .and. xxx < 0. .and. &
!!!!                  xxx > - 33.0 .and. &
!!!                  abs(yyy)< 18. .and. zzz < 3. .and. zzz > - 9.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. xxx < -6. .and. &
                  xxx > - 45.0 .and. &
                   abs(yyy)< 18. .and.  zzz < 3. .and. zzz > - 9.) DoRefine = .true.
!!!!                  abs(yyy)< 21. .and.  zzz < 3. .and. zzz > - 9.) DoRefine = .true.

!!!!             if (dx_BLK(iBLK)>0.125 .and. xxx < -40. .and. &
!!!!                  xxx > - 57.0 .and. &
!!!!                  abs(yyy)< 18. .and.  zzz < 3. .and. zzz > -9.) DoRefine = .true.




             !!                  if (dx_BLK(iBLK)>0.0625 .and. xxx < -8. .and. &
             !!                     xxx > - 24.0 .and. &
             !!                    abs(yyy)<12. .and. abs(zzz)<3.) DoRefine = .true.
             !!
             !!                      if (dx_BLK(iBLK)>0.0625 .and. xxx < -8. .and. &
             !!                     xxx > - 36.0 .and. &
             !!                    abs(yyy)<9. .and. abs(zzz)<1.5) DoRefine = .true.
             !!
             !!
             !!                  if (dx_BLK(iBLK)>0.03125 .and. xxx < -10. .and. &
             !!                     xxx > - 16.0 .and. &
             !!                    abs(yyy)<9. .and. abs(zzz)<0.75) DoRefine = .true.
          end if
       end if


    case ('ror6_mp_105')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<48.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -9. .and. xxx < 33.0 .and. &
                  abs(yyy)<36..and. abs(zzz)<36.) DoRefine = .true.


             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<36..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 18. .and. xxx> -87.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. xxx> -48. .and. &
                  abs(yyy)<18. .and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -3.0 .and. &
                  abs(yyy)<18. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 18. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 12.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.
             !
             ! IMF angle = 105
             RRR = sqrt(xxx*xxx*1.03*1.03 + 0.83*0.83*yyy*yyy + zzz*zzz)
             RRR1 = 0.5*xxx+sqrt(0.83*0.83*yyy*yyy + zzz*zzz)
             !          
             if (dx_BLK(iBLK)>0.125 .and. RRR1>12. .and. &
                  RRR1<18.  .and. xxx > -4.5 .and. xxx < 7.5 .and. &
                  yyy > 1.5 .and. yyy<21 .and. zzz > 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR1>12. .and. &
                  RRR1<18.  .and. xxx > -4.5 .and. xxx < 7.5 .and. &
                  yyy > -21. .and. yyy<-1.5 .and. zzz < 0.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.125 .and. RRR>12. .and. &
                  RRR<16.  .and. xxx > 7. .and. xxx < 10.5 .and. &
                  yyy > 0. .and. zzz > 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR>12. .and. &
                  RRR<16.  .and. xxx > 7. .and. xxx < 10.5 .and. &
                  yyy < 0. .and. zzz < 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR>12. .and. &
                  RRR<16.5  .and. xxx > 10.   &
                  ) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. RRR1>14 .and. &
                  RRR1<17.  .and. xxx > -3. .and. xxx < 6. .and. &
                  yyy > 2. .and. yyy<18. .and. zzz > 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. RRR1>14. .and. &
                  RRR1<17.  .and. xxx > -3. .and. xxx < 6. .and. &
                  yyy > -18. .and. yyy<-2. .and. zzz < 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. RRR>12.5 .and. &
                  RRR<15.  .and. xxx > 11.   &
                  ) DoRefine = .true.


             !
             !                  if (dx_BLK(iBLK)>0.125 .and.  abs(xxx) < 8. .and. &
!!!                    abs(yyy) < 8. .and. abs(zzz)< 8. ) DoRefine = .true.
!!!
!!!                   if (dx_BLK(iBLK)>0.125 .and.  xxx > 0. .and. xxx < 13.5 .and. &
!!!                    abs(yyy) < 12. .and. abs(zzz)< 12. ) DoRefine = .true.
!!!
!!!                    if (dx_BLK(iBLK)>0.125 .and.  xxx > 0. .and. xxx < 15. .and. &
!!!                    abs(yyy) < 6. .and. abs(zzz)< 6. ) DoRefine = .true.

             !!                    if (dx_BLK(iBLK)>0.125 .and. xxx < 0. .and. &
             !!                     xxx > - 24.0 .and. &
             !!                    abs(yyy)< 15. .and. abs(zzz)<3.) DoRefine = .true.

             !!                  if (dx_BLK(iBLK)>0.0625 .and. xxx < -6. .and. &
             !!                     xxx > - 18.0 .and. &
             !!                    abs(yyy)<12. .and. abs(zzz)<1.5) DoRefine = .true.
             !!
             !!                  if (dx_BLK(iBLK)>0.03125 .and. xxx < -10. .and. &
             !!                     xxx > - 16.0 .and. &
             !!                    abs(yyy)<9. .and. abs(zzz)<0.75) DoRefine = .true.
          end if
       end if
       !
    case ('ror6_mp_120')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<48.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -9. .and. xxx < 33.0 .and. &
                  abs(yyy)<36..and. abs(zzz)<36.) DoRefine = .true.


             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<36..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 18. .and. xxx> -87.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. xxx> -48. .and. &
                  abs(yyy)<18. .and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -3.0 .and. &
                  abs(yyy)<18. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 18. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 12.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.
             !
             ! IMF angle = 120
             !
             if (dx_BLK(iBLK)>0.25 .and. xxx < 0 .and. xxx > -6. .and. &
                  yyy > 0. .and. yyy < 21. .and. &
                  zzz > 0. .and. zzz< 9.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0 .and. xxx > -6. .and. &
                  yyy < 0. .and. yyy > -21. .and. &
                  zzz < 0. .and. zzz>  -9.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0 .and. xxx > -6. .and. &
                  yyy > 0. .and. yyy < 12. .and. &
                  zzz > 0. .and. zzz< 18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0 .and. xxx > -6. .and. &
                  yyy < 0. .and. yyy > -12. .and. &
                  zzz < 0. .and. zzz>  -18.) DoRefine = .true.

             RRR = sqrt(xxx*xxx*1.03*1.03 + 0.83*0.83*yyy*yyy + zzz*zzz)
             RRR1 = 0.5*xxx+sqrt(0.83*0.83*yyy*yyy + zzz*zzz)
             RRR2 = 0.5*xxx+sqrt(0.85*0.85*yyy*yyy + zzz*zzz)
             RRR3 = 0.6*xxx+sqrt(0.83*0.83*yyy*yyy + zzz*zzz)

             if (dx_BLK(iBLK)>0.125 .and. RRR1>12. .and. &
                  RRR1<18.  .and. xxx > -6. .and. xxx < 3. .and. &
                  yyy > 9.  .and. zzz > 3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR1>12. .and. &
                  RRR1<18.  .and. xxx > -6. .and. xxx < 3. .and. &
                  yyy<-9. .and. zzz < -3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR1>12. .and. &
                  RRR1<18.  .and. xxx > 3. .and. xxx < 6. .and. &
                  yyy > 1.5  .and. zzz > 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR1>12. .and. &
                  RRR1<18.  .and. xxx > 3. .and. xxx < 6. .and. &
                  yyy<-1.5 .and. zzz < 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR3>12. .and. &
                  RRR3<18.  .and. xxx > 6. .and. xxx < 9.  .and. &
                  yyy > 0.  .and. zzz > 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR3>12. .and. &
                  RRR3<18.  .and. xxx > 6. .and. xxx < 9. .and. &
                  yyy< 0. .and. zzz < 0.) DoRefine = .true.
             !

             if (dx_BLK(iBLK)>0.125 .and. RRR>12. .and. &
                  RRR<16.5  .and. xxx > 9.  .and. &
                  abs(zzz)<9. .and. abs(yyy)< 12. &
                  ) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR>12. .and. &
                  RRR<16.5  .and. xxx > 12.  .and. &
                  abs(zzz)<9. .and. abs(yyy)< 12. &
                  ) DoRefine = .true.

             !
             if (dx_BLK(iBLK)>0.0625 .and. RRR1>14. .and. &
                  RRR1<17.  .and. xxx > -3 .and. xxx < 3. .and. &
                  yyy > 7.5  .and. zzz > 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. RRR1>14. .and. &
                  RRR1<17.  .and. xxx > -3. .and. xxx < 3. .and. &
                  yyy<-7.5 .and. zzz < 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. RRR1>14. .and. &
                  RRR1<17.  .and. xxx > 3. .and. xxx < 6. .and. &
                  yyy > 3.  .and. zzz > 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. RRR1>14. .and. &
                  RRR1<17.  .and. xxx > 3. .and. xxx < 6. .and. &
                  yyy<-3. .and. zzz < 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. RRR3>14. .and. &
                  RRR3<17.  .and. xxx > 6. .and. xxx < 9.  .and. &
                  yyy > 0.  .and. zzz > 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. RRR3>14. .and. &
                  RRR3<17.  .and. xxx > 6. .and. xxx < 9. .and. &
                  yyy< 0. .and. zzz < 0.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.0625 .and. RRR>12.5 .and. &
                  RRR<15.  .and. xxx > 9.  .and. &
                  abs(zzz)<7.5 .and. abs(yyy)< 9. &
                  ) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. RRR>12.5 .and. &
                  RRR<16.  .and. xxx > 12.  .and. &
                  abs(zzz)<7.5 .and. abs(yyy)< 9. &
                  ) DoRefine = .true.

             !
             !
             !                  if (dx_BLK(iBLK)>0.125 .and.  abs(xxx) < 8. .and. &
!!!                    abs(yyy) < 8. .and. abs(zzz)< 8. ) DoRefine = .true.
!!!
!!!                   if (dx_BLK(iBLK)>0.125 .and.  xxx > 0. .and. xxx < 13.5 .and. &
!!!                    abs(yyy) < 12. .and. abs(zzz)< 12. ) DoRefine = .true.
!!!
!!!                    if (dx_BLK(iBLK)>0.125 .and.  xxx > 0. .and. xxx < 15. .and. &
!!!                    abs(yyy) < 6. .and. abs(zzz)< 6. ) DoRefine = .true.

             !!                    if (dx_BLK(iBLK)>0.125 .and. xxx < 0. .and. &
             !!                     xxx > - 24.0 .and. &
             !!                    abs(yyy)< 15. .and. abs(zzz)<3.) DoRefine = .true.

             !!                  if (dx_BLK(iBLK)>0.0625 .and. xxx < -6. .and. &
             !!                     xxx > - 18.0 .and. &
             !!                    abs(yyy)<12. .and. abs(zzz)<1.5) DoRefine = .true.
             !!
             !!                  if (dx_BLK(iBLK)>0.03125 .and. xxx < -10. .and. &
             !!                     xxx > - 16.0 .and. &
             !!                    abs(yyy)<9. .and. abs(zzz)<0.75) DoRefine = .true.
          end if
       end if
       !
    case ('ror6_mp_225')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<48.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -9. .and. xxx < 33.0 .and. &
                  abs(yyy)<36..and. abs(zzz)<36.) DoRefine = .true.


             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<36..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 18. .and. xxx> -87.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. xxx> -48. .and. &
                  abs(yyy)<18. .and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -3.0 .and. &
                  abs(yyy)<18. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 18. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 12.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.
             !
             ! IMF angle = 225
             if (dx_BLK(iBLK)>0.25 .and. xxx < 0 .and. xxx > -6. .and. &
                  yyy < 0. .and. yyy > -21. .and. &
                  zzz > 0. .and. zzz< 9.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0 .and. xxx > -6. .and. &
                  yyy > 0. .and. yyy < 21. .and. &
                  zzz < 0. .and. zzz>  -9.) DoRefine = .true.

             RRR = sqrt(xxx*xxx*1.03*1.03 + 0.83*0.83*yyy*yyy + zzz*zzz)
             RRR1 = 0.5*xxx+sqrt(0.83*0.83*yyy*yyy + zzz*zzz)
             RRR2 = 0.5*xxx+sqrt(0.85*0.85*yyy*yyy + zzz*zzz)
             RRR3 = 0.6*xxx+sqrt(0.83*0.83*yyy*yyy + zzz*zzz)

             if (dx_BLK(iBLK)>0.125 .and. RRR2>12. .and. &
                  RRR2<18.  .and. xxx > -3. .and. xxx < 3. .and. &
                  yyy < -6.  .and. zzz > 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR2>12. .and. &
                  RRR2<18.  .and. xxx > -3. .and. xxx < 3. .and. &
                  yyy > 6. .and. zzz < 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR1>12. .and. &
                  RRR1<18.  .and. xxx > 3. .and. xxx < 6. .and. &
                  yyy <  -4.5  .and. zzz > 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR1>12. .and. &
                  RRR1<18.  .and. xxx > 3. .and. xxx < 6. .and. &
                  yyy > 4.5 .and. zzz < 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR3>12. .and. &
                  RRR3<18.  .and. xxx > 6. .and. xxx < 9.  .and. &
                  yyy < 0.  .and. zzz > 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR3>12. .and. &
                  RRR3<18.  .and. xxx > 6. .and. xxx < 9. .and. &
                  yyy >  0. .and. zzz < 0.) DoRefine = .true.
             !

             if (dx_BLK(iBLK)>0.125 .and. RRR>12. .and. &
                  RRR<16.5  .and. xxx > 9.  .and. &
                  abs(zzz)<9. .and. abs(yyy)< 12. &
                  ) DoRefine = .true.


             if (dx_BLK(iBLK)>0.125 .and. RRR>12. .and. &
                  RRR<16.  .and. xxx > 9.  .and. &
                  abs(zzz)<6. .and. abs(yyy)< 13.5 &
                  ) DoRefine = .true.
             !
             if (dx_BLK(iBLK)>0.0625 .and. RRR2>14. .and. &
                  RRR2<18.  .and. xxx > -1.5 .and. xxx < 3. .and. &
                  yyy < -9.  .and. zzz > 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. RRR2>14. .and. &
                  RRR2<18.  .and. xxx > -1.5 .and. xxx < 3. .and. &
                  yyy > 9. .and. zzz < 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. RRR1>14. .and. &
                  RRR1<17.  .and. xxx > 3. .and. xxx < 6. .and. &
                  yyy < -6.  .and. zzz > 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. RRR1>14. .and. &
                  RRR1<17.  .and. xxx > 3. .and. xxx < 6. .and. &
                  yyy > 6. .and. zzz < 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. RRR3>14. .and. &
                  RRR3<17.  .and. xxx > 6. .and. xxx < 9.  .and. &
                  yyy < 0.  .and. zzz > 0.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. RRR3>14. .and. &
                  RRR3<17.  .and. xxx > 6. .and. xxx < 9. .and. &
                  yyy >  0. .and. zzz < 0.) DoRefine = .true.
             !

             if (dx_BLK(iBLK)>0.0625 .and. RRR>12.5 .and. &
                  RRR<15.  .and. xxx > 9.  .and. &
                  abs(zzz)<7.5 .and. abs(yyy)< 9. &
                  ) DoRefine = .true.
             !
             !
             !                  if (dx_BLK(iBLK)>0.125 .and.  abs(xxx) < 8. .and. &
!!!                    abs(yyy) < 8. .and. abs(zzz)< 8. ) DoRefine = .true.
!!!
!!!                   if (dx_BLK(iBLK)>0.125 .and.  xxx > 0. .and. xxx < 13.5 .and. &
!!!                    abs(yyy) < 12. .and. abs(zzz)< 12. ) DoRefine = .true.
!!!
!!!                    if (dx_BLK(iBLK)>0.125 .and.  xxx > 0. .and. xxx < 15. .and. &
!!!                    abs(yyy) < 6. .and. abs(zzz)< 6. ) DoRefine = .true.

             !!                    if (dx_BLK(iBLK)>0.125 .and. xxx < 0. .and. &
             !!                     xxx > - 24.0 .and. &
             !!                    abs(yyy)< 15. .and. abs(zzz)<3.) DoRefine = .true.

             !!                  if (dx_BLK(iBLK)>0.0625 .and. xxx < -6. .and. &
             !!                     xxx > - 18.0 .and. &
             !!                    abs(yyy)<12. .and. abs(zzz)<1.5) DoRefine = .true.
             !!
             !!                  if (dx_BLK(iBLK)>0.03125 .and. xxx < -10. .and. &
             !!                     xxx > - 16.0 .and. &
             !!                    abs(yyy)<9. .and. abs(zzz)<0.75) DoRefine = .true.
          end if
       end if

    case ('ror6_mp_135')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<48.)  DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -15. .and. xxx < 33.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>1. .and. xxx > -9. .and. xxx < 33.0 .and. &
                  abs(yyy)<36..and. abs(zzz)<36.) DoRefine = .true.


             if (dx_BLK(iBLK)>1. .and. xxx < 0. .and. xxx > -111. .and. &
                  abs(yyy)<36..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 15.0 .and. xxx > -9. .and. &
                  abs(yyy)<18..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx > -9. .and. &
                  abs(yyy)<24..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 18. .and. xxx> -87.0 .and. &
                  abs(yyy)<24..and. abs(zzz)<18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0. .and. xxx> -48. .and. &
                  abs(yyy)<18. .and. abs(zzz)<6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> -3.0 .and. &
                  abs(yyy)<18. .and. abs(zzz)<15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 18. .and. xxx> -6.0 .and. &
                  abs(yyy)<9. .and. abs(zzz)<9.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. abs(xxx) < 12.0 .and. &
                  abs(yyy)<12..and. abs(zzz)<12.) DoRefine = .true.
             !
             ! IMF angle = 135
             if (dx_BLK(iBLK)>0.25 .and. xxx < 0 .and. xxx > -6. .and. &
                  yyy > 0. .and. yyy < 21. .and. &
                  zzz > 0. .and. zzz< 9.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 0 .and. xxx > -6. .and. &
                  yyy < 0. .and. yyy > -21. .and. &
                  zzz < 0. .and. zzz>  -9.) DoRefine = .true.

             RRR = sqrt(xxx*xxx*1.03*1.03 + 0.83*0.83*yyy*yyy + zzz*zzz)
             RRR1 = 0.5*xxx+sqrt(0.83*0.83*yyy*yyy + zzz*zzz)
             RRR2 = 0.5*xxx+sqrt(0.85*0.85*yyy*yyy + zzz*zzz)
             RRR3 = 0.6*xxx+sqrt(0.83*0.83*yyy*yyy + zzz*zzz)

             if (dx_BLK(iBLK)>0.125 .and. RRR2>12. .and. &
                  RRR2<18.  .and. xxx > -3. .and. xxx < 3. .and. &
                  yyy > 6.  .and. zzz > -6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR2>12. .and. &
                  RRR2<18.  .and. xxx > -3. .and. xxx < 3. .and. &
                  yyy<-6. .and. zzz < 6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR1>12. .and. &
                  RRR1<18.  .and. xxx > 3. .and. xxx < 6. .and. &
                  yyy > 4.5  .and. zzz > -6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR1>12. .and. &
                  RRR1<18.  .and. xxx > 3. .and. xxx < 6. .and. &
                  yyy<-4.5 .and. zzz < 6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR3>12. .and. &
                  RRR3<18.  .and. xxx > 6. .and. xxx < 9.  .and. &
                  yyy > 0.  .and. zzz > -6.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. RRR3>12. .and. &
                  RRR3<18.  .and. xxx > 6. .and. xxx < 9. .and. &
                  yyy< 0. .and. zzz < 6.) DoRefine = .true.
             !

             if (dx_BLK(iBLK)>0.125 .and. RRR>12. .and. &
                  RRR<16.5  .and. xxx > 9.  .and. &
                  abs(zzz)<9. .and. abs(yyy)< 12. &
                  ) DoRefine = .true.


             if (dx_BLK(iBLK)>0.125 .and. RRR>12. .and. &
                  RRR<16.  .and. xxx > 9.  .and. &
                  abs(zzz)<6. .and. abs(yyy)< 13.5 &
                  ) DoRefine = .true.
             !
             if (dx_BLK(iBLK)>0.0625 .and. RRR2>14. .and. &
                  RRR2<18.  .and. xxx > -1.5 .and. xxx < 3. .and. &
                  yyy > 9.  .and. zzz > -3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. RRR2>14. .and. &
                  RRR2<18.  .and. xxx > -1.5 .and. xxx < 3. .and. &
                  yyy<-9. .and. zzz < 3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. RRR1>14. .and. &
                  RRR1<17.  .and. xxx > 3. .and. xxx < 6. .and. &
                  yyy > 6.  .and. zzz > -3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. RRR1>14. .and. &
                  RRR1<17.  .and. xxx > 3. .and. xxx < 6. .and. &
                  yyy<-6. .and. zzz < 3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. RRR3>14. .and. &
                  RRR3<17.  .and. xxx > 6. .and. xxx < 9.  .and. &
                  yyy > 0.  .and. zzz > -3.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. RRR3>14. .and. &
                  RRR3<17.  .and. xxx > 6. .and. xxx < 9. .and. &
                  yyy< 0. .and. zzz < 3.) DoRefine = .true.
!!!!

             if (dx_BLK(iBLK)>0.0625 .and. RRR>12.5 .and. &
                  RRR<15.  .and. xxx > 9.  .and. &
                  abs(zzz)<9. .and. abs(yyy)< 9. &
                  ) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. RRR>12.5 .and. &
                  RRR<15.  .and. xxx > 9.  .and. &
                  abs(zzz)<3. .and. abs(yyy)< 15. &
                  ) DoRefine = .true.


             !
             !
             !                  if (dx_BLK(iBLK)>0.125 .and.  abs(xxx) < 8. .and. &
!!!                    abs(yyy) < 8. .and. abs(zzz)< 8. ) DoRefine = .true.
!!!
!!!                   if (dx_BLK(iBLK)>0.125 .and.  xxx > 0. .and. xxx < 13.5 .and. &
!!!                    abs(yyy) < 12. .and. abs(zzz)< 12. ) DoRefine = .true.
!!!
!!!                    if (dx_BLK(iBLK)>0.125 .and.  xxx > 0. .and. xxx < 15. .and. &
!!!                    abs(yyy) < 6. .and. abs(zzz)< 6. ) DoRefine = .true.

             !!                    if (dx_BLK(iBLK)>0.125 .and. xxx < 0. .and. &
             !!                     xxx > - 24.0 .and. &
             !!                    abs(yyy)< 15. .and. abs(zzz)<3.) DoRefine = .true.

             !!                  if (dx_BLK(iBLK)>0.0625 .and. xxx < -6. .and. &
             !!                     xxx > - 18.0 .and. &
             !!                    abs(yyy)<12. .and. abs(zzz)<1.5) DoRefine = .true.
             !!
             !!                  if (dx_BLK(iBLK)>0.03125 .and. xxx < -10. .and. &
             !!                     xxx > - 16.0 .and. &
             !!                    abs(yyy)<9. .and. abs(zzz)<0.75) DoRefine = .true.
          end if
       end if
       !
    case ('ror6_cluster1')
       ! 6 x 6 x 6 block for High resolution Runs On Request
       ! Box sizes: -255<X<33, -48<Z<48, -48<Y<48
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             if (dx_BLK(iBLK) > 4. .and. xxx > -159.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -39.0)  DoRefine = .true.

             if (dx_BLK(iBLK) > 2. .and. xxx > -135.0 .and. &
                  abs(yyy)<48..and. abs(zzz)<48.)  DoRefine = .true.


             if (dx_BLK(iBLK)>1. .and. xxx < 33. .and. xxx > -63. .and. &
                  abs(yyy)<36..and. abs(zzz)<36.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.5 .and. abs(xxx) < 15.  .and. &
                  abs(yyy)<30..and. abs(zzz)<24.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx> -15.0 .and. &
                  abs(yyy)<24..and. abs(zzz) < 18. ) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < 21. .and. xxx> -27.0 .and. &
                  abs(yyy)<18..and. abs(zzz )< 12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.5 .and. xxx < -15. .and. xxx> -39.0 .and. &
                  yyy<30. .and. yyy > 6. .and. abs(zzz )< 18.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.25 .and. xxx < 9. .and. xxx> -12. .and. &
                  abs(yyy)<24. .and. zzz<18. .and. zzz> -21.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 12. .and. xxx> 9. .and. &
                  abs(yyy)<24. .and. zzz<15. .and. zzz> -18.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 15. .and. xxx> 12. .and. &
                  abs(yyy)<24. .and. zzz<12. .and. zzz> -15.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < 18. .and. xxx> 15. .and. &
                  abs(yyy)<24. .and. zzz<9. .and. zzz> -12.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.25 .and. xxx < -12. .and. xxx> -36.0 .and. &
                  yyy<27. .and. yyy > 12. .and. abs(zzz )< 15.) DoRefine = .true.


             RRR1 = xxx+0.035*(yyy*yyy + zzz*zzz)
             RRR2 = xxx+0.070*(yyy*yyy + zzz*zzz)
             RRR3 = xxx+0.040*(yyy*yyy + zzz*zzz)
             RRR4 = xxx+0.060*(yyy*yyy + zzz*zzz)
             RRR5 = xxx+0.080*(yyy*yyy + zzz*zzz)
             !
             ! Don Fairfield 2a
             !!
             !!                   if (dx_BLK(iBLK)>0.125 .and. xxx > -10.5 .and. yyy> -9. .and. &
             !!                    yyy<22.5 .and. zzz<9. .and. zzz> -9. .and. &
             !!                   RRR1 < 15. .and. RRR2 > 8.) DoRefine = .true.
             !!
             !!                  if (dx_BLK(iBLK)>0.125 .and. xxx < -6. .and. xxx > -30. .and. &
             !!                    yyy> 13.5 .and. &
             !!                    yyy<22.5 .and. zzz<9. .and. zzz> -9.) DoRefine = .true. 
             !!                  
             !!                   if (dx_BLK(iBLK)>0.125 .and. xxx > -10.5 .and. yyy> 10. .and. &
             !!                    yyy<22.5 .and. zzz<9. .and. zzz> -9. .and. &
             !!                   RRR1 < 15. .and. RRR5 > 8.) DoRefine = .true.
             !!
             !
             !!                    if (dx_BLK(iBLK)>0.0625 .and. xxx > -7. .and. yyy> -4. .and. &
             !!                    yyy<19.5 .and. zzz<7.5 .and. zzz> -6. .and. &
             !!                   RRR3 < 13. .and. RRR2 > 10.) DoRefine = .true.
             !!
             !!                    if (dx_BLK(iBLK)>0.0625 .and. xxx < -6.5 .and. xxx > -24. .and. &
             !!                     yyy> 16.5 .and. &
             !!                    yyy<20.25 .and. zzz<6. .and. zzz> -4.) DoRefine = .true.
             !!
             !
             !
             ! Don Fairfield 3a 

             if (dx_BLK(iBLK)>0.125 .and. xxx > -10.5 .and. yyy> -6. .and. &
                  yyy<22.5 .and. zzz<9. .and. zzz> -9. .and. &
                  RRR1 < 15. .and. RRR2 > 9.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. xxx < -6. .and. xxx > -12. .and. &
                  yyy> 13.5 .and. &
                  yyy<22.5 .and. zzz<9. .and. zzz> -9.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.125 .and. xxx < -12. .and. xxx > -30. .and. &
                  yyy> 15. .and. &
                  yyy<24. .and. zzz<9. .and. zzz> -9.) DoRefine = .true.


             if (dx_BLK(iBLK)>0.125 .and. xxx > -10.5 .and. yyy> 10. .and. &
                  yyy<22.5 .and. zzz<9. .and. zzz> -9. .and. &
                  RRR1 < 15. .and. RRR5 > 9.) DoRefine = .true.

             !
             if (dx_BLK(iBLK)>0.0625 .and. xxx > -7. .and. yyy> -3. .and. &
                  yyy<19.5 .and. zzz<5.25 .and. zzz> -5.25 .and. &
                  RRR3 < 13. .and. RRR2 > 11.) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. xxx < -6.5 .and. xxx > -12. .and. &
                  yyy> 16.5 .and. &
                  yyy<20.25 .and. zzz<5.25 .and. zzz> -5.25) DoRefine = .true.

             if (dx_BLK(iBLK)>0.0625 .and. xxx < -12. .and. xxx > -24. .and. &
                  yyy> 18. .and. &
                  yyy<21.75 .and. zzz<5.25 .and. zzz> -5.25) DoRefine = .true.


          end if
       end if


    case ('mag_new')
       ! Refine for generic magnetosphere
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             ! Refine all blocks intersecting body
             if (minRblk <= 1.5*Rcurrents) DoRefine = .true.

             ! Refine all blocks intersecting body
             if (SizeMax > 1.0 .and. minRblk <= 5.0*Rcurrents) DoRefine = .true.

             ! Refine magnetopause/shock
             if (SizeMax > 0.5) then
                if (xxx > 0.) then

                   if (nJ == 4) then

                      x_off = 14.0 ; curve_rad = 20.0 ; mach = 20.0 ; frac = .3
                      x_off1 = x_off*(1.0-frac)
                      x_off2 = x_off*(1.0+frac)
                      curve_rad1 = curve_rad*(1.0-0.5*frac)
                      curve_rad2 = curve_rad*(1.0+0.5*frac)
                      tmpminRblk = sqrt((min(abs(yy1),abs(yy2)))**2 + &
                           (min(abs(zz1),abs(zz2)))**2)
                      tmpmaxRblk = sqrt((max(abs(yy1),abs(yy2)))**2 + &
                           (max(abs(zz1),abs(zz2)))**2)
                      bound1 = 0.0
                      if (xx2 < x_off1) then
                         bound1 =  sqrt((xx2-x_off1)**2/(mach**2-1) &
                              - 2*curve_rad1*(xx2-x_off1));
                      end if
                      bound2 = 0.0
                      if (xx1 < x_off2) then
                         bound2 =  sqrt((xx1-x_off2)**2/(mach**2-1) &
                              - 2*curve_rad2*(xx1-x_off2));
                      end if
                      if (xx1 < x_off2 .and. &
                           (tmpminRblk < bound2 .and. tmpmaxRblk > bound1) ) then
                         if (SizeMax >= 2.0) DoRefine = .true.
                         if (SizeMax >= 1.0 .and. xx1 >=  10.0) DoRefine = .true.
                         if (SizeMax >= 0.5 .and. xx1 >=  2.0) DoRefine = .true.
                      end if

                      xxx = 0.0
                      yyy = 0.50*(y_BLK(nI,nJ,nK,iBLK)+y_BLK(1,1,1,iBLK))
                      zzz = 0.50*(z_BLK(nI,nJ,nK,iBLK)+z_BLK(1,1,1,iBLK))
                      RR = sqrt( xxx*xxx + yyy*yyy + zzz*zzz )

                      if (RR < 13.0 .and. SizeMax > 1.0) then
                         DoRefine = .true.
                      endif

                      xxx = 0.50*(x_BLK(nI,nJ,nK,iBLK)+x_BLK(1,1,1,iBLK))
                      RR = sqrt( xxx*xxx + yyy*yyy + zzz*zzz )

                   else
                      if (((abs(zz1) == 0) .or. &
                           (abs(zz2) == 0)) .and. &
                           ((abs(yy1) == 0) .or. &
                           (abs(yy2) == 0))) then
                         DoRefine = .true.
                      endif
                   endif

                end if
             end if

             ! Refine tail sheet

             ! This will make it:
             !    2   Re at lt 120
             !    1   Re at lt  60
             !    0.5 Re at lt  30

             ! Thickness
             !   5 Re in Z
             !  12 Re in Y

             if (SizeMax > 0.5 .and. &
                  (TypeGeometry=='cartesian'.or.xxx>x1+(x2-x1)*cRefinedTailCutoff).and.&
                  (xxx < 0. .and. xxx > -SizeMax*28.0)) then

                if (nJ == 4) then
                   tmpminRblk = sqrt((min(abs(yy1),abs(yy2)))**2 + &
                        25.*((min(abs(zz1),abs(zz2)))**2))
                   if (tmpminRblk < 12.) DoRefine = .true.
                else
                   if (((abs(zz1) == 0) .or. &
                        (abs(zz2) == 0)) .and. &
                        ((abs(yy1) == 0) .or. &
                        (abs(yy2) == 0))) then
                      DoRefine = .true.
                   endif
                endif

             end if

             if (SizeMax >= 0.5 .and. &
                  (xxx < 10.0 .and. xxx > -32.0)) then
                if (((abs(zz1) == 0) .or. &
                     (abs(zz2) == 0)) .and. &
                     ((yyy >= -8.0) .and. &
                     (yyy <= 8.0))) then
                   DoRefine = .true.
                endif
             endif

             ! Refine tail
             !                 if (SizeMax > 4 .and. &
             !                      (xxx < 0. .and. xxx > -150.)) DoRefine = .true.


          end if

       endif

    case ('coupledhelio')
       !refine to have resolution not worse 4.0 and
       !refine the body intersecting blocks
       DoRefine=minRblk<=rBody.or.dx_BLK(iBLK)>4.01

    case ('magnetosphere')
       ! Refine for generic magnetosphere
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then                                        
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else

             ! Refine inner blocks for corotation
             if ((UseRotatingBc).and.(minRblk <= 6.0)) then
                select case(TypeGeometry)
                case('cartesian')
                   if (min(abs(zz1),abs(zz2)) < Rcurrents) DoRefine = .true.   
                case('spherical','spherical_lnr')
                   if(minRblk*min(abs(cos(XyzStart_BLK(3,iBLK)-0.5*&
                        dz_BLK(iBLK))),abs(cos(XyzStart_BLK(3,iBLK)+(nK-0.5)&
                        *dz_BLK(iBLK)))) < Rcurrents)DoRefine = .true.         
                case default
                   call stop_mpi('Unknown TypeGeometry = '//TypeGeometry)
                end select
             endif

             ! Refine magnetopause/shock
             if (SizeMax > 0.5) then
                if ((xxx > 0.).and.(minRblk < 12.0)) then

                   if (nJ == 4) then
                      DoRefine = .true.
                   else
                      if  (((abs(zz1) == 0) .or. (abs(zz2) == 0)) .and. &
                           ((abs(yy1) == 0) .or. (abs(yy2) == 0))) then
                         DoRefine = .true.
                      endif
                   endif

                end if
             end if

             ! Refine tail sheet

             ! This will make it:
             !    2   Re at lt 120
             !    1   Re at lt  60
             !    0.5 Re at lt  30

             ! Thickness
             !   5 Re in Z
             !  12 Re in Y

             if (SizeMax>0.5 .and.&
                  (TypeGeometry=='cartesian' &
                  .or. xxx>x1+(x2-x1)*cRefinedTailCutoff).and. &
                  (xxx<0. .and. xxx>-SizeMax*28.0)) then
                if (nJ == 4) then
                   select case(TypeGeometry)
                   case('cartesian')
                      tmpminRblk = sqrt((min(abs(yy1),abs(yy2)))**2 + &
                           25.*((min(abs(zz1),abs(zz2)))**2))
                   case('spherical','spherical_lnr')
                      tmpminRblk=min(cos(XyzStart_BLK(3,iBLK)-0.5*&          
                           dz_BLK(iBLK))**2,cos(XyzStart_BLK(3,iBLK)+&  
                           (nK-0.5)*dz_BLK(iBLK))**2)
                      tmpminRblk=min(sin(XyzStart_BLK(2,iBLK)-0.5*&          
                           dy_BLK(iBLK))**2,sin(XyzStart_BLK(2,iBLK)+&          
                           (nJ-0.5)*dy_BLK(iBLK))**2)*(1.0- tmpminRblk)&     
                           + 25.*tmpminRblk                                         
                      tmpminRblk=minRblk*sqrt(tmpminRblk)
                   case default
                      call stop_mpi('Unknown TypeGeometry = '//TypeGeometry)               
                   end select
                   if (tmpminRblk < 12.) DoRefine = .true.
                else
                   if  (((abs(zz1) == 0) .or. (abs(zz2) == 0)) .and. &
                        ((abs(yy1) == 0) .or. (abs(yy2) == 0))) then
                      DoRefine = .true.
                   endif
                endif

             end if

             ! Refine tail
             !                 if (SizeMax > 4 .and. &
             !                      (xxx < 0. .and. xxx > -150.)) DoRefine = .true.

             ! Refine all blocks intersecting body
             if (minRblk <= Rcurrents) then
                if(TypeGeometry=='cartesian')&
                     DoRefine = .true.
                if(TypeGeometry=='spherical')&
                     DoRefine = dy_BLK(iBLK)>cPi/128+cTiny.or.&
                     dz_BLK(iBLK)>cPi/128+cTiny
             end if
          end if
       end if

    case ('magneto12')
       ! Refine for generic magnetosphere
       if (maxRblk > Rbody) then

          if (dx_BLK(iBLK) > 8.) then
             ! Refine all blocks with dx greater than 8
             DoRefine = .true.
          else
             select case(TypeGeometry)
             case('cartesian')
                minx = minval(abs(x_BLK(1:nI, 1, 1,iBLK)))
                miny = minval(abs(y_BLK(1, 1:nJ, 1,iBLK)))
                minz = minval(abs(z_BLK(1, 1, 1:nK,iBLK)))

                minRblk = sqrt(minx*minx + miny*miny + minz*minz)
                minx = minval(x_BLK(1:nI, 1, 1,iBLK))

             case('spherical','spherical_lnr')
                minx = minval(x_BLK(1:nI, 1:nJ, 1:nK,iBLK))
                miny = minval(abs(y_BLK(1:nI, 1:nJ, 1:nK,iBLK)))
                minz = minval(abs(z_BLK(1:nI, 1:nJ, 1:nK,iBLK)))  
             case default
                call stop_mpi('Specify Refinement: Unknown TypeGeometry'//TypeGeometry)
             end select

             ! With 8x8x8 blocks
             ! levels = 7 (i.e. 1/4 Re) 495616 cells rbody+corotation
             ! levels = 8 (i.e. 1/8 Re) 1671168 cells rbody+corotation
             ! Refine inner blocks for corotation
             if ((UseRotatingBc).and.(minRblk <= 6.0).and.SizeMax>=0.25) then

                if (minz < Rcurrents+ SizeMax) DoRefine = .true.

             endif
             ! With 8x8x8 blocks
             ! levels = 7 (i.e. 1/4 Re) 811008 cells rbody+corotation+mp
             ! levels = 8 (i.e. 1/8 Re) 1671168 cells rbody+corotation+mp

             ! Refine magnetopause/shock
             ! Refine magnetopause/shock
             if (SizeMax > 0.5) then
                if ((minx > 0.).and.(RR < 20.0)) then
                   DoRefine = .true.
                end if
             end if

             ! Refine tail sheet

             ! This will make it:
             !    2   Re at lt 120
             !    1   Re at lt  60
             !    0.5 Re at lt  30

             ! Thickness
             !   5 Re in Z
             !  12 Re in Y

             ! With 8x8x8 blocks
             ! levels = 7 (i.e. 1/4 Re) 1398784 cells for rbody+corotation+mp+tail
             ! levels = 8 (i.e. 1/8 Re) 2201600 cells for rbody+corotation+mp+tail

             ! With 4x4x4 blocks
             ! levels = 8 (i.e. 1/4 Re) 742912  cells for rbody+corotation+mp+tail
             ! levels = 9 (i.e. 1/8 Re) 1214208 cells for rbody+corotation+mp+tail
             if ( SizeMax >0.4.and.&                               
                  (TypeGeometry=='cartesian'.or.&
                  minx>x1+(x2-x1)*cRefinedTailCutoff).and.&
                  (minx<0. .and. &
                  minx > -SizeMax*50.0)) then
                if (miny < 15.0*dx_BLK(iBLK) .and. minz < SizeMax) DoRefine = .true.
             end if

             ! With 8x8x8 blocks
             ! levels = 7 (i.e. 1/4 Re) 208896 cells  (rbody only)
             ! levels = 8 (i.e. 1/8 Re) 1441792 cells (rbody only)

             ! Refine all blocks intersecting body
             ! This line didn't work in previous versions, since the
             ! cell centers were further away from 0.0 than RCurrents.
             ! Now it checks to see if it is within 1 cell of rcurrents.
             if (minRblk - Rcurrents <= dx_BLK(iBLK)) then
                if(TypeGeometry=='cartesian')&
                     DoRefine = .true.
                if(TypeGeometry=='spherical'.or.&
                     TypeGeometry=='spherical_lnr')&
                     DoRefine = dy_BLK(iBLK)>cPi/128+cTiny.or.&
                     dz_BLK(iBLK)>cPi/128+cTiny
             end if

             minRblk = sqrt((min(abs(xx1),abs(xx2)))**2 + &
                  (min(abs(yy1),abs(yy2)))**2 + &
                  (min(abs(zz1),abs(zz2)))**2)

          end if
       end if
    case ('magneto_fine')

       ! Refine for generic magnetosphere
       if (maxRblk > Rbody) then
          if (dx_BLK(iBLK) > 4.) then
             ! Refine all blocks with dx greater than 4
             DoRefine = .true.
          else
             ! Refine all blocks intersecting body
             if (minRblk <= Rcurrents) DoRefine = .true.

             ! Refine magnetopause/shock
             if (SizeMax > 0.25) then
                if (xxx > 0.) then

                   if ((RR < 20.0).and.(RR > 6.0)) DoRefine = .true.
                end if
             end if

             ! Refine tail sheet

             ! This will make it:
             !    2   Re at lt 207
             !    1   Re at lt 107
             !    0.5 Re at lt  57
             !    0.25   at lt  32
             !    0.125  at lt  19.5

             if (SizeMax>0.125 .and.&
                  (TypeGeometry=='cartesian'.or.xxx>x1+(x2-x1)*cRefinedTailCutoff).and.&
                  (xxx<-7. .and. xxx>-SizeMax*50.0-5.0)) then
                select case(TypeGeometry)
                case('cartesian')
                   tmpminRblk = sqrt((min(abs(yy1),abs(yy2)))**2 + &
                        25.*((min(abs(zz1),abs(zz2)))**2))
                case('spherical')
                   tmpminRblk=min(cos(XyzStart_BLK(3,iBLK)-0.5*&          
                        dz_BLK(iBLK))**2,cos(XyzStart_BLK(3,iBLK)+&          
                        (nK-0.5)*dz_BLK(iBLK))**2)                         
                   tmpminRblk=min(sin(XyzStart_BLK(2,iBLK)-0.5*&          
                        dy_BLK(iBLK))**2,sin(XyzStart_BLK(2,iBLK)+&          
                        (nJ-0.5)*dy_BLK(iBLK))**2)*(1.0- tmpminRblk)&     
                        + 25.*tmpminRblk
                   tmpminRblk=minRblk*sqrt(tmpminRblk)
                case default
                   call stop_mpi('Unknown TypeGeometry = '//TypeGeometry)
                end select

                if (tmpminRblk < 12.) DoRefine = .true.
             end if

             ! Refine tail
             if (SizeMax>2 .and. (xxx<0. .and. xxx>-150.)) DoRefine = .true.

          end if
       end if

    case ('carrington')
       ! Refine for Carrington Event grid

       select case(TypeGeometry)
       case('spherical')
          call stop_mpi('No spherical definitions for carrington refinement')
       end select

       ! Cycle if the whole block is inside body
       if(maxRblk < Rbody) then
          RETURN
       end if

       ! Refine all blocks with dx greater than 8
       if(dx_BLK(iBLK) > 8.) then
          DoRefine = .true.
          RETURN
       end if

       ! Refine all blocks inside Rcurrents
       if(minRblk <= Rcurrents) then
          DoRefine = .true.
          RETURN
       end if

       ! Refine magnetopause/shock
       if(SizeMax > 0.25 .and. xxx > 0.) then
          tmpminRblk = sqrt( &
               1.00*((min(abs(xx1),abs(xx2)))**2) + &
               0.75*((min(abs(yy1),abs(yy2)))**2) + &
               0.75*((min(abs(zz1),abs(zz2)))**2) )
          if(tmpminRblk < 12.) then
             DoRefine = .true.
             RETURN
          end if
       end if

       ! Refine tail
       if (SizeMax>2 .and. (xxx<0. .and. xxx>-100.)) then
          tmpminRblk = sqrt((min(abs(yy1),abs(yy2)))**2 + &
               1.*((min(abs(zz1),abs(zz2)))**2))
          if (tmpminRblk < 50.) then
             DoRefine = .true.
             RETURN
          end if
       end if

       ! Refine tail sheet
       if (SizeMax>0.25 .and.&
            (xxx<0. .and. xxx>-SizeMax*32.0-3.1)) then
          tmpminRblk = sqrt((min(abs(yy1),abs(yy2)))**2 + &
               25.*((min(abs(zz1),abs(zz2)))**2))
          if (tmpminRblk < 12.) then
             DoRefine = .true.
             RETURN
          end if
       end if

    case default
       call stop_mpi(NameSub//' ERROR: unknown grid name='//trim(NameGrid))
    end select

  contains
    !======================================================================
    real function minmod(x,y)
      real, intent(in) :: x,y
      minmod = max(0.0,min(abs(x),sign(1.0,x)*y))
    end function minmod

  end subroutine user_specify_initial_refinement
  !======================================================================

end module ModUser

