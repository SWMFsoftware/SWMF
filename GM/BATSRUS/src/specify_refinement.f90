!^CFG COPYRIGHT UM
subroutine specify_initial_refinement(refb, lev)
  use ModSize
  use ModMain, ONLY : &
       UseUserSpecifyRefinement,&           !^CFG IF USERFILES
       body1,UseRotatingBc,UseMassLoading,unusedBLK
  use ModGeometry, ONLY : XyzMin_D,XyzMax_D,XyzStart_BLK,&
       dy_BLK,dz_BLK,TypeGeometry,x1,x2,&                !^CFG IF NOT CARTESIAN
       x_BLK,y_BLK,z_BLK,dx_BLK
  use ModPhysics, ONLY : Rbody,Rcurrents,Qprod
  use ModAMR, ONLY : InitialRefineType
  use ModNumConst
  implicit none

  logical, intent(out) :: refb(nBLK)
  integer, intent(in) :: lev

  integer :: iBLK
  real :: xxx,yyy,zzz,RR, xxPoint,yyPoint,zzPoint
  real :: minRblk,maxRblk,xx1,xx2,yy1,yy2,zz1,zz2, tmpminRblk,tmpmaxRblk
  real :: critx, critDX
  real :: x_off,x_off1,x_off2,curve_rad,curve_rad1,curve_rad2
  real :: mach,frac,bound1,bound2
  real :: lscale,rcyl, minx, miny, minz, maxx, maxy, maxz
  real :: SizeMax

  real,parameter::cRefinedTailCutoff=(cOne-cQuarter)*(cOne-cEighth)
  logical::IsFound
  character(len=*), parameter :: NameSub='specify_initial_refinement'
  logical :: oktest, oktest_me

  !----------------------------------------------------------------------------
  call set_oktest('initial_refinement',oktest,oktest_me)

  refb = .false.

  do iBLK = 1,nBLK
     if (.not. unusedBLK(iBLK)) then
        ! Block min, max radius values
        xx1 = cHalf*(x_BLK( 0, 0, 0,iBLK)+x_BLK(   1,   1  , 1,iBLK))
        xx2 = cHalf*(x_BLK(nI,nJ,nK,iBLK)+x_BLK(nI+1,nJ+1,nK+1,iBLK))
        yy1 = cHalf*(y_BLK( 0, 0, 0,iBLK)+y_BLK(   1,   1,   1,iBLK))
        yy2 = cHalf*(y_BLK(nI,nJ,nK,iBLK)+y_BLK(nI+1,nJ+1,nK+1,iBLK))
        zz1 = cHalf*(z_BLK( 0, 0, 0,iBLK)+z_BLK(   1,   1,   1,iBLK))
        zz2 = cHalf*(z_BLK(nI,nJ,nK,iBLK)+z_BLK(nI+1,nJ+1,nK+1,iBLK))
        
        select case(TypeGeometry)                              !^CFG IF NOT CARTESIAN
        case('cartesian')                                      !^CFG IF NOT CARTESIAN
           SizeMax=dx_BLK(iBLK)
           ! Block center coordinates
           xxx = cHalf*(x_BLK(nI,nJ,nK,iBLK)+x_BLK(1,1,1,iBLK)) 
           yyy = cHalf*(y_BLK(nI,nJ,nK,iBLK)+y_BLK(1,1,1,iBLK))
           zzz = cHalf*(z_BLK(nI,nJ,nK,iBLK)+z_BLK(1,1,1,iBLK))         
           RR = sqrt( xxx*xxx + yyy*yyy + zzz*zzz )
           minRblk = sqrt(&
             minmod(xx1,xx2)**2 + minmod(yy1,yy2)**2 + minmod(zz1,zz2)**2)
!^CFG IF NOT SIMPLE BEGIN
           maxRblk = sqrt((max(abs(xx1),abs(xx2)))**2 + &
                (max(abs(yy1),abs(yy2)))**2 + &
                (max(abs(zz1),abs(zz2)))**2)
           if(body1.and.maxRblk<rBody)CYCLE
        case('spherical')                                          !^CFG IF NOT CARTESIAN BEGIN
           SizeMax=max(dx_BLK(iBLK),dy_BLK(iBLK)*maxRblk)
           minRblk = XyzStart_BLK(1,iBLK)-cHalf*dx_BLK(iBLK)            
           maxRblk = minRblk+nI*dx_BLK(iBLK)                            
           RR=cHalf*(minRblk+maxRblk)                                   
           xxx=RR*sin(XyzStart_BLK(3,iBLK)+cHalf*(nK-1)*dz_BLK(iBLK))*& 
                cos(XyzStart_BLK(2,iBLK)+cHalf*(nJ-1)*dy_BLK(iBLK))       
           yyy=RR*sin(XyzStart_BLK(3,iBLK)+cHalf*(nK-1)*dz_BLK(iBLK))*& 
                sin(XyzStart_BLK(2,iBLK)+cHalf*(nJ-1)*dy_BLK(iBLK))       
           zzz=RR*cos(XyzStart_BLK(3,iBLK)+cHalf*(nK-1)*dz_BLK(iBLK))
        case('spherical_lnr')                                          
           SizeMax=max(dx_BLK(iBLK)*maxRblk,dy_BLK(iBLK)*maxRblk)
           minRblk =exp( XyzStart_BLK(1,iBLK)-cHalf*dx_BLK(iBLK))            
           maxRblk =exp(nI*dx_BLK(iBLK))*minRblk                            
           RR=cHalf*(minRblk+maxRblk)                                   
           xxx=RR*sin(XyzStart_BLK(3,iBLK)+cHalf*(nK-1)*dz_BLK(iBLK))*& 
                cos(XyzStart_BLK(2,iBLK)+cHalf*(nJ-1)*dy_BLK(iBLK))       
           yyy=RR*sin(XyzStart_BLK(3,iBLK)+cHalf*(nK-1)*dz_BLK(iBLK))*& 
                sin(XyzStart_BLK(2,iBLK)+cHalf*(nJ-1)*dy_BLK(iBLK))       
           zzz=RR*cos(XyzStart_BLK(3,iBLK)+cHalf*(nK-1)*dz_BLK(iBLK))
        case default
           call stop_mpi('Unknown TypeGeometry = '//TypeGeometry)
        end select                                                 !^CFG END CARTESIAN
!^CFG END SIMPLE

        select case (InitialRefineType)
        case ('none')
           ! Refine no blocks
!^CFG IF NOT SIMPLE BEGIN
        case ('all')
           ! Refine all used blocks
           refb(iBLK) = .true.

        case ('2Dbodyfocus')
           ! 

        case ('3Dbodyfocus')
           ! Refine, focusing on body
           if (maxRblk > Rbody) then
              if (lev <= 2) then
                 ! Refine all blocks first two times through
                 refb(iBLK) = .true.
              else
                 ! Refine blocks intersecting body
                 if (minRblk < 1.5*Rbody) refb(iBLK) = .true.
              end if
           end if

        case('xplanefocus')
           ! concentrate refinement around x=0 plane
           if (lev <= 2) then
              ! Refine all blocks first two times through
              refb(iBLK) = .true.
           else
              ! Refine blocks intersecting x=0 plane
              if(xx1<=0.and.0<=xx2)refb(iBLK) = .true.
           end if
        case ('spherefocus')
           ! Refine, focusing spherically - like
           ! for a body, but can be used when there is
           ! no body
           if (lev <= 2) then
              ! Refine all blocks first two times through
              refb(iBLK) = .true.
           else
              ! Refine blocks intersecting R=4.5 sphere
              if (minRblk < 4.5) refb(iBLK) = .true.
           end if

        case ('helio_init')
           ! Refine all blocks time through (starting with 1 block)
           if (lev <= 4) then
              refb(iBLK) = .true.
           else
              critx=(XyzMax_D(1)-XyzMin_D(1))/(2.0**real(lev-2))
              if ( RR < 1.10*rBody + critx ) then
                 refb(iBLK) = .true.
              else
                 refb(iBLK) = .false.
              end if
           endif

        case('helio_z=4')
           ! Refine all blocks time through (starting with 1 block)
           if (lev <= 2) then
              refb(iBLK) = .true.
           else
              !         if ((xxx < 60.0).and.(xxx >-60.0)) then
              !           if ((yyy < 60.0).and.(yyy >-60.0)) then
              if ((SizeMax <=2.0).and.(SizeMax >1.0) &
                   .and.(zzz<=10.0).and.(zzz>=-10.0)) then
                 refb(iBLK) =.true.
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
              refb(iBLK) = .true.
           else
              if (lev <= 14) then
                 critx=(XyzMax_D(1)-XyzMin_D(1))/(cTwo**real(lev-1))
                 if ( RR < 1.10*rBody + critx ) then
                    refb(iBLK) = .true.
                 else
                    refb(iBLK) = .false.
                 end if
              end if
           endif

        case ('heat_init')
           ! Refine all blocks time through (starting with 1 block)
           if (lev.le.2) then
              refb(iBLK) = .true.
           else
              critx=(XyzMax_D(1)-XyzMin_D(1))/(cTwo**real(lev-1))
              if (RR.le.(cOne+cFour)/cHundred+critx)then
                 refb(iBLK) = .true.
              else
                 refb(iBLK) = .false.
              end if
           endif

        case('cme')
           !\
           ! Initial refinement criteria for self-similar cme
           !/
           if (lev <= 2) then
              ! Refine all blocks first times through
              refb(iBLK) = .true.
           else
              if (lev <= 10) then
                 critx=(XyzMax_D(1)-XyzMin_D(1))/(cTwo**real(lev-1))
                 if ( RR < 1.10*rBody + critx ) then
                    refb(iBLK) = .true.
                 else
                    refb(iBLK) = .false.
                 end if
              end if
           endif
           if (lev .eq. 11) then
              if( (abs(yyy) < 34.0) .and. (SizeMax > 2.0) )then
                 refb(iBLK) = .true.
              endif
           endif
           if (lev .eq. 12) then
              if( (abs(yyy) < 17.0) .and. (SizeMax > 3.0) )then
                 refb(iBLK) = .true.
              endif
           endif
           if (lev .eq. 13) then
              if( (abs(yyy) < 2.0) .and. (SizeMax > 0.3) )then
                 refb(iBLK) = .true.
              endif
           endif
           if (lev .eq. 14) then
              if ((-0.8 .lt. xxx) .and. (xxx .lt. 0.8) .and. &
                   (-0.8 .lt. yyy) .and. (yyy .lt. 0.8) .and. &
                   (-0.7 .lt. zzz) .and. (zzz .lt. 2.0)) then
                 refb(iBLK) = .true.          
              endif
           endif

        case ('points')
           ! Refine around given points
           if(lev<=3)then
              ! Refine all blocks first 3 times through
              refb(iBLK) = .true.
           else
              ! List as many points as desired
              xxPoint =  0.
              yyPoint =  0.
              zzPoint =  0.
              RR = sqrt( (xxx-xxPoint)**2 + (yyy-yyPoint)**2 + (zzz-zzPoint)**2 )
              if (RR < 2.*real(nI)*dx_BLK(iBLK)) refb(iBLK) = .true.
           end if

        case ('mag_new')
           ! Refine for generic magnetosphere
           if (maxRblk > Rbody) then

              if (dx_BLK(iBLK) > 8.) then
                 ! Refine all blocks with dx greater than 8
                 refb(iBLK) = .true.
              else

                 ! Refine all blocks intersecting body
                 if (minRblk <= 1.5*Rcurrents) refb(iBLK) = .true.

                 ! Refine all blocks intersecting body
                 if (SizeMax > 1.0 .and. minRblk <= 5.0*Rcurrents) refb(iBLK) = .true.

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
                             if (SizeMax >= 2.0) refb(iBLK) = .true.
                             if (SizeMax >= 1.0 .and. xx1 >=  10.0) refb(iBLK) = .true.
                             if (SizeMax >= 0.5 .and. xx1 >=  2.0) refb(iBLK) = .true.
                          end if

                          xxx = 0.0
                          yyy = 0.50*(y_BLK(nI,nJ,nK,iBLK)+y_BLK(1,1,1,iBLK))
                          zzz = 0.50*(z_BLK(nI,nJ,nK,iBLK)+z_BLK(1,1,1,iBLK))
                          RR = sqrt( xxx*xxx + yyy*yyy + zzz*zzz )

                          if (RR < 13.0 .and. SizeMax > 1.0) then
                             refb(iBLK) = .true.
                          endif

                          xxx = 0.50*(x_BLK(nI,nJ,nK,iBLK)+x_BLK(1,1,1,iBLK))
                          RR = sqrt( xxx*xxx + yyy*yyy + zzz*zzz )

                       else
                          if (((abs(zz1) == 0) .or. &
                               (abs(zz2) == 0)) .and. &
                               ((abs(yy1) == 0) .or. &
                               (abs(yy2) == 0))) then
                             refb(iBLK) = .true.
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
                      (TypeGeometry=='cartesian'.or.xxx>x1+(x2-x1)*cRefinedTailCutoff).and.& !^CFG IF NOT CARTESIAN
                      (xxx < 0. .and. xxx > -SizeMax*28.0)) then

                    if (nJ == 4) then
                       tmpminRblk = sqrt((min(abs(yy1),abs(yy2)))**2 + &
                            25.*((min(abs(zz1),abs(zz2)))**2))
                       if (tmpminRblk < 12.) refb(iBLK) = .true.
                    else
                       if (((abs(zz1) == 0) .or. &
                            (abs(zz2) == 0)) .and. &
                            ((abs(yy1) == 0) .or. &
                            (abs(yy2) == 0))) then
                          refb(iBLK) = .true.
                       endif
                    endif

                 end if

                 if (SizeMax >= 0.5 .and. &
                      (xxx < 10.0 .and. xxx > -32.0)) then
                    if (((abs(zz1) == 0) .or. &
                         (abs(zz2) == 0)) .and. &
                         ((yyy >= -8.0) .and. &
                         (yyy <= 8.0))) then
                       refb(iBLK) = .true.
                    endif
                 endif

                 ! Refine tail
                 !                 if (SizeMax > 4 .and. &
                 !                      (xxx < 0. .and. xxx > -150.)) refb(iBLK) = .true.


              end if

           endif

        case ('coupledhelio')
           !refine to have resolution not worse 4.0 and
           !refine the body intersecting blocks
           refb(iBLK)=minRblk<=rBody.or.dx_BLK(iBLK)>4.01

        case ('magnetosphere')
           ! Refine for generic magnetosphere
           if (maxRblk > Rbody) then

              if (dx_BLK(iBLK) > 8.) then                                        
                 ! Refine all blocks with dx greater than 8
                 refb(iBLK) = .true.
              else
   
                 ! Refine inner blocks for corotation
                 if ((UseRotatingBc).and.(minRblk <= 6.0)) then
                    select case(TypeGeometry)                                     !^CFG IF NOT CARTESIAN
                    case('cartesian')                                             !^CFG IF NOT CARTESIAN
                       if (min(abs(zz1),abs(zz2)) < Rcurrents) refb(iBLK) = .true.   
                    case('spherical','spherical_lnr')                                        !^CFG IF NOT CARTESIAN BEGIN  
                       if(minRblk*min(abs(cos(XyzStart_BLK(3,iBLK)-cHalf*&           
                            dz_BLK(iBLK))),abs(cos(XyzStart_BLK(3,iBLK)+(nK-cHalf)&  
                            *dz_BLK(iBLK)))) < Rcurrents)refb(iBLK) = .true.         
                    case default
                       call stop_mpi('Unknown TypeGeometry = '//TypeGeometry)
                    end select                                                       !^CFG END CARTESIAN
                 endif

                 ! Refine magnetopause/shock
                 if (SizeMax > 0.5) then
                    if ((xxx > 0.).and.(minRblk < 12.0)) then

                       if (nJ == 4) then
                          refb(iBLK) = .true.
                       else
                          if  (((abs(zz1) == 0) .or. (abs(zz2) == 0)) .and. &
                               ((abs(yy1) == 0) .or. (abs(yy2) == 0))) then
                             refb(iBLK) = .true.
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
                      (TypeGeometry=='cartesian'.or.xxx>x1+(x2-x1)*cRefinedTailCutoff).and.& !^CFG IF NOT CARTESIAN
                      (xxx<0. .and. xxx>-SizeMax*28.0)) then
                    if (nJ == 4) then
                       select case(TypeGeometry)                           !^CFG IF NOT CARTESIAN
                       case('cartesian')                                   !^CFG IF NOT CARTESIAN
                          tmpminRblk = sqrt((min(abs(yy1),abs(yy2)))**2 + &
                               25.*((min(abs(zz1),abs(zz2)))**2))
                       case('spherical','spherical_lnr')                             !^CFG IF NOT CARTESIAN BEGIN
                          tmpminRblk=min(cos(XyzStart_BLK(3,iBLK)-cHalf*&          
                               dz_BLK(iBLK))**2,cos(XyzStart_BLK(3,iBLK)+&          
                               (nK-cHalf)*dz_BLK(iBLK))**2)                         
                          tmpminRblk=min(sin(XyzStart_BLK(2,iBLK)-cHalf*&          
                               dy_BLK(iBLK))**2,sin(XyzStart_BLK(2,iBLK)+&          
                               (nJ-cHalf)*dy_BLK(iBLK))**2)*(cOne- tmpminRblk)&     
                               + 25.*tmpminRblk                                         
                          tmpminRblk=minRblk*sqrt(tmpminRblk)
                       case default
                          call stop_mpi('Unknown TypeGeometry = '//TypeGeometry)               
                       end select                                    !^CFG END CARTESIAN
                       if (tmpminRblk < 12.) refb(iBLK) = .true.
                    else
                       if  (((abs(zz1) == 0) .or. (abs(zz2) == 0)) .and. &
                            ((abs(yy1) == 0) .or. (abs(yy2) == 0))) then
                          refb(iBLK) = .true.
                       endif
                    endif

                 end if

                 ! Refine tail
                 !                 if (SizeMax > 4 .and. &
                 !                      (xxx < 0. .and. xxx > -150.)) refb(iBLK) = .true.

                 ! Refine all blocks intersecting body
                 if (minRblk <= Rcurrents) then
                    if(TypeGeometry=='cartesian')&                !^CFG IF NOT CARTESIAN
                         refb(iBLK) = .true.
                    if(TypeGeometry=='spherical')&                !^CFG IF NOT CARTESIAN BEGIN
                         refb(iBLK) = dy_BLK(iBLK)>cPi/128+cTiny.or.&
                         dz_BLK(iBLK)>cPi/128+cTiny               !^CFG END CARTESIAN
                 end if
              end if
           end if

!^CFG END SIMPLE

        case ('magneto12')
           ! Refine for generic magnetosphere
           if (maxRblk > Rbody) then

              if (dx_BLK(iBLK) > 8.) then
                 ! Refine all blocks with dx greater than 8
                 refb(iBLK) = .true.
              else
                 select case(TypeGeometry)                       !^CFG IF NOT CARTESIAN
                 case('cartesian')                               !^CFG IF NOT CARTESIAN
                    minx = minval(abs(x_BLK(1:nI, 1, 1,iBLK)))
                    miny = minval(abs(y_BLK(1, 1:nJ, 1,iBLK)))
                    minz = minval(abs(z_BLK(1, 1, 1:nK,iBLK)))

                    minRblk = sqrt(minx*minx + miny*miny + minz*minz)
                    minx = minval(x_BLK(1:nI, 1, 1,iBLK))

                 case('spherical','spherical_lnr')              !^CFG IF NOT CARTESIAN BEGIN
                    minx = minval(x_BLK(1:nI, 1:nJ, 1:nK,iBLK))
                    miny = minval(abs(y_BLK(1:nI, 1:nJ, 1:nK,iBLK)))
                    minz = minval(abs(z_BLK(1:nI, 1:nJ, 1:nK,iBLK)))                     
                 case default
                    call stop_mpi('Specify Refinement: Unknown TypeGeometry'//TypeGeometry)
                 end select                                      !^CFG END CARTESIAN

                 ! With 8x8x8 blocks
                 ! levels = 7 (i.e. 1/4 Re) 495616 cells rbody+corotation
                 ! levels = 8 (i.e. 1/8 Re) 1671168 cells rbody+corotation                 
                 ! Refine inner blocks for corotation
                 if ((UseRotatingBc).and.(minRblk <= 6.0).and.SizeMax>=cQuarter) then
                       
                    if (minz < Rcurrents+ SizeMax) refb(iBLK) = .true.
                       
                 endif
                 ! With 8x8x8 blocks
                 ! levels = 7 (i.e. 1/4 Re) 811008 cells rbody+corotation+mp
                 ! levels = 8 (i.e. 1/8 Re) 1671168 cells rbody+corotation+mp

                 ! Refine magnetopause/shock
                    ! Refine magnetopause/shock
                 if (SizeMax > 0.5) then
                    if ((minx > 0.).and.(RR < 20.0)) then
                       refb(iBLK) = .true.
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
                      (TypeGeometry=='cartesian'.or.&           !^CFG IF NOT CARTESIAN
                      minx>x1+(x2-x1)*cRefinedTailCutoff).and.& !^CFG IF NOT CARTESIAN
                      (minx<0. .and. &
                      minx > -SizeMax*50.0)) then
                    if (miny < 15.0*dx_BLK(iBLK) .and. minz < SizeMax) refb(iBLK) = .true.
                 end if

                 ! With 8x8x8 blocks
                 ! levels = 7 (i.e. 1/4 Re) 208896 cells  (rbody only)
                 ! levels = 8 (i.e. 1/8 Re) 1441792 cells (rbody only)

                 ! Refine all blocks intersecting body
                 ! This line didn't work in previous versions, since the
                 ! cell centers were further away from 0.0 than RCurrents.
                 ! Now it checks to see if it is within 1 cell of rcurrents.
                 if (minRblk - Rcurrents <= dx_BLK(iBLK)) then
                    if(TypeGeometry=='cartesian')&                !^CFG IF NOT CARTESIAN
                         refb(iBLK) = .true.
                    if(TypeGeometry=='spherical'.or.&                !^CFG IF NOT CARTESIAN BEGIN
                         TypeGeometry=='spherical_lnr')&
                         refb(iBLK) = dy_BLK(iBLK)>cPi/128+cTiny.or.&
                         dz_BLK(iBLK)>cPi/128+cTiny               !^CFG END CARTESIAN
                 end if

!^CFG IF NOT SIMPLE BEGIN
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
                 refb(iBLK) = .true.
              else
                 ! Refine all blocks intersecting body
                 if (minRblk <= Rcurrents) refb(iBLK) = .true.

                 ! Refine magnetopause/shock
                 if (SizeMax > 0.25) then
                    if (xxx > 0.) then

                       if ((RR < 20.0).and.(RR > 6.0)) refb(iBLK) = .true.
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
                      (TypeGeometry=='cartesian'.or.xxx>x1+(x2-x1)*cRefinedTailCutoff).and.& !^CFG IF NOT CARTESIAN
                      (xxx<-7. .and. xxx>-SizeMax*50.0-5.0)) then
                       select case(TypeGeometry)                           !^CFG IF NOT CARTESIAN
                       case('cartesian')                                   !^CFG IF NOT CARTESIAN
                          tmpminRblk = sqrt((min(abs(yy1),abs(yy2)))**2 + &
                               25.*((min(abs(zz1),abs(zz2)))**2))
                       case('spherical')                             !^CFG IF NOT CARTESIAN BEGIN
                          tmpminRblk=min(cos(XyzStart_BLK(3,iBLK)-cHalf*&          
                               dz_BLK(iBLK))**2,cos(XyzStart_BLK(3,iBLK)+&          
                               (nK-cHalf)*dz_BLK(iBLK))**2)                         
                          tmpminRblk=min(sin(XyzStart_BLK(2,iBLK)-cHalf*&          
                               dy_BLK(iBLK))**2,sin(XyzStart_BLK(2,iBLK)+&          
                               (nJ-cHalf)*dy_BLK(iBLK))**2)*(cOne- tmpminRblk)&     
                               + 25.*tmpminRblk                                         
                          tmpminRblk=minRblk*sqrt(tmpminRblk)
                       case default
                          call stop_mpi('Unknown TypeGeometry = '//TypeGeometry)               
                       end select                                    !^CFG END CARTESIAN

                    if (tmpminRblk < 12.) refb(iBLK) = .true.
                 end if

                 ! Refine tail
                 if (SizeMax>2 .and. (xxx<0. .and. xxx>-150.)) refb(iBLK) = .true.

              end if
           end if

        case ('carrington')
           ! Refine for Carrington Event grid

           select case(TypeGeometry)                     !^CFG IF NOT CARTESIAN BEGIN
           case('spherical')
              call stop_mpi('No spherical definitions for carrington refinement')               
           end select                                    !^CFG END CARTESIAN

           ! Cycle if the whole block is inside body
           if(maxRblk < Rbody) then
              CYCLE
           end if

           ! Refine all blocks with dx greater than 8
           if(dx_BLK(iBLK) > 8.) then
              refb(iBLK) = .true.
              CYCLE
           end if

           ! Refine all blocks inside Rcurrents
           if(minRblk <= Rcurrents) then
              refb(iBLK) = .true.
              CYCLE
           end if

           ! Refine magnetopause/shock
           if(SizeMax > 0.25 .and. xxx > 0.) then
              tmpminRblk = sqrt( &
                   1.00*((min(abs(xx1),abs(xx2)))**2) + &
                   0.75*((min(abs(yy1),abs(yy2)))**2) + &
                   0.75*((min(abs(zz1),abs(zz2)))**2) )
              if(tmpminRblk < 12.) then
                 refb(iBLK) = .true.
                 CYCLE
              end if
           end if

           ! Refine tail
           if (SizeMax>2 .and. (xxx<0. .and. xxx>-100.)) then
              tmpminRblk = sqrt((min(abs(yy1),abs(yy2)))**2 + &
                   1.*((min(abs(zz1),abs(zz2)))**2))
              if (tmpminRblk < 50.) then
                 refb(iBLK) = .true.
                 CYCLE
              end if
           end if

           ! Refine tail sheet
           if (SizeMax>0.25 .and.&
                (xxx<0. .and. xxx>-SizeMax*32.0-3.1)) then
              tmpminRblk = sqrt((min(abs(yy1),abs(yy2)))**2 + &
                   25.*((min(abs(zz1),abs(zz2)))**2))
              if (tmpminRblk < 12.) then
                 refb(iBLK) = .true.
                 CYCLE
              end if
           end if

        case ('magnetosaturn')
           ! Refine for Saturnian magnetosphere
           if (maxRblk > Rbody) then

              if (dx_BLK(iBLK) > 8.) then
                 ! Refine all blocks with dx greater than 8
                 refb(iBLK) = .true.
              else
                 ! Refine all blocks intersecting Rcurrents
                 if (minRblk <= Rcurrents .and. SizeMax > 0.5) &
                      refb(iBLK) = .true.

                 ! Refine all blocks intersecting body
                 if (minRblk <= Rbody .and. SizeMax > 0.25) refb(iBLK) = .true.

                 ! Refine bowshock
                 if (SizeMax >= 1.0 .and. (.not. UseMassLoading)) then
                    !with mass loading    x_off = 28.0 ; curve_rad = 50.0 ; mach = 8.0 ; frac = .2
                    x_off = 23.0 ; curve_rad = 35.0 ; mach = 8.0 ; frac = .1
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
                       if (SizeMax >= 8.0 .and. xx1 > -100.0) refb(iBLK) = .true.
                       if (SizeMax >= 4.0 .and. xx1 >  -75.0) refb(iBLK) = .true.
                       if (SizeMax >= 2.0 .and. xx1 >=   0.0) refb(iBLK) = .true.
                       if (SizeMax >= 1.0 .and. xx1 >=  13.0) refb(iBLK) = .true.
                    end if
                 end if

                 ! Refine magnetopause
                 if (SizeMax >= 1.0 .and. (.not. UseMassLoading)) then
                    !with mass loading    x_off = 27.0 ; curve_rad = 20.0 ; mach = 8.0 ; frac = .25
                    x_off = 18.0 ; curve_rad = 14.0 ; mach = 8.0 ; frac = .10
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
                       if (SizeMax >= 8.0 .and. xx1 > -100.0) refb(iBLK) = .true.
                       if (SizeMax >= 4.0 .and. xx1 >  -75.0) refb(iBLK) = .true.
                       if (SizeMax >= 2.0 .and. xx1 >=   0.0) refb(iBLK) = .true.
                       if (SizeMax >= 1.0 .and. xx1 >=   0.0) refb(iBLK) = .true.
                    end if
                 end if

                 ! Refine tail sheet
                 ! This will make it:
                 !    3.125    Re at lt 218
                 !    1.5625   Re at lt 109
                 !    0.78125  Re at lt  55
                 ! Thickness
                 !   2 blocks thick in Z
                 !  18.0 Re in Y
                 if (SizeMax > 1.5 .and. &
                      (TypeGeometry=='cartesian'.or.xxx>x1+(x2-x1)*cRefinedTailCutoff).and.& !^CFG IF NOT CARTESIAN
                      (xxx < 0. .and. xxx > -SizeMax*30.0) .and. &
                      (.not. UseMassLoading)) then
                    tmpminRblk = sqrt((min(abs(yy1),abs(yy2)))**2 + &
                         25.*((min(abs(zz1),abs(zz2)))**2))
                    if (tmpminRblk < 18.) refb(iBLK) = .true.
                 end if

                 if (UseMassLoading) then
                    ! Refine the Titan torus
                    if (SizeMax > 1.5 .and. zz1 < 8. .and. &
                         zz2 > -8.0 .and. minRblk < 28.0) refb(iBLK)=.true.

                    ! Refine the icy satellite torus
                    if (SizeMax > .5 .and. zz1 < .5 .and. &
                         zz2 > -.5 .and. minRblk < 15.0) refb(iBLK)=.true.
                    if (SizeMax > .25 .and. zz1 < .5 .and. &
                         zz2 > -.5 .and. minRblk < 10.0) refb(iBLK)=.true.
                 end if
              end if
           end if

        case ('magnetojupiter')
           ! Refine for Jovian magnetosphere
           if (maxRblk > Rbody) then

              if (dx_BLK(iBLK) > 32.) then
                 ! Refine all blocks with dx greater than 32
                 refb(iBLK) = .true.
              else
                 ! Refine all blocks intersecting Rcurrents
                 if (minRblk <= Rcurrents .and. SizeMax > 0.5) &
                      refb(iBLK) = .true.

                 ! Refine all blocks intersecting body
                 if (minRblk <= Rbody .and. SizeMax > 0.25) refb(iBLK) = .true.

	         ! Refine all blocks that encounter the nominal magnetosphere
                 if (xx2>-512.0 .and. SizeMax > 16.) then 

                    tmpminRblk = sqrt((min(abs(yy1),abs(yy2)))**2 + &
                         (min(abs(zz1),abs(zz2)))**2)

                    if (tmpminRblk < 256.) refB(iBLK) = .true.

		 end if

                 ! Refine bowshock
                 if (SizeMax >= 2.0 .and. (.not. UseMassLoading)) then
                    !with mass loading    x_off = 65.0 ; curve_rad = 85.0 ; mach = 8.0 ; frac = .1
                    x_off = 56.0 ; curve_rad = 80.0 ; mach = 8.0 ; frac = .1
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
                       if (SizeMax > 16.0 .and. xx1 > -300.0) refb(iBLK) = .true.
                       if (SizeMax >  8.0 .and. xx1 > -150.0) refb(iBLK) = .true.
                       if (SizeMax >  4.0 .and. xx1 >= -16.0) refb(iBLK) = .true.
                       if (SizeMax >  2.0 .and. xx1 >=  16.0) refb(iBLK) = .true.
                    end if
                 end if

                 ! Refine magnetopause
                 if (SizeMax > 2.0) then
                    x_off = 45.0 ; curve_rad = 28.0 ; mach = 8.0 ; frac = .25
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
                       if (SizeMax >  8.0 .and. xx1 >  -300.0) refb(iBLK) = .true.
                       if (SizeMax >  4.0 .and. xx1 >= -150.0) refb(iBLK) = .true.
                       if (SizeMax >  2.0 .and. xx1 >=  -25.0) refb(iBLK) = .true.
                    end if
                 end if

                 ! Refine tail sheet
                 ! This will make it:
                 !    3.125    Re at lt 218
                 !    1.5625   Re at lt 109
                 !    0.78125  Re at lt  55
                 ! Thickness
                 !   2 blocks thick in Z
                 !  18.0 Re in Y
                 if (SizeMax > 2.0 .and. &
                      (TypeGeometry=='cartesian'.or.xxx>x1+(x2-x1)*cRefinedTailCutoff).and.& !^CFG IF NOT CARTESIAN
                      (xxx < 0. .and. xxx > -SizeMax*35.0)) then

                    tmpminRblk = sqrt((min(abs(yy1),abs(yy2)))**2 + &
                         50.*((min(abs(zz1),abs(zz2)))**2))
                    if (tmpminRblk < 7.*dx_BLK(iBLK)) refb(iBLK) = .true.
                 end if

                 if (UseMassLoading) then
                    ! Refine the IO torus
                    if (SizeMax >= 2.0 .and. zz1 < 2.0 .and. &
                         zz2 > -2.0 .and. minRblk < 25.0) refb(iBLK)=.true.
                    if (SizeMax >= 1.0 .and. zz1 < 1.0 .and. &
                         zz2 > -1.0 .and. minRblk < 15.0) refb(iBLK)=.true.
                    if (SizeMax >= 0.5 .and. zz1 < 0.5 .and. &
                         zz2 > -0.5 .and. minRblk < 8.0) refb(iBLK)=.true.
                 end if
              end if
           end if


        case ('paleo')
           ! Refine for a quadrupole paleo magnetosphere
           if (maxRblk > Rbody) then
              if (dx_BLK(iBLK) > 2.) then
                 ! Refine all blocks with dx greater than 2
                 refb(iBLK) = .true.
              else
                 ! Refine all blocks intersecting body
                 if (minRblk < 1.5*Rbody) refb(iBLK) = .true.

                 ! Refine magnetopause/shock
                 if (SizeMax > 0.25) then
                    x_off = 5.5 ; curve_rad = 9.0 ; mach = 8.0 ; frac = .2
                    x_off1 = x_off*(1.0-frac)
                    x_off2 = x_off*(1.0+frac)
                    curve_rad1 = curve_rad*(1.0-frac)
                    curve_rad2 = curve_rad*(1.0+frac)
                    tmpminRblk = sqrt((min(abs(yy1),abs(yy2)))**2 + &
                         (min(abs(zz1),abs(zz2)))**2)
                    tmpmaxRblk = sqrt((max(abs(yy1),abs(yy2)))**2 + &
                         (max(abs(zz1),abs(zz2)))**2)
                    bound1 = 0.0
                    bound2 = 0.0
                    if (xx2 < x_off1 .and. xx1 < x_off2) then
                       bound1 =  sqrt((xx2-x_off1)**2/(mach**2-1) &
                            - 2*curve_rad1*(xx2-x_off1));
                       bound2 =  sqrt((xx1-x_off2)**2/(mach**2-1) &
                            - 2*curve_rad2*(xx1-x_off2));
                       if (SizeMax > 1.0 .and. &
                            (tmpminRblk < bound2 .and. tmpmaxRblk > bound1) .and. &
                            (xx1 > -16.*dx_BLK(iBLK))) refb(iBLK) = .true.
                       if (SizeMax > 0.5 .and. &
                            (tmpminRblk < bound2 .and. tmpmaxRblk > bound1) .and. &
                            (xx1 > -8.*dx_BLK(iBLK))) refb(iBLK) = .true.
                    else
                       if (minRBLK < SizeMax*curve_rad2 .and. &
                            xx1 <= x_off2) refb(iBLK) = .true.
                    end if
                 end if


		 ! Refine tail
		 if (SizeMax > 0.5 .and. &
		      (xxx < 0. .and. xxx > -25.)) then
		    if (sqrt(yyy*yyy+zzz*zzz) < 8.*Rbody ) refb(iBLK) = .true.
		 end if
                 if (SizeMax > 0.25 .and. &
                      (xxx < 0. .and. xxx > -8.)) then
                    if (sqrt(yyy*yyy+zzz*zzz) < abs(xxx)+2 ) refb(iBLK) = .true.
                 end if

              end if
           end if

        case ('comet')
           ! Refine for comet interaction
           if (maxRblk > Rbody) then
              lscale = 0.4*Qprod/7E29
              ! use Halley/linear scaling as characteristic length

              if (SizeMax > 2*lscale) then
                 refb(iBLK) = .true.
                 ! Refine all blocks with dx greater than 4*lscale
              else
                 ! Refine near bowshock if dx greater than 0.05 lscale
                 if (xxx >= 0 .and. maxRblk >= 0.75*lscale .and. &
                      minRblk <= 1.25*lscale .and. &
                      SizeMax >= 0.05*lscale) then 
                    refb(iBLK) = .true.
                 endif
                 rcyl=sqrt( min(abs(zz1),abs(zz2))**2 + &
                      min(abs(yy1),abs(yy2))**2 )
                 !Refine in downstream, tail region
                 if (xxx <=0 .and. rcyl <= 0.5*lscale .and. &
                      SizeMax > 0.1*lscale) then
                    refb(iBLK) = .true.
                 endif
                 !Refine to make grid proportional to r
                 if (SizeMax > 0.1*maxRblk .and. &
                      minRblk < 0.25*lscale ) then
                    refb(iBLK) = .true.
                 endif
                 if (SizeMax > 0.05*maxRblk .and. &
                      minRblk < 0.1*lscale ) then
                    refb(iBLK) = .true.
                 endif
!^CFG END SIMPLE
              endif
           endif
        case default
           IsFound=.false.                       !^CFG IF USERFILES BEGIN
           if (UseUserSpecifyRefinement) &
                call user_specify_initial_refinement(iBLK,refb(iBLK),lev,dx_BLK(iBLK), &
                xxx,yyy,zzz,RR,minx,miny,minz,minRblk,maxx,maxy,maxz,maxRblk,IsFound)
           if(.not.IsFound) &                    !^CFG IF USERFILES END
           call stop_mpi(NameSub//' ERROR: unknown InitialRefineType='// &
                trim(InitialRefineType)//'!!!')
        end select
     end if
  end do

  call specify_area_refinement(refb)

contains

  real function minmod(x,y)
    real, intent(in) :: x,y
    minmod = max(cZero,min(abs(x),sign(cOne,x)*y))
  end function minmod

end subroutine specify_initial_refinement

!==============================================================================

subroutine specify_area_refinement(DoRefine_B)

  !DESCRIPTION:
  ! Set DoRefine_B to .true. for blocks touching the predefined areas
  ! if the area has a finer resolution than the block

  use ModProcMH,   ONLY: iProc
  use ModMain,     ONLY: MaxBlock, nBlock, nBlockMax, nI, nJ, nK, UnusedBlk
  use ModAMR,      ONLY: nArea, Area_I
  use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, dx_BLK

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  logical, intent(inout) :: DoRefine_B(MaxBlock)

  !LOCAL VARIABLES:
  integer, parameter :: nCorner = 8
  real     :: CornerOrig_DI(3, nCorner), Corner_DI(3, nCorner)
  real     :: DistMin_D(3), DistMax_D(3)
  integer  :: i, j, k, iDim, iArea, iBlock, iCorner
  real     :: CurrentResolution

  character(len=*), parameter :: NameSub = 'specify_area_refinement'

  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------
  if(nArea <= 0) RETURN

  call set_oktest(NameSub,DoTest,DoTestMe)
  if(DoTestMe)write(*,*)NameSub,' nArea, nBlock, nBlockMax, MaxBlock=',&
       nArea, nBlock, nBlockMax, MaxBlock

  BLOCK: do iBlock = 1, nBlockMax

     if( UnusedBlk(iBlock) ) CYCLE BLOCK

     ! No need to check block if it is to be refined already
     if(DoRefine_B(iBlock)) CYCLE BLOCK

     CurrentResolution = dx_BLK(iBlock)

     iCorner = 0
     do k=1,nK+1,nK; do j=1,nJ+1,nJ; do i=1,nI+1,nI
        iCorner = iCorner+1
        CornerOrig_DI(1,iCorner) = 0.125 * sum(x_BLK(i-1:i,j-1:j,k-1:k,iBlock))
        CornerOrig_DI(2,iCorner) = 0.125 * sum(y_BLK(i-1:i,j-1:j,k-1:k,iBlock))
        CornerOrig_DI(3,iCorner) = 0.125 * sum(z_BLK(i-1:i,j-1:j,k-1:k,iBlock))
     end do; end do; end do

     AREA: do iArea = 1, nArea

        ! No need to check area if block is finer than area resolution
        if(Area_I(iArea) % Resolution >= CurrentResolution) CYCLE AREA

        ! Check if area refines the whole domain
        if(Area_I(iArea) % Name == 'all')then
           DoRefine_B(iBlock) = .true.
           CYCLE AREA
        endif
           
        ! Shift corner coordinates to the center of area
        do iCorner = 1, nCorner
           Corner_DI(:,iCorner) = &
                CornerOrig_DI(:,iCorner) - Area_I(iArea) % Center_D
        end do

        ! Rotate corners into the orientation of the area if required
        if(Area_I(iArea) % DoRotate) &
             Corner_DI = matmul(Area_I(iArea) % Rotate_DD, Corner_DI)

        ! Normalize coordinates to the size of the area in all 3 directions
        do iCorner = 1, nCorner
           Corner_DI(:,iCorner) = Corner_DI(:,iCorner) / Area_I(iArea) % Size_D
        end do

        ! Calculate maximum and minimum distances in all 3 directions
        do iDim = 1, 3
           DistMax_D(iDim) = maxval(abs(Corner_DI(iDim,:)))

           if( maxval(Corner_DI(iDim,:))*minval(Corner_DI(iDim,:)) <= 0.0)then
              ! The block covers the center point in this dimension
              DistMin_D(iDim) = 0.0
           else
              ! Select the point that is closer in this dimension
              DistMin_D(iDim) = minval(abs(Corner_DI(iDim,:)))
           end if
        end do

        ! Check if this area is intersecting with the block
        select case( Area_I(iArea) % Name)
        case('brick')
           if( all( DistMin_D < 1.0 ) ) DoRefine_B(iBlock) = .true.

        case('sphere')
           if( sum(DistMin_D**2) < 1.0 ) DoRefine_B(iBlock) = .true.

        case('shell')
           ! Check if block intersects with the enclosing sphere
           ! but it is not fully inside the inner sphere
           if(  sum(DistMin_D**2) < 1.0 .and. &
                sum(DistMax_D**2) > Area_I(iArea) % Radius1**2 ) &
                DoRefine_B(iBlock) = .true.

        case('cylinderx')
           if( DistMin_D(1) < 1.0 .and. sum(DistMin_D(2:3)**2) < 1.0 ) &
                DoRefine_B(iBlock) = .true.

        case('cylindery')
           if( DistMin_D(2) < 1.0 .and. sum(DistMin_D(1:3:2)**2) < 1.0 ) &
                DoRefine_B(iBlock) = .true.

        case('cylinderz')
           if( DistMin_D(3) < 1.0 .and. sum(DistMin_D(1:2)**2) < 1.0 ) &
                DoRefine_B(iBlock) = .true.

        case('ringx')
           ! Check if block intersects with the enclosing cylinder
           ! but it is not fully inside the inner cylinder
           if( DistMin_D(1) < 1.0 .and. sum(DistMin_D(2:3)**2) < 1.0    &
                .and. sum(DistMax_D(2:3)**2) > Area_I(iArea) % Radius1**2 ) &
                DoRefine_B(iBlock) = .true.

        case('ringy')
           ! Check if block intersects with the enclosing cylinder
           ! but it is not fully inside the inner cylinder
           if( DistMin_D(2) < 1.0 .and. sum(DistMin_D(1:3:2)**2) < 1.0 &
                .and. sum(DistMax_D(1:3:2)**2) > Area_I(iArea) % Radius1**2 ) &
                DoRefine_B(iBlock) = .true.

        case('ringz')
           ! Check if block intersects with the enclosing cylinder
           ! but it is not fully inside the inner cylinder
           if( DistMin_D(3) < 1.0 .and. sum(DistMin_D(1:2)**2) < 1.0 &
                .and. sum(DistMax_D(1:2)**2) > Area_I(iArea) % Radius1**2 ) &
                DoRefine_B(iBlock) = .true.

        case default
           call stop_mpi(NameSub // &
                ' ERROR: Unknown NameArea = ',Area_I(iArea) % Name)

        end select

        ! No need to check more areas if block is to be refined already
        if(DoRefine_B(iBlock)) EXIT AREA

     end do AREA

  end do BLOCK

  if(DoTest)write(*,*)NameSub,' on iProc=',iProc,&
       ' number of selected blocks=',count(DoRefine_B)

end subroutine specify_area_refinement
