!^CFG COPYRIGHT UM
subroutine amr_criteria(ref_criteria)
  use ModMain
  use ModGeometry, ONLY:x_BLK,y_BLK,z_BLK,R_BLK,dx_BLK,dy_BLK,dz_BLK,true_cell
  use ModAdvance
  use ModAMR,      ONLY:nRefineCrit,RefineCrit
  use ModPhysics,  ONLY:cosTHETAtilt,sinTHETAtilt,Rcurrents
  use ModPhysics,  ONLY:UseSunEarth
  use ModConst
  use ModUser, ONLY: user_amr_criteria
  implicit none

  real, intent(out) :: ref_criteria(4,nBLK)
  real :: userCriteria

  logical :: UseSwitchAMR,IsFound
  integer :: iBLK, iCrit, iVar, i, j, k
  real :: xxx,yyy,zzz,RR, dsMIN,dsMAX,ds2, RcritAMR,AMRsort_1,AMRsort_2
  real :: XTilt, YTilt, ZTilt, Zmax, ZTiltmin, ZTiltmax

  real, dimension(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn) :: outVAR,&
       Rho_G, RhoUx_G, RhoUy_G, RhoUz_G, Bx_G, By_G, Bz_G, P_G

  ! initialize all criteria to zero
  ref_criteria = cZero
  do iBLK=1,nBLK
     if (unusedBLK(iBLK)) CYCLE
     ! set 4th criteria to block radius, used in amr_physics to preserve symmetry
     ref_criteria(4,iBLK) = maxval(R_BLK(1:nI,1:nJ,1:nK,iBLK))
     ! Initialize values to use below for criteria
     dsMIN = min(dx_BLK(iBLK), dy_BLK(iBLK), dz_BLK(iBLK))
     dsMAX = max(dx_BLK(iBLK), dy_BLK(iBLK), dz_BLK(iBLK))
     ds2 = dsMIN*dsMIN
     if (UseSunEarth) then
        RcritAMR = cOne+(cOne+cHalf)*cSqrtTwo*dsMAX
     else
        RcritAMR = cZero
     end if
    
     do k=1-gcn,nK+gcn; do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
        Rho_G(i,j,k)  = State_VGB(rho_,i,j,k,iBLK)
        RhoUx_G(i,j,k)= State_VGB(rhoUx_,i,j,k,iBLK)
        RhoUy_G(i,j,k)= State_VGB(rhoUy_,i,j,k,iBLK)
        RhoUz_G(i,j,k)= State_VGB(rhoUz_,i,j,k,iBLK)
        Bx_G(i,j,k)   = State_VGB(Bx_,i,j,k,iBLK)
        By_G(i,j,k)   = State_VGB(By_,i,j,k,iBLK)
        Bz_G(i,j,k)   = State_VGB(Bz_,i,j,k,iBLK)
        P_G(i,j,k)    = State_VGB(P_,i,j,k,iBLK)
     end do; end do; end do

     do iCrit=1,nRefineCrit
        select case(RefineCrit(iCrit))
        case('gradT','gradt')
           ! Temperature gradient.
           outVAR = P_G/Rho_G
           call grad1D(1,iBLK,outVAR,gradX_VAR,gradY_VAR,gradZ_VAR,"none",0)
           ref_criteria(iCrit,iBLK) = maxval(ds2*sqrt( &
                gradX_VAR(1:nI,1:nJ,1:nK)**2 + &
                gradY_VAR(1:nI,1:nJ,1:nK)**2 + &
                gradZ_VAR(1:nI,1:nJ,1:nK)**2))
        case('gradlogrho')
           ! Log of density gradient.
           outVAR = log10(Rho_G)
           call grad1D(1, iBLK, outVAR, gradX_VAR,gradY_VAR,gradZ_VAR, "none",0)
           ref_criteria(iCrit,iBLK) = maxval(ds2*sqrt( &
                gradX_VAR(1:nI,1:nJ,1:nK)**2 + &
                gradY_VAR(1:nI,1:nJ,1:nK)**2 + &
                gradZ_VAR(1:nI,1:nJ,1:nK)**2))
        case('gradlogP','gradlogp')
           ! Log of pressure gradient
           outVAR = log10(P_G)
           call grad1D(1, iBLK, outVAR, gradX_VAR,gradY_VAR,gradZ_VAR, "none",0)
           ref_criteria(iCrit,iBLK) = maxval(ds2*sqrt( &
                gradX_VAR(1:nI,1:nJ,1:nK)**2 + &
                gradY_VAR(1:nI,1:nJ,1:nK)**2 + &
                gradZ_VAR(1:nI,1:nJ,1:nK)**2))
        case('gradP','gradp')
           ! Pressure gradient 2.
           call grad1D(1, iBLK, P_G, gradX_VAR,gradY_VAR,gradZ_VAR, "none",0)
           ref_criteria(iCrit,iBLK) = maxval(ds2*sqrt( &
                gradX_VAR(1:nI,1:nJ,1:nK)**2 + &
                gradY_VAR(1:nI,1:nJ,1:nK)**2 + &
                gradZ_VAR(1:nI,1:nJ,1:nK)**2))
        case('gradE')
           ! Electric field gradient.
           outVAR = sqrt( &
                ( -((RhoUy_G/Rho_G)* &
                   (Bz_G+B0zCell_BLK(:,:,:,iBLK)) - &
                   (RhoUz_G/Rho_G)* &
                   (By_G+B0yCell_BLK(:,:,:,iBLK))) )**2 + &
                ( -((RhoUz_G/Rho_G)* &
                   (Bx_G+B0xCell_BLK(:,:,:,iBLK)) - &
                   (RhoUx_G/Rho_G)* &
                   (Bz_G+B0zCell_BLK(:,:,:,iBLK))) )**2 + &
                ( -((RhoUx_G/Rho_G)* &
                   (By_G+B0yCell_BLK(:,:,:,iBLK)) - &
                   (RhoUy_G/Rho_G)* &
                   (Bx_G+B0xCell_BLK(:,:,:,iBLK))) )**2 )
           call grad1D(1, iBLK, outVAR, gradX_VAR,gradY_VAR,gradZ_VAR, "none",0)
           ref_criteria(iCrit,iBLK) = maxval(ds2*sqrt( &
                gradX_VAR(1:nI,1:nJ,1:nK)**2 + &
                gradY_VAR(1:nI,1:nJ,1:nK)**2 + &
                gradZ_VAR(1:nI,1:nJ,1:nK)**2))
        case('curlV','curlv','curlU','curlu')
           ! Curl of velocity
           call grad1D(1, iBLK, RhoUx_G/Rho_G, &
                gradX_Ux,gradY_Ux,gradZ_Ux, "none",Ux_)
           call grad1D(1, iBLK, RhoUy_G/Rho_G, &
                gradX_Uy,gradY_Uy,gradZ_Uy, "none",Uy_)
           call grad1D(1, iBLK, RhoUz_G/Rho_G, &
                gradX_Uz,gradY_Uz,gradZ_Uz, "none",Uz_)
           ref_criteria(iCrit,iBLK) = maxval(ds2*sqrt( &
                (gradY_Uz(1:nI,1:nJ,1:nK)-gradZ_Uy(1:nI,1:nJ,1:nK))**2 + &
                (gradZ_Ux(1:nI,1:nJ,1:nK)-gradX_Uz(1:nI,1:nJ,1:nK))**2 + &
                (gradX_Uy(1:nI,1:nJ,1:nK)-gradY_Ux(1:nI,1:nJ,1:nK))**2))
        case('curlB','curlb')
           ! Curl of magnetic field (current)
           call grad1D(1, iBLK, Bx_G, &
                gradX_Bx,gradY_Bx,gradZ_Bx,"none",Bx_)
           call grad1D(1, iBLK, By_G, &
                gradX_By,gradY_By,gradZ_By, "none",By_)
           call grad1D(1, iBLK, Bz_G, &
                gradX_Bz,gradY_Bz,gradZ_Bz, "none",Bz_)
           ref_criteria(iCrit,iBLK) = maxval(ds2*sqrt( &
                (gradY_Bz(1:nI,1:nJ,1:nK)-gradZ_By(1:nI,1:nJ,1:nK))**2 + &
                (gradZ_Bx(1:nI,1:nJ,1:nK)-gradX_Bz(1:nI,1:nJ,1:nK))**2 + &
                (gradX_By(1:nI,1:nJ,1:nK)-gradY_Bx(1:nI,1:nJ,1:nK))**2))
        case('divU','divu','divV','divv')
           ! Divergence of velocity.
           call grad1D(1, iBLK, RhoUx_G/Rho_G, &
                gradX_Ux,gradY_Ux,gradZ_Ux, "none",Ux_)
           call grad1D(1, iBLK, RhoUy_G/Rho_G, &
                gradX_Uy,gradY_Uy,gradZ_Uy, "none",Uy_)
           call grad1D(1, iBLK, RhoUz_G/Rho_G, &
                gradX_Uz,gradY_Uz,gradZ_Uz, "none",Uz_)
    
           ref_criteria(iCrit,iBLK) = maxval(ds2*abs( &
                gradX_Ux(1:nI,1:nJ,1:nK) + &
                gradY_Uy(1:nI,1:nJ,1:nK) + &
                gradZ_Uz(1:nI,1:nJ,1:nK)))
        case('divB','divb')
           ! Divergence of magnetic field.
           ref_criteria(iCrit,iBLK) = maxval(abs(DivB1_GB(1:nI,1:nJ,1:nK,iBLK)),&
                MASK=true_cell(1:nI,1:nJ,1:nK,iBLK))
        case('Rcurrents','rcurrents')	
           ! Inverse distance from Rcurrents, squared
           outVAR(1:nI,1:nJ,1:nK) = cOne/((max(cTiny, &
                abs(Rcurrents-R_BLK(1:nI,1:nJ,1:nK,iBLK))))**2)
           ref_criteria(iCrit,iBLK) = maxval(outVAR(1:nI,1:nJ,1:nK),&
                MASK=true_cell(1:nI,1:nJ,1:nK,iBLK))
        case('Transient','transient')	
           xxx = cHalf*(x_BLK(nI,nJ,nK,iBLK)+x_BLK(1,1,1,iBLK))
           yyy = cHalf*(y_BLK(nI,nJ,nK,iBLK)+y_BLK(1,1,1,iBLK))
           zzz = cHalf*(z_BLK(nI,nJ,nK,iBLK)+z_BLK(1,1,1,iBLK))
           RR = sqrt(xxx**2+yyy**2+zzz**2 )
           if (UseSunEarth) then
              UseSwitchAMR = (RR.gt.RcritAMR)
           else
              UseSwitchAMR = (RR.gt.RcritAMR).and.(abs(zzz).le.dz_BLK(iBLK))
           end if
           if (UseSwitchAMR) then
              !\
              ! Use dynamic refinement if there is a transient event 
              !/
              call trace_transient(iCrit,iBLK,AMRsort_1)
              ref_criteria(iCrit,iBLK) = AMRsort_1*ds2
              !\              
              ! Restrict the refinement to the particular ray Sun-Earth only
              ! Only if UseSunEarth == .true.
              !/
              if (UseSunEarth) then
                 call refine_sun_earth_cyl(iBLK,xxx,yyy,zzz,AMRsort_2)
              else
                 AMRsort_2 = cOne
              end if
              ref_criteria(iCrit,iBLK) = AMRsort_2*ref_criteria(iCrit,iBLK)
           end if
        case default
           if (UseUserAMR) then
              IsFound=.false.
              call user_amr_criteria(iBLK, userCriteria, RefineCrit(iCrit), IsFound)
              if (IsFound) then
                 ref_criteria(iCrit,iBLK) = userCriteria
              else            
                 write(*,*) 'User refinement criteria not found in user_amr_criteria:', &
                            RefineCrit(iCrit)
                 call stop_mpi('Fix user_amr_criteria or PARAM.in!')
              end if
           else
              call stop_mpi('Unknown RefineCrit='//RefineCrit(iCrit))
           end if
        end select
     end do ! iCrit
  end do ! iBLK
contains
  subroutine trace_transient(iCrit,iBLK,refine_crit)
    use ModAMR,      ONLY:TypeTransient_I,nRefineCrit,RefineCrit


    integer, intent(in) :: iBLK,iCrit
    real, intent(out) :: refine_crit
    real :: AMRsort
    real, dimension(1:nI,1:nJ,1:nK) :: scrARR

    real, dimension(1:nI, 1:nJ, 1:nK) :: RhoOld_C, RhoUxOld_C, &
       RhoUyOld_C, RhoUzOld_C, BxOld_C, ByOld_C, BzOld_C, POld_C    

     do k=1,nK; do j=1,nJ; do i=1,nI
        RhoOld_C(i,j,k)  = StateOld_VCB(rho_,i,j,k,iBLK)
        RhoUxOld_C(i,j,k)= StateOld_VCB(rhoUx_,i,j,k,iBLK)
        RhoUyOld_C(i,j,k)= StateOld_VCB(rhoUy_,i,j,k,iBLK)
        RhoUzOld_C(i,j,k)= StateOld_VCB(rhoUz_,i,j,k,iBLK)
        BxOld_C(i,j,k)   = StateOld_VCB(Bx_,i,j,k,iBLK)
        ByOld_C(i,j,k)   = StateOld_VCB(By_,i,j,k,iBLK)
        BzOld_C(i,j,k)   = StateOld_VCB(Bz_,i,j,k,iBLK)
        POld_C(i,j,k)    = StateOld_VCB(P_,i,j,k,iBLK)
     end do; end do; end do
    select case(TypeTransient_I(iCrit))
    case('P_dot','p_dot')
       !\
       ! refine_crit = abs(|p|-|p|_o)/max(|p|,|p|_o,cTiny)
       ! over all the cells of block iBLK
       !/
       scrARR(1:nI,1:nJ,1:nK) = abs(P_G(1:nI,1:nJ,1:nK) - POld_C(1:nI,1:nJ,1:nK))
       scrARR(1:nI,1:nJ,1:nK) = scrARR(1:nI,1:nJ,1:nK) / max(cTiny,P_G(1:nI,1:nJ,1:nK), &
            POld_C(1:nI,1:nJ,1:nK))
       refine_crit = maxval(scrARR)
    case('T_dot','t_dot')
       !\
       ! refine_crit = abs(|T|-|T|_o)/max(|T|,|T|_o,cTiny)
       ! over all the cells of block iBLK
       !/
       scrARR(1:nI,1:nJ,1:nK) = abs(P_G(1:nI,1:nJ,1:nK)/Rho_G(1:nI,1:nJ,1:nK)  - &
            POld_C(1:nI,1:nJ,1:nK)/RhoOld_C(1:nI,1:nJ,1:nK))
       scrARR(1:nI,1:nJ,1:nK) = scrARR(1:nI,1:nJ,1:nK) / max(cTiny,P_G(1:nI,1:nJ,1:nK)/ &
            Rho_G(1:nI,1:nJ,1:nK),POld_C(1:nI,1:nJ,1:nK)/ &
            RhoOld_C(1:nI,1:nJ,1:nK))
       refine_crit = maxval(scrARR)
    case('Rho_dot','rho_dot')
       !\
       ! refine_crit = abs(|rho|-|rho|_o)/max(|rho|,|rho|_o,cTiny)
       ! over all the cells of block iBLK
       !/
       scrARR(1:nI,1:nJ,1:nK) = abs(Rho_G(1:nI,1:nJ,1:nK) - RhoOld_C(1:nI,1:nJ,1:nK))
       scrARR(1:nI,1:nJ,1:nK) = scrARR(1:nI,1:nJ,1:nK) / max(cTiny,Rho_G(1:nI,1:nJ,1:nK), &
            RhoOld_C(1:nI,1:nJ,1:nK))
       refine_crit = maxval(scrARR)
    case('RhoU_dot','rhoU_dot','rhou_dot')
       !\
       ! refine_crit = abs(|rhoU|-|rhoU|_o)/max(|rhoU|,|rhoU|_o,cTiny)
       ! over all the cells of block iBLK
       !/
       scrARR(1:nI,1:nJ,1:nK) = abs(sqrt(RhoUx_G(1:nI,1:nJ,1:nK)**2              + &
            RhoUy_G(1:nI,1:nJ,1:nK)**2 + RhoUz_G(1:nI,1:nJ,1:nK)**2)      - &
            sqrt(RhoUxOld_C(1:nI,1:nJ,1:nK)**2                                    + &
            RhoUyOld_C(1:nI,1:nJ,1:nK)**2 + RhoUzOld_C(1:nI,1:nJ,1:nK)**2))
       scrARR(1:nI,1:nJ,1:nK) = scrARR(1:nI,1:nJ,1:nK) / max(cTiny, &
            sqrt(RhoUx_G(1:nI,1:nJ,1:nK)**2 + RhoUy_G(1:nI,1:nJ,1:nK)**2  + &
            RhoUz_G(1:nI,1:nJ,1:nK)**2),sqrt(RhoUxOld_C(1:nI,1:nJ,1:nK)**2 + &
            RhoUyOld_C(1:nI,1:nJ,1:nK)**2 + RhoUzOld_C(1:nI,1:nJ,1:nK)**2))
       refine_crit = maxval(scrARR)
    case('B_dot','b_dot')
       !\
       ! refine_crit = abs(|B|-|B|_o)/max(|B|,|B|_o,cTiny)
       ! over all the cells of block iBLK
       !/
       scrARR(1:nI,1:nJ,1:nK) = abs(sqrt(Bx_G(1:nI,1:nJ,1:nK)**2           + &
            By_G(1:nI,1:nJ,1:nK)**2 + Bz_G(1:nI,1:nJ,1:nK)**2)      - &
            sqrt(BxOld_C(1:nI,1:nJ,1:nK)**2                                 + &
            ByOld_C(1:nI,1:nJ,1:nK)**2 + BzOld_C(1:nI,1:nJ,1:nK)**2))
       scrARR(1:nI,1:nJ,1:nK) = scrARR(1:nI,1:nJ,1:nK) / max(cTiny, &
            sqrt(Bx_G(1:nI,1:nJ,1:nK)**2 + By_G(1:nI,1:nJ,1:nK)**2  + &
            Bz_G(1:nI,1:nJ,1:nK)**2),sqrt(BxOld_C(1:nI,1:nJ,1:nK)**2 + &
            ByOld_C(1:nI,1:nJ,1:nK)**2 + BzOld_C(1:nI,1:nJ,1:nK)**2))
       refine_crit = maxval(scrARR)
    case('AK47')
       scrARR(1:nI,1:nJ,1:nK) = abs(sqrt(RhoUx_G(1:nI,1:nJ,1:nK)**2              + &
            RhoUy_G(1:nI,1:nJ,1:nK)**2 + RhoUz_G(1:nI,1:nJ,1:nK)**2)      - &
            sqrt(RhoUxOld_C(1:nI,1:nJ,1:nK)**2                                    + &
            RhoUyOld_C(1:nI,1:nJ,1:nK)**2 + RhoUzOld_C(1:nI,1:nJ,1:nK)**2))
       scrARR(1:nI,1:nJ,1:nK) = scrARR(1:nI,1:nJ,1:nK) / max(cTiny, &
            sqrt(RhoUx_G(1:nI,1:nJ,1:nK)**2 + RhoUy_G(1:nI,1:nJ,1:nK)**2  + &
            RhoUz_G(1:nI,1:nJ,1:nK)**2),sqrt(RhoUxOld_C(1:nI,1:nJ,1:nK)**2 + &
            RhoUyOld_C(1:nI,1:nJ,1:nK)**2 + RhoUzOld_C(1:nI,1:nJ,1:nK)**2))
       AMRsort = maxval(scrARR)
       
       scrARR(1:nI,1:nJ,1:nK) = abs(sqrt(Bx_G(1:nI,1:nJ,1:nK)**2           + &
            By_G(1:nI,1:nJ,1:nK)**2 + Bz_G(1:nI,1:nJ,1:nK)**2)      - &
            sqrt(BxOld_C(1:nI,1:nJ,1:nK)**2                                 + &
            ByOld_C(1:nI,1:nJ,1:nK)**2 + BzOld_C(1:nI,1:nJ,1:nK)**2))
       scrARR(1:nI,1:nJ,1:nK) = scrARR(1:nI,1:nJ,1:nK) / max(cTiny, &
            sqrt(Bx_G(1:nI,1:nJ,1:nK)**2 + By_G(1:nI,1:nJ,1:nK)**2  + &
            Bz_G(1:nI,1:nJ,1:nK)**2),sqrt(BxOld_C(1:nI,1:nJ,1:nK)**2 + &
            ByOld_C(1:nI,1:nJ,1:nK)**2 + BzOld_C(1:nI,1:nJ,1:nK)**2))
       refine_crit = AMRsort*maxval(scrARR)
    case('Rho_2nd_1')
       scrARR(1:nI,1:nJ,1:nK) = ( & 
            abs(Rho_G(0:nI-1,1:nJ,1:nK) + Rho_G(2:nI+1,1:nJ,1:nK) - &
            2 * Rho_G(1:nI,1:nJ,1:nK))                                + &
            abs(Rho_G(1:nI,0:nJ-1,1:nK) + Rho_G(1:nI,2:nJ+1,1:nK) - &
            2 * Rho_G(1:nI,1:nJ,1:nK))                                + &
            abs(Rho_G(1:nI,1:nJ,0:nK-1) + Rho_G(1:nI,1:nJ,2:nK+1) - &
            2 * Rho_G(1:nI,1:nJ,1:nK)))                               / &
            Rho_G(1:nI,1:nJ,1:nK)
       refine_crit = maxval(scrARR)
    case('Rho_2nd_2')
       scrARR(1:nI,1:nJ,1:nK) = abs( & 
            (Rho_G(0:nI-1,1:nJ,1:nK) + Rho_G(2:nI+1,1:nJ,1:nK) - &
            2 * Rho_G(1:nI,1:nJ,1:nK))                             + &
            (Rho_G(1:nI,0:nJ-1,1:nK) + Rho_G(1:nI,2:nJ+1,1:nK) - &
            2 * Rho_G(1:nI,1:nJ,1:nK))                             + &
            (Rho_G(1:nI,1:nJ,0:nK-1) + Rho_G(1:nI,1:nJ,2:nK+1) - &
            2 * Rho_G(1:nI,1:nJ,1:nK)))                            / &
            Rho_G(1:nI,1:nJ,1:nK)
       refine_crit = maxval(scrARR)
    case default
       call stop_mpi('Unknown TypeTransient_I='//TypeTransient_I(iCrit))
    end select
  
  end subroutine trace_transient
end subroutine amr_criteria


subroutine refine_sun_earth_cone(iBLK,xBLK,yBLK,zBLK,refine_profile)
  use ModMain,     ONLY:BLKtest
  use ModProcMH,   ONLY:iProc
  use ModPhysics,  ONLY:Rbody,xEarth,yEarth,zEarth,InvD2Ray
  use ModNumConst, ONLY:cOne,cHalf,cZero,cTiny,cPi
  implicit none
  
  integer, intent(in) :: iBLK
  real, intent(in) :: xBLK,yBLK,zBLK
  real, intent(out) :: refine_profile
  real, parameter:: cOneighty=180.0000000000000000000000000000000
  real :: rBLK, xxx,yyy,zzz, cutFACT
  real :: signY,cosPHI,cosTHETA
  real :: signY_BLK,cosPHI_BLK,cosTHETA_BLK
  !\
  ! This module aims to restrict the refinement mainly along the ray 
  ! Sun-Earth, in a cone with a user-defined opening angle.
  !/
  cutFACT = InvD2Ray*2*cOneighty/cPi  
  !\
  ! For InvD2Ray = 1. ==> refine_profile = 0.174587 at angle 5deg around the ray.
  ! For InvD2Ray = 2. ==> refine_profile = 0.030481 at angle 5deg around the ray.
  !/
  xxx = xEarth
  yyy = yEarth
  zzz = zEarth
  
  if (yyy == cZero) then
     signY = cOne
  else
     signY = abs(yyy)/yyy
  end if
  cosTHETA = zzz/sqrt(xxx**2+yyy**2+zzz**2)
  cosPHI   = xxx/sqrt(xxx**2+yyy**2)
  if ((iProc==0).and.(iBLK==BLKtest)) then
     write(*,*) ''
     write(*,*) '>>>>>>>>>>>>>>>>>>>                  <<<<<<<<<<<<<<<<<<<<'
     write(*,*) '                 Position of the Earth'
     write(*,*) '' 
     write(*,*) 'cosPHI   =',cosPHI
     write(*,*) 'PHI      =',acos(cosPHI)*cOneighty/cPi
     write(*,*) 'cosTHETA =',cosTHETA
     write(*,*) 'THETA    =',acos(cosTHETA)*cOneighty/cPi
     write(*,*) '' 
     write(*,*) '>>>>>>>>>>>>>>>>>>>                  <<<<<<<<<<<<<<<<<<<<<'
     write(*,*) '' 
  end if
  rBLK = sqrt(xBLK**2+yBLK**2+zBLK**2)
  if (rBLK.gt.Rbody) then
     if (yBLK == cZero) then
        signY_BLK = cOne
     else
        signY_BLK = abs(yBLK)/yBLK
     end if
     cosTHETA_BLK = zBLK/sqrt(xBLK**2+yBLK**2+zBLK**2)
     cosPHI_BLK   = xBLK/sqrt(xBLK**2+yBLK**2)
     
     refine_profile = abs(0.5*(signY+signY_BLK)* &
          exp(-cutFACT*(acos(cosPHI_BLK)-acos(cosPHI))**2)* &
          exp(-cutFACT*(acos(cosTHETA_BLK)-acos(cosTHETA))**2))
  else
     refine_profile = cZero
  end if
  
end subroutine refine_sun_earth_cone

subroutine refine_sun_earth_cyl(iBLK,xBLK,yBLK,zBLK,refine_profile)
  use ModPhysics,  ONLY:Rbody,xEarth,yEarth,zEarth,InvD2Ray
  use ModNumConst, ONLY:cOne,cHalf,cZero,cTiny
  implicit none
  
  integer, intent(in) :: iBLK
  real, intent(in) :: xBLK,yBLK,zBLK
  real, intent(out) :: refine_profile
  real :: rBLK, xxx,yyy,zzz, cutFact
  real :: cosPHI,sinPHI,cosTHETA,sinTHETA
  real :: dist2BLK,yPrimeBLK
  !\
  ! This module aims to restrict the refinement mainly along the ray 
  ! Sun-Earth, in a cylinder with user-defined profile across.
  !/
  cutFact = (InvD2Ray*10)**2
  !\
  ! For InvD2Ray = 1. ==> refine_profile = 0.3679 at distance 0.1*Rsun from the ray
  ! For InvD2Ray = 2. ==> refine_profile = 0.0183 at distance 0.1*Rsun from the ray
  ! For InvD2Ray = 3. ==> refine_profile = 0.0001 at distance 0.1*Rsun from the ray
  !/
  xxx = xEarth
  yyy = yEarth
  zzz = zEarth
  cosTHETA = zzz/ &
             sqrt(xxx**2+yyy**2+zzz**2)
  sinTHETA = sqrt(xxx**2+yyy**2)/ &
             sqrt(xxx**2+yyy**2+zzz**2)
  cosPHI   = xxx/sqrt(xxx**2+yyy**2)
  sinPHI   = yyy/sqrt(xxx**2+yyy**2)

  dist2BLK  = (xBLK*sinPHI-yBLK*cosPHI)**2                + &
              (-xBLK*cosTHETA*cosPHI-yBLK*cosTHETA*sinPHI + &
              zBLK*sinTHETA)**2
  yPrimeBLK = (xBLK*sinTHETA*cosPHI+yBLK*sinTHETA*sinPHI  + &
              zBLK*cosTHETA)
  
  rBLK = sqrt(xBLK**2+yBLK**2+zBLK**2)
  if ((rBLK.gt.Rbody).and.(yPrimeBLK >= 0.0)) then
     if (cutFact*dist2BLK <= 150.0) then
        refine_profile = exp(-cutFact*dist2BLK)
     else
        refine_profile = cZero
     end if
  else
     refine_profile = cZero
  end if
  
end subroutine refine_sun_earth_cyl
