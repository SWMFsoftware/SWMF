!^CFG COPYRIGHT UM
!General routines used by the Godunov schemes with an exact Riemann solver
module ModExactRS
  implicit none
  real,private::GammaHere=1.6666666666666666
  real        ::PStar,UnStar !pressure and velocity in the Star Region
  real        ::RhoL, UnL, PL, RhoR, UnR, PR !Rho,Un,P: (L)eft and (R)ight
  real,private::CL,CR       !Sound speeds
  real::WL,WR       !Total perturbation speed 
                            !(advection+speed of sound for weak wave)
                            !Shock wave speed for a shock.
contains
  !========================================================================!
  subroutine get_gamma(GammaIn)
    real,intent(in)::GammaIn
    !Should be applied only if: 
    !FIRST, Gamma=const, and, 
    !SECOND, Gamma=const /= (5/3)
    GammaHere=GammaIn
  end subroutine get_gamma
  !========================================================================!
  subroutine  pu_star(GammaLIn,GammaRIn)
    implicit none
    !     Programer: E. F. Toro                                            *
    !                                                                      *
    !     Last revision: 31st May 1999                                     *
    !                                                                      *
    !     Theory is found in Ref. 1, Chapt. 4 and in original              *
    !     references therein                                               *
    !                                                                      *
    !     1. Toro, E. F., "Riemann Solvers and Numerical                   *
    !                      Methods for Fluid Dynamics"                     *
    !                      Springer-Verlag, 1997                           *
    !                      Second Edition, 1999                            *
    !                                                                      *
    !     This program is part of                                          *
    !                                                                      *
    !     NUMERICA                                                         *
    !     A Library of Source Codes for Teaching,                          *
    !     Research and Applications,                                       *
    !     by E. F. Toro                                                    *
    !     Published by NUMERITEK LTD, 1999                                 *
    !     Website: www.numeritek.com                 
    !     Revision history: re-wrote in f90, allow GammaL/=GammaR
    !     Purpose: to compute the solution for pressure and
    !              velocity in the Star Region
    !
    !
    real,optional,intent(in)::GammaLIn,GammaRIn

    real :: FL, FLD, FR, FRD !Pressure (F)unctions and its (D)erivatives 
    real :: POld             !Iterative values for pressure
    real :: UDiff            !Function to be nullified
    integer,parameter::nIterMax=10
    integer::iIter
    real,parameter::TolP=0.0010,ChangeStart=2.0*TolP 
    real::Change,UVac,GammaL,GammaR
    real,parameter::cSafetyFactor=0.999 !Limits the speed of expansion
    if(present(GammaLIn))then
       GammaL=GammaLIn
    else
       GammaL=GammaHere
    end if
    CL=sqrt(GammaL*PL/RhoL)
    
    if(present(GammaRIn))then
       GammaR=GammaRIn
    else
       GammaR=GammaHere
    end if
    CR=sqrt(GammaR*PR/RhoR)

    UVac=2.0*(CL/(GammaL-1.0)+CR/(GammaR-1.0))

    UDiff = UnR - UnL

    !Prevents the formation of vacuum regions:
    if(UDiff > cSafetyFactor*UVac)call CON_stop('Vacuum in RS')

    !
    !     Guessed value for PStar  is computed
    !
    call guess_p(POld)

    !Initialize loop:
    Change=ChangeStart; iIter=0
    do while(Change > TolP.and.iIter<nIterMax)

       call pressure_function(FL, FLD, POLD, RhoL, PL, CL,WL, GammaLIn)
       call pressure_function(FR, FRD, POLD, RhoR, PR, CR,WR, GammaRIn)
       PStar      = POld - (FL + FR + UDiff)/(FLD + FRD)
       Change= 2.0*abs((PStar - POld)/(PStar + POld))
       POLD  = PStar; iIter=iIter+1
    end do

    !     Compute velocity in Star Region

    UnStar = 0.50*(UnL + UnR + FR - FL)
    WR= UnR + WR
    WL= UnL - WL
  contains
    !====================================================!
    subroutine pressure_function(F,FD,P,RhoK,PK,CK,WK,GammaIn)
      !
      !     Purpose: to evaluate the pressure functions
      !              FL and FR in exact Riemann solver
      !              and their first derivatives 
      !
      real,intent(in)::P,RhoK,PK,CK
      real,intent(out)::F,FD,WK
      real,optional,intent(in)::GammaIn
      real::Gamma
      !Misc:
      REAL::    BK, PRatio,PRatioPowG1, QRT
      !-------------------------------------------!
      if(present(GammaIn))then
         Gamma=GammaIn
      else
         Gamma=GammaHere
      end if
      IF(P.LE.PK)THEN
         !
         !        Rarefaction wave
         !
         PRatio = P/PK 
         PRatioPowG1= PRatio**((GAMMA - 1.0)/(2.0*GAMMA))
         F    = CK*(PRatioPowG1 - 1.0)*2.0/(GAMMA - 1.0)
         FD   = PRatioPowG1/(PRatio*RhoK*CK)
         WK = CK
      ELSE
         !
         !        Shock wave
         !
         BK = PK*(GAMMA - 1.0)/(GAMMA + 1.0)+P
         QRT = SQRT(2.0/((GAMMA + 1.0)*RhoK*BK))
         F   = (P - PK)*QRT
         FD  = (1.0 - 0.5*(P - PK)/BK)*QRT
         WK=0.50*(GAMMA + 1.0)*BK*QRT
      ENDIF
    end subroutine pressure_function
    !===============================================!
    subroutine guess_p(PGuess)

      !     Purpose: to provide a guessed value for pressure
      !              PM in the Star Region. 
      !     Author of the version : P.Voinovich, 1992
      real,intent(out)::PGuess
      real::ImpedanceL,ImpedanceR,PAvr,ImpedanceTotalInv
      real,parameter::cTinyFractionOf=1.0e-8
      !----------------------------------------------------!
      ImpedanceL = RhoL*CL; ImpedanceR = RhoR*CR

      ImpedanceTotalInv=1.0/(ImpedanceL+ImpedanceR)

      PAvr=ImpedanceTotalInv*(PL*ImpedanceR+PR*ImpedanceL)

      PGuess=max(cTinyFractionOf*PAvr,&
           PAvr+(UnL-UnR)*ImpedanceTotalInv*ImpedanceL*ImpedanceR)
    end subroutine guess_p
  end subroutine pu_star
  !==================================================!
  subroutine sample(S, Rho, Un, P, GammaLIn,GammaRIn)
    !
    !     Purpose: to sample the solution throughout the wave
    !             pattern. Pressure and velocity in the
    !              Star Region are known. Sampling is performed
    !              in terms of the 'speed' S = X/T. Sampled
    !              values are Rho, U, P
    !
    real,intent(in)::S
    real,intent(out):: Rho, Un, P
    real,optional,intent(in)::GammaLIn,GammaRIn

    !
    IF(S >  WR)THEN
       !
       !              Sampled point is data state
       !
       Rho = RhoR
       Un = UnR
       P = PR
    ELSEIF(S < WL)THEN
       !
       !              Sampled point is data state
       !
       Rho = RhoL
       Un = UnL
       P = PL
    ELSEIF(S.LE.UnStar)THEN
       !
       !        Sampling point lies to the left of the contact
       !        discontinuity
       call simple_wave(-S+UnStar,RhoL,-UnL+UnStar,PL,CL,GammaLIn)
       Un = UnStar - Un
    ELSE
       !
       !        Sampling point lies to the right of the contact
       !        discontinuity
       !    
       call simple_wave(S-UnStar,RhoR,UnR-UnStar,PR,CR, GammaRIn)
       Un = Un+UnStar
    end IF
  contains
    subroutine simple_wave(SRel,RhoK,UnK,PK,CK, GammaIn)
      real,intent(in)::SRel,RhoK,UnK,PK,CK
      real,optional,intent(in)::GammaIn
      real::Gamma
      !Misc:
      real::C,G7
      if(present(GammaIn))then
         Gamma=GammaIn
      else
         Gamma=GammaHere
      end if
      IF(PStar.LE.PK)THEN
         !
         !           rarefaction
         !
         IF(SRel< CK*(PStar/PK)**((Gamma - 1.0)/(2.0*Gamma)))THEN
            !
            !                 Sampled point is Star  state
            !

            Un = 0.0
            P = PStar
         ELSE
            !
            !                 Sampled point is inside left fan
            !

            G7 = (Gamma - 1.0)/2.0

            !Solve the system of equations:
            !
            ! CK - G7 * UnK = C - G7* Un
            ! SRel = C + Un

            C = (CK - G7*(UnK - SRel))/(1.0+G7)
            Un = SRel - C
            P = PK*(C/CK)**(2.0*Gamma/(Gamma - 1.0))
         ENDIF
         Rho = RhoK*(P/PK)**(1.0/Gamma)
      ELSE
         !
         !          Shock
         ! 
         !
         !              Sampled point is Star state
         !
         Rho= RhoK*((Gamma + 1.0)*PStar + PK*(Gamma - 1.0))/&
              ((Gamma - 1.0)*PStar + PK*(Gamma + 1.0))
         Un = 0.0
         P = PStar
      ENDIF
    end  subroutine simple_wave
  end subroutine sample
  !=====================================================================!
  !Calculates the Godunov flux with the exact Riemann solver. The
  !interface is the most generic (nDim may be 1,2,3), the state vector
  !of nVar .ge. 2+nDim may inslude an arbitrary number of extra state 
  !variable, for all of the the continuity equation being implied. At 
  !the same time, the interface is the least convenient: it requires to
  !pass nDim and nVar and the unity vector normal to the face which is 
  !(/1,0,0/), (/0,1,0/), (0,0,1/) for 3D cartesian grid while calculate 
  !the fluxes along x,y,z axis respectively. 
  subroutine get_godunov_flux(&
                 nDim,Normal_D, & !Unity vector normal to face 
                 nVar,LeftState_V,RightState_V,&!Input states
                 Flux_V, CMax, &  !Output flux and max perturbation speed
                 GammaL,GammaR,&  !Optional Gamma, if needed
                 Energy0L,Energy0R, & !Optional Energy0, if needed
                 DoTest)
    integer,intent(in)::nDim,nVar
    real,intent(in),dimension(nDim)::Normal_D !Unity vector normal to the face
    real,intent(in),dimension(nVar)::LeftState_V,RightState_V
    real,intent(out),dimension(nVar)::Flux_V
    real,intent(out)::CMax
    real,intent(in),optional::GammaL,GammaR !Effective gamma, left and right
    real,intent(in),optional::Energy0L,Energy0R
    logical,optional::DoTest
    real::Rho, Un, P, StateStar_V(nVar)
    real::RhoSide,UnSide,GammaSideM1Inv,E0=0.0
    integer::iVar
    !----------------------------------------
    !Take scalars

    RhoL=LeftState_V(1)   ;RhoR=RightState_V(1)
    PL  =LeftState_V(nVar);PR  =RightState_V(nVar)
    UnL =sum( LeftState_V(2:1+nDim) * Normal_D)
    UnR =sum(RightState_V(2:1+nDim) * Normal_D)

    !Take the parameters at the Contact Discontinuity (CD)

    call pu_star(GammaL,GammaR)

    if(UnStar>0.0)then
       !The CD is to the right from the face
       !The Left gas passes through the face
       RhoSide = RhoL ; UnSide = UnL
       StateStar_V=LeftState_V
       if(present(GammaL))then
          GammaSideM1Inv=1.0/(GammaL-1.0)
       else
          GammaSideM1Inv=1.0/(GammaHere-1.0)
       end if
       if(present(Energy0L))then
          E0=Energy0L
       else
          E0=0.0
       end if
    else
       !The CD is to the left from the face
       !The Right gas passes through the face
       RhoSide = RhoR ; UnSide = UnR
       StateStar_V=RightState_V
       if(present(GammaR))then
          GammaSideM1Inv=1.0/(GammaR-1.0)
       else
          GammaSideM1Inv=1.0/(GammaHere-1.0)
       end if
       if(present(Energy0R))then
          E0=Energy0R
       else
          E0=0.0
       end if
    end if

    !Take the parameters at the face

    call sample(0.0,Rho, Un, P, GammaL,GammaR)
    StateStar_V(1)=Rho
    StateStar_V(2:1+nDim) = StateStar_V(2:1+nDim)+(Un-UnSide)*Normal_D
    StateStar_V(nVar) = P
    do iVar=2+nDim,nVar-1
       StateStar_V(iVar) = StateStar_V(iVar)*(Rho/RhoSide)
    end do

    !Calculate flux (1) take conservative variable

    StateStar_V(nVar) = (1.0 +  GammaSideM1Inv) *  StateStar_V(nVar) &
         + 0.5 * Rho *sum(StateStar_V(2:1+nDim)**2) + Rho*E0
    StateStar_V(2:1+nDim) = StateStar_V(2:1+nDim) * Rho

    ! (2) take advective part of the flux

    Flux_V = StateStar_V * Un 

    ! (3) add the pressure tensor

    Flux_V(2:1+nDim) = Flux_V(2:1+nDim) + P * Normal_D

    CMax = max( WR, -WL)
    if(present(DoTest))then
       if(DoTest)then
          write(*,'(a)')'Test for get_godunov_flux:'
          write(*,*)'Left state:', LeftState_V
          write(*,*)'Right state:',RightState_V
          write(*,*)'Star state:',StateStar_V
          write(*,*)'PStar,UStar=',PStar,UnStar
          write(*,*)'Flux:',Flux_V
          write(*,*)'CMax',CMax
          if(present(GammaL))write(*,*)'GammaL,GammaR=',GammaL,GammaR
       end if
    end if
  end subroutine get_godunov_flux
  !==================================================================!
end module ModExactRS
