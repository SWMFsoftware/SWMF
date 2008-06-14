!^CFG COPYRIGHT UM
!General routines used by the Godunov schemes with an exact Riemann solver
module ModExactRS
  implicit none
  real,private::GammaHere=1.6666666666666666
  real        ::PStar,UStar !pressure and velocity in the Star Region
  real,private::RhoL, UnL, PL, RhoR, UnR, PR !Rho,Un,P: (L)eft and (R)ight
  real,private::GammaL,GammaR
  real,private::CL,CR       !Sound speeds
contains
  !========================================================================!
  subroutine get_gamma(GammaIn)
    real,intent(in)::GammaIn
    !Should be applied only if: 
    !FIRST, Gamma=const, and, 
    !SECOND, Gamma /= (5/3)
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
    real::Change
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
    !
    !     Guessed value for PStar  is computed
    !
    call guess_p(POld)

    UDiff = UnR - UnL
    !Initialize loop:
    Change=ChangeStart; iIter=0
    do while(Change > TolP.and.iIter<nIterMax)

       call pressure_function(FL, FLD, POLD, RhoL, PL, CL, GammaLIn)
       call pressure_function(FR, FRD, POLD, RhoR, PR, CR, GammaRIn)
       PStar      = POld - (FL + FR + UDiff)/(FLD + FRD)
       Change= 2.0*abs((PStar - POld)/(PStar + POld))
       POLD  = PStar; iIter=iIter+1
    end do

    !     Compute velocity in Star Region

    UStar = 0.50*(UnL + UnR + FR - FL)
  contains
    !====================================================!
    subroutine pressure_function(F,FD,P,RhoK,PK,CK,GammaIn)
      !
      !     Purpose: to evaluate the pressure functions
      !              FL and FR in exact Riemann solver
      !              and their first derivatives
      !
      real,intent(in)::P,RhoK,PK,CK
      real,intent(out)::F,FD
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
      ELSE
         !
         !        Shock wave
         !
         BK = PK*(GAMMA - 1.0)/(GAMMA + 1.0)
         QRT = SQRT(2.0/((GAMMA + 1.0)*RhoK*(BK + P)))
         F   = (P - PK)*QRT
         FD  = (1.0 - 0.5*(P - PK)/(BK + P))*QRT
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
  SUBROUTINE SAMPLE(S, Rho, Un, P)
    !
    !     Purpose: to sample the solution throughout the wave
    !             pattern. Pressure and velocity in the
    !              Star Region are known. Sampling is performed
    !              in terms of the 'speed' S = X/T. Sampled
    !              values are Rho, U, P
    !
    real,intent(in)::S
    real,intent(out):: Rho, Un, P
    !Misc:
    real::SHL,CML,STL,SHR,CMR,STR,C,SL,SR,G3,G4,G5,G6,G7
    real::PStarL,PStarR
    !
    IF(S.LE.UStar)THEN
       !
       !        Sampling point lies to the left of the contact
       !        discontinuity

       IF(PStar.LE.PL)THEN
          !
          !           Left rarefaction
          !
          SHL = UnL - CL
          !
          IF(S.LE.SHL)THEN
             !
             !              Sampled point is left data state
             !
             Rho = RhoL
             Un = UnL
             P = PL
          ELSE
             CML = CL*(PStar/PL)**((GammaL - 1.0)/(2.0*GammaL))
             STL = UStar - CML
             IF(S.GT.STL)THEN
                !
                !                 Sampled point is Star Left state
                !
                Rho = RhoL*(PStar/PL)**(1.0/GammaL)
                Un = UStar
                P = PStar
             ELSE
                !
                !                 Sampled point is inside left fan
                !
                G4 = 2.0/(GammaL - 1.0)
                G3 = G4*GammaL
                
                G5 = 2.0/(GammaL + 1.0)    
                G7 = (GammaL - 1.0)/2.0
 
                Un = G5*(CL + G7*UnL + S)
                C = G5*(CL + G7*(UnL - S))
                Rho = RhoL*(C/CL)**G4
                P = PL*(C/CL)**G3
             ENDIF
          ENDIF
       ELSE
          !
          !           Left shock
          !
          PStarL = PStar/PL
          SL  = UnL - CL*SQRT(((GammaL + 1.0)*PStarL + (GammaL - 1.0))/&
                           (2.0*GammaL))  
          !
          IF(S.LE.SL)THEN
             !
             !              Sampled point is left data state
             !
             Rho = RhoL
             Un = UnL
             P = PL
             !
          ELSE
             !
             !              Sampled point is Star Left state
             !
             G6 = (GammaL - 1.0)/(GammaL + 1.0)
             
             Rho = RhoL*(PStarL + G6)/(PStarL*G6 + 1.0)
             Un = UStar
             P = PStar
          ENDIF
       ENDIF
    ELSE
       !
       !        Sampling point lies to the right of the contact
       !        discontinuity
       !
       IF(PStar.GT.PR)THEN
          !
          !           Right shock
          !
          PStarR = PStar/PR
          SR  = UnR + CR*SQRT(((GammaR + 1.0)*PStarR + (GammaR - 1.0))/&
                           (2.0*GammaR))  
          !
          IF(S.GE.SR)THEN
             !
             !              Sampled point is right data state
             !
             Rho = RhoR
             Un = UnR
             P = PR
          ELSE
             !
             !              Sampled point is Star Right state
             !
             G6 = (GammaR - 1.0)/(GammaR + 1.0)
             Rho= RhoR*(PStarR + G6)/(PStarR*G6 + 1.0)
             Un = UStar
             P = PStar
          ENDIF
       ELSE
          !
          !           Right rarefaction
          !
          SHR = UnR + CR
          !
          IF(S.GE.SHR)THEN
             !
             !              Sampled point is right data state
             !
             Rho = RhoR
             Un = UnR
             P = PR
          ELSE
             CMR = CR*(PStar/PR)**((GammaR - 1.0)/(2.0*GammaR))
             STR = UStar + CMR
             !
             IF(S.LE.STR)THEN
                !
                !                 Sampled point is Star Right state
                !
                Rho = RhoR*(PStar/PR)**(1.0/GammaR)
                Un = UStar
                P = PStar
             ELSE
                !
                !                 Sampled point is inside left fan
                !
                G4 = 2.0/(GammaR - 1.0)
                G3 = G4*GammaR
                
                G5 = 2.0/(GammaR + 1.0)    
                G7 = (GammaR - 1.0)/2.0

                Un = G5*(-CR + G7*UnR + S)
                C = G5*(CR - G7*(UnR - S))
                Rho = RhoR*(C/CR)**G4
                P = PR*(C/CR)**G3
             ENDIF
          ENDIF
       ENDIF
    ENDIF
  end SUBROUTINE SAMPLE
end module ModExactRS
