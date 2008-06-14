!^CFG COPYRIGHT UM
!General routines used by the Godunov schemes with an exact Riemann solver
module ModGodunov
  implicit none
  real,private::GammaHere=1.6666666666666666
  real,private::PStar,UStar !pressure and velocity in the Star Region
  real,private::RhoL, UnL, PL, RhoR, UnR, PR !Rho,Un,P: (L)eft and (R)ight
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
  subroutine  pu_star(GammaL,GammaR)
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
    real,optional,intent(in)::GammaL,GammaR

    real :: FL, FLD, FR, FRD !Pressure (F)unctions and its (D)erivatives 
    real :: POld             !Iterative values for pressure
    real :: UDiff            !Function to be nullified
    integer,parameter::nIterMax=10
    integer::iIter
    real,parameter::TolP=0.0010,ChangeStart=2.0*TolP 
    real::Change
    if(present(GammaL))then
       CL=sqrt(GammaL*PL/RhoL)
    else
       CL=sqrt(GammaHere*PL/RhoL)
    end if
    if(present(GammaR))then
       CR=sqrt(GammaR*PR/RhoR)
    else
       CR=sqrt(GammaHere*PR/RhoR)
    end if
    !
    !     Guessed value for PStar  is computed
    !
    call guess_p(POld)

    UDiff = UnR - UnL
    !Initialize loop:
    Change=ChangeStart; iIter=0
    do while(Change > TolP.and.iIter<nIterMax)

       call pressure_function(FL, FLD, POLD, RhoL, PL, CL, GammaL)
       call pressure_function(FR, FRD, POLD, RhoR, PR, CR, GammaR)
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
  
end module ModGodunov
