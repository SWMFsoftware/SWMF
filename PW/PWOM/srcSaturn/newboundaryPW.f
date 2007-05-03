 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE NEWBGD
      use ModCommonVariables
      use ModGmPressure
      use ModPWOM, ONLY: iLine
      implicit none

      real :: ScaleHeightO,ScaleHeightH
      real :: xhtm
C     
C     
C     
C     
!******************************************************************************
!     Set the Temperature in the boundary using a isothermal assumption
!******************************************************************************

      TBGNDO=TOXYG(NDIM)
      TBGNHE=THEL(NDIM)
      TBGNDH=THYD(NDIM)

      TBGNDO2=TBGNDO
      TBGNHE2=TBGNHE
      TBGNDH2=TBGNDH
C     
      XHTM=1.+ELFXIN*EXP(-(TIME-300.)**2/2./150./150.)
C     
C     Top boundary electron temperature is calculated from
C     Te(nDim) + Tgradient*Dx where Tgradient = Etop/kappa
C     and Etop is the electron heat flux or proportional to that
      
      TBGNDE=TELECT(NDIM)+XHTM*ETOP/TELECT(NDIM)**2.5

C**********************************************************************
C Set the pressure boundary conditions
C**********************************************************************

! Use Plasma ScaleHeight=k(Ti+Te)/mg
      ScaleHeightO=1.380658e-16*(TOXYG(nDim)+Telect(nDim))
     &     /(abs(Gravty(nDim))*xmso)
      ScaleHeightH=1.380658e-16*(THYD(nDim)+Telect(nDim))
     &     /(abs(Gravty(nDim))*xmsh)

!      ScaleHeightO=2.0*1.380658e-16*(TOXYG(nDim)+Telect(nDim))
!     &     /(abs(Gravty(nDim))*xmso)
!      ScaleHeightH=2.0*1.380658e-16*(THYD(nDim)+Telect(nDim))
!     &     /(abs(Gravty(nDim))*xmsh)
!      
!      ScaleHeightO=1.380658e-16*(TOXYG(nDim)+Telect(nDim))
!     &     /(abs(Gravty(nDim)+4.8e-10*Efield(nDim-100)/xmso)*xmso)
!      ScaleHeightH=1.380658e-16*(THYD(nDim)+Telect(nDim))
!     &     /(abs(Gravty(nDim)+4.8e-10*Efield(nDim-100)/xmsh)*xmsh)
!      
      PBGNDO=POXYG(NDIM)*exp(-DrBnd/ScaleHeightO)
      PBGNDH=PHYD(NDIM) *exp(-DrBnd/ScaleHeightH)

      PBGNDO2=PBGNDO
      PBGNDH2=PBGNDH

C******************************************************************************
C      Set derived quantities:
C     (Density is set from the ideal gas law using the pressure and temperature
C     in the ghost cell. Velocity is set using a steady state continuity eq.)
C******************************************************************************

      DBGNDO=PBGNDO/RGASO/TBGNDO
      DBGNDH=PBGNDH/RGASH/TBGNDH

      DBGNDO2=PBGNDO2/RGASO/TBGNDO2
      DBGNDH2=PBGNDH2/RGASH/TBGNDH2

      UBGNDO=UOXYG(NDIM)
      UBGNDH=UHYD(NDIM)

      UBGNDO2=UBGNDO
      UBGNDH2=UBGNDH

!      UBGNDO=DOXYG(nDim)/DBGNDO*AR12TOP(1)/AR23(nDim)*UOXYG(NDIM)
!      UBGNDH=DHYD(nDim) /DBGNDH*AR12TOP(1)/AR23(nDim)*UHYD(NDIM)
!
!      UBGNDO2=DBGNDO/DBGNDO2*AR12Top(2)/AR23Top(1)*UBGNDO
!      UBGNDH2=DBGNDH/DBGNDH2*AR12Top(2)/AR23Top(1)*UBGNDH

C n_e   = n_H+  + n_O+ + n_He+
C rho_e = rho_H*(m_e/m_H) + rho_O*(m_e/m_O) + rho_He*(m_e/m_He)

      DBGNDE=RTHDEL*DBGNDH+RTOXEL*DBGNDO+RTHEEL*DBGNHE

cAlex      write(*,*)'DBGNDE=RTHDEL,DBGNDH,RTOXEL,DBGNDO,RTHEEL,DBGNHE=',
cAlex     &     DBGNDE,RTHDEL,DBGNDH,RTOXEL,DBGNDO,RTHEEL,DBGNHE


      PBGNDE=RGASE*TBGNDE*DBGNDE

      if (UseGmToPw) then
         PBGNDO = (TBGNDO / (TBGNDO+TBGNDE))*p_I(iLine)
         DBGNDO = PBGNDO/RGASO/TBGNDO
      endif
cAlex      write(*,*)'PBGNDE=RGASE*TBGNDE*DBGNDE:',
cAlex     &     PBGNDE,RGASE,TBGNDE,DBGNDE

cGabor      CURHLP=EXP(-0.5*(TIME-CURTIM0)**2/CURTIM**2)

      UBGNDE=(RTHDEL*DBGNDH*UBGNDH+RTOXEL*DBGNDO*UBGNDO+
     &     RTHEEL*DBGNHE*UBGNHE-1.8965E-18*CURRMX)/DBGNDE

CALEX these are sound speeds      
CALEX      write(*,*) dbgndo,pbgndo,POXYG(NDIM)
      
      WBGNDO=SQRT(GAMMA*PBGNDO/DBGNDO)
CALEX      WBGNHE=SQRT(GAMMA*PBGNHE/DBGNHE)
      WBGNDH=SQRT(GAMMA*PBGNDH/DBGNDH)
      WBGNDE=SQRT(GAMMA*PBGNDE/DBGNDE)

cAlex write(*,*) 'WBGNDE=SQRT(GAMMA*PBGNDE/DBGNDE):',
cAlex     &     WBGNDE, GAMMA, PBGNDE, DBGNDE
C!      TCBGO=HLPO*TBGNDO**2.5
C!      TCBGE=HLPE*TBGNDE**2.5
C!      TCBGH=HLPH*TBGNDH**2.5
     
      TCBGO=HLPO*(DBGNDO/DBGNDE)*TBGNDO**2.5
      TCBGE=HLPE*TBGNDE**2.5
      TCBGH=HLPH*(DBGNDH/DBGNDE)*TBGNDH**2.5
      TCBGHE=HLPHE*(DBGNHE/XMSHE)*TBGNHE**2.5

      RETURN
      END
