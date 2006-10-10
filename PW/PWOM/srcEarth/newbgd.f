
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE NEWBGD
      use ModCommonVariables



      if (TypeSolver == 'Godunov') then

         UBGNDO=UOXYG(NDIM)
         UBGNHE=UHEL(NDIM)
         UBGNDH=UHYD(NDIM)
         
         TBGNDO=TOXYG(NDIM)
         TBGNHE=THEL(NDIM)
         TBGNDH=THYD(NDIM)

         PBGNDO = (POXYG(nDim-1)+ POXYG(nDim)- POXYG(nDim-2))*.8
         PBGNDH = (PHYD(nDim-1) + PHYD(nDim) - PHYD(nDim-2)) *.8
         PBGNHE = (PHEL(nDim-1) + PHEL(nDim) - PHEL(nDim-2)) *.8
      
         DBGNDO=PBGNDO/RGASO/TBGNDO
         DBGNHE=PBGNHE/RGASHE/TBGNHE
         DBGNDH=PBGNDH/RGASH/TBGNDH
      ENDIF

      if (TypeSolver == 'Rusanov') then

         TBGNDO=TOXYG(NDIM)
         TBGNHE=THEL(NDIM)
         TBGNDH=THYD(NDIM)
         
         TBGNDO2=TBGNDO
         TBGNHE2=TBGNHE
         TBGNDH2=TBGNDH

!alternative pressure condition using free expansion
!      PBGNDO=POXYG(NDIM)*(CellVolume_C(nDim)/CellVolumeTop(1)  
!     &     *                   (Rad(nDim)/(Rad(nDim)+DrBnd))**3.0)
!      PBGNHE=PHEL(NDIM)*(CellVolume_C(nDim)/CellVolumeTop(1)  
!     &     *                   (Rad(nDim)/(Rad(nDim)+DrBnd))**3.0)
!      PBGNDH=PHYD(NDIM)*(CellVolume_C(nDim)/CellVolumeTop(1)  
!     &     *                   (Rad(nDim)/(Rad(nDim)+DrBnd))**3.0)

!      PBGNDO2=PBGNDO
!      PBGNHE2=PBGNHE
!      PBGNDH2=PBGNDH

!      PBGNDO2=PBGNDO*(CellVolumeTop(1)/CellVolumeTop(2)  
!     &     *        ((Rad(nDim)+DrBnd)/(Rad(nDim)+2.0*DrBnd))**3.0)
!      PBGNHE2=PBGNHE*(CellVolumeTop(1)/CellVolumeTop(2)  
!     &     *        ((Rad(nDim)+DrBnd)/(Rad(nDim)+2.0*DrBnd))**3.0)
!      PBGNDH2=PBGNDH*(CellVolumeTop(1)/CellVolumeTop(2)  
!     &     *        ((Rad(nDim)+DrBnd)/(Rad(nDim)+2.0*DrBnd))**3.0)

      ScaleHeightO=2.0*1.380658e-16*(TOXYG(nDim)+Telect(nDim))
     &     /(abs(Gravty(nDim))*xmso)
      ScaleHeightH=2.0*1.380658e-16*(THYD(nDim)+Telect(nDim))
     &     /(abs(Gravty(nDim))*xmsh)
      ScaleHeightHe=2.0*1.380658e-16*(THEL(nDim)+Telect(nDim))
     &     /(abs(Gravty(nDim))*xmshe)

!      ScaleHeightO=1.380658e-16*(TOXYG(nDim)+Telect(nDim))
!     &     /(abs(Gravty(nDim)+4.8e-10*Efield(nDim-1)/xmso)*xmso)
!      ScaleHeightH=1.380658e-16*(THYD(nDim)+Telect(nDim))
!     &     /(abs(Gravty(nDim)+4.8e-10*Efield(nDim-1)/xmsh)*xmsh)
!      ScaleHeightHe=1.380658e-16*(THEL(nDim)+Telect(nDim))
!     &     /(abs(Gravty(nDim)+4.8e-10*Efield(nDim-1)/xmshe)*xmshe)
!      
      PBGNDO=POXYG(NDIM)*exp(-DrBnd/ScaleHeightO)
      PBGNDH=PHYD(NDIM) *exp(-DrBnd/ScaleHeightH)
      PBGNHE=PHEL(NDIM) *exp(-DrBnd/ScaleHeightHe)
      
!      PBGNDO2=PBGNDO*exp(-DrBnd/ScaleHeightO)
!      PBGNDH2=PBGNDH*exp(-DrBnd/ScaleHeightH)
!      PBGNHE2=PBGNHE*exp(-DrBnd/ScaleHeightHe)

      PBGNDO2=PBGNDO
      PBGNHE2=PBGNHE
      PBGNDH2=PBGNDH

      DBGNDO=PBGNDO/RGASO/TBGNDO
      DBGNHE=PBGNHE/RGASHE/TBGNHE
      DBGNDH=PBGNDH/RGASH/TBGNDH
      
      DBGNDO2=PBGNDO2/RGASO/TBGNDO2
      DBGNHE2=PBGNHE2/RGASHE/TBGNHE2
      DBGNDH2=PBGNDH2/RGASH/TBGNDH2         

      UBGNDO=DOXYG(nDim)/DBGNDO*AR12TOP(1)/AR23(nDim)*UOXYG(NDIM)
      UBGNHE=DHEL(nDim) /DBGNHE*AR12TOP(1)/AR23(nDim)*UHEL(NDIM)
      UBGNDH=DHYD(nDim) /DBGNDH*AR12TOP(1)/AR23(nDim)*UHYD(NDIM)
      
      UBGNDO2=DBGNDO/DBGNDO2*AR12Top(2)/AR23Top(1)*UBGNDO
      UBGNHE2=DBGNHE/DBGNHE2*AR12Top(2)/AR23Top(1)*UBGNHE
      UBGNDH2=DBGNDH/DBGNDH2*AR12Top(2)/AR23Top(1)*UBGNDH
      
!      UBGNDO2=UBGNDO
!      UBGNHE2=UBGNHE
!      UBGNDH2=UBGNDH
     
      endif
      


C     
      XHTM=1.+ELFXIN*EXP(-(TIME-300.)**2/2./150./150.)
C     
C     Top boundary electron temperature is calculated from
C     Te(nDim) + Tgradient*Dx where Tgradient = Etop/kappa
C     and Etop is the electron heat flux or proportional to that
      TBGNDE=TELECT(NDIM)+XHTM*ETOP/TELECT(NDIM)**2.5



!      write(*,*) 'O'
!      write(*,*) 'd', doxyg(nDim),dbgndo,dbgndo2
!      write(*,*) 'u', uoxyg(nDim),ubgndo,ubgndo2
!      write(*,*) 'p', poxyg(nDim),pbgndo,pbgndo2
!      write(*,*) 'B', ubgndo, doxyg(nDim)/DBGNDO, AR12Top(1)/AR23(nDim),uoxyg(nDim) 
!      write(*,*) AR12Top(1)/AR23(nDim),AR12Top(1),AR23(nDim)
!     
!      write(*,*) 'H'
!      write(*,*) 'd', dhyd(nDim),dbgndh,dbgndh2
!      write(*,*) 'u', uhyd(nDim),ubgndh,ubgndh2
!      write(*,*) 'p', phyd(nDim),pbgndh,pbgndh2
!      write(*,*) 'B', ubgndh, dhyd(nDim)/DBGNDh, AR12Top(1)/AR23(nDim),uhyd(nDim)
!
!      write(*,*) 'He'
!      write(*,*) 'd', dhel(nDim),dbgnhe,dbgnhe2
!      write(*,*) 'u', uhel(nDim),ubgnhe,ubgnhe2
!      write(*,*) 'p', phel(nDim),pbgnhe,pbgnhe2
!      write(*,*) 'B', ubgnhe, dhel(nDim)/DBGNDO, AR12Top(1)/AR23(nDim),uhel(nDim)




C n_e   = n_H+  + n_O+ + n_He+
C rho_e = rho_H*(m_e/m_H) + rho_O*(m_e/m_O) + rho_He*(m_e/m_He)



      DBGNDE=RTHDEL*DBGNDH+RTOXEL*DBGNDO+RTHEEL*DBGNHE
      PBGNDE=RGASE*TBGNDE*DBGNDE

cAlex      write(*,*)'PBGNDE=RGASE*TBGNDE*DBGNDE:',
cAlex     &     PBGNDE,RGASE,TBGNDE,DBGNDE

!      CURHLP=EXP(-0.5*(TIME-CURTIM0)**2/CURTIM**2)
      UBGNDE=(RTHDEL*DBGNDH*UBGNDH+RTOXEL*DBGNDO*UBGNDO+
     &     RTHEEL*DBGNHE*UBGNHE-1.8965E-18*CURRMX)/DBGNDE

CALEX these are sound speeds      
      WBGNDO=SQRT(GAMMA*PBGNDO/DBGNDO)
      WBGNHE=SQRT(GAMMA*PBGNHE/DBGNHE)
      WBGNDH=SQRT(GAMMA*PBGNDH/DBGNDH)
      WBGNDE=SQRT(GAMMA*PBGNDE/DBGNDE)

cAlex write(*,*) 'WBGNDE=SQRT(GAMMA*PBGNDE/DBGNDE):',
cAlex     &     WBGNDE, GAMMA, PBGNDE, DBGNDE

      TCBGO=HLPO*(DBGNDO/DBGNDE)*TBGNDO**2.5
      TCBGE=HLPE*TBGNDE**2.5
      TCBGH=HLPH*(DBGNDH/DBGNDE)*TBGNDH**2.5
      THLP=SQRT(TBGNDH)*(1.-0.047*ALOG10(TBGNDH))**2
      TKLP=DBGNDO/(TBGNDO+16.*TBGNDH)**1.5
      TCBGH=TCBGH/(1.+(0.7692*(CZHN2+CZHO2)+1.0962*CZHO*THLP)/CZHOX
     $     /TKLP)
      THLP=0.5*(XTNMAX+TBGNHE)
      THLP=SQRT(THLP)*(1.-0.093*ALOG10(THLP))**2
      TKLP=DBGNDO/(TBGNDO+4.*TBGNHE)**1.5
      TLLP=DBGNDH/(TBGNHE+4.*TBGNDH)**1.5
      TCBGHE=HLPHE*(DBGNHE/XMSHE)*TBGNHE
      TCBGHE=TCBGHE/(0.99*CZHEN2+0.99*CZHEO2+1.02*CZHEO+1.48*CZHEHE*
     $     THLP+2.22*CZHEH+1.21*CZHEOX*TKLP+2.23*CZHEHD*TLLP)      

      RETURN
      END
