      SUBROUTINE ELFLDW
      use ModCommonVariables
      REAL PRO(0:MaxGrid),PRHE(0:MaxGrid),PRH(0:MaxGrid),PRE(0:MaxGrid),
     &     PR(0:MaxGrid),EZ(0:MaxGrid)
 
C
C
      PRO(0) =USURFO**2*DSURFO
      PRHE(0)=USURHE**2*DSURHE
      PRH(0) =USURFH**2*DSURFH
      PRE(0) =USURFE**2*DSURFE
      PR(0)  =PSURFE+PRE(0)-RTHDEL*(PSURFH+PRH(0))-RTOXEL*
     &     (PSURFO+PRO(0))-RTHEEL*(PSURHE+PRHE(0))
 1    EZ(0)=0.
      
      EZ(NDIM1)  =0.
      PRO(NDIM1) =UBGNDO**2*DBGNDO
      PRHE(NDIM1)=UBGNHE**2*DBGNHE
      PRH(NDIM1) =UBGNDH**2*DBGNDH
      PRE(NDIM1) =UBGNDE**2*DBGNDE
      PR(NDIM1)  =PBGNDE+PRE(NDIM1)-RTHDEL*(PBGNDH+PRH(NDIM1))
     &     -RTOXEL*(PBGNDO+PRO(NDIM1))-RTHEEL*(PBGNHE+PRHE(NDIM1))
      DO K=1,NDIM
         !PR_species = rho * u^2
         PRO(K) =UOXYG(K)**2*DOXYG(K)
         PRHE(K)=UHEL(K)**2*DHEL(K)
         PRH(K) =UHYD(K)**2*DHYD(K)
         PRE(K) =UELECT(K)**2*DELECT(K)
         !PR = sum( M_e/M_species *(P_thermal(k) + rho(k) * u(k)^2))
         PR(K)  =PELECT(K)+PRE(K)-RTHDEL*(PHYD(K)+PRH(K))
     $        -RTOXEL*(POXYG(K)+PRO(K))-RTHEEL*(PHEL(K)+PRHE(K))

C     dM_e/dt-sum(dM_i/dt * m_e/m_i) :: dM_i/dt = momentum collision term
         EZ(k)=FCLSNE(K)-RTHDEL*FCLSNH(K)-RTOXEL*FCLSNO(K)
     $        -RTHEEL*FCLSHE(K)

C     sum(m_e/m_i * (u_e - u_i) * S_i :: S_i is the mass density source
         EZ(k)=EZ(k)+RTHDEL*(UELECT(K)-UHYD(K))*ADMSH(K)+
     ;        RTOXEL*(UELECT(K)-UOXYG(K))*ADMSO(K)+
     ;        RTHEEL*(UELECT(K)-UHEL(K))*ADMSHE(K)

      enddo
      DZ=0.5/DRBND
      DO k=1,NDIM
         
         !dnom = q*n_e
         DNOM=4.803242E-10*(DHYD(k)/XMSH+DOXYG(k)/XMSO+DHEL(k)/XMSHE)

C      A'/A * [u_e^2*rho_e - sum(m_e/m_i * u_i^2*rho_i)]
         EZ1=-DAREA(k)*(PRE(k)-RTHDEL*PRH(k)-RTOXEL*PRO(k)
     $        -RTHEEL*PRHE(k))

C     E_z = -1/(e*n_e) * d/dz [ sum(m_e/m_s(p_s + u_s^2*rho_s)) 
C           -sum_ions ( m_e/m_i*( (u_e-u_i)*S_i - dM_i/dt ) ) - dM_e/dt ] 
C           -1/(e*n_e) * A'/A * ( sum( m_e/m_s * u_s^2 * rho_s ) )
C note: derivatives are central difference

         EFIELD(K) = ( EZ1 + DZ*(EZ(k+1)-EZ(k-1)) - DZ*(PR(k+1)-PR(k-1)) )
     &        /DNOM
c         EFIELD(K) = -DZ*(PR(k+1)-PR(k-1))/DNOM

 
C upwind stencil
c         EFIELD(K) = ( EZ1 + (EZ(k)-EZ(k-1))/DrBnd - (PR(k)-PR(k-1))/DrBnd )
c     &        /DNOM


      enddo
      EFIELD(NDIM)=2.13E-7*CURRMX/TELECT(NDIM)**1.5

c      Efield(nDim)=Efield(nDim-1)
c      if (time .eq. 0.) then
c         EfieldConstant=efield(nDim-201)
c      endif

c      EfieldConstant=2.85e-12
!      EfieldConstant=0.0
c      
!      do k=1, nDim 
!        
!         EFIELD(k) = EfieldConstant
!         Gravty(k)=0.
!      enddo  

      

      IF (EFIELD(NDIM).LT.0.) EFIELD(NDIM)=0.
      
c     Update the momentum and energy source terms now that the 
c     Efield is calculated.


!      Centrifugal(:)=0.0
      
40    DO 50 K=1,NDIM
         FCLSNO(K)=FCLSNO(K)+ADMSO(K)*UOXYG(K)+DOXYG(K)*
     ;        (Centrifugal(K)+GRAVTY(K)+4.803242E-10*EFIELD(K)/XMSO)
         If (NamePlanet .eq. 'Earth ') then
            FCLSHE(K)=FCLSHE(K)+ADMSHE(K)*UHEL(K)+DHEL(K)*
     ;           (Centrifugal(K)+GRAVTY(K)+4.803242E-10*EFIELD(K)/XMSHE)
         endif
         If (NamePlanet .eq. 'Saturn') then 
            FCLSHE(K)=0.
         endif
         FCLSNH(K)=FCLSNH(K)+ADMSH(K)*UHYD(K)+DHYD(K)*
     ;        (Centrifugal(K)+GRAVTY(K)+4.803242E-10*EFIELD(K)/XMSH)
         ECLSNO(K)=ECLSNO(K)+UOXYG(K)*(FCLSNO(K)-0.5*UOXYG(K)*ADMSO(K))
         ECLSHE(K)=ECLSHE(K)+UHEL(K)*(FCLSHE(K)-0.5*UHEL(K)*ADMSHE(K))
         ECLSNH(K)=ECLSNH(K)+UHYD(K)*(FCLSNH(K)-0.5*UHYD(K)*ADMSH(K))

!!!!!Temporary to turn off sources and forces
!         FCLSNO(K)=0.
!         FCLSNH(K)=0.
!         ECLSNO(K)=0.
!         ECLSNH(K)=0.
 50   CONTINUE
      RETURN
      END
