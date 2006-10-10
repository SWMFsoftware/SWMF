



      SUBROUTINE Solver (NameVariable, uSpecies,dSpecies,pSpecies,
     &     uSur,dSur,pSur,dBgnd,pBgnd,uBgnd,
     &     RhoSource, eCollision, fCollision, qSpecies,
     &     RgasSpecies,
     &     Density, Velocity, Pressure, Temperature)

      use ModCommonVariables

C     Input: state variables in the grid, at bottom and top boundaries, 
C            mass, energy, external heating and force source terms,
C            gas constant normalized for the mass of this species
      character(len=*), intent(in) :: NameVariable
      REAL uSpecies(MAXGRID),dSpecies(MAXGRID),pSpecies(MAXGRID),
     &     uSur, dSur, pSur, dBgnd, pBgnd, uBgnd,
     &     RhoSource(MAXGRID), eCollision(MAXGRID),
     &     fCollision(MAXGRID), qSpecies(MAXGRID),
     &     RgasSpecies

C      Output
      REAL Density(MaxGrid),Velocity(MaxGrid),Pressure(MaxGrid),
     &     Temperature(MaxGrid) 

C     Local:
      REAL DBN1(MaxGrid),UBN1(MaxGrid),PBN1(MaxGrid),
     &     U12(MaxGrid),P12(MaxGrid),D12(MaxGrid),
     &     WL(MaxGrid),WR(MaxGrid),PDIFF(MaxGrid),PSUM(MaxGrid)

      logical :: DoTest
c      DoTest = index(NameVariable,'O') > 0
c      write(*,*)'godunov solver starting for variable ',NameVariable, DoTest

CLFMO
C
C
CALEX work on bottom of grid

      USUM=USUR+USPECIES(1)
      UDIFF=USUR-USPECIES(1)
      PDIFF(1)=PSUR-PSPECIES(1)
      PSUM(1)=PSUR+PSPECIES(1)
      DSUM=DSUR+DSPECIES(1)
      A0=SQRT(GAMMA*PSUM(1)*DSUM/4.)
      P12(1)=(PSUM(1)+A0*UDIFF)/2.
      U12(1)=(USUM+PDIFF(1)/A0)/2.
      DD=CVMGP(DSUR,DSPECIES(1),U12(1))
      PP=CVMGP(PSUR,PSPECIES(1),U12(1))
      D12(1)=DD*(GPL1*P12(1)+GMIN1*PP)/(GMIN1*P12(1)+GPL1*PP)
      WL(1)=USUR-A0/DSUR
      WR(1)=USPECIES(1)+A0/DSPECIES(1)
      D12(1)=CVMGP(DSUR,D12(1),WL(1))
      U12(1)=CVMGP(USUR,U12(1),WL(1))
      P12(1)=CVMGP(PSUR,P12(1),WL(1))
      U12(1)=CVMGP(U12(1),USPECIES(1),WR(1))
      P12(1)=CVMGP(P12(1),PSPECIES(1),WR(1))
      D12(1)=CVMGP(D12(1),DSPECIES(1),WR(1))

CALEX work on rest of grid     
      DO K=2,NDIM
         L=K-1
         USUM=USPECIES(L)+USPECIES(K)
         UDIFF=USPECIES(L)-USPECIES(K)
         PDIFF(K)=PSPECIES(L)-PSPECIES(K)
         PSUM(K)=PSPECIES(L)+PSPECIES(K)
         DSUM=DSPECIES(L)+DSPECIES(K)
         A0=SQRT(GAMMA*PSUM(K)*DSUM/4.)
         P12(K)=(PSUM(K)+A0*UDIFF)/2.
c         write(50,*) USUM,PDIFF(NDIM),AO
         U12(K)=(USUM+PDIFF(K)/A0)/2.
         DD=CVMGP(DSPECIES(L),DSPECIES(K),U12(K))
         PP=CVMGP(PSPECIES(L),PSPECIES(K),U12(K))
         D12(K)=DD*(GPL1*P12(K)+GMIN1*PP)/(GMIN1*P12(K)+GPL1*PP)
         WL(K)=USPECIES(L)-A0/DSPECIES(L)
         WR(K)=USPECIES(K)+A0/DSPECIES(K)
         D12(K)=CVMGP(DSPECIES(L),D12(K),WL(K))
         U12(K)=CVMGP(USPECIES(L),U12(K),WL(K))
         P12(K)=CVMGP(PSPECIES(L),P12(K),WL(K))
         U12(K)=CVMGP(U12(K),USPECIES(K),WR(K))
         P12(K)=CVMGP(P12(K),PSPECIES(K),WR(K))
         D12(K)=CVMGP(D12(K),DSPECIES(K),WR(K))
      end do
 
CALEX work on top of grid
      USUM=USPECIES(NDIM)+UBGND
      UDIFF=USPECIES(NDIM)-UBGND
      PDIFF(NDIM1)=PSPECIES(NDIM)-PBGND
      PSUM(NDIM1)=PSPECIES(NDIM)+PBGND
      DSUM=DSPECIES(NDIM)+DBGND
      A0=SQRT(GAMMA*PSUM(NDIM1)*DSUM/4.)
      P12(NDIM1)=(PSUM(NDIM1)+A0*UDIFF)/2.
      U12(NDIM1)=(USUM+PDIFF(NDIM1)/A0)/2.
      DD=CVMGP(DSPECIES(NDIM),DBGND,U12(NDIM1))
      PP=CVMGP(PSPECIES(NDIM),PBGND,U12(NDIM1))
      D12(NDIM1)=DD*(GPL1*P12(NDIM1)+GMIN1*PP)/(GMIN1*P12(NDIM1)
     $     +GPL1*PP)
      WL(NDIM1)=USPECIES(NDIM)-A0/DSPECIES(NDIM)
      WR(NDIM1)=UBGND+A0/DBGND
      D12(NDIM1)=CVMGP(DSPECIES(NDIM),D12(NDIM1),WL(NDIM1))
      U12(NDIM1)=CVMGP(USPECIES(NDIM),U12(NDIM1),WL(NDIM1))
      P12(NDIM1)=CVMGP(PSPECIES(NDIM),P12(NDIM1),WL(NDIM1))
      U12(NDIM1)=CVMGP(U12(NDIM1),UBGND,WR(NDIM1))
      P12(NDIM1)=CVMGP(P12(NDIM1),PBGND,WR(NDIM1))
      D12(NDIM1)=CVMGP(D12(NDIM1),DBGND,WR(NDIM1))
      DO K=1,NDIM1
         DBN1(K)=U12(K)*D12(K)
         UFL=U12(K)*DBN1(K)
         PFL=GMIN2*U12(K)*UFL
         UBN1(K)=UFL+P12(K)
         PBN1(K)=PFL+GAMMA*U12(K)*P12(K)
      end do
c     
      XHTM=EXP(-(TIME-300.)**2/2./150./150.)
c     

      if(DoTest)then
!         write(*,*)'OldRho    =',Dspecies(nDim-1:nDim)
!         write(*,*)'DHYD      =',DHYD(nDim-1:nDim)
!         write(*,*)'OldP      =',Pspecies(nDim-1:nDim)
!         write(*,*)'PHYD      =',PHYD(nDim-1:nDim)
!         write(*,*)'THYD      =',THYD(nDim-1:nDim)
!         write(*,*)'Rhoflux   =',DBN1(nDim-1:nDim+1)
!         write(*,*)'RhoUflux  =',UBN1(nDim-1:nDim+1)
!         write(*,*)'eflux     =',PBN1(nDim-1:nDim+1)/gmin1

!         write(*,*)'OldRho    =',Dspecies(1:3)
!         write(*,*)'DHYD      =',DHYD(1:3)
!         write(*,*)'OldP      =',Pspecies(1:3)
!         write(*,*)'PHYD      =',PHYD(1:3)
!         write(*,*)'THYD      =',THYD(1:3)
!         write(*,*)'Rhoflux   =',DBN1(1:3)
!         write(*,*)'RhoUflux  =',UBN1(1:3)
!         write(*,*)'eflux     =',PBN1(1:3)/gmin1

      end if

CALEX Update
      DO K=1,NDIM
         KK=K+1
         DENSITY(K)=DSPECIES(K)-DTR1*(AR23(K)*DBN1(KK)-AR12(K)*DBN1(K))
     $        +DTX1*RhoSource(K)
         VELOCITY(K)=(USPECIES(K)*DSPECIES(K)-DTR1*(AR23(K)*UBN1(KK)-AR12(K)*
     $        UBN1(K))+DTX1*(DAREA(K)*PSPECIES(K)+FCOLLISION(K)))/DENSITY(K)

         PRESSURE(K)=PSPECIES(K)-GMIN2*(VELOCITY(K)**2*DENSITY(K)-
     $     USPECIES(K)**2*DSPECIES(K))-DTR1*(AR23(K)*PBN1(KK)-AR12(K)*PBN1(K))
     $     +DTX1*GMIN1*(ECOLLISION(K)+XHTM*QSPECIES(K)*DSPECIES(K))
         TEMPERATURE(K)=PRESSURE(K)/RGASSPECIES/DENSITY(K)
         
c         write(30,*) PSPECIES(K),GMIN2*(CELLNW(1,K)**2*CELLNW(3,K)-
c     $     USPECIES(K)**2*DSPECIES(K)),DTR1*(AR23(K)*PBN1(KK)-AR12(K)*PBN1(K))
c     $     ,DTX1*GMIN1*(ECOLLISIONO(K)+XHTM*QSPECIES(K)*DSPECIES(K)),CELLNW(2,K)

   
      end do
      RETURN
      END
