

      SUBROUTINE Solver (iIon, nCell, DtIn,InitialState_GV,
     &     RhoSource, eCollision, fCollision,
     &     RgasSpecies,
     &     HeatCon_G,OutputState_GV)

      use ModCommonVariables
      use ModPWOM, ONLY: IsFullyImplicit, UseIonHeat
C     Input: state variables in the grid, at bottom and top boundaries, 
C            mass, energy, external heating and force source terms,
C            gas constant normalized for the mass of this species
      integer, intent(in) :: iIon

      real, intent(in):: DtIn, InitialState_GV(-1:nCell+2,4) 
      real, intent(in):: RhoSource(nCell), eCollision(nCell),
     &     fCollision(nCell),
     &     RgasSpecies
      real, intent(in)  :: HeatCon_G(0:nCell+1)
C     Output
      real, intent(out):: OutputState_GV(-1:nCell+2,4)

C     Local:
      integer :: iCell
      REAL uSpecies(MAXGRID),dSpecies(MAXGRID),pSpecies(MAXGRID)
      REAL GradT_F(0:MAXGRID),Conduction_F(0:MAXGRID),Diffusion_C(MAXGRID)
      REAL T_G(0:MAXGRID),eSource_C(MaxGrid)
      
      real ::uSur, dSur, pSur, dBgnd, pBgnd, uBgnd

      REAL DBN1(MaxGrid),UBN1(MaxGrid),PBN1(MaxGrid),
     &     U12(MaxGrid),P12(MaxGrid),D12(MaxGrid),
     &     WL(MaxGrid),WR(MaxGrid),PDIFF(MaxGrid),PSUM(MaxGrid)
      
      logical :: DoTest
      !------------------------------------------------------------------------
      eSource_C(1:nCell) = eCollision(1:nCell)
      
      dSpecies(1:nCell) = InitialState_GV(1:nCell,1)
      uSpecies(1:nCell) = InitialState_GV(1:nCell,2)
      pSpecies(1:nCell) = InitialState_GV(1:nCell,3)
      T_G(0:nCell+1)      = InitialState_GV(0:nCell+1,4)

      dSur = InitialState_GV(0,1)
      uSur = InitialState_GV(0,2)
      pSur = InitialState_GV(0,3)

      dBgnd = InitialState_GV(nCell+1,1)
      uBgnd = InitialState_GV(nCell+1,2)
      pBgnd = InitialState_GV(nCell+1,3)
      
      ! Get heat flow in fully implicit case
      if (IsFullyImplicit .and. UseIonHeat)then
      ! Get Temperature Gradient
         do iCell = 0,nCell
            GradT_F(iCell)=(T_G(iCell+1)-T_G(iCell)) / DrBnd
         enddo
      ! Get facevalues of Kappa*Grad(T)
         do iCell = 0,nCell
            Conduction_F(iCell) = 
     &           (HeatCon_G(iCell+1)+HeatCon_G(iCell)) / 2.0 
     &           * GradT_F(iCell)
         enddo
      ! Get heat diffusion term
         do iCell = 1,nCell
            Diffusion_C(iCell) = RgasSpecies/Gmin1 *
     &      (Ar23(iCell)*Conduction_F(iCell)-Ar12(iCell)*Conduction_F(iCell-1))
     &           / (CellVolume_C(iCell))
         enddo
      ! Add heat flux to energy source
         eSource_C(1:nCell) = eCollision(1:nCell)+Diffusion_C(1:nCell)
      
      endif

!     work on bottom of grid

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

! work on rest of grid     
      DO K=2,nCell
         L=K-1
         USUM=USPECIES(L)+USPECIES(K)
         UDIFF=USPECIES(L)-USPECIES(K)
         PDIFF(K)=PSPECIES(L)-PSPECIES(K)
         PSUM(K)=PSPECIES(L)+PSPECIES(K)
         DSUM=DSPECIES(L)+DSPECIES(K)
         A0=SQRT(GAMMA*PSUM(K)*DSUM/4.)
         P12(K)=(PSUM(K)+A0*UDIFF)/2.
c         write(50,*) USUM,PDIFF(nCell),AO
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
      USUM=USPECIES(nCell)+UBGND
      UDIFF=USPECIES(nCell)-UBGND
      PDIFF(nCell+1)=PSPECIES(nCell)-PBGND
      PSUM(nCell+1)=PSPECIES(nCell)+PBGND
      DSUM=DSPECIES(nCell)+DBGND
      A0=SQRT(GAMMA*PSUM(nCell+1)*DSUM/4.)
      P12(nCell+1)=(PSUM(nCell+1)+A0*UDIFF)/2.
      U12(nCell+1)=(USUM+PDIFF(nCell+1)/A0)/2.
      DD=CVMGP(DSPECIES(nCell),DBGND,U12(nCell+1))
      PP=CVMGP(PSPECIES(nCell),PBGND,U12(nCell+1))
      D12(nCell+1)=DD*(GPL1*P12(nCell+1)+GMIN1*PP)/(GMIN1*P12(nCell+1)
     $     +GPL1*PP)
      WL(nCell+1)=USPECIES(nCell)-A0/DSPECIES(nCell)
      WR(nCell+1)=UBGND+A0/DBGND
      D12(nCell+1)=CVMGP(DSPECIES(nCell),D12(nCell+1),WL(nCell+1))
      U12(nCell+1)=CVMGP(USPECIES(nCell),U12(nCell+1),WL(nCell+1))
      P12(nCell+1)=CVMGP(PSPECIES(nCell),P12(nCell+1),WL(nCell+1))
      U12(nCell+1)=CVMGP(U12(nCell+1),UBGND,WR(nCell+1))
      P12(nCell+1)=CVMGP(P12(nCell+1),PBGND,WR(nCell+1))
      D12(nCell+1)=CVMGP(D12(nCell+1),DBGND,WR(nCell+1))
      DO K=1,nCell+1
         DBN1(K)=U12(K)*D12(K)
         UFL=U12(K)*DBN1(K)
         PFL=GMIN2*U12(K)*UFL
         UBN1(K)=UFL+P12(K)
         PBN1(K)=PFL+GAMMA*U12(K)*P12(K)
      end do
c     
      XHTM=EXP(-(TIME-300.)**2/2./150./150.)
c     

!     Update
      
      OutputState_GV(-1:0,1)=InitialState_GV(-1:0,1)
      OutputState_GV(-1:0,2)=InitialState_GV(-1:0,2)
      OutputState_GV(-1:0,3)=InitialState_GV(-1:0,3)
      OutputState_GV(-1:0,4)=InitialState_GV(-1:0,4)
      
      DO K=1,nCell
         KK=K+1
         OutputState_GV(K,1)=DSPECIES(K)-DtIn/DRBND*(AR23(K)*DBN1(KK)-AR12(K)*DBN1(K))
     $        +DtIn*RhoSource(K)
         OutputState_GV(K,2)=(USPECIES(K)*DSPECIES(K)-DtIn/DRBND*(AR23(K)*UBN1(KK)-AR12(K)*
     $        UBN1(K))+DtIn*(DAREA(K)*PSPECIES(K)+FCOLLISION(K)))/OutputState_GV(K,1)

         OutputState_GV(K,3)=PSPECIES(K)-GMIN2*(OutputState_GV(K,2)**2*OutputState_GV(K,1)-
     $     USPECIES(K)**2*DSPECIES(K))-DtIn/DRBND*(AR23(K)*PBN1(KK)-AR12(K)*PBN1(K))
     $     +DtIn*GMIN1*(eSource_C(K))
         OutputState_GV(K,4)=OutputState_GV(K,3)/RGASSPECIES/OutputState_GV(K,1)
      end do

      OutputState_GV(nCell+1:nCell+2,1)=InitialState_GV(nCell+1:nCell+2,1)
      OutputState_GV(nCell+1:nCell+2,2)=InitialState_GV(nCell+1:nCell+2,2)
      OutputState_GV(nCell+1:nCell+2,3)=InitialState_GV(nCell+1:nCell+2,3)
      OutputState_GV(nCell+1:nCell+2,4)=InitialState_GV(nCell+1:nCell+2,4)

      RETURN
      END
