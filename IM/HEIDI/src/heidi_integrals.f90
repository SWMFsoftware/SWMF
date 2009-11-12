!\
! Contains the bounced averaged numeric integrals
!/
!==================================================================================

subroutine get_IntegralH(IntegralH_III)

  use ModHeidiSize,  ONLY: nPoint, nPa, nT, nR, io, jo,lo
  use ModConst,      ONLY: cPi,  cTiny
  use ModHeidiInput, ONLY: TypeBField
  use ModHeidiMain,  ONLY: Phi, LZ, mu
  use ModNumConst,   ONLY: cDegToRad
  use ModHeidiBField
  
  implicit none 

  real :: y,x,alpha,beta,a1,a2,a3,a4

  real                 :: HalfPathLength,Sb
  real                 :: bMirror_I(nPa),bMirror
  integer              :: iMirror_I(2)
  real                 :: bFieldMagnitude_III(nPoint,nR,nT)! Magnitude of magnetic field 
  real                 :: RadialDistance_III(nPoint,nR,nT)
  real                 :: dLength_III(nPoint-1,nR,nT)      ! Length interval between i and i+1  
  real                 :: Length_III(nPoint,nR,nT) 
  real                 :: PitchAngle_I(nPa)   
  real, intent(out)    :: IntegralH_III(nPa,nR,nT)
  integer              :: iPhi, iR,iPitch
  !----------------------------------------------------------------------------------

  if (TypeBField == 'analytic') then
     alpha=1.+alog(2.+sqrt(3.))/2./sqrt(3.)
     beta=alpha/2.-cPi*sqrt(2.)/12.
     a1=0.055
     a2=-0.037
     a3=-0.074
     a4=0.056

     do iPhi =1 , nT 
        do iR =1, nR    
           do iPitch =1 , nPa
              PitchAngle_I(iPitch) = acos(mu(iPitch))
              x = cos(PitchAngle_I(iPitch))
              y=sqrt(1-x*x)
              
              IntegralH_III(iPitch,iR, iPhi) =alpha-beta*(y+sqrt(y))+ &
                   a1*y**(1./3.)+a2*y**(2./3.)+ a3*y+a4*y**(4./3.)
           end do
        end do
     end do
  endif


  if (TypeBField == 'numeric') then

     call initialize_b_field(LZ, Phi, nPoint, nR, nT, bFieldMagnitude_III, &
          RadialDistance_III,Length_III, dLength_III)

     do iPhi = 1, nT
        do iR =1, nR
           do iPitch =1, nPa
              PitchAngle_I(iPitch) = acos(mu(iPitch))
              call find_mirror_points (nPoint,  PitchAngle_I(iPitch), bFieldMagnitude_III(:,iR,iPhi), &
                      bMirror_I(iPitch),iMirror_I)
              
              call half_bounce_path_length(nPoint, iMirror_I(:),bMirror_I(iPitch),&
                   bFieldMagnitude_III(:,iR,iPhi), dLength_III(:,iR,iPhi), LZ(iR), HalfPathLength,Sb)

              if (HalfPathLength==0.0) HalfPathLength = cTiny
              IntegralH_III(iPitch, iR, iPhi) = HalfPathLength
           end do
        end do
     end do

  endif

end subroutine get_IntegralH

!============================================================
subroutine get_IntegralI(IntegralI_III)

  use ModHeidiSize,  ONLY: nPoint, nPa, nT, nR
  use ModConst,      ONLY: cPi,  cTiny
  use ModHeidiInput, ONLY: TypeBField
  use ModHeidiMain,  ONLY: Phi, LZ, mu
  use ModHeidiBField
 
  implicit none

  real    :: y,x,alpha,beta,a1,a2,a3,a4
  
  real                 :: SecondAdiabInv
  real                 :: bMirror_I(nPa),bMirror
  integer              :: iMirror_I(2)
  real                 :: bFieldMagnitude_III(nPoint,nR,nT) ! Magnitude of magnetic field 
  real                 :: RadialDistance_III(nPoint,nR,nT)
  real                 :: dLength_III(nPoint-1,nR,nT)      ! Length interval between i and i+1  
  real                 :: Length_III(nPoint,nR,nT) 
  real                 :: PitchAngle_I(nPa)
  real, intent(out)    :: IntegralI_III(nPa,nR,nT)
  integer              :: iPhi, iR,iPitch

  !----------------------------------------------------------------------------------

  if (TypeBField == 'analytic') then
     alpha=1.+alog(2.+sqrt(3.))/2./sqrt(3.)
     beta=alpha/2.-cPi*sqrt(2.)/12.
     a1=0.055
     a2=-0.037
     a3=-0.074
     a4=0.056

     do iPhi = 1, nT
        do iR = 1, nR
           do iPitch = 1, nPa
              PitchAngle_I(iPitch) = acos(mu(iPitch))
              if (PitchAngle_I(iPitch)==0.0) PitchAngle_I(iPitch)=cTiny
              x = cos(PitchAngle_I(iPitch))
              y=sqrt(1.-x*x)
                                          
              IntegralI_III(iPitch,iR,iPhi)=2.*alpha*(1.-y)+2.*beta*y*alog(y)+4.*beta*(y-sqrt(y))+  &
                   3.*a1*(y**(1./3.)-y)+6.*a2*(y**(2./3.)-y)+6.*a4*(y-y**(4./3.))  &
                   -2.*a3*y*alog(y)
           end do
        end do
     end do

  endif


  if (TypeBField == 'numeric') then

     call initialize_b_field(LZ, Phi, nPoint, nR, nT, bFieldMagnitude_III, &
          RadialDistance_III,Length_III, dLength_III)

     do iPhi = 1, nT
        do iR =1, nR
           do iPitch =1, nPa
              PitchAngle_I(iPitch) = acos(mu(iPitch))
              
              call find_mirror_points (nPoint,  PitchAngle_I(iPitch), bFieldMagnitude_III(:,iR,iPhi), &
                   bMirror_I(iPitch),iMirror_I)
              
              call second_adiabatic_invariant(nPoint, iMirror_I(:), bMirror_I(iPitch), &
                   bFieldMagnitude_III(:,iR,iPhi),dLength_III(:,iR,iPhi), LZ(iR), SecondAdiabInv)
              !if (SecondAdiabInv ==0.0) SecondAdiabInv = cTiny
              IntegralI_III(iPitch,iR,iPhi) = SecondAdiabInv           
           end do
        end do
     end do

  endif


end subroutine get_IntegralI
!============================================================
subroutine get_neutral_hydrogen(NeutralHydrogen_III,PitchAngle_I)

  use ModHeidiSize,  ONLY: nPoint, nPa, nT, nR
  use ModConst,      ONLY: cPi,  cTiny
  use ModHeidiInput, ONLY: TypeBField
  use ModHeidiMain,  ONLY: Phi, LZ, mu
  use ModHeidiBField
  use ModHeidiHydrogenGeo

  implicit none

  real                 :: AvgHDensity
  real                 :: bMirror_I(nPa)
  integer              :: iMirror_I(2)
  real                 :: bFieldMagnitude_III(nPoint,nR,nT) ! Magnitude of magnetic field 
  real                 :: RadialDistance_III(nPoint,nR,nT)
  real                 :: dLength_III(nPoint-1,nR,nT)      ! Length interval between i and i+1  
  real                 :: Length_III(nPoint,nR,nT) 
  real                 :: PitchAngle_I(nPa)
  real, intent(out)    :: NeutralHydrogen_III(nPa,nR,nT)
  real                 :: Rho_II(nR,nPoint)
  integer              :: iPhi, iR,iPitch
  real                 :: HalfPathLength,Sb

  !----------------------------------------------------------------------------------
  if (TypeBField == 'numeric') then
     
     call initialize_b_field(LZ, Phi, nPoint, nR, nT, bFieldMagnitude_III, &
          RadialDistance_III,Length_III, dLength_III)
     
     do iPhi = 1, nT
        do iR =1, nR
           do iPitch =1, nPa
              PitchAngle_I(iPitch) = acos(mu(iPitch))
              
              call find_mirror_points (nPoint,  PitchAngle_I(iPitch), bFieldMagnitude_III(:,iR,iPhi), &
                   bMirror_I(iPitch),iMirror_I)

              
              call get_rairden_density(nPoint, nR, LZ(iR), Rho_II(iR,:))
              call get_hydrogen_density(nPoint, LZ(iR), bFieldMagnitude_III(:,iR,iPhi), bMirror_I(iPitch)&
                   ,iMirror_I(:),dLength_III(:,iR,iPhi),Rho_II(iR,:),AvgHDensity)
              
              call half_bounce_path_length(nPoint, iMirror_I(:),bMirror_I(iPitch),&
                   bFieldMagnitude_III(:,iR,iPhi), dLength_III(:,iR,iPhi), LZ(iR), HalfPathLength,Sb)

              if (Sb==0.0) Sb = cTiny
              if (AvgHDensity ==0.0) AvgHDensity = cTiny
              NeutralHydrogen_III(iPitch,iR,iPhi) = AvgHDensity/Sb
           end do
        end do
     end do


  end if

end subroutine get_neutral_hydrogen
!============================================================
subroutine get_B_field(bFieldMagnitude_III)
  
  use ModHeidiSize,  ONLY: nPoint, nPa, nT, nR
  use ModHeidiBField
  use ModHeidiMain,  ONLY: Phi, LZ
  implicit none

  real                 :: bFieldMagnitude_III(nPoint,nR,nT) ! Magnitude of magnetic field 
  real                 :: RadialDistance_III(nPoint,nR,nT)
  real                 :: dLength_III(nPoint-1,nR,nT)       ! Length interval between i and i+1  
  real                 :: Length_III(nPoint,nR,nT) 
  integer              :: iPhi, iR,iPitch

  !----------------------------------------------------------------------------------
  
  call initialize_b_field(LZ, Phi, nPoint, nR, nT, bFieldMagnitude_III, &
       RadialDistance_III,Length_III, dLength_III)
  
  
end subroutine get_B_field

!============================================================
subroutine get_coef
  
  use ModHeidiSize,  ONLY: nPoint, nPa, nT, nR,nE
  use ModConst,      ONLY: cTiny  
  use ModHeidiMain,  ONLY: Phi, LZ, mu,EKEV
  use ModHeidiBACoefficients

  use ModHeidiBField
  
  implicit none

  real, dimension(nPa)              :: bMirror_I,PitchAngle_I
  real, dimension(nPoint,nR,nT)     :: bFieldMagnitude_III,RadialDistance_III, Length_III, b4_III
  real, dimension(nR,nT,nE,nPA)     :: dRdt_IIII ,dPhiDt_IIII, InvRdRdt_IIII
  real, dimension(nPoint,nR,nT,NE)  :: DriftR_IIII, DriftPhi_IIII
  real, dimension(3,nPoint,nR,nT)   :: GradB2over2_VIII
  real, dimension(3,nPoint,nR,nT,NE):: VDrift_VIII
  real, dimension(nPoint-1,nR,nT)   :: dLength_III
  real                              :: BouncedDriftR,BouncedDriftPhi,BouncedInvRdRdt
  real                              :: HalfPathLength,Sb
  integer                           :: iMirror_I(2)
  integer                           :: iPhi, iR, iPitch,iE

  !----------------------------------------------------------------------------------

  call initialize_b_field(LZ, Phi, nPoint, nR, nT, bFieldMagnitude_III, &
       RadialDistance_III,Length_III, dLength_III,GradB2over2_VIII)
  
 
  do iE = 1, nE
     call get_drift_velocity(nPoint,nR,nT, EKEV(iE), b4_III,VDrift_VIII(:,:,:,:,iE))
     
     do iPhi = 1, nT
        do iR =1, nR
           do iPitch =1, nPa
              PitchAngle_I(iPitch) = acos(mu(iPitch))
              
              call find_mirror_points (nPoint,  PitchAngle_I(iPitch), bFieldMagnitude_III(:,iR,iPhi), &
                   bMirror_I(iPitch),iMirror_I)
              
              call half_bounce_path_length(nPoint, iMirror_I(:),bMirror_I(iPitch),&
                   bFieldMagnitude_III(:,iR,iPhi), dLength_III(:,iR,iPhi), LZ(iR), HalfPathLength,Sb)
              
              DriftR_IIII(:,:,:,:) =  VDrift_VIII(1,:,:,:,:)  ! Radial Drift
              DriftPhi_IIII(:,:,:,:) = VDrift_VIII(2,:,:,:,:) ! Azimuthal Drift
              
              ! \
              !                 Calculates :
              !              BouncedDrift = <dR/dt>
              !              BouncedInvRdRdt = <1/R dR/dt>
              !/
              call get_bounced_drift(nPoint, LZ(iR),bFieldMagnitude_III(:,iR,iPhi),&
                   bMirror_I(iPitch), iMirror_I(:),dLength_III(:,iR,iPhi),&
                   DriftR_IIII(:,iR,iPhi,iE),BouncedDriftR,RadialDistance_III(:,iR,iPhi),BouncedInvRdRdt ) 
           
              call get_bounced_drift(nPoint, LZ(iR),bFieldMagnitude_III(:,iR,iPhi),&
                   bMirror_I(iPitch), iMirror_I(:),dLength_III(:,iR,iPhi),&
                   DriftPhi_IIII(:,iR,iPhi,iE),BouncedDriftPhi,RadialDistance_III(:,iR,iPhi)) 
              
              if (Sb==0.0) Sb = cTiny
              dRdt_IIII  (iR,iPhi,iE,iPitch)   = BouncedDriftR/Sb
              dPhiDt_IIII(iR,iPhi,iE,iPitch)   = BouncedDriftPhi/Sb
              InvRdRdt_IIII(iR,iPhi,iE,iPitch) = BouncedInvRdRdt/Sb
                           

           end do
        end do
     end do
  end do

end subroutine get_coef
