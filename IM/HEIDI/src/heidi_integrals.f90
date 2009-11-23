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
  real                 :: GradBCrossB_VIII(3,nPoint,nR,nT)
  real                 :: GradB_VIII(3,nPoint,nR,nT)
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
          RadialDistance_III,Length_III, dLength_III,GradBCrossB_VIII,GradB_VIII)

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
  real                 :: GradBCrossB_VIII(3,nPoint,nR,nT)
  real                 :: GradB_VIII(3,nPoint,nR,nT)
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
          RadialDistance_III,Length_III, dLength_III,GradBCrossB_VIII,GradB_VIII)

!     write(*,*) GradB_VIII
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
  real                 :: GradBCrossB_VIII(3,nPoint,nR,nT)
  real                 :: GradB_VIII(3,nPoint,nR,nT)
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
          RadialDistance_III,Length_III, dLength_III,GradBCrossB_VIII,GradB_VIII)
     
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
  real                 :: GradBCrossB_VIII(3,nPoint,nR,nT)
  real                 :: GradB_VIII(3,nPoint,nR,nT)
  real                 :: Length_III(nPoint,nR,nT) 
  integer              :: iPhi, iR,iPitch

  !----------------------------------------------------------------------------------
  
  call initialize_b_field(LZ, Phi, nPoint, nR, nT, bFieldMagnitude_III, &
       RadialDistance_III,Length_III, dLength_III,GradBCrossB_VIII,GradB_VIII)
  
  
end subroutine get_B_field

!============================================================
subroutine get_coef(dRdt_IIII ,dPhiDt_IIII, InvRdRdt_IIII,dEdt_IIII,dMudt_III)
  
  use ModHeidiSize,  ONLY: nPoint, nPa, nT, nR,nE
  use ModConst,      ONLY: cTiny  
  use ModHeidiMain,  ONLY: Phi, LZ, mu, EKEV, EBND
  use ModHeidiBACoefficients
  use ModHeidiIO,    ONLY: Kp
  use ModHeidiBField
  
  implicit none

  real, dimension(nPa)              :: bMirror_I,PitchAngle_I
  real, dimension(nPoint,nR,nT)     :: bFieldMagnitude_III,RadialDistance_III, Length_III, b4_III
  real, dimension(nR,nT,nE,nPA)     :: dRdt_IIII ,dPhiDt_IIII, InvRdRdt_IIII,dEdt_IIII
  real, dimension(nR,nT,nPA)        :: dMudt_III
  real, dimension(nPoint,nR,nT,NE)  :: DriftR_IIII, DriftTheta_IIII,DriftPhi_IIII
  real, dimension(3,nPoint,nR,nT)   :: GradBCrossB_VIII, GradB_VIII
  real, dimension(3,nPoint,nR,nT)   :: VDrift_VIII
  real, dimension(nPoint-1,nR,nT)   :: dLength_III
  real                              :: BouncedDriftR,BouncedDriftPhi,BouncedInvRdRdt,BouncedInvRdRdt1 
  real                              :: HalfPathLength,Sb,SecondAdiabInv
  integer                           :: iMirror_I(2)
  integer                           :: iPhi, iR, iPitch,iE
  real                              :: bField_VIII(3,nPoint,nR,nT)
  real                              :: Energy
  real                              :: t2,t3,t4,t5,sinPitch,sin2Pitch,cosPitch,cos2Pitch
  real                              :: I2, CoeffE, CoeffMu
  integer                           :: iPoint
  real, parameter                   :: Me = 7.9e15 ,Re = 6.371e6 
  real                              :: VrConv(nR,nT), A
  !----------------------------------------------------------------------------------

  A=7.05E-6/(1.-0.159*Kp + 0.0093 * Kp**2)**3 

  call initialize_b_field(LZ, Phi, nPoint, nR, nT, bFieldMagnitude_III, &
       RadialDistance_III, Length_III, dLength_III, GradBCrossB_VIII,GradB_VIII)
  
  do iE = 1, nE
     Energy = 1000. * EKEV(iE)

     call get_drift_velocity(nPoint,nR,nT, Energy, bFieldMagnitude_III,VDrift_VIII,GradBCrossB_VIII)

     do iPhi = 1,nT
        do iR =1, nR
           do iPitch =1,nPa
             
              PitchAngle_I(iPitch) = acos(mu(iPitch))
              sinPitch = sin(PitchAngle_I(iPitch))
              sin2Pitch = sinPitch * sinPitch

              cosPitch = cos(PitchAngle_I(iPitch))
              cos2Pitch = cosPitch * cosPitch
              
              call find_mirror_points (nPoint,  PitchAngle_I(iPitch), bFieldMagnitude_III(:,iR,iPhi), &
                   bMirror_I(iPitch),iMirror_I)

              call second_adiabatic_invariant(nPoint, iMirror_I(:), bMirror_I(iPitch), &
                   bFieldMagnitude_III(:,iR,iPhi),dLength_III(:,iR,iPhi), LZ(iR), SecondAdiabInv) 

              call half_bounce_path_length(nPoint, iMirror_I(:),bMirror_I(iPitch),&
                   bFieldMagnitude_III(:,iR,iPhi), dLength_III(:,iR,iPhi), LZ(iR), HalfPathLength,Sb)
                            

              DriftR_IIII(:,iR,iPhi,iE)   = VDrift_VIII(1,:,iR,iPhi) ! Radial Drift
              DriftTheta_IIII(:,iR,iPhi,iE) = VDrift_VIII(2,:,iR,iPhi) ! Theta Drift
              DriftPhi_IIII(:,iR,iPhi,iE) = VDrift_VIII(3,:,iR,iPhi) ! Azimuthal Drift

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
                   DriftPhi_IIII(:,iR,iPhi,iE),BouncedDriftPhi,RadialDistance_III(:,iR,iPhi),BouncedInvRdRdt1 ) 
            
              if (Sb==0.0) Sb = cTiny
              dRdt_IIII  (iR,iPhi,iE,iPitch)   = BouncedDriftR/Sb
              dPhiDt_IIII(iR,iPhi,iE,iPitch)   = BouncedDriftPhi/Sb
              InvRdRdt_IIII(iR,iPhi,iE,iPitch) = BouncedInvRdRdt/Sb
              if (SecondAdiabInv ==0.0) SecondAdiabInv = cTiny

              call get_bounced_VGradB(nPoint, LZ(iR), bFieldMagnitude_III(:,iR,iPhi),&
                   GradB_VIII(1,:,iR,iPhi),bMirror_I(iPitch),iMirror_I(:),&
                   dLength_III(:,iR,iPhi),DriftR_IIII(:,iR,iPhi,iE),&
                   RadialDistance_III(:,iR,iPhi),t2)
              
              call get_bounced_VGradB(nPoint, LZ(iR), bFieldMagnitude_III(:,iR,iPhi),&
                   GradB_VIII(2,:,iR,iPhi),bMirror_I(iPitch),iMirror_I(:),&
                   dLength_III(:,iR,iPhi),DriftTheta_IIII(:,iR,iPhi,iE),&
                   RadialDistance_III(:,iR,iPhi),t4)


              call get_bounced_VGradB(nPoint, LZ(iR), bFieldMagnitude_III(:,iR,iPhi),&
                   GradB_VIII(3,:,iR,iPhi),bMirror_I(iPitch),iMirror_I(:),&
                   dLength_III(:,iR,iPhi),DriftPhi_IIII(:,iR,iPhi,iE),&
                   RadialDistance_III(:,iR,iPhi),t5)


              
              VrConv(iR, iPhi) = -A* cos(Phi(iPhi))*(LZ(iR))**4/&              
                   (RE*RE/ME)        
              
              I2 = SecondAdiabInv*SecondAdiabInv
              CoeffE = Energy*(1.- I2/(I2+sin2Pitch))
              ! Needs the dBdt term too
!              dEdt_IIII(iR,iPhi,iE,iPitch)  = CoeffE * t2/Sb -(&
!                   1./RadialDistance_III(nPoint/2,iR,iPhi))* (2.* Energy*I2)/(I2+sin2Pitch)*(BouncedDriftR/Sb+VrConv(iR,iPhi))+&
!                   CoeffE * t4/Sb + CoeffE * t5/Sb
 

              dEdt_IIII(iR,iPhi,iE,iPitch)  = CoeffE * (t2/Sb+VrConv(iR,iPhi)) -(&
                   1./RadialDistance_III(nPoint/2,iR,iPhi))* (2.* Energy*I2)/(I2+sin2Pitch)*(BouncedDriftR/Sb + VrConv(iR,iPhi))+&
                   CoeffE * t4/Sb + CoeffE * t5/Sb
             
              CoeffMu = -(1.-cos2Pitch)*I2/(cosPitch*(I2 + 1. - cos2Pitch))
 

              !Needs the dBdt term too
             dMudt_III(iR,iPhi,iPitch) = CoeffMu * &
                  ((1./RadialDistance_III(nPoint/2,iR,iPhi))*(BouncedDriftR/Sb+VrConv(iR,iPhi))+&
                   0.5 * (t2/Sb+VrConv(iR,iPhi)) + 0.5 * t4/Sb +0.5* t5/Sb) 

 
           end do
        end do
     end do
     
  end do




end subroutine get_coef
