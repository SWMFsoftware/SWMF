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
  real                 :: dBdt_III(nPoint,nR,nT)
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
              y = sqrt(1-x*x)
              
              IntegralH_III(iPitch,iR, iPhi) =alpha-beta*(y+sqrt(y))+ &
                   a1*y**(1./3.)+a2*y**(2./3.)+ a3*y+a4*y**(4./3.)
           end do
        end do
     end do
  endif


  if (TypeBField == 'numeric') then

     call initialize_b_field(LZ, Phi, nPoint, nR, nT, bFieldMagnitude_III, &
          RadialDistance_III,Length_III, dLength_III,GradBCrossB_VIII,GradB_VIII,dBdt_III)

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
           IntegralH_III(1,iR,iPhi) = IntegralH_III(2,iR,iPhi)

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
  real                 :: dBdt_III(nPoint, nR, nT)
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
          RadialDistance_III,Length_III, dLength_III,GradBCrossB_VIII,GradB_VIII,dBdt_III)

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
  real                 :: dBdt_III(nPoint,nR,nT)

  !----------------------------------------------------------------------------------
  if (TypeBField == 'numeric') then
     
     call initialize_b_field(LZ, Phi, nPoint, nR, nT, bFieldMagnitude_III, &
          RadialDistance_III,Length_III, dLength_III,GradBCrossB_VIII,GradB_VIII,dBdt_III)
     
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
           NeutralHydrogen_III(1,iR,iPhi) = NeutralHydrogen_III(2,iR,iPhi)
        end do
     end do


  end if

end subroutine get_neutral_hydrogen
!============================================================
subroutine get_B_field(bFieldMagnitude_III)
  
  use ModHeidiSize,  ONLY: nPoint, nPa, nT, nR
  use ModHeidiBField
  use ModHeidiMain,  ONLY: Phi, LZ,Z
  implicit none

  real                 :: bFieldMagnitude_III(nPoint,nR,nT) ! Magnitude of magnetic field 
  real                 :: RadialDistance_III(nPoint,nR,nT)
  real                 :: dLength_III(nPoint-1,nR,nT)       ! Length interval between i and i+1  
  real                 :: GradBCrossB_VIII(3,nPoint,nR,nT)
  real                 :: GradB_VIII(3,nPoint,nR,nT)
  real                 :: Length_III(nPoint,nR,nT) 
  integer              :: iPhi, iR,iPitch
  real                 :: dBdt_III(nPoint,nR,nT)

  !----------------------------------------------------------------------------------
  
  call initialize_b_field(Z, Phi, nPoint, nR, nT, bFieldMagnitude_III, &
       RadialDistance_III,Length_III, dLength_III,GradBCrossB_VIII,GradB_VIII,dBdt_III)
  
  
end subroutine get_B_field

!============================================================
subroutine get_coef(dEdt_IIII,dMudt_III)
  
  use ModHeidiSize,   ONLY: nPoint,nPointEq, nPa, nT, nR,nE,DT
  use ModConst,       ONLY: cTiny  
  use ModHeidiMain,   ONLY: Phi, LZ, mu, EKEV, EBND,Z, funi, funt
  use ModHeidiDrifts, ONLY: VrConv, P1,P2
  use ModHeidiBField
  use ModHeidiBACoefficients

   use ModPlotFile, only: save_plot_file

  implicit none

  real, dimension(nPa)              :: bMirror_I, PitchAngle_I
  real, dimension(nPoint,nR,nT)     :: bFieldMagnitude_III, RadialDistance_III, Length_III, b4_III
  real, dimension(nR,nT,nE,nPA)     :: dEdt_IIII
  real, dimension(nR,nT,nE,nPA)     :: dMudt_III
  real, dimension(nPoint)           :: DriftR_I, DriftLambda_I,DriftPhi_I
  real, dimension(3,nPoint,nR,nT)   :: GradBCrossB_VIII, GradB_VIII
  real, dimension(3,nPoint)         :: VDrift_II
  real, dimension(nPoint-1,nR,nT)   :: dLength_III
  real, dimension(nPa,nR,nT)        :: h,s
  real, dimension(nR,nT,nE,nPA)     :: VPhi_IIII,VR_IIII
  real                              :: BouncedDriftR, BouncedDriftPhi, BouncedDriftLambda
  real                              :: HalfPathLength,Sb,SecondAdiabInv
  integer                           :: iMirror_I(2)
  real                              :: Energy
  real                              :: sinPitch,sin2Pitch,cosPitch,cos2Pitch
  real                              :: I2, CoeffE, CoeffMu
  real                              :: InvB, InvR, TermER, TermELambda,TermEPhi
  real                              :: TermER1, TermER2, TermEB, TermMuB
  real                              :: TermMuR, TermMuLambda,TermMuPhi
  real                              :: TermMuR1, TermMuR2
  real                              :: GradEqBR, GradEqBLambda, GradEqBPhi
  real, parameter                   :: Me = 7.9e15 ,Re = 6.371e6 
  integer                           :: iPhi, iR, iPitch,iE
  real                              :: dBdt_III(nPoint,nR,nT)
  real                              :: BouncedDBdt
  real, dimension(nPoint)           :: TR1, TR2, TLambda,TPhi
  real                              :: BouncedDriftR1,BouncedDriftR2
  integer                           :: iPoint

  character(LEN=500):: StringVarName, StringHeader, NameFile                              
  character(len=20) :: TypePosition                                                        
  character(len=20) :: TypeFile = 'ascii'  
  !----------------------------------------------------------------------------------

  call initialize_b_field(LZ, Phi, nPoint, nR, nT, bFieldMagnitude_III, &
       RadialDistance_III, Length_III, dLength_III, GradBCrossB_VIII,GradB_VIII,dBdt_III)

  do iE = 1, nE
     Energy = 1000. * EKEV(iE)
     
     do iPhi = 1, nT
        do iR = 1, nR
           do iPitch = 1, nPa
              
              PitchAngle_I(iPitch) = acos(mu(iPitch))
              sinPitch = sin(PitchAngle_I(iPitch))
              sin2Pitch = sinPitch * sinPitch

              cosPitch = cos(PitchAngle_I(iPitch))
              cos2Pitch = cosPitch * cosPitch
              
              call find_mirror_points (nPoint,  PitchAngle_I(iPitch), bFieldMagnitude_III(:,iR,iPhi), &
                   bMirror_I(iPitch),iMirror_I)


              call half_bounce_path_length(nPoint, iMirror_I(:),bMirror_I(iPitch),&
                   bFieldMagnitude_III(:,iR,iPhi), dLength_III(:,iR,iPhi), LZ(iR), HalfPathLength,Sb)
              

              call get_drift_velocity(nPoint, Energy,PitchAngle_I(iPitch), bFieldMagnitude_III(:,iR,iPhi),&
                   GradBCrossB_VIII(:,:,iR,iPhi), VDrift_II)
              

              DriftR_I(:)      = VDrift_II(1,:) ! Radial Drift
              DriftLambda_I(:) = VDrift_II(2,:) ! Theta Drift
              DriftPhi_I(:)    = VDrift_II(3,:) ! Azimuthal Drift
              
!              do iPoint = 1, nPoint
!                 TR1(iPoint) = 1./LZ(iR)*(DriftR_I(iPoint)+VrConv(iR,iPhi,iE,iPitch))
!                 TR2(iPoint) = (1./bFieldMagnitude_III(iPoint,iR,iPhi)) * (DriftR_I(iPoint)+&
!                      VrConv(iR,iPhi,iE,iPitch)) * GradB_VIII(1,iPoint,iR,iPhi)
!                 TLambda(iPoint) =(1./bFieldMagnitude_III(iPoint,iR,iPhi)) * DriftLambda_I(iPoint) *&
!                      GradB_VIII(2,iPoint,iR,iPhi) 
!                 TPhi(iPoint) =(1./bFieldMagnitude_III(iPoint,iR,iPhi)) * DriftPhi_I(iPoint) *&
!                      GradB_VIII(3,iPoint,iR,iPhi) 
!              end do

!              call get_bounced_drift(nPoint, LZ(iR),bFieldMagnitude_III(:,iR,iPhi),&
!                   bMirror_I(iPitch), iMirror_I(:),dLength_III(:,iR,iPhi),&
!                   TR1(:),BouncedDriftR1,RadialDistance_III(:,iR,iPhi)) 
!
!              call get_bounced_drift(nPoint, LZ(iR),bFieldMagnitude_III(:,iR,iPhi),&
!                   bMirror_I(iPitch), iMirror_I(:),dLength_III(:,iR,iPhi),&
!                   TR2(:),BouncedDriftR2,RadialDistance_III(:,iR,iPhi)) 
!            
!
!              call get_bounced_drift(nPoint, LZ(iR),bFieldMagnitude_III(:,iR,iPhi),&
!                   bMirror_I(iPitch), iMirror_I(:),dLength_III(:,iR,iPhi),&
!                   TLambda(:),BouncedDriftLambda,RadialDistance_III(:,iR,iPhi)) 
!              
!              call get_bounced_drift(nPoint, LZ(iR),bFieldMagnitude_III(:,iR,iPhi),&
!                   bMirror_I(iPitch), iMirror_I(:),dLength_III(:,iR,iPhi),&
!                   TPhi(:),BouncedDriftPhi,RadialDistance_III(:,iR,iPhi)) 



              call get_bounced_drift(nPoint, LZ(iR),bFieldMagnitude_III(:,iR,iPhi),&
                   bMirror_I(iPitch), iMirror_I(:),dLength_III(:,iR,iPhi),&
                   DriftR_I(:),BouncedDriftR,RadialDistance_III(:,iR,iPhi)) 

              call get_bounced_drift(nPoint, LZ(iR),bFieldMagnitude_III(:,iR,iPhi),&
                   bMirror_I(iPitch), iMirror_I(:),dLength_III(:,iR,iPhi),&
                   DriftLambda_I(:),BouncedDriftLambda,RadialDistance_III(:,iR,iPhi)) 
            

              call get_bounced_drift(nPoint, LZ(iR),bFieldMagnitude_III(:,iR,iPhi),&
                   bMirror_I(iPitch), iMirror_I(:),dLength_III(:,iR,iPhi),&
                   DriftPhi_I(:),BouncedDriftPhi,RadialDistance_III(:,iR,iPhi)) 
              
              call get_bounced_drift(nPoint, LZ(iR),bFieldMagnitude_III(:,iR,iPhi),&
                   bMirror_I(iPitch), iMirror_I(:),dLength_III(:,iR,iPhi),&
                   dBdt_III(:,iR,iPhi),BouncedDBdt,RadialDistance_III(:,iR,iPhi)) 
              
              if (Sb==0.0) Sb = cTiny
              I2 = funi(iPitch,iR,iPhi) * funi(iPitch,iR,iPhi)    

              GradEqBR      = GradB_VIII(1,nPointEq,iR,iPhi)
              GradEqBLambda = GradB_VIII(2,nPointEq,iR,iPhi)
              GradEqBPhi    = GradB_VIII(3,nPointEq,iR,iPhi)

              !\
              ! dE/dt bounce average coefficient
              !/

              CoeffE = EBND(iE)*(1.- I2/(I2+sin2Pitch))

              InvB = 1./ bFieldMagnitude_III(nPointEq,iR,iPhi)
              InvR =  1./LZ(iR)

              
              TermEB = CoeffE * InvB * BouncedDBdt/Sb
              TermER1 = CoeffE * InvB * GradEqBR * (BouncedDriftR/Sb + VrConv(iR,iPhi,iE,iPitch))
              TermER2 = -(2*EBND(iE)*I2)/(I2+sin2Pitch) * InvR * (BouncedDriftR/Sb+ VrConv(iR,iPhi,iE,iPitch)) 
              TermER  = TermER1 + TermER2
              TermELambda = CoeffE * InvB * GradEqBLambda * BouncedDriftLambda/Sb
              TermEPhi = CoeffE * InvB * GradEqBPhi* &
                   (BouncedDriftPhi/Sb + P1(iR,iPhi) + P2(iR,iPhi,iE,iPitch))


!              TermEB = CoeffE * InvB * BouncedDBdt/Sb
!              TermER1 = CoeffE * BouncedDriftR2/Sb
!              TermER2 = -(2*EBND(iE)*I2)/(I2+sin2Pitch) *BouncedDriftR1/Sb
!              TermER  = TermER1 + TermER2
!              TermELambda = CoeffE * BouncedDriftLambda/Sb
!              TermEPhi = CoeffE * BouncedDriftPhi/Sb

              dEdt_IIII(iR,iPhi,iE,iPitch)  = (TermEB + TermER + TermELambda + TermEPhi)!*6.371e6


              !\
              !  dMu/dT bounce averaged coefficient 
              !/


              CoeffMu = -(funi(iPitch,iR,iPhi)/(2.* funt(iPitch,iR,iPhi)))
!              CoeffMu = -(1.-cos2Pitch)*(SecondAdiabInv/(2.*HalfPathLength))/&
!                   (cosPitch*(I2 + 1. - cos2Pitch))
 

              TermMuB      = 0.5 * InvB*BouncedDBdt/Sb! * 1.e6
              TermMuR1     = 1./LZ(iR) * (BouncedDriftR/Sb + VrConv(iR,iPhi,iE,iPitch)) 
              TermMuR2     = 0.5 * InvB * GradEqBR * (BouncedDriftR/Sb+ VrConv(iR,iPhi,iE,iPitch))!*1e6 
              TermMuR      = TermMuR1 + TermMuR2
              TermMuLambda = -0.5 * InvB * GradEqBLambda * BouncedDriftLambda/Sb!*1e6
              TermMuPhi    = 0.5 * InvB * GradEqBPhi *&
                   (BouncedDriftPhi/Sb + P1(iR,iPhi) + P2(iR,iPhi,iE,iPitch))!*1.e6   ! in m



!              TermMuB      = 0.5 * InvB*BouncedDBdt/Sb! * 1.e6
!              TermMuR1     = BouncedDriftR1/Sb
!              TermMuR2     = 0.5 * BouncedDriftR2/Sb
!              TermMuR      = TermMuR1 + TermMuR2
!              TermMuLambda = -0.5 *BouncedDriftLambda/Sb
!              TermMuPhi    = 0.5 * BouncedDriftPhi/Sb



              dMudt_III(iR,iPhi,iE,iPitch) =  CoeffMu*(TermMuB + TermMuR + TermMuLambda + TermMuPhi)
                            

          end do
          dEdt_IIII(iR,iPhi,iE,1) = dEdt_IIII(iR,iPhi,iE,2) 
          dMudt_III(iR,iPhi,iE,1) = dMudt_III(iR,iPhi,iE,2)
       end do
    end do
     
 end do


 NameFile = 'EDOT.out'                                                                 
 StringHeader = 'EDOT values'
 StringVarName = 'R MLT Edot E PA'                                                               
 TypePosition = 'rewind'                                                                 
                                                                                          
  do iPitch =1, nPa
     do iE = 1, nE

        call save_plot_file(NameFile, &                                                   
             TypePositionIn = TypePosition,&                                              
             TypeFileIn     = TypeFile,&                                                  
             StringHeaderIn = StringHeader, &                                             
             nStepIn = 0, &                                                               
             TimeIn = 0.0, &                                                              
             ParamIn_I = (/ EKEV(iE), acos(mu(iPitch))*180./3.14159265, real(nR), real(NT)/), & 
             NameVarIn = StringVarName, &                                                 
             nDimIn = 2, &                                                                
             CoordMinIn_D = (/1.75, 0.0/),&                                               
             CoordMaxIn_D = (/6.5, 24.0/),&                                               
             VarIn_IIV = dEdt_IIII(:,:,iE,iPitch:iPitch))                      
        TypePosition = 'append'   
     end do
  end do


 NameFile = 'MUDOT.out'                                                                                     
 StringHeader = 'EDOT values'                                                                              
 StringVarName = 'R MLT Mudot E PA'                                                                         
 TypePosition = 'rewind'                                                                                   
                                                                                                           
  do iPitch =1, nPa                                                                                        
     do iE = 1, nE                                                                                         
        call save_plot_file(NameFile, &                                                                    
             TypePositionIn = TypePosition,&                                                               
             TypeFileIn     = TypeFile,&                                                                   
             StringHeaderIn = StringHeader, &                                                              
             nStepIn = 0, &                                                                                
             TimeIn = 0.0, &                                                                               
             ParamIn_I = (/ EKEV(iE), acos(mu(iPitch))*180./3.14159265, real(nR), real(NT)/), &            
             NameVarIn = StringVarName, &                                                                  
             nDimIn = 2, &                                                                                 
             CoordMinIn_D = (/1.75, 0.0/),&                                                                
             CoordMaxIn_D = (/6.5, 24.0/),&                                                                
             VarIn_IIV = dMudt_III(:,:,iE,iPitch:iPitch))                                                  
        TypePosition = 'append'                                                                            
     end do                                                                                                
  end do             


end subroutine get_coef

!============================================================

subroutine get_grad_curv_drift(VPhi_IIII,VR_IIII)
  
  use ModHeidiSize,   ONLY: nPoint,nPointEq, nPa, nT, nR,nE,DT
  use ModConst,       ONLY: cTiny  
  use ModHeidiMain,   ONLY: Phi, LZ, mu, EKEV, EBND,wmu,DPHI,Z
  use ModHeidiBField
  use ModHeidiBACoefficients
  use ModHeidiDrifts, ONLY: VrConv, P1,P2
  implicit none

  real, dimension(nPa)              :: bMirror_I, PitchAngle_I
  real, dimension(nPoint,nR,nT)     :: bFieldMagnitude_III, RadialDistance_III, Length_III, b4_III
  real, dimension(nPoint)           :: DriftR_I, DriftLambda_I,DriftPhi_I
  real, dimension(3,nPoint,nR,nT)   :: GradBCrossB_VIII, GradB_VIII
  real, dimension(3,nPoint)         :: VDrift_II
  real, dimension(nPoint-1,nR,nT)   :: dLength_III
  real, dimension(nPa,nR,nT)        :: h,s
  real, dimension(nR,nT,nE,nPA)     :: VPhi_IIII,VR_IIII
  real                              :: BouncedDriftR, BouncedDriftPhi, BouncedDriftLambda
  real                              :: HalfPathLength,Sb,SecondAdiabInv
  integer                           :: iMirror_I(2)
  real                              :: Energy
  real                              :: sinPitch,sin2Pitch,cosPitch,cos2Pitch
  real                              :: I2, CoeffE, CoeffMu
  real                              :: InvB, InvR, TermER, TermELambda,TermEPhi
  real                              :: TermER1, TermER2
  real                              :: TermMuR, TermMuLambda,TermMuPhi
  real                              :: TermMuR1, TermMuR2
  real                              :: GradEqBR, GradEqBLambda, GradEqBPhi
  real, parameter                   :: Me = 7.9e15 ,Re = 6.371e6 
  integer                           :: iPhi, iR, iPitch,iE
  real                              :: dBdt_III(nPoint,nR,nT)

  !----------------------------------------------------------------------------------

  call initialize_b_field(Z, Phi, nPoint, nR, nT, bFieldMagnitude_III, &
       RadialDistance_III, Length_III, dLength_III, GradBCrossB_VIII,GradB_VIII,dBdt_III)

  do iE = 1, nE
     Energy = 1000. * EKEV(iE)
     do iPhi = 1, nT
        do iR = 1, nR
           do iPitch = 1, nPa
              
              PitchAngle_I(iPitch) = acos(mu(iPitch))
              sinPitch = sin(PitchAngle_I(iPitch))
              sin2Pitch = sinPitch * sinPitch

              cosPitch = cos(PitchAngle_I(iPitch))
              cos2Pitch = cosPitch * cosPitch
              
              call find_mirror_points (nPoint,  PitchAngle_I(iPitch), bFieldMagnitude_III(:,iR,iPhi), &
                   bMirror_I(iPitch),iMirror_I)

              call half_bounce_path_length(nPoint, iMirror_I(:),bMirror_I(iPitch),&
                   bFieldMagnitude_III(:,iR,iPhi), dLength_III(:,iR,iPhi), Z(iR), HalfPathLength,Sb)
              

              call get_drift_velocity(nPoint, Energy,PitchAngle_I(iPitch), bFieldMagnitude_III(:,iR,iPhi),&
                   GradBCrossB_VIII(:,:,iR,iPhi), VDrift_II)
              

              DriftR_I(:)      = VDrift_II(1,:) ! Radial Drift
              DriftLambda_I(:) = VDrift_II(2,:) ! Theta Drift
              DriftPhi_I(:)    = VDrift_II(3,:) ! Azimuthal Drift


              call get_bounced_drift(nPoint, Z(iR),bFieldMagnitude_III(:,iR,iPhi),&
                   bMirror_I(iPitch), iMirror_I(:),dLength_III(:,iR,iPhi),&
                   DriftR_I(:),BouncedDriftR,RadialDistance_III(:,iR,iPhi)) 

              call get_bounced_drift(nPoint, Z(iR),bFieldMagnitude_III(:,iR,iPhi),&
                   bMirror_I(iPitch), iMirror_I(:),dLength_III(:,iR,iPhi),&
                   DriftLambda_I(:),BouncedDriftLambda,RadialDistance_III(:,iR,iPhi)) 
            

              call get_bounced_drift(nPoint, Z(iR),bFieldMagnitude_III(:,iR,iPhi),&
                   bMirror_I(iPitch), iMirror_I(:),dLength_III(:,iR,iPhi),&
                   DriftPhi_I(:),BouncedDriftPhi,RadialDistance_III(:,iR,iPhi)) 
              
              if (Sb==0.0) Sb = cTiny

              VPhi_IIII(iR,iPhi,iE,iPitch) = BouncedDriftPhi/Sb
!              VR_IIII(iR,iPhi,iE,iPitch) = sqrt( (BouncedDriftR/Sb)**2 + (BouncedDriftLambda/Sb)**2)
              VR_IIII(iR,iPhi,iE,iPitch) = BouncedDriftR/Sb

          end do
          VPhi_IIII(iR,iPhi,iE,1) = VPhi_IIII(iR,iPhi,iE,2)
          VR_IIII(iR,iPhi,iE,1)   =  VR_IIII(iR,iPhi,iE,2) 
      end do
    end do
 end do

!stop

 

end subroutine get_grad_curv_drift
