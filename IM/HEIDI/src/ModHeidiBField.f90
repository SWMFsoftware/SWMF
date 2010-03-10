module ModHeidiBField

  implicit none

contains

  subroutine initialize_b_field (L_I, Phi_I, nPoint, nR, nPhi, bFieldMagnitude_III, &
       RadialDistance_III, Length_III, dLength_III,GradBCrossB_VIII,GradB_VIII,dBdt_III)

    use ModHeidiInput, ONLY: TypeBfieldGrid
    use ModNumConst,   ONLY: cTiny, cPi
    use ModConst,      ONLY: cMu
    use ModHeidiIO,    ONLY: time
    use ModHeidiMain,  ONLY: LZ, BHeidi_III, SHeidi_III, RHeidi_III!,t

    integer, intent(in)    :: nPoint                              ! Number of points along the field line
    integer, intent(in)    :: nR                                  ! Number of points in the readial direction
    integer, intent(in)    :: nPhi                                ! Number of points in the azimuthal direction
    real,    intent(in)    :: L_I(nR)                             ! L shell value
    real,    intent(in)    :: Phi_I(nPhi)                         ! Phi values
    real,    intent(inout) :: bFieldMagnitude_III(nPoint,nR,nPhi) ! Magnitude of magnetic field 
    real,    intent(out)   :: Length_III(nPoint,nR,nPhi)          ! Length of the field line
    real,    intent(out)   :: RadialDistance_III(nPoint,nR,nPhi)
    real,    intent(out)   :: dLength_III(nPoint-1,nR,nPhi)       ! Length interval between i and i+1
    !Local Variables    
    real                   :: LatMin                              ! Minimum Latitude 
    real                   :: LatMax                              ! Maximum Latitude 
    real                   :: Lat                                 ! Latitude
    real                   :: dLat                                ! Latitude cell size
    real                   :: SinLat2, CosLat2, SinLat4 
    real                   :: LatMin2, LatMax2, LatMinMax
    real                   :: BetaF,dLatNew,beta1,beta2
    real                   :: alpha2, alpha4
    real                   :: bField_VIII(3,nPoint,nR,nPhi)
    real                   :: GradBCrossB_VIII(3,nPoint,nR,nPhi) 
    real                   :: GradB_VIII(3,nPoint,nR,nPhi)
    double precision       :: r,x,y,Vr,Vtheta,Vphi,mag,Br,Btheta,Bphi,GradR, GradTheta, GradPhi,dBdt
    real                   :: DipoleFactor,delta, dSdTheta
    real                   :: a,b,c,gamma,cos2Lat,cosLat,L
    integer                :: iR, iPhi,iPoint     

    real                   :: cos2LatMin,cos2LatMax,cos2Lat1,cos2Lat2,cos2Lat3
    complex                :: root(3)
    integer                :: nroot, i
    real                   :: dd
    real                   :: VrConv(nR,nPhi)

    !Parameters
    real, parameter        :: DipoleStrength =  0.32   ! nTm^-3
    real, parameter        :: alpha0 =1.0              ! alpha is the stretching factor in z direction
    real, parameter        :: Me = 7.19e15

    real :: dBdt_III(nPoint,nR,nPhi), dBdtTemp,p,w,c1
     real, dimension(nPoint,nR,nPhi) :: GradB0R_III, GradB0Theta_III, GradB0Phi_III
    real :: beta,alpha,t
    !----------------------------------------------------------------------------------

    DipoleFactor= Me

    !\
    ! Dipole magnetic field with uniform number of points along the field line
    !/

    select case(TypeBFieldGrid)

    case('mhd')
       bFieldMagnitude_III = BHeidi_III
       RadialDistance_III  = RHeidi_III 
       Length_III          = SHeidi_III 
!       write(*,*) bFieldMagnitude_III(50,:,1)

       call get_gradB0(nPoint,nR,nPhi,L_I,Phi_I,GradB0R_III, GradB0Theta_III, GradB0Phi_III)


    case('uniform')
       do iPhi =1, nPhi
          do iR =1, nR 
             LatMax =  acos(sqrt(1./L_I(iR)))
             LatMin = -LatMax
             dLat   = (LatMax-LatMin)/(nPoint-1)
             Lat = LatMin

             do iPoint = 1, nPoint
                SinLat2 = sin(Lat)**2
                CosLat2 = 1.0 - SinLat2  
                bFieldMagnitude_III(iPoint,iR,iPhi) = DipoleStrength*sqrt(1.0+3.0*SinLat2)/(L_I(iR)*CosLat2)**3
                RadialDistance_III(iPoint,iR,iPhi) = L_I(iR)*CosLat2
                Length_III(iPoint,iR,iPhi) = dipole_length(L_I(ir),LatMin,Lat) 
                Lat = Lat + dLat
             end do
          end do
       end do

       !\
       ! Dipole magnetic field with non-uniform number of points along the field line. 
       ! More refined at the equator, coarser towards the poles
       !/

    case('nonuniform')
       do iPhi =1, nPhi
          do iR =1, nR 
             LatMax =  acos(sqrt(1./L_I(iR)))
             LatMin = -LatMax
             dLat   = (LatMax-LatMin)/(nPoint-1)
             Lat = LatMin

             LatMin2 = LatMin**2
             LatMax2 = LatMax**2
             LatMinMax = LatMin*LatMax

             beta1 = 6.0*(-nPoint+cTiny*nPoint+1.0)*(nPoint-1.0) 
             beta2 = nPoint*(2.0*LatMax2*nPoint + 2.0*LatMinMax*nPoint + 2.0*LatMin2*nPoint &
                  - 4.0*LatMinMax - LatMin2 - LatMax2)
             BetaF = - beta1/beta2  

             do iPoint = 1, nPoint
                Lat = LatMin+(iPoint-1)*dLat
                SinLat2 = sin(Lat)**2
                CosLat2 = 1.0 - SinLat2
                bFieldMagnitude_III(iPoint,iR,iPhi) = DipoleStrength*sqrt(1.0+3.0*SinLat2)/(L_I(iR)*CosLat2)**3
                RadialDistance_III(iPoint,iR,iPhi) = L_I(iR)*CosLat2
                Length_III(iPoint,iR,iPhi) = dipole_length(L_I(iR),LatMin,Lat) 
                dLatNew =(BetaF*(LatMin+(iPoint-1)*dLat)**2+cTiny)*dLat
                Lat = Lat+dLatNew
             end do
          end do
       end do

       !\
       ! Stretched dipole magnetic field, with azimuthal symmetry. 
       !/


    case('stretched')
       t = 0.0!t/3600.
       w = 2*cPi/50.
       ! Time dependent B field
       !       alpha = alpha0 + 1.1*sin(w*t/3600.)

       alpha = alpha0

       beta = 1./alpha

       alpha2 = alpha * alpha
       alpha4 = alpha2 * alpha2
       beta2 = beta* beta

       dd = 0.0
       do iPhi =1, nPhi
          y = cos(Phi_I(iPhi))

          !\
          ! Calculate the maximum latitude for this case (from the equation of the field line)
          ! Yields a cubic equation in latitude
          !/ 

          gamma = y**2 + beta2 * (1-y**2)
          a = gamma - alpha2
          b = alpha2

          do iR =1, nR 
             c = -1./(LZ(iR)*LZ(iR))
             call get_cubic_root(a,b,dd,c,root,nroot)
             do i =1, nroot
                if ((aimag(root(i)) <= 1.e-5).and. (real(root(i))<=1.0) .and. (real(root(i))>=0.0)) &
                     cos2Lat1 = real(root(i))
             end do

             LatMax = acos(sqrt(cos2Lat1))
             LatMin = -LatMax
             dLat   = (LatMax-LatMin)/(nPoint-1)
             Lat = LatMin


             do iPoint = 1, nPoint
                x = sin(Lat)





                !\
                ! Radial distance for the stretched dipole is calclulated as a function of the 
                ! dipole radial distance; Easy calculation for this case.
                !/

                RadialDistance_III(iPoint,iR,iPhi) = L_I(iR)*(1. - x**2)*&
                     sqrt((1. - x ** 2) * (y**2 + beta**2 *(1-y**2)) + alpha ** 2 * x ** 2)

                r =  RadialDistance_III(iPoint,iR,iPhi)

                ! \
                !  Magnetic field components for the uniformly stretched dipole in y and z.
                !/


                Br = -3 * r ** 2 * x * (1 - x ** 2) * y ** 2 * alpha * (r ** 2 * (1 - x ** 2) * &
                     y ** 2 + r ** 2 * (1 - x ** 2) * (1 - y ** 2) * beta ** 2 + r ** 2 * x ** 2 *&
                     alpha ** 2) ** (-0.5D1 / 0.2D1) - 3 * r ** 2 * x * (1 - x ** 2) * (1 - y ** 2) *&
                     alpha * (r ** 2 * (1 - x ** 2) * y ** 2 + r ** 2 * (1 - x ** 2) * (1 - y ** 2) *&
                     beta ** 2 + r ** 2 * x ** 2 * alpha ** 2) ** (-0.5D1 / 0.2D1) + 1 / alpha * &
                     (-2 * r ** 2 * x ** 2 * alpha ** 2 * (r ** 2 * (1 - x ** 2) * y **2 + &
                     r ** 2 * (1 - x ** 2) * (1 - y ** 2) * beta ** 2 + r ** 2 * x** 2 * &
                     alpha ** 2) ** (-0.5D1 / 0.2D1) - (-r ** 2 * (1 - x ** 2) * y ** 2 - &
                     r ** 2 * (1 - x ** 2) * (1 - y ** 2) * beta ** 2) * (r ** 2 * (1 - x ** 2) *&
                     y ** 2 + r ** 2 * (1 - x ** 2) * (1 - y ** 2) * beta ** 2 + r ** 2 * x ** 2 * &
                     alpha ** 2) ** (-0.5D1 / 0.2D1)) * x


                Btheta = -0.3D1 * dble(r ** 2) * dble(x ** 2) * sqrt(dble(1 - x ** 2)) * &
                     dble(y ** 2) * dble(alpha) * dble((r ** 2 * (1 - x ** 2) * y ** 2 + &
                     r ** 2 * (1 - x ** 2) * (1 - y ** 2) * beta ** 2 + r ** 2 * x ** 2 * &
                     alpha ** 2) ** (-0.5D1 / 0.2D1)) - 0.3D1 * dble(r ** 2) *dble(x ** 2) * &
                     sqrt(dble(1 - x ** 2)) * dble(1 - y ** 2) * dble(alpha) * dble((r ** 2 * &
                     (1 - x ** 2) * y ** 2 + r ** 2 * (1 - x ** 2) * (1 - y ** 2) * &
                     beta ** 2 + r ** 2 * x ** 2 * alpha ** 2) ** (-0.5D1 / 0.2D1)) - &
                     0.1D1 / dble(alpha) * dble(-2 * r ** 2 * x ** 2 * alpha ** 2 *&
                     (r ** 2 * (1 - x ** 2) * y ** 2 + r ** 2 * (1 - x ** 2) * (1 - y ** 2) * &
                     beta ** 2 + r ** 2 * x ** 2 * alpha ** 2) ** (-0.5D1 / 0.2D1) -&
                     (-r ** 2 * (1 - x ** 2) * y ** 2 - r ** 2 * (1 - x ** 2) * (1 - y ** 2) *&
                     beta ** 2) * (r ** 2 * (1 - x ** 2) * y ** 2 + r ** 2 * (1 - x ** 2) *&
                     (1 - y ** 2) * beta ** 2 + r ** 2 * x ** 2 * alpha ** 2) ** (-0.5D1 / 0.2D1)) *&
                     sqrt(dble(1 - x ** 2))


                Bphi = 0.0


                !\
                ! Gradient drift components (Vr,Vteta,Vphi)
                !/



                Vr = -0.3D1 / dble(r ** 9) * dble(-1 + x ** 2) * dble(y) * &
                     sqrt(dble(1 - y ** 2)) * dble(-1 + beta ** 2) * dble(-x ** 4 * y ** 4 - 15 &
                     * alpha ** 4 * x ** 2 + 7 * x ** 4 * alpha ** 4 + 2 * y ** 4 * x** 2 + 2 * &
                     beta ** 4 * x ** 2 - beta ** 4 * x ** 4 + 4 * beta ** 2 * x ** 2 * &
                     y ** 2 - y ** 4 * beta ** 4 - 4 * y ** 4 * beta ** 2 * x ** 2 - 6 * &
                     x ** 4 * alpha ** 2 * y ** 2 - 6 * x ** 4 * alpha ** 2 * beta ** 2 - &
                     2 * x ** 4 * y ** 2 * beta ** 2 + 2 * x ** 4 * y** 4 * beta ** 2 - &
                     4 * beta ** 4 * x ** 2 * y ** 2 + 2 * y ** 4 *beta ** 4 * x ** 2 + &
                     2 * beta ** 4 * x ** 4 * y ** 2 - beta ** 4 * x ** 4 * y ** 4 + 2 * &
                     y ** 4 * beta ** 2 + 2 * beta ** 4 * y ** 2 - 2 * y ** 2 * beta ** 2 - &
                     beta ** 4 - y ** 4 - 6 * alpha ** 2 *beta ** 2 * x ** 2 * y ** 2 + &
                     6 * x ** 4 * alpha ** 2 * beta ** 2 * y ** 2 + 6 * alpha ** 2 * &
                     beta ** 2 * x ** 2 + 6 * y ** 2 * alpha ** 2 * x ** 2) / dble((y ** 2 - &
                     x ** 2 * y ** 2 + beta ** 2 - y** 2 * beta ** 2 - beta ** 2 * x ** 2 &
                     + beta ** 2 * x ** 2 * y ** 2 + x ** 2 * alpha ** 2) ** 7) * &
                     dble((r ** 2 * (y ** 2 - x ** 2 * y ** 2 + beta ** 2 - y ** 2 * &
                     beta ** 2 - beta ** 2 * x ** 2 + beta ** 2 * x ** 2 * y ** 2 + &
                     x ** 2 * alpha ** 2)) ** (-0.1D1 / 0.2D1)) / dble(alpha ** 3)




                Vtheta = dble(30 * alpha ** 2 * y ** 4 * beta ** 2 * x ** 6 - 9 * beta ** 4 * alpha ** 2 - &
                     15 * alpha ** 2 * y ** 4 * beta ** 4 * x ** 6 - 99 * x ** 2 * alpha ** 4 * y ** 2 &
                     * beta ** 2 - 30 * alpha ** 2 * y ** 2 * beta ** 2 * x ** 6 + 9 * y ** 6 * beta ** 4 - &
                     6 * alpha ** 2 * y ** 4 * beta ** 2 * x ** 2 + 21 * alpha ** 2 * y ** 4 *beta ** 4 &
                     * x ** 4 + 27 * y ** 6 * beta ** 2 * x ** 2 + 3 * beta ** 6 - 6 * alpha ** 2 * y ** 2 &
                     * beta ** 4 * x ** 2 - 3 * y ** 6 *beta ** 6 - 9 * y ** 6 * beta ** 4 * x ** 6 - &
                     27 * y ** 2 * beta ** 4 * x ** 2 + 39 * alpha ** 4 * beta ** 2 * x ** 6 + &
                     18 * y ** 2* beta ** 4 * alpha ** 2 - 138 * alpha ** 4 * beta ** 2 * x ** 4 + &
                     27 * y ** 2 * beta ** 4 * x ** 4 + 3 * y ** 6 * beta ** 6 * x **6 - 9 * y ** 4 * &
                     beta ** 4 * alpha ** 2 - 9 * y ** 4 * beta ** 2 * x ** 6 + 27 * y ** 4 * beta ** 2 * x ** 4 + &
                     9 * y ** 6 * beta ** 2 * x ** 6 - 42 * alpha ** 2 * y ** 2 * beta ** 4 * x ** 4 + &
                     3 * alpha ** 2 * y ** 4 * beta ** 4 * x ** 2 + 9 * y ** 4 * beta ** 6 - 9 * y ** 2 * &
                     beta ** 6 - 18 * y ** 4 * beta ** 4 - 9 * y ** 6 * beta ** 2 + 108 * x ** 4 * alpha ** 6 - &
                     135 * x ** 2 * alpha ** 6 - 9 * x ** 2 * y ** 6 + 9 * x ** 4 * y ** 6 - 3 * x ** 6 * y ** 6 - &
                     3 * x ** 6 * beta ** 6 + 9 * x ** 4 * beta ** 6 + 9 * y ** 4 * beta ** 2 + 9 * y ** 2 * &
                     beta ** 4 - 9 * x ** 2 * beta ** 6 - 9 * y ** 4 * alpha ** 2 - 21 * alpha ** 6 * x ** 6 - &
                     39 * alpha ** 4 * y ** 2 * beta ** 2 * x ** 6 + 138 * alpha ** 4 * y ** 2 * &
                     beta ** 2 * x ** 4 - 18 * y ** 2 * beta ** 2 * alpha ** 2 + 54 * y ** 4 * beta ** 4 * &
                     x ** 2 + 39 * alpha ** 4 * y ** 2 * x ** 6 - 9 * y ** 4 * beta ** 6 * x ** 6 + &
                     18 * y ** 4 * beta ** 4 * x ** 6 + 27 * y** 4 * beta ** 6 * x ** 4 - 27 * y ** 6 * &
                     beta ** 4 * x ** 2 - 9 * y ** 2 * beta ** 4 * x ** 6 + 3 * alpha ** 2 * y ** 4 * x ** 2 - &
                     138 * alpha ** 4 * y ** 2 * x ** 4 + 27 * y ** 6 * beta ** 4 * x ** 4 + 27 * y ** 2 * &
                     beta ** 6 * x ** 2 + 99 * y ** 2 * x ** 2 * alpha ** 4 - 27 * y ** 4 * beta ** 6 * x ** 2 + &
                     3 * y ** 6 + 42 * beta ** 2 * y ** 2 * x ** 4 * alpha ** 2 - 15 * alpha ** 2 * &
                     beta ** 4 * x ** 6 - 27 * y ** 4 * beta ** 2 * x ** 2 + 99 * x ** 2 * beta ** 2 * &
                     alpha ** 4 + 9 * y ** 2 * beta ** 6 * x ** 6 - 27 * y **2 * beta ** 6 * x ** 4 - &
                     9 * y ** 6 * beta ** 6 * x ** 4 - 27 * y** 6 * beta ** 2 * x ** 4 - 54 * y ** 4 *&
                     beta ** 4 * x ** 4 + 3 * alpha ** 2 * beta ** 4 * x ** 2 + 18 * y ** 4 * &
                     beta ** 2 * alpha ** 2 + 21 * alpha ** 2 * y ** 4 * x ** 4 + 6 * alpha ** 2 * &
                     beta ** 2 * x ** 2 * y ** 2 + 30 * alpha ** 2 * y ** 2 * beta ** 4 * x ** 6 - &
                     42 * alpha ** 2 * y ** 4 * beta ** 2 * x ** 4 + 9 * y ** 6 * beta ** 6 * x ** 2 + &
                     21 * alpha ** 2 * beta ** 4 * x ** 4 - 15 * alpha ** 2 * y ** 4 * x ** 6) * &
                     dble(x) * dble(-1 + beta ** 2) *sqrt(dble(1 - y ** 2)) * dble(y) *&
                     sqrt(dble(1 - x ** 2)) / dble(r ** 9) / dble((y ** 2 - y ** 2 * x ** 2 + beta ** 2 - &
                     beta ** 2 * y ** 2 - beta ** 2 * x ** 2 + beta ** 2 * y ** 2 * x ** 2 + x ** 2 * &
                     alpha ** 2) ** 8) * dble((r ** 2 * (y ** 2 - y ** 2 * x ** 2 + beta ** 2 -&
                     beta ** 2 * y ** 2 - beta ** 2 * x ** 2 + beta ** 2 *y ** 2 * x ** 2 + &
                     x ** 2 * alpha ** 2)) ** (-0.1D1 / 0.2D1)) / dble(alpha ** 3)


                Vphi = dble(144*alpha**2*y**6*beta**2*x**4+72*alpha**2*y**4*beta**2*x**6-3*alpha**4*beta**4*x**6&
                     *y**4-99*x**4*alpha**6*beta**2*y**2+126*x**2*alpha**4*y**4*beta**2+ &
                     72*alpha**2*y**6*beta**4*x**6-144*alpha**2*y**4*beta**4*x**6-126*x**2*alpha**4*y**2*beta**2-9*x**4*y**8+ &
                     36 * y ** 6 * beta ** 4 - 18 * y ** 8 * beta ** 4 + 3 * x ** 6 * y ** 8 + 3 * beta ** 8 * &
                     x ** 6 - 18 * beta** 8 * y ** 4 + 72 * alpha ** 2 * y ** 4 * beta ** 2 * x ** 2 + &
                     288 * alpha ** 2 * y ** 4 * beta ** 4 * x ** 4 + 36 * y ** 6 * beta ** 2 * x ** 2 - &
                     3 * beta ** 8 + 48 * alpha ** 2 * beta ** 6 * x ** 4 * y ** 6 - 72 * alpha ** 2 * &
                     beta ** 6 * x ** 6 * y ** 2 + 72 * alpha ** 2 * y ** 2 * beta ** 4 * x ** 2 + &
                     12 * beta ** 8 * y **2 - 36 * y ** 6 * beta ** 6 - 27 * x ** 2 * alpha ** 6 * &
                     beta ** 2 * y ** 2 - 48 * alpha ** 2 * beta ** 6 * x ** 4 - 3 * alpha ** 4* y ** 4 * &
                     x ** 6 - 78 * alpha ** 6 * y ** 2 * x ** 6 + 18 * y **8 * beta ** 4 * x ** 6 - &
                     36 * y ** 6 * beta ** 4 * x ** 6 - 63 * x ** 2 * alpha ** 4 * beta ** 4 - 78 * alpha ** 6 * &
                     beta ** 2 * x ** 6 + 36 * beta ** 8 * x ** 4 * y ** 2 - 36 * beta ** 8 * x ** 2 * y ** 2 - &
                     12 * beta ** 8 * x ** 6 * y ** 2 + 36 * y ** 6 * beta ** 6 * x ** 6 + 36 * beta ** 8 * &
                     x ** 4 * y ** 6 + 27 * x ** 2 * alpha ** 6 * beta ** 2 + 99 * x ** 4 * alpha ** 6 * &
                     beta ** 2 + 12 *y ** 6 * beta ** 2 * x ** 6 - 144 * alpha ** 2 * y ** 2 * beta **4 * x ** 4 - &
                     144 * alpha ** 2 * y ** 6 * beta ** 4 * x ** 4 + 72 * alpha ** 2 * y ** 6 * beta ** 4 * &
                     x ** 2 - 144 * alpha ** 2 * y ** 4 * beta ** 4 * x ** 2 - 24 * alpha ** 2 * beta ** 6 * &
                     x ** 2 * y ** 6 - 72 * alpha ** 2 * beta ** 6 * x ** 2 * y ** 2 + 9 * x **2 * y ** 8 - &
                     3 * beta ** 8 * y ** 8 + 36 * y ** 4 * beta ** 6 + 144 * alpha ** 2 * beta ** 6 * x ** 4 *&
                     y ** 2 + 72 * alpha ** 2 * beta ** 6 * x ** 2 * y ** 4 - 144 * alpha ** 2 * beta ** 6 *&
                     x**4*y**4-63*x**2*alpha**4*beta**4*y**4+12*y**8*beta**2-12*y**2*beta**6-18*y**4*beta**4-12*y**6*beta**2+ &
                     54*x**6*alpha**8-108*x**4*alpha**8-6*alpha**4*y**2*beta**2*x**6+132*alpha**4*y**2*beta**2*x**4-63*&
                     x**2*alpha**4*y**4+ &
                     3*beta**8*y**8*x**6+54*y**4*beta**4*x**2+27*x**2*alpha**6*y**2+36*y**8*beta**2*x**4+18*beta**8*x**6*y**4- &
                     12*y**8*beta**2*x**6+24*alpha**2*beta**6*x**2+66*alpha**4*y**4*x**4-54*y**8*beta**4*x**4-36*y**4*beta**6*x**6+&
                     18*y**4*beta**4*x**6+108*y**4*beta**6*x**4-108*y**6*beta**4*x**2+24*alpha**2*beta**6*x**6+108*y**6*beta**4* &
                     x**4+36*y**2*beta**6*x**2-108*y**4*beta**6*x**2-3*y**8+9*beta**8*x**2+6*alpha**4*&
                     y**4*beta**2*x**6-132*alpha**4* &
                     y**4*beta**2*x**4-9*beta**8*y**8*x**4+66*alpha**4*beta**4*x**4+54*beta**8*x**2*y**4-48*alpha**2*y**6*x**4- &
                     3*alpha**4*beta**4*x**6+9*beta**8*y**8*x**2+54*y**8*beta**4*x**2+24*alpha**2*y**6*x**6+36*y**8*beta**6*x**4+&
                     12*y**2*beta**6*x**6-36*y**8*beta**6*x**2+99*x**4*alpha**6*y**2-36*y**2*beta**6*x**4-54*beta**8 *&
                     x**4*y**4-108*y**6*beta**6*x**4-36*y**6*beta**2*x**4-54*y**4*beta**4*x**4-12*y**8* beta**6*x**6-36*y**8* &
                     beta**2*x**2-12*beta**8*x**6*y**6+72*alpha**2*beta**6*x**6*y**4+126*x**2*alpha**4*beta**4*y**2+6*alpha**4*&
                     beta**4*x**6*y**2+12*y**8*beta**6+12*beta**8*y**6-9*beta**8*x**4+72*alpha**2*y**2*beta**4*x**6+ &
                     78*alpha**6*beta**2*x**6*y**2-144*alpha**2*y**4*beta**2*x**4+108*y**6*beta**6*x**2+24*alpha**2*y**6*x**2- &
                     36*beta**8*x**2*y**6-72*alpha**2*y**6*beta**2*x**2-24*alpha**2*beta**6*x**6*y**6-132*alpha**4*beta**4*x**4* &
                     y**2-72*alpha**2*y**6*beta**2*x**6+66*alpha**4*beta**4*x**4*y**4)*sqrt(dble(1-x**2))/dble(r**9)/&
                     dble((y**2-y**2*x**2+beta**2-beta**2*y**2-beta**2*x**2+beta**2*y**2*x**2+x**2*alpha**2)**8)*dble((r**2*(y**2- &
                     y**2*x**2+beta**2- beta**2*y**2-beta**2*x**2+beta**2*y**2*x**2+x**2*alpha**2))**(-0.1D1/0.2D1))/dble(alpha**3)

                GradR = 3 * (4 * beta ** 2 * x ** 2 * y ** 2 - beta ** 4 * y ** 4 + 2* beta ** 4 * y ** 2 + &
                     2 * beta ** 4 * x ** 2 + 5 * x ** 4 * alpha ** 4 - 4 * x ** 2 * alpha ** 2 * beta ** 2 *&
                     y ** 2 + 4 * beta ** 2 * x ** 4 * y ** 2 * alpha ** 2 + 2 * beta ** 4 * x ** 2 * &
                     y **4 - 4 * x ** 4 * alpha ** 2 * beta ** 2 + 2 * beta ** 4 * x ** 4 * y ** 2 - &
                     beta ** 4 * x ** 4 * y ** 4 - 4 * beta ** 4 * x ** 2 * y ** 2 + 4 * beta ** 2 * x ** 2 * &
                     alpha ** 2 - beta ** 4 * x ** 4 + 2 * y ** 4 * x ** 2 + 2 * y ** 4 * beta ** 2 - x ** 4 * &
                     y ** 4 - y ** 4 - 2 * beta ** 2 * y ** 2 - beta ** 4 - 9 * alpha ** 4 * x ** 2 + 2 * &
                     x ** 4 * y ** 4 * beta ** 2 - 4 * x ** 4 * y ** 2 * alpha ** 2 - 2 * x ** 4 * beta ** 2 * &
                     y ** 2 - 4 * y ** 4 * beta ** 2 * x ** 2 + 4 * y ** 2 * x ** 2 * alpha ** 2) / &
                     r ** 7 * (-(4 * beta ** 2 * x ** 2 * y ** 2 - beta ** 4 * y ** 4 + 2 * beta ** 4 * &
                     y** 2 + 2 * beta ** 4 * x ** 2 + 5 * x ** 4 * alpha ** 4 - 4 * x ** 2 * alpha ** 2 * beta ** 2 * &
                     y ** 2 + 4 * beta ** 2 * x ** 4 * y ** 2 * alpha ** 2 + 2 * beta ** 4 * x ** 2 * y ** 4 - &
                     4 * x ** 4 * alpha ** 2 * beta ** 2 + 2 * beta ** 4 * x ** 4 * y ** 2 - beta ** 4 * x ** 4 * &
                     y ** 4 - 4 * beta ** 4 * x ** 2 * y ** 2 + 4 * beta ** 2 * x ** 2 * alpha ** 2 - beta ** 4 * &
                     x ** 4 + 2 * y ** 4 * x** 2 + 2 * y ** 4 * beta ** 2 - x ** 4 * y ** 4 - y ** 4 - &
                     2 * beta ** 2 * y ** 2 - beta ** 4 - 9 * alpha ** 4 * x ** 2 + 2 * x ** 4 * y ** 4 * beta ** 2 - &
                     4 * x ** 4 * y ** 2 * alpha ** 2 - 2 * x ** 4 * beta ** 2 * y ** 2 - 4 * y ** 4 * beta ** 2 * &
                     x ** 2 + 4 * y ** 2 * x ** 2 * alpha ** 2) / r ** 6 / alpha ** 2 / (y ** 2 - x ** 2 * y ** 2 + &
                     beta ** 2 - beta ** 2 * y ** 2 - beta ** 2 * x ** 2 + beta ** 2 * x ** 2 * y ** 2 + x ** 2 * &
                     alpha ** 2) ** 5) ** (-0.1D1 / 0.2D1) / (y ** 2 - x ** 2 * y ** 2 + beta ** 2 - beta ** 2 * &
                     y ** 2 - beta ** 2 * x ** 2 + beta ** 2 * x ** 2 * y ** 2 + x ** 2 * alpha ** 2) ** 5 / alpha ** 2

                GradTheta = -0.3D1 * dble(-3 * y ** 4 * alpha ** 2 - 3 * beta ** 4 * alpha ** 2 + 6 * &
                     y ** 4 * alpha ** 2 * beta ** 2 + 14 * x ** 2 * alpha ** 4 * y ** 2 - 6 * y ** 2 * &
                     beta ** 2 * alpha ** 2 - 14 * x ** 2 * alpha ** 4 * beta ** 2 * y ** 2 + beta ** 6 - &
                     3 * beta ** 4 * y ** 4 * alpha ** 2 - 9 * beta ** 2 * x ** 4 * alpha ** 4 - 6 * beta ** 6 * &
                     x ** 2 * y ** 4 + 3 * beta ** 4 * x ** 4 * alpha ** 2 - 3 * beta ** 6 * x ** 4 * y ** 2 + 3 * &
                     beta ** 6 * x ** 4 * y ** 4 + 6 * beta ** 6 * x ** 2 * y ** 2 - 6 * beta ** 4 * y ** 4 + 3 * &
                     beta ** 4 * y ** 2 + beta ** 6 * x ** 4 - beta ** 6 * y ** 6 + 6 * y ** 6 * beta ** 2 * x ** 2 + &
                     2 * y ** 6 * beta ** 6 * x ** 2 + 3 * y ** 4 * alpha ** 2 * x ** 4 + 6 * beta ** 2 * x ** 4 * &
                     y ** 2 *alpha ** 2 + 12 * beta ** 4 * x ** 2 * y ** 4 + 3 * beta ** 4 * x** 4 * y ** 2 - 6 * &
                     beta ** 4 * x ** 4 * y ** 4 - 6 * beta ** 4 *x ** 2 * y ** 2 + 3 * y ** 4 * beta ** 2 + y ** 6 -&
                     y ** 6 * beta** 6 * x ** 4 + 6 * beta ** 4 * alpha ** 2 * y ** 2 + 14 * x ** 2 * alpha ** 4 *&
                     beta ** 2 + 3 * x ** 4 * y ** 4 * beta ** 2 - 6 * y ** 4 * beta ** 2 * x ** 2 + 3 * y ** 6 * &
                     beta ** 4 - 2 * y ** 6 * x ** 2 - 3 * y ** 6 * beta ** 2 - 6 * y ** 4 * beta ** 2 * x ** 4 * &
                     alpha ** 2 - 6 * beta ** 4 * x ** 4 * y ** 2 * alpha ** 2 + 9 * beta ** 2 * y ** 2 * x ** 4 * &
                     alpha ** 4 + 3 * beta ** 4 * y ** 4 * x ** 4 * alpha ** 2 + y ** 6 * x ** 4 + 3 * beta ** 6 * &
                     y ** 4 - 3 * beta ** 6 * y ** 2 - 2 * beta ** 6 * x ** 2 - 9 * y ** 2 * x ** 4 * alpha ** 4 - &
                     6 * y ** 6 * beta ** 4 * x ** 2 + 3 * y ** 6 * beta ** 4 * x ** 4 - 3 * y ** 6 * x ** 4 * &
                     beta ** 2 - 3 * alpha ** 4 * beta ** 2 * y ** 2 + 5 * alpha ** 6 * x ** 4 - 12 * alpha** 6 * &
                     x ** 2 + 3 * alpha ** 4 * y ** 2 + 3 * alpha ** 4 * beta ** 2) / dble(r ** 7) * dble(x) * &
                     sqrt(dble(1 - x ** 2)) / dble(alpha ** 2) / dble((y ** 2 - x ** 2 * y ** 2 + beta ** 2 - &
                     beta ** 2 *y ** 2 - beta ** 2 * x ** 2 + beta ** 2 * x ** 2 * y ** 2 + x ** 2 * alpha ** 2) ** 6) * &
                     dble((-(4 * beta ** 2 * x ** 2 * y ** 2 - beta ** 4 * y ** 4 + 2 * beta ** 4 * y ** 2 + 2 * &
                     beta ** 4 * x **2 + 5 * x ** 4 * alpha ** 4 - 4 * x ** 2 * alpha ** 2 * beta ** 2* y ** 2 + 4 * &
                     beta ** 2 * x ** 4 * y ** 2 * alpha ** 2 + 2 * beta ** 4 * x ** 2 * y ** 4 - 4 * x ** 4 * &
                     alpha ** 2 * beta ** 2 + 2* beta ** 4 * x ** 4 * y ** 2 - beta ** 4 * x ** 4 * y ** 4 - 4 * &
                     beta ** 4 * x ** 2 * y ** 2 + 4 * beta ** 2 * x ** 2 * alpha ** 2 - beta ** 4 * x ** 4 + 2 * &
                     y ** 4 * x ** 2 + 2 * y ** 4 * beta **2 - x ** 4 * y ** 4 - y ** 4 - 2 * beta ** 2 * y ** 2 - &
                     beta ** 4 - 9 * alpha ** 4 * x ** 2 + 2 * x ** 4 * y ** 4 * beta ** 2 - 4 *x ** 4 * y ** 2 * &
                     alpha ** 2 - 2 * x ** 4 * beta ** 2 * y ** 2 - 4 * y ** 4 * beta ** 2 * x ** 2 + 4 * y ** 2 * &
                     x ** 2 * alpha ** 2) / r ** 6 / alpha ** 2 / (y ** 2 - x ** 2 * y ** 2 + beta ** 2 - &
                     beta ** 2 * y ** 2 - beta ** 2 * x ** 2 + beta ** 2 * x ** 2 * y ** 2 + x ** 2 * &
                     alpha ** 2) ** 5) ** (-0.1D1 / 0.2D1))


                GradPhi = 0.3D1 * dble(4 * beta ** 2 * x ** 2 * y ** 2 - beta ** 4 * y** 4 + 2 * beta ** 4 * y ** 2 + &
                     2 * beta ** 4 * x ** 2 + 7 * x **4 * alpha ** 4 - 6 * x ** 2 * alpha ** 2 * beta ** 2 * y ** 2 + &
                     6* beta ** 2 * x ** 4 * y ** 2 * alpha ** 2 + 2 * beta ** 4 * x **2 * y ** 4 - 6 * x ** 4 * &
                     alpha ** 2 * beta ** 2 + 2 * beta ** 4 * x ** 4 * y ** 2 - beta ** 4 * x ** 4 * y ** 4 - 4 * &
                     beta ** 4 * x ** 2 * y ** 2 + 6 * beta ** 2 * x ** 2 * alpha ** 2 - beta ** 4 * x ** 4 + 2 * &
                     y ** 4 * x ** 2 + 2 * y ** 4 * beta ** 2 - x ** 4 *y ** 4 - y ** 4 - 2 * beta ** 2 * y ** 2 - &
                     beta ** 4 - 15 * alpha ** 4 * x ** 2 + 2 * x ** 4 * y ** 4 * beta ** 2 - 6 * x ** 4 * y ** 2 * &
                     alpha ** 2 - 2 * x ** 4 * beta ** 2 * y ** 2 - 4 * y ** 4 *beta ** 2 * x ** 2 + 6 * y ** 2 * &
                     x ** 2 * alpha ** 2) / dble(r ** 7) * dble(-1 + beta ** 2) * sqrt(dble(1 - y ** 2)) * &
                     dble(y) / dble(alpha ** 2) / dble((y ** 2 - x ** 2 * y ** 2 + beta ** 2 - beta ** 2 * y ** 2 - &
                     beta ** 2 * x ** 2 + beta ** 2 * x ** 2 * y ** 2 + x ** 2 * alpha ** 2) ** 6) * &
                     dble((-(4 * beta ** 2 * x ** 2 * y** 2 - beta ** 4 * y ** 4 + 2 * beta ** 4 * y ** 2 + &
                     2 * beta ** 4 * x ** 2 + 5 * x ** 4 * alpha ** 4 - 4 * x ** 2 * alpha ** 2 * beta ** 2 * y ** 2 + &
                     4 * beta ** 2 * x ** 4 * y ** 2 * alpha ** 2 + 2 * beta ** 4 * x ** 2 * y ** 4 - 4 * x ** 4 * &
                     alpha ** 2 * beta ** 2 + 2 * beta ** 4 * x ** 4 * y ** 2 - beta ** 4 * x ** 4 * y **4 - 4 * &
                     beta ** 4 * x ** 2 * y ** 2 + 4 * beta ** 2 * x ** 2 * alpha ** 2 - beta ** 4 * x ** 4 + 2 * &
                     y ** 4 * x ** 2 + 2 * y ** 4 *beta ** 2 - x ** 4 * y ** 4 - y ** 4 - 2 * beta ** 2 * y ** 2 - &
                     beta ** 4 - 9 * alpha ** 4 * x ** 2 + 2 * x ** 4 * y ** 4 * beta **2 - 4 * x ** 4 * y ** 2 * &
                     alpha ** 2 - 2 * x ** 4 * beta ** 2 * y** 2 - 4 * y ** 4 * beta ** 2 * x ** 2 + 4 * y ** 2 * &
                     x ** 2 * alpha ** 2) / r ** 6 / alpha ** 2 / (y ** 2 - x ** 2 * y ** 2 + beta** 2 - &
                     beta ** 2 * y ** 2 - beta ** 2 * x ** 2 + beta ** 2 * x **2 * y ** 2 + x ** 2 * &
                     alpha ** 2) ** 5) ** (-0.1D1 / 0.2D1)) * sqrt(dble(1 - x ** 2))


                ! dB/dt for the stretched dipole

                p = alpha

                dBdtTemp = alpha0 * cos(w * t) * w * p ** 3 * (0.9D1 * x ** 4 * y ** 6 * p **2 - 0.9D1 * x ** 2 * y ** 6 * p ** 2 &
                     - 0.3D1 * x ** 6 * y ** 6 * p ** 2 + 0.6D1 * x ** 6 * y ** 4 * p ** 2 - 0.96D2 * x ** 4 * y **2 * &
                     p ** 4 + 0.9D1 * x ** 4 * p ** 2 * y ** 2 + 0.60D2 * x ** 4 *p ** 6 * y ** 2 + 0.3D1 * y ** 2 * p ** 2 &
                     - 0.24D2 * x ** 6 * p ** 4 * y ** 4 - 0.3D1 * x ** 6 * p ** 2 * y ** 2 + 0.48D2 * x ** 6 * y ** 2 * &
                     p ** 4 + 0.60D2 * x ** 4 * p ** 8 * y ** 2 - 0.54D2 * x** 2 * p ** 8 * y ** 2 + 0.9D1 * x ** 2 * &
                     p ** 10 * y ** 2 - 0.6D1 * x ** 6 * p ** 8 * y ** 2 + 0.18D2 * x ** 2 * y ** 4 * p ** 2 -0.18D2 * &
                     x ** 4 * y ** 4 * p ** 2 + 0.3D1 * y ** 6 * p ** 6 * x ** 2 + 0.54D2 * x ** 2 * p ** 8 - 0.60D2 * &
                     x ** 4 * p ** 8 - 0.24D2* x ** 2 * p ** 4 + 0.48D2 * x ** 4 * p ** 4 - 0.2D1 * y ** 6 - 0.2D1 * x ** 6 &
                     - 0.6D1 * x ** 2 - 0.60D2 * x ** 4 * p ** 6 * y ** 4+ 0.12D2 * x ** 4 * p ** 8 * y ** 4 + 0.6D1 * &
                     x ** 4 * p ** 10 * y ** 2 + 0.30D2 * x ** 2 * p ** 6 * y ** 4 - 0.6D1 * x ** 2 * p **8 * y ** 4 + &
                     0.30D2 * x ** 6 * p ** 6 * y ** 4 - 0.6D1 * x ** 6 *p ** 8 * y ** 4 - 0.30D2 * x ** 6 * p ** 6 * &
                     y ** 2 - 0.15D2 * x ** 6 * p ** 10 * y ** 2 - 0.6D1 * x ** 4 * y ** 6 + 0.6D1 * x ** 2* y ** 6 - &
                     0.6D1 * x ** 6 * y ** 4 + 0.6D1 * x ** 6 * y ** 2 + 0.2D1 * x ** 6 * y ** 6 - 0.36D2 * x ** 4 * &
                     p ** 12 + 0.6D1 * x ** 6* p ** 8 + 0.20D2 * x ** 6 * p ** 12 + 0.3D1 * y ** 6 * p ** 2 - 0.6D1 * &
                     y ** 2 - 0.3D1 * x ** 4 * y ** 6 * p ** 6 + x ** 6 * y ** 6 * p ** 6 - 0.24D2 * x ** 6 * p ** 4 - &
                     0.1D1 * y ** 6 * p ** 6 + 0.48D2 * x ** 2 * y ** 2 * p ** 4 - 0.6D1 * y ** 4 * p ** 2 - 0.9D1 * &
                     x ** 2 * p ** 2 * y ** 2 - 0.30D2 * x ** 2 * p ** 6 * y ** 2 -0.24D2 * x ** 2 * y ** 4 * p ** 4 + &
                     0.48D2 * x ** 4 * y ** 4 * p ** 4 - 0.18D2 * x ** 2 * y ** 4 + 0.18D2 * x ** 4 * y ** 4 + 0.6D1* &
                     y ** 4 + 0.6D1 * x ** 4 - 0.18D2 * x ** 4 * y ** 2 + 0.18D2 * x** 2 * y ** 2 + 0.2D1) * (-0.1D1 * &
                     (-0.1D1 + 0.4D1 * x ** 4 * y ** 2 * p ** 4 - 0.2D1 * x ** 4 * p ** 2 * y ** 2 - 0.4D1 * x ** 4 *&
                     p ** 6 * y ** 2 - 0.2D1 * y ** 2 * p ** 2 - 0.4D1 * x ** 2 * y **4 * p ** 2 + 0.2D1 * x ** 4 * y ** 4 &
                     * p ** 2 - 0.9D1 * x ** 2 * p ** 8 + 0.5D1 * x ** 4 * p ** 8 - 0.1D1 * y ** 4 * p ** 4 + 0.4D1* &
                     x ** 2 * p ** 4 - 0.4D1 * x ** 4 * p ** 4 + 0.2D1 * x ** 2 + 0.2D1 * y ** 2 - 0.4D1 * x ** 2 * y ** 2 &
                     * p ** 4 + 0.2D1 * y ** 4 *p ** 2 + 0.4D1 * x ** 2 * p ** 2 * y ** 2 + 0.4D1 * x ** 2 * p **6 * y ** 2 &
                     + 0.2D1 * x ** 2 * y ** 4 * p ** 4 - 0.1D1 * x ** 4 * y ** 4 * p ** 4 + 0.2D1 * x ** 2 * y ** 4 - &
                     0.1D1 * x ** 4 * y ** 4 - 0.1D1 * y ** 4 - 0.1D1 * x ** 4 + 0.2D1 * x ** 4 * y ** 2 - 0.4D1 * x ** 2 &
                     * y ** 2) * p ** 4 / r ** 6 / (x ** 2 * y ** 2 - 0.1D1 * x ** 2 * p ** 2 * y ** 2 - 0.1D1 * x ** 2 +&
                     x ** 2 * p ** 4 - 0.1D1 * y ** 2 + y ** 2 * p ** 2 + 0.1D1) ** 5) ** (-0.1D1 / 0.2D1) / r ** 6 / &
                     (x ** 2 * y ** 2 - 0.1D1 * x ** 2 * p ** 2 * y ** 2 -0.1D1 * x ** 2 + x ** 2 * p ** 4 - 0.1D1 * &
                     y ** 2 + y ** 2 * p **2 + 0.1D1) ** 6


                dBdt_III(iPoint,iR,iPhi) = dBdtTemp
                !                write(*,*) dBdtTemp

                ! Gradient B
                GradB_VIII(1,iPoint, iR, iPhi) = DipoleFactor*GradR
                GradB_VIII(2,iPoint, iR, iPhi) = DipoleFactor*GradTheta
                GradB_VIII(3,iPoint, iR, iPhi) = DipoleFactor*GradPhi               


                ! drift Velocity components
                GradBCrossB_VIII(1,iPoint,iR,iPhi) = Vr * (DipoleFactor)**3
                GradBCrossB_VIII(2,iPoint,iR,iPhi) = 1./(L_I(iR))*Vtheta * (DipoleFactor)**3
                GradBCrossB_VIII(3,iPoint,iR,iPhi) = 1./(L_I(iR)* cos(Lat)) *Vphi * (DipoleFactor)**3


                mag = DipoleFactor*(sqrt(Br**2+Btheta**2+Bphi**2))


                bFieldMagnitude_III(iPoint,iR,iPhi) = mag
                Length_III(iPoint,iR,iPhi) = stretched_dipole_length(L_I(iR), LatMin,Lat,Phi_I(iPhi), alpha, beta)  
                Lat = Lat + dLat 


             end do
          end do
       end do

    end select

    do iPhi =1, nPhi
       do iR =1, nR 
          do iPoint = 1, nPoint-1
             dLength_III(iPoint,iR,iPhi) = Length_III(iPoint+1,iR,iPhi) - Length_III(iPoint,iR,iPhi)
          end do
       end do
    end do

  end subroutine initialize_b_field
  !==================================================================================
  subroutine find_mirror_points (nPoint, PitchAngle, bField_I, bMirror,iMirror_II)

    integer, intent(in) :: nPoint                ! Number of points along the field line
    real, intent(in)    :: PitchAngle            ! Pitch angle values
    real, intent(in)    :: bField_I(nPoint)      ! Magnetic field values
    real, intent(out)   :: bMirror               ! B magnitude at mirror points    
    real                :: bMin,bMax                  ! Minimum value of magnetic field along a field line
    integer             :: iMinB,iMaxB           ! Location of minimum B
    integer,intent(out) :: iMirror_II(2)  ! Location of each mirror point for all pitch angles
    integer             :: i_I(1),j_I(1)
    integer             :: iPoint, iPitch
    real, parameter     :: cTiny = 1.e-6
    real                :: PALossCone
    !----------------------------------------------------------------------------------
    i_I  = minloc(bField_I)
    iMinB= i_I(1)
    bMin = bField_I(iMinB)

    j_I = maxloc(bField_I)
    iMaxB = j_I(1)
    bMax = bField_I(iMaxB)

    PALossCone = asin(sqrt( bField_I(nPoint/2)/bMax))

    if (PitchAngle <= PALossCone) then
       bMirror = bMin/(sin(PALossCone))**2
    else if (PitchAngle == 0.0) then   
       bMirror = bMin/(sin(PitchAngle+cTiny))**2 
    else 
       bMirror = bMin/(sin(PitchAngle))**2 
    end if

    do iPoint = iMinB,1, -1
       if (bField_I(iPoint) >= bMirror ) then 
          iMirror_II(1) = (iPoint + 1)
          EXIT
       end if
    end do

    do iPoint = iMinB, nPoint
       if (bField_I(iPoint) >= bMirror ) then 
          iMirror_II(2) = (iPoint - 1)
          EXIT
       end if
    end do

  end subroutine find_mirror_points

  !==================================================================================

  subroutine second_adiabatic_invariant(nPoint, iMirror_I, bMirror, bField_I, dLength_I,L, SecondAdiabInv)
    !\
    ! Calculate integral of sqrt((B-Bm)/Bm) between the mirror points using the
    ! trapezoidal rule
    !/

    integer             :: nPoint
    integer, intent(in) :: iMirror_I(2)
    real, intent(in)    :: bMirror  
    real, intent(in)    :: dLength_I(nPoint-1)
    real, intent(in)    :: bField_I(nPoint) 
    real, intent(out)   :: SecondAdiabInv
    real, intent(in)    :: L
    real                :: InvL
    integer             :: iPoint, iFirst, iLast
    real                :: DeltaS1, DeltaS2, b1, b2, Coeff
    !----------------------------------------------------------------------------------
    iFirst = iMirror_I(1)
    iLast  = iMirror_I(2)
    SecondAdiabInv = 0.0


    if (iFirst > iLast) RETURN

    InvL = 1.0/L
    Coeff = InvL/sqrt(bMirror)

    DeltaS1 = abs((bMirror-bField_I(iFirst))*(dLength_I(iFirst-1))/(bField_I(iFirst-1)-bField_I(iFirst)))


    SecondAdiabInv= SecondAdiabInv + Coeff*(2./3.)*DeltaS1*sqrt(bMirror-bField_I(iFirst))

    do iPoint = iFirst, iLast-1
       b1 = bField_I(iPoint)
       b2 =  bField_I(iPoint+1)
       SecondAdiabInv = SecondAdiabInv + Coeff*(2./3.)*dLength_I(iPoint)/(b1 - b2) &
            *( sqrt(bMirror  - b2)**3 - sqrt(bMirror  - b1)**3 )  

    end do

    DeltaS2 = abs((bMirror-bField_I(iLast))*(dLength_I(iLast))/(bField_I(iLast+1)-bField_I(iLast)))
    SecondAdiabInv= SecondAdiabInv + Coeff*(2./3.)*DeltaS2*(sqrt(bMirror-bField_I(iLast)))

  end subroutine second_adiabatic_invariant

  !==================================================================================
  subroutine half_bounce_path_length(nPoint, iMirror_I, bMirror, bField_I, dLength_I,L, HalfPathLength,Sb)

    !\
    ! Calculate integral of ds/sqrt((B-Bm)/Bm) between the mirror points using 
    !/

    integer             :: nPoint
    integer, intent(in) :: iMirror_I(2)
    real, intent(in)    :: bMirror  
    real, intent(in)    :: dLength_I(nPoint-1)
    real, intent(in)    :: bField_I(nPoint) 
    real, intent(out)   :: HalfPathLength,Sb
    real, intent(in)    :: L
    real                :: Inv2L
    integer             :: iPoint, iFirst, iLast
    real                :: DeltaS1, DeltaS2,b1,b2, Coeff,CoeffSb
    !----------------------------------------------------------------------------------
    iFirst = iMirror_I(1)
    iLast  = iMirror_I(2)
    HalfPathLength = 0.0
    Sb = 0.0


    if (iFirst >iLast) RETURN

    Inv2L = 1.0/(2.*L)
    Coeff = Inv2L*sqrt(bMirror)
    CoeffSb = sqrt(bMirror)



    DeltaS1 = abs((bMirror-bField_I(iFirst))*(dLength_I(iFirst-1))/(bField_I(iFirst-1)-bField_I(iFirst)))
    HalfPathLength = HalfPathLength + Coeff*2.*DeltaS1/(sqrt(bMirror-bField_I(iFirst)))
    Sb = Sb + CoeffSb*2.*DeltaS1/(sqrt(bMirror-bField_I(iFirst)))

    do iPoint = iFirst, iLast-1
       b1 = bField_I(iPoint)
       b2 = bField_I(iPoint+1)


       HalfPathLength = HalfPathLength + Coeff*2.*dLength_I(iPoint)/(b1 - b2) &
            *( sqrt(bMirror  - b2) - sqrt(bMirror  - b1) )
       Sb = Sb + CoeffSb*2.*dLength_I(iPoint)/(b1 - b2) &
            *( sqrt(bMirror  - b2) - sqrt(bMirror  - b1) )

    end do

    DeltaS2 = abs((bMirror-bField_I(iLast))*(dLength_I(iLast))/(bField_I(iLast+1)-bField_I(iLast)))
    HalfPathLength= HalfPathLength + Coeff*2.*DeltaS2/(sqrt(bMirror-bField_I(iLast)))
    Sb = Sb + CoeffSb*2.*DeltaS2/(sqrt(bMirror-bField_I(iLast)))


  end subroutine half_bounce_path_length

  !================================================================================== 
  subroutine get_gradB0(nPoint,nR,nPhi,L_I,Phi_I,GradB0R_III, GradB0Theta_III, GradB0Phi_III)

    implicit none

    integer, intent(in)    :: nPoint                              ! Number of points along the field line
    integer, intent(in)    :: nR                                  ! Number of points in the readial direction
    integer, intent(in)    :: nPhi                                ! Number of points in the azimuthal direction
    real,    intent(in)    :: L_I(nR) 
    real,    intent(in)    :: Phi_I(nPhi) 
    real, dimension(nPoint,nR,nPhi), intent(out) :: GradB0R_III, GradB0Theta_III, GradB0Phi_III

    !Local Variables
    real    :: LatMax, LatMin,dLat, Lat
    real    :: y, x, r
    integer :: iPhi, iR, iPoint
    

    !----------------------------------------------------------------------------------
    
    
    do iPhi =1, nPhi
       y = cos(Phi_I(iPhi))
       do iR =1, nR 
          LatMax =  acos(sqrt(1./L_I(iR)))
          LatMin = -LatMax
          dLat   = (LatMax-LatMin)/(nPoint-1)
          Lat = LatMin
          do iPoint = 1, nPoint
             x = sin(Lat)
             r = L_I(iR) * (cos(Lat))**2
                          
             GradB0R_III(iPoint,iR,iPhi) = dble((4/r**6 * x**2 + 1/r**6 * (1 - x**2))**(-0.1D1/0.2D1) * &
                  (-24/r**7 * x**2 - 6/r**7 * (1 - x**2)))/0.2D1
             GradB0Theta_III(iPoint,iR,iPhi) = -0.3D1/dble(r**7) * dble((4/r**6 * x**2 + 1 / r **6 * &
                  (1 - x**2))**(-0.1D1/0.2D1)) * dble(x) * sqrt(dble(1 - x**2))
             GradB0Phi_III(iPoint,iR,iPhi) = 0.0
             
          end do
       end do
    end do
    
  end subroutine get_gradB0
  !================================================================================== 
  real function dipole_length(L, LatMin, LatMax)

    implicit none

    real               :: L                       ! L shell value
    real               :: LatMin                  ! Minimum Latitude
    real               :: LatMax                  ! Maximum Latitude
    real               :: x, y
    !----------------------------------------------------------------------------------

    x = sin(LatMax)
    y = sin(LatMin)

    dipole_length = abs(L*(0.5*Sqrt(3.0*x**2.0+1.0)*x + asinh(sqrt(3.0)*x)/2.0/sqrt(3.0)) &
         - L*(0.5*Sqrt(3.0*y**2.0+1.0)*y + asinh(sqrt(3.0)*y)/2.0/sqrt(3.0))) 

  end function dipole_length

  !==================================================================================

  real function stretched_dipole_length(L,LatMin,LatMax,Phi,alpha,beta)


    use ModHeidiSize, only: nPoint

    implicit none
    real               :: L                       ! L shell value   
    real               :: LatMin                  ! Minimum Latitude                                               
    real               :: LatMax                  ! Maximum Latitude
    real               :: dLat, dSdTheta(nPoint)
    real               :: x(nPoint), y, Phi,Lat
    real               :: alpha, beta              ! Stretching parameters
    real               :: alpha2,alpha4,gamma
    real               :: cos2Lat,cosLat,a,b,c
    integer            :: i
    !----------------------------------------------------------------------------------
    y = cos(Phi)
    stretched_dipole_length = 0.0
    dLat   = (LatMax-LatMin)/(nPoint-1)
    Lat = LatMin

    x(1) = sin(LatMin)

    dSdTheta(1) = (sqrt(-dble((-9 * x(1) ** 4 * beta ** 2 + 9 * beta ** 2 * &
         x(1) ** 2 + alpha ** 2 + 9 * x(1) ** 4 * y ** 2 * beta ** 2 - 9 * x(1) ** 4 * y **2 - &
         9 * beta ** 2 * y ** 2 * x(1) ** 2 + 9 * y ** 2 * x(1) ** 2 - 6 * alpha ** 2 * & 
         x(1) ** 2 + 9 * x(1) ** 4 * alpha ** 2) * (y ** 2 - y ** 2 * x(1) ** 2 + beta ** 2 - & 
         beta ** 2 * y ** 2 - beta ** 2 * x(1) ** 2 + beta ** 2 * y ** 2 * x(1) ** 2 + &
         alpha ** 2 * x(1) ** 2) ** 2 * (1 - x(1) ** 2) * L ** 2 / (-y ** 2 + beta ** 2 * &
         y ** 2 - beta ** 2) / alpha ** 2)))

    do i = 2, nPoint
       Lat = Lat + dLat
       x(i) = sin(Lat)
       dSdTheta(i) = (sqrt(-dble((-9 * x(i) ** 4 * beta ** 2 + 9 * beta ** 2 * &
            x(i) ** 2 + alpha ** 2 + 9 * x(i) ** 4 * y ** 2 * beta ** 2 - 9 * x(i) ** 4 * y **2 - &
            9 * beta ** 2 * y ** 2 * x(i) ** 2 + 9 * y ** 2 * x(i) ** 2 - 6 * alpha ** 2 * & 
            x(i) ** 2 + 9 * x(i) ** 4 * alpha ** 2) * (y ** 2 - y ** 2 * x(i) ** 2 + beta ** 2 - &
            beta ** 2 * y ** 2 - beta ** 2 * x(i) ** 2 + beta ** 2 * y ** 2 * x(i) ** 2 + &
            alpha ** 2 * x(i) ** 2) ** 2 * (1 - x(i) ** 2) * L ** 2 / (-y ** 2 + beta ** 2 * &
            y ** 2 - beta ** 2) / alpha ** 2)))

       stretched_dipole_length = stretched_dipole_length + 0.5*(dSdTheta(i)+dSdTheta(i-1))*dLat

    end do
  end function stretched_dipole_length

  !==================================================================================

  real function asinh(x)              

    implicit none

    real, intent(in) :: x
    !----------------------------------------------------------------------------------

    asinh = log(x+sqrt(x**2+1.0))

  end function asinh

  !==================================================================================

  real function second_adiab_invariant(x)		
    !\
    ! This is an approximate formula for second adiabatic invariant 
    ! = the integral of sqrt(1 - B/Bm)ds for a dipole B field. 
    !	function I(mu) taken from Ejiri, JGR,1978
    !/

    real, intent(in)     :: x ! cosine of the equatorial pitch angle
    real                 :: y, alpha, beta, a1, a2, a3, a4
    real, parameter      :: Pi = 3.141592654
    !----------------------------------------------------------------------------------

    y=sqrt(1.-x*x)
    alpha=1.+alog(2.+sqrt(3.))/2./sqrt(3.)
    beta=alpha/2.-Pi*sqrt(2.)/12.
    a1=0.055
    a2=-0.037
    a3=-0.074
    a4=0.056
    second_adiab_invariant = &
         2.*alpha*(1.-y)+2.*beta*y*alog(y)+4.*beta*(y-sqrt(y))+  &
         3.*a1*(y**(1./3.)-y)+6.*a2*(y**(2./3.)-y)+6.*a4*(y-y**(4./3.))  &
         -2.*a3*y*alog(y)

  end function second_adiab_invariant

  !==================================================================================

  real function  analytic_h(x)		
    !\
    ! This is an approximate formula for second adiabatic invariant 
    ! = the integral of ds/sqrt(1 - B/Bm) for a dipole B field. 
    !	function f(y) taken from Ejiri, JGR,1978
    !/

    real, intent(in)     :: x ! cosine of the equatorial pitch angle 
    real                 :: y, alpha, beta, a1, a2, a3, a4
    real, parameter      :: Pi = 3.141592654
    !----------------------------------------------------------------------------------

    y=sqrt(1-X*X)
    alpha=1.+alog(2.+sqrt(3.))/2./sqrt(3.)
    beta=alpha/2.-Pi*sqrt(2.)/12.
    a1=0.055
    a2=-0.037
    a3=-0.074
    a4=0.056
    analytic_h = alpha-beta*(y+sqrt(y))+a1*y**(1./3.)+a2*y**(2./3.)+  &
         a3*y+a4*y**(4./3.)

  end function analytic_h

  !=====================================================================                      
  subroutine get_quadratic_root(a,b,c,x,nroot)                                                
    !\                                                                                        
    ! Solve the  quadratic equation: ax**2 + bx +c = 0                                        
    !/                                                                                        

    implicit none                                                                             

    real    :: a, b, c            ! coefficients                                              
    real    :: discriminant       ! discriminant                                              
    complex :: x(2)               ! roots                                                     
    integer :: nroot              ! number of roots                                           

    !--------------------------------------------------------------------                     

    if (a == 0.0) then                                                                        
       if(b == 0.0) then                                                                      
          nroot = 0               ! Nothing to solve                                          
       else                                                                                   
          nroot =1                                                                            
          x(1) = cmplx(-c/b, 0.0) ! Linear equation                                           
       end if
    else                                                                                      

       !Start solving the quadratic equation                                                  
       nroot  = 2                                                                             
       discriminant = b*b -4.*a*c                                                             
       if (discriminant >= 0.0) then                                                          
          x(1) = cmplx((-b + sqrt(discriminant))/(2.*a), 0.0)                                 
          x(2) = cmplx((-b - sqrt(discriminant))/(2.*a), 0.0)                                 
       else                                                                                   
          x(1) = cmplx(-b/(2.*a), sqrt(-discriminant)/(2.*a))                                  
          x(2) = cmplx(-b/(2.*a), -sqrt(-discriminant)/(2.*a))                                 
       end if
    end if


  end subroutine get_quadratic_root
  !=====================================================================                      

  subroutine get_cubic_root(a,b,c,d,x,nroot)                                                  

    !\                                                                                        
    ! Solve a cubic equation ax**3 + bx**2 + cx + d = 0                                       
    !/                                                                                        

    use ModNumConst,   ONLY: cPi

    implicit none

    real    :: a,b,c,d       ! coefficients                                                   
    complex :: x(3)          ! roots                                                          
    integer :: nroot         ! number of roots                                                
    real    :: p,q,discriminant                                                               
    real    :: alpha, beta, phi                                                               
    real    :: y1,y2,y3,y2real,y2imag                                                         
    real    :: factor                                                                         
    real    :: u,v                                                                            

    !---------------------------------------------------------------------                      

    y1 = 0.0
    y2 = 0.0
    y3 = 0.0
    y2real = 0.0

    if (a == 0.0) then                                                                        
       call get_quadratic_root(b,c,d,x,nroot)                                                 
       return                                                                                 
    end if

    nroot = 3                                                                                 

    !\                                                                                        
    ! define p and q and the discriminant                                                     
    ! definitions in Tuma and Walsh 'Engineering mathematics handbook'                        
    ! x = y -b/3a ; y**3+py+q = 0                                                             
    ! If discriminant > 0 => one real and two complex roots                                   
    !    discriminant = 0 =>  there are three real roots of which at least two are equal      
    !    discriminant < 0 =>  there are three real unequal roots                              
    !/                                                                                        

    p = (1./3.)*(3.*(c/a)-(b/a)*(b/a))
    q = (1./27.)*(2.*(b/a)*(b/a)*(b/a) - 9.*(b/a)*(c/a) + 27.*(d/a))
    discriminant = (p*p*p)/27. + (q*q)/4.

    if (discriminant < 0.0) then 
       phi = acos((-q/2.)/(sqrt(abs(p*p*p)/ 27.)))
       alpha =  2.* sqrt(abs(p)/3.)                                                           
       y1 = alpha*cos(phi/3.)                                                                 
       y2 = -alpha*cos((phi + cPi)/3.)                                                        
       y3 = -alpha*cos((phi - cPi)/3.)
    else                                                                                      
       alpha = -(q/2.) + sqrt(discriminant)                                                    
       beta  = -(q/2.) - sqrt(discriminant)                                                    
       u = abs(alpha)**(1./3.)                                                                
       v = abs(beta)**(1./3.)                                                                 
       if (alpha < 0.0) u = -u                                                                
       if (beta < 0.0)  v = -v                                                                
       y1 = u + v                                                                             
       y2real = -(u + v)/2.   
       y2imag = (u -v) * sqrt(3.)/2                                                           
    end if

    factor = b/(a*3.)
    y1 = y1 - factor
    y2 = y2 -factor
    y3 = y3 -factor
    y2real = y2real -factor

    if (discriminant < 0.0) then                                                              
       x(1) = cmplx(y1, 0.0)                                                                  
       x(2) = cmplx(y2, 0.0)                                                                  
       x(3) = cmplx(y3, 0.0)                                                                  
    else if (discriminant == 0.0) then                                                             
       x(1) = cmplx(y1, 0.0)                                                                  
       x(2) = cmplx(y2real, 0.0)                                                              
       x(3) = cmplx(y2real, 0.0)                                                              
    else                                                                                      
       x(1) = cmplx(y1, 0.0)                                                                  
       x(2) = cmplx(y2real, y2imag)                                                           
       x(3) = cmplx(y2real, -y2imag)                                                          
    end if


  end subroutine get_cubic_root

  !==================================================================================

  subroutine test_general_b

    integer, parameter   :: nPoint = 10001
    integer, parameter   :: nPitch = 1
    integer, parameter   :: nR = 1
    integer, parameter   :: nPhi =1
    real                 :: L_I(nR)
    real                 :: Phi_I(nPhi)                      ! Phi values
    real                 :: bFieldMagnitude_III(nPoint,nR,nPhi) ! Magnitude of magnetic field 
    real                 :: RadialDistance_III(nPoint,nR,nPhi)
    real                 :: dLength_III(nPoint-1,nR,nPhi)      ! Length interval between i and i+1  
    real                 :: Length_III(nPoint,nR,nPhi) 
    real                 :: GradBCrossB_VIII(3,nPoint,nR,nPhi),GradB_VIII(3,nPoint,nR,nPhi)
    real                 :: PitchAngle
    real                 :: Ds_I(nPoint)
    real                 :: SecondAdiabInv, IntegralBAnalytic
    real                 :: HalfPathLength, IntegralHAnalytic,Sb
    real                 :: bMirror
    integer              :: iMirror_I(2)
    real, parameter      :: Pi = 3.141592654
    real                 :: Percent1, Percent2
    integer              :: iPoint
    real                 :: dBdt_III(nPoint,nR,nPhi)
    !----------------------------------------------------------------------------------
    !open (unit = 2, file = "Convergence_nonuniform.dat")
    !write (2,*)'Numerical values for the second adiabatic invariant integration btw a mirror point and eq'
    !write (2,*)'nPoint  IntegralBAnalytic   2nd_orderB  IntegralHAnalytic  2nd_orderH Percent1, Percent2'

    L_I(1) = 10. 
    Phi_I(1) = 1.0
    PitchAngle = Pi/10.

    do iPoint = 101, nPoint,100

       call initialize_b_field(L_I(1), Phi_I, nPoint, nR, nPhi, bFieldMagnitude_III, &
            RadialDistance_III,Length_III, dLength_III,GradBCrossB_VIII,GradB_VIII,dBdt_III)

       IntegralBAnalytic = second_adiab_invariant(cos(PitchAngle))
       IntegralHAnalytic = analytic_h(cos(PitchAngle))
       call find_mirror_points (iPoint,  PitchAngle, bFieldMagnitude_III, &
            bMirror,iMirror_I)

       call second_adiabatic_invariant(iPoint, iMirror_I, bMirror, bFieldMagnitude_III, Ds_I,L_I(1), SecondAdiabInv)

       Percent1 = 100*abs(IntegralBAnalytic - SecondAdiabInv)/IntegralBAnalytic

       call half_bounce_path_length(iPoint, iMirror_I, bMirror,  bFieldMagnitude_III, Ds_I,L_I(1), HalfPathLength,Sb)


       Percent2 = 200*abs(IntegralHAnalytic - HalfPathLength)/(IntegralHAnalytic+HalfPathLength)

       !       write (2,*) iPoint, IntegralBAnalytic, SecondAdiabInv,IntegralHAnalytic,&
       !            HalfPathLength ,Percent1, Percent2

       !       write (*,*) iPoint, IntegralBAnalytic, SecondAdiabInv,IntegralHAnalytic,&
       !            HalfPathLength ,Percent1, Percent2


    end do

    !close(2)

  end subroutine test_general_b
  !==================================================================================

end module ModHeidiBField

