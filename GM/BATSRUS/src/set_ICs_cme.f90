!^CFG COPYRIGHT UM
!^CFG FILE NOT SIMPLE
!---------------------------------------------------------------------------
!Subroutine             set_ICs_cme
!---------------------------------------------------------------------------
!\
! Calculates magnetic field, pressure and density for a coronal flux rope 
! capable of self-similar expansion and erupting in a CME.
! The analytical solution is taken from Gibson and Low 
! Astrophysical Journal, Vol 493, p. 460.
!
! Written by Chip Manchester Jan 18 2001
!/
!   Bug fixes
!   March 18       dpres_1dr1 = cme_a1*(dA2dr*dr2dr1 + dA2dth*dth2dr1)
!   change to..... dpres_1dr1 = a1scl*(dA2dr*dr2dr1 + dA2dth*dth2dr1)
!   without fix, used for runs 12, 13, 14
!   Feb   07, 2002       Br1 is changed to Br1**2 (line 323 thanks to Ilia)
!===========================================================================
subroutine set_ICs_cme
  use ModMain, ONLY : nI,nJ,nK,gcn,globalBLK
  use ModVarIndexes,ONLY:rho_,rhoUx_,rhoUy_,rhoUz_,&
       Bx_,By_,Bz_,P_
  use ModAdvance, ONLY : &
       State_VGB
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK
  use ModPhysics
  implicit none

  integer :: i, j, k
  real :: x, y, z,   x_1, y_1, z_1,   x_2, y_2, z_2
  real :: r,   cos_theta,  sin_theta,  cos_phi,  sin_phi
  real :: r_1, cos_theta1, sin_theta1, cos_phi1, sin_phi1, lambda
  real :: r_2, cos_theta2, sin_theta2, cos_phi2, sin_phi2
  real :: dr2dr1, dth2dr1, cos_thmax, sin_thmax, dsin_thmaxdr
  real :: Br,  Btheta,  Bphi, Br1, Btheta1, Bphi1 
  real :: Br2, Btheta2, Bphi2, Bx_1, By_1, Bz_1
  real :: Br_r0, Btheta_r0, Psurface, dPsurfdr
  real :: dBr1dr, dBtheta1dr, dBphi1dr, dBr_r0dr, dBtheta_r0dr 
  real :: dBr2dr2, dBr2dth2, dBth2dr2, dBth2dth2, dBphi2dr2, dBphi2dth2
  real :: dA1dr, dA1dth,  d2A1dr2, d2A1drdth 
  real :: dA1dr_r0, dA1dth_r0, d2A1dr_r0dr, d2A1dth_r0dr
  real :: A2, dA2dr, dA2dth,  d2A2dr2, d2A2drdth, d2A2dth2
  real :: pres_1, dpres_1dr1, F_grav, alpha0, ga0r0 
  real :: denom, denom_r0, vr, delta
  !
  !  PARAMETER LIST: cme_a, cme_r1, cme_r0, cme_a1, cme_rho1, cme_rho2, B1_dim,
  !                  RHOsun,Vscl
  !  Definition of Parameters used for the initial state
  !   cme_a    = contraction distance as in   r --> r -a
  !   cme_r1   = distance of flux rope from sun center = 1.2
  !   cme_r0   = radius of flux rope
  !   cme_a1   = constant for setting pressure in flux rope
  !   Rscl     = 1.0  scaled radius of the sun
  !   RHOscl   = 1.0  scaled density of RHOsun
  !   SSPscl   = 1.0  scaled soundspeed of the sun
  !   rho1scl  = uniform backround density of the solution before contraction
  !   rho2scl  = background powerlaw density added to after contraction
  !   B1scl    = magnetic field strength parameter of the flux rope
  !   Gsun     = gravitational acceleration at sun's surface = 2.734e4 cm/s**2
  !   Vscl     = V/SSPsun     
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////
  !=====================================================================
  !  if(iProc==0)then
  !     write(6,*) 'cme_init called'
  !     write(6,*) 'Scale height = ',P_height
  !     write(6,*) 'B1_dim = ',B1_dim
  !  endif
  !
  delta = 0.1
  do i=1-gcn,nI+gcn
     do j=1-gcn,nJ+gcn
        do k=1-gcn,nK+gcn
           !CALCULATE CELL CENTER FOR GLOBAL Cartesian COORDINATES 
           x = x_BLK(i,j,k,globalBLK)
           y = y_BLK(i,j,k,globalBLK)
           z = z_BLK(i,j,k,globalBLK)
           !CALCULATE CELL CENTER FOR GLOBAL SPHERICAL COORDINATES 
           r = sqrt(x**2 + y**2 + z**2)
           cos_theta = z/r
           sin_theta = sqrt(x**2 + y**2)/r
           cos_phi   = x/sqrt(x**2 + y**2)
           sin_phi   = y/sqrt(x**2 + y**2)
           if (abs(r) .lt. delta) then 
              r =  delta
              x =  delta*sin_theta*cos_phi
              y =  delta*sin_theta*sin_phi
              z =  delta*cos_theta
           endif
           !CALCULATE CELL CENTER FOR TRANSFORMED SPHERICAL COORDINATES 
           !self-similar transformation of variables r --> r - a
           lambda = r + cme_a
           r_1 = lambda
           cos_theta1 = cos_theta
           sin_theta1 = sin_theta 
           cos_phi1   = cos_phi
           sin_phi1   = sin_phi
           !CALCULATE CELL CENTER FOR TRANSFORMED Cartesian COORDINATES 
           x_1 = lambda*sin_theta1*cos_phi1
           y_1 = lambda*sin_theta1*sin_phi1
           z_1 = lambda*cos_theta1
           !             
           !POTENTIAL BIPOLE MAGNETIC FIELD
           denom = ( (cme_r0**2 -cme_r1**2)**2 + (cme_r1*r_1)**2 +  & 
                2.0*r_1*cme_r1*(cme_r0**2 -cme_r1**2)*cos_theta1 )

           dA1dth = -sin_theta1 +  & 
                r_1*cme_r1*sin_theta1/ &
                (cme_r0*sqrt(r_1**2 + cme_r1**2 - 2.0*r_1*cme_r1*cos_theta1))   + &
                !
           r_1*(cme_r0**2 - 2.0*cme_r1**2)*sin_theta1/(cme_r0*sqrt(denom)) - &
                !
           r_1*cme_r1*(cme_r0**2-cme_r1**2)*sin_theta1 * &
                ( cme_r1*(r_1**2 + cme_r1**2 - cme_r0**2)   + &
                r_1*(cme_r0**2-2.0*cme_r1**2)*cos_theta1 )/ & 
                (cme_r0*denom*sqrt(denom))     

           dA1dr = (r_1-cme_r1*cos_theta1) / &
                ( cme_r0*sqrt(r_1**2 + cme_r1**2 -2.0*r_1*cme_r1*cos_theta1) )  - &
                !
           (2.0*r_1*cme_r1 + (cme_r0**2 -2.0*cme_r1**2)*cos_theta1) / &
                (cme_r0*sqrt(denom))                                            + &
                !
           ( cme_r1*(r_1**2 +cme_r1**2 -cme_r0**2) + & 
                r_1*(cme_r0**2 -2.0*cme_r1**2)*cos_theta1 ) * & 
                (r_1*cme_r1**2 + cme_r1*(cme_r0**2 -cme_r1**2)*cos_theta1) / &
                (cme_r0*denom*sqrt(denom))
           !
           d2A1dr2 = 1.0/(cme_r0*sqrt(r_1**2 +cme_r1**2 -2.0*r_1*cme_r1*cos_theta1))-&
                !
           (r_1 - cme_r1*cos_theta1)**2 / &
                (cme_r0*(r_1**2 +cme_r1**2 -2.0*r_1*cme_r1*cos_theta1)**(1.5))  - &
                !
                2.0*cme_r1/(cme_r0*sqrt(denom))                                 + &
                !
           (2.0*r_1*cme_r1 + (cme_r0**2-2.0*cme_r1**2)*cos_theta1) * &
                (r_1*cme_r1**2 + cme_r1*(cme_r0**2-cme_r1**2)*cos_theta1) / &
                (cme_r0*denom*sqrt(denom))                                      + &
                !
           ( (2.0*r_1*cme_r1 + (cme_r0**2 -2.0*cme_r1**2)*cos_theta1) * &
                (r_1*cme_r1**2 + cme_r1*(cme_r0**2 -cme_r1**2)*cos_theta1)+ &
                r_1*(cme_r1*(r_1**2 +cme_r1**2 -cme_r0**2) + &
                r_1*(cme_r0**2 -2.0*cme_r1**2)*cos_theta1) ) / &
                (cme_r0*denom*sqrt(denom))                                      - &
                !
                3.0*(cme_r1*(r_1**2 +cme_r1**2 -cme_r0**2) + &
                r_1*(cme_r0**2 -2.0*cme_r1**2)*cos_theta1)* &
                (r_1*cme_r1**2 + cme_r1*(cme_r0**2 -cme_r1**2)*cos_theta1)**2 / &
                (cme_r0*sqrt(denom)*denom**2) 
           !
           d2A1drdth = cme_r1*sin_theta1 / & 
                (cme_r0*sqrt(r_1**2 +cme_r1**2 -2.0*r_1*cme_r1*cos_theta1))     - &
                !
           r_1*cme_r1*sin_theta1*(r_1 -cme_r1*cos_theta1) / &
                (cme_r0*(r_1**2 +cme_r1**2 -2.0*r_1*cme_r1*cos_theta1)**(1.5))  + &
                !
           (cme_r0**2 -2.0*cme_r1**2)*sin_theta1/(cme_r0*sqrt(denom))      - &
                !
           r_1*sin_theta1*(cme_r0**2 -2.0*cme_r1**2) * & 
                (r_1*cme_r1**2 +cme_r1*(cme_r0**2 -cme_r1**2)*cos_theta1) / &
                (cme_r0*denom*sqrt(denom))                                      - &
                !
           sin_theta1*cme_r1*(cme_r0**2 -cme_r1**2) * &
                ( (cme_r1*(r_1**2 +cme_r1**2 -cme_r0**2) + & 
                r_1*(cme_r0**2 -2.0*cme_r1**2)*cos_theta1) - &
                r_1*(2.0*r_1*cme_r1 +(cme_r0**2 -2.0*cme_r1**2)*cos_theta1) ) / &
                (cme_r0*denom*sqrt(denom))                                      + &
                !
                3.0*sin_theta1*r_1*cme_r1*(cme_r0**2 -cme_r1**2) * &
                (cme_r1*(r_1**2 +cme_r1**2 -cme_r0**2) + &
                r_1*(cme_r0**2 -2.0*cme_r1**2)*cos_theta1)     * &
                (r_1*cme_r1**2 + cme_r1*(cme_r0**2 -cme_r1**2)*cos_theta1) / &
                (cme_r0*sqrt(denom)*denom**2)
           !
           !simple dipole test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !                  Br1 = 2.0*cos_theta1/r_1**3
           !              Btheta1 = sin_theta1/r_1**3
           !                Bphi1 = 0.0
           ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !magnetic field components in global spherical coordinates
           Br1 =  B1scl*(1.0/(sin_theta1*r_1**2))*dA1dth
           Btheta1 = -B1scl*(1.0/(sin_theta1*r_1   ))*dA1dr
           Bphi1 =  0.0
           !make the field dipolar
           if (y .lt. 0.0) then
              Br1 = -Br1
              Btheta1 = -Btheta1
              Bphi1 = -Bphi1
           endif
           dBr1dr = -1.0/(sin_theta1*r_1**2)*((2.0/r_1)*dA1dth - d2A1drdth)
           dBtheta1dr =  1.0/(sin_theta1*r_1   )*((1.0/r_1)*dA1dr  -  d2A1dr2)
           dBphi1dr =  0.0
           !PRESSURE
           pres_1 =  inv_g*rho1scl*SSPscl**2
           dpres_1dr1 =  0.0
           !---------------------------FLUX ROPE REGION--------------------------------
           !CALCULATE CELL CENTER Cartesian COORDINATES for CME FLUX ROPE
           !self-similar transformation r = r --> r -a
           x_2 = x_1
           z_2 = y_1
           y_2 = z_1 - cme_r1
           !CALCULATE CELL CENTER SPHERICAL COORDINATES for CME FLUX ROPE
           r_2 = sqrt(x_2**2 + y_2**2 + z_2**2)
           cos_theta2 = z_2/r_2
           sin_theta2 = sqrt(x_2**2 + y_2**2)/r_2
           cos_phi2   = x_2/sqrt(x_2**2 + y_2**2)
           sin_phi2   = y_2/sqrt(x_2**2 + y_2**2)
           if (abs(r_2) .lt. delta) then
              r_2 =  delta
              x_2 =  delta*sin_theta2*cos_phi2
              y_2 =  delta*sin_theta2*sin_phi2
              z_2 =  delta*cos_theta2
           endif

           alpha0 = 5.763854/cme_r0
           ga0r0 = sin(alpha0*cme_r0)/alpha0*cme_r0 - cos(alpha0*cme_r0)

           A2 = (4.0*cPi*a1scl/alpha0**2)*( (cme_r0**2/ga0r0)               * &
                (sin(alpha0*r_2)/(alpha0*r_2) - cos(alpha0*r_2)) - r_2**2) * &
                sin_theta2**2

           dA2dr = ((4.0*cPi*a1scl/alpha0**2)*( (cme_r0**2/ga0r0)              * & 
                (cos(alpha0*r_2)/r_2 - sin(alpha0*r_2)/(alpha0*r_2**2)     + & 
                alpha0*sin(alpha0*r_2)) -  2.0*r_2 ))* sin_theta2**2

           dA2dth = (8.0*cPi*a1scl/alpha0**2)*( (cme_r0**2/ga0r0)               * & 
                (sin(alpha0*r_2)/(alpha0*r_2) - cos(alpha0*r_2)) - r_2**2) * &
                sin_theta2*cos_theta2 

           d2A2dr2 = (4.0*cPi*a1scl/alpha0**2)*sin_theta2**2                     * &
                ( (cme_r0**2/ga0r0)*(2.0*sin(alpha0*r_2)/(alpha0*r_2**3)   - &
                2.0*cos(alpha0*r_2)/(r_2**2) -alpha0*sin(alpha0*r_2)/r_2 + &
                (alpha0**2)*cos(alpha0*r_2)) - 2.0 )  

           d2A2drdth = (8.0*cPi*a1scl/alpha0**2)*sin_theta2*cos_theta2             * &
                ( (cme_r0**2/ga0r0)*(cos(alpha0*r_2)/r_2                   - & 
                sin(alpha0*r_2)/(alpha0*r_2**2)                         + &
                alpha0*sin(alpha0*r_2)) -2.0*r_2 ) 

           d2A2dth2 = (8.0*cPi*a1scl/alpha0**2)*( (cme_r0**2/ga0r0)               * &
                (sin(alpha0*r_2)/(alpha0*r_2) - cos(alpha0*r_2)) - r_2**2) * &
                (cos_theta2**2 - sin_theta2**2)

           dr2dr1   = (r_1 -cme_r1*cos_theta1)                                   / &
                sqrt(r_1**2 -2.0*r_1*cme_r1*cos_theta1 + cme_r1**2)

           dth2dr1 = (cme_r1/(r_1**2))*(cos_theta1 - cme_r1/r_1)                / &
                ( (sin_theta1*sin_phi1/r_1)                                * &
                sqrt(r_1**2*sin_theta1**2*cos_phi1**2                      + &
                r_1**2*cos_theta1**2 -2.0*r_1*cme_r1*cos_theta1 + cme_r1**2) + &
                (r_1**2*sin_theta1**2*cos_phi1**2 +  r_1**2*cos_theta1**2  - &
                2.0*r_1*cme_r1*cos_theta1 + cme_r1**2)**(1.5)             / &
                (sin_theta1*sin_phi1*r_1**3) )

           !derivative of field components in flux rope spherical coordinates
           dBr2dr2 = -2.0*dA2dth/(sin_theta2*r_2**3) + d2A2drdth/(sin_theta2*r_2**2)

           dBr2dth2 = -cos_theta2*dA2dth/(r_2**2*sin_theta2**2)                  + & 
                d2A2dth2/(sin_theta2*r_2**2)

           dBth2dr2 =  dA2dr/(sin_theta2*r_2**2) - d2A2dr2/(r_2 *sin_theta2)

           dBth2dth2 =  cos_theta2*dA2dr/(r_2*sin_theta2**2) -d2A2drdth/(r_2*sin_theta2)

           dBphi2dr2 =  alpha0*dA2dr /(r_2*sin_theta2) - alpha0*A2/(sin_theta2*r_2**2)

           dBphi2dth2 =  alpha0*dA2dth/(r_2*sin_theta2)                            - & 
                alpha0*cos_theta2*A2/(r_2*sin_theta2**2)

           !magnetic field components in the flux rope spherical coordinates
           Br2 = (1.0/(sin_theta2*r_2**2))*dA2dth
           Btheta2 =-(1.0/(sin_theta2*r_2   ))*dA2dr
           Bphi2 = (1.0/(sin_theta2*r_2   ))*alpha0*A2
           !magnetic field components in global cartesian coordinates
           !X COMPONENT OF MAGNETIC FIELD
           Bx_1 = Br2    *sin_theta2*cos_phi2 + &
                Btheta2*cos_theta2*cos_phi2 - &
                Bphi2  *sin_phi2
           !z COMPONENT OF MAGNETIC FIELD
           Bz_1 = Br2    *sin_theta2*sin_phi2 + &
                Btheta2*cos_theta2*sin_phi2 + &
                Bphi2  *cos_phi2
           !Y COMPONENT OF MAGNETIC FIELD
           By_1 = Br2    *cos_theta2          - &
                Btheta2*sin_theta2
           !INSIDE THE MAGNETIC FLUX ROPE
           if( sqrt(x_1**2 + y_1**2 + (z_1-cme_r1)**2) .lt. cme_r0 ) then
              !magnetic field components in global sperical coordinates
              Br1 = Bx_1   *sin_theta1*cos_phi1 + &
                   By_1   *sin_theta1*sin_phi1 + &
                   Bz_1   *cos_theta1
              Btheta1 = Bx_1   *cos_theta1*cos_phi1 + &
                   By_1   *cos_theta1*sin_phi1 - &
                   Bz_1   *sin_theta1
              Bphi1 = Bx_1   *(-1.0)*sin_phi1     + &
                   By_1   *cos_phi1 
              dBr1dr = dBr2dr2  *dr2dr1 + dBr2dth2  *dth2dr1 
              dBtheta1dr = dBth2dr2 *dr2dr1 + dBth2dth2 *dth2dr1
              dBphi1dr = dBphi2dr2*dr2dr1 + dBphi2dth2*dth2dr1
              !PRESSURE 
              pres_1 = inv_g*rho1scl*SSPscl**2 + a1scl*A2
              dpres_1dr1 = a1scl*(dA2dr*dr2dr1 + dA2dth*dth2dr1)
           endif

           ! !!!!!!!!!!INCLUSIVE OF BOTH MAGNETIC FIELD REGIONS!!!!!!!!!!!!!!!!!!!!!!
           !MAGNETIC FIELD with self-similar transformation
           Br = Br1    *(lambda/r)**2
           Btheta = Btheta1*(lambda/r)
           Bphi = Bphi1  *(lambda/r)
           !magnetic field components in global cartesian coordinates
           !X COMPONENT OF MAGNETIC FIELD
           State_VGB(Bx_,i,j,k,globalBLK) = Br    *sin_theta1*cos_phi1 + &
                Btheta*cos_theta1*cos_phi1 - &
                Bphi  *sin_phi1
           !Y COMPONENT OF MAGNETIC FIELD
           State_VGB(By_,i,j,k,globalBLK) = Br    *sin_theta1*sin_phi1 + &
                Btheta*cos_theta1*sin_phi1 + &
                Bphi  *cos_phi1
           !Z COMPONENT OF MAGNETIC FIELD
           State_VGB(Bz_,i,j,k,globalBLK) = Br    *cos_theta1          - &
                Btheta*sin_theta1
           !PLASMA DENSITY with self-similar transformation
           !     F_grav = (Ggrav*Msun/(r**2) + cme_alpha*r) 
           F_grav = (abs(Gbody)/(r**2) + cme_alpha*r) 
           State_VGB(rho_,i,j,k,globalBLK) = (1.0/ F_grav)             * &
                ( ((lambda/r)**2)*( (lambda/r)**2 - 1.0)           * &   
                (dpres_1dr1 + (1.0/(4.0*cPi))                    * &
                (Br1*dBr1dr + Btheta1*dBtheta1dr + Bphi1* dBphi1dr))     + &
                2.0*lambda*cme_a*pres_1/r**3                    + &
                cme_a*lambda/(4.0*cPi*r**3)*(1.0 - 2.0*(lambda/r)**2)*Br1**2 +&
                ((lambda/r)**2)*((cme_a/r)**2 +2.0*cme_a/r)     * &
                (Btheta1**2 + Bphi1**2)/(4.0*cPi*lambda) )

           State_VGB(rho_,i,j,k,globalBLK) = &
                State_VGB(rho_,i,j,k,globalBLK)  + &
                rho2scl/r**3

           if( (cme_a .eq. 0.0) .and. (rho2scl .eq. 0.0)) then
              State_VGB(rho_,i,j,k,globalBLK) = rho1scl
           endif

           !PLASMA PRESSURE  with self-similar transformation
           State_VGB(P_,i,j,k,globalBLK) = pres_1*(lambda/r)**2        - &
                (1.0/(8.0*cPi))*((lambda/r)**2)*((lambda/r)**2 - 1.0)*Br1**2 

           State_VGB(P_,i,j,k,globalBLK) =  &
                State_VGB(P_,i,j,k,globalBLK)     + &
                abs(Gbody)*rho2scl/(4.0*r**4)

           !MAGNETIC PRESSURE AT BUBBLE INTERFACE
           if( sqrt(x_1**2 + y_1**2 + (z_1-cme_r1)**2) .lt. cme_r0 ) then
              cos_thmax = (cme_r1**2 + r_1**2 -cme_r0**2)/(2.0*cme_r1*r_1)

              sin_thmax = sqrt(1.0 - cos_thmax**2)

              dsin_thmaxdr = -(cos_thmax/sin_thmax)* &
                   (r_1**2 +cme_r0**2 -cme_r1**2)/(2.0*cme_r1*r_1**2)

              denom_r0 = ( (cme_r0**2 -cme_r1**2)**2 + (cme_r1*r_1)**2                + &
                   2.0*r_1*cme_r1*(cme_r0**2 -cme_r1**2)*cos_thmax )

              dA1dth_r0 = -sin_thmax +  &
                   r_1*cme_r1*sin_thmax/ &
                   (cme_r0*sqrt(r_1**2 + cme_r1**2 - 2.0*r_1*cme_r1*cos_thmax))   + &
                   !
              r_1*(cme_r0**2 - 2.0*cme_r1**2)*sin_thmax/ & 
                   (cme_r0*sqrt(denom_r0))                                        - &
                   !
              r_1*cme_r1*(cme_r0**2-cme_r1**2)*sin_thmax * &
                   ( cme_r1*(r_1**2 + cme_r1**2 - cme_r0**2)   + &
                   r_1*(cme_r0**2-2.0*cme_r1**2)*cos_thmax )/ &
                   (cme_r0*denom_r0*sqrt(denom_r0))

              dA1dr_r0 = (r_1-cme_r1*cos_thmax) / &
                   ( cme_r0*sqrt(r_1**2 + cme_r1**2 -2.0*r_1*cme_r1*cos_thmax) )  - &
                   !
              (2.0*r_1*cme_r1 + (cme_r0**2 -2.0*cme_r1**2)*cos_thmax) / &
                   (cme_r0*sqrt(denom_r0))                                        + &
                   !
              ( cme_r1*(r_1**2 +cme_r1**2 -cme_r0**2) + &
                   r_1*(cme_r0**2 -2.0*cme_r1**2)*cos_thmax ) * &
                   (r_1*cme_r1**2 + cme_r1*(cme_r0**2 -cme_r1**2)*cos_thmax) / &
                   (cme_r0*denom_r0*sqrt(denom_r0))

              d2A1dth_r0dr = cme_r1*sin_thmax / &
                   (cme_r0*sqrt(r_1**2 +cme_r1**2 -2.0*r_1*cme_r1*cos_thmax))     - &
                   !
              r_1*cme_r1*sin_thmax*(r_1 -cme_r1*cos_thmax) / &
                   (cme_r0*(r_1**2 +cme_r1**2 -2.0*r_1*cme_r1*cos_thmax)**(1.5))  + &
                   !
              (cme_r0**2 -2.0*cme_r1**2)*sin_thmax/(cme_r0*sqrt(denom_r0))   - &
                   !
              r_1*sin_thmax*(cme_r0**2 -2.0*cme_r1**2) * &
                   (r_1*cme_r1**2 +cme_r1*(cme_r0**2 -cme_r1**2)*cos_thmax) / &
                   (cme_r0*denom_r0*sqrt(denom_r0))                               - &
                   !
              sin_thmax*cme_r1*(cme_r0**2 -cme_r1**2) * &
                   ( (cme_r1*(r_1**2 +cme_r1**2 -cme_r0**2) + &
                   r_1*(cme_r0**2 -2.0*cme_r1**2)*cos_thmax) - &
                   r_1*(2.0*r_1*cme_r1 +(cme_r0**2 -2.0*cme_r1**2)*cos_thmax) )/ &
                   (cme_r0*denom_r0*sqrt(denom_r0))                               + &
                   !
                   3.0*sin_thmax*r_1*cme_r1*(cme_r0**2 -cme_r1**2) * &
                   (cme_r1*(r_1**2 +cme_r1**2 -cme_r0**2) + &
                   r_1*(cme_r0**2 -2.0*cme_r1**2)*cos_thmax)     * &
                   (r_1*cme_r1**2 + cme_r1*(cme_r0**2 -cme_r1**2)*cos_thmax) / &
                   (cme_r0*sqrt(denom_r0)*denom_r0**2)

              d2A1dr_r0dr = 1.0/(cme_r0*sqrt(r_1**2 +cme_r1**2 -2.0*r_1*cme_r1*cos_thmax))-&
                   !
              (r_1 - cme_r1*cos_thmax)**2 / &
                   (cme_r0*(r_1**2 +cme_r1**2 -2.0*r_1*cme_r1*cos_thmax)**(1.5))  - &
                   !
                   2.0*cme_r1/(cme_r0*sqrt(denom_r0))                             + &
                   !
              (2.0*r_1*cme_r1 + (cme_r0**2-2.0*cme_r1**2)*cos_thmax) * &
                   (r_1*cme_r1**2 + cme_r1*(cme_r0**2-cme_r1**2)*cos_thmax) / &
                   (cme_r0*denom_r0*sqrt(denom_r0))                               + &
                   !
              ( (2.0*r_1*cme_r1 + (cme_r0**2 -2.0*cme_r1**2)*cos_thmax) * &
                   (r_1*cme_r1**2 + cme_r1*(cme_r0**2 -cme_r1**2)*cos_thmax)+ &
                   r_1*(cme_r1*(r_1**2 +cme_r1**2 -cme_r0**2) + &
                   r_1*(cme_r0**2 -2.0*cme_r1**2)*cos_thmax) ) / &
                   (cme_r0*denom_r0*sqrt(denom_r0))                               - &
                   !
                   3.0*(cme_r1*(r_1**2 +cme_r1**2 -cme_r0**2) + &
                   r_1*(cme_r0**2 -2.0*cme_r1**2)*cos_thmax)* &
                   (r_1*cme_r1**2 + cme_r1*(cme_r0**2 -cme_r1**2)*cos_thmax)**2 / &
                   (cme_r0*sqrt(denom_r0)*denom_r0**2) 

              Br_r0 = ((lambda/r)**2)*(1.0/(r_1*sin_thmax))*(dA1dth_r0/r_1) 

              Btheta_r0 = -(lambda/r)*    (1.0/(r_1*sin_thmax))* dA1dr_r0

              dBr_r0dr = r_1/(sin_thmax*r**4)*dA1dth_r0* &
                   (2.0 - 4.0*r_1/r - (r_1/sin_thmax)*dsin_thmaxdr) + &
                   r_1**2/(sin_thmax*r**4)*d2A1dth_r0dr

              dBtheta_r0dr = 1.0/(sin_thmax*r**2)*dA1dr_r0 * &
                   (1.0 - 2.0*r_1/r - (r_1/sin_thmax)*dsin_thmaxdr) + &
                   r_1   /(sin_thmax*r**2)*d2A1dr_r0dr

              Psurface = (1.0/(8.0*cPi))*(Br_r0**2 + Btheta_r0**2)

              dPsurfdr = (1.0/(4.0*cPi))*(Br_r0*dBr_r0dr + Btheta_r0*dBtheta_r0dr)

              State_VGB(P_,i,j,k,globalBLK) =  &
                   State_VGB(P_,i,j,k,globalBLK) + Psurface

              State_VGB(rho_,i,j,k,globalBLK) =&
                   State_VGB(rho_,i,j,k,globalBLK) - &
                   ((r**2)/abs(Gbody))*dPsurfdr 
           endif

           !VELOCITY FIELD of SELF-SIMILAR EXPANSION (radial flow only)
           vr = (Vscl/Rscl)*r
           State_VGB(rhoUx_,i,j,k,globalBLK) = vr*sin_theta1*cos_phi1
           State_VGB(rhoUy_,i,j,k,globalBLK) = vr*sin_theta1*sin_phi1
           State_VGB(rhoUz_,i,j,k,globalBLK) = vr*cos_theta1

           !NEGATIVE PRESSURE and DENSIY TEST
           if(State_VGB(P_,i,j,k,globalBLK) .le. 0.0) then
              write(*,*) 'cme_init: negative pressure at', i,j,k,globalBLK
              call stop_mpi('ERROR in set_ICs_cme')
           end if
           if(State_VGB(rho_,i,j,k,globalBLK) .le. 0.0) then
              write(*,*) 'cme_init: negative density at', i,j,k,globalBLK
              call stop_mpi('ERROR in set_ICs_cme')
           end if
        end do
     end do
  end do

end subroutine set_ICs_cme
