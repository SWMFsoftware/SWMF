
subroutine calc_electron_temperature(iBlock)

  use ModSizeGitm
  use ModGITM
  use ModPlanet
  use ModRates
  use ModEUV
  use ModSources, only: eEuvHeating,JouleHeating
  use ModConstants
  use ModTime, only: CurrentTime

  implicit none

  integer, external :: jday

  integer, intent(in) :: iBlock
  integer :: i,j,k,N,iLon,iLat,iAlt

  real :: NA,ddalt,h,kbc,omega,temp_A, temp_B

  real,dimension(nAlts) :: a,b,c,r,u,dz2

  real,dimension(nLons,nLats) :: flux

  integer, dimension(7) :: iTime

  integer :: nSteps, isteps
  real :: DtSub

  integer :: iError,DoY

  real, dimension(nLons, nLats, nAlts) :: &
       TOld, Conduction, Te_Advection,        &
       Heating, cooling, expansion,        &
       qD_N2,qD_O2,qD_O,qD_Total,f,g,hh,T0,T1,  &
       ZZ,Dx1,Dx2,Dx3,Ex1,Ex2,Ex3,ABT,d,Ki2,Ke2,   &
       lnA, Temp_T,   &
       Temp_vel,Grad_Vx,Grad_Vy,Grad_Vz,          &
       Grad_Qx,Grad_Qy,Grad_Qz,Temp_cooling, tmp,source_last,&
       IJouleHeating

!!!!!!!!!!!!!!!!!!!!!!!!!!

  call report("Calc_electron_temperature",1)

  !  eTemperature(:,:,:,iBlock)= Temperature(:,:,:,iBlock)*TempUnit * 2.0
  !  ITemperature(:,:,:,iBlock)= Temperature(:,:,:,iBlock)*TempUnit * 1.5
  !  return

  h=6.626e-34
  kbc=1.381e-23
  NA=6.0221e23

  call report("Calc_etemp_sources",5)
  call calc_etemp_sources(Heating,cooling,iBlock)

  call report("electron t coeff conductivity",5)
  call Coeff_Conductivity(iBlock)

  call report("electron t boundary conditions",5)
  call Boundary_Conditions(flux,iBlock)

  do iLon=1,nlons
     do iLat=1,nLats
        do iAlt = 1, nAlts
           dz2(iAlt) = dAlt_GB(iLon,iLat,iAlt,iBlock)**2

!!$              a(iAlt) = Ke(iLon,iLat,iAlt,iBlock)- &
!!$                   dKe(iLon,iLat,iAlt,iBlock)/2.                
!!$
!!$              c(iAlt) = Ke(iLon,iLat,iAlt,iBlock)+ &
!!$                   dKe(iLon,iLat,iAlt,iBlock)/2.                
!!$
!!$              b(iAlt) = -2*Ke(iLon,iLat,iAlt,iBlock)
!!$              r(iAlt) = (Cooling(iLon,iLat,iAlt)- &
!!$                         Heating(iLon,iLat,iAlt))*&
!!$                         dz2(iAlt)

           a(iAlt)=1.
           b(iAlt)=-2.
           c(iAlt)=1.
           r(iAlt)=(Cooling(iLon,iLat,iAlt)- &
                Heating(iLon,iLat,iAlt))*&
                dz2(iAlt)/Ke(iLon,iLat,iAlt,iBlock)/1.
           source_last(iLon,iLat,iAlt)=r(iAlt)

        enddo

        ! Boundary Conditions:

        a(1)=0
!!$           b(1)=1
!!$           c(1)=0
        r(1)=-Temperature(iLon,iLat,1,iBlock)*TempUnit(iLon,iLat,1)

        c(nAlts)=0

        !!! Why is the altitude difference face centered ???

        r(nAlts)=&!-0.00005* &
             -flux(ilon,ilat)/Ke(ilon,ilat,k,iBlock) &
             *dAlt_GB(iLon,iLat,nAlts,iBlock) &
             -eTemperature(ilon,ilat,nalts,iblock)

        call tridag(a,b,c,r,u)

        do iAlt=1,nAlts

           eTemperature(iLon,iLat,iAlt,iBlock)=(u(iAlt) + &
                eTemperature(iLon,iLat,iAlt,iBlock))/2.0

        enddo

     enddo
  enddo




  call report("Calc_etemp_sources (again)",5)
  call calc_etemp_sources(Heating,cooling,iBlock)

  call report("electron t coeff conductivity (again)",5)
  call Coeff_Conductivity(iBlock)

  call report("electron t boundary conditions (again)",5)
  call Boundary_Conditions(flux,iBlock)


  call report("electron t loop (again)",5)

  do iLon=1,nlons
     do iLat=1,nLats

        do iAlt = 1, nAlts
           r(iAlt)=(Cooling(iLon,iLat,iAlt)- &
                Heating(iLon,iLat,iAlt))*&
                dz2(iAlt)/Ke(iLon,iLat,iAlt,iBlock)/1.

           r(iAlt)=(r(iAlt)+source_last(iLon,iLat,iAlt))/2.

        enddo


        ! Boundary Conditions:

        a(1)=0
!!$           b(1)=1
!!$           c(1)=0
        r(1)=-Temperature(iLon,iLat,1,iBlock)*TempUnit(iLon,iLat,1)

        c(nAlts)=0

        !!! Why is this repeated ???
        r(nAlts)=&!-0.00005* &
             -flux(ilon,ilat)/Ke(iLon,iLat,k,iBlock) &
             *dAlt_GB(iLon,iLat,nAlts,iBlock) &
             -eTemperature(iLon,iLat,nAlts,iBlock)


        call tridag(a,b,c,r,u)



        do iAlt=1,nAlts

           eTemperature(iLon,iLat,iAlt,iBlock)=(u(iAlt) + &
                eTemperature(iLon,iLat,iAlt,iBlock))/2.0

           eTemperature(iLon,iLat,iAlt,iBlock)=  &
                max(eTemperature(iLon,iLat,iAlt,iBlock),&
                1.1*Temperature(iLon,iLat,iAlt,iBlock)*&
                TempUnit(iLon,iLat,iAlt))

           if (eTemperature(iLon,iLat,iAlt,iBlock)<10.0) then
              write(*,*) "eTemperature(i,j,k,iBlock)<10.0!!!"
              write(*,*)"eTemperature=",eTemperature(ilon,ilat,ialt,iBlock)
              write(*,*)"i=",ilon,"j=",ilat,"k=",ialt,&
                   'heating=',Heating(iLon,iLat,iAlt), &
                   'cooling=',cooling(iLon,iLat,iAlt)
              pause
           endif

        enddo




!!!!!!!!!!Calculate Ion temperature!!!!!!!!!!!!!!!!!!!!!!
        do iAlt=1,nAlts
           temp_A = 15 *3.2e-8*IDensityS(iLon,iLat,iAlt,ie_,iBlock)/ &
                (eTemperature(iLon,iLat,iAlt,iBlock))**1.5 *  &
                (IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock) + &
                0.5 *IDensityS(iLon,iLat,iAlt,iO2P_,iBlock) + & 
                0.53 *IDensityS(iLon,iLat,iAlt,iNOP_,iBlock))

           temp_B = (6.6e-14*NDensityS(iLon,iLat,iAlt,iN2_,iBlock) + &
                5.8e-14*NDensityS(iLon,iLat,iAlt,iO2_,iBlock) + &
                0.21e-14*NDensityS(iLon,iLat,iAlt,iO_3P_,iBlock)*  &
                sqrt(2*Temperature(iLon,iLat,iAlt,iBlock)*&
                TempUnit(iLon,iLat,iAlt))) *  &
                IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock)

           temp_B = temp_B + &
                (5.9e-14*NDensityS(iLon,iLat,iAlt,iN2_,iBlock) + &
                5.45e-14*NDensityS(iLon,iLat,iAlt,iO2_,iBlock) + &
                4.5e-14*NDensityS(iLon,iLat,iAlt,iO_3P_,iBlock)) *  &
                IDensityS(iLon,iLat,iAlt,iNOP_,iBlock)


           temp_B = temp_B + &
                (5.8e-14*NDensityS(iLon,iLat,iAlt,iN2_,iBlock) + &
                0.14e-14*NDensityS(iLon,iLat,iAlt,iO2_,iBlock)*  &
                sqrt(Temperature(iLon,iLat,iAlt,iBlock)* &
                TempUnit(iLon,iLat,iAlt)) + &
                4.4e-14*NDensityS(iLon,iLat,iAlt,iO_3P_,iBlock)) *  &
                IDensityS(iLon,iLat,iAlt,iO2P_,iBlock)


           IJouleHeating(iLon,iLat,iAlt)=3./2.* Rho(iLon,iLat,iAlt,iBlock) *  &
                (MeanMajorMass(iLon,iLat,iAlt)/AMU) /   &
                (MeanIonMass(iLon,iLat,iAlt)/AMU)*  &
                JouleHeating(iLon,iLat,iAlt) 

           ITemperature(iLon,iLat,iAlt,iBlock) = &
                (temp_A * eTemperature(iLon,iLat,iAlt,iBlock) +         &
                temp_B * Temperature(iLon,iLat,iAlt,iBlock)*&
                TempUnit(iLon,iLat,iAlt) + &
                IJouleHeating(iLon,iLat,iAlt))/ &
                (temp_A + temp_B)

        enddo

        !!           ITemperature(iLon,iLat,0,iBlock) = &
        !!                ITemperature(iLon,iLat,1,iBlock)
        !!           ITemperature(iLon,iLat,-1,iBlock) = &
        !!                ITemperature(iLon,iLat,1,iBlock)
        !
        !           ITemperature(iLon,iLat,nAlts+1,iBlock) = &
        !                ITemperature(iLon,iLat,nAlts,iBlock)
        !           ITemperature(iLon,iLat,nAlts+2,iBlock) = &
        !                ITemperature(iLon,iLat,nAlts,iBlock)
        !
        !!           eTemperature(iLon,iLat,0,iBlock) = &
        !!                eTemperature(iLon,iLat,1,iBlock)
        !!           eTemperature(iLon,iLat,-1,iBlock) = &
        !!                eTemperature(iLon,iLat,1,iBlock)
        !
        !           eTemperature(iLon,iLat,nAlts+1,iBlock) = &
        !                eTemperature(iLon,iLat,nAlts,iBlock)
        !           eTemperature(iLon,iLat,nAlts+2,iBlock) = &
        !                eTemperature(iLon,iLat,nAlts,iBlock)

     enddo
  enddo

contains


  !-------------------------------------------------
  !
  !-------------------------------------------------

  subroutine Boundary_Conditions(flux,iBlock)

    use ModSizeGitm
    implicit none

    real, intent(out)::flux(nlons,nlats)
    integer:: iBlock
    real :: Dst,Kp,bb1,kk1,P1,P2,bb2,kk2
    real :: x1,y1,z1,aa,bb,cc,temp,mlon,mlat
    real :: geo_lat,geo_lon

!!!! Lower Boundary Conditions!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    k=1
    eTemperature(:,:,k,iBlock)=Temperature(:,:,k,iBlock)*TempUnit(:,:,k)



!!!! Upper Boundary Conditions!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  k=nAlts  
    !  NewTemperature(:,:,k,Electron,iBlock)=NewTemperature(:,:,k-1,Electron,iBlock)&
    !  +0.00005*(Altitude_GB(:,:,k,iBlock)-Altitude_GB(:,:,k-1,iBlock))


!!!!!![Reference: 1. J.T. Hastings and R.G. Roble, Planet. Space Sci., Vol.25, pp.209, 1977.  Equation(31)
!!!!!!!           2. M.W. Liemohn et al, J. Geohys. Res., Vol.105, pp.27,767, 2000   Plate 1.            ]


    Kp=2.0
    Dst=-Kp*20.0

    call time_real_to_int(CurrentTime, iTime)
    DoY = jday(iTime(1), &
         iTime(2), &
         iTime(3)) 



    bb1=2365-508*exp(-((DoY-jday(1998,6,12))/180.)**2)
    kk1=7.8+6.2*sin(DoY*pi/365.)
    P1=1200.+ 500.*(1+cos(DoY*4*pi/365.))
    P2=4200.+ 650.*(1+cos(DoY*4*pi/365.))
    bb2=64.+3.*sin(DoY*3*pi/365.)
    kk2=1.6+0.5*sin(DoY*3*pi/365.)


    x1=16.0
    ! y1=bb2-kk2*Kp-5.0
    y1=50
    z1=7.0
    aa=10.0
    bb=20.0
    ! cc=2.5+(bb1-kk1*Dst-P1)/(P2-P1)
    cc=2.5

    k=nAlts  
    do j=1,nlats
       do i=1,nlons

          if (MLT(i,j,k) < 0.0) then
             MLT(i,j,k) = MLT(i,j,k) + 24.0
          endif

          if (abs(MLT(i,j,k)-x1)<12) then
             temp=1-(MLT(i,j,k)-x1)**2/aa**2- &
                  (abs(MLatitude(i,j,k,iBlock))-y1)**2/bb**2
          else
             temp=1-(24-abs(MLT(i,j,k)-x1))**2/aa**2- &
                  (abs(MLatitude(i,j,k,iBlock))-y1)**2/bb**2
          endif


          if (temp>0.0) then
             flux(i,j)=z1+cc*sqrt(temp)
          else
             flux(i,j)=z1
          endif

          flux(i,j)=exp(log(10.)*flux(i,j))


          !         flux(i,j)=1e9
          flux(i,j)=flux(i,j) * 1.602e-19 * 1e4 * 2.

          eTemperature(i,j,k,iBlock)= &
               eTemperature(i,j,k-1,iBlock)+&
               flux(i,j)/Ke(i,j,k,iBlock)*&
               (Altitude_GB(i,j,k,  iBlock) &
               -Altitude_GB(i,j,k-1,iBlock))

       enddo
    enddo


  end subroutine Boundary_Conditions





  !----------------------------------------------
  !
  !
  !----------------------------------------------


  Subroutine tridag(a,b,c,r,u)

    use ModSizeGitm
    implicit none

    real, intent(in)::a(nAlts),b(nAlts),c(nAlts),r(nAlts)
    real, intent(out)::u(nAlts)
    real :: bet,temp_u
    real :: gam(nAlts)

    !C-----------------------------------------------------------------------
    !      Subroutine tridag(a,b,c,r,u,n)
    !      Implicit Double Precision (A-z)
    !      Integer j, n, NMAX
    !      real :: a(n), b(n), c(n), r(n), u(n)
    !      real :: gam(n)
    !      Parameter (NMAX=500)
    !      !
    !      If (b(1).eq.0.D0) Pause 'tridag: rewrite equations'
    !      bet = b(1)
    !      u(1) = r(1)/bet
    !      Do 11 j = 2,n
    !         gam(j) = c(j-1)/bet
    !         bet = b(j)-a(j)*gam(j)
    !         If (bet.eq.0.D0) Pause 'tridag failed'
    !         u(j) = (r(j)-a(j)*u(j-1))/bet
    !   11 Continue
    !      Do 12 j = n-1,1,-1
    !         u(j) = u(j)-gam(j+1)*u(j+1)
    !   12 Continue
    !      Return
    !      End
    !C-----------------------------------------------------------------------


    ! This solver comes out of Numerical Recipies, page 48.


    bet = b(1)
    u(1) = r(1)/bet
    do iAlt = 2,nAlts
       gam(iAlt) = c(iAlt-1)/bet
       bet = b(iAlt)-a(iAlt)*gam(iAlt)
       If (bet.eq.0.0) then
          call stop_gitm("Error in tridiaginal solver in calc_cond")
       endif
       u(iAlt) = (r(iAlt)-a(iAlt)*u(iAlt-1))/bet
    enddo
    Do iAlt = nAlts-1,1,-1
       temp_u=u(iAlt)
       u(iAlt) = u(iAlt)-gam(iAlt+1)*u(iAlt+1)

    enddo

  end subroutine tridag



  !--------------------------------------------------------------------
  !
  !
  !--------------------------------------------------------------------

  subroutine Coeff_Conductivity(iBlock)

    use ModGITM
    implicit None

    integer::i,j,k,iBlock

    real, dimension(-1:nLons+2,                 &
         -1:nLats+2,                  &
         1:nAlts) ::  qD_N2,qD_O2,qD_O,qD_Total

    real, dimension(-1:nLons+2,-1:nLats+2,-1:nAlts+2, &
         nSpeciesTotal) :: Temp_NDensity

    real, dimension(-1:nLons+2,-1:nLats+2,-1:nAlts+2, &
         nIons) :: Temp_IDensity             

    Temp_NDensity(:,:,:,:)=NDensityS(:,:,:,:,iBlock)*1e-6
    Temp_IDensity(:,:,:,:)=IDensityS(:,:,:,:,iBlock)*1e-6



!!!!  Conduction !!!!!!!!!!!!!

    ! qD-N2,qD-O2,qD-O are the Maxwellian-averaged momentum transfer cross sections
    ! Ke is the electron thermal conductivity !


    do j=1,nLats
       do i=1,nLons
          do k = 1, nAlts

             qD_N2(i,j,k) = &
                  2.82e-17 * (1 - 1.21e-4 *                         &
                  eTemperature(i,j,k,iBlock)) *                         &
                  sqrt(eTemperature(i,j,k,iBlock))

             qD_O2(i,j,k)= &
                  2.2e-16 * (1 + 3.6e-2 *                            &
                  sqrt(eTemperature(i,j,k,iBlock)))

             qD_O(i,j,k)= &
                  1.1e-16 * (1 + 5.7e-4 * eTemperature(i,j,k,iBlock))  

             qD_Total(i,j,k)= &
                  !              Temp_NDensity(i,j,k,iN2_)*qD_N2(i,j,k)+            &
                  Temp_NDensity(i,j,k,iO2_)*qD_O2(i,j,k)+            &
                  Temp_NDensity(i,j,k,iO_3P_)*qD_O(i,j,k)

!!$         Ke(i,j,k,iBlock) = (7.7e5 * eTemperature(i,j,k,iBlock)**2.5) / &
!!$              (1 + 3.22e4 * eTemperature(i,j,k,iBlock)**2 *           &
!!$              qD_Total(i,j,k) / Temp_IDensity(i,j,k,ie_))

             Ke(i,j,k,iBlock) = 7.7e5 * eTemperature(i,j,k,iBlock)**2.5


             Ke(i,j,k,iBlock) = Ke(i,j,k,iBlock) * 1.602e-19 * 1e2 * 2.   


!!!!!! change the units eV cm-1 s-1 deg-1 to Joule m-1 s-1 deg-1 !!!!!!!!!1

          enddo
          ke(i,j,0,iBlock) = ke(i,j,1,iBlock)
          ke(i,j,nAlts+1,iBlock) = ke(i,j,nAlts,iBlock)

          do k=1,nAlts
!!$         dKe(i,j,k,iBlock)=  &
!!$                   (Ke(i,j,k+1,iBlock) - &
!!$                   Ke(i,j,k-1,iBlock))/2.
!!$
!!$              dke(i,j,k,iBlock) = 0.0

          enddo

       enddo
    enddo

  end subroutine Coeff_Conductivity
  !----------------------------------------------------
  !
  !
  !----------------------------------------------------

end subroutine calc_electron_temperature

!===============================================================================

subroutine calc_etemp_sources(Heating,Cooling,iBlock)

  use ModSizeGitm
  use ModGITM
  use ModPlanet
  use ModRates
  use ModEUV
  use ModSources, only: eEuvHeating
  use ModConstants

  implicit none

  integer, external :: jday

  integer, intent(in) :: iBlock
  real, intent(out)::Heating(nLons, nLats, nAlts)
  real, intent(out)::Cooling(nLons, nLats, nAlts)

  integer :: i,j,k,N

  real :: R,NA,ddalt,c,h,kbc,omega
  real :: c1,c2,c3,E1,E2,E3,AA1,AA2,AA3,B1,B2,B3 

  integer :: nSteps, isteps
  real :: DtSub

  integer :: iError,DoY

  real, dimension(nLons, nLats, nAlts) :: TOld,        &
       SolarHeating,eJouleHeating,expansion,        &
       f,g,hh,T0,T1,  &
       ZZ,Dx1,Dx2,Dx3,Ex1,Ex2,Ex3,ABT,d,Ki2,Ke2,   &
       lnA,Lrot_E_N2,Lrot_E_O2,Lvib_e_N2,         &
       Lvib_e_O2,Lf_e_O, L_e_O1D, L_e_i,Temp_T

  real, dimension(-1:nLons+2,-1:nLats+2,-1:nAlts+2, nSpeciesTotal) :: Temp_NDensity

  real, dimension(-1:nLons+2,-1:nLats+2,-1:nAlts+2, nIons) :: Temp_IDensity

!!!!!!!!!!!!!!!!!!!!!!!!!!

  !  call report("Heating Terms",1)

  c=2.998e8
  h=6.626e-34
  kbc=1.381e-23
  NA=6.0221e23
  R= 8.3145 

  Temp_NDensity(:,:,:,:)=NDensityS(:,:,:,:,iBlock)*1e-6
  Temp_IDensity(:,:,:,:)=IDensityS(:,:,:,:,iBlock)*1e-6


!!$ Ke(i,j,k,iBlock) = 7.7e5 * eTemperature(i,j,k,iBlock)**2.5 

!!!!!! change the units eV cm-1 s-1 deg-1 to Joule m-1 s-1 deg-1 !!!!!!!!!1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  do i=1,nLons
     do j=1,nLats
        do k = 1, nAlts

           !     write(*,*)'It is Ok----------1' 


           Lrot_E_N2(i,j,k)=2.9e-14 * Temp_IDensity(i,j,k,ie_) *       &
                Temp_NDensity(i,j,k,iN2_) *               &
                (eTemperature(i,j,k,iBlock) -        &
                (Temperature(i,j,k,iBlock)*TempUnit(i,j,k)))/         &
                sqrt(eTemperature(i,j,k,iBlock))

           Lrot_E_O2(i,j,k) =6.9e-14 * Temp_IDensity(i,j,k,ie_) *       &
                Temp_NDensity(i,j,k,iO2_) *               &
                (eTemperature(i,j,k,iBlock) -        &
                (Temperature(i,j,k,iBlock)*TempUnit(i,j,k)))/         &
                sqrt(eTemperature(i,j,k,iBlock))


           f(i,j,k) = 1.06e4 + 7.51e3 * tanh(1.10e-3 *                     &
                (eTemperature(i,j,k,iBlock)-1800.))

           g(i,j,k) = 3300. + 1.233 * (eTemperature(i,j,k,iBlock) - 1000.) -   &
                2.056e-4 * (eTemperature(i,j,k,iBlock) - 1000.) *     &
                (eTemperature(i,j,k,iBlock) - 4000.)

           hh(i,j,k) = 3300. - 839. * sin(1.91e-4 *                       &
                (eTemperature(i,j,k,iBlock)-2700.))


           Lvib_e_N2(i,j,k) = -2.99e-12 *Temp_IDensity(i,j,k,ie_) *            &
                Temp_NDensity(i,j,k,iN2_) *                    &
                exp(f(i,j,k)*(eTemperature(i,j,k,iBlock) - 2000)/   & 
                (2000.*eTemperature(i,j,k,iBlock)))*     &
                (exp(-g(i,j,k)*(eTemperature(i,j,k,iBlock)-          &
                Temperature(i,j,k,iBlock)*TempUnit(i,j,k))/ &
                (eTemperature(i,j,k,iBlock)*          &
                Temperature(i,j,k,iBlock)*TempUnit(i,j,k)))-1)


           !     write(*,*)'It is Ok----------1bb'

           Lvib_e_O2(i,j,k) = -5.196e-13 *Temp_IDensity(i,j,k,ie_) *           &
                Temp_NDensity(i,j,k,iO2_) *                    &
                exp(hh(i,j,k)*(eTemperature(i,j,k,iBlock) - 700)/    & 
                (700.*eTemperature(i,j,k,iBlock)))*      &
                (exp(-2770*(eTemperature(i,j,k,iBlock)-       &
                (Temperature(i,j,k,iBlock)*TempUnit(i,j,k)))/ &
                (eTemperature(i,j,k,iBlock)*          &
                (Temperature(i,j,k,iBlock)*TempUnit(i,j,k))))-1)                              
           !     write(*,*)'It is Ok----------1b'
           T0(i,j,k) =  Temperature(i,j,k,iBlock)*TempUnit(i,j,k)
           T1(i,j,k) =  Temperature(i,j,k,iBlock)*TempUnit(i,j,k)
           ZZ(i,j,k) = 5. + 3.* exp(-228./T1(i,j,k)) + exp(-326./T0(i,j,k))
           c1 = 0.02
           c2 = 0.028
           c3 = 0.008
           Dx1(i,j,k) = exp(-228./T1(i,j,k))
           Dx2(i,j,k) = exp(-326./T0(i,j,k))
           Dx3(i,j,k) = exp(-326./T0(i,j,k))
           Ex1(i,j,k) = exp(-228./eTemperature(i,j,k,iBlock))
           Ex2(i,j,k) = exp(-326./eTemperature(i,j,k,iBlock))
           Ex3(i,j,k) = exp(-98./eTemperature(i,j,k,iBlock)-228./T1(i,j,k))
           E1 = 228.
           E2 = 326.
           E3 = 98.
           AA1 = 8.58e-6
           AA2 = 7.201e-6
           AA3 = 2.463e-7
           B1 = 1.008
           B2 = 0.9617
           B3 = 1.1448


           ABT(i,j,k) =AA1*B1*eTemperature(i,j,k,iBlock)**(B1-0.5)*          &
                ( c1*(Dx1(i,j,k)-Ex1(i,j,k))+5.91e-9*(TempUnit(i,j,k)* &
                Temperature(i,j,k,iBlock)-     &
                eTemperature(i,j,k,iBlock))*((1+B1)*Dx1(i,j,k)+               &
                (E1/eTemperature(i,j,k,iBlock)+1+B1)*Ex1(i,j,k)))

           ABT(i,j,k) =ABT(i,j,k)+AA2*B2*eTemperature(i,j,k,iBlock)**(B2-0.5)*&
                ( c2*(Dx2(i,j,k)-Ex2(i,j,k))+5.91e-9*(TempUnit(i,j,k)* &
                Temperature(i,j,k,iBlock)-     &
                eTemperature(i,j,k,iBlock))*((1+B2)*Dx2(i,j,k)+               &
                (E2/eTemperature(i,j,k,iBlock)+1+B2)*Ex2(i,j,k)))

           ABT(i,j,k) =ABT(i,j,k)+AA3*B3*eTemperature(i,j,k,iBlock)**(B3-0.5)*&
                ( c3*(Dx3(i,j,k)-Ex3(i,j,k))+5.91e-9*(TempUnit(i,j,k)*Temperature(i,j,k,iBlock)-     &
                eTemperature(i,j,k,iBlock))*((1+B3)*Dx3(i,j,k)+               &
                (E3/eTemperature(i,j,k,iBlock)+1+B3)*Ex3(i,j,k)))


           Lf_e_O(i,j,k) = -8.629e-6 * Temp_IDensity(i,j,k,ie_) *           &
                Temp_NDensity(i,j,k,iO_3P_) *            &
                ABT(i,j,k)/ZZ(i,j,k)

           Lf_e_O(i,j,k) = 3.4e-12*(1-7e-5*eTemperature(i,j,k,iBlock))*    &
                Temp_IDensity(i,j,k,ie_) *           &
                Temp_NDensity(i,j,k,iO_3P_) *            &
                (150./eTemperature(i,j,k,iBlock)+0.4)* &
                (eTemperature(i,j,k,iBlock)- &
                Temperature(i,j,k,iBlock)*TempUnit(i,j,k))/ &
                (Temperature(i,j,k,iBlock)*TempUnit(i,j,k))


           !     write(*,*)'It is Ok----------1c'


           d(i,j,k) = 2.4e4 + 0.3 * (eTemperature(i,j,k,iBlock)-1500.) -      &
                1.947e-5 * (eTemperature(i,j,k,iBlock)-1500.) *         &
                (eTemperature(i,j,k,iBlock)-4000.)

           L_e_O1D(i,j,k) = -1.57e-12 *Temp_IDensity(i,j,k,ie_) *              &
                Temp_NDensity(i,j,k,iO_3P_) *                     &
                exp(d(i,j,k)*(eTemperature(i,j,k,iBlock) -3000)/    & 
                (3000.*eTemperature(i,j,k,iBlock)))*     &
                (exp(-22713*(eTemperature(i,j,k,iBlock)-      &
                TempUnit(i,j,k)*Temperature(i,j,k,iBlock))/     &
                (eTemperature(i,j,k,iBlock)*          &
                TempUnit(i,j,k)*Temperature(i,j,k,iBlock)))-1)
           !     write(*,*)'It is Ok----------1d'

           Ki2(i,j,k) = 4 * 3.14159 * Temp_IDensity(i,j,k,ie_) *              &
                (1.602e-19)**2/(kbc*ITemperature(i,j,k,iBlock))
           Ke2(i,j,k) = 4 * 3.14159 * Temp_IDensity(i,j,k,ie_) *              &
                (1.602e-19)**2/(kbc*eTemperature(i,j,k,iBlock))

           lnA(i,j,k) = LOG(4. * kbc * eTemperature(i,j,k,iBlock) /                &
                ((1.602e-19)**2 * sqrt(Ke2(i,j,k)))) - 2 * 0.577 -     &
                (Ki2(i,j,k) + Ke2(i,j,k)) / Ki2(i,j,k) *               &
                LOG(sqrt((Ki2(i,j,k) + Ke2(i,j,k)) / Ke2(i,j,k)))

           lnA(i,j,k) = 15

           L_e_i(i,j,k) = 3.2e-8/3. * Temp_IDensity(i,j,k,ie_) *                 &
                (eTemperature(i,j,k,iBlock) -                  &
                iTemperature(i,j,k,iBlock)) /                      &
                (eTemperature(i,j,k,iBlock))**1.5*            &
                lnA(i,j,k) * (Temp_IDensity(i,j,k,iO_4SP_)+         &
                !                          4.* Temp_IDensity(i,j,k,iHeP_)+                    &
                !                          16.*Temp_IDensity(i,j,k,iHP_) +                    &
                0.5 * Temp_IDensity(i,j,k,iO2P_) +                 &
                0.53 * Temp_IDensity(i,j,k,iNOP_))
           !     write(*,*)'It is Ok----------1e'

           TOld(i,j,k) = eTemperature(i,j,k,iBlock)

           !     write(*,*)'It is Ok----------2'           


!!!!  Heating !!!!!!!!!!!!!!!!!

           SolarHeating(i,j,k)= eEuvHeating(i,j,k,iBlock)       
           eJouleHeating(i,j,k)=0.0
           heating(i,j,k) = SolarHeating(i,j,k) + &
                eJouleHeating(i,j,k)
           heating(i,j,k) = heating(i,j,k) * 0.5
           !     write(*,*)'It is Ok----------3'

!!!!  cooling !!!!!!!!!!!!!!!!


           cooling(i,j,k) = Lvib_E_O2(i,j,k)             &
                +Lvib_E_N2(i,j,k)            &
                + L_e_O1D(i,j,k)             &
                !                           +  L_e_i(i,j,k)              &
                +  Lrot_e_N2(i,j,k)          &
                + Lrot_e_O2(i,j,k)!           &
           !                           + Lf_e_O(i,j,k) 


           cooling(i,j,k) = cooling(i,j,k) * 1.602e-19 * 1e6 * 0.5
!!!!!! change the units eV cm-3 s-1  to Joule m-3 s-1 !!!!!!!!!

           if (Altitude_GB(i,j,k,iBlock)/1000. <250.) &
                cooling(i,j,k) = 0.99*Heating(i,j,k)

        enddo


        do k = 1, nAlts           

           if (cooling(i,j,k)> 1.5*Heating(i,j,k)) &
                cooling(i,j,k)= 1.5*Heating(i,j,k)
           if (cooling(i,j,k)< 0.6*Heating(i,j,k)) &
                cooling(i,j,k)= 0.6*Heating(i,j,k)


        enddo

     enddo
  enddo

end subroutine calc_etemp_sources
