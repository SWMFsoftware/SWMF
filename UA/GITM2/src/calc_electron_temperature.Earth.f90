
subroutine calc_electron_temperature(iBlock)

  use ModSize
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
  integer :: i,j,k,N, CLAWIter,iLon,iLat,iAlt

  real, dimension(1:2) :: msis_temp
  real, dimension(1:8) :: msis_dens

  real :: NA,ddalt,h,kbc,omega,temp_A, temp_B

  real,dimension(nAlts) :: a,b,c,r,u,dz2

  real,dimension(nLons,nLats) ::flux

  integer, dimension(7) :: iTime

  integer :: nSteps, isteps
  real :: DtSub

  integer :: iError,DoY

  real, dimension(nLons,                 &
                  nLats,                  &
                  1:nAlts) :: TOld,Conduction,Te_Advection,        &
                   Heating, cooling,expansion,        &
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
 
call calc_etemp_sources(Heating,cooling,iBlock)

call Coeff_Conductivity(iBlock)

call Boundary_Conditions(flux,iBlock)




  do iLon=1,nlons
     do iLat=1,nLats

           do iAlt = 1, nAlts

          

              dz2(iAlt) = ((Altitude(iAlt+1) - Altitude(iAlt-1))/2.)**2


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
           r(1)=-Temperature(iLon,iLat,1,iBlock)*TempUnit



           c(nAlts)=0

           r(nAlts)=&!-0.00005* &
	-flux(ilon,ilat)/Ke(ilon,ilat,k,iBlock)* &
                    (altitude(nAlts)-altitude(nAlts-1))-     &
                    eTemperature(ilon,ilat,nalts,iblock)


           call tridag(a,b,c,r,u)

           do iAlt=1,nAlts

           eTemperature(iLon,iLat,iAlt,iBlock)=(u(iAlt) + &
                eTemperature(iLon,iLat,iAlt,iBlock))/2.0

           enddo

     enddo
  enddo




call calc_etemp_sources(Heating,cooling,iBlock)

call Coeff_Conductivity(iBlock)

call Boundary_Conditions(flux,iBlock)



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
           r(1)=-Temperature(iLon,iLat,1,iBlock)*TempUnit



           c(nAlts)=0
           r(nAlts)=&!-0.00005* &
		-flux(ilon,ilat)/Ke(ilon,ilat,k,iBlock)* &
                     (altitude(nAlts)-altitude(nAlts-1))-    &
                    eTemperature(ilon,ilat,nalts,iblock)
  
       
           call tridag(a,b,c,r,u)

 

              do iAlt=1,nAlts

                eTemperature(iLon,iLat,iAlt,iBlock)=(u(iAlt) + &
                eTemperature(iLon,iLat,iAlt,iBlock))/2.0

                eTemperature(iLon,iLat,iAlt,iBlock)=  &
                     max(eTemperature(iLon,iLat,iAlt,iBlock),&
                         1.1*Temperature(iLon,iLat,iAlt,iBlock)*TempUnit)

                if (eTemperature(iLon,iLat,iAlt,iBlock)<100.0) then
                   write(*,*) "Temperature(i,j,k,iBlock)<100.0!!!"
                   write(*,*)"Temperature=",eTemperature(ilon,ilat,ialt,iBlock)
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
                    0.21e-14*NDensityS(iLon,iLat,iAlt,iO_,iBlock)*  &
                    sqrt(2*Temperature(iLon,iLat,iAlt,iBlock)*TempUnit)) *  &
                    IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock)

           temp_B = temp_B + &
                    (5.9e-14*NDensityS(iLon,iLat,iAlt,iN2_,iBlock) + &
                    5.45e-14*NDensityS(iLon,iLat,iAlt,iO2_,iBlock) + &
                    4.5e-14*NDensityS(iLon,iLat,iAlt,iO_,iBlock)) *  &
                    IDensityS(iLon,iLat,iAlt,iNOP_,iBlock)

           
           temp_B = temp_B + &
                    (5.8e-14*NDensityS(iLon,iLat,iAlt,iN2_,iBlock) + &
                    0.14e-14*NDensityS(iLon,iLat,iAlt,iO2_,iBlock)*  &
                      sqrt(Temperature(iLon,iLat,iAlt,iBlock)*TempUnit) + &
                    4.4e-14*NDensityS(iLon,iLat,iAlt,iO_,iBlock)) *  &
                    IDensityS(iLon,iLat,iAlt,iO2P_,iBlock)


           IJouleHeating(iLon,iLat,iAlt)=3./2.* Rho(iLon,iLat,iAlt,iBlock) *  &
                    (MeanMajorMass(iLon,iLat,iAlt)/AMU) /   &
                    (MeanIonMass(iLon,iLat,iAlt)/AMU)*  &
                    JouleHeating(iLon,iLat,iAlt) 

           ITemperature(iLon,iLat,iAlt,iBlock) = &
               (temp_A * eTemperature(iLon,iLat,iAlt,iBlock) +         &
                temp_B * Temperature(iLon,iLat,iAlt,iBlock)*TempUnit + &
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
   
      use ModSize
      implicit none

      real, intent(out)::flux(nlons,nlats)
      integer:: iBlock
      real :: Dst,Kp,bb1,kk1,P1,P2,bb2,kk2
      real :: x1,y1,z1,aa,bb,cc,temp,mlon,mlat
      real :: geo_lat,geo_lon

!!!! Lower Boundary Conditions!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  k=1
  eTemperature(:,:,k,iBlock)=Temperature(:,:,k,iBlock)*TempUnit



!!!! Upper Boundary Conditions!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  k=nAlts  
!  NewTemperature(:,:,k,Electron,iBlock)=NewTemperature(:,:,k-1,Electron,iBlock)&
!  +0.00005*(altitude(:,:,k,iBlock)-altitude(:,:,k-1,iBlock))


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
              flux(i,j)/Ke(i,j,k,iBlock)*(altitude(k)-altitude(k-1))

     enddo
  enddo


end subroutine Boundary_Conditions


 


!----------------------------------------------
!
!
!----------------------------------------------


    Subroutine tridag(a,b,c,r,u)
      
      use ModSize
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
              Temp_NDensity(i,j,k,iO_)*qD_O(i,j,k)

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







