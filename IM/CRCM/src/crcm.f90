
subroutine crcm_run(delta_t)
  use ModConst,       ONLY: cLightSpeed, cElectronCharge
  use ModCrcmInitialize
  use ModCrcm,        ONLY: f2, dt, Time, phot, Pressure_IC,FAC_C
  use ModCrcmPlanet,  ONLY: re_m, dipmom, Hiono, nspec, amu_I, &
                            dFactor_I,tFactor_I
  use ModFieldTrace,  ONLY: fieldpara, brad=>ro, ftv=>volume, xo,yo,rb,irm,&
                            ekev,iba,bo,pp,Have, sinA, vel, alscone, iw2
  use ModGmCrcm,      ONLY: Den_IC,Temp_IC,StateIntegral_IIV,AveP_,AveDens_, &
                            iLatMin,DoMultiFluidGMCoupling,AveDen_I,AveP_I
  use ModIeCrcm,      ONLY: pot
  use ModCrcmPlot,    ONLY: Crcm_plot, Crcm_plot_fls, DtOutput, DoSavePlot,&
                            DoSaveFlux
  use ModCrcmRestart, ONLY: IsRestart
  implicit none


  integer n,nstep,ib0(nt)
  real delta_t, FactorTotalDens
  real flux(nspec,np,nt,neng,npit)
  real achar(nspec,np,nt,nm,nk)
  real vl(nspec,0:np,nt,nm,nk),vp(nspec,np,nt,nm,nk),fb(nspec,nt,nm,nk),rc
  integer iLat, iLon, iSpecies
  logical, save :: IsFirstCall =.true.
  !----------------------------------------------------------------------------

  if (dt==0) then
     nstep = 0
     dt = 0.0
  else
     nstep=nint(delta_t/dt)
     dt=delta_t/nstep         ! new dt
  endif


  
  ! do field line integration and determine vel, ekev, momentum (pp), etc.
  rc=(re_m+Hiono*1000.)/re_m        ! ionosphere distance in RE`
  call fieldpara(Time,dt,cLightSpeed,cElectronCharge,rc,re_m,xlat,xmlt,phi,xk,&
                 dipmom)
  
  !set boundary density and temperature inside iba
  if (.not. DoMultiFluidGMCoupling) then
     ! When not Multifluid we get the total density from Rho_MHD as follows:
     ! Rho = (m1n1+m2n2+...) = n * sum(m_i*dFactor_i)
     ! where sum(dFactor_i)=1 (over ions) and n_i=dFactor_i*n 
     ! n_i = dFactor_i*Rho_MHD/(sum(m_i*dFactor_i))
     ! n_total = Rho_MHD/sum(m_i*dfactor_i)
     ! FactorTotalDens = sum(m_i*dfactor_i)
     FactorTotalDens = sum(dFactor_I(1:nspec-1)*amu_I(1:nspec-1))
     do iSpecies = 1, nspec
        do iLon=1,nt
           do iLat=1,irm(iLon) 
              if (iLat < iLatMin) then
                 !Inside MHD boundary set den and temp to value at boundary
                 Den_IC(iSpecies,iLat,iLon) = dFactor_I(iSpecies) * &
                      StateIntegral_IIV(iLatMin,iLon,AveDens_)/FactorTotalDens
                 Temp_IC(iSpecies,iLat,iLon) = tFactor_I(iSpecies) * &
                      StateIntegral_IIV(iLatMin,iLon,AveP_) * FactorTotalDens &
                      / StateIntegral_IIV(iLatMin,iLon,AveDens_) &
                      * 6.2415e18 !J-->eV
!                 Den_IC(iSpecies,iLat,iLon) = dFactor_I(iSpecies) * 1.0e6
!                 Temp_IC(iSpecies,iLat,iLon) = tFactor_I(iSpecies)* 5000.0
              else
                 !Outside MHD boundary set den and temp from MHD
                 Den_IC(iSpecies,iLat,iLon) = dFactor_I(iSpecies) * &
                      StateIntegral_IIV(iLat,iLon,AveDens_)/FactorTotalDens
                 Temp_IC(iSpecies,iLat,iLon) = tFactor_I(iSpecies) * &
                      StateIntegral_IIV(iLat,iLon,AveP_) * FactorTotalDens &
                      / StateIntegral_IIV(iLat,iLon,AveDens_) &
                      * 6.2415e18 !J-->eV  
              endif
           end do
        end do
     end do
  else
     !Multifluid Case
     !Set Ion density and temperature
     do iSpecies = 1, nspec-1
        do iLon=1,nt
           do iLat=1,irm(iLon) 
              if (iLat < iLatMin) then
                 !Inside MHD boundary set den and temp to value at boundary
                 Den_IC(iSpecies,iLat,iLon) = &
                      StateIntegral_IIV(iLatMin,iLon,AveDen_I(iSpecies))&
                      / amu_I(iSpecies)
                 Temp_IC(iSpecies,iLat,iLon) = &
                      StateIntegral_IIV(iLatMin,iLon,AveP_I(iSpecies))&
                      /(Den_IC(iSpecies,iLat,iLon)) &
                        * 6.2415e18 !J-->eV
!                 Den_IC(iSpecies,iLat,iLon) = dFactor_I(iSpecies) * 1.0e6
!                 Temp_IC(iSpecies,iLat,iLon) = tFactor_I(iSpecies)* 5000.0
              else
                 !Outside MHD boundary set den and temp from MHD
                 Den_IC(iSpecies,iLat,iLon) = &
                      StateIntegral_IIV(iLat,iLon,AveDen_I(iSpecies))&
                      / amu_I(iSpecies)
                 Temp_IC(iSpecies,iLat,iLon) = &
                      StateIntegral_IIV(iLat,iLon,AveP_I(iSpecies))&
                      /(Den_IC(iSpecies,iLat,iLon)) &
                         * 6.2415e18 !J-->eV  
              endif
           end do
        end do
     end do
     !Set Electron density and temperature
     do iLon=1,nt
        do iLat=1,irm(iLon) 
           ! Density set by quasineutrality
           Den_IC(nspec,iLat,iLon)  = sum(Den_IC(1:nspec-1,iLat,iLon))
           ! Temp is set by 1/7 of weighted sum of ion temperatures
           Temp_IC(nspec,iLat,iLon) = 0.128205 * sum( &
                Den_IC(1:nspec-1,iLat,iLon)*Temp_IC(1:nspec-1,iLat,iLon)) &
                / Den_IC(nspec,iLat,iLon)
        end do
     end do
    !call CON_STOP('CRCM not set to use multifluid')
  endif


  ! setup initial distribution
  if (IsFirstCall .and. .not.IsRestart) then
     !set initial state when no restarting
     call initial_f2(nspec,np,nt,iba,amu_I,vel,xjac,ib0)
     IsFirstCall=.false.
  elseif(IsFirstCall .and. IsRestart) then
     ib0=iba
     IsFirstCall=.false.
  endif

  ! calculate boundary flux (fb) at the CRCM outer boundary at the equator
  call boundaryIM(nspec,np,nt,nm,nk,iba,irm,amu_I,xjac,vel,fb)
  
  ! calculate the drift velocity
  call driftV(nspec,np,nt,nm,nk,irm,re_m,Hiono,dipmom,dphi,xlat, &
       dlat,ekev,pot,vl,vp) 
  
  ! calculate the depreciation factor, achar, due to charge exchange loss
  call ceparaIM(nspec,np,nt,nm,nk,irm,dt,vel,ekev,Have,achar)
  
  ! Calculate the strong diffusion lifetime for electrons
  call StDiTime(dt,vel,ftv,rc,re_m,dipmom,iba)

  ! time loop
  do n=1,nstep
     call driftIM(iw2,nspec,np,nt,nm,nk,iba,dt,dlat,dphi,brad,rb,vl,vp, &
          fb,f2,ib0)
     call charexchangeIM(np,nt,nm,nk,nspec,iba,achar,f2)
     call lossconeIM(np,nt,nm,nk,nspec,iba,alscone,f2)
     call StrongDiff(iba)                               
     Time = Time+dt
 enddo
  
  ! Calculate CRCM output: flux, fac, phot
  call crcm_output(np,nt,nm,nk,nspec,neng,npit,iba,ftv,f2,ekev, &
       sinA,energy,sinAo,delE,dmu,amu_I,xjac,pp,xmm, &
       dmm,dk,xlat,dphi,re_m,Hiono,flux,FAC_C,phot,Pressure_IC)

  if (DoSavePlot.and.&
       (floor((Time+1.0e-5)/DtOutput))/=floor((Time+1.0e-5-delta_t)/DtOutput))&
       then
     call Crcm_plot(np,nt,xo,yo,Pressure_IC,phot,Den_IC,bo,ftv,pot,FAC_C,Time,dt)
     if (DoSaveFlux) call Crcm_plot_fls(rc,flux,time)
  endif
end subroutine Crcm_run

!-----------------------------------------------------------------------------
subroutine crcm_init
  !---------------------------------------------------------------------------
  ! Routine does CRCM initialization: fill arrays
  !
  ! Input: np,nt,neng,npit,nspec,re_m,dipmom,Hiono
  ! Output: xlat,xmlt,energy,sinAo (through augments)
  !         xmm1,xk1,phi1,dlat1,dphi1,dmm1,dk1,delE1,dmu1,xjac,amu (through 
  !         common block cinitialization

  use ModPlanetConst, ONLY: Earth_,DipoleStrengthPlanet_I,rPlanet_I
  use ModConst,       ONLY: cElectronCharge
  use ModNumConst,    ONLY: cDegToRad,cRadToDeg,cPi
  use ModCrcmPlanet,  ONLY: re_m, dipmom, Hiono, amu_I
  use ModCrcmInitialize
  use ModCrcmRestart, ONLY: IsRestart, crcm_read_restart

  implicit none

  integer i,n,k
  
  real rw,rsi,rs1
  real xjac1,sqrtm

  ! Define constants
  re_m = rPlanet_I(Earth_)                            ! earth's radius (m)
  dipmom=abs(DipoleStrengthPlanet_I(Earth_)*re_m**3)  ! earth's dipole moment
  

  ! CRCM xlat and xmlt grids
  do i=1,np
     xlat(i)=xlat_data(i)
     dlat(i)=0.5*(xlat_data(i+1)-xlat_data(i-1))*cDegToRad    ! dlat in radian
  enddo
  xlatr=xlat*cDegToRad  
  dphi=2.*cPi/nt
  do i=1,nt
     phi(i)=(i-1)*dphi
     xmlt(i)=mod(phi(i)*12.0/cPi + 12.0,24.0)   
  enddo

  ! CRCM output grids: energy, sinAo, delE1, dmu1
  energy=(/1.0000,1.6795,2.8209,4.7378,7.9574,13.365, &
       22.447,37.701,63.320,106.35,178.62,300.00/)
  delE=0.5243*energy
  sinAo=(/0.010021,0.030708,0.062026,0.086108,0.16073,0.27682, &
       0.430830,0.601490,0.753790,0.863790,0.94890,0.98827/)
  dmu=(/0.000207365,0.000868320,0.00167125,0.00489855,0.0165792,0.0404637, &
       0.078819500,0.121098000,0.14729600,0.16555900,0.1738560,0.2486830/)
  do k=2,neng
     Ebound(k)=sqrt(energy(k-1)*energy(k))
  enddo
  Ebound(1) = energy(1)**2.0/Ebound(2)
  Ebound(neng+1)=energy(neng)**2.0/Ebound(neng)

  ! CRCM magnetic moment, xmm1
  xmm(1)=energy(1)*cElectronCharge/(dipmom/(2*re_m)**3.0)
  dmm(1)=xmm(1)*2.              
  rw=1.55                       
  do i=2,nm                    
     dmm(i)=dmm(1)*rw**(i-1)           
     xmm(i)=xmm(i-1)+0.5*(dmm(i-1)+dmm(i))
  enddo

  ! CRCM K, xk
  rsi=1.47
  xk(1)=40.*rsi
  rs1=(rsi-1.)/sqrt(rsi) ! in following sutup: xk(i+0.5)=sqrt(xk(i)*xk(i+1))
  do i=1,nk
     if (i.gt.1) xk(i)=xk(i-1)*rsi
     dk(i)=xk(i)*rs1                 
  enddo

  ! Calculate Jacobian, xjac
  do n=1,nspec 
     xjac1=4.*sqrt(2.)*cPi*(1.673e-27*amu_I(n))*dipmom/(re_m+Hiono*1000.)
     sqrtm=sqrt(1.673e-27*amu_I(n))
     do i=1,np
        do k=1,nm
           xjac(n,i,k)=xjac1*sin(2.*xlatr(i))*sqrt(xmm(k))*sqrtm
        enddo
     enddo
  enddo

  if(IsRestart) then
     !set initial state when restarting
     call crcm_read_restart
  endif

end subroutine crcm_init

!-------------------------------------------------------------------------------
subroutine initial_f2(nspec,np,nt,iba,amu_I,vel,xjac,ib0)
  !-----------------------------------------------------------------------------
  ! Routine setup initial distribution.
  ! 
  ! Input: nspec,np,nt,iba,Den_IC,Temp_IC,amu,vel,xjac
  ! Output: ib0,f2 (through common block cinitial_f2)
  Use ModGmCrcm, ONLY: Den_IC, Temp_IC
  use ModCrcm,   ONLY: f2
  use ModCrcmInitialize,   ONLY: IsEmptyInitial

  implicit none

  integer,parameter :: np1=51,nt1=48,nspec1=1  
  integer,parameter :: nm=35,nk=28 ! dimension of CRCM magnetic moment and K
  integer nspec,np,nt,iba(nt),ib0(nt),n,j,i,k,m
  real amu_I(nspec),vel(nspec,np,nt,nm,nk)
  real xjac(nspec,np,nm),pi,xmass,chmass,f21,vtchm

  pi=acos(-1.)

  ib0=iba
  f2=0.

  if (IsEmptyInitial) then
     ! Set initial f2 to a small number
     f2(:,:,:,:,:)=1.0e-40
  else
     ! Set initial f2 based on Maxwellian
     do n=1,nspec
        xmass=amu_I(n)*1.673e-27
        chmass=1.6e-19/xmass
        do j=1,nt
           do i=1,iba(j)
              f21=Den_IC(n,i,j)/(2.*pi*xmass*Temp_IC(n,i,j)*1.6e-19)**1.5
              do k=1,nm
                 do m=1,nk
                    vtchm=&
                         -vel(n,i,j,k,m)*vel(n,i,j,k,m)/2./Temp_IC(n,i,j)/chmass
                    f2(n,i,j,k,m)=xjac(n,i,k)*f21*exp(vtchm)
                 end do
              end do
           end do
        end do
     end do
  end if
end subroutine initial_f2


!-------------------------------------------------------------------------------
subroutine boundaryIM(nspec,np,nt,nm,nk,iba,irm,amu_I,xjac,vel,fb)
  !-----------------------------------------------------------------------------
  ! Routine setup the boundary distribution for the CRCM. Distribution at the
  ! boundary is assumed to be Maxwellian. Boundary temperature and density are
  ! from MHD.
  !
  ! Input: nspec,np,nt,nm,nk,iba,irm,amu,xjac,Den_IC,Temp_IC,vel
  ! Output: fb
  Use ModGmCrcm, ONLY: Den_IC, Temp_IC
  implicit none

  integer nspec,np,nt,nm,nk,iba(nt),irm(nt),j,n,k,m,ib1
  real amu_I(nspec),xjac(nspec,np,nm)
  real vel(nspec,np,nt,nm,nk),fb(nspec,nt,nm,nk),pi,xmass,chmass,fb1,vtchm

  pi=acos(-1.)

  do n=1,nspec
     xmass=amu_I(n)*1.673e-27
     chmass=1.6e-19/xmass
     do j=1,nt
        ib1=iba(j)+1
        if (ib1.gt.irm(j)) ib1=irm(j)
        fb1=Den_IC(n,ib1,j)/(2.*pi*xmass*Temp_IC(n,ib1,j)*1.6e-19)**1.5
        do k=1,nm
           do m=1,nk
              vtchm=-vel(n,ib1,j,k,m)*vel(n,ib1,j,k,m)/2./Temp_IC(n,ib1,j)/chmass
              fb(n,j,k,m)=xjac(n,ib1,k)*fb1*exp(vtchm)
           enddo
        enddo
     enddo
  enddo

end subroutine boundaryIM


!-------------------------------------------------------------------------------
subroutine ceparaIM(nspec,np,nt,nm,nk,irm,dt,vel,ekev,Have,achar)
  !-----------------------------------------------------------------------------
  ! Routine calculates the depreciation factor of H+, achar, due to charge
  ! exchange loss
  !
  ! Input: irm,nspec,np,nt,nm,nk,dt,vel,ekev,Have     ! Have: bounce-ave [H]
  ! Output: achar
  use ModCrcmPlanet,  ONLY: a0_I,a1_I,a2_I,a3_I,a4_I
  
  implicit none

  integer np,nt,nspec,nk,irm(nt),nm,i,j,k,m,n
  real vel(nspec,np,nt,nm,nk),ekev(np,nt,nm,nk),Have(np,nt,nk)
  real achar(nspec,np,nt,nm,nk),dt,Havedt,x,d,sigma,alpha

  do n=1,nspec-1
     do j=1,nt
        do i=1,irm(j)
           do m=1,nk
              Havedt=Have(i,j,m)*dt
              do k=1,nm
                 x=log10(ekev(i,j,k,m))
                 if (x.lt.-2.) x=-2.
                 d=a0_I(n)+a1_I(n)*x+a2_I(n)*x**2+a3_I(n)*x**3+a4_I(n)*x**4
                 sigma=10.**d        ! charge exchange cross section of H+ in m2
                 alpha=vel(n,i,j,k,m)*sigma*Havedt
                 achar(n,i,j,k,m)=exp(-alpha) ! charge. exchange decay rate
              enddo
           enddo
        enddo
     enddo
  enddo

end subroutine ceparaIM


!-------------------------------------------------------------------------------
subroutine driftV(nspec,np,nt,nm,nk,irm,re_m,Hiono,dipmom,dphi,xlat, &
     dlat,ekev,pot,vl,vp)
  !-----------------------------------------------------------------------------
  ! Routine calculates the drift velocities
  !
  ! Input: re_m,Hiono,dipmom,dphi,xlat,dlat,ekev,pot,nspec,np,nt,nm,nk,irm
  ! Output: vl,vp

  implicit none

  integer nspec,np,nt,nm,nk,irm(nt),n,i,ii,j,k,m,i0,i2,j0,j2,icharge
  real kfactor,xlat(np),xlatr(np),dlat(np),ekev(np,nt,nm,nk),pot(np,nt)
  real ksai,ksai1,xlat1,sf0,sf2,dlat2,re_m,Hiono,dipmom,dphi,pi,dphi2,cor
  real ham(np,nt),vl(nspec,0:np,nt,nm,nk),vp(nspec,np,nt,nm,nk)

  pi=acos(-1.)
  dphi2=dphi*2.
  kfactor=dipmom/(re_m+Hiono*1000.)
  cor=2.*pi/86400.                        ! corotation speed in rad/s
  xlatr=xlat*pi/180.

  nloop: do n=1,nspec
     if (n < nspec) then
        icharge=1
     else
        icharge=-1
     endif

     mloop: do m=1,nk
        kloop: do k=1,nm  

           ! ham: Hamiltonian/q
           ham(1:np,1:nt)=icharge*ekev(1:np,1:nt,k,m)*1000.+pot(1:np,1:nt)

           ! calculate drift velocities vl and vp
           iloop: do i=0,np
              ii=i
              if (i.eq.0) ii=1
              if (i.ge.1) ksai=kfactor*sin(2.*xlatr(i))
              if (i.lt.np) xlat1=0.5*(xlatr(ii)+xlatr(i+1))    ! xlat(i+0.5)
              ksai1=kfactor*sin(2.*xlat1)                   ! ksai at i+0.5
              jloop: do j=1,nt
                 j0=j-1
                 if (j0.lt.1) j0=j0+nt
                 j2=j+1
                 if (j2.gt.nt) j2=j2-nt

                 ! calculate vl
                 if (irm(j0).gt.i.and.irm(j2).gt.i) then
                    sf0=0.5*ham(ii,j0)+0.5*ham(i+1,j0)
                    sf2=0.5*ham(ii,j2)+0.5*ham(i+1,j2)
                    vl(n,i,j,k,m)=-(sf2-sf0)/dphi2/ksai1   ! vl at (i+0.5,j)
                 else
                    vl(n,i,j,k,m)=vl(n,i-1,j,k,m)
                 endif

                 ! calculate vp
                 if (i.ge.1) then
                    if (irm(j2).gt.i) then
                       i0=i-1
                       if (i.eq.1) i0=1
                       i2=i+1
                       if (i.eq.np) i2=np
                       dlat2=xlatr(i2)-xlatr(i0)
                       sf0=0.5*(ham(i0,j2)+ham(i0,j))
                       sf2=0.5*(ham(i2,j2)+ham(i2,j))
                       vp(n,i,j,k,m)=cor+(sf2-sf0)/dlat2/ksai  ! vp at (i,j+0.5)
                    else
                       vp(n,i,j,k,m)=vp(n,i-1,j,k,m)
                    endif
                 endif

              enddo jloop
           enddo iloop
        enddo kloop
     enddo mloop
  enddo nloop

end subroutine driftV


!-------------------------------------------------------------------------------
subroutine driftIM(iw2,nspec,np,nt,nm,nk,iba,dt,dlat,dphi,brad,rb,vl,vp, &
     fb,f2,ib0)
  !-----------------------------------------------------------------------------
  ! Routine updates f2 due to drift
  !
  ! Input: iw2,nspec,np,nt,nm,nk,iba,dt,dlat,dphi,brad,rb,vl,vp,fb 
  ! Input/Output: f2,ib0

  implicit none

  integer nk,iw2(nk),nspec,np,nt,nm,iba(nt),ib0(nt)
  integer n,i,j,k,m,j1,j_1,ibaj,ib,ibo,nrun,nn
  real dt,dlat(np),dphi,brad(np,nt),vl(nspec,0:np,nt,nm,nk),vp(nspec,np,nt,nm,nk)
  real rb,fb(nspec,nt,nm,nk),f2(nspec,np,nt,nm,nk)
  real f2d(np,nt),cmax,cl1,cp1,cmx,dt1,fb0(nt),fb1(nt),fo_log,fb_log,f_log
  real slope,cl(np,nt),cp(np,nt),fal(0:np,nt),fap(np,nt),fupl(0:np,nt),fupp(np,nt)
  logical :: UseUpwind=.false.
  nloop: do n=1,nspec
     mloop: do m=1,nk
        kloop: do k=1,iw2(m)
           f2d(1:np,1:nt)=f2(n,1:np,1:nt,k,m)         ! initial f2

           ! find nrun and new dt (dt1)
           cmax=0.
           do j=1,nt
              j1=j+1
              if (j1.gt.nt) j1=j1-nt
              ibaj=max(iba(j),iba(j1))
              do i=1,ibaj
                 cl1=dt/dlat(i)*vl(n,i,j,k,m)
                 cp1=dt/dphi*vp(n,i,j,k,m)
                 cmx=max(abs(cl1),abs(cp1))
                 cmax=max(cmx,cmax)
              enddo
           enddo
           nrun=ifix(cmax/0.50)+1     ! nrun to limit the Courant number
           dt1=dt/nrun                ! new dt
           ! Setup boundary fluxes and Courant numbers
           do j=1,nt
              ib=iba(j)
              ibo=ib0(j)
              fb0(j)=f2d(1,j)                   ! psd at inner boundary
              fb1(j)=fb(n,j,k,m)              ! psd at outer boundary
              if (ib.gt.ibo) then             ! during dipolarization
                 fo_log=-50.
                 if (f2d(ibo,j).gt.1.e-50) fo_log=log10(f2d(ibo,j))
                 fb_log=-50.
                 if (fb1(j).gt.1.e-50) fb_log=log10(fb1(j))
                 slope=(fo_log-fb_log)/(brad(ibo,j)-rb)
                 do i=ibo+1,ib
                    f_log=fo_log+slope*(brad(i,j)-brad(ibo,j))
                    f2d(i,j)=10.**f_log
                 enddo
              endif
              do i=1,np
                 cl(i,j)=dt1/dlat(i)*vl(n,i,j,k,m)
                 cp(i,j)=dt1/dphi*vp(n,i,j,k,m)
              enddo
           enddo

           ! run drift nrun times
           do nn=1,nrun
              UseUpwind=.false.
              call FLS_2D(np,nt,iba,fb0,fb1,cl,cp,f2d,fal,fap,fupl,fupp)
              fal(0,1:nt)=f2d(1,1:nt)
              jloop: do j=1,nt
                 j_1=j-1
                 if (j_1.lt.1) j_1=j_1+nt
                 iloop: do i=1,iba(j)
                    f2d(i,j)=f2d(i,j)+dt1/dlat(i)* &
                         (vl(n,i-1,j,k,m)*fal(i-1,j)-vl(n,i,j,k,m)*fal(i,j))+ &
                         cp(i,j_1)*fap(i,j_1)-cp(i,j)*fap(i,j)
                    if (f2d(i,j).lt.0.) then
                       if (f2d(i,j).gt.-1.e-30) then
                          f2d(i,j)=0.
                       else
                          write(*,*)'IM WARNING: f2d < 0 in drift ',n,i,j,k,m
                          write(*,*)'IM WARNING: Retrying step with upwind scheme'
                          UseUpwind=.true.
                          exit jloop
                       endif
                    endif
                 enddo iloop
              enddo jloop
              ! When regular scheme fails, try again with upwind scheme before 
              ! returning an error
              if (UseUpwind) then
                 fupl(0,1:nt)=f2d(1,1:nt)
                 do j=1,nt
                    j_1=j-1
                    if (j_1.lt.1) j_1=j_1+nt
                    do i=1,iba(j)
                       f2d(i,j)=f2d(i,j)+dt1/dlat(i)* &
                        (vl(n,i-1,j,k,m)*fupl(i-1,j)-vl(n,i,j,k,m)*fupl(i,j))+ &
                        cp(i,j_1)*fupp(i,j_1)-cp(i,j)*fupp(i,j)
                       if (f2d(i,j).lt.0.) then
                          if (f2d(i,j).gt.-1.e-30) then
                             f2d(i,j)=0.
                          else
                             write(*,*)'IM ERROR: f2d < 0 in drift ',n,i,j,k,m
                             call CON_STOP('CRCM dies in driftIM')
                          endif
                       endif
                    enddo
                 enddo
              endif
           enddo
           f2(n,1:np,1:nt,k,m)=f2d(1:np,1:nt)
        enddo kloop
     enddo mloop
  enddo nloop
  ! Update ib0
  ib0(1:nt)=iba(1:nt)

end subroutine driftIM


!-------------------------------------------------------------------------------
subroutine charexchangeIM(np,nt,nm,nk,nspec,iba,achar,f2)
  !-----------------------------------------------------------------------------
  ! Routine updates f2 due to charge exchange loss
  !
  ! Input: np,nt,nm,nk,nspec,achar   ! charge exchange depreciation of H+ 
  ! Input/Output: f2

  implicit none

  integer np,nt,nm,nk,nspec,iba(nt),n,i,j
  real achar(nspec,np,nt,nm,nk),f2(nspec,np,nt,nm,nk)

  do n=1,nspec-1             
     do j=1,nt
        do i=1,iba(j)
           f2(n,i,j,1:nm,1:nk)=f2(n,i,j,1:nm,1:nk)*achar(n,i,j,1:nm,1:nk)
        enddo
     enddo
  enddo
 
end subroutine charexchangeIM

!******************************************************************************
!                                StDiTime                                      
!  Routine calculate the strong diffusion lifetime for electrons.     
!*****************************************************************************
subroutine StDiTime(dt,vel,volume,rc,re_m,xme,iba)
  use ModCrcm,       ONLY: SDtime
  use ModCrcmGrid,   ONLY: np,nt,nm,nk, xlatr
  use ModCrcmPlanet, ONLY: nspec
  real vel(nspec,np,nt,nm,nk),volume(np,nt)
  integer iba(nt)
  
  
  eb=0.25                         ! fraction of back scatter e-   
  xmer3=xme/(rc*re_m)**3
  
  do j=1,nt
     do i=1,iba(j)
        !              xlat2=xlati(i)*xlati(i)!- from M.-Ch., Aug 1 2007  
        !              Bi=xmer3*sqrt(3.*xlat2+1.)     
        sinlat2=sin(xlatr(i))*sin(xlatr(i))
        Bi=xmer3*sqrt(3.*sinlat2+1.)      ! magnetic field at ionosphere 
        
        vBe=2.*volume(i,j)*Bi/(1.-eb)
        do k=1,nm
           do m=1,nk
              SDtime1=vBe/vel(nspec,i,j,k,m) !strong diff T,(gamma*mo/p = 1/v)
              SDtime(i,j,k,m)=exp(-dt/SDtime1)
           enddo
        enddo
     enddo
  enddo
  
  return
end subroutine StDiTime


!***********************************************************************   
!                            StrongDiff                                     
!  Routine calculate the change of electron psd (f2) by strong diffusion  
!***********************************************************************        
subroutine StrongDiff(iba)                               
  use ModCrcm,       ONLY: SDtime,f2
  use ModCrcmGrid,   ONLY: np,nt,nm,nk
  use ModCrcmPlanet, ONLY: nspec  
  implicit none
  integer iba(nt),i,j,k,m
  
  do j=1,nt
     do i=2,iba(j)
        do m=1,nk
           do k=1,nm
              f2(nspec,i,j,k,m)=f2(nspec,i,j,k,m)*SDtime(i,j,k,m)
           enddo
        enddo
     enddo
  enddo
  
  return
end subroutine StrongDiff



!-------------------------------------------------------------------------------
subroutine lossconeIM(np,nt,nm,nk,nspec,iba,alscone,f2)
  !-----------------------------------------------------------------------------
  ! Routine calculate the change of f2 due to lossconeIM loss
  ! 
  ! Input: np,nt,nm,nk,nspec,iba,alscone
  ! Input/Output: f2

  implicit none

  integer np,nt,nm,nk,nspec,iba(nt),n,i,j,k,m
  real alscone(nspec,np,nt,nm,nk),f2(nspec,np,nt,nm,nk)

  do n=1,nspec
     do j=1,nt
        do i=1,iba(j)
           do k=1,nm
              do m=1,nk
                 if (alscone(n,i,j,k,m).lt.1.) &
                      f2(n,i,j,k,m)=f2(n,i,j,k,m)*alscone(n,i,j,k,m)
              enddo
           enddo
        enddo
     enddo
  enddo

end subroutine lossconeIM


!-------------------------------------------------------------------------------
subroutine crcm_output(np,nt,nm,nk,nspec,neng,npit,iba,ftv,f2,ekev, &
     sinA,energy,sinAo,delE,dmu,amu_I,xjac,pp,xmm, &
     dmm,dk,xlat,dphi,re_m,Hiono,flux,fac,phot,Pressure_IC)
  !-----------------------------------------------------------------------------
  ! Routine calculates CRCM output, flux, fac and phot from f2
  !
  ! Input: np,nt,nm,nk,nspec,neng,npit,iba,ftv,f2,ekev,sinA,energy,sinAo,xjac
  !        delE,dmu,amu_I,xjac,pp,xmm,dmm,dk,xlat,dphi,re_m,Hiono
  ! Output: flux,fac,phot,Den_IC,Temp_IC
  Use ModGmCrcm, ONLY: Den_IC, Temp_IC
  use ModConst,   ONLY: cProtonMass
  use ModNumConst,ONLY: cPi, cDegToRad
  implicit none

  integer np,nt,nm,nk,nspec,neng,npit,iba(nt),i,j,k,m,n,j1,j_1
  real f2(nspec,np,nt,nm,nk),ekev(np,nt,nm,nk),sinA(np,nt,nk),re_m,Hiono,rion
  real ftv(np,nt),ftv1,energy(neng),sinAo(npit),delE(neng),dmu(npit),aloge(neng)
  real flux2D(nm,nk),pp(nspec,np,nt,nm,nk),xjac(nspec,np,nm)
  real sinA1D(nk),flx,ekev2D(nm,nk),flx_lo,pf(nspec),delEE(neng),pi
  real amu_I(nspec),amu1,psd1,psd(nspec,np,nt,nm,nk),fave(nspec,np,nt,neng)
  real xmm(nm),dmm(nm),dk(nk),xlat(np),xlatr(np),dphi,eta(nspec,np,nt,nm,nk)
  real flux(nspec,np,nt,neng,npit),detadi,detadj,dwkdi,dwkdj
  real fac(np,nt),phot(nspec,np,nt),Pressure_IC(nspec,np,nt)
  real Pressure1

  flux=0.
  fac=0.
  phot=0.

  ! Some constants for pressure, fac calculations
  rion=re_m+Hiono*1000.                      ! ionosphere distance in meter
  do n=1,nspec
     pf(n)=4.*cPi*1.e4/3.*sqrt(2.*cProtonMass*amu_I(n))*sqrt(1.6e-16)*1.e9  ! phot(nPa)
  enddo
  delEE=delE*sqrt(energy)
  xlatr=xlat*cDegToRad

  ! Calculate CRCM ion density (m^-3), Den_IC, and flux (cm^-2 s^-1 keV^-1 sr^-1)
  ! at fixed energy & pitch-angle grids 
  aloge=log10(energy)
  jloop1: do j=1,nt
     iloop1: do i=1,iba(j)
        ftv1=ftv(i,j)     ! ftv1: flux tube volume in m^3/Wb
        Pressure1=0.0
        nloop: do n=1,nspec
           Pressure1=0.0
           Den_IC(n,i,j)=0.0
           amu1=amu_I(n)**1.5
!!!! Calculate Den_IC, and 2D flux, fl2D(log), ekev2D(log) and sinA1D
           do m=1,nk
              sinA1D(m)=sinA(i,j,m)
              do k=1,nm
                 !write(*,*) 'n,i,k,xjac(n,i,k)',n,i,k,xjac(n,i,k)
                 psd1=f2(n,i,j,k,m)/1.e20/1.e19/xjac(n,i,k)  ! mug^-3cm^-6s^3
                 flx=psd1*(1.6e19*pp(n,i,j,k,m))*pp(n,i,j,k,m)
                 flux2D(k,m)=-50.
                 if (flx.gt.1.e-50) flux2D(k,m)=log10(flx)
                 ekev2D(k,m)=log10(ekev(i,j,k,m))
                 eta(n,i,j,k,m)=amu1*1.209*psd1*sqrt(xmm(k))*dmm(k)*dk(m)
                 psd(n,i,j,k,m)=psd1
                 Den_IC(n,i,j)=Den_IC(n,i,j)+eta(n,i,j,k,m)/ftv1
                 Pressure1=Pressure1+eta(n,i,j,k,m)*ekev(i,j,k,m)/ftv1
              enddo
           enddo
           Pressure1=Pressure1*1.6e-16*2./3.      ! pressure in Pa
           Pressure_IC(n,i,j)=Pressure1*1.e9           ! pressure in nPa
!!!! Map flux to fixed energy and pitch-angle grids (energy, sinAo)
           do k=1,neng
              do m=1,npit
                 call lintp2aIM(ekev2D,sinA1D,flux2D,nm,nk,aloge(k),sinAo(m),flx_lo)
                 flux(n,i,j,k,m)=10.**flx_lo
              enddo
           enddo
        enddo nloop
     enddo iloop1
  enddo jloop1

  ! Calculate pressure of the 'hot' ring current, phot, and temperature, Temp_IC
  jloop2: do j=1,nt
     iloop2: do i=1,iba(j)
!!!! calculate pitch-angle averaged flux
        do n=1,nspec
           do k=1,neng
              fave(n,i,j,k)=0.
              do m=1,npit
                 fave(n,i,j,k)=fave(n,i,j,k)+flux(n,i,j,k,m)*dmu(m)
              enddo
           enddo
        enddo
!!!! calculate pressure and temperature
        do n=1,nspec
           do k=1,neng
              phot(n,i,j)=phot(n,i,j)+fave(n,i,j,k)*delEE(k)*pf(n) ! phot in nPa
           enddo
           Temp_IC(n,i,j)=0.
           if (Den_IC(n,i,j).gt.0.) &
                Temp_IC(n,i,j)=phot(n,i,j)*1.e-9/Den_IC(n,i,j)/1.6e-19   ! eV
        enddo
     enddo iloop2
  enddo jloop2

  ! Calculate field aligned current, fac
  jloop3: do j=1,nt
     j1=j+1
     j_1=j-1
     if (j1.gt.nt) j1=j1-nt          !! periodic boundary condition
     if (j_1.lt.1) j_1=j_1+nt        !!
     iloop3: do i=2,iba(j)-1
        do k=1,nm
           do m=1,nk
              dwkdi=(ekev(i+1,j,k,m)-ekev(i-1,j,k,m))/(xlatr(i+1)-xlatr(i-1))
              dwkdj=(ekev(i,j1,k,m)-ekev(i,j_1,k,m))/(2.*dphi)
              do n=1,nspec
                 detadi=(eta(n,i+1,j,k,m)-eta(n,i-1,j,k,m))/(xlatr(i+1)-xlatr(i-1))
                 detadj=(eta(n,i,j1,k,m)-eta(n,i,j_1,k,m))/(2.*dphi)
                 fac(i,j)=fac(i,j)+(detadi*dwkdj-detadj*dwkdi)
              enddo
           enddo
        enddo
        fac(i,j)=1.6e-16*fac(i,j)/cos(xlatr(i))/rion**2    ! fac in Amp/m^2
     enddo iloop3
  enddo jloop3

end subroutine crcm_output


!-------------------------------------------------------------------------------
subroutine FLS_2D(np,nt,iba,fb0,fb1,cl,cp,f2d,fal,fap,fupl,fupp)
!-------------------------------------------------------------------------------
  !  Routine calculates the inter-flux, fal(i+0.5,j) and fap(i,j+0.5), using
  !  2nd order flux limited scheme with super-bee flux limiter method
  !
  !  Input: np,nt,iba,fb0,fb1,cl,cp,f2d
  !  Output: fal,fap
  use ModCrcm, ONLY: UseMcLimiter, BetaLimiter

  implicit none

  integer np,nt,iba(nt),i,j,j_1,j1,j2,ib
  real cl(np,nt),cp(np,nt),f2d(np,nt),fal(0:np,nt),fap(np,nt),fwbc(0:np+2,nt)
  real fb0(nt),fb1(nt),x,fup,flw,xsign,corr,xlimiter,r

  real,intent(out) :: fupl(0:np,nt), fupp(np,nt)
  fwbc(1:np,1:nt)=f2d(1:np,1:nt)        ! fwbc is f2d with boundary condition

  ! Set up boundary condition
  fwbc(0,1:nt)=fb0(1:nt)
  do j=1,nt
     ib=iba(j)
     fwbc(ib+1:np+2,j)=fb1(j)
  enddo

  ! find fal and fap
  jloop: do j=1,nt
     j_1=j-1
     j1=j+1
     j2=j+2
     if (j_1.lt.1) j_1=j_1+nt
     if (j1.gt.nt) j1=j1-nt
     if (j2.gt.nt) j2=j2-nt
     iloop: do i=1,np
        ! find fal
        xsign=sign(1.,cl(i,j))
        fupl(i,j)=0.5*(1.+xsign)*fwbc(i,j)+0.5*(1.-xsign)*fwbc(i+1,j) ! upwind
        flw=0.5*(1.+cl(i,j))*fwbc(i,j)+0.5*(1.-cl(i,j))*fwbc(i+1,j)   ! LW
        x=fwbc(i+1,j)-fwbc(i,j)
        if (abs(x).le.1.e-27) fal(i,j)=fupl(i,j)
        if (abs(x).gt.1.e-27) then
           if (xsign.eq.1.) r=(fwbc(i,j)-fwbc(i-1,j))/x
           if (xsign.eq.-1.) r=(fwbc(i+2,j)-fwbc(i+1,j))/x
           if (r.le.0.) fal(i,j)=fupl(i,j)
           if (r.gt.0.) then
              if(UseMcLimiter)then
                 xlimiter = min(BetaLimiter*r, BetaLimiter, 0.5*(1+r))
              else
                 xlimiter = max(min(2.*r,1.),min(r,2.))
              end if
              corr=flw-fupl(i,j)
              fal(i,j)=fupl(i,j)+xlimiter*corr
           endif
        endif
        ! find fap
        xsign=sign(1.,cp(i,j))
        fupp(i,j)=0.5*(1.+xsign)*fwbc(i,j)+0.5*(1.-xsign)*fwbc(i,j1) ! upwind
        flw=0.5*(1.+cp(i,j))*fwbc(i,j)+0.5*(1.-cp(i,j))*fwbc(i,j1)   ! LW
        x=fwbc(i,j1)-fwbc(i,j)
        if (abs(x).le.1.e-27) fap(i,j)=fupp(i,j)
        if (abs(x).gt.1.e-27) then
           if (xsign.eq.1.) r=(fwbc(i,j)-fwbc(i,j_1))/x
           if (xsign.eq.-1.) r=(fwbc(i,j2)-fwbc(i,j1))/x
           if (r.le.0.) fap(i,j)=fupp(i,j)
           if (r.gt.0.) then
              if(UseMcLimiter)then
                 xlimiter = min(BetaLimiter*r, BetaLimiter, 0.5*(1+r))
              else
                 xlimiter = max(min(2.*r,1.),min(r,2.))
              end if
              corr=flw-fupp(i,j)
              fap(i,j)=fupp(i,j)+xlimiter*corr
           endif
        endif
     enddo iloop
  enddo jloop

end subroutine FLS_2D

! OLD LINTP
!!-----------------------------------------------------------------------
!subroutine lintp(xx,yy,n,x,y,ier)
!  !-----------------------------------------------------------------------
!  !  Routine does 1-D interpolation.  xx must be increasing or decreasing
!  !  monotonically.  x is between xx(1) and xx(n)
!  ! 
!  !  input: xx,yy,n,x
!  !  output: y,ier
!
!  implicit none
!
!  integer n,ier,i,jl,ju,jm,j
!  real xx(n),yy(n),x,y,d
!
!  ier = 0
!
!  ! Make sure xx is increasing or decreasing monotonically
!  do i=2,n
!     if (xx(n).gt.xx(1).and.xx(i).lt.xx(i-1)) then
!        write(*,*) ' lintp: xx is not increasing monotonically '
!        write(*,*) n,xx
!        stop
!     endif
!     if (xx(n).lt.xx(1).and.xx(i).gt.xx(i-1)) then
!        write(*,*) ' lintp: xx is not decreasing monotonically '
!        write(*,*) n,xx
!        stop
!     endif
!  enddo
!
!  ! Set ier=1 if out of range
!  if (xx(n).gt.xx(1)) then
!     if (x.lt.xx(1).or.x.gt.xx(n)) ier=1
!  else
!     if (x.gt.xx(1).or.x.lt.xx(n)) ier=1
!  endif
!  if (ier.eq.1) then
!     write(*,*) ' Error: ier.eq.1'
!     print *,'n,x ',n,x
!     print *,'xx(1:n) ',xx(1:n)
!     stop
!  endif
!
!  ! initialize lower and upper values
!  jl=1
!  ju=n
!
!  ! if not done compute a midpoint
!10 if (ju-jl.gt.1) then
!     jm=(ju+jl)/2
!     ! now replace lower or upper limit
!     if ((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm))) then
!        jl=jm
!     else
!        ju=jm
!     endif
!     ! try again
!     go to 10
!  endif
!
!  ! this is the j
!  j=jl      ! if x.le.xx(1) then j=1
!  ! if x.gt.x(j).and.x.le.x(j+1) then j=j
!  ! if x.gt.x(n) then j=n-1
!  d=xx(j+1)-xx(j)
!  y=(yy(j)*(xx(j+1)-x)+yy(j+1)*(x-xx(j)))/d
!
!end subroutine lintp


!-------------------------------------------------------------------------------
subroutine lintp2aIM(x,y,v,nx,ny,x1,y1,v1)
  !-----------------------------------------------------------------------------
  !  This sub program takes 2-d interplation. x is 2-D and y is 1-D.
  !
  !  Input: x,y,v,nx,ny,x1,y1
  !  Output: v1

  implicit none               

  integer nx,ny,j,j1,i,i1,i2,i3
  real x(nx,ny),y(ny),v(nx,ny),x1,y1,v1,a,a1,b,x1d(1000)   ! max(nx)=1000
  real q00,q01,q10,q11

  call locate1IM(y,ny,y1,j)
  j1=j+1
  if (j.eq.0.or.j1.gt.ny) then
     b=1.
     if (j.eq.0) j=j1
     if (j1.gt.ny) j1=j
  else
     b=(y1-y(j))/(y(j+1)-y(j))
  endif

  x1d(1:nx)=x(1:nx,j)
  call locate1IM(x1d,nx,x1,i)
  i1=i+1
  if (i.eq.0.or.i1.gt.nx) then
     a=1.
     if (i.eq.0) i=i1
     if (i1.gt.nx) i1=i
  else
     a=(x1-x1d(i))/(x1d(i+1)-x1d(i))
  endif

  x1d(1:nx)=x(1:nx,j1)
  call locate1IM(x1d,nx,x1,i2)
  i3=i2+1
  if (i2.eq.0.or.i3.gt.nx) then
     a1=1.
     if (i2.eq.0) i2=i3
     if (i3.gt.nx) i3=i2
  else
     a1=(x1-x1d(i2))/(x1d(i2+1)-x1d(i2))
  endif

  q00=(1.-a)*(1.-b)
  q01=(1.-a1)*b
  q10=a*(1.-b)
  q11=a1*b
  v1=q00*v(i,j)+q01*v(i2,j1)+q10*v(i1,j)+q11*v(i3,j1)

end subroutine lintp2aIM


!--------------------------------------------------------------------------
subroutine locate1IM(xx,n,x,j)
  !--------------------------------------------------------------------------
  !  Routine return a value of j such that x is between xx(j) and xx(j+1).
  !  xx must be increasing or decreasing monotonically.
  !  If xx is increasing:
  !     If x=xx(m), j=m-1 so if x=xx(1), j=0  and if x=xx(n), j=n-1
  !     If x < xx(1), j=0  and if x > xx(n), j=n
  !  If xx is decreasing:
  !     If x=xx(m), j=m so if x=xx(1), j=1  and if x=xx(n), j=n
  !     If x > xx(1), j=0  and if x < xx(n), j=n
  !
  !  Input: xx,n,x
  !  Output: j

  implicit none

  integer n,j,i,jl,ju,jm
  real xx(n),x

  ! Make sure xx is increasing or decreasing monotonically
  do i=2,n
     if (xx(n).gt.xx(1).and.xx(i).lt.xx(i-1)) then
        write(*,*) ' locate1IM: xx is not increasing monotonically '
        write(*,*) n, (xx(j),j=1,n)
        call CON_STOP('CRCM stopped in locate1IM')
     endif
     if (xx(n).lt.xx(1).and.xx(i).gt.xx(i-1)) then
        write(*,*) ' locate1IM: xx is not decreasing monotonically '
        write(*,*) ' n, xx  ',n,xx
        call CON_STOP('CRCM stopped in locate1IM')
     endif
  enddo

  jl=0
  ju=n+1
  test: do
     if (ju-jl.le.1) exit test
     jm=(ju+jl)/2
     if ((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm))) then
        jl=jm
     else
        ju=jm
     endif
  end do test
  j=jl

end subroutine locate1IM


!Old CLOSED SUBROUTINE
!!--------------------------------------------------------------------------
!subroutine closed(n1,n2,yy,dx,ss)
!  !--------------------------------------------------------------------------
!  ! Routine does numerical integration using closed form.
!  ! 
!  ! Input: n1,n2,yy,dx
!  ! Output: ss
!
!  implicit none
!
!  integer n1,n2,i
!  real yy(n2),dx(n2),ss
!
!  ss=0.
!  do i=n1,n2
!     ss=ss+yy(i)*dx(i)
!  enddo
!
!end subroutine closed

