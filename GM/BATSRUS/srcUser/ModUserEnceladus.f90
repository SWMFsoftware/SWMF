!^CFG COPYRIGHT UM
!this file contains the ModUser foe a 1 species model of inosphereic mgnetospheric interaction in Enceladus, it contains Hptoionization of H2O, recombination, charge exchange and electron impact. 
!==============================================================================

module ModUser
  ! This is the user module for Enceladus
  use ModSize
  use ModUserEmpty,                      &
       IMPLEMENTED1 => user_read_inputs, &
       IMPLEMENTED2 => user_calc_sources

  include 'user_module.h' !list of public methods

  ! This is a temporary solution
  real, public:: BoundaryCoeffB=2.0, BoundaryCoeffU=2.0

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real, parameter :: VersionUserModule = 1.1
  integer,parameter :: MaxNuSpecies=1
  integer::iBlock
  real::NumDenNeutral_VC(nI,nJ,nK)

  real,dimension(0:26)::Impact,Te
  character (len=*), parameter :: NameUserModule = &
       'Enceladus 1 species MHD code, Dalal Najib'

contains

  subroutine user_read_inputs

    use ModReadParam
    character (len=100) :: NameCommand

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       
       case('#BOUNDARYCOEFF')
          call read_var('BoundaryCoeffB',BoundaryCoeffB)
          call read_var('BoundaryCoeffU',BoundaryCoeffU)
          !write(*,*)'BoundaryCoeffB=',BoundaryCoeffB
          !write(*,*)'BoundaryCoeffU=',BoundaryCoeffU
       case('#USERINPUTEND')
          EXIT
       case default
          call stop_mpi('Unknown user command:'//trim(NameCommand))
       end select
    end do

  end subroutine user_read_inputs

  !===========================================================================

  subroutine user_calc_sources(iBlock)

    use ModVarIndexes
    use BATL_size,  ONLY: nI, nJ, nK, nIJK
    use ModAdvance, ONLY: State_VGB, Source_VC, Energy_
    use ModProcMH,   ONLY: iProc
    use ModMain, ONLY: ProcTest, BlkTest, iTest, jTest, kTest 
    use ModPhysics, ONLY: Rbody, inv_gm1, gm1, &
         No2Si_V, Si2No_V, No2Io_V, Io2No_V,UnitN_, UnitT_,UnitTemperature_,&
         UnitX_, UnitRhoU_, UnitU_
    use ModProcMH,   ONLY: iProc
    use ModGeometry, ONLY: Xyz_DGB, R_BLK
    use ModBlockData,ONLY: use_block_data, put_block_data, get_block_data, &
         MaxBlockData

    integer, intent(in) :: iBlock

    character (len=*), parameter :: Name='user_calc_sources'

    integer::i,j,k,h,iBlockLast=-1
    real:: R0, R0_Si

    ! to convert from eV to Kelvins 1eV=1.1604e+4 K
    real,parameter::kTn_dim = 180 !in K
    ! density of hot electron that cause imapact ionization
    real, parameter::ne_hot=0.2 
    ! electron impact rate from Burger only for electron of Te=12.5eV
    real::kappa 
    !real,parameter:: PhotoIonRate_H2O=9.1e-9 !Solarmax condition
    !real,parameter::PhotoIonRate_H2O=3.6e-9 !Solarmin condition
    real,parameter::PhotoIonRate_H2O=3.6e-7 !Solarmin condition test

    ! cross section for charge exchange in cm2
    real, parameter:: sigma_exchange=8.1e-16
    real::kTi, kTi_dim, kTn,rhotimeU,rhotimeU_Si ,nion 
    real::inv_rho, inv_rho2, uu2, U_Si, U_No, Rho_No
    real::totalNumRho, totalSourceNumRho
    real:: RhoDot,RhoNumDot, RhoDotPhoto,RhoDotImpact,RhoDotRecomb, &
         RhoDotchargeX, RhoDotL, RhoNumDotL, RhoDotLx, RhoNumDotLx
    real:: ReactionRate_H2O_dim, ReactionRate_H2O,ImpactRate_dim, &
         RecombRate_dim, RecombRate 
    real:: Nu_C,Nu_C_dim, Nu_exchange!Collision frequency
    real,parameter::nu0=5.4745e-10!this is in cm^3/s

    !------------------------------------------------------------------
    !write(*,*)'iBlock',iBlock

!!$
!!$    write(*,*)'before update'
!!$    write(*,*)' Source_VC(Rho)=',Source_VC(rho_,iTest,jTest,kTest)
!!$    write(*,*)' Source_VC(Ux,Uy,Uz)=',Source_VC(RhoUx_:RhoUz_,iTest,jTest,kTest)
!!$    write(*,*)' Source_VC(Energy)=',Source_VC(Energy_,iTest,jTest,kTest)
!!$    write(*,*)' Source_VC(P)=',Source_VC(p_,iTest,jTest,kTest)

    ! Max number of reals saved in ModBlockData (NumDensNeutral is not saved!)
    MaxBlockData = nIJK

    !put the interpolated data in ModBlockData     
    if(iBlock /= iBlockLast)then       
       iBlockLast = iBlock       
       if(use_block_data(iBlock))then          
          !call get_block_data(iBlock, MaxNuSpecies, nI, nJ, nK, NumDenNeutral_VC)         
          call get_block_data(iBlock, nI, nJ, nK, NumDenNeutral_VC)

       else         
          call enceladus_input(iBlock)
          !write(*,*) 'enceladus_input has been read'
          !write(*,*)'Te(22)=',Te(22)

          !call put_block_data(iBlock, MaxNuSpecies, nI, nJ, nK, NumDenNeutral_VC)         
          call put_block_data(iBlock, nI, nJ, nK, NumDenNeutral_VC)

       end if
    end if

    !End of interpolation of Data

    do k = 1, nK ;   do j = 1, nJ ;  do i = 1, nI  
       R0=R_BLK(i,j,k,iBlock)
       R0_Si=R_BLK(i,j,k,iBlock)*No2Si_V(UnitX_)
       ! write(*,*)'R_block_Si=',R0_Si
       !write(*,*)'Density_Valiery=',NumDenNeutral_VC(i,j,k)
       !write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       totalNumRho=State_VGB(rho_,i,j,k,iBlock)/(MassFluid_I(1))

       ReactionRate_H2O_dim=PhotoIonRate_H2O*NumDenNeutral_VC(i,j,k)!neutral_density(R0)

       ! Now normalize the Reaction Rate
       ReactionRate_H2O= ReactionRate_H2O_dim*Io2No_V(UnitN_)*No2Io_V(UnitT_)

       ! Calculating the source term due to photoionization
       RhoDotPhoto= ReactionRate_H2O*(MassFluid_I(1)) !in this case the mass of the fluid is 18 amu      
       !the ion and electron temperature
       kTi = State_VGB(p_,i,j,k,iBlock)/ totalNumRho/2.0
       kTi_dim=kTi*No2Io_V(UnitTemperature_)

       !Calculating recombination rates, it is actually dissociative recombination rates are in Nagy's book p.228/229

       if( kTi_dim.LT.800)then
          RecombRate_dim=1.57e-5*(kTi_dim)**(-0.569)
       elseif(kTi_dim.LT.4000)then
          RecombRate_dim=4.73e-5*(kTi_dim)**(-0.74)
       else
          RecombRate_dim=1.03e-3*(kTi_dim)**(-1.111)
       end if
       RecombRate=RecombRate_dim*No2Io_V(UnitN_)*No2Io_V(UnitT_)

       !Calculating the loss term due to recombination
       RhoDotRecomb=RecombRate*totalNumRho*State_VGB(rho_,i,j,k,iBlock)              

       !Calculating the RhoDot due to charge exchange   
       !sigma_exchange=sigma_exchange_dim*Si2No_V(UnitX_)*Si2No_V(UnitX_)

       ! nion=State_VGB(rho_,i,j,k,iBlock)*No2Io_V(UnitN_)/(MassFluid_I(1))    
       rhotimeU=sqrt(State_VGB(Ux_,i,j,k,iBlock)*State_VGB(Ux_,i,j,k,iBlock)  &
            +State_VGB(Uy_,i,j,k,iBlock)*State_VGB(Uy_,i,j,k,iBlock)  &
            +State_VGB(Uz_,i,j,k,iBlock)*State_VGB(Uz_,i,j,k,iBlock))
       !rhotimeU_Si=rhotimeU*No2Si_V(UnitRhoU_)

       inv_rho = 1.00/State_VGB(rho_,i,j,k,iBlock)
       inv_rho2 = inv_rho**2
       uu2 =(State_VGB(Ux_,i,j,k,iBlock)*State_VGB(Ux_,i,j,k,iBlock)  &
            +State_VGB(Uy_,i,j,k,iBlock)*State_VGB(Uy_,i,j,k,iBlock)  &
            +State_VGB(Uz_,i,j,k,iBlock)*State_VGB(Uz_,i,j,k,iBlock)) &
            *inv_rho2

       ! U_Si=sqrt(uu2)*No2Si_V(UnitU_)
       !U_No=sqrt(uu2)
       !Rho_No=State_VGB(rho_,i,j,k,iBlock)

       Nu_exchange=sigma_exchange*rhotimeU*No2Si_V(UnitU_)*100 
       !Nu_exchange=sigma_exchange*sqrt(uu2)*No2Si_V(UnitU_)*100*State_VGB(rho_,i,j,k,iBlock)

       RhoDotchargeX=Nu_exchange*NumDenNeutral_VC(i,j,k)
       RhoDotchargeX=RhoDotchargeX*Io2No_V(UnitN_)

       !Calculating impact ionization
       !write(*,*)'Ti_dim=',kTi_dim
       if (kTi_dim.ge.1.0e4.and.kTi_dim.le.1.0e7)then
          if(kTi_dim.le.1.0e5)then
             h=int(kTi_dim/1.04 -1.0)
          else if(kTi_dim.le.1.0e6)then
             h=int(kTi_dim/1.0e5+9.0)
          else if(kTi_dim.le.5.0e6)then
             h=int(kTi_dim/1.0e6+19.0)
          else if(kTi_dim.le.7.0e6)then
             h=int(kTi_dim/1.0e6+18.0)
          else 
             h=25
          end if
          !write(*,*)'before calculating kappa'
          !write(*,*)'h',h
          !write(*,*)'Te(22)=',Te(22)
          !write(*,*)'Te(23)=',Te(23)
          !write(*,*)'(Te(h+1)-Te(h)',(Te(h+1)-Te(h))
          kappa = Impact(h)+&
               (Impact(h+1)-Impact(h))*(kTi_dim-Te(h))/(Te(h+1)-Te(h))
          !write(*,*)'after calculating kappa'
       else
          kappa=0.0
       end if

       ImpactRate_dim=kappa*NumDenNeutral_VC(i,j,k)*ne_hot
       RhoDotImpact=ImpactRate_dim*Io2No_V(UnitN_)*No2Io_V(UnitT_)*(MassFluid_I(1))

       !Calculating the general source and Loss term
       !RhoDotImpact=0.0
       !RhoDotchargeX=0.0
       RhoDot=RhoDotPhoto+RhoDotchargeX+RhoDotImpact
       RhoDotL=RhoDotRecomb+RhoDotchargeX
       RhoNumDot=RhoDot/(MassFluid_I(1))
       RhoNumDotL=RhoDotL/(MassFluid_I(1))
       RhoDotLx=RhoDotL/State_VGB(rho_,i,j,k,iBlock)

       !Now calculate the collision frequency
       !nu0=5.4745e-10
       Nu_C_dim=nu0*NumDenNeutral_VC(i,j,k)!the units are s-1

       !Normalizing the collision frequency
       Nu_C=Nu_C_dim*No2Io_V(UnitT_)

       !Calculating the source terms due to photoionization and collisions 

!!$       inv_rho = 1.00/State_VGB(rho_,i,j,k,iBlock)
!!$       inv_rho2 = inv_rho**2
!!$       uu2 =(State_VGB(Ux_,i,j,k,iBlock)*State_VGB(Ux_,i,j,k,iBlock)  &
!!$            +State_VGB(Uy_,i,j,k,iBlock)*State_VGB(Uy_,i,j,k,iBlock)  &
!!$            +State_VGB(Uz_,i,j,k,iBlock)*State_VGB(Uz_,i,j,k,iBlock)) &
!!$            *inv_rho2
!!$       U_Si=sqrt(uu2)*No2Si_V(UnitU_)

       Source_VC(rho_,i,j,k) = Source_VC(rho_,i,j,k) + RhoDot-RhoDotL
       Source_VC(RhoUx_:RhoUz_,i,j,k) = Source_VC(RhoUx_:RhoUz_,i,j,k) &
            -Nu_C*State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock)&
            -RhoDotLx*State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock)

       kTn = kTn_dim*Io2No_V(UnitTemperature_)


       Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
            -0.5*State_VGB(rho_,i,j,k,iBlock)*uu2*Nu_C     &
            +inv_gm1*(RhoNumDot*kTn-RhoNumDotL*kTi) &
            +1.5*totalNumRho*(kTn-kTi)*Nu_C&
            -0.50*uu2*RhoDotL


       Source_VC(p_,i,j,k) = Source_VC(p_,i,j,k) &
            +0.5*gm1*State_VGB(rho_,i,j,k,iBlock)*uu2*Nu_C  &
            +(RhoNumDot*kTn-RhoNumDotL*kTi) &
            +0.50*(gm1)*uu2*(RhoDot) &
            +1.5*gm1*totalNumRho*(kTn-KTi)*Nu_C 


       !write(*,*)'RhoDotPhoto=',RhoDotPhoto   
!!$    write(*,*)'RhoDotChargeX=',RhoDotChargeX
!!$    !write(*,*)'charge exchangecrosssection=',sigma_exchange
!!$    write(*,*)'Rho*U=',rhotimeU      
       !write(*,*)'RhoDotImpact=',RhoDotImpact
       !write(*,*)'RhoDotRecomb=',RhoDotRecomb
!!$    write(*,*)'RhoDottotal=',RhoDot   
    end do; end do; end do


    !if(oktest_me)then

!!$    write(*,*)'Blocknumber= ', iBlock
!!$       write(*,*)'Block_radius= ', R_BLK(iTest,jTest,kTest,iBlock)

!!$    write(*,*) 'neutral_density= ',neutral_density(R_BLK(iTest,jTest,kTest,iBlock))
!!$    uu2 = sum(( State_VGB(Ux_:Uz_,itest,jtest,ktest,iBlock) &
!!$         / State_VGB(rho_,itest,jtest,ktest,iBlock))**2)
!!$    !write(*,*)'uu2=', uu2
!!$    write(*,*)'kTi=   ',kTi
!!$    write(*,*)'kTn=   ',kTn
!!$   
!!$    write(*,*)'After update'
!!$    write(*,*)' Source_VC(Rho)=',Source_VC(rho_,iTest,jTest,kTest)
!!$    write(*,*)' Source_VC(Ux,Uy,Uz)=',Source_VC(RhoUx_:RhoUz_,iTest,jTest,kTest)
!!$    write(*,*)' Source_VC(Energy)=',Source_VC(Energy_,iTest,jTest,kTest)
!!$    write(*,*)' Source_VC(P)=',Source_VC(p_,iTest,jTest,kTest)
!!$  
    !  end if


    !Now I am making a test to see if the source terms
    !write(*,*) 'I just went through the loop in user_source for the block number',iBlock
    !Do I need to do stop_user?
    ! call stop_user(Name)

  end subroutine user_calc_sources

  !==========================================================================
  subroutine enceladus_input(iBlock)

    use ModInterpolate, ONLY: trilinear
    use ModPhysics
    use ModGeometry, ONLY: Xyz_DGB, R_BLK

    integer, intent(in) :: iBlock

    integer,parameter::imax=49
    integer,parameter::jmax=49
    integer,parameter::kmax=49
    character(2000)::line
    real,dimension(0:imax,0:jmax,0:kmax)::density
    !real,dimension(0:26)::Impact,Te
    real::Xyz_D(3)
    real::x,y,z,dx,x_BLK_Si,y_BLK_Si,z_BLK_Si, NumDenNeutral_VC_Si
    real,parameter::dx_dim=120.0e3!the step in m

    integer::m,n,p,i,j,k,h
    !-------------------------------------------------------------------------

    open(130, file='DensityH2O.dat ', status='old')
    do m=0,1
       read(130,'(a)')line
    end do

    do m=0,imax
       do n=0,jmax
          do p=0,kmax
             read(130,*)x,y,z,density(m,n,p)
          end do
       end do
    end do
    dx=dx_dim*Si2No_V(UnitX_)
    do k=1,nK; do j=1,nJ;do i=1,nI
       Xyz_D(1)=Xyz_DGB(x_,i,j,k,iBlock)/dx+25.0
       Xyz_D(2)=Xyz_DGB(y_,i,j,k,iBlock)/dx+25.0
       Xyz_D(3)=Xyz_DGB(z_,i,j,k,iBlock)/dx+25.0
       !x_BLK_Si=Xyz_DGB(x_,i,j,k,iBlock)*No2Si_V(UnitX_)
       !y_BLK_Si=Xyz_DGB(y_,i,j,k,iBlock)*No2Si_V(UnitX_)
       !z_BLK_Si=Xyz_DGB(z_,i,j,k,iBlock)*No2Si_V(UnitX_)
       if (any(Xyz_D<(/0, 0, 0/)) .or. &
            any(Xyz_D>(/iMax,jMax,kMax/)))then
          NumDenNeutral_VC(i,j,k)=0.0
       else
          NumDenNeutral_VC(i,j,k)=trilinear(density,0,imax,0,jmax,0,kmax,Xyz_D)
          NumDenNeutral_VC_Si=NumDenNeutral_VC(i,j,k)
          !write(*,*)'Si2No_V(UnitN_)',Si2No_V(UnitN_)
          !write(*,*)'x_BLK_Si=',x_BLK_Si
          !write(*,*)'Xyz_D(1)=',Xyz_D(1)
          !write(*,*)'y_BLK_Si=',y_BLK_Si
          !write(*,*)'Xyz_D(2)=',Xyz_D(2)
          !write(*,*)'z_BLK_Si=',z_BLK_Si
          !write(*,*)'Xyz_D(3)=',Xyz_D(3)
          !write(*,*)'NumDenValieryinterpolated',NumDenNeutral_VC_Si!But I want it in Si to be able to read it
       end if
    end do; end do; end do
    NumDenNeutral_VC=NumDenNeutral_VC*Si2No_V(UnitN_)
    close (130)

!!$write(*,*)'density(-3000,-3000, -3000)',density(0,0,0)
!!$write(*,*)'density(0,0,0)',density(25,25,25)
!!$write(*,*)'density(2880,2880,2880)', density(49,49,49)

    open(150,file='Impact_rates.dat')
    !write(*,*)'Impact_rates has been opened'
    read(150,'(a)')line
    do h=0,26
       read(150,*)Te(h),Impact(h)
    end do
    !write(*,*)'Te(22)=',Te(22)
    !write(*,*)'Te(23)=',Te(23)
    close(150)

    !write(*,*)'Impact_Rates has been read'
  end subroutine enceladus_input
  !========================================================================

  !=========================================================================
  real function neutral_density(R0)
    use ModPhysics, ONLY :Rbody,cZero

    real, intent(in) :: R0
    real, parameter:: n0_H2O=2.8143082e+7

    !This is obtained by interpolating Mike combi's number at closest approach
    real, parameter:: HNuH2O=0.2

    !-----------------------------------------------------------------------
    neutral_density = 0.0

    if( R0 >= 0.9*Rbody .and. R0< 3.0*Rbody ) &
         neutral_density = n0_H2O* exp(-(R0-Rbody)/HNuH2O)

  end function neutral_density


end module ModUser
!==============================================================================
