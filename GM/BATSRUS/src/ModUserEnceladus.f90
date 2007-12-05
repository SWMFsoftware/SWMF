!^CFG COPYRIGHT UM
!this file contains the ModUser foe a 1 species model of inosphereic mgnetospheric interaction in Enceladus, it contains Hptoionization of H2O, recombination, charge exchange and electron impact. 
!==============================================================================

module ModUser
  ! This is the user module for Enceladus
  use ModSize
  use ModUserEmpty,               &
        IMPLEMENTED1 => user_calc_sources


  include 'user_module.h' !list of public methods

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real, parameter :: VersionUserModule = 1.1
  real::NumDenNeutral_VC(nI,nJ,nK)
  character (len=*), parameter :: NameUserModule = &
       'Enceladus 1 species MHD code, Dalal Najib'

contains

  subroutine user_calc_sources

    use ModVarIndexes
    use ModAdvance, ONLY: State_VGB, Source_VC, Energy_
    use ModProcMH,   ONLY: iProc
    use ModMain, ONLY: PROCTEST,GLOBALBLK,BLKTEST, nI, nJ,nK,iTest,jTest,kTest 
    use ModPhysics, ONLY: Rbody, inv_gm1, gm1, Si2No_V,No2Io_V, Io2No_V,UnitN_, UnitT_,UnitTemperature_,UnitX_
    use ModProcMH,   ONLY: iProc
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, R_BLK, vInv_CB 
    use ModBlockData,ONLY: use_block_data, put_block_data, get_block_data

    character (len=*), parameter :: Name='user_calc_sources'

    integer::iBlock,i,j,k,iBlockLast=-1
    real:: R0
   
    ! to convert from eV to Kelvins 1eV=1.1604e+4 K
    real,parameter::kTn_dim = 180 !in K
    real, parameter::ne_hot=0.2!density of hot electron that cause imapact ionization
    real,parameter::kappa=1.2e-8 !electron impact rate from Burger only for electron of Te=12.5eV
    !real,parameter:: PhotoIonRate_H2O=9.1e-9 !Solarmax condition
    real,parameter::PhotoIonRate_H2O=3.6e-9 !Solarmin condition
    !real,parameter::PhotoIonRate_H2O=3.6e-7 !Solarmin condition test
    real, parameter:: sigma_exchange_dim=8.1e-20!cross section for charge exchange
    real::kTi, kTi_dim, kTn,sigma_exchange 
    real::inv_rho, inv_rho2, uu2
    real::totalNumRho, totalSourceNumRho
    real:: RhoDot,RhoNumDot, RhoDotPhoto,RhoDotImpact,RhoDotRecomb,RhoDotchargeX, RhoDotL, RhoNumDotL, RhoDotLx, RhoNumDotLx
    real:: ReactionRate_H2O_dim, ReactionRate_H2O,ImpactRate_dim, RecombRate_dim, RecombRate 
    real:: Nu_C,Nu_C_dim, Nu_exchange!Collision frequency
    real,parameter::nu0=5.4745e-10!this is in cm^3/s
   

    !------------------------------------------------------------------
    iBlock = globalBlk
    write(*,*)'iBlock',iBlock
    
!!$
!!$    write(*,*)'before update'
!!$    write(*,*)' Source_VC(Rho)=',Source_VC(rho_,iTest,jTest,kTest)
!!$    write(*,*)' Source_VC(Ux,Uy,Uz)=',Source_VC(RhoUx_:RhoUz_,iTest,jTest,kTest)
!!$    write(*,*)' Source_VC(Energy)=',Source_VC(Energy_,iTest,jTest,kTest)
!!$    write(*,*)' Source_VC(P)=',Source_VC(p_,iTest,jTest,kTest)
    
    do k = 1, nK ;   do j = 1, nJ ;  do i = 1, nI

!!$ !put the interpolated data in ModBlockData
!!$       if(iBlock /= iBlockLast)then
!!$       iBlockLast = iBlock
!!$       if(use_block_data(iBlock))then
!!$          !call get_block_data(iBlock, nI, nJ, nK, Nu_C)
!!$          call get_block_data(iBlock, MaxNuSpecies, nI, nJ, nK, NumDenNeutral_VC)
!!$          !call get_block_data(iBlock, MaxSpecies, nI, nJ, nK, PhotoIonRate_VC)
!!$          !call get_block_data(iBlock, MaxSpecies, nI, nJ, nK, RecombRate_VC)
!!$       else
!!$          call enceladus_input(iBlock)
!!$          !call put_block_data(iBlock, nI, nJ, nK, Nu_C)
!!$          call put_block_data(iBlock, MaxNuSpecies, nI, nJ, nK, NumDenNeutral_VC)
!!$          !call put_block_data(iBlock, MaxSpecies, nI, nJ, nK, PhotoIonRate_VC)
!!$          !call put_block_data(iBlock, MaxSpecies, nI, nJ, nK, RecombRate_VC)
!!$       end if
!!$    end if

 !End of interpolation of Data
       
       R0=R_BLK(i,j,k,iBlock)
       write(*,*)'R_block',R0 
       totalNumRho=State_VGB(rho_,i,j,k,iBlock)/(MassFluid_I(1))
      
       ReactionRate_H2O_dim=PhotoIonRate_H2O*neutral_density(R0)

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
         sigma_exchange=sigma_exchange_dim*Si2No_V(UnitX_)*Si2No_V(UnitX_) 
         Nu_exchange=sigma_exchange*sqrt(State_VGB(Ux_,i,j,k,iBlock)*State_VGB(Ux_,i,j,k,iBlock)  &
            +State_VGB(Uy_,i,j,k,iBlock)*State_VGB(Uy_,i,j,k,iBlock)  &
            +State_VGB(Uz_,i,j,k,iBlock)*State_VGB(Uz_,i,j,k,iBlock)) 
         RhoDotchargeX=Nu_exchange*neutral_density(R0)
         RhoDotchargeX=RhoDotchargeX*Io2No_V(UnitN_)
         
       !Calculating impact ionization
         ImpactRate_dim=kappa*neutral_density(R0)*ne_hot
         RhoDotImpact=ImpactRate_dim*Io2No_V(UnitN_)*No2Io_V(UnitT_)
   
    !Calculating the general source and Loss term
       RhoDot=RhoDotPhoto+RhoDotchargeX+RhoDotImpact
       RhoDotL=RhoDotRecomb+RhoDotchargeX
       RhoNumDot=RhoDot/(MassFluid_I(1))
       RhoNumDotL=RhoDotL/(MassFluid_I(1))
       RhoDotLx=RhoDotL/State_VGB(rho_,i,j,k,iBlock)

       !Now calculate the collision frequency
       !nu0=5.4745e-10
       Nu_C_dim=nu0*neutral_density(R0)!the units are s-1

       !Normalizing the collision frequency
       Nu_C=Nu_C_dim*No2Io_V(UnitT_)

       !Calculating the source terms due to photoionization and collisions 

       inv_rho = 1.00/State_VGB(rho_,i,j,k,iBlock)
       inv_rho2 = inv_rho**2
       uu2 =(State_VGB(Ux_,i,j,k,iBlock)*State_VGB(Ux_,i,j,k,iBlock)  &
            +State_VGB(Uy_,i,j,k,iBlock)*State_VGB(Uy_,i,j,k,iBlock)  &
            +State_VGB(Uz_,i,j,k,iBlock)*State_VGB(Uz_,i,j,k,iBlock)) &
            *inv_rho2
      
       
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
            +totalNumRho*(kTn-KTi)*Nu_C 
       
    write(*,*)'RhoDotPhoto=',RhoDotPhoto
    write(*,*)'RhoDotChargeX=',RhoDotChargeX
    write(*,*)'RhoDotImpact=',RhoDotImpact
    write(*,*)'RhoDotRecomb=',RhoDotRecomb
    write(*,*)'RhoDottotal=',RhoDot   
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
    !write(*,*) 'I just went through the loop in user_source for the block number',GLOBALBLK
    !Do I need to do stop_user?
    ! call stop_user(Name)

  end subroutine user_calc_sources

  !==========================================================================
!!$subroutine enceladus_input
!!$use ModInterpolate, ONLY: trilinear
!!$use ModPhysics
!!$use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, R_BLK
!!$!-------------------------------------------------------------------------
!!$integer,parameter::imax=49
!!$integer,parameter::jmax=49
!!$integer,parameter::kmax=49
!!$character(2000)::line
!!$real,dimension(0:imax,0:jmax,0:kmax)::density
!!$real::Xyz_D(3)
!!$real::x,y,z,dx
!!$real::parameter:dx_dim=120.0e3!the step in m
!!$integer::m,n,p,i,j,k
!!$!===================================================================
!!$open(130, file='DensityH2O.dat ', status='old')
!!$do m=0,1
!!$read(130,'(a)')line
!!$end do
!!$
!!$do m=0,imax
!!$   do n=0,jmax
!!$      do p=0,kmax
!!$         read(130,*)x,y,z,density(m,n,p)
!!$      end do
!!$   end do
!!$end do
!!$dx=dx_dim*Si2No_V(UnitX_)
!!$Xyz_D(1)=x_BLK(i,j,k,iBlock)/dx
!!$Xyz_D(2)=y_BLK(i,j,k,iBlock)/dx
!!$Xyz_D(3)=z_BLK(i,j,k,iBlock)/dx
!!$NumDenNeutral_VC(i,j,k)=trilinear(density,0,imax,0,jmax,0,kmax,Xyz_D)   
!!$close (130)
!!$
!write(*,*)'density(-3000,-3000, -3000)',density(0,0,0)
!write(*,*)'density(0,0,0)',density(25,25,25)
!write(*,*)'density(2880,2880,2880)', density(49,49,49)
!!$
!!$
!!$end subroutine enceladus_input
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
