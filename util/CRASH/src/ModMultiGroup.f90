module CRASH_ModMultiGroup
  use CRASH_ModIonMix
  use CRASH_ModOpacityVoigt, ONLY : line_width, voigt_profile, UseVoigt
  use CRASH_ModAtomicDataMix,ONLY : nMix, nZ_I, nExcitation, nMixMax
  use CRASH_ModAtomicDataMix,ONLY : Concentration_I !(1:nMixMax)
  use CRASH_ModAtomicDataMix,ONLY : IonizPotential_II !(1:nZMax,1:nMixMax)
  use CRASH_ModAtomicMass,   ONLY : cAtomicMass_I !(1:nZMax)
  use CRASH_ModStatSumMix,   ONLY : Population_II !(0:nZMax,1:nMixMax)
  use CRASH_ModExcitationData, ONLY:n_ground !(iZ,nZ)
  use CRASH_ModExcitationData, ONLY:n_screened !(iZ,nZ)
  use CRASH_ModExcitation,   ONLY : Partition_III !(nExcitation,0:nZMax  ,nMixMax)
  use CRASH_ModExcitation,   ONLY : ExcitationEnergy_III !nExcitation,0:nZMax-1,nMixMax)
  use CRASH_ModStatSumMix,   ONLY : Na, Te, zAv
  use CRASH_ModStatSumMix,   ONLY : iZMin_I !(1:nMixMax)
  use CRASH_ModStatSumMix,   ONLY : iZMax_I !(1:nMixMax)
  implicit none
  SAVE
  PRIVATE !Except

  !Public members
  public :: meshhv !Creates the grid of photon energies
  public :: lines  !Calculates the absorbtion and emission in lines


  !       nPhotonMax  - photon energy mesh points                                    
  !       nGroupMax  - opacity groups     
  !       nfrqbb  - number of photon energy points near a line center   
  !       at which absorption coefficients will be computed  

  integer,parameter:: nPhotonMax = 5000, nGroupMax = 100, nfrqbb = 9

  !Radiation frequency groups:
  integer :: nGrups = nGroupMax 
  real :: EnGrup(0:nGroupMax)

  integer,parameter::nptspg = 30

  ! Meshs to evaluate the absorbtion coefficients
  integer::nPhoton
  real ::PhotonEnergy_I(nPhotonMax)

  ! The absorbtion coefficients
  !       abscfs  -  array of absorption coefficients (cm**-1)            
  !       emscfs  -  array of emission coefficients (cm**-1)         
  !       sctcfs  -  array of scattering coefficients (cm**-1) 

  real,dimension(nPhotonMax) ::abscfs,sctcfs,emscfs

  !\
  ! LOGICALS
  !/
  !Determine if we account for the optical transitions with no change 
  !in the principal quantum number

  logical,parameter :: UseDeltaNEq0Transition = .false.
  
  !Determine if we account for bound-free transition, in which the
  !bound electron is photo-ionized from inner shells ("core electron")
  logical,parameter :: UseCoreElectron = .false.

  !Determine if we account for the corrections used in HYADES for
  !the oscillator strength

  logical :: UseHYADESCorrection4Strength = .true.
 
  !Switch to add or not to add the controbutions from the line core if
  !it is added somewhere else
  logical,public :: DoNotAddLineCore = .true.
  
  logical :: UseBremsstrahlung = .true.
  logical :: UsePhotoionization = .true.
  logical :: UseScattering      = .false.
contains
  !======================================================================
  real function oscillator_strength(nI,nF)
    integer,intent(in):: nI,nF

    ! ... oscillator strength taken from Zeldovich      
    !     & Raizer for the upward transition from "ni" to "nf")             
    !      oscstr(ni,nf) = 1.96 / ni**5 / nf**3 / (1./ni**2-1./nf**2)**3
    !Here 1.96 = 32/(3\pi\sqrt{3})
    
    !Corrections used in HYADES 
    !f(1,2) = 0.4246 
    !f(1,3) = 0.0808 
    !f(1,4) = 0.0296 
    !f(1,5) = 0.1418 
    !f(2,3) = 0.6500 
    !f(2,4) = 0.1214 
    !f(2,5) = 0.0452 
    !f(3,4) = 0.8580 
    !f(3,5) = 0.1530 
    !f(4,5) = 1.058 
    real,parameter,dimension(1:4,2:5):: HYADESCorrection_II = reshape( (/&
         0.4246 , 0.0    , 0.0    , 0.0,   &  !nF=2
         0.0808 , 0.6500 , 0.0    , 0.0,   &  !nF=3
         0.0296 , 0.1214 , 0.8580 , 0.0,   &  !nF=4
         0.1418 , 0.0452 , 0.1530 , 1.058  &  !nF=5
                                         /),(/4,4/))
    !    nI =1  ! nI=2   ! nI=3   ! nI=4   !
    
    real :: rNI, rNF, rNI2, rNF2
    !---------------------------
    if(UseHYADESCorrection4Strength &
         .and. nI<=4 .and. nF<=5) then
       oscillator_strength = HYADESCorrection_II(nI,nF)
    else
       rNI = real(nI)
       rNF = real(nF)
       rNI2 = rNI * rNI
       rNF2 = rNF * rNF

       oscillator_strength = 1.96 * rNI * (rNF * rNF2)/ (rNF2 - rNI2)**3 
    end if
  end function oscillator_strength

  !======================================================================
  subroutine lines (ephot,abstot,emstot )       

    !                                                                       
    ! ... compute the contribution to the absorption coefficient from       
    !     all lines (bound-bound transitions)                               
    !                                                                       
    ! ... input variables:                                                 
    !       Te      -  plasma temperature (eV)                             
    !       densnn  -  number density of all nuclei (cm**-3)               
    !       densne  -  electron density (cm**-3)                            
    !       ephot   -  photon energy (eV)                                   
    real,intent(in) :: ephot
    real :: densnn,densne
    !                                                                       
    ! ... output variables:                                                 
    !       abstot  -  absorption coefficient due to lines (cm**-1)         
    !       emstot  -  emission coefficient due to lines (cm**-1)           
    !   

    real,intent(out)::absTot,emsTot                            

    integer:: iSav !Integer to count the total number of lines involved
    integer,parameter:: nSavMax = 200 !The upper bound for iSav

    !\
    ! Loop variables
    !/
    integer :: iMix, & !runs over the mixture component
         iZ,   & !runs over the charge number
         iN,   & ! runs over the quntum principal number for the lower level
         iNUpper,& !the same, for upper level 
         nBound, &!Number of bound electrons
         nGround  !For a given iZ, the principal number of the ground state.

    real:: denlq  !Density of ions of a given sort, [cm-3]
    real:: denlqn !The same, for the ion in the lower state, for a given transition
    real:: dnlqnp !The same, for the ion in the upper state.

    !The controbutions to the total absorbtion and the total emission,
    !calculated by the 'abslin' subroutine
    real:: abscof, emscof 
    !-------------------



    iSav = 0                                                          
    abstot = 0.                                                       
    emstot = 0.                      

    DensNN = Na * 1.0e-6  !To convert to cm-3
    DensNe = DensNN * zAv !Electron density, in cm-3

    ! ... loop over gas species                                             
    IMIXLOOP: do iMix =1,nMix                                              
       if (Concentration_I(iMix) .lt. con(2) ) CYCLE IMIXLOOP

       ! ...    loop over ionization states                         
       IZLOOP: do  iZ=iZMin_I(iMix),iZMax_I(iMix)                                                                            
          if ( Concentration_I(iMix)*Population_II(iZ,iMix) .lt. con(3))&
               CYCLE IZLOOP

          ! ... find the principal quantum number of the valence electrons  
          !           in their ground state                                       
                            
          nGround = n_ground(iZ,nZ_I(iMix))                                   

          denlq  = densnn * Concentration_I(iMix) * Population_II(iZ,iMix)          

          ! ...  loop over initial quantum states                            
          INLOOP: do iN=nGround,nExcitation-1                              

             if ( Concentration_I(iMix) * &
                  Population_II(iZ,iMix)* &
                  Partition_III(iN,iZ,iMix)  &
                  <  con(4) ) CYCLE  INLOOP                           

             ! ... compute the density of ions with an electron in          
             !             in state "n"                                             
             denlqn = denlq * Partition_III(iN,iZ,iMix)                  

             !         loop over final quantum states                           
             do iNUpper=iN,nExcitation                                 

                dnlqnp = denlq * Partition_III(iN,iZ,iMix)                 
                call abslin   

                if ( abscof .ne. 0. ) then                            
                   abstot = abstot + abscof                           
                   emstot = emstot + emscof                           
                   isav = isav + 1                                    
                   if ( isav==nSavMax ) then                            
                      call CON_stop(' you are keeping track of too many lines') 
                   endif
                endif
             end do
          end do INLOOP
       end do IZLOOP
    end do IMIXLOOP
  contains
    !=============================================================
    subroutine abslin   
      use ModConst,ONLY: cPlanckH,cEV

      ! ... computes the absorption coefficient from a particular            
      !     bound-bound transition from quantum state "n" to "np".           
      !                                                                      
      ! ... input variables                                                  
      !       Te      =  plasma temperature (eV)                              
      !       densnn  =  number density of all nuclei (cm**-3)                
      !       densne  =  electron density (cm**-3)                            
      !       denlqn  =  number density of nuclei of the "l"th gas, in        
      !                  the "q"th ionization state, with an electron in     
      !                  the "n"th quantum state                             
      !       dnlqnp  =  number density of nuclei of the "l"th gas, in       
      !                  the "q"th ionization state, with an electron in     
      !                  the "np"th quantum state                             
      !       ephot   =  photon energy (eV)                                   
      !       potiz   =  ionization potential for the ion in its ground      
      !                  state (eV)                                          
      !       atomwt  =  atomic weight of the ion (amu)                       
      !       iN       =  initial principal quantum number for transition     
      !       iNUpper  =  final principal quantum number                      
      !       nprin0  =  principal quantum number of the valence electrons   
      !                  in their ground state                              
      !       nbound  =  number of electrons bound to the ion                
      !       izgas   =  atomic number of the ion                           
      !                                                                     
      ! ... output variable                                                  
      !       abscof  =  absorption coefficient (cm**-1)                     
      !       emscof  =  emission coefficient (cm**-1)                       
  
     
      real:: ennp !The transition energy
      real:: gamma,avoigt,dnudop,vvoigt  !Line width parameters

      !Difference between the line center and the photon energy, 
      !for which to calculate the absorbtion

      real:: DeltaNu 
      
      !       hplank  -  4.136e-15   Planck's constant (eV sec)
      real,parameter:: cHPlanckEV = cPlanckH/cEV

      !The oscillator strength
      real :: OscillatorStrength

      !The line height at a given frequency (the line shape)
      real:: Shape

      !The absorbtion coefficient, to be corrected for stimulated emission
      real:: Alpha

      !\
      !Coefficients to account for stimulated emission
      !/
      real::corsea,corsee,corrse

      !Mics:
      real::ex1

      !Misc, are not used at the time
      real:: fnn,gnn
      !-------------------------



      ! ... compute the transition energy                                     
      if ( iN.eq.iNUpper ) then    
         nBound = nZ_I(iMix) - iZ                                            
         call bbneq0 (iZ,nbound, &                                
              ennp,fnn,gnn )                 
         ! ...    some transitions are not allowed                               
         if ( ennp.le.0) return                     
      else                                                              
         ennp = ExcitationEnergy_III(iNUpper,iZ,iMix) - ExcitationEnergy_III(iN,iZ,iMix)                  
      endif

      ! ... compute the line widths for natural, Doppler, and pressure        
      !     broadening to be used with Lorentzian line shape                  
      call line_width ( Te,densnn,ennp,cAtomicMass_I(nZ_I(iMix)),  &                            
           gamma,avoigt,dnudop )        

      ! ... compute this contribution if the photon energy is not far         
      !     from the line center                                              
      DeltaNu = ( ennp-ephot ) / cHPlanckEV                                  
      if ( abs( deltaNu ) > con(5)*gamma ) then 
         !Line is too far, ignore it
         abscof = 0.                                                    
         emscof = 0.                                                    
         return
      endif

      ! ...    if the contribution from line "cores" will be added            
      !        elsewhere (-opacbb-), then use the value at the                
      !        line core boundary.                                            
      if (DoNotAddLineCore .and. abs(DeltaNu).lt.con(6)*gamma ) then     
         DeltaNu = con(6)*gamma                                       
      endif

      ! isw(14):Voigt or Lorentzian line profile (0=>V
      ! if ( isw(14).eq.0 ) then     

      if(UseVoigt)then
         ! ...       use Voigt profile                                           
         vvoigt = deltaNu / dnudop                                    
         shape = voigt_profile ( avoigt,vvoigt ) / dnudop / 1.7725           
      else                                                           
         ! ...       use Lorentzian profile                                      
         shape = (gamma/39.48) / (DeltaNu**2 + (gamma/12.57)**2)      
      endif

      !Commented out lines below are relevant to non-LTE
      !case only

      ex1 = exp( -ephot/Te )                                         
      !non-LTE: ex2 = exp( - ennp/Te )                                         

      if ( iN==iNUpper) then                                            

         alpha = 2.65e-2 * fnn * shape * denlqn * ( 1.-ex1 )         

         ! ...       correct for the "effective" stimulated emission to          
         !           give the proper form for the cooling rate                   
         !non-LTE:if ( isw(6) .ne. 1 ) then                                   
            ! ...          general form                                             
         !non-LTE:   xnjoni = gnn * densne * ex2 /    &                        
         !non-LTE        ( gnn * densne + ennp**3 * sqrt(Te) * 2.74e12 )   
         !non-LTE   corsee = gnn * densne /          &                       
         !non-LTE        ( gnn * densne + ennp**3 * sqrt(Te) * 2.74e12 )   
         !non-LTE   corsea = ( 1.-xnjoni ) / ( 1.-ex2 )  
         !non-LTE else                                                        
            ! ...          LTE assumed                                              
            corsea = 1.0                                            
            corsee = 1.0                                             
         !non-LTE: endif

                       
         abscof = alpha * corsea                                     
         emscof = alpha * corsee                                     

      else                                                          

         OscillatorStrength = oscillator_strength ( iN, iNUpper) 
                       
         !non-LTE: dum = ( iN * iN * dnlqnp ) / ( np*np*denlqn )
         !in LTE dum=ex2                     

         ! ...       correct for stimulated emission  
                           
         !non-LTE: corrse = ( 1.-dum ) * ( 1.-ex1 ) / ( 1.-ex2 ) 

         corrse = 1.0 - ex1
              
         alpha  = 2.65e-2 * OscillatorStrength * shape * denlqn                    
         abscof = alpha * corrse   
                                  
         !non-LTE: if ( ex2 .gt. 0. ) then                                     
         !non-LTE:   emscof = alpha * dum * ( 1.-ex1 ) / ex2                   
         !non-LTE: else                                                        
            emscof = abscof                                           
         !non-LTE: endif

      endif


    end subroutine abslin
    !===================
  end subroutine lines
  !=============================================
  subroutine bbneq0 (iZGas,nBound, enn,fnn,gnn )              

    ! ... calculates the effective energy, oscillator strength, and         
    !     gaunt factor for bound-bound transitions in which the             
    !     principal quantum number does not change.  This approximation     
    !     was taken from Post, et.al., Atomic & Nuclear Data Tables,       
    !     vol. 20, p.397 (1977).                                            

    ! ... input variables                                            
    !       Te      -  plasma temperature (eV)                              
    !       izgas   -  nuclear charge of the ion                            
    !       nbound  -  number of electrons bound to the ion

    integer,intent(in) :: iZGas !nZ
    integer,intent(in) :: nBound

    !                                                                       
    ! ... output variable                                                   
    !       enn  -  effective energy of the transition (eV)                 
    !       fnn  -  effective oscillator strength                          
    !       gnn  -  effective gaunt factor                                  

    real, intent(out) :: enn, fnn, gnn

    !       bbnton  -  constants for bound-bound transitions in which       
    !                  the principal quantum number does not change         

    ! ... parameters for bound-bound transitions when "delta n" = 0.      
    !     these values are taken from Post, et.al., Atomic and Nuclear     
    !     Data Tables, Vol. 20, p.397 (1977).                              

    real,parameter:: bbnton(55,6) = reshape( (/&                                       
         0.,  0.,-.02,-.01,-.04, .04,-.01,-.01, .01,  0., &           
         .01, .10, .25, .17, .08,-.02,-.14,-.27,-.29,-.30, &           
         -.30,-.29,-.27,-.24,-.20,-.14,-.08,  0., .97,1.96, &            
         1.92,1.89,1.86,1.83,1.78,1.73,1.41,1.05, .67, .26, &          
         -.17,-.64,-1.14,-1.67,-2.26,-2.88,-2.90,-2.83,-2.72,-2.61,&
         -2.45,-2.27,-2.05,-1.81,-1.55,&
         0.,  0.,2.00,4.33,4.68,3.60,2.85,2.40,1.30,  0., &  
         10.5,20.0,16.4,24.7,34.1,44.8,56.8,70.3,68.6,65.8, & 
         61.9,56.9,50.7,45.5,34.4,24.2,12.8,  0.,-25.8,-53.1,&  
         -28.1,-2.97,23.0,50.3,79.1,109.,142.,176.,214.,257.,&         
         300.,348.,401.,457.,518.,584.,572.,551.,525.,497.,  &          
         463.,427.,385.,340.,291.,&
         0.,  0.,2.04,4.49,6.80,8.16,6.80,10.6,13.1,  0.,  &         
         2.30,3.88,5.71,5.44,8.16,6.80,8.30,4.32,5.11,6.04,  &         
         7.08,8.12,9.69,11.3,13.7,14.9,17.1,  0.,.015,.019,  &         
         .056,.120,.213,.334,.466,.666,.910,1.25,1.69,1.98,  &         
         3.03,3.92,5.19,6.40,8.05,9.80,10.8,11.9,13.2,14.7,  &         
         15.9,18.1,19.2,21.0,23.3, &
         1.,  1., 1.0, .93, .86, .80, .87, .77, .77,  1.,  &         
         .99, .90, .83, .87, .77, .82, .79,1.06,1.03,1.00,  &         
         .07, .94, .91, .88, .85, .83, .80,  1.,2.46,2.40,  &         
         2.14,1.96,1.83,1.72,1.64,1.57,1.49,1.41,1.34,1.30,  &         
         1.19,1.13,1.06,1.01, .95, .91, .89, .87, .85, .83,  &         
         .82, .80, .78, .76, .74, &
         0.,  0., 1.2, .91, .89, .92, .86, .87, .85,  0.,  &         
         .72, .78, .76, .78, .75, .70, .68, .65, .65, .70,  &         
         .7,  .7,  .7,  .7,  .7,  .7,  .7,  .7,  .7,  .7,  &         
         .7,  .7,  .7,  .7,  .7,  .7,  .7,  .7,  .7,  .7,  &         
         .7,  .7,  .7,  .7,  .7,  .7,  .7,  .7,  .7,  .7,  &         
         .7,  .7,  .7,  .7,  .7, & 
         0.,  0., .54, .77, .58, .57, .41, .63, .66,  0.,  &         
         .97, .96, .88, .86, .87, .85, .83, .81, .80, .80,  &         
         .8,  .8,  .8,  .8,  .8,  .8,  .8,  .8,  .8,  .8,  &         
         .8,  .8,  .8,  .8,  .8,  .8,  .8,  .8,  .8,  .8,  &        
         .8,  .8,  .8,  .8,  .8,  .8,  .8,  .8,  .8,  .8,  &         
         .8,  .8,  .8,  .8,  .8 /),(/55,6/))                                   

    integer:: nb !Limited from above nBound
    integer:: iZ 

    !Misc
    real :: c, ExpInt, xnn
    !----------------------------------
    if(.not.UseDeltaNEq0Transition)then
       enn = -1.0
       return
    end if

    iZ =izgas - nbound                                            
    nb = min( nbound , 55 )                                           

    fnn = bbnton(nb,1) + bbnton(nb,2) / izgas                         
    enn = bbnton(nb,3) * (iZ+1)**bbnton(nb,4)                             

    if ( enn .le. 0. ) return                                         

    xnn = enn / Te                                                   
    if ( iZ == 0 ) then                                               
       c = 0.06 * (sqrt( xnn )-2.) / (1.+xnn)                         
    else                                                              
       c = bbnton(nb,5) * ( 1. - bbnton(nb,6)/iZ )                
    endif

    ! ... compute the exponential integral (see Abramowitz & Stegun), and   
    !     multiply by exp(xnn)                                              
    if ( xnn .ge. 1. ) then                                           
       expint = (xnn**2+2.3347*xnn+0.25062) / &                        
            (xnn**2+3.3306*xnn+1.6815) / xnn                      
    else                                                              
       expint = -log( xnn ) - 0.57722 + 0.99999*xnn - 0.24991*xnn**2 &
            + 0.05520*xnn**3 - 0.00976*xnn**4 + 0.00108*xnn**5    
       expint = expint * exp( xnn )                                   
    endif

    gnn = c + 0.276 * expint                                          

  end subroutine bbneq0
  !==========================================
  subroutine meshhv  
    !                                                                       
    ! ... set a mesh of photon energies at which 
    !     we would like to evaluate  
    !     absorption coefficients                                           
    !   


    ! ... input variables:                                                  
    !       Te      -  plasma temperature (eV)                              
    !       densnn  -  number density of all nuclei (cm**-3)                
    !       densne  -  electron density (cm**-3)                            
    !       engrup  -  photon energy group boundaries (ngrups+1) (eV)       
    !       ngrups  -  number of photon energy groups                       
    !       nptspg  -  minimum number of mesh points per energy group       
    !

    real :: DensNN, DensNE
    
    real,parameter:: dPlus =  1.001, dMinus = 0.999
    !\
    !Loop variables
    !/
    integer :: iGroup    !Enumerates energy groups
    integer :: iSubGrid  !Photon energy grid for each energy group
    integer :: nMax      !

    integer :: iMix      !Over the mixture components
    integer :: iZ        !Ion charge 
    integer :: iN        !Principal quantum number
    integer :: iNUpper   !Principal quantum number for the excited state

    integer :: iPointPerLine

    !Enumerate the energies of the photoionization transitions
    !from inner shells
    integer :: nCore 

    !\
    !Grid parameters, for each group
    !/
    real :: HvMin, HvMax, dLogHv, LogHv

    !CutOff energy (for h\nu=\hbar\omega_{pe}
    real :: hvCut


    ! Averaged energy for transitions between the states with no
    ! change in the principal quantum number

    real :: ennp

    real::    Edge !Cutoff energy for photoionization from the ground state
    integer:: nBound !The number of bound electrons

    !A gap between the line centered frequency and two neighboring 
    !frequency points used to represent its profile.
    real:: dHNu
    !A width of the line
    real :: Gamma

    !The principal quantum number of the ion in the
    !ground state
    integer :: nGround 

    !Misc
    real::Dum1,Dum2,avoigt,dnudop !non-used output parameters
    real :: KK
    !------------------------------------                                                                 


    ! ... initialize variables     

    nPhoton = 0                                                        

    PhotonEnergy_I( 1:nPhotonMax ) = 0.0

    DensNN = Na * 1.0e-6
    DensNe = DensNN * zAv


    ! ... set up initial grid within each photon energy group
    do iGroup = 1,ngrups                                                

       hvmin = engrup(iGroup)                                             
       hvmax = engrup(iGroup+1)                                           
       dloghv = log( hvmax/hvmin ) / nptspg                           
       LogHv = log( hvmin ) - dloghv                                  
       if ( iGroup.lt.ngrups ) then                                       
          nmax = nptspg                                               
       else                                                           
          nmax = nptspg + 1                                           
       endif

       do  iSubGrid = 1, nmax                                                
          LogHv = LogHv + dloghv                                      
          PhotonEnergy_I(nPhoton+1) = exp( LogHv )                              
          nPhoton = nPhoton + 1                                           
       end do

    end do

    ! ... add 2 points near the plasma cutoff frequency                     

    hvcut = sqrt( densne / 7.25e20 )                                  
    if ( hvcut.gt.engrup(0) .and. hvcut.lt.engrup(ngrups) ) then    
       PhotonEnergy_I(nPhoton+1) = hvcut * dminus                               
       PhotonEnergy_I(nPhoton+2) = hvcut * dplus                                
       nPhoton = nPhoton + 2                                              
    endif

    ! ... add points near line centers and ionization edges                 

    do iMix=1,nMix                                              
       if ( Concentration_I(iMix) .lt. con(2) ) CYCLE                     

       do  iZ=iZMin_I(iMix),iZMax_I(iMix)                                           
          if ( Concentration_I(iMix)*Population_II(iZ,iMix) .lt. con(3) )CYCLE

          ! ...       find the principal quantum number 
          !           of the valence electrons  
          !           in their ground state                                       
          nbound = nZ_I(iMix) - iZ                               
          nGround = n_ground(iZ,nZ_I(iMix))                                     
          
          do iN = nGround, nExcitation-1                                
             if (Concentration_I(iMix) * &
                  Population_II(iZ,iMix) * &
                  Partition_III(iN, iZ, iMix)& 
                  < con(4) )CYCLE
             if(.not.DoNotAddLineCore)then        

                do  iNUpper = iN, nExcitation                              
                   ! ... calculate the energy of the transition             
                   if ( iN == iNUpper ) then                               
                      call bbneq0 ( nZ_I(iMix),nbound, &            
                           ennp,dum1,dum2 )  
                      if ( ennp.le.0) CYCLE
                   else                                                
                      ennp = ExcitationEnergy_III(iNUpper, iZ,iMix)- &
                           ExcitationEnergy_III(iN, iZ,iMix)
                   endif
                   if ( ennp .gt. engrup(0) .and. &                     
                        ennp .lt. engrup(ngrups) ) then              
                      
                      call line_width ( Te,densnn,ennp,cAtomicMass_I(nZ_I(iMix)), &      
                           gamma,avoigt,dnudop ) 
                      
                      ! ...set points at and near the line center           
                      dhnu = con(9) * 4.14e-15 * gamma / 12.57         
                      PhotonEnergy_I(nPhoton+1) = ennp                           
                      nPhoton = nPhoton + 1                                
                      
                      kk = 1.0
                      do  iPointPerLine = 1, nfrqbb / 2 - 1                                                                  
                         PhotonEnergy_I(nPhoton+1) = ennp - kk * dhnu             
                         PhotonEnergy_I(nPhoton+2) = ennp + kk * dhnu             
                         nPhoton = nPhoton + 2                              
                         kk = kk * 2.0
                      end do
                      
                      ! ... stop if too many mesh pts are requested          
                      if ( nPhoton .gt. nPhotonMax - 5 - nfrqbb )&           
                           stop ' too many pts in -meshhv-'              
                   end if
                end do !iNUpper

             end if
             ! ... add 2 points near the ionization edge                    
             edge = IonizPotential_II(iZ+1,iMix) - ExcitationEnergy_III(iN,iZ,iMix)   
             
             if ( edge.gt.engrup(0) .and. edge.lt.engrup(ngrups) )&  
                  then                                                  
                PhotonEnergy_I(nPhoton+1) = edge * dminus                       
                PhotonEnergy_I(nPhoton+2) = edge * dplus                        
                nPhoton = nPhoton + 2                                     
             endif
          
          end do   !Over iN 

       end do    !Over iZ

       ! ... add 2 points near photoionization edges of core electrons     
       if(.not.UseCoreElectron)CYCLE
       do ncore=1,nZ_I(iMix)-1                                   
          edge = IonizPotential_II(nZ_I(iMix)+1-ncore,iMix)                        
          if ( edge.gt.engrup(0) .and. &                              
               edge.lt.engrup(ngrups) ) then                        
             PhotonEnergy_I(nPhoton+1) = edge * dminus                          
             PhotonEnergy_I(nPhoton+2) = edge * dplus                           
             nPhoton = nPhoton + 2                                        
          endif
       end do
    end do !Over iMix

    ! ... now, arrange the photon energies in monotonically increasing      
    !     order                                                             
    call sort 
  contains
    !==============
    subroutine sort 
      integer:: nCycle, j, i, iMinus1
      real :: rSave
      logical:: DoRepeat    
      !-----------------
      ncycle = 0                                                        

      ! ... first, throw out those values < or = zero                         

      j = 0                                                             
      do  i=1,nPhoton                                                    
         if ( PhotonEnergy_I(i)<= 0.0 ) CYCLE                                     
         j = j + 1                                                   
         PhotonEnergy_I(j) = PhotonEnergy_I(i)                                             
      end do
      nPhoton = j                                                          
      DoRepeat = .true.                                                                
      do while(DoRepeat)
         DoRepeat = .false.
         do i = 2,nPhoton                                                 
            iMinus1 = i - 1                                                 
            if ( PhotonEnergy_I(iMinus1) .gt. PhotonEnergy_I(i) ) then                            
               rSave = PhotonEnergy_I(iMinus1)                                          
               PhotonEnergy_I(iMinus1) = PhotonEnergy_I(i)                                        
               PhotonEnergy_I(i) = rSave                                            
               DoRepeat = .true.                                           
            endif
         end do
         nCycle = nCycle + 1                                           

         ! ...    check for error                                                

         if ( ncycle .gt. nPhoton ) call CON_stop( ' error in -sort-')                                                                                                  
         ! ...    if any of the elements have been swapped, try again
      end do
    end subroutine sort
    !===============
  end subroutine meshhv
  !====================

  subroutine abscon  
    use CRASH_ModStatSumMix, ONLY: Z2_I
                                                                        
    ! ... this routine calculates the absorption, emission, and scattering 
    !     coefficients for an array of photon energies                      
    !                                                                      
    ! ... input variables:                                                  
    !       Te      -  plasma temperature (eV)                              
    !       densnn  -  number density of all nuclei (cm**-3)                
    !       densne  -  electron density (cm**-3)                            
    !       PhotonEnergy_I  -  array of photon energies at which the absorption
    !                  coefficient will be evaluated (eV)                   
    !       nPhoton   -  number of elements in the "PhotonEnergy_I" array   
    !                                                                      
    ! ... output variables:                                                 
    !       abscfs  -  array of absorption coefficients (cm**-1)            
    !       emscfs  -  array of emission coefficients (cm**-1)         
    !       sctcfs  -  array of scattering coefficients (cm**-1)          
    !                                                                      
                                                                       
    real:: DensNN, DensNe    ![cm -3]

    real,dimension(nMixMax) :: brems = 0.0, fotiza = 0.0, fotize = 0.0
    
    
    !\
    ! Loop variables
    !/
    integer :: iMix, & !runs over the mixture component
         iZ,   & !runs over the charge number
         iN,   & ! runs over the quntum principal number for the lower level
         nBound, &!Number of bound electrons
         nGround  !For a given iZ, the principal number of the ground state.
 
    integer::iPhot

    real :: eTransition, ennpot

    !\
    !Partial sums
    !/
    real:: sum1a, sum2a, sum1e, sum2e

    !\
    ! Bound-bound contributions
    !/
    real :: abslns,emslns

    !\
    ! Cutoff energy for the photon with the frequency correspondent to
    ! the plasma electron frequency
    !/
    real :: HNuCut

    !Misc: functions of the photon energy
    real :: exhvot,  photn3

    !Misc: Density times Z2
    real :: densi2

    !Misc: coefficients used in free-free gaunt factor
    real :: gam2lg, gBrem

    !Number of screened and valence electrons, as well as iZEff
    integer :: nScreened, nValence, iZEff, iQ

    ! densities to calculate photoionization
    real :: deni, eqdeni, dumden, conpi

    !\
    !Coefficients to account for stimulated emission
    !/
    real :: stcorr

    !Misc: Scattering parameters
    real :: ScatNe, tScatt, pScatt

    !Misc: Photoionization parameters
    real:: Const, XSec
    !---------------------------------------------------

    !Initialization
    
    abscfs  = 0.0  !-  array of absorption coefficients (cm**-1)            
    emscfs  = 0.0  !-  array of emission coefficients (cm**-1)         
    sctcfs  = 0.0  !-  array of scattering coefficients (cm**-1)     
    
    DensNN = Na * 1.0e-6
    DensNe = DensNN * zAv
    
    ! ... loop over photon energies
    do iphot=1,nPhoton                                              
                                                                        
       exhvot = exp( -PhotonEnergy_I(iphot) / Te )                            
       photn3 = PhotonEnergy_I(iphot)**3                                   
                                                                        
       ! ...    loop over gas species                                          
                                                                        
       IMIXLOOP: do iMix = 1,nMix                                           
                                                                        
          brems(iMix) = 0.0                                            
          fotiza(iMix) = 0.0                                           
          fotize(iMix) = 0.0                                           
          
          ! ...       Bremsstrahlung                                              
          !           --------------                                              
          
          
                                                                        
          densi2 = Z2_I(iMix) * Concentration_I(iMix) * densnn * 1.e-16               
          
          ! ...       the free-free gaunt factor is a simple fit to the results   
          !           of Karzas and Latter (Ap. J. Suppl., 6, 167 (1961))         

          gam2lg = log10( 13.6 * Z2_I(iMix) / Te )                           
          gBrem  = 1. + 0.44 * exp( -0.25*(gam2lg+0.25)**2 )          
                                                                        
          brems(iMix) = 2.4e-21 * densi2 * gbrem * densne * &          
               (1.-exhvot) / ( sqrt( Te ) * photn3 )                   
                                                                        

          ! ...       photoionization                                             
          !           ----------------                                            
                                                                        
          ! ...       sum over ionization levels for the transition from          
          !           "iZ" to "iZ+1"                                      
          
          sum1a = 0.0                                                  
          sum1e = 0.0                                                  
          conpi = 1.66e-22 * densne / Te**1.5                         
                                                                        
          IZLOOP: do iZ = iZMin_I(iMix), min( iZMax_I(iMix), nZ_I(iMix) - 1)                               
                                                                       
                                                                        
             ! ...          find the principal quantum number of the valence electron
             !              in the ground state for the ion before ("nprin0");       
             !              "nbound" is the number of electrons bound to the ion     
             !              after another is captured (or, before ionization).       
             
             nBound = nZ_I(iMix) - iZ                            
             nGround = n_ground( iZ, nZ_I(iMix) )                                  
               
             ! ...          first, consider the contibution from valence shell       
             !              electrons                                                
               
             sum2a = 0.0                                               
             sum2e = 0.0
                                               
             if ( Concentration_I(iMix) * Population_II(iZ,iMix)< con(3) ) &
                  CYCLE IZLOOP     
                                                                        
             ! ...            sum over quantum states                                
             do  iN = nGround, nExcitation                       
                                                                        
                !  calculate the energy to excite the electron into     
                !  the continuum 
                                       
                eTransition = IonizPotential_II(iZ+1, iMix) - &
                     ExcitationEnergy_III(iN, iZ, iMix)         
                                                                        
                ! ...  the photon energy must exceed the binding energy     
               
                if ( PhotonEnergy_I(iphot) <  eTransition ) CYCLE            
                                                                        
                ! ... find the number of "screening" electrons "nscren"    
                !     and the number of electrons in the outermost shell   
                !     "nvalen"                                             
                
                if ( iN == nGround ) then                          
                   ! ...ground state ion                                   
                   nScreened = n_screened(iZ, nZ_I(iMix))                           
                   nValence = nBound - nScreened                           
                else                                                 
                   ! ...                ion is excited                                     
                   nScreened = nbound - 1                                
                   nValence = 1                                         
                endif
                                                                        
                ennpot = eTransition / Te                                   
                                                      
                                                                        
                ! ... use an "effective charge" seen by the electron,      
                !     corrected for screening                              
                iZEff = nZ_I(iMix) - nScreened                        
                                                                        
                deni = Population_II(iZ,iMix) * Partition_III(iN, iZ, iMix)     
                if ( ennpot .lt. 50. ) then                          
                   eqdeni = iN**2 * conpi * exp( ennpot ) *  &     
                        Population_II(iZ+1,iMix)*Partition_III(nGround,iZ+1,iMix)  
                else                                                 
                   eqdeni = 0.0                                       
                endif
                                                                        
                ! ...  the degeneracy level of the fully stripped ion       
                !      is 1, while a value of 2 is used for other ions.     
                if ( iZ +1== nZ_I(iMix) ) eqdeni = eqdeni * 2.    
                                                                        
                ! ... correction for stimulated emission; do not allow     
                !       this to be < 0 for the abs. coef.
                !

                !\
                ! Below lines for non_LTE are commented out  
                !/                  
                !non-LTE:if ( isw(6) .ne. 1 ) then                            
                   ! ...                 non-LTE correction                                
                !non-LTE:   dumden = deni                                     
                !non-LTE: else                                                 
                   ! ... LTE correction                                    
                   dumden = eqdeni                                   
                !non-LTE: endif

                stcorr = max ( 0., dumden - eqdeni * exhvot )        
                                                                        
                xsec = nValence * iZEff**4 /iN**5                        
                sum2a = sum2a + stcorr * xsec                        
                sum2e = sum2e + eqdeni*(1.-exhvot) * xsec            
                
             end do
                                                                                                                            
                                                                        
             ! ...          now, add core electron photoionization cross-sections    
             !              to the absorption term                                   
             
             ! ...          loop over inner shells (K,L,M,...); each inner shell     
             !              is assume to be full                                     
             
             !sum1ac = 0.0                                              
             ! nshels = nGround - 1                    
             if ( & !nshels > 0 .and. 
                  UseCoreElectron ) then
                call CON_stop('UseCoreElectron should be set to .false.')
                !  do ishell=1,nshels                                 
                
                ! ... determine the photoionization cutoff energy (eV); use
                !     the ionization potential of the outermost bound      
                !     electron; "nocc" is the number of electrons occupying
                !     shell "ishell", "nscren" is the number of electrons  
                !     screening shell "ishell"                             
                
                !    nocc = 2*ishell*ishell                               
                !    nscren = nscrsh( ishell )                            
                !    izeff = izgas(iMix) - nscren                         
                !    enpi = pot(iMix,izgas(iMix)+1-nscrsh(ishell+1))      
                !    if ( PhotonEnergy_I(iphot) .ge. enpi ) then                  
                !       sum1ac = sum1ac + nocc * izeff**4 / ishell**5      
                !    endif
                
                ! end do
                
                ! sum1a = sum1a + sum1ac * Population_II(iZ,iMix)             
                
             endif
             
             sum1a = sum1a + sum2a                                    
             sum1e = sum1e + sum2e                                    
             
          end do IZLOOP
          
          const = (1.99e-14*densnn) * Concentration_I(iMix) / photn3           
          fotiza(iMix) = const * sum1a                                
          fotize(iMix) = const * sum1e                                
       end do IMIXLOOP   !Over iMix                                                                    

                                                                        
! ...    Scattering contributions                                       
!        ------------------------                                       
                                                                        
! ...    Thomson scattering contribution                                
!        first, find the "effective" electron density for scattering;   
!        i.e., if the photon energy is greater than the binding energy  
!        of bound electrons, include bound electrons to the density.    
         scatne = 0.0                                                  
         do  iMix=1,nMix                                           
            iq = 0                                                      
            do  iZ = 0,nZ_I(iMix) -1                              
                                                      
               if ( PhotonEnergy_I(iphot) > IonizPotential_II(iZ+1,iMix)) then            
                  iq = iq + 1                                           
                  scatne = scatne + Concentration_I(iMix)
               else
                  EXIT
               endif                                                    
            end do
            !iQ is the last state at which the ionization potential 
            !is less than the photon energy
            !All the other state give only the free electron contribution to
            !
                                                                        
            do iZ= max(iZMin_I(iMix),iQ+1), iZMax_I(iMix)
               scatne = scatne + (iZ-iq) * Population_II(iZ,iMix)*&                    
                         Concentration_I(iMix) 
            end do
         end do
         scatne = scatne * densnn                                       
         tscatt = 6.66e-25 * scatne !Thomson cros-section                                     
                                                                        
         ! ...    contribution from plasma oscillations                          
         hnucut = sqrt( densne / 7.25e20 )      
                        
         if ( PhotonEnergy_I(iphot) .lt. hnucut ) then                          
            pscatt = 5.05e4 * sqrt( hnucut**2 - PhotonEnergy_I(iphot)**2 )      
         else                                                           
            pscatt = 0.                                                 
         endif                                                          
                                                                       
         if ( UseScattering) sctcfs(iphot) = tscatt + pscatt            
                                                                        
         !\
         !Sum up al the contributions                                                               
         abscfs(iphot) = 0.0                                            
         emscfs(iphot) = 0.0                                             
         if(UseBremsstrahlung)then                                         
            abscfs(iphot) = abscfs(iphot)+sum(brems(1:nMix)) 
            emscfs(iphot) = emscfs(iphot)+sum(brems(1:nMix))
         end if
         if(UsePhotoionization)then
            
            abscfs(iphot) = abscfs(iphot)+sum(fotiza(1:nMix)) 
            emscfs(iphot) = emscfs(iphot)+sum(fotize(1:nMix)) 
        end if
                                                      
                                                                        
           ! ...    line contributions                                             
           !        ------------------                                             

           ! ...    add in the contribution from bound-bound transitions           
                   
            call lines ( PhotonEnergy_I(iphot), abslns,emslns )      
                                                        
            abscfs(iphot) = abscfs(iphot) + abslns              
            emscfs(iphot) = emscfs(iphot) + emslns                                                             
                                                         
                                                                        
  end do  !over iphot                                                         
end subroutine abscon

  !====================
end module CRASH_ModMultiGroup

