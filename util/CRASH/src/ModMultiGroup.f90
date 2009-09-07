
  module CRASH_ModMultiGroup
  use CRASH_ModIonMix
  use CRASH_ModOpacityVoigt, ONLY : line_width, voigt_profile, UseVoigt
  use CRASH_ModAtomicDataMix,ONLY : nMix, nZ_I, nExcitation, nMixMax
  use CRASH_ModAtomicDataMix,ONLY : Concentration_I !(1:nMixMax)
  use CRASH_ModAtomicDataMix,ONLY : IonizPotential_II !(1:nZMax,1:nMixMax)
  use CRASH_ModAtomicMass,   ONLY : cAtomicMass_I !(1:nZMax)
  use CRASH_ModPartition,   ONLY : Population_II !(0:nZMax,1:nMixMax)
  use CRASH_ModExcitationData, ONLY:n_ground !(iZ,nZ)
  use CRASH_ModExcitationData, ONLY:n_screened !(iZ,nZ)
  use CRASH_ModExcitation,   ONLY : Partition_III !(nExcitation,0:nZMax  ,nMixMax)
  use CRASH_ModExcitation,   ONLY : ExcitationEnergy_III !nExcitation,0:nZMax-1,nMixMax)
  use CRASH_ModPartition,   ONLY : Na, Te, zAv
  use CRASH_ModPartition,   ONLY : iZMin_I !(1:nMixMax)
  use CRASH_ModPartition,   ONLY : iZMax_I !(1:nMixMax)
  use CRASH_ModFermiGas,    ONLY : LogGe
  use ModConst
  implicit none
  SAVE
  PRIVATE !Except


  !Public members
  public :: meshhv !Creates the grid of photon energies
  public :: lines  !Calculates the absorption and emission in lines
  public :: abscon !Calculates the absorption, emission, and scattering 
  public :: opacys !Calculates opacities
  public :: nGroup, OpacityPlanck_I, OpacityRosseland_I

  !For test:
  public :: PhotonEnergy_I, AbsorptionCoefficient_I, nPhoton, set_multigroup,EnergyGroup_I


  !       nPhotonMax  - photon energy mesh points                                    
  !       nGroupMax  - opacity groups     
  !       nfrqbb  - number of photon energy points near a line center   
  !       at which absorption coefficients will be computed  

  integer,parameter:: nPhotonMax = 10000, nGroupMax = 100, nfrqbb = 9

  !Radiation frequency groups:
  integer :: nGroup = nGroupMax 
  real :: EnergyGroup_I(0:nGroupMax)
  real :: DeltaLogFrequency

  public:: get_energy_g_from_temperature, get_temperature_from_energy_g


  integer,parameter::nptspg = 30

  ! Meshs to evaluate the absorption coefficients
  integer::nPhoton
  real ::PhotonEnergy_I(nPhotonMax)

  ! The absorption coefficients
  !       AbsorptionCoefficient_I  -  array of absorption coefficients (cm**-1)            
  !       emscfs  -  array of emission coefficients (cm**-1)         
  !       ScatteringCoefficient_I  -  array of scattering coefficients (cm**-1) 


  real,dimension(nPhotonMax) ::AbsorptionCoefficient_I,ScatteringCoefficient_I

  !\
  !The opacities averaged over photon groups
  !/
  real :: OpacityPlanck_I(nGroupMax),OpacityRosseland_I(nGroupMax) 

  !\
  !the same, averaged over the whole frequency range
  !/
  real :: OpacityPlanckTotal, OpacityRosselandTotal  


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

  logical,public :: UseBremsstrahlung = .true.
  logical,public :: UsePhotoionization = .true.
  logical,public :: UseScattering      = .false.
contains
  !======================================================================
  subroutine get_energy_g_from_temperature(iGroup, TgSI, EgSI, CgSI)
    !\
    !Input parameters
    !/
    integer,intent(in):: iGroup
    real,   intent(in):: TgSI    !Group temperature [K]

    !\
    !Output parameters
    !/
    real, optional, intent(out) :: EgSI    !Radiation energy per group, J/m3
    real, optional, intent(out) :: CgSI    !Radiation specific heat per group, J/(K.m3)
    real :: xMin, xMax
    !---------------------------------
    xMin = EnergyGroup_I(iGroup - 1)/(TgSI * cKToEV)
    xMax = EnergyGroup_I(iGroup    )/(TgSI * cKToEV)
    
    if(present(EgSI))EgSI = cNormG5 * gint(5,xMin,xMax) * (      cRadiation * TgSI**4)
    if(present(CgSI))CgSI = cNormG6 * gint(6,xMin,xMax) * (4.0 * cRadiation * TgSI**3)
    
  end subroutine get_energy_g_from_temperature
  !======================================================================
  subroutine get_temperature_from_energy_g(iGroup, EgSI, TgSIOut, CgSIOut)
    !\
    !Input parameters
    !/
    integer,intent(in):: iGroup
    real,   intent(in):: EgSI    !Group temperature [K]

    !\
    !Output parameters
    !/
    real, optional, intent(out) :: TgSIOut    !Radiation energy per group, J/m3
    real, optional, intent(out) :: CgSIOut    !Radiation specific heat per group, J/(K.m3)


    real :: xMin, xMax, FreqMin, FreqMax
    real :: TgSI, CgSI, ToleranceEg, DeltaEg

    real, parameter:: cTolerance = 1.0E-3
    
    integer, parameter :: nIter = 10
    integer :: iIter
    !--------------------------------------------------
    FreqMin = EnergyGroup_I(iGroup - 1) * cEVToK 
    FreqMax = EnergyGroup_I(iGroup    ) * cEVToK
    !\
    !Approximation to start:
    !/
    TgSI = sqrt(FreqMin*FreqMax)&
         /log(1.0 + cRadiation * FreqMin**2 * FreqMax**2 * cNormG5 * &
         DeltaLogFrequency / EgSI)
   
    
    ToleranceEg = cTolerance * EgSI

    iIter = 0
    DeltaEg = 2.0 * ToleranceEg !To start Newton-Rapson iterations
    do while (abs(DeltaEg) > ToleranceEg.and.iIter < nIter)
       xMin = FreqMin /TgSI
       xMax = FreqMax /TgSI

       DeltaEg = EgSI -  cNormG5 * gint(5, xMin, xMax) * cRadiation * TgSI**4
       CgSI    =         cNormG6 * gint(6,xMin,xMax) * (4.0 * cRadiation * TgSI**3) 
       TgSI = TgSI + DeltaEg/CgSI

       iIter = iIter + 1
    end do

    
    if(present(TgSIOut))TgSIOut = TgSI
    if(present(CgSIOut))CgSIOut = CgSI
  end subroutine get_temperature_from_energy_g
  !======================================================================
  subroutine set_multigroup(nGroupIn, FreqMinSI, FreqMaxSI)
    !\
    !Set the values of PHOTON ENERGY grid
    !/
    integer, intent(in) :: nGroupIn  !The number of photon energy groups

    real,    intent(in) :: FreqMinSI, FreqMaxSI  !Min and max FREQUENCIES [Hz]
    real:: elnmin,elnmax,elog
    integer:: iGroup
    !-----------------------
    nGroup = nGroupIn
    EnergyGroup_I(0) = FreqMinSI * cHPlanckEV !Photon energy

    EnergyGroup_I(nGroup) = FreqMaxSI * cHPlanckEV  !Photon energy

    if ( nGroup <=1 ) return                              
    elnmin = log( EnergyGroup_I(0) )                            
    elnmax = log( EnergyGroup_I(nGroup) )                    
    DeltaLogFrequency = ( elnmax-elnmin ) / nGroup                   
    elog = elnmin                                        
    do iGroup=1,nGroup-1                                  
       elog = elog + DeltaLogFrequency                               
       EnergyGroup_I(iGroup) = exp( elog )                         
    end do

  end subroutine set_multigroup
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
  subroutine lines (ephot,abstot)       

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

    real,intent(out) :: absTot   !!$ ,emsTot                            

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

    !The controbutions to the total absorption,
    !calculated by the 'abslin' subroutine
    real:: abscof  
    !-------------------


    iSav = 0                                                          
    abstot = 0.                                                       
                         

    DensNN = Na * 1.0e-6  !To convert to cm-3
    DensNe = DensNN * zAv !Electron density, in cm-3

    ! ... loop over gas species                                             
    IMIXLOOP: do iMix =1,nMix                                              
       if (Concentration_I(iMix) .lt. con(2) ) CYCLE IMIXLOOP

       ! ...    loop over ionization states                         
       IZLOOP: do  iZ=iZMin_I(iMix), min(iZMax_I(iMix), nZ_I(iMix) - 1)

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

                dnlqnp = denlq * Partition_III(iNUpper,iZ,iMix)                 
                call abslin   
                if ( abscof .ne. 0. ) then                            
                   abstot = abstot + abscof                           
                           
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
      !for which to calculate the absorption

      real:: DeltaNu 



      !The oscillator strength
      real :: OscillatorStrength

      !The line height at a given frequency (the line shape)
      real:: Shape

      !The absorption coefficient, to be corrected for stimulated emission
      real:: Alpha

      !\
      !Coefficients to account for stimulated emission
      !/
      real::corsea 
      real::corrse

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
         if ( ennp.le.0) then
            abscof = 0.0
            return
         end if
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
                                           

      if ( iN==iNUpper) then   
                                         
         call CON_stop('Inappropriate')

         alpha = 2.65e-2 * fnn * shape * denlqn * ( 1.-ex1 )         

                                                      
         corsea = 1.0                                            
         abscof = alpha * corsea                                     
                                     

      else                                                          

         OscillatorStrength = oscillator_strength ( iN, iNUpper) 

                         

         ! ...       correct for stimulated emission  

         corrse = 1.0 - ex1

         alpha  = 2.65e-2 * OscillatorStrength * shape * denlqn                    
         abscof = alpha * corrse   

        

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
    !       EnergyGroup_I  -  photon energy group boundaries (nGroup+1) (eV)       
    !       nGroup  -  number of photon energy groups                       
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
    do iGroup = 0,nGroup-1                                                

       hvmin = EnergyGroup_I(iGroup)                                             
       hvmax = EnergyGroup_I(iGroup+1)                                           
       dloghv = log( hvmax/hvmin ) / nptspg                           
       LogHv = log( hvmin ) - dloghv                                  
       if ( iGroup.lt.nGroup - 1 ) then                                       
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
    if ( hvcut.gt.EnergyGroup_I(0) .and. hvcut.lt.EnergyGroup_I(nGroup) ) then    
       PhotonEnergy_I(nPhoton+1) = hvcut * dminus                               
       PhotonEnergy_I(nPhoton+2) = hvcut * dplus                                
       nPhoton = nPhoton + 2                                              
    endif

    ! ... add points near line centers and ionization edges                 

    do iMix=1,nMix                                              
       if ( Concentration_I(iMix) .lt. con(2) ) CYCLE                     

       do  iZ=iZMin_I(iMix), min(iZMax_I(iMix), nZ_I(iMix) - 1)                                           
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
                   if ( ennp .gt. EnergyGroup_I(0) .and. &                     
                        ennp .lt. EnergyGroup_I(nGroup) ) then              

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

             if ( edge.gt.EnergyGroup_I(0) .and. edge.lt.EnergyGroup_I(nGroup) )&  
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
          if ( edge.gt.EnergyGroup_I(0) .and. &                              
               edge.lt.EnergyGroup_I(nGroup) ) then                        
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
    use CRASH_ModPartition, ONLY: Z2_I

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
    !       AbsorptionCoefficient_I  -  array of absorption coefficients (cm**-1)            
    !       emscfs  -  array of emission coefficients (cm**-1)         
    !       ScatteringCoefficient_I  -  array of scattering coefficients (cm**-1)          
    !                                                                      

    real:: DensNN, DensNe    ![cm -3]

    !\
    !Arrays to collect contributions from different effects for a given component
    real,dimension(nMixMax) :: BremsStrahlung_I = 0.0
    real,dimension(nMixMax) :: PhotoIonizationAbs_I = 0.0


    !\
    ! Loop variables
    !/
    integer :: iMix, & !runs over the mixture component
         iZ,   & !runs over the charge number
         iN,   & ! runs over the quntum principal number for the lower level
         nBound, &!Number of bound electrons
         nGround,&!For a given iZ, the principal number of the ground state.
         iPhoton

    real :: eTransition, ETransitionPerTe

    !\
    !Partial sums
    !/
    real:: SumOverZ, SumOverN 

    !\
    ! Bound-bound contributions
    !/
    real :: abslns

    !\
    ! Cutoff energy for the photon with the frequency correspondent to
    ! the plasma electron frequency
    !/
    real :: HNuCut

    !Misc: functions of the photon energy
    real :: ExpMinusHNuPerT,  PhotonEnergyCubed

    !Misc: Density times Z2
    real :: DensityZ2

    !Misc: coefficients used in free-free gaunt factor
    real :: Log10OfGamma2, GauntFactorBrems

    !Number of screened and valence electrons, as well as iZEff
    integer :: nScreened, nValence, iZEff, iQ

    ! densities to calculate photoionization
    real :: AbsorberDensity, eqdeni, dumden, ExpMuPerT

    !\
    !Coefficients to account for stimulated emission
    !/
    real :: stcorr

    !Misc: Scattering parameters
    real :: ScatNe, tScatt, pScatt

    !Misc: Photoionization parameters
    real :: PhotoIonizationConst
    
    !---------------------------------------------------

    !Initialization

    AbsorptionCoefficient_I  = 0.0  !-  array of absorption coefficients (cm**-1)                    
    ScatteringCoefficient_I  = 0.0  !-  array of scattering coefficients (cm**-1)     

    DensNN = Na * 1.0e-6
    DensNe = DensNN * zAv

    ! ... loop over photon energies
    do iPhoton = 1,nPhoton                                              
       !\
       !The dependence of the cross-sections \propto (Photon Energy)^{-3} is typical
       !Therefore, introduce:
       !/                           
       PhotonEnergyCubed = PhotonEnergy_I(iPhoton)**3

       !\
       !The factor needed to account for the stimulated emission:
       !/
       ExpMinusHNuPerT = exp( -PhotonEnergy_I(iPhoton) / Te )                                    

       ! ...    loop over gas species                                          

       IMIXLOOP: do iMix = 1,nMix                                           

          BremsStrahlung_I(iMix) = 0.0                                            
                                                  

          ! ...       Bremsstrahlung                                              
          !           --------------                                              



          DensityZ2 = Z2_I(iMix) * Concentration_I(iMix) * densnn * 1.e-16               

          ! ...       the free-free gaunt factor is a simple fit to the results   
          !           of Karzas and Latter (Ap. J. Suppl., 6, 167 (1961))         

          !
          Log10OfGamma2 = log10( 13.6 * Z2_I(iMix) / Te )                           
          GauntFactorBrems  = 1. + 0.44 * exp( -0.25*(Log10OfGamma2+0.25)**2 )          

          BremsStrahlung_I(iMix) = 2.4e-21 * DensityZ2 * GauntFactorBrems * densne * &          
               (1.-ExpMinusHNuPerT) / ( sqrt( Te ) * PhotonEnergyCubed )
          !!NOTE: there is a more accurate calculation of the Gaunt factor in 
          !ggff.f file from HYADES


          ! ...       photoionization                                             
          !           ----------------                                            

          ! ...       sum over ionization levels for the transition from          
          !           "iZ" to "iZ+1"                                      

          PhotoIonizationAbs_I(iMix) = 0.0                                            
          SumOverZ = 0.0                                                  
                                                  
          ExpMuPerT = 1.66e-22 * densne / Te**1.5                         

          IZLOOP: do iZ = iZMin_I(iMix), min( iZMax_I(iMix), nZ_I(iMix) - 1)                               


             ! ...          find the principal quantum number of the valence electron
             !              in the ground state for the ion before ("nprin0");       
             !              "nbound" is the number of electrons bound to the ion     
             !              after another is captured (or, before ionization).       

             nBound = nZ_I(iMix) - iZ                            
             nGround = n_ground( iZ, nZ_I(iMix) )                                  

             ! ...          first, consider the contibution from valence shell       
             !              electrons                                                

             SumOverN = 0.0                                               

             if ( Concentration_I(iMix) * Population_II(iZ,iMix)< con(3) ) &
                  CYCLE IZLOOP     

             ! ...            sum over quantum states                                
             do  iN = nGround, nExcitation                       

                !  calculate the energy to excite the electron into     
                !  the continuum 

                eTransition = IonizPotential_II(iZ+1, iMix) - &
                     ExcitationEnergy_III(iN, iZ, iMix)         

                ! ...  the photon energy must exceed the binding energy     

                if ( PhotonEnergy_I(iPhoton) <  eTransition ) CYCLE  
          

                ! ... find the number of "screening" electrons     
                !     and the number of electrons in the outermost shell   
                !                 

                if ( iN == nGround ) then                          
                   ! ...ground state ion                                   
                   nScreened = n_screened(iZ, nZ_I(iMix))                           
                   nValence = nBound - nScreened                           
                else                                                 
                   ! ...ion is excited                                     
                   nScreened = nbound - 1                                
                   nValence = 1                                         
                endif

                ETransitionPerTe = eTransition / Te                                   


                AbsorberDensity = Population_II(iZ,iMix) * Partition_III(iN, iZ, iMix)     
               

                                        
                !\
                ! The version more close to that implemented in
                ! Emilio Minguez et al. With this approach we substitute
                ! eTransition for the combination Ry * Z_eff^3/iN^2, which
                ! is implied in the ionmix.f.
                ! By this account, our factor PhotoIonizationConst
                ! differs from the version in the ionmix by a foctor of
                ! (1/13.6)**3 where 13.6 = cRyEV.
                ! The last multiplier accounts for effects of the Fermi
                ! statistics in an electron gas. The formula is available in
                ! PhotoIonization.pdf document. Note, that our LogGe = -\mu/T
                !/         
                SumOverN = SumOverN + nValence * iN/(real(iZ+1)**2) * AbsorberDensity&
                      * (1.-ExpMinusHNuPerT) * (eTransition**3)/&
                      (1.0 + ExpMinusHNuPerT*exp(ETransitionPerTe - LogGe))

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
                !    if ( PhotonEnergy_I(iPhoton) .ge. enpi ) then                  
                !       sum1ac = sum1ac + nocc * izeff**4 / ishell**5      
                !    endif

                ! end do

                ! sum1a = sum1a + sum1ac * Population_II(iZ,iMix)             

             endif

             SumOverZ = SumOverZ + SumOverN                                                                       

          end do IZLOOP
          !\
          ! Commented out version from the ionmix:
          ! const = (1.99e-14*densnn) * Concentration_I(iMix) / PhotonEnergyCubed
          ! Insetad we use here the expression from Zel'dovich and Raizer, Eq.5.34,
          ! which differs by a factor of (1/13.6)**3, es we explained above.
          ! The exact espression for the photoionization cross-section is:
          ! 7.9e-18 cm^2= \frac{64\pi}{3\sqrt{3}}\alpha a_0^2, where \alpha is the fine
          ! structure constant and a_0 is the Bohr radius
          !/
          PhotoIonizationConst = 7.9e-18 * Concentration_I(iMix) / PhotonEnergyCubed   
          PhotoIonizationAbs_I(iMix) = PhotoIonizationConst * SumOverZ
                 
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

             if ( PhotonEnergy_I(iPhoton) > IonizPotential_II(iZ+1,iMix)) then            
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

       if ( PhotonEnergy_I(iPhoton) .lt. hnucut ) then                          
          pscatt = 5.05e4 * sqrt( hnucut**2 - PhotonEnergy_I(iPhoton)**2 )      
       else                                                           
          pscatt = 0.                                                 
       endif

       if ( UseScattering) ScatteringCoefficient_I(iPhoton) = tscatt + pscatt            

       !\
       !Sum up al the contributions                                                               
       AbsorptionCoefficient_I(iPhoton) = 0.0                                                                                      
       if(UseBremsstrahlung)then                                         
          AbsorptionCoefficient_I(iPhoton) = AbsorptionCoefficient_I(iPhoton)+sum(BremsStrahlung_I(1:nMix)) 
       end if
       if(UsePhotoionization)then

          AbsorptionCoefficient_I(iPhoton) = AbsorptionCoefficient_I(iPhoton)+sum(PhotoIonizationAbs_I(1:nMix))  
       end if



       ! ...    line contributions                                             
       !        ------------------                                             

       ! ...    add in the contribution from bound-bound transitions           

       call lines ( PhotonEnergy_I(iPhoton), abslns)      

       AbsorptionCoefficient_I(iPhoton) = AbsorptionCoefficient_I(iPhoton) + abslns                                                            


    end do  !over iPhoton                                                         
  end subroutine abscon
  !====================



  subroutine opacys ( TRadIn )                    

    ! ... this routine computes the Planck and Rosseland opacities for      
    !     the photon energy group "EnergyGroup_I".  "opacpm" and "OpacityRosseland_Im" are the  
    !     Planck and Rosseland mean opacities (integrated over frequency).  


    real, optional, intent(in) :: TRadIn

    ! ... input variables:                                                  
    !       Te     - plasma temperature (eV)                                
    !       densnn - number density of all nuclei (cm**-3)                 
    !       densne - electron density (cm**-3)                              
    !       photen - photon energies (eV)                                   
    !       nphot  - number of elements in "photen" array                   
    !       AbsorptionCoefficient_I - array of absorption coefficients (cm**-1)              
    !       emscfs - array of emission coefficients (cm**-1)                
    !       ScatteringCoefficient_I - array of scattering coefficients (cm**-1)              
    !       EnergyGroup_I - photon energy group boundaries (eV) (nGroup+1)  
    !       nGroup - number of energy bins                                  
    !       ntrad  - number of radiation temperatures                       
    !       trad   - array of radiation temperatures (eV)                   


    ! ... output variables:                                                 
    !       OpacityPlanck_I - Planck group opacities for absorption (cm-1)        
    !       OpacityRosseland_I  - Rosseland group opacity (cm-1)                      
    !       OpacityPlanckTotal - Planck mean opacity for absorption (cm-1)                        
    !       OpacityRosselandTotal  - Rosseland mean opacity (cm-1)                       
    !       culrat - plasma cooling rate (erg cm**3/sec)                    
    !       op2tp  - Planck 2-temperature opacities (cm**2/g)               
    !       op2tr  - Rosseland 2-temperature opacities (cm**2/g)            
    !                                                                       

    real :: DensNE, DensNN
    real :: TRad
    !\
    !Loop Variables
    !/
    integer:: iGroup

    real :: XGroupMin, XGroupMax, LineCoreOpacity
    !---------------------
    DensNN = Na * 1.0e-6  !To convert to cm-3
    DensNe = DensNN * zAv !Electron density, in cm-3
    if(present(TRadIn))then
       TRad = TRadIn
    else
       TRad = Te
    end if



    ! ... compute the Planck and Rosseland opacities for each photon        
    !     energy group                                                      

    OpacityPlanckTotal = 0.0                                                                                                           
    OpacityRosselandTotal = 0.0                                                        
    do iGroup = 1, nGroup                                                 
       XGroupMin = EnergyGroup_I(iGroup-1) / TRad                                
       XGroupMax = EnergyGroup_I(iGroup  ) / TRad                                  

       call opacgp (OpacityPlanck_I(iGroup),OpacityRosseland_I(iGroup) )  
       if ( DoNotAddLineCore ) then                                     
          ! ...       use analytic solution to bound-bound opacities 
          call opacbb (LineCoreOpacity)  
          OpacityPlanck_I(iGroup) = OpacityPlanck_I(iGroup) + LineCoreOpacity                                                            
       endif

       OpacityPlanckTotal = OpacityPlanckTotal + gint(5,XGroupMin,XGroupMax) * OpacityPlanck_I(iGroup)                         
       OpacityRosselandTotal  = OpacityRosselandTotal  + gint(6,XGroupMin,XGroupMax) / OpacityRosseland_I(iGroup)               
    end do
    ! ... compute the total opacities based on the group opacities         
    OpacityPlanckTotal = OpacityPlanckTotal * cNormG5                                                                                 
    OpacityRosselandTotal = 1.0/ (OpacityRosselandTotal * cNormG6)                                            


  contains

    subroutine opacgp ( OpacityPlanck, OpacityRosseland ) 
      !... this routine computes the Planck and Rosseland opacities for
      !     the photon energy group from "XGroupMin" to "XGroupMax".
      !
      ! ... input variables:
      !       Te     - plasma temperature (eV)                         
      !       rho    - mass density (g/cm**3)                                 
      !       trad   - radiation temperature (eV)                             
      !       PhotonEnergy_I - photon energies at which the absorption coefficients   
      !                were calculated                                        
      !       nphot  - number of photon energy mesh points                    
      !       AbsorptionCoefficient_I - array of absorption coefficients (cm**-1) (w/o scatt)  
      !       emscfs - array of emission coefficients (cm**-1) (w/o scatt)    
      !       ScatteringCoefficient_I - array of scattering coefficients (cm**-1)             
      !       XGroupMin  - minimum photon energy (in units of kT)                
      !       XGroupMax  - maximum photon energy (in units of kT)   



      !                                                                      
      ! ... output variables:                                                 
      !       OpacityPlanck - Planck opacity for absorption (cm**-1)                                
      !       OpacityRosseland  - Rosseland opacity (cm**-1)                             
      !                                                                       

      real, intent(out) :: OpacityPlanck, OpacityRosseland
      !\
      !Loop variables
      !/
      integer :: iPhot, iPhot0

      !\
      !Partial sums
      !/
      real :: dsumpa, dsmpa0, dsumr, dsumr0

      !Misc:
      real :: xg1,xg2, emx, fpa2,fpa20,fr2, fr20,sumpa,sumpa0,sumr,sumr0
      real :: fpa1, fpa10,fr1, fr10, dxg, g5, g6
      !-----------------------------                                     


      ! ... find the index where the photon energy equals the lower group     
      !     boundary (note: mesh pts should be placed at each boundary)      
      do iPhot=1,nPhoton                                                   
         if ( PhotonEnergy_I(iPhot) >= XGroupMin * Trad * 0.99999 ) then
            iPhot0 = iPhot
            EXIT
         end if
      end do

      ! ... integrate to get the group opacities                              
      !     ------------------------------------                              

      ! ... initialize values for the first point                             

      xg2 = PhotonEnergy_I(iphot0) / trad                                       
      emx = exp( -xg2 )                                                 

      fpa2 = xg2**3 * emx * AbsorptionCoefficient_I(iphot0) / ( 1.-emx )                 
      fpa20 = AbsorptionCoefficient_I(iphot0)

      fr2 = xg2**4 * emx / (1.-emx)**2 / &
           ( AbsorptionCoefficient_I(iphot0) + ScatteringCoefficient_I(iphot0) )                         
      fr20 = 1. / ( AbsorptionCoefficient_I(iphot0) + ScatteringCoefficient_I(iphot0) )                   

      ! ... loop over photon energy mesh points                              

      sumpa  = 0.0                                                      
      sumpa0 = 0.0                                                      

      sumr   = 0.0                                                      
      sumr0  = 0.0                                                      

      do iPhot = iPhot0 + 1, nPhoton

         if ( PhotonEnergy_I(iPhot) > XGroupMax * trad  * 1.00001 ) EXIT

         xg1 = xg2                                                      
         xg2 = PhotonEnergy_I(iPhot) / TRad                                     

         fpa1  = fpa2                                                   
         fpa10 = fpa20                                                  

         fr1   = fr2                                                    
         fr10  = fr20                                                   

         if ( xg2 .lt. 1.e20 ) then                                     

            emx = exp( -xg2 )                                           

            fpa2  = xg2**3 * emx * AbsorptionCoefficient_I(iphot) / ( 1.-emx )           
            fpa20 = AbsorptionCoefficient_I(iphot)                                      

            fr2   = xg2**4 * emx / (1.-emx)**2 / &                       
                 ( AbsorptionCoefficient_I(iPhot) + ScatteringCoefficient_I(iPhot) )                   
            fr20  = 1. / ( AbsorptionCoefficient_I(iphot) + ScatteringCoefficient_I(iphot) )              

            dxg = xg2 - xg1                                             

            ! ...       integrate using a logarithmic interpolation scheme   

            ! ...       Planck mean absorption opacities                     
            if ( abs( fpa1-fpa2 ) .gt. 1.e-3*fpa2 ) then                
               if ( fpa2 .ne. 0. ) then                                  
                  dsumpa = dxg * ( fpa2-fpa1 ) / log( fpa2/fpa1 )         
               else                                                      
                  dsumpa = 0.                                             
               endif
            else                                                       
               dsumpa = dxg * fpa1                                       
            endif

            ! ...       mean (non-weighted) absorption opacities                 
            if ( abs( fpa10-fpa20 ) .gt. 1.e-3*fpa20 ) then             
               if ( fpa20 .ne. 0. ) then                                 
                  dsmpa0 = dxg * ( fpa20-fpa10 ) / log( fpa20/fpa10 )     
               else                                                      
                  dsmpa0 = 0.                                             
               endif
            else                                                        
               dsmpa0 = dxg * fpa10                                      
            endif

            !NON-LTE part 
            ! ...       Planck mean emission opacities                              
            !if ( abs( fpe1-fpe2 ) .gt. 1.e-3*fpe2 ) then                
            !  if ( fpe2 .ne. 0. ) then                                  
            !    dsumpe = dxg * ( fpe2-fpe1 ) / log( fpe2/fpe1 )         
            !  else                                                      
            !    dsumpe = 0.                                             
            !  endif                                                     
            !else                                                        
            !  dsumpe = dxg * fpe1                                      
            !endif                                                      

            ! ...       Rosseland mean opacities                                 
            if ( abs( fr1-fr2 ) .gt. 1.e-3*fr2 ) then                   
               if ( fr2 .ne. 0. ) then                                   
                  dsumr = dxg * ( fr2-fr1 ) / log( fr2/fr1 )              
               else                                                      
                  dsumr = 0.                                              
               endif
            else                                                        
               dsumr = dxg * fr1                                         
            endif

            ! ...       mean (non-weighted) "transport" opacities                   
            if ( abs( fr10-fr20 ) .gt. 1.e-3*fr20 ) then                
               if ( fr20 .ne. 0. ) then                                  
                  dsumr0 = dxg * ( fr20-fr10 ) / log( fr20/fr10 )        
               else                                                      
                  dsumr0 = 0.                                             
               endif
            else                                                        
               dsumr0 = dxg * fr10                                       
            endif

            sumpa  = sumpa + dsumpa                                    
            sumpa0 = sumpa0 + dsmpa0                                    

            sumr   = sumr + dsumr                                       
            sumr0  = sumr0 + dsumr0                                     

         endif
      end do
      ! ... normalize to get the group opacities                  

      g5 = gint(5,XGroupMin,XGroupMax)                                          
      g6 = gint(6,XGroupMin,XGroupMax)                                          
      if ( XGroupMax.ge.1./con(7) .and. XGroupMin.le.con(7) .and. &              
           g5.gt.0. ) then                                              
         ! ...    weight absorption coef. by Planck function          
         OpacityPlanck = sumpa  / g5                                                                            
      else                                                             
         ! ...    use straight average of absorption coef.                       
         OpacityPlanck = sumpa0  / ( XGroupMax-XGroupMin )                        
      endif

      if ( XGroupMax.ge.1./con(7) .and. XGroupMin.le.con(7) .and. &              
           sumr.gt.0. .and. g6.gt.0. ) then                             
         ! ...    weight absorption and scattering coefs. by Planck function 
         !        derivative wrt temperature                                     
         OpacityRosseland = g6  / sumr                                        
      else if ( sumr0.gt.0. ) then                                      
         ! ...    use straight average of absorption plus scattering coefs.      
         OpacityRosseland = ( XGroupMax-XGroupMin )  / sumr0                          
      else                                                              
         ! ...    in the very high PhotonEnergy_I energy limit, this just becomes equal 
         !        to Thomson scattering contribution (assuming the absorption  
         !        is much less)                                            
         OpacityRosseland = ( ScatteringCoefficient_I(nPhoton) + AbsorptionCoefficient_I(nPhoton) )                 
      endif


    end subroutine opacgp
    !=============================
    subroutine opacbb (OpacityPlanck)  
      !                                                                       
      ! ... compute the contribution to the group opacity from                
      !     all lines (bound-bound transitions)                               
      !                                                                       
      ! ... input variables:                                             
      !       Te      -  plasma temperature (eV)                    
      !       densnn  -  number density of all nuclei (cm**-3)       
      !       densne  -  electron density (cm**-3)                
      !       trad    -  radiation temperature (eV)                  
      !       XGroupMin   -  minimum photon energy of group (in units of kT) 
      !       XGroupMax   -  maximum photon energy of group (in units of kT) 



      !                                                                
      ! ... output variables:                                  
      !       OpacityPlanck   -  Planck absorption due to lines (cm^-1)                   
      !       OpacityRosseland  -  Rosseland opacity due to lines (cm^)            
      !---------------------------------------------                                                                

      real, intent(out) :: OpacityPlanck           


      !\
      ! Loop variables
      !/
      integer :: iMix, & !runs over the mixture component
           iZ,   & !runs over the charge number
           iN,   & ! runs over the quntum principal number for the lower level
           iNUpper,& !the same, for upper level 
           nBound, &!Number of bound electrons
           nGround  !For a given iZ, the principal number of the ground state.

      real :: TransitionEnergy, Energy2TRadRatio 
      real :: ExpOfEnergy2TRadRatio, ExpOfEnergy2TeRatio
      !Misc:
      real :: g5, const, const1, opacij, fnnp


      !----------------------------


      OpacityPlanck = 0.0                                                                                                           

      g5  = gint(5,XGroupMin, XGroupMax )                                      

      if ( g5 .eq. 0.0 ) return                                         

      const = 1.10e-16 * DensNN / ( Te * g5 )                     

      ! ... loop over gas species                                             
      IMIXLOOP: do iMix = 1,nMix                                              
         if ( Concentration_I(iMix) .lt. con(2) ) CYCLE IMIXLOOP

         ! ...    loop over ionization states                                  
         IZLOOP:   do iZ = iZMin_I(iMix), min( iZMax_I(iMix), nZ_I(iMix) - 1)                                                                              
            if ( Concentration_I(iMix) * Population_II(iZ,iMix) .lt. con(3)) &
                 CYCLE IZLOOP  

            ! ...       find the principal quantum number of the valence electrons  
            !           in their ground state                                                                            

            const1 = const * Concentration_I(iMix) * Population_II(iZ,iMix)       


            ! ...       loop over initial quantum states 

            INITIAL: do iN = n_ground(iZ, nZ_I(iMix)) , nExcitation -1
               if ( Concentration_I(iMix) * Population_II(iZ,iMix) * Partition_III(iN, iZ, iMix)&
                    .lt. con(4) ) CYCLE INITIAL                             

               ! ... calculate opacity due to "delta n" = 0 transitions       

               !call bbneq0 ( Te,izgas(lgas),nbound,&
               !     enn,fnn,gnn )        
               !if ( isw(11).ne.0 ) enn = 0.                             
               !Energy2TRadRatio = enn / trad                                          
               !if ( Energy2TRadRatio.gt.XGroupMin .and. Energy2TRadRatio.le.XGroupMax ) then                
               !   ExpOfEnergy2TRadRatio = exp( -Energy2TRadRatio )                                      
               !   opacnn = const1 * fnn * fraclv(lgas,izp1,n) *&
               !        ExpOfEnergy2TRadRatio * Energy2TRadRatio**3
               ! ...  correct for the "effective" stimulated emission to    
               !      give the proper form for the cooling rate             
               !xnjoni = gnn * densne * ExpOfEnergy2TRadRatio / &
               !   ( gnn * densne + enn**3 * sqrt(Te) * 2.74e12 ) 
               !   corsee = xnjoni / ExpOfEnergy2TRadRatio                                 
               !   corsea = ( 1.-xnjoni ) / ( 1.-ExpOfEnergy2TRadRatio )                   
               !   OpacityPlanck = OpacityPlanck + opacnn * corsea                       
               !   opems = opems + opacnn * corsee                      
               !endif
               ! ...          loop over final quantum states

               FINAL: do iNUpper = iN+1, nExcitation                               

                  ! ...            compute the transition energy     
                  TransitionEnergy = ExcitationEnergy_III(iNUpper, iZ,iMix)- &
                       ExcitationEnergy_III(iN, iZ,iMix)
                  Energy2TRadRatio = TransitionEnergy / TRad                                       
                  if ( Energy2TRadRatio .le. XGroupMin .or. Energy2TRadRatio .gt. XGroupMax ) CYCLE FINAL    

                  ExpOfEnergy2TRadRatio = exp( -Energy2TRadRatio )                                      
                  fnnp = oscillator_strength(iN, iNUpper)

                  opacij = const1 * fnnp * Partition_III(iN, iZ, iMix)  * &
                       ExpOfEnergy2TRadRatio * Energy2TRadRatio**3                                  

                  !\
                  !Account for the stimulated emission
                  !/
                  ! ...            note that in LTE, "ExpOfEnergy2TeRatio" = "ExpOfEEnergy2TRadRatio0"      

                  ExpOfEnergy2TeRatio = iN * iN * Partition_III(iNUpper, iZ, iMix)  / &                 
                       ( iNUpper * iNUpper * Partition_III(iN, iZ, iMix)  )                    


                  OpacityPlanck = OpacityPlanck  +  &
                       opacij * ( 1.-ExpOfEnergy2TeRatio ) / ( 1.-ExpOfEnergy2TRadRatio )

               end do FINAL
            end do INITIAL
         end do IZLOOP
      end do IMIXLOOP
    end subroutine opacbb
  end subroutine opacys
  !====================
end module CRASH_ModMultiGroup

