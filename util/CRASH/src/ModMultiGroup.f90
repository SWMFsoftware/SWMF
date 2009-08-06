module CRASH_ModMultiGroup
  use CRASH_ModIomix
  use CRASH_ModOpacityVoigt, ONLY : linwid
  use CRASH_ModAtomicDataMix,ONLY : nMix, nZ_I, nExcitation
  use CRASH_ModAtomicDataMix,ONLY : Concentration_I
  use CRASH_ModAtomicDataMix,ONLY : IonizPotential_II
  use CRASH_ModAtomicMass,   ONLY : cAtomicMass_I
  use CRASH_ModStatSumMix,   ONLY : Population_II
  use CRASH_ModExcitationData, ONLY:n_ground
  use CRASH_ModExcitation,   ONLY : Partition_III, ExcitationEnergy_III
  implicit none
  SAVE
  PRIVATE !Except
  integer,parameter:: mxphot = 3000, mxgrps = 100, nfrqbb = 9
  !       mxphot  - photon energy mesh points                                    
  !       mxgrps  - opacity groups     
  !       nfrqbb  - number of photon energy points near a line center   
  !       at which absorption coefficients will be computed   



contains
  !=============================================
  subroutine bbneq0 ( TP,iZGas,nBound, enn,fnn,gnn )              

    ! ... calculates the effective energy, oscillator strength, and         
    !     gaunt factor for bound-bound transitions in which the             
    !     principal quantum number does not change.  This approximation     
    !     was taken from Post, et.al., Atomic & Nuclear Data Tables,       
    !     vol. 20, p.397 (1977).                                            

    ! ... input variables                                            
    !       tp      -  plasma temperature (eV)                              
    !       izgas   -  nuclear charge of the ion                            
    !       nbound  -  number of electrons bound to the ion

    real,intent(in) :: TP
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

    iZ =izgas - nbound                                            
    nb = min( nbound , 55 )                                           

    fnn = bbnton(nb,1) + bbnton(nb,2) / izgas                         
    enn = bbnton(nb,3) * (iZ+1)**bbnton(nb,4)                             

    if ( enn .le. 0. ) return                                         

    xnn = enn / tp                                                    
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
  subroutine meshhv ( TP,&
       DensNN,&
       DensNE,&
       EnGrup,&
       nGrups,&
       nptspg,&        
       PhotEn,&
       nphot )  
    !                                                                       
    ! ... set a mesh of photon energies at which 
    !     we would like to evaluate  
    !     absorption coefficients                                           
    !   


    ! ... input variables:                                                  
    !       tp      -  plasma temperature (eV)                              
    !       densnn  -  number density of all nuclei (cm**-3)                
    !       densne  -  electron density (cm**-3)                            
    !       engrup  -  photon energy group boundaries (ngrups+1) (eV)       
    !       ngrups  -  number of photon energy groups                       
    !       nptspg  -  minimum number of mesh points per energy group       
    !

    real, intent(in) :: TP, DensNN, DensNE
    integer,intent(in) :: nGrups, NptsPG !Inputs per group

    real, intent(in) :: EnGrup(nGrups+1)


    ! ... output variables:                                                 
    !       photen  -  array of photon energies (eV)                        
    !       nphot   -  number of mesh points including contrbutions from    
    !                  lines and photoionization edges                      
    !                                                                       
    integer,intent(out)::nPhot
    real,intent(out)::PhotEn(mxphot)

    ! ... set the maximum number of:                                        
    !       gases (mxgass), photen energy mesh points (mxphot),             
    !       temperatures (mxtemp), densities (mxdens),                      
    !       opacity groups (mxgrps), and gas atomic number (= mxatom-1)     

    !parameter ( mxtemp = 20, mxdens = 10 )                            
    !parameter ( mxphot = 3000, mxgrps = 100 )                         
    !parameter ( mxgass = 10, mxatom = 55 )                            


    !common / contrl / isw(30),con(30),iplot(30),dtheat,critsc         
    ! ................................................................      
    !                                                                  
    !common / gases / ngases,izgas(mxgass),atomwt(mxgass),             
    !                avgatw,avgatn,                                  
    !                fracsp(mxgass),fraciz(mxgass,mxatom),            
    !                fraclv(mxgass,mxatom,20),pot(mxgass,mxatom),     
    !                numocc(mxgass,mxatom)                            
    ! ................................................................

    !common / params / npqmax(mxgass), npmaxp, nfrqbb,                
    !                 npring(100), neopen(5,100), bbnton(6,55),       
    !                 defpot(27,55), noccdf(27,55)                    
    !................................................................     


    logical:: DoTest                                                  

    real,parameter:: dPlus =  1.001, dMinus = 0.999
    !\
    !Loop variables
    !/
    integer :: iGroup    !Enumerates energy groups
    integer :: iSubGrid  !Photon energy grid for each energy group
    integer :: nMax      !

    integer :: iMix !Over the mixture components
    integer :: iZ   !Ion charge
    integer :: iZPlus1 
    integer :: n    !Principal quantum number
    integer :: nP   !Principal quantum number for the excited state

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

    nphot = 0                                                        

    photen( 1:mxphot ) = 0.0


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
          photen(nphot+1) = exp( LogHv )                              
          nphot = nphot + 1                                           
       end do

    end do

    ! ... add 2 points near the plasma cutoff frequency                     

    hvcut = sqrt( densne / 7.25e20 )                                  
    if ( hvcut.gt.engrup(1) .and. hvcut.lt.engrup(ngrups+1) ) then    
       photen(nphot+1) = hvcut * dminus                               
       photen(nphot+2) = hvcut * dplus                                
       nphot = nphot + 2                                              
    endif

    ! ... add points near line centers and ionization edges                 

    do iMix=1,nMix                                              
       if ( Concentration_I(iMix) .lt. con(2) ) CYCLE                     

       do  iZ=0,nZ_I(iMix)-1                                  
          iZPlus1 = iZ + 1                                           
          if ( Concentration_I(iMix)*Population_II(iZ,iMix) .lt. con(3) )CYCLE

          ! ...       find the principal quantum number 
          !           of the valence electrons  
          !           in their ground state                                       
          nbound = nZ_I(iMix) - iZ                               
          nGround = n_ground(iZ,nZ_I(iMix))                                     

          do n = nGround, nExcitation                                
             if (.not. (Concentration_I(iMix) * &
                  Population_II(iZ,iMix) * &
                  Partition_III(n, iZ, iMix)& 
                  .lt. con(4) ).and.&
                  n.lt.nExcitation) then         

                do  np=n,nExcitation                              
                   ! ... calculate the energy of the transition             
                   if ( n .eq. np ) then                               
                      call bbneq0 ( TP,nZ_I(iMix),nbound, &            
                           ennp,dum1,dum2 )  
                      if ( ennp.le.0) CYCLE
                   else                                                
                      ennp = nGround**2 * IonizPotential_II(iZPlus1,iMix) * &             
                           (1./n**2-1./np**2)                        
                   endif
                   if ( ennp .gt. engrup(1) .and. &                     
                        ennp .lt. engrup(ngrups+1) ) then              

                      call linwid ( tp,densnn,ennp,cAtomicMass_I(nZ_I(iMix)), &      
                           gamma,avoigt,dnudop ) 

                      ! ...set points at and near the line center           
                      dhnu = con(9) * 4.14e-15 * gamma / 12.57         
                      photen(nphot+1) = ennp                           
                      nphot = nphot + 1                                

                      kk = 1.0
                      do  iPointPerLine = 1, nfrqbb / 2 - 1                                                                  
                         photen(nphot+1) = ennp - kk * dhnu             
                         photen(nphot+2) = ennp + kk * dhnu             
                         nphot = nphot + 2                              
                         kk = kk * 2.0
                      end do

                      ! ... stop if too many mesh pts are requested          
                      if ( nphot .gt. mxphot - 5 - nfrqbb )&           
                           stop ' too many pts in -meshhv-'              
                   endif
                end do

             endif

             ! ... add 2 points near the ionization edge                    
             edge = nGround**2 * IonizPotential_II(iZPlus1,iMix) / n**2                 
             if ( edge.gt.engrup(1) .and. edge.lt.engrup(ngrups+1) )&  
                  then                                                  
                photen(nphot+1) = edge * dminus                       
                photen(nphot+2) = edge * dplus                        
                nphot = nphot + 2                                     
             endif

          end do   !Over n 

       end do    !Over iZ

       ! ... add 2 points near photoionization edges of core electrons      
       do ncore=1,nZ_I(iMix)-1                                   
          edge = IonizPotential_II(nZ_I(iMix)+1-ncore,iMix)                        
          if ( edge.gt.engrup(1) .and. &                              
               edge.lt.engrup(ngrups+1) ) then                        
             photen(nphot+1) = edge * dminus                          
             photen(nphot+2) = edge * dplus                           
             nphot = nphot + 2                                        
          endif
       end do
    end do

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
      do  i=1,nPhot                                                    
         if ( PhotEn(i)<= 0.0 ) CYCLE                                     
         j = j + 1                                                   
         PhotEn(j) = PhotEn(i)                                             
      end do
      nPhot = j                                                          
      DoRepeat = .true.                                                                
      do while(DoRepeat)
         DoRepeat = .false.
         do i = 2,nPhot                                                 
            iMinus1 = i - 1                                                 
            if ( PhotEn(iMinus1) .gt. PhotEn(i) ) then                            
               rSave = PhotEn(iMinus1)                                          
               PhotEn(iMinus1) = PhotEn(i)                                        
               PhotEn(i) = rSave                                            
               DoRepeat = .true.                                           
            endif
         end do
         nCycle = nCycle + 1                                           

         ! ...    check for error                                                

         if ( ncycle .gt. nPhot ) call CON_stop( ' error in -sort-')                                                                                                  
         ! ...    if any of the elements have been swapped, try again
      end do
    end subroutine sort
    !===============
  end subroutine meshhv
  !====================
end module CRASH_ModMultiGroup
