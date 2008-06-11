!^CFG COPYRIGHT UM
! In the initial CRASH treatment of materials,
!
! 1. We can treat any given material as having a single average Z.
!
! 2. Mixtures can be treated either as an average material or as mixtures.
!
! 3. If mixtures are treated as mixtures, then collisional rates should be 
!    calculated using the "effective Z", which is the average of Z squared 
!    divided by the average Z.
!
! 4. The average Z can be determined from the Saha equation, but must not 
!    exceed the nuclear charge of the material in question.
!
! 5. We do not need to account for electron degeneracy in the initial model.
!
! 6. In our regime of interest, the electrons behave as an ideal gas in an 
!    ion-sphere environment within which Coulomb interactions do affect the 
!    electron pressure and internal energy. The electron pressure and internal 
!    energy are best calculated using equations 3.47 through 3.50 in 
!    R. P. Drake, High Energy Density Physics
!
! 7. The ion pressure is the ideal gas pressure. The ion internal energy 
!    includes the particle energy of random motion and the energy of ionization. 
!    The model of eqs. 3.74 through 3.76 in the mentioned book is acceptable. 
!    Alternatively, a more complex model using actual ionization energies would 
!    be acceptable.
!
! 8. The materials that matter are
!    - Beryllium
!    - Xenon
!    - Polyimide (C_22 H_10 N_2 O_5)
!
!***********************************************************************
!    calculation of ionization equilibrium, single material
!    for given - concentration of  atomic particles Co+Ciz [cm^{-3}]
!              - electron temperature Te[eV]
!***********************************************************************

module ModSaha
  use ModConst
  use ModStatSum  ! separated computations of StatSum

  implicit NONE
  include 'CRASH_definitions.h'
  real,  parameter ::  &
       cToleranceHere =  cTiny**2 ,     &  !012, 009 givs up to 4-5 digits
       cTwoPiInv = cOne /cTwoPi   ,     &
       Hmax  = cPlankH  /cErg     ,     &  !   erg*s Planck
       HmaxT = cPlankHBar                  !   the same/2Pi

  PUBLIC  proGON;  



  PRIVATE
  real:: &  
       vNe        , &   ! [cm^{-3}]
       vTe        , &   ! [eV]
       vNatomII   , &   ! a sum over concentrations of atoms+ions, cm^{-3}
       Z          , &   ! vNe/sum(C_I) <= only ions  
       C_I(0:99)  , &   ! C_I[0]-neutrals resided after ionization
       U_I(1:99)  , &   ! U_I(iZ) - the energy to make an ion iZ+ from (iZ-1)+    
       CUITotal         ! sum(Uiz*Ciz)/Ce  <~> ionization energy per 1 electron


  integer :: nCi = 54   ! for Xenon, must be a var for other elements

! end module ModSaha variables
!=============================================================================
 

  contains

 
 !=========================!  IONIZation , based on ModStatSaha, only
   subroutine proGON       !  DO the whole cycle for IONIZ
   real, parameter  :: &
                  Geo = 6.06e21 , &  ! cm-3 ev-3/2
                  Nao = 1.00e19 , &  ! 25 , &  ! 25, &     ! cm-3
                  Zo  = 1.00

   real :: Ziter = Zo, ZoLd
   real :: Ge  , Na , Ne,  nn, dnn, dzz,  NeInv
   real :: T32 , dTe
   real :: xx1 ,xx2 ,xx3
  
   integer :: iT, iter 

     dTe = 50.0 
     call set_element( 54 )
     call mod_init


  razTe:    do  iT = 1,7
           
           Na  = Nao
           vTe = dTe * (iT-1) +5.
           T32 = vTe * sqrt(vTe)
            Ge =  Geo* T32
 
           dzz =  1.0
! assigned above Ziter =  Zo
! Ziter would be saved between steps with different {Na,Te} values 
! and will used as a
! good initial value 
           iter = 0
      
       write (*,*) " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       write (*,"( a,f6.1,a,1pg12.4)") " ~~~~~~ Te=", vTe, " No=", Na



    razZ:  do while ( (abs(dzz) >=  cToleranceHere ) .AND. (iter<22)     )       
      iter = iter+1
      vNe =  Na * Ziter

      xx1 =  1.00/vTe          
      xx2 =  Ge /vNe ;  xx2 =  log(xx2) 

      write(*,"(a,f6.1,a,f6.3,a,1pg12.4)") "_gonStatS__", vTe," eStatW/Ne=",xx1," Ne:",vNe
      call  set_population( xx1 , xx2  )  ! TeInv, logGe)
        write(*,*) "_gon__<>Stat:", "iZmin=",iZmin, " iZdom=", iZdominant, " iZmax=",iZmax
      xx1 = z_averaged() ;       xx2 = z2_averaged()/xx1
        write(*,*) "_gon__<<Stat:", " <Z>=",xx1, " <Z^2>/<Z>=", xx2

      Zold = Ziter
      Ziter= zNew(Zold)
       
       dzz = (Zold-Ziter);!/(Zold +Ziter)
         write(*,*)'iTer=',Iter, " dZz=", dzz

     end do razZ

         write(*,*)'============================================='

     call set_ionization_equilibrium(Na*1.0e6, vTe*11610.)
     write(*,*)'Pavels <z>, <z^2>/z=', z_averaged(), z2_averaged()/z_averaged()
   end do razTe
 contains
  !=======================================!
  real function zNew( zOLD)! return (zITER)
  
   real, intent(in) :: zOLD
   real           :: z1, z2, zITER 
  
   z1    = z_averaged()
   z2    = z2_averaged()
   zITER = zOLD - ( zOLD - z1)/( cOne +( z2 -z1*z1)/zOLD )

   !vNe   = NatomII * zIter
   zNEW  = zITER
  end function zNew
  !=========================================
 end subroutine proGON

   

   !=======================================================================
   !=======================================================================
end module ModSaha
