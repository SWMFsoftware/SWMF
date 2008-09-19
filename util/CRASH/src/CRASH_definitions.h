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
         integer,parameter:: Xe_=0      !Xenon
         integer,parameter:: Be_=1      !Beryllium
	 integer,parameter:: Plastic_=2 ! Chemical formula:C_22 H_10 N_2 O_5
         integer,parameter:: C_=3, H_=4, N_=5, O_=6 !Composition elements
         integer,parameter, dimension(0:2) :: nZ_I=(/54, 4, 6 /)
