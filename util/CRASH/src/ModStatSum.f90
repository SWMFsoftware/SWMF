module ModStatSum
   use ModIonizPotential
   implicit none  
   integer,parameter::nZMax=54 !Atomic number of element in question
   integer:: nZ,iZMin,iZMax
   real,dimension(nZMax)::IonizPotential_I
   real,dimension(nZMax)::Population_I

Contains
   !Set the element and its Ionization Potentials
   subroutine set_element(nZIn)
        integer,intent(in)::nZIn
        nZ=nZIn
        call get_ioniz_potential(nZ,IonizPotential_I(1:nZ))
   end subroutine set_element

   !Finding the populations of the ions
   subroutine set_population(TeInv,GeLog)
         real,intent(in)::TeInv !the inverse of the electron temperature
         real,intent(in)::GeLog ! Logarithm of the electron stat weight
                                !log(1/Ne\lambda^3)
         !The program sets iZMin and iZMax, such that the ionization state
         !population is negligibly small for iZ<iZMin or iZ>iZMax. It also
         !sets the ionization state populations, Population_I(iZMin:iZMax),
         !such that:
         !log(Population_I(I)/Population_I(I-1))=
         !=-IonizPotential_I(I)*TeInv + GeLog
         !and sum(Pupulation_I(iZMin:iZMax)=1
   end subroutine set_population

   !Calculating the Z average values
   real function z_averaged()
     integer::iLoop
     z_averaged=0
     do iLoop=iZMin,iZMax
       z_averaged = z_averaged + Population_I(iLoop)*iLoop
     end do
   end function z_averaged

   !Calculating the Z^2 average values
   real function z2_averaged()
     integer::iLoop
     z2_averaged=0
     do iLoop=iZMin,iZMax
       z2_averaged = z2_averaged + Population_I(iLoop) * real(iLoop)**2
     end do
   end function z2_averaged

end module ModStatSum			 
