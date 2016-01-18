!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModInductiveE

  ! Calculate inductive electric field Eind 
  ! from the total electric field Efield
  ! minus the potential electric field Epot.
  ! Epot = -grad Potential, and Potential satisfies the Poisson equation
  !
  ! Laplace Potential = - div(E)
  
  implicit none
  SAVE

  private ! except

  public:: calc_inductive_e

  ! Total and inductive electric fields
  real, public, allocatable:: Efield_DGB(:,:,:,:,:)
  real, public, allocatable:: Eind_DGB(:,:,:,:,:)

  ! Electric potential
  real, public, allocatable:: Potential_GB(:,:,:,:)

  ! local variables

  contains
    subroutine calc_inductive_e

      call get_electric_field(Efield_DGB)
      call get_div_e

      ! Calculate Potential_GB from the Poisson equation 
      ! fill ghost cells at the end
      call calc_potential

      ! Calculate the inductive E field
      ! fill ghost cells at the end
      call calc_inductive_field

    end subroutine calc_inductive_e
    !=========================================================================
    subroutine get_div_e

      rInside = rBody - 0.4

      ! Fill in ghost cells
      call message_pass_cell(3, Efield_DGB)

      ! Fill in inner boundary cells
      do iBlock = 1, nBlock
         if(Unused_B(iBlock)) CYCLE
         if(body_BLK(iBLOCK))then
            do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
               if(r_BLK(i,j,k,iBlock) > rBody) CYCLE
               if(r_BLK(i,j,k,iBlock) < rInside) CYCLE
               call get_ie_potential(Xyz_DGB(:,i,j,k,iBlock), &
                    Potential_GB(i,j,k,iBlock))
            end do; end do; end do
         end do

         if(far_field_bcs_blk(iBlock))then
            call set_cell_boundary(nG, iBlock, 3, 

      ! Fill outer ghost cells with floating values
      do iBlock = 1, nBlock
         
         
      end do

      if(IsCartesianGrid)then
         
      else
      endif

    end subroutine get_div_e

    subroutine
      call message_pass_cell(Potential_GB)
    end module ModInductiveE

end module ModInductiveE

