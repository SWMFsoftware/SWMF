   SUBROUTINE RCM_Compute_ionospheric_potential
      USE rcm_variables
      implicit none

      ! a separate RCM routine to handle solving for
      ! ionospheric potential including setting the
      ! boundary condition (at high-L boundary)

      INTEGER :: j, info_bicgstab
    
  INTERFACE
     SUBROUTINE Gmresm_Matvec (x, y, n)
       USE rcm_variables, ONLY : iprec, rprec
       IMPLICIT NONE
       INTEGER (iprec), INTENT (IN) :: n
       REAL (rprec), INTENT (IN) :: x(n)
       REAL (rprec), INTENT (OUT) :: y(n)
     END SUBROUTINE Gmresm_Matvec
  END INTERFACE


      IF (ipot == 4) then    ! GMRES

         !  Note: new_cfive thinks c5w is the old c5w without d
         !  denominator. Before activating winds, check this.  5/29/99

         if (IsPartofFramework) then

            DO j = 1-n_gc, jsize+n_gc
               vbnd(j) = v(imin_j(j),j)
            END DO

         else

            CALL Get_potential_on_bndy ()

         end if

         CALL Comput_c5_wind
         CALL New_cfive
         CALL New_coeff
         CALL Comput_lowlat_boundary
         CALL Define_pde_matrix ()
         CALL Gmresm (Gmresm_Matvec, nij_pde, x0_pde, b_mtrx_pde, &
              tol_gmres, iter_max_gmres)
         CALL Gmresm_unwrap_pde_solution ()
         !
      ELSE IF(ipot == 5) then    ! BICGSTAB
         !
         CALL Comput_c5_wind
         CALL New_cfive
         CALL New_coeff
         CALL Comput_lowlat_boundary
         CALL Define_pde_matrix ()
         iter_bicgstab = iter_max_bicgstab
         tol_bicgstab = 3.0E-4
         CALL Impl_bicgstab (Gmresm_Matvec, b_mtrx_pde, x0_pde, SIZE(x0_pde), &
                tol_bicgstab, 'rel', iter_bicgstab, info_bicgstab)
         CALL Gmresm_unwrap_pde_solution ()
         if (info_bicgstab /= 0 ) then
            print*,info_bicgstab
            call CON_stop('IM_ERROR: bicgstab failed')
         end if
         !
      ELSE IF (ipot == 6) THEN
         !
         ! Potential should have been set elsewhere,
         ! here we assume that array V holds the necesssary
         ! values.
         !
      ELSE
         !
         call CON_STOP('IM_ERROR: invalid value for ipot')
         !
      END IF

  
   RETURN
   END SUBROUTINE RCM_Compute_ionospheric_potential





      SUBROUTINE Get_potential_on_bndy
      USE Rcm_variables
      IMPLICIT NONE
!
!
      INTEGER (iprec) :: j
      REAL    (rprec) :: r_eq, p_eq, Gntrp_2d, Gntrp_2d_ang
!
      IF (ipcp_type == 11) THEN
          DO j = 1, jsize
             vbnd (j) = -vdrop * SIN(aloct(1,j)-vdrop_phase*HTR ) / &
                         2.0_rprec * 1.0E+3_rprec
          END DO
!
      ELSE IF (ipcp_type == 13) THEN
!
!       Maynard and Chen [JGR, 1975]:
!
        DO j = 1, jsize
           r_eq = Gntrp_2d (rmin, isize, jsize, bndloc(j), REAL(j,rprec))
           p_eq = Gntrp_2d_ang (pmin, isize, jsize, bndloc(j), REAL(j,rprec))
           vbnd (j) = (92.4_rprec / r_eq - &
              A_coeff_MC(Kp)*r_eq**2*SIN(p_eq)) * 1000.0_rprec
        END DO
!
      ELSE
         call CON_stop('ERROR in IM/RCM2/src/rcm_comput.f90:'// &
             'VBOUND: IPCP_TYPE NOT IMPLEMENTED')
      END IF
      vbnd (0) = vbnd (jsize)
      vbnd (-1) = vbnd (jsize-1)
      vbnd (jsize+1) = vbnd (1)
      vbnd (jsize+2) = vbnd (2)
      RETURN
!
      CONTAINS
      FUNCTION A_coeff_MC (Kp)
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: Kp
      REAL (rprec) :: A_coeff_MC
      A_coeff_MC = 0.045_rprec / (1.0_rprec-0.159_rprec*Kp+0.0093_rprec*Kp**2)**3
      RETURN
      END FUNCTION A_coeff_MC
!
      END SUBROUTINE Get_potential_on_bndy
