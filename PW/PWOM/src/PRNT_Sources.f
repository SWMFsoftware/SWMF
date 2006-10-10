
C ============================================================================

      SUBROUTINE PRNT_Sources
      use ModCommonVariables 

      write (iUnitSourceGraphics,"(a79)") 'Polarwind output_var11'
      write (iUnitSourceGraphics,"(i7,1pe13.5,3i3)") nint(time/dt),time,1,1,14
      write (iUnitSourceGraphics,"(3i4)") nDim
      write (iUnitSourceGraphics,"(100(1pe13.5))") Gamma
      write (iUnitSourceGraphics,"(a100)") 
     &     'alt Mass_O Mass_H Mass_He Mass_e Mom_O Mom_H Mom_He Mom_e E_O E_H E_He E_e Efield Gravity gamma'
      
      DO K=1,NDIM


         WRITE (iUnitSourceGraphics,"(100(1pe18.10))") 
     &        ALTD(K),ADMSO(K),ADMSH(K),ADMSHE(K),ADMSE(K),FCLSNO(K),FCLSNH(K),
     &        FCLSHE(K),FCLSNE(K),ECLSNO(K),ECLSNE(K),ECLSHE(K),ECLSNE(K),
     &        EFIELD(K), GRAVTY(K)
      enddo

      RETURN
      END
