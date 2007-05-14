
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
     &        ALTD(K),Source_GV(K,RhoO_),Source_GV(K,RhoH_),Source_GV(K,RhoHe_),Source_GV(K,RhoE_),Source_GV(K,uO_),Source_GV(K,uH_),
     &        Source_GV(K,uHe_),Source_GV(K,uE_),Source_GV(K,pO_),Source_GV(K,pE_),Source_GV(K,pHe_),Source_GV(K,pE_),
     &        EFIELD(K), GRAVTY(K)
      enddo

      RETURN
      END
