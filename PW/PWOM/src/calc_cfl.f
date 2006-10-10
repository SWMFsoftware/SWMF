CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Calculates the Max speeds and cfl condition
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      subroutine calc_cfl
      use ModCommonVariables
      
      !calculate the max speed 
      cMax_O = SQRT(GAMMA * POXYG  /DOXYG) + abs(UOXYG)
      cMax_H = SQRT(GAMMA * PHYD   /DHYD)  + abs(UHYD)
      cMax_e = SQRT(GAMMA * PELECT /DELECT)+ abs(UELECT)
      
      !calculate the cfl
      cfl_H  = dt / drbnd * cMax_H
      cfl_O  = dt / drbnd * cMax_O
      cfl_e  = dt / drbnd * cMax_O
      
      cfl = max(cfl_H, cfl_O, cfl_e)
      
      MaxCfl = maxval(cfl)
      
      cfl_dt = dt/MaxCfl  

c      write(*,*) Maxcfl, cfl_dt

      return
      end
