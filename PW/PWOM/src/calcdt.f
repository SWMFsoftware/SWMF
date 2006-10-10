CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Set the adaptive timestep, based on the pc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine calcdt(pc)
      use ModCommonVariables
      real pc
      real MaxPc
      parameter (MaxPc =.1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     If the percent change of the pressure is large cut the timestep in
C     half. If the percent change of the pressure is small increase 
C     the timestep by .1% never allow to exceed cfl.The sonic speed is 
C     used to calculate cfl
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      !refine dt
c      if (IsVariableDt) then         
c         if (pc > MaxPc) then
c            dt = min(dt/2.,cfl_dt)
c            DTX1=DT
c            DTR1=DTX1/DRBND
c            DTX2=DT*NTS
c            DTR2=DTX2/DRBND
c            
c            write(*,*) "percent change, dt : ",pc, dt, time
c         else
c            dt = min(dt * 1.01, cfl_dt)
c            DTX1=DT
c            DTR1=DTX1/DRBND
c            DTX2=DT*NTS
c            DTR2=DTX2/DRBND
c            endif
c         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         
         
         
         return
         end
