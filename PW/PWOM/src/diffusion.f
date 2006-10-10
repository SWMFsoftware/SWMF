CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This subroutine adds diffusion to the code. An array representing
C     any physical quantity can be inputed. A delta array is calculated 
C     that is then added back to the origional array. The delta array
C     is the effect of diffusion.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine diffusion (Var_C,VarBottom,VarTop,cMax_C)
      use ModCommonVariables
      real, intent(inout) :: Var_C(MaxGrid)
      
      real, intent(in)    :: VarBottom, VarTop, cMax_C(MaxGrid)
      
      real :: DiffVar_C(MaxGrid),DiffusionCoef(MaxGrid)
      
      integer iMin, iMax

      logical IsGlobal
      
c      IsGlobal=.true.
      IsGlobal=.false.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Set the range of the Diffusion region
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c      iMin=2
c      iMax=ndim-1
c      iMax=nDim-20
c      iMin=800

c      iMin=465
c      iMax=840

      iMin=1
      iMax=nDim     

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Set the diffusion coeficient. The coeficient must be between 
C     0 and 1/2 for stability. The coef. is multiplied by the courant
C     number. This helps to keep the amount of diffusion from depending
C     on the length of the timestep. This can be done with the global or
C     local courant number
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      !set coeficient in first order scheme
c      DiffusionCoef= 5.0e-2
c      if (IsGlobal) then
c         CflMax = maxval(Cfl)
c         DiffusionCoef = DiffusionCoef * MaxCfl
c      else
c         DiffusionCoef = DiffusionCoef * Cfl
c      endif

c      write(*,*) cfl

c      write(*,*) diffusioncoef


      
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C     Calculate the diffusion. dU(i)=diffcoef*[U(i+1)-2U(i)+U(i-1)]
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      ! Calculate diffusion
c      do i=iMin, iMax
c         DiffVar_C(i) = DiffusionCoef(i)*(Var_C(i+1)-2*Var_C(i)+Var_C(i-1))
c      end do
      
      call calc_rusanov(nDim, Var_C, VarBottom, VarTop, cMax_C, DiffVar_C)

      
      ! Apply diffusion
      do i=iMin, iMax
c     local time stepping
c         Var_C(i)=Var_C(i) + Dt/DrBnd*DiffVar_C(i)

c     Global time stepping

         Var_C(i)=Var_C(i) + 0.5*DiffVar_C(i)




      end do

      return
      end
      

      
