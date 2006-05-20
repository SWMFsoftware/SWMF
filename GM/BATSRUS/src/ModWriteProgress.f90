!^CFG COPYRIGHT UM
!==============================================================================

subroutine write_progress(inopt)
  use ModProcMH
  use ModMain
  use ModIO, ONLY: iUnitOut, write_prefix
  implicit none

  integer, intent(in) :: inopt
  !----------------------------------------------------------------------------

  if (iProc /= 0 .or. lVerbose<=0) return

  select case(inopt)
  case (0)
     call write_prefix; write(iUnitOut,&
          '(1X,''BATSRUS: Block Adaptive Tree Solar-Wind '', &
          & ''Roe Upwind Scheme'')')
     call write_prefix; write(iUnitOut,&
          '(1X,''         for 3D Heliospheric Flows,'')')
     call write_prefix; write(iUnitOut,&
          '(1X,''University of Michigan, 1995-2004.'')')
     call write_prefix; write(iUnitOut,*)
     if(IsStandAlone)then
        write(*,'(a,f4.2,a,i4,a)') &
             ' BATSRUS version ',CodeVersion,&
             ' is running as '//NameThisComp//' on ',nProc,' PE(s)'
        write(*,*)
     end if
     call option_list
  case (1)
     call write_prefix; write(iUnitOut,*)
     call write_prefix
     if (time_accurate) then
        write(iUnitOut,*) 'Restarted run from N = ',n_step,' steps ', &
             'and T = ',Time_Simulation, &
             ' (',Time_Simulation/60.00, &
             ' min, ',Time_Simulation/3600.00,' hrs).'
     else
        write(iUnitOut,*) 'Restarted run from N = ',n_step,' steps.'
     end if
  end select
end subroutine write_progress

!==============================================================================

subroutine write_runtime_values()

  use ModIO, ONLY: iUnitOut, write_prefix
  use ModProcMH
  use ModMain
  use ModFaceValue, ONLY: TypeLimiter, BetaLimiter
  use ModAdvance,  ONLY : FluxType
  use ModGeometry, ONLY : x1,x2,y1,y2,z1,z2,minDXvalue,maxDXvalue,dx_BLK
  use ModParallel, ONLY : proc_dims
  use ModPhysics
  use ModMpi
  use ModUser, ONLY: user_write_progress
  implicit none

  integer :: iError
  real :: minDX,maxDX    
  !------------------------------------------------------------------------

  ! Find new min and max dx
  minDX=minval(dx_BLK, MASK=(.not.unusedBLK))
  maxDX=maxval(dx_BLK, MASK=(.not.unusedBLK))
  call MPI_allreduce(minDX,minDXvalue,1,MPI_REAL,MPI_MIN,iComm,iError)
  call MPI_allreduce(maxDX,maxDXvalue,1,MPI_REAL,MPI_MAX,iComm,iError)

  if (iProc /= 0 .or. lVerbose<=0) return

  call write_prefix; write(iUnitOut,*)
  call write_prefix; write(iUnitOut,*)'   Begin Numerical Simulation'
  call write_prefix; write(iUnitOut,*)'   --------------------------'
  call write_prefix; write(iUnitOut,*)
  call write_prefix; write(iUnitOut,*)
  call write_prefix; write(iUnitOut,*)'   Problem Type'
  call write_prefix; write(iUnitOut,*)'   ------------'
  call write_prefix; write(iUnitOut,*)
  call write_prefix; write(iUnitOut,'(10X,A14,I8)') 'problem_type: ',problem_type
  call write_prefix; write(iUnitOut,'(10X,A)') trim(StringProblemType_I(problem_type))
  call write_prefix; write(iUnitOut,*)
  call write_prefix; write(iUnitOut,*) '   Physical Model Input Solution Parameters'
  call write_prefix; write(iUnitOut,*) '   ----------------------------------------'
  call write_prefix; write(iUnitOut,*)
  if(problem_type==problem_heliosphere .or. problem_type==problem_cme)then
     call write_prefix; write(iUnitOut,'(10X,2(A13,E13.5))') &
          'Rhosun:      ',Rhosun,   ', Presun:    ',Presun
     call write_prefix; write(iUnitOut,'(10X,2(A13,E13.5))') &
          'SSPsun:      ',SSPsun,   ', Velsun:    ',Velsun
     call write_prefix; write(iUnitOut,'(10X,2(A13,E13.5))') &
          'Tsunrot:     ',Tsunrot
  endif
  if(problem_type==problem_heliosphere)then
     call write_prefix; write(iUnitOut,'(10X,2(A13,E13.5))') &
          'Qsun:        ',Qsun     ,', Theat:     ',Theat
     call write_prefix; write(iUnitOut,'(10X,2(A13,E13.5))') &
          'SIGMAheat:   ',SIGMAheat,', Rheat:     ',Rheat
  endif
  if(body1)then
     call write_prefix; write(iUnitOut,'(10X,2(A13,E13.5))') &
          'rBody:       ', rBody,      ', rPlanet:   ',unitSI_x
     call write_prefix; write(iUnitOut,'(10X,2(A13,E13.5))') &
          'Body_rho_dim:',Body_rho_dim,', Body_T_dim:',Body_T_dim
     call write_prefix; write(iUnitOut,'(10X,2(A13,E13.5))') &
          'Bdp:         ',Bdp      ,', Tilt:      ',ThetaTilt
     if(UseRotatingBc)then
        call write_prefix; write(iUnitOut,'(10X,a)') 'Corotation is used'
     end if
     if(UseGravity)then
        call write_prefix; write(iUnitOut,'(10X,a)') 'Gravity is used'
     end if
  else
     call write_prefix; write(iUnitOut,'(10X,''body1: .false.'')')
  end if
  call write_prefix; write(iUnitOut,*)
  !^CFG IF SECONDBODY BEGIN
  if (UseBody2) then
     call write_prefix
     write(iUnitOut,'(10X,''body2: .true.'')')
     call write_prefix; write(iUnitOut,'(10X,2(A13,3E13.5))') &
          'rBody2:      ',RBody2   ,', BdpDimBody2: ',&
          sqrt(sum(BdpDimBody2_D**2))
     call write_prefix; write(iUnitOut,'(10X,2(A13,E13.5))') &
          'xBody2:      ',xBody2   ,', Bdp2Body_x:    ',BdpBody2_D(1)
     call write_prefix; write(iUnitOut,'(10X,2(A13,E13.5))') &
          'yBody2:      ',yBody2   ,', Bdp2Body_y:    ',BdpBody2_D(2)
     call write_prefix; write(iUnitOut,'(10X,2(A13,E13.5))') &
          'zBody2:      ',zBody2   ,', Bdp2Body_z:    ',BdpBody2_D(3)
     call write_prefix
     write(iUnitOut,'(10X,2(A13,E13.5))') &
          'RhoDimBody2: ',RhoDimBody2,', TDimBody2:   ',TDimBody2
  else
     call write_prefix
     write(iUnitOut,'(10X,''body2: .false.'')')
  end if
  !^CFG END SECONDBODY

  call write_prefix; write(iUnitOut,*)
  call write_prefix; write(iUnitOut,'(10X,2(A13,E13.5))')&
       'cLIGHTfactor:',boris_cLIGHT_factor,', cLIGHT:    ',cLIGHT
  call write_prefix; write(iUnitOut,*)
  select case(problem_type)
  case(problem_shocktube, problem_uniform, problem_heliosphere, problem_cme)
     call write_prefix; write(iUnitOut,*)
  case default
     call write_prefix; write(iUnitOut,*)
     call write_prefix
     write(iUnitOut,'(10X,A19,F15.6,A11,F15.6)') 'SW_rho_dim [n/cc]: ',SW_rho_dim,',  SW_rho: ',SW_rho
     call write_prefix
     write(iUnitOut,'(10X,A19,F15.6,A11,F15.6)') 'SW_Ux_dim  [km/s]: ',SW_Ux_dim ,',  SW_Ux:  ',SW_Ux 
     call write_prefix
     write(iUnitOut,'(10X,A19,F15.6,A11,F15.6)') 'SW_Uy_dim  [km/s]: ',SW_Uy_dim ,',  SW_Uy:  ',SW_Uy 
     call write_prefix
     write(iUnitOut,'(10X,A19,F15.6,A11,F15.6)') 'SW_Uz_dim  [km/s]: ',SW_Uz_dim ,',  SW_Uz:  ',SW_Uz 
     call write_prefix
     write(iUnitOut,'(10X,A19,F15.6,A11,F15.6)') 'SW_p_dim   [ nPa]: ',SW_p_dim  ,',  SW_p:   ',SW_p  
     call write_prefix
     write(iUnitOut,'(10X,A19,F15.6,A11,F15.6)') 'SW_Bx_dim  [  nT]: ',SW_Bx_dim ,',  SW_Bx:  ',SW_Bx 
     call write_prefix
     write(iUnitOut,'(10X,A19,F15.6,A11,F15.6)') 'SW_By_dim  [  nT]: ',SW_By_dim ,',  SW_By:  ',SW_By 
     call write_prefix
     write(iUnitOut,'(10X,A19,F15.6,A11,F15.6)') 'SW_Bz_dim  [  nT]: ',SW_Bz_dim ,',  SW_Bz:  ',SW_Bz 
     call write_prefix
     write(iUnitOut,'(10X,A19,F15.6)')           'SW_a_dim   [km/s]: ',SW_a_dim
     call write_prefix
     write(iUnitOut,'(10X,A19,F15.6)')           'SW_T_dim   [   K]: ',SW_T_dim
  end select
  call write_prefix; write(iUnitOut,*)
  call write_prefix; write(iUnitOut,*)'   MHD Numerical Solution Parameters'
  call write_prefix; write(iUnitOut,*)'   ---------------------------------'
  call write_prefix; write(iUnitOut,*)
  select case (nORDER)
  case (1)
     call write_prefix
     write(iUnitOut,'(10X,''1st-Order Scheme'')')
  case (2)
     call write_prefix
     write(iUnitOut,'(10X,''2nd-Order Scheme'')')
     call write_prefix
     write(iUnitOut,'(10x,a,a)')'with limiter ',TypeLimiter
     if(TypeLimiter /= 'minmod') then
        call write_prefix
        write(iUnitOut,'(10x,a,e13.5)')'beta=',BetaLimiter
     end if
  end select
  call write_prefix
  if (time_accurate) then
     write(iUnitOut,'(10X,''Time accurate calculation'')')
  else
     write(iUnitOut,'(10X,''Steady state calculation'')')
  end if

  call write_prefix
  write(iUnitOut,'(10X,a,a)') FluxType,' Flux Function'

  call write_prefix
  if (UseImplicit) then                            !^CFG IF IMPLICIT BEGIN
     write(iUnitOut,'(10X,''Implicit Time Stepping'')')
  else                                             !^CFG END IMPLICIT
     write(iUnitOut,'(10X,''Explicit Time Stepping'')')
  end if                                           !^CFG IF IMPLICIT
  if(boris_correction)then                         !^CFG IF BORISCORR BEGIN 
     call write_prefix     
     write(iUnitOut,'(10X,''With Boris Correction, factor ='',f10.4)') &
          boris_cLIGHT_factor
  end if                                           !^CFG END BORISCORR
  if(UseBorisSimple)then                           !^CFG IF SIMPLEBORIS BEGIN
     call write_prefix
     write(iUnitOut,'(10X,''With Simple Boris Correction, factor ='',f10.4)') &
          boris_cLIGHT_factor                         
  end if                                           !^CFG END SIMPLEBORIS

  call write_prefix; write(iUnitOut,*)
  call write_prefix; write(iUnitOut,*)'   Other Runtime Parameters'
  call write_prefix; write(iUnitOut,*)'   ------------------------'
  call write_prefix; write(iUnitOut,*)
  call write_prefix; write(iUnitOut,*)'Available processors: nProc = ',nProc
  call write_prefix; write(iUnitOut,*)
  call write_prefix; write(iUnitOut,*)'After initial grid setup:'
  call write_prefix; write(iUnitOut,*)'  nBlockMax = ',nBlockMax,' nBLK = ',nBLK
  call write_prefix; write(iUnitOut,*)'  Total number of blocks used = ',nBlockALL
  call write_prefix; write(iUnitOut,*)'  Total number of cells = ',nBlockALL*nI*nJ*nK
  call write_prefix; write(iUnitOut,*)'  Smallest cell dx: ',minDXvalue,'  Largest cell dx: ',maxDXvalue
  call write_prefix; write(iUnitOut,*)
  call write_prefix
  write(iUnitOut,'(1x,a,3i8)')    'proc_dims:   ',proc_dims(1:3)
  call write_prefix
  write(iUnitOut,'(1x,a,3i8)')    'nCells:      ',nCells(1:3)
  call write_prefix
  write(iUnitOut,'(1x,a,2f16.8)') 'x:           ',x1,x2
  call write_prefix
  write(iUnitOut,'(1x,a,2f16.8)') 'y:           ',y1,y2
  call write_prefix
  write(iUnitOut,'(1x,a,2f16.8)') 'z:           ',z1,z2
  call write_prefix
  write(iUnitOut,'(1x,a,1i8)')    'multistage:  ',nSTAGE
  call write_prefix
  write(iUnitOut,'(1x,a,1f16.8)') 'cfl:         ',cfl
  call write_prefix
  write(iUnitOut,'(1x,a,4f16.8)') 'gamma:       ',g,gm1,gm2,gp1
  call write_prefix; write(iUnitOut,*)
  if(UseUserEchoInput) call user_write_progress
end subroutine write_runtime_values

!==============================================================================

subroutine write_timeaccurate

  use ModMain, ONLY : Time_Simulation
  implicit none

  write(*, '(a,e13.5,a,f12.6,a,f12.6,a)') &
       '   Simulated Time T = ',Time_Simulation, &
       ' (',Time_Simulation/60.00, &
       ' min, ',Time_Simulation/3600.00,' hrs)'

end subroutine write_timeaccurate

!===========================================================================

subroutine option_list

  use ModIO, ONLY: iUnitOut, write_prefix
  use ModUser, ONLY: NameUserModule, VersionUserModule
  implicit none

  logical            :: on
  character (len=40) :: name
  real               :: number
  !---------------------------------------------------------------------------
  call write_prefix; write(iUnitOut,'(a)') &
       '#=================================================================#'
  call OPTION_CONSTRAIN_B(on,name);  call write_option     !^CFG IF CONSTRAINB
  call OPTION_PROJECTION(on,name);   call write_option     !^CFG IF PROJECTION
  call OPTION_RAYTRACING(on,name);   call write_option     !^CFG IF RAYTRACE
  call OPTION_IMPLICIT(on,name);     call write_option     !^CFG IF IMPLICIT
  call write_prefix; write(iUnitOut,'(a)') &
       '#=================================================================#'

  call write_prefix; write(iUnitOut,'(a,f5.2)') &
       '# USER MODULE: '//NameUserModule,VersionUserModule
  call write_prefix; write(iUnitOut,'(a)') &
       '#=================================================================#'


  call write_prefix; write(iUnitOut,*)

contains

  subroutine write_option
    character(len=9) :: String
    if(on)then
       String=' is ON  #'
    else
       String=' is OFF #'
    end if
    call write_prefix; write(iUnitOut,'(a,a,a9,a)')'# OPTION ',name,' ',String

  end subroutine write_option

end subroutine option_list
