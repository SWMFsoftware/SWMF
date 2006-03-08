!^CFG COPYRIGHT UM
!========================================================================
module ModUser
  ! This is the default user module which contains empty methods defined
  ! in ModUserEmpty.f90
  !
  ! Please see the documentation, and the files ModUserEmpty.f90 and 
  ! srcUser/ModUserExamples.f90 for information about what the different user
  ! subroutines do and how to implement them for your specific problem.

  use ModUserEmpty, ONLY:               &
!!!       user_read_inputs,                &
       user_init_session,               &
!!!       user_set_ics,                    &
       user_initial_perturbation,       &
       user_set_boundary_cells,         &
       user_face_bcs,                   &
       user_set_outerbcs,               &
       user_specify_initial_refinement, &
       user_amr_criteria,               &
       user_write_progress,             &
       user_get_log_var,                &
       user_calc_sources,               &
       user_heat_source,                &
       user_get_b0,                     &
       user_update_states

  include 'user_module.h' !list of public methods

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'Hall MHD test, Yingjuan Ma'

  real :: &
       AmplBx=0., AmplBy=0., AmplBz=0., &
       AmplUx=0., AmplUy=0., AmplUz=0., &
       WidthB=1,  lamdax=4.0,kwave
  integer :: iDim             
       
contains

  subroutine user_read_inputs
    use ModMain
    use ModProcMH,    ONLY: iProc
    use ModReadParam
    use ModNumConst, ONLY : cTwo,cPi
    character (len=100) :: NameCommand
    !-------------------------------------------------------------------------

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#PERTUBATION")
          call read_var('iDim',iDim)
          call read_var('WidthB',WidthB)
          call read_var('lamdax',lamdax)
          kwave = cTwo*cPi/lamdax 
          write(*,*)'lamda,k=',lamdax, kwave
          select case(iDim)
             case(x_)
                call read_var('AmplBy',AmplBy)
                call read_var('AmplBz',AmplBz)
                call read_var('AmplUy',AmplUy)                            
                call read_var('AmplUz',AmplUz)
             case(y_)
                call read_var('AmplBz',AmplBz)
                call read_var('AmplBx',AmplBx)
                call read_var('AmplUz',AmplUz)                            
                call read_var('AmplUx',AmplUx)
             case(z_)
                call read_var('AmplBx',AmplBx)
                call read_var('AmplBy',AmplBy)
                call read_var('AmplUx',AmplUx)                            
                call read_var('AmplUy',AmplUy)
             case default
                write(*,*)'wrong direction'
             end select
       case('#USERINPUTEND')
          if(iProc==0) write(*,*)'USERINPUTEND'
          EXIT
       case default
          if(iProc==0) call stop_mpi( &
               'read_inputs: unrecognized command: '//NameCommand)
       end select
    end do
  end subroutine user_read_inputs

  subroutine user_set_ics
    use ModMain,     ONLY: nI, nJ, nK, globalBLK, nBlock, ProcTest,x_,y_,z_
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, dx_BLK, dy_BLK
    use ModAdvance,  ONLY: State_VGB, Rho_,  RhoUx_, RhoUy_, RhoUz_, P_, Bx_, By_, Bz_
    use ModProcMH,   ONLY: iProc

    real, parameter :: pPerturb = 1.1
    integer :: i, j, k, iBlock
    !--------------------------------------------------------------------------
    iBlock = globalBLK
    
    if(iProc==PROCtest .and. iBlock == 1)then
       write(*,*)'Initializing HALL MHD TEST problem'
       write(*,*)'Parameters:'
       write(*,*)'WidthB=',WidthB,'lamda=',lamdax, 'kwave =',kwave
       write(*,*)'AmplBx =',AmplBx, 'AmplBy =',AmplBy,'AmplBz =',AmplBz
       write(*,*)'AmplUx =',AmplUx, 'AmplUy =',AmplUy,'AmplUz =',AmplUz
    end if

    select case(iDim)
    case(x_)    
       where(abs(x_BLK(:,:,:,iBlock))<WidthB)          
          State_VGB(By_,:,:,:,iBlock)=                    &
               AmplBy*cos(kwave*x_BLK(:,:,:,iBlock))          
          State_VGB(Bz_,:,:,:,iBlock)=                    &
               AmplBz*sin(kwave*x_BLK(:,:,:,iBlock))
          State_VGB(RhoUy_,:,:,:,iBlock)=                 &
               AmplUy*State_VGB(Rho_,:,:,:,iBlock)        & 
               *cos(kwave*x_BLK(:,:,:,iBlock))
          State_VGB(RhoUz_,:,:,:,iBlock)=                 &
               AmplUz*State_VGB(Rho_,:,:,:,iBlock)        & 
               *sin(kwave*x_BLK(:,:,:,iBlock))
       end where
    case(y_)    
       where(abs(y_BLK(:,:,:,iBlock))<WidthB)          
          State_VGB(Bz_,:,:,:,iBlock)=                    &
               AmplBz*cos(kwave*y_BLK(:,:,:,iBlock))          
          State_VGB(Bx_,:,:,:,iBlock)=                    &
               AmplBx*sin(kwave*y_BLK(:,:,:,iBlock))
          State_VGB(RhoUz_,:,:,:,iBlock)=                 &
               AmplUz*State_VGB(Rho_,:,:,:,iBlock)        & 
               *cos(kwave*y_BLK(:,:,:,iBlock))
          State_VGB(RhoUx_,:,:,:,iBlock)=                 &
               AmplUx*State_VGB(Rho_,:,:,:,iBlock)        & 
               *sin(kwave*y_BLK(:,:,:,iBlock))
       end where
    case(z_)    
       where(abs(z_BLK(:,:,:,iBlock))<WidthB)          
          State_VGB(Bx_,:,:,:,iBlock)=                    &
               AmplBx*cos(kwave*z_BLK(:,:,:,iBlock))          
          State_VGB(By_,:,:,:,iBlock)=                    &
               AmplBy*sin(kwave*z_BLK(:,:,:,iBlock))
          State_VGB(RhoUx_,:,:,:,iBlock)=                 &
               AmplUx*State_VGB(Rho_,:,:,:,iBlock)        & 
               *cos(kwave*z_BLK(:,:,:,iBlock))
          State_VGB(RhoUy_,:,:,:,iBlock)=                 &
               AmplUy*State_VGB(Rho_,:,:,:,iBlock)        & 
               *sin(kwave*z_BLK(:,:,:,iBlock))
       end where
    case default
       write(*,*)'wrong direction'                       
    end select

  end subroutine user_set_ics

end module ModUser
