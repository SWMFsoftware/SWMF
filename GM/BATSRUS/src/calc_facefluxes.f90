!^CFG COPYRIGHT UM
subroutine calc_facefluxes(DoResChangeOnly)

  use ModProcMH
  use ModMain, ONLY : &
       boris_correction,&               !^CFG IF BORISCORR
       BlkTest,PROCtest,globalBLK
  use ModAdvance, ONLY : FluxType

  implicit none

  logical, intent (in) :: DoResChangeOnly

  logical :: oktest, oktest_me, oktest_row
  integer::iVar
  !--------------------------------------------------------------------------

  if(iProc==PROCtest .and. globalBLK==BlkTest)then
     call set_oktest('calc_facefluxes',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if

  if(oktest_me) call print_values

  select case (FluxType)
  case ('Roe')                                  !^CFG IF ROEFLUX
     call calc_flux_Roe(DoResChangeOnly)        !^CFG IF ROEFLUX
  case ('Rusanov')                              !^CFG IF RUSANOVFLUX
     call calc_flux_Rusanov(DoResChangeOnly)    !^CFG IF RUSANOVFLUX
  case('Linde')                                 !^CFG IF LINDEFLUX
     call calc_flux_Linde(DoResChangeOnly)      !^CFG IF LINDEFLUX
  case ('Sokolov')                              !^CFG IF AWFLUX BEGIN
     if(boris_correction)then                      !^CFG IF BORISCORR BEGIN
        call calc_flux_AWboris(DoResChangeOnly)
     else                                          !^CFG END BORISCORR
        call calc_flux_AW(DoResChangeOnly)
     end if                                        !^CFG IF BORISCORR
                                                !^CFG END AWFLUX
  case default
     call stop_mpi("Invalid flux function in calc_facefluxes")
  end select

contains

  subroutine print_values

    use ModVarIndexes
    use ModMain, ONLY : iTest,jTest,kTest,DimTest,x_,y_,z_
    use ModAdvance

    if(DoResChangeOnly)then
       write(*,*)'calc_facefluxes for DoResChangeOnly'
       RETURN
    end if

    if(DimTest==x_ .or. DimTest==0)then
       write(*,*)&
            'Calc_facefluxes, left and right states at i-1/2 and i+1/2:'

       do iVar=1,nVar
          write(*,'(a,4(1pe13.5))')NameVar_V(iVar),'=',&
               LeftState_VX(iVar,iTest,jTest,kTest),&
               RightState_VX(iVar,iTest,  jTest,kTest),&
               LeftState_VX(iVar,iTest+1,jTest,kTest),&
               RightState_VX(iVar,iTest+1,jTest,kTest)
       end do
       write(*,'(a,1pe13.5,a13,1pe13.5)')'B0x:',&
            B0xFace_x_BLK(iTest,jTest,kTest,BlkTest),' ',&
            B0xFace_x_BLK(iTest+1,jTest,kTest,BlkTest)
       write(*,'(a,1pe13.5,a13,1pe13.5)')'B0y:',&
            B0yFace_x_BLK(iTest,jTest,kTest,BlkTest),' ',&
            B0yFace_x_BLK(iTest+1,jTest,kTest,BlkTest)
       write(*,'(a,1pe13.5,a13,1pe13.5)')'B0z:',&
            B0zFace_x_BLK(iTest,jTest,kTest,BlkTest),' ',&
            B0zFace_x_BLK(iTest+1,jTest,kTest,BlkTest)
    end if

    if(DimTest==y_ .or. DimTest==0)then
       write(*,*)&
            'Calc_facefluxes, left and right states at j-1/2 and j+1/2:'

       do iVar=1,nVar
          write(*,'(a,4(1pe13.5))')NameVar_V(iVar),'=',&
               LeftState_VY(iVar,iTest,jTest,kTest),&
               RightState_VY(iVar,iTest,  jTest,kTest),&
               LeftState_VY(iVar,iTest,jTest+1,kTest),&
               RightState_VY(iVar,iTest,jTest+1,kTest)
       end do
       write(*,'(a,1pe13.5,a13,1pe13.5)')'B0x:',&
            B0xFace_y_BLK(iTest,jTest,kTest,BlkTest),' ',&
            B0xFace_y_BLK(iTest,jTest+1,kTest,BlkTest)
       write(*,'(a,1pe13.5,a13,1pe13.5)')'B0y:',&
            B0yFace_y_BLK(iTest,jTest,kTest,BlkTest),' ',&
            B0yFace_y_BLK(iTest,jTest+1,kTest,BlkTest)
       write(*,'(a,1pe13.5,a13,1pe13.5)')'B0z:',&
            B0zFace_y_BLK(iTest,jTest,kTest,BlkTest),' ',&
            B0zFace_y_BLK(iTest,jTest+1,kTest,BlkTest)
    end if

    if(DimTest==z_ .or. DimTest==0)then
       do iVar=1,nVar
          write(*,'(a,4(1pe13.5))')NameVar_V(iVar),'=',&
               LeftState_VZ(iVar,iTest,jTest,kTest),&
               RightState_VZ(iVar,iTest,  jTest,kTest),&
               LeftState_VZ(iVar,iTest,jTest,kTest+1),&
               RightState_VZ(iVar,iTest,jTest,kTest+1)
       end do
       write(*,'(a,1pe13.5,a13,1pe13.5)')'B0x:',&
            B0xFace_z_BLK(iTest,jTest,kTest,BlkTest),' ',&
            B0xFace_z_BLK(iTest,jTest,kTest+1,BlkTest)
       write(*,'(a,1pe13.5,a13,1pe13.5)')'B0y:',&
            B0yFace_z_BLK(iTest,jTest,kTest,BlkTest),' ',&
            B0yFace_z_BLK(iTest,jTest,kTest+1,BlkTest)
       write(*,'(a,1pe13.5,a13,1pe13.5)')'B0z:',&
            B0zFace_z_BLK(iTest,jTest,kTest,BlkTest),' ',&
            B0zFace_z_BLK(iTest,jTest,kTest+1,BlkTest)
    end if

  end subroutine print_values

end subroutine calc_facefluxes
