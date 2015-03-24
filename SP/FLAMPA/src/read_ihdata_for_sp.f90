!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModReadMhData
  use SP_ModMain, ONLY : State_VI, X_DI, nX, iStdOut, prefix
  implicit none
  SAVE
  !\
  ! Directory with input data fiels from BATSRUS
  !/
  character(LEN=*),parameter:: IH_dir='./SP/MHDATA/'
contains
  !============================================================================!
  subroutine read_ihdata_for_sp(nLineHeader,nLineBad)
    use SP_ModMain, ONLY : SP_allocate, DataInputTime, iDataSet, DInner_I 
    use ModIOUnit, ONLY  : io_unit_new
    !-----------------------------------------------!
    !\
    ! Input parameters: how many lines should be ignored at the beginning and 
    ! at the end 
    !/
    integer,intent(in):: nLineHeader
    integer,intent(in):: nLineBad
    !\
    ! File name to be composed, for the given file number
    !/
    character(LEN=80):: NameFile
    !\
    !String variable, to read the header line(s)
    !/
    character(LEN=80):: NameHeader
    !\
    ! Misc
    !/
    integer:: iLine,nPoint,iFile,iError=0
    !\
    ! State vecor to read
    real :: State_V(11)
    !--------------------------------------------------------!
    write(NameFile,'(a,i4.4,a)')IH_dir//'mhdata_',iDataSet,'.dat'
    write(iStdout,*)prefix,'Reading IH data from ',trim(NameFile)
    !\
    ! First reading of file sets the value of nX::
    !/
    iFile=io_unit_new()
    open(iFile,file=NameFile,status='old',iostat=iError)
    if(iError>0)call CON_stop('File does not exist')
    if (nLineHeader>0) then
       do iLine=1,nLineHeader
          read(iFile,'(a)') NameHeader
          read(NameHeader(20:40),*) DataInputTime
       end do
    end if
    nPoint=0
    do
       read(iFile,*,iostat=iError) State_V(:)
       if(iError/=0)exit
       nPoint=nPoint+1
    end do
    close(iFile)
    nX=nPoint-nLineBad
    write(iStdOut,*)prefix,' nX=',nX
    if(.not.allocated(DInner_I))call SP_allocate
    nX=min(nX,ubound(DInner_I,1))
    
    !\
    ! Second reading of file with data from IH_::
    !/
    iFile=io_unit_new()
    open(iFile,file=NameFile,status='old',iostat=iError)
    if (nLineHeader>0) then
       do iLine=1,nLineHeader
          read(iFile,'(a)') NameHeader
          write(iStdout,*)prefix,'at ',NameHeader
       end do
    end if

    do iLine=1,nX
       read(iFile,*)X_DI(:,iLine), State_VI(:,iLine)
    end do
    close(iFile)
    call mh_transform_for_flampa
  end subroutine read_ihdata_for_sp
  !==================================
  subroutine mh_transform_for_flampa
    use ModConst,   ONLY: cProtonMass, rSun
    use SP_ModMain, ONLY: BSmoothOld_I, BSmooth_I, RhoSmoothOld_I, RhoSmooth_I,&
                          i_shock, T_I, iShock, iShockOld, U_I, NameEnergyUnit
    implicit none
    integer        ::  iPoint
    real, external ::  energy_in
    !\
    ! Named indexes
    !/
    integer,parameter:: Rho_=1, Ux_=2, Uz_=4, Bx_=5, Bz_=7, p_=8 
    integer,parameter:: x_=1, y_=2, z_=3
    
    RhoSmoothOld_I(1:nX)= RhoSmooth_I(1:nX)
    BSmoothOld_I(1:nX)= BSmooth_I(1:nX)
    iShockOld=iShock
    
    do iPoint=1,nX
       X_DI(:,iPoint)=X_DI(x_:z_,iPoint)*Rsun                     ! in [m]
       RhoSmooth_I(iPoint) =State_VI(Rho_,iPoint)/cProtonMass     ! in [m^-3]
       T_I(iPoint)   =State_VI(p_,iPoint)/RhoSmooth_I(iPoint)/&   
            energy_in(NameEnergyUnit)                ! in NameEnergyUnit [KeV]
       BSmooth_I(iPoint)   =sqrt(sum(State_VI(Bx_:Bz_,iPoint)**2))! in [Tesla]
       U_I(iPoint)   =sqrt(sum(State_VI(Ux_:Uz_,iPoint)**2))      ! in [m/s]
    end do
    !\
    ! Write min/max values of MHD variables being read from file::
    !/
    write(iStdout,*)prefix//'nX is set to: ',nX 
    write(iStdout,*)prefix,'min/max value of number density in [m^-3]    = ',&
         minval(RhoSmooth_I),maxval(RhoSmooth_I)
    write(iStdout,*)prefix,'min/max value of temperature in ['//             &
         trim(NameEnergyUnit)//']       = ',minval(T_I),maxval(T_I)
    write(iStdout,*)prefix,'min/max value of radial distance in [Rsun]   = ',&
         sqrt(minval(X_DI(x_,:)**2+X_DI(y_,:)**2+X_DI(z_,:)**2))/Rsun, &
         sqrt(maxval(X_DI(x_,:)**2+X_DI(y_,:)**2+X_DI(z_,:)**2))/Rsun
    write(iStdout,*)prefix,'min/max value of field strength in [Tesla]   = ',&
         minval(BSmooth_I),maxval(BSmooth_I)
    iShock=i_shock(RhoSmooth_I(1:nX))
  end subroutine mh_transform_for_flampa
  !======================================================
end module SP_ModReadMhData
