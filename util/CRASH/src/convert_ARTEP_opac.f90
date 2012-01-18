module ModReadArtep
  use ModIoUnit, ONLY: io_unit_new
  implicit none

  integer:: nRho = 201, nTe = 201, nGroup = 30

  SAVE
  real, allocatable:: Value_VII(:,:,:) 
  real, allocatable:: EGroup_I(:)
  integer::line
  real:: TMin, TMax, RhoMin, RhoMax
  character(LEN=15)NameEmpty
contains
  !============================================================================
  subroutine read_artep
    integer:: iFile, iRho, iTe
    real:: Aux_V(7), Aux
    !-----------------------------------------
    iFile = io_unit_new()
    open(iFile, file='header.C5O2H8',status = 'old')
    do line = 1,17
       read(iFile,*)
    end do
    read(iFile,*)NameEmpty, TMin, TMax, nTe
    write(*,*)'NameEmpty, TMin, TMax, nTe=',NameEmpty, TMin, TMax, nTe
    read(iFile,*)NameEmpty, RhoMin, RhoMax, nRho
    write(*,*)'NameEmpty, RhoMin, RhoMax, nRho=', &
         NameEmpty, RhoMin, RhoMax, nRho
    read(iFile,*)
    !read(iFile,*)
    !read(iFile,*)
    read(iFile,*)NameEmpty,NameEmpty,NameEmpty,nGroup
    write(*,*)'nGroup=', nGroup
    allocate(Value_VII(2*nGroup,nRho,nTe), EGroup_I(1+nGroup))
   
    read(iFile,*)EGroup_I
    close(iFile)
    write(*,*)EGroup_I

    iFile = io_unit_new()
    open(iFile, file='opa.C5O2H8', form='formatted',status='old')
    line=0
    do iTe = 1, nTe
       do iRho = 1, nRho
          line = line+1;write(*,*)line
          read(iFile,*) Aux_V, Value_VII(1:nGroup,iRho,iTe), &
               Aux,            Value_VII(nGroup+1:2*nGroup,iRho,iTe)
       end do
    end do
    close(iFile)

  end subroutine read_artep

end Module ModReadArtep
!==============================================================================
program save_eos_table

  use CRASH_ModStatSum
  use CRASH_ModPartition, ONLY: UseCoulombCorrection
  use CRASH_ModPolyimide
  use CRASH_ModFermiGas
  use CRASH_ModMultiGroup
  use CRASH_ModExcitation
  use CRASH_ModAtomicDataMix
  use CRASH_ModExcitationData,ONLY : n_ground, cExcitationN_III
  use CRASH_ModIonMix
  use ModReadArtep
  use ModPlotFile
  use ModMpi
  use CRASH_ModAtomicMass,ONLY: cAtomicMass_I
  use ModConst
  implicit none
  !----------------------------------------------------------------------------
  call read_artep
  Value_VII = 0.10 * Value_VII  !cm^2/g = 1e-4 m^2/(1e-3 kg) = 0.1 m2/kg
  call save_plot_file( &
         'Ay_opac.dat',                                                      &
         TypeFileIn     = 'real8',                                           &
         StringHeaderIn = 'ARTEP Opacity for Ay',                            &
         NameVarIn      = 'logRho logTe Planck(30) Ross(30) Ev(31)',         &
         CoordMinIn_D   = (/log10(1.0e3*RhoMin), log10(TMin)/),              &
         CoordMaxIn_D   = (/log10(1.0e3*RhoMax), log10(TMax)/),              &
         ParamIn_I      = EGroup_I,                                          &
         VarIn_VII      = Value_VII)

end program save_eos_table
!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  implicit none
  character (len=*), intent(in) :: StringError
  write(*,*)StringError
  stop
end subroutine CON_stop
!============================================================================
subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe
end subroutine CON_set_do_test
