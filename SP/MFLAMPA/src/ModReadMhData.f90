module SP_ModReadMhData

  ! This module contains methods for reading input MH data

  use SP_ModSize, ONLY: &
       nDim, nLat, nLon, nNode, &
       iParticleMin, iParticleMax, nParticle,&
       nMomentumBin, &
       Particle_, OriginLat_, OriginLon_

  use SP_ModGrid, ONLY: &
       get_node_indexes, &
       iComm, &
       nVar, nBlock, State_VIB, iGridLocal_IB, iNode_B, &
       Distribution_IIB, LogEnergyScale_I, LogMomentumScale_I, &
       DMomentumOverDEnergy_I, &
       Proc_, Begin_, End_, X_, Y_, Z_, Bx_, By_, Bz_, &
       B_, Ux_, Uy_, Uz_, U_, Rho_, T_, S_, EFlux_, &
       NameVar_V

  use ModPlotFile, ONLY: read_plot_file

  implicit none

  SAVE

  private ! except
 
  public:: &
       set_read_mh_data_param, read_mh_data, &
       DoReadMhData


  !\
  !----------------------------------------------------------------------------
  ! Format of output files
  integer, parameter:: &
       Tec_ = 0, &
       Idl_ = 1
  !----------------------------------------------------------------------------
  ! the input directory
  character (len=100):: NameInputDir=""
  ! the input file name base
  character (len=100):: NameFileBase="MH_data"
  character (len=4)  :: NameFormat
  character (len=20) :: TypeFile

  real,allocatable:: Buffer_II(:,:)

  real:: TimeRead, TimeReadStart, TimeReadMax, DtRead
  integer:: iIterRead, iIterReadStart, DnRead

  logical:: DoReadMhData = .false.
  !/
contains
  
  subroutine set_read_mh_data_param
    use ModReadParam, ONLY: read_var
    ! set parameters of output files: file format, kind of output etc.
    character(len=300):: StringPlot
    ! loop variables
    integer:: iFile, iNode
    character(len=*), parameter :: NameSub='SP:set_read_mh_data_param'
    !--------------------------------------------------------------------------
    !
    call read_var('DoReadMhData', DoReadMhData)
    if(.not. DoReadMhData)&
         RETURN
    ! the input directory
    call read_var('NameInputDir', NameInputDir)
    ! ADD "/" IF NOT PRESENT
        
    !
    call read_var('TypeFile', TypeFile)

    ! time and iteration of the first file to be read
    call read_var('TimeReadStart',  TimeReadStart)
    call read_var('iIterReadStart', iIterReadStart)

    TimeRead = TimeReadStart
    iIterRead= iIterReadStart

    ! max time
    call read_var('TimeReadMax', TimeReadMax)

    ! time step
    call read_var('DtRead', DtRead)
    call read_var('DnRead', DnRead)

    ! the format of output file must be set
    select case(trim(TypeFile))
    case('tec')
       NameFormat='.dat'
    case('idl')
       NameFormat='.out'
    case default
       call CON_stop(NameSub//': input format was not set in PARAM.in')
    end select

    allocate(Buffer_II(14,nParticle))
  end subroutine set_read_mh_data_param

  !============================================================================

  subroutine read_mh_data(TimeOut, IsLast)
    real,    intent(out):: TimeOut
    logical, intent(out):: IsLast
    ! write the output data
!    real,    intent(in):: Time ! current time
!    integer, intent(in):: iIter! current iteration

    
      ! write output with 1D MH data in the format to be read by IDL/TECPLOT;
      ! separate file is created for each field line, name format is
      ! MH_data_<iLon>_<iLat>_n<iIter>.{out/dat}
      !------------------------------------------------------------------------
      ! name of the output file
      character(len=100):: NameFile
      ! loop variables
      integer:: iBlock, iParticle, iVarPlot
      ! indexes of corresponding node, latitude and longitude
      integer:: iNode, iLat, iLon, n
!!      ! index of first/last particle on the field line
!!      integer:: iFirst, iLast
!!      ! for better readability
!!      integer:: nVarPlot
character(len=300):: NameVar
      !------------------------------------------------------------------------
      do iBlock = 1, nBlock
         iNode = iNode_B(iBlock)
         call get_node_indexes(iNode, iLon, iLat)

         ! set the file name
         write(NameFile,'(a,i3.3,a,i3.3,a,i8.8,a,i6.6,a)') &
              trim(NameInputDir)//trim(NameFileBase)//'_',iLon,'_',iLat,&
              '_t',floor(TimeRead),'_n',iIterRead, NameFormat

         ! read the header first
         call read_plot_file(&
              NameFile = NameFile,&
              TypeFileIn = TypeFile,&
              NameVarOut = NameVar,&
              n1out = n)

         !\
         ! DETERMINE ORDER
         !/

         ! read the data itself
         call read_plot_file(&
              NameFile = NameFile,&
              TypeFileIn = TypeFile,&
              VarOut_VI = Buffer_II)

         State_VIB(&
              (/X_,Y_,Z_,S_,Rho_,T_,Ux_,Uy_,Uz_,U_,Bx_,By_,Bz_,B_/),&
              :,iBlock) = Buffer_II(:,1:n)

         iGridLocal_IB(Begin_,iBlock) = 1
         iGridLocal_IB(End_,iBlock) = n
      end do

      ! advance read time and iteration
      TimeRead = TimeRead + DtRead
      iIterRead = iIterRead + DnRead

      TimeOut = TimeRead

      IsLast = TimeRead > TimeReadMax

    end subroutine read_mh_data

end module SP_ModReadMhData
