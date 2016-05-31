module ModWrite

  ! This module contains methods for writing output files

  use ModSize, ONLY: &
       nDim, nVar, nLat, nLon, nNode, &
       RMin, iParticleMin, iParticleMax, nParticle,&
       Particle_, OriginLat_, OriginLon_, R_, Lat_, Lon_

  use ModGrid, ONLY: &
       get_node_indexes, &
       nBlock, State_VIB, iGrid_IA, iNode_B, &
       Proc_, Begin_, End_

  use ModPlotFile, ONLY: save_plot_file

  use ModIoUnit, ONLY: UnitTmp_

  implicit none

  SAVE

  private ! except
 
  public:: set_write_param, write_output, NamePlotDir

  !\
  !----------------------------------------------------------------------------
  ! Write modes:
  ! Data mode, i.e. what data to include in the output file
  integer:: OutputData
  integer, parameter:: &
       All_  = 0,  &
       Mesh_ = 1
  ! Write mode, i.e. the format of the output file
  integer:: OutputFormat
  integer, parameter:: &
       Tec_ = 0, &
       Idl_ = 1
  !----------------------------------------------------------------------------
  ! the output directory
  character (len=100) :: NamePlotDir="SP/IO2/"
  !/
contains
  
  subroutine set_write_param(StringPlot)
    ! set how output files are written:
    ! the data to be printed, the format of files
    character(len=*), intent(in):: StringPlot

    character(len=*), parameter :: NameSub='SP:set_write_param'
    !--------------------------------------------------------------------------
    ! the output data must be set in the input file
    if(    index(StringPlot,'msh') > 0)then
       OutputData = Mesh_
    elseif(index(StringPlot,'all') > 0)then
       OutputData = All_
    else
       OutputData = All_
!       call CON_stop(NameSub//': output data mode was not set in PARAM.in')
    end if

    ! the type of output file must be set
    if(    index(StringPlot,'tec') > 0)then
       OutputFormat = Tec_
    elseif(index(StringPlot,'idl') > 0)then
       OutputFormat = Idl_
    else
       call CON_stop(NameSub//': output format was not set in PARAM.in')
    end if
  end subroutine set_write_param

  !============================================================================

  subroutine write_output(Time, iIter)
    ! write the output data
    use ModNumConst, ONLY: cHalfPi
    use ModCoordTransform, ONLY: sph_to_xyz
    ! current time and iteration
    real,    intent(in):: Time
    integer, intent(in):: iIter

    ! file name
    character(len=100) :: NameFile, NameVar
    ! file counter
    integer, save:: iPlotFile = 1
    !    integer, parameter:: iVarMesh_I(3) = (/1,2,3/)

    ! loop variables
    integer:: iBlock, iParticle
    ! status for file open 
    integer:: iError

    ! cartesian and HGI coordinates
    real:: Xyz_D(nDim), Coord_D(nDim)

    character(len=*), parameter:: NameVarMesh = '"R", "Lat", "Lon"'
    character(len=*), parameter:: NameSub = 'SP:write_output'
    !--------------------------------------------------------------------------
    if(OutputFormat == Idl_)then
       call write_idl
       RETURN
    end if

    ! name of output file
    write(NameFile,'(a,i6.6)') &
         trim(NamePlotDir)//'lines_n',iPlotFile

    select case(OutputFormat)
    case(Tec_)
       NameFile = trim(NameFile) // '.dat'
    case(Idl_)
       NameFile = trim(NameFile) // '.out'
    end select

    ! open a new file
    open(UnitTmp_, file=NameFile, status="replace", iostat=iError)    

    ! check if it has been opened successfully
    if(iError /= 0) &
         call CON_stop(NameSub//': could not open '//trim(NameFile))

    ! file contents depends on the format
    select case(OutputFormat)
    case (Tec_)
       ! header for the file
       write(UnitTmp_,'(a)') 'Temporary header'
       ! list the data names
       select case(OutputData)
       case(Mesh_)
          write(UnitTmp_,'(a)') 'VARIABLES = "R", "Lat", "Lon"'
       case(All_)
          call CON_stop(NameSub//': not implemented')
       end select
       ! the format is scattered points
       write(UnitTmp_,'(a)')'ZONE F=POINT'
       ! print the data itself
       do iBlock = 1, nBlock
          do iParticle = iParticleMin, iParticleMax
             if(State_VIB(R_,iParticle, iBlock) < 0) CYCLE
             ! first, convert HGI coordinates to the cartesian
             Coord_D = State_VIB((/R_,Lat_,Lon_/), iParticle, iBlock)
             call sph_to_xyz(&
                  Coord_D(R_),cHalfPi-Coord_D(Lat_),Coord_D(Lon_), &
                  Xyz_D)
             ! write data to the output file
             write(UnitTmp_,'(100es18.10)') Xyz_D
             !             write(UnitTmp_,'(100es18.10)')State_VIB(iVarMesh_I,iParticle, iBlock)
          end do
       end do
    end select
    
    ! output is finished, close the file
    close(UnitTmp_)

    ! update the file counter
    iPlotFile = iPlotFile + 1

  contains

    subroutine write_idl
      ! write output file in the format to be read by IDL;
      ! separate fiel is created for each field line, name format is
      ! MH_data_<iLat>_<iLon>_n<iIter>.out
      !------------------------------------------------------------------------
      ! name of the output file
      character(len=100):: NameFile
      ! loop variables
      integer:: iBlock, iParticle
      ! indexes of corresponding node, latitude and longitude
      integer:: iNode, iLat, iLon
      ! index of first/last particle on the field line
      real:: PFirst
      real:: PLast
      ! coordinates in spherical coordinates
      real:: Coord_D(nDim)
      ! storage for output data
      real:: DataOut_VI(nDim,nParticle)
      !------------------------------------------------------------------------
      do iBlock = 1, nBlock
         iNode = iNode_B(iBlock)
         call get_node_indexes(iNode, iLon, iLat)

         ! set the file name
         write(NameFile,'(a,i3.3,a,i3.3,a,i6.6,a)') &
              trim(NamePlotDir)//'MH_data_',iLon,'_',iLat,'_n',iIter,'.out'

         ! get min and max particle indexes on this field line
         PFirst = iGrid_IA(Begin_, iNode)
         PLast  = iGrid_IA(End_,   iNode)

         ! convert coordinates to cartesian before output
         do iParticle = int(PFirst), int(PLast)
            Coord_D = State_VIB((/R_,Lat_,Lon_/), iParticle, iBlock)
            call sph_to_xyz(&
                 Coord_D(R_),cHalfPi-Coord_D(Lat_),Coord_D(Lon_), &
                 DataOut_VI(1:nDim,iParticle))
         end do

         ! print data to file
         call save_plot_file(&
              NameFile     = NameFile, &
              nDimIn       = 1, &
              TimeIn       = Time, &
              nStepIn      = iIter, &
              CoordMinIn_D = (/PFirst/), &
              CoordMaxIn_D = (/PLast/), &
              NameVarIn    = 'ParticleIndex X Y Z', &
              VarIn_VI     = DataOut_VI(1:nDim, int(PFirst):int(PLast))&
              )
      end do
      
    end subroutine write_idl

  end subroutine write_output

end module ModWrite
