module ModWrite

  ! This module contains methods for writing output files

  use ModSize, ONLY: &
       nDim, nLat, nLon, nNode, &
       iParticleMin, iParticleMax, nParticle,&
       nMomentumBin, &
       Particle_, OriginLat_, OriginLon_

  use ModGrid, ONLY: &
       get_node_indexes, &
       nVar, nBlock, State_VIB, iGridLocal_IB, iNode_B, &
       Distribution_IIIB, MomentumScale_I, LogMomentumScale_I, &
       Proc_, Begin_, End_, R_, Lat_, Lon_, Bx_, By_, Bz_, &
       B_, Ux_, Uy_, Uz_, U_, Rho_, T_, S_, &
       NameVar_V

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
  character(len=4):: NameFormat
  character(len=20):: TypeFile
  !----------------------------------------------------------------------------
  ! List of MH variables that will be outputed
  logical:: DoPlot_V(nVar) = .false.
  ! Number of the variables to be outputed
  integer:: nVarPlot = 0
  ! Output Buffer
  real, allocatable:: DataOut_VI(:,:)
  ! Indices of the output variables in the state vector
  integer, allocatable:: iVarPlot_V(:)
  ! Coordinates are always fixed
  integer, parameter:: PlotX_ = 1, PlotY_ = 2, PlotZ_ = 3
  ! Names of variables
  character(len=300) :: NameVarPlot = ''
  !----------------------------------------------------------------------------
  ! the output directory
  character (len=100) :: NamePlotDir="SP/IO2/"
  !/
contains
  
  subroutine set_write_param(StringPlot)
    ! set how output files are written:
    ! the data to be printed, the format of files
    character(len=*), intent(in):: StringPlot

    integer:: iVar, iVarPlot
    character(len=*), parameter :: NameSub='SP:set_write_param'
    !--------------------------------------------------------------------------
    ! the type of output file must be set
    if(    index(StringPlot,'tec') > 0)then
       OutputFormat = Tec_
       NameFormat   ='.dat'
       TypeFile     ='tec'
    elseif(index(StringPlot,'idl') > 0)then
       OutputFormat = Idl_
       NameFormat   ='.out'
       TypeFile     ='ascii'
    else
       call CON_stop(NameSub//': output format was not set in PARAM.in')
    end if

    ! set variables to plot
    !----------------------
    ! coordinates are always printed
    DoPlot_V((/R_, Lon_, Lat_, S_/)) = .true.
    nVarPlot = nVarPlot + 4
    ! plasma density ----------
    if(index(StringPlot,'rho') > 0)then
       DoPlot_V(Rho_) = .true.
       nVarPlot = nVarPlot + 1
    end if
    ! temperature -------------
    if(index(StringPlot,'temp')> 0)then
       DoPlot_V(T_) = .true.
       nVarPlot = nVarPlot + 1
    end if
    ! velocity ----------------
    if(index(StringPlot,'ux')> 0 .or. index(StringPlot,'uy')> 0 .or. &
         index(StringPlot,'uz')> 0)then
       DoPlot_V((/Ux_, Uy_, Uz_/)) = .true.
       nVarPlot = nVarPlot + 3
    end if
    if(index(StringPlot,'|u|')> 0)then
       DoPlot_V(U_) = .true.
       nVarPlot = nVarPlot + 1
    end if
    ! magnetic field ----------
    if(index(StringPlot,'bx')> 0 .or. index(StringPlot,'by')> 0 .or. &
         index(StringPlot,'bz')> 0)then
       DoPlot_V((/Bx_, By_, Bz_/)) = .true.
       nVarPlot = nVarPlot + 3
    end if
    if(index(StringPlot,'|b|')> 0)then
       DoPlot_V(B_) = .true.
       nVarPlot = nVarPlot + 1
    end if
    !--------------------------
    ! prepare the output data container
    allocate(DataOut_VI(nVarPlot, iParticleMin:iParticleMax))
    ! indices in the state vector
    allocate(iVarPlot_V(nVarPlot))
    ! coordinates
    iVarPlot_V(PlotX_:PlotZ_) = (/1,2,3/)
    NameVarPlot  = 'X Y Z'
    ! the rest of variables
    iVarPlot = nDim + 1
    do iVar = nDim + 1, nVar
       if(.not.DoPlot_V(iVar)) CYCLE
       NameVarPlot = trim(NameVarPlot)//' '//trim(NameVar_V(iVar))
       iVarPlot_V(iVarPlot) = iVar
       iVarPlot = iVarPlot + 1
    end do
  end subroutine set_write_param

  !============================================================================

  subroutine write_output(Time, iIter)
    ! write the output data
    use ModNumConst, ONLY: cHalfPi
    use ModCoordTransform, ONLY: sph_to_xyz, rlonlat_to_xyz
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
    call write_mh_data
!    call write_distribution

  contains

    subroutine write_mh_data
      ! write output file with MH data in the format to be read by IDL/TECPLOT;
      ! separate file is created for each field line, name format is
      ! MH_data_<iLon>_<iLat>_n<iIter>.{out/dat}
      !------------------------------------------------------------------------
      ! name of the output file
      character(len=100):: NameFile
      ! loop variables
      integer:: iBlock, iParticle, iVarPlot
      ! indexes of corresponding node, latitude and longitude
      integer:: iNode, iLat, iLon
      ! index of first/last particle on the field line
      integer:: iFirst
      integer:: iLast
      ! coordinates in spherical coordinates
      real:: Coord_D(nDim)
      !------------------------------------------------------------------------
      do iBlock = 1, nBlock
         iNode = iNode_B(iBlock)
         call get_node_indexes(iNode, iLon, iLat)

         ! set the file name
         write(NameFile,'(a,i3.3,a,i3.3,a,i6.6,a)') &
              trim(NamePlotDir)//'MH_data_',iLon,'_',iLat,'_n',iIter,NameFormat

         ! get min and max particle indexes on this field line
         iFirst = iGridLocal_IB(Begin_, iBlock)
         iLast  = iGridLocal_IB(End_,   iBlock)

         do iParticle = iFirst, iLast
            ! convert coordinates to cartesian before output
            Coord_D = State_VIB((/R_,Lon_,Lat_/), iParticle, iBlock)
            call rlonlat_to_xyz(Coord_D, DataOut_VI(PlotX_:PlotZ_,iParticle))
            ! the rest of variables
            do iVarPlot = nDim + 1, nVarPlot

               DataOut_VI(iVarPlot, iParticle) = &
                    State_VIB(iVarPlot_V(iVarPlot), iParticle, iBlock)
            end do
         end do

         ! print data to file
         call save_plot_file(&
              NameFile     = NameFile, &
              TypeFileIn   = TypeFile, &
              nDimIn       = 1, &
              TimeIn       = Time, &
              nStepIn      = iIter, &
              CoordMinIn_D = (/real(iFirst)/), &
              CoordMaxIn_D = (/real(iLast)/), &
              NameVarIn    = 'ParticleIndex '//trim(NameVarPlot), &
              VarIn_VI     = DataOut_VI(1:nVarPlot, iFirst:iLast)&
              )
      end do     
    end subroutine write_mh_data

    !=========================================================================
    
    subroutine write_distribution
      ! write file with distribution in the format to be read by IDL/TECPLOT;
      ! separate fiel is created for each field line, name format is
      ! Distribution_<iLon>_<iLat>_n<iIter>.{out/dat}
      !------------------------------------------------------------------------
      ! name of the output file
      character(len=100):: NameFile
      ! loop variables
      integer:: iBlock, iParticle, iVarPlot
      ! indexes of corresponding node, latitude and longitude
      integer:: iNode, iLat, iLon
      ! index of first/last particle on the field line
      integer:: iFirst
      integer:: iLast
      ! coordinates in spherical coordinates
      real:: Coord_D(nDim)
      ! distribution function
      real:: F_II(iParticleMin:iParticleMax, nMomentumBin)
      ! auxliray array for scale with respect to particle index
      real, save:: IndexScale_I(iParticleMin:iParticleMax) = -1.0
      logical, save:: IsFirstCall = .true.
      !------------------------------------------------------------------------
      if(IsFirstCall)then
         IsFirstCall = .false.
         do iParticle = iParticleMin, iParticleMax
            IndexScale_I(iParticle) = real(iParticle)
         end do
      end if
      do iBlock = 1, nBlock
         iNode = iNode_B(iBlock)
         call get_node_indexes(iNode, iLon, iLat)

         ! set the file name
         write(NameFile,'(a,i3.3,a,i3.3,a,i6.6,a)') &
              trim(NamePlotDir)//&
              'Distribution_',iLon,'_',iLat,'_n',iIter,NameFormat

         ! get min and max particle indexes on this field line
         iFirst = iGridLocal_IB(Begin_, iBlock)
         iLast  = iGridLocal_IB(End_,   iBlock)
         
         do iParticle = iParticleMin, iParticleMax
            ! reset values outside the line's range
            if(iParticle < iFirst .or. iParticle > iLast)then
               F_II(iParticle,:) = 0.0
               CYCLE
            end if
            ! the actual distribution
            F_II(iParticle,:) = log(Distribution_IIIB(:,1,iParticle,iBlock))
         end do

         ! print data to file
         call save_plot_file(&
              NameFile     = NameFile, &
              TypeFileIn   = TypeFile, &
              nDimIn       = 2, &
              TimeIn       = Time, &
              nStepIn      = iIter, &
              Coord1In_I   = State_VIB(S_,iFirst:iLast,iBlock), &!IndexScale_I, &
              Coord2In_I   = LogMomentumScale_I, &
              NameVarIn    = 'Distance LogMomentum LogDistribution', &
              VarIn_II     = F_II(iFirst:iLast,:) &
              )
      end do     
    end subroutine write_distribution

  end subroutine write_output

end module ModWrite
