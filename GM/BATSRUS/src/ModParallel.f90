!^CFG COPYRIGHT UM
module ModParallel
  use ModSize
  implicit none
  save

  !\
  ! Domain Decomposition
  !/
  integer, dimension(3) :: proc_dims       ! Initial layout of processors
  logical :: periodic3D(3) = .false.

  !\
  ! Code timing values.
  !/
  real, dimension(:), allocatable :: PE_float_list, PE_time_list, PE_mflop_list

  !\
  ! Index arrays to convert global block number to local block and proc index,
  ! and the block number used for the restart file name.
  !/
  integer, allocatable, dimension(:) :: iBlock_A, iProc_A, iBlockRestartALL_A

  !\                                               !^CFG IF IMPLICIT BEGIN
  ! logical for explicit vs. implicit blocks
  !/
  logical, allocatable :: implicitBlock_BP(:,:)    !^CFG END IMPLICIT

  !\
  ! Neighbor solution block refinement levels
  ! ( 0=neighbors at same level, 
  !  -1=neighbors at lower level,
  !  +1=neighbors at higher level,
  !  NOBLK=no neighbors).
  !/
  integer, parameter :: NOBLK=-100

  integer, dimension(nBLK) :: &
       neiLtop, neiLbot, neiLeast, neiLwest, neiLnorth, neiLsouth

  integer, dimension(east_:top_,nBLK):: neiLEV

  !\
  ! Neighbor processor and block numbers (a value of NOBLK
  ! means not used).  As only one level change is permitted
  ! between neighboring solution blocks, there are either 1 or 4 
  ! neighboring blocks in each of the six directions.
  !/
  integer, dimension(4,nBLK) :: &
       neiPtop, neiPbot, neiPeast, neiPwest, neiPnorth, neiPsouth, &
       neiBtop, neiBbot, neiBeast, neiBwest, neiBnorth, neiBsouth

  integer, dimension(4,east_:top_,nBLK) :: neiPE, neiBLK

  integer, dimension( -1:1, -1:1, -1:1, 4, nBLK) :: &
       BLKneighborPE, BLKneighborBLK, BLKneighborCHILD
  integer, dimension( -1:1, -1:1, -1:1, nBLK) :: BLKneighborLEV

  !\
  ! Variables for nonblocking sends used in fix_refine and fix_coarsen.
  !/
  integer, parameter :: NO_BLK=-1, NUM_NB_BUFFS=12

  integer :: nb_buffs(3,NUM_NB_BUFFS,0:nBLK)
  integer :: curr_nb_buff(0:nBLK)
  integer :: nb_buff_req(NUM_NB_BUFFS,0:nBLK)

  ! Use special options in exchange_messages if creating plot files
  logical :: UsePlotMessageOptions = .false.

end module ModParallel
