!^CFG COPYRIGHT UM
module ModInterface

  interface
     subroutine message_pass_dir(&
          idirmin,idirmax,width,sendghostcells,prolongorder,nVar,&
          sol_BLK,sol2_BLK,sol3_BLK,sol4_BLK,sol5_BLK,sol6_BLK,sol7_BLK,&
          sol8_BLK,sol9_BLK,Sol_VGB,restrictface)

       use ModSize
       implicit none

       integer, intent(in) :: idirmin,idirmax, width
       logical, intent(in) :: sendghostcells
       integer, intent(in) :: prolongorder, nvar

       real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), optional,intent(inout)::&
            sol_BLK,sol2_BLK,sol3_BLK,sol4_BLK,sol5_BLK,sol6_BLK,sol7_BLK,&
            sol8_BLK,sol9_BLK
       real, optional, intent(inout) :: &
            Sol_VGB(nVar,-1:nI+2,-1:nJ+2,-1:nK+2,nBLK)
       logical, optional, intent(in):: restrictface

     end subroutine message_pass_dir
  end interface

end module ModInterface
