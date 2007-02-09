module ModPwGrid
  
  real, allocatable :: CoordPw_DI(:,:), StatePw_VI(:,:)
  real, allocatable :: CoordXY_DI(:,:)
  integer, parameter :: Theta_=1, Phi_=2, Rho_=1, RhoUb_=2
  integer :: nVarPw=-1, nLinePw=-1

  ! variables for triangulation
  integer, allocatable   :: triangle_node(:,:)
  real   , allocatable   :: table(:,:)
  integer :: node_num
  integer triangle_num
contains
  !==================================================

  subroutine set_PwGrid_buffer
    
    logical :: DoInitialize
    integer :: X_=1, Y_=2,nCoord=2
    !-----------------------------------------------

    if (DoInitialize)then
       allocate(Table(1:nCoord,1:nLinePw))
       allocate(CoordXY_DI(1:nCoord,1:nLinePw))
       allocate(table(1:nCoord,1:nLinePw))
       DoInitialize = .false.
    endif

    CoordXY_DI(X_,:) =  &
         rBoundary*sin(CoordPw_DI(Theta_,:)) * cos(CoordPw_DI(Phi_,:))
    CoordXY_DI(Y_,:) =  &
         rBoundary*sin(CoordPw_DI(Theta_,:)) * sin(CoordPw_DI(Phi_,:))

    
    Table(:,:) = CoordXY_DI(:,:)
    
    ! create the triangulation
    call table_delaunay
    
    
  end subroutine set_PwGrid_buffer
  
  !===========================================================
  
  subroutine get_from_PwGrid_buffer(CoordIN_D)
    real, intent(in)    :: CoordIN_D(2)
    integer, parameter  :: X_=1, Y_=2
    real                :: Node_DI(2,3)
    real                :: Triangle1_DI(2,3), &
                           Triangle2_DI(2,3), Triangle3_DI(2,3)
    logical             :: IsTriangleFound
    integer    :: node1,node2,node3
    real       :: Area1,Area2,Area3
    !--------------------------------------------------------

        
    !Find triangle containing point
    call find_triangle(&
         CoordIN_D, node1,node2,node3,IsTriangleFound)

    

    Node_DI(:,1)=CoordXY_DI(:,node1)
    Node_DI(:,2)=CoordXY_DI(:,node2)
    Node_DI(:,3)=CoordXY_DI(:,node3)
         
    
    

    
    ! interpolate values
    If (IsTriangleFound) then
       Triangle1_DI(:,1)   = CoordIN_D(:)
       Triangle1_DI(:,2:3) = Node_DI(:,1:2)

       Triangle2_DI(:,1)   = CoordIN_D(:)
       Triangle2_DI(:,2:3) = Node_DI(:,2:3)

       Triangle3_DI(:,1)   = CoordIN_D(:)
       Triangle3_DI(:,2) = Node_DI(:,2)
       Triangle3_DI(:,3) = Node_DI(:,3)
    
       Area1=Area(Triangle1_DI)
       Area2=Area(Triangle2_DI)
       Area3=Area(Triangle3_DI)

       RhoOut= (&
            Area1*StatePw_VI(Rho_,node3) &
            + Area2*StatePw_VI(Rho_,node1)&
            + Area3*StatePw_VI(Rho_,node2))&
            /(Area1+Area2+Area3)
       RhoUOut= (&
            Area1*StatePw_VI(RhoUb_,node3) &
            + Area2*StatePw_VI(RhoUb_,node1)&
            + Area3*StatePw_VI(RhoUb_,node2))&
            /(Area1+Area2+Area3)
    else
       !don't interpolate since not in triangle
    end If
  end subroutine get_from_PwGrid_buffer
  
  !==========================================================

  function Area(Node_DI)
    !use heron's formula to get triangle area
    real,intent(in)       :: Node_DI(2,3)
    integer, parameter  :: X_=1, Y_=2

    !-------------------------------------------------------

    Side1 = sqrt( &
         (Node_DI(X_,1)-Node_DI(X_,2))**2.0 &
         +(Node_DI(Y_,1)-Node_DI(Y_,2))**2.0)
    Side2 = sqrt( &
         (Node_DI(X_,2)-Node_DI(X_,3))**2.0 &
         +(Node_DI(Y_,2)-Node_DI(Y_,3))**2.0)
    Side3 = sqrt( &
         (Node_DI(X_,1)-Node_DI(X_,3))**2.0 &
         +(Node_DI(Y_,1)-Node_DI(Y_,3))**2.0)
    
    s = 0.5*(Side1+Side2+Side3)
    
    Area = sqrt(s*(s-Side1)*(s-Side2)*(s-Side3))
    return
  end function Area
end module ModPwGrid
