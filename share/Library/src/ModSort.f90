module ModSort

  implicit none

contains

  subroutine sort_quick(n,arr,indx)

    ! Quick sort algorithm: sorts indx array according to 'arr'
    ! so that arr(indx(i)) <= arr(indx(i+1)) for i=1..n

    ! Based on the F77 code indexx from Numerical Recipes

    integer, intent(in)  :: n
    integer, intent(out) :: indx(n)
    real, intent(in)     :: arr(n)

    integer, parameter   :: M=7, NSTACK=50

    integer :: i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
    real    :: a

    !-------------------------------------------------------------------
    do j=1,n
       indx(j)=j
    end do
    jstack=0
    l=1
    ir=n
1   if(ir-l < M)then
       do j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do i=j-1,1,-1
             if(arr(indx(i)).le.a)goto 2
             indx(i+1)=indx(i)
          end do
          i=0
2         indx(i+1)=indxt
       end do
       if(jstack.eq.0)return
       ir=istack(jstack)
       l=istack(jstack-1)
       jstack=jstack-2
    else
       k=(l+ir)/2
       itemp=indx(k)
       indx(k)=indx(l+1)
       indx(l+1)=itemp
       if(arr(indx(l+1)) > arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
       endif
       if(arr(indx(l)) > arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
       endif
       if(arr(indx(l+1)) > arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
       endif
       i=l+1
       j=ir
       indxt=indx(l)
       a=arr(indxt)
3      continue
       i=i+1
       if(arr(indx(i)) < a)goto 3
4      continue
       j=j-1
       if(arr(indx(j)) > a)goto 4
       if(j < i)goto 5
       itemp=indx(i)
       indx(i)=indx(j)
       indx(j)=itemp
       goto 3
5      indx(l)=indx(j)
       indx(j)=indxt
       jstack=jstack+2
       if(jstack > NSTACK) &
            stop 'ERROR in ModSort::sort_quick: NSTACK too small'
       if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
       else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
       endif
    endif
    goto 1

  end subroutine sort_quick
  
  !==========================================================================
  subroutine sort_test

    integer, parameter :: n=5
    real    :: a_I(n)
    integer :: i_I(n)

    a_I = (/0.0, 2.0, 1.0, 4.0, 0.0/)

    call sort_quick(n,a_I,i_I)

    write(*,*)'a      =',a_I
    write(*,*)'indx   =',i_I
    write(*,*)'a(indx)=',a_I(i_I)

  end subroutine sort_test

end module ModSort
