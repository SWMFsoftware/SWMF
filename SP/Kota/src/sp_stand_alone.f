      program sep   
      implicit none
      integer iLine         !The # of the magnetic filed line to process
      integer iSnapshot     !The # of coupling 
      real TimeToRead,tSimulation
      common/contime/TimeToRead,tSimulation
      real tFinal
      real tCoupleHr
      common /couple/tCoupleHr
      include 'stdout.h'
      logical UseSelfSimilarity,UseRefresh
      common/log/UseSelfSimilarity,UseRefresh

      call sp_init


      if(DoWriteAll)write(iStdout,*) prefix, 
     1     'Select line (1-2-3): '
CCC   if(DoWriteAll)write(iStdout,*) prefix, 'PRESS 1'
CCC   read(*,*)  lineno

      if(UseSelfSimilarity)then
         iLine=1
         if(DoWriteAll)write(iStdout,*) prefix, 
     1        'lineno=',iLine
CCC   if(DoWriteAll)write(iStdout,*) prefix, 'PRESS 9'
CCC   read(*,*)  iSnapshot
         iSnapshot=9            


         call sp_get_from_ih(iLine,iSnapshot)
      
! Physical time to start the update
         tSimulation=tCoupleHr*3.60d3 
! The start time is read from the snapshot
         if(DoWriteAll)write(iStdout,*)'tSimulation=',tSimulation

! Final value of physical time is not defined
         tFinal=1.0e8

         call sp_run(tSimulation,tFinal)
      else
         tSimulation=-1.0
         do iLine=114,197
            call read_ihdata_for_sp(iLine,1,5)
            call sp_sharpen_and_run(tSimulation,TimeToRead)
         end do
      end if
      write(iStdout,*)'tSimulation=',tSimulation 
      
      call closetime  !Finalize
      stop
      end 
!=============================================================!
      subroutine sp_read_inputs(nLine,NameList_I)
      implicit none
      integer::nLine
      integer::iLine,Misc
      character*80::NameList_I(nLine)
      call get_io_unit_new(Misc)
      open(Misc,file = 'violet.in')
      do iLine=1,nLine
         read(Misc,'(a)')Namelist_I(iLine)
      end do
      close(Misc)
      end
!===========================================================
      subroutine sp_get_from_ih(iLine,iSnapshot)
      implicit none
      include 'stdout.h'
      include 'coupler.h'
      real tCoupleHr
      common /couple/tCoupleHr
      integer iLine
      integer iFile
      integer iLoop
      integer iLoopSnapshot
      integer iSnapshot
      integer iMisc
C ======= START COMMUNICATION WITH IH COMPONENT HERE =====================

      call get_io_unit_new(iFile)

      if (iLine.eq.1) open(iFile,file = 'evolv1.cme' )
      if (iLine.eq.2) open(iFile,file = 'evolv2.cme' )
      if (iLine.eq.3) open(iFile,file = 'evolv3.cme' )

      read(iFile,*)
      read(iFile,*)
      read(iFile,*) iLine,iMisc  
      read(iFile,*) 
      write(*,*)DoWriteAll
      if(DoWriteAll)write(iStdout,*) prefix,
     1           'energy in megavolts'
      if(DoWriteAll)write(iStdout,*) prefix,
     1         'lineno,kmax :  ',iLine,iMisc
      if(DoWriteAll)write(iStdout,*) prefix   
      if(DoWriteAll)write(iStdout,*) prefix,
     1           'selected line: ',iLine
      if(DoWriteAll)write(iStdout,*) prefix,
     1       'Which snapshot to take (1-20): '
   
      do  iLoopSnapshot=1,iSnapshot
         read(iFile,*)  iMisc,iMax,tCoupleHr
      

         do iLoop=1,iMax 
            read(iFile,*)  iMisc,rx(iLoop),ry(iLoop),rz(iLoop)
            read(iFile,*)  iMisc,vx(iLoop),vy(iLoop),vz(iLoop)
            read(iFile,*)  iMisc,bx(iLoop),by(iLoop),bz(iLoop)
            read(iFile,*)  iMisc,dd(iLoop),pp(iLoop)
         end do

         read(iFile,*) 
      end do
      close(iFile)
      if(DoWriteAll)write(iStdout,*) prefix,
     1      'consider snapshot : ',iSnapshot
      end subroutine sp_get_from_ih

C ======= END COMMUNICATION WITH IH COMPONENT HERE =====================
