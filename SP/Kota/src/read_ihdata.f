!      call read_ihdata_for_sp(
!     1     114     ,1,1)
!      stop
!      end
      subroutine read_ihdata_for_sp(
     1     NumFileIn,nLineHeader,nLineBad)
      implicit none
      include 'stdout.h'
      include 'coupler.h'
      include 'param.h'   
      integer nr,nmu,nw
      real dim1
      common /size/ nr,nmu,nw, dim1
      !--------------------------------------------------!
      integer  NumFileIn
      integer nLineHeader
      integer nLineBad
      !----------------------------------------!
      character*30 NameFile
      character*80 TextHeader
      integer iLine,nPoint,iFile,iError
      real Factor
      real TimeToRead,tSimulation
      common/contime/ TimeToRead,tSimulation
      integer iV
      real cSecondPerHour
      PARAMETER(cSecondPerHour=3600.00000000000000000)
      real cE1
      PARAMETER(cE1=10.000000000000000000000000000000)
      real cE6
      PARAMETER(cE6=cE1**6)
      real cE9
      PARAMETER(cE9=cE1**9)
      real Rsun
      PARAMETER(Rsun=0.696*cE9)
      real cAU
      PARAMETER(cAU= 149.59787000000000000*cE9)
      real cProtonMass
      PARAMETER( cProtonMass=1.6726/cE9**3)
      
      !------------------------------------------!
      write(NameFile,11)'./IO_IH/IHDATA/ihdata',
     1     NumFileIn,
     2     '.dat'
 11   FORMAT(a,i5.5,a)

      write(iStdout,*)prefix,'Reading IH data from ',NameFile
      write(*,*)NameFile
 !     stop
      !\
      ! First reading of file sets the value of nX::
      !/
      call get_io_unit_new(iFile)
      open(iFile,file=NameFile,status='old',iostat=iError)
      if (nLineHeader>0) then
         do iLine=1,nLineHeader
            read(iFile,'(a)') TextHeader
            read(TextHeader(20:40),*) TimeToRead
         end do
      end if
      nPoint=0
      iError=0
      do while(iError.eq.0)
         read(iFile,*,iostat=iError)(Smooth_VII(iV,iLine,Old_),iV=1,11)
         nPoint=nPoint+1
      end do
      close(iFile)
      iMax=min(nPoint-nLineBad,nRMax+1)
      nr=iMax-1
      !\
      ! Second reading of file with data from IH_::
      !/
      call get_io_unit_new(iFile)
      open(iFile,file=NameFile,status='old',iostat=iError)
      if (nLineHeader>0) then
         do iLine=1,nLineHeader
            read(iFile,'(a)') TextHeader
            write(iStdout,*)prefix,'at ',TextHeader
         end do
      end if
      do iLine=1,iMax
         do iV=1,11
            Smooth_VII(iV,iLine,Old_)=Smooth_VII(iV,iLine,New_)
         end do
         read(iFile,*,iostat=iError) 
     1         (Smooth_VII(iV,iLine,New_),iV=1,11)            
         Factor=Rsun/cAU   ! in AU
         do iV=1,3
            Smooth_VII(iV,iLine,New_)=
     1           Factor*Smooth_VII(iV,iLine,New_)
         end do
         rx(iLine)=Smooth_VII(1,iLine,New_)
         ry(iLine)=Smooth_VII(2,iLine,New_)
         rz(iLine)=Smooth_VII(3,iLine,New_)
         
         Smooth_VII(4,iLine,New_)=Smooth_VII(4,iLine,New_)/
     1        (cProtonMass*cE6) ! in [cm^-3]
         
         dd(iLine) = Smooth_VII(4,iLine,New_)     ! in [cm^-3]

         Factor=cSecondPerHour/cAU     ! in [AU/Hr]
         do iV=5,7
            Smooth_VII(iV,iLine,New_)=
     1           Factor*Smooth_VII(iV,iLine,New_)
         end do
         vx(iLine)=Smooth_VII(5,iLine,New_)
         vy(iLine)=Smooth_VII(6,iLine,New_)
         vz(iLine)=Smooth_VII(7,iLine,New_)

         Factor=cE9     ! in [nT]
         do iV=8,10
            Smooth_VII(iV,iLine,New_)=
     1        Factor*Smooth_VII(iV,iLine,New_)
         end do
         bx(iLine)=Smooth_VII(8,iLine,New_)
         by(iLine)=Smooth_VII(9,iLine,New_)
         bz(iLine)=Smooth_VII(10,iLine,New_)
         
         Smooth_VII(11,iLine,New_)=cE1*
     1        Smooth_VII(11,iLine,New_)
         pp(iLine)   = Smooth_VII(11,iLine,New_) ! in [erg/cm^3]
      end do
      close(iFile)
      end subroutine read_ihdata_for_sp
!=============================================================================


  
