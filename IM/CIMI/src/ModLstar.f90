!***********************************************************************
!                              calc_Lstar1  
!
! Routine modified from calc_Lstar. This routine calculates L * at given 
! (xlati,xmlt) in a given magnetic field configuration for the CIMI model 
! of Mei-Ching H. Fok at (Code 673) NASA/GSFC.
! Calculation of magnetic flux and L* is made at ionophere
!
! 
! In put: xlati,xmlt,bo,iba (field parameters)
! Output: L*, max L*
!
! Creaded by Suk-Bin Kang (Code 673) at NASA/GSFC on November 12, 2015
!
! Subroutines
! - calc_lstar1: calculates L* for K = 0 (or 90deg) particles.
! - calc_lstar2: calculates L* with K dependence.
! - locateB: is modified from locate1 and calculates latitudinal index
!            of bo or bm excluding magnetic island.
! - dBdMLAT: calculates magnetic flux at the ionosphere
!
! Parameters
! - BE: magnetic field strength at the Earth surface at magnetic equator.
! - Bflux: magnetic flux enclosed by drift path at each grid.
! - dBdlat: analytic integration of B * dMLAT
! - Lstar: out L*
! - Lstar_max: L* of last closed drift shell.
! - dmlt: MLT grid size at ionosphere.
! - mlt1,mlt2,dp,dmlt1,j1,j2: variables to be used in dmlt.
! - xlat_i, dBdlat_i: interpolated MLAT and dBdlat, respectively.
! - b1,b2,mlat2: variables to be used in xlat_i and dBdlat_i.
! - logBo_i: bo from i=1 to iba(j).
! - jdawn,jdusk,jnoon: MLT indices for dawn,dusk, and noon, respectively.
! - B_island: contains information if bo(i,j) is magnetic island
!           if 0, then normal, if 1, then magnetic island.
! - dBdMLAT: function that calculates integral B * dMLAT
!
! Modification history
!  * November 16, 2015 - correct interpolation of xlat_i between ii and ii+1
!                        to ii and ii1, skipping magnetic island.
!  * November 19, 2015 - To avoid magnetic island
!                        locateB: ju=j(ii)+1 -> j(ii)-1
!                                 jl=j(1)    -> j(1)+1
!  * November 20. 215  - direct linear interpolation -> log linear interpol.
!  * February 11, 2016 - allocable logBo_i(:) -> real logBo_i(ir)
!                      - remove allocate and deallocate logBo_i
!  * March 21, 2016    - Add calc_lstar2.f90
!**************************************************************************
  Module ModLstar

contains

  subroutine calc_Lstar1(Lstar,Lstar_max,rc)

  use ModNumConst, only:pi=>cPi
  use ModCimiPlanet, only:re_m,xme=>dipmom
  use ModFieldTrace,only: bo,ro,iba,rb
  use ModCimiGrid, only:xlati=>xlatr,ir=>np,ip=>nt,xmlt
  implicit none

  real BE,Bflux(ir,ip),logBo,dBdlat(ir),Lstar(ir,ip),Lstar_max,rc,&
       dmlt(ip),mlt1,mlt2,dp,dmlt1,xlat_i,dBdlat_i,b1,b2,mlat2,logBo_i(ir)
  integer i,j,j0,ib0,ii,ii1,jdawn,jdusk,jnoon,j1,j2,B_island(ir,ip)
  !real,external :: dBdMLAT

  jdawn=nint(float(ip)*0.75)+1     ! MLT index at dawn
  jdusk=nint(float(ip)*0.25)+1     ! MLT index at dusk
  jnoon=1      ! MLT index at noon
  dp=24./ip                        ! MLT grid size
  ib0=11          ! ro(ib0,*)=2.1, inner boundary of L* calculation

  BE=xme/re_m**3                    ! mangetic field at L* = 1 in


! Calulates dmlt at MLT grid 
  dmlt(1:ip)=24./ip

! Find magnetic island
  B_island(1:ir,1:ip)=0    ! 0=normal, 1=island
  do j=1,ip
     ii=1000
     do i=iba(j)-1,1,-1
        if (ii.lt.i) go to 20 
        if (bo(i,j).lt.bo(i+1,j)) then
           B_island(i,j)=1          
           ii=i
           do while ( bo(ii-1,j).lt.bo(i+1,j) .and. ii.gt.1)
             ii=ii-1
             B_island(ii,j)=1
           enddo
           write(*,'(" bo(i,j) has a island at, j = ",i3,"  i = ",i3," -",i3 )') j,ii,i
        endif
20      continue
     enddo
  enddo
   
! Caculates L* at each grid
  Bflux=0.  ! Initialize magnetic fluxes
  Lstar(1:ib0-1,1:ip)=ro(1:ib0-1,1:ip)  ! set L*=ro at i<ib0
  do j0=1,ip                     ! MLT grid
     do i=ib0,iba(j0)
        Bflux(i,j0)=0.             ! initialization of magnetic flux  
        if (B_island(i,j0).eq.1) then 
           lstar(i,j0)=-1.*rb      
           go to 40                ! skipping magnetic island
        endif
        logBo=log10(bo(i,j0))
        do j=1,ip
           logBo_i(1:iba(j))=log10(bo(1:iba(j),j))
           if (j.ne.j0.and.logBo.lt.logBo_i(1).and.logBo.ge.logBo_i(iba(j))) &
              call locateB(logBo_i(1:iba(j)),iba(j),logBo,ii) ! Find the MLAT index for the same bo
           if (logBo.ge.logBo_i(1)) ii=1
           if (logBo.lt.logBo_i(iba(j))) ii=iba(j)
           if (j.eq.j0) ii=i
           if (ii.lt.iba(j)) then
              ii1=ii+1
              do while (ii1.lt.iba(j) .and. B_island(ii1,j).eq.1 ) 
                 ii1=ii1+1
              enddo
              b2=logBo_i(ii1)-logBo
              b1=logBo-logBo_i(ii)
              xlat_i=(xlati(ii)*b2+xlati(ii1)*b1)/(b2+b1) ! interpolates xlati
           else
              xlat_i=xlati(iba(j))                   
           endif
           dBdlat_i=dBdMLAT(BE,rc,xlat_i)
           Bflux(i,j0)=Bflux(i,j0)+dBdlat_i*dmlt(j)*pi/12. ! analytic intergration of magnetic flux
        enddo          ! end of j
        lstar(i,j0)=2.*pi*BE/Bflux(i,j0)
40      continue
     enddo             ! end of i
     if (iba(j0).lt.ir) lstar(iba(j0)+1:ir,j0)=rb     ! arbitrary big L* 
  enddo                ! end of j0

  ! Determines L* max
  Lstar_max=Lstar(iba(jnoon),jnoon) ! L* at noon magnetopause
  do j=jdawn,jdusk                  ! find L* max from dawn to dusk
     if (Lstar_max.gt.Lstar(iba(j),j)) Lstar_max=Lstar(iba(j),j)
  enddo

write(*,'("L*max =",f8.2)') Lstar_max
write(*,'("L* =",10f8.2)') (Lstar(i,1),i=1,ir)
  end subroutine calc_Lstar1

!**************************************************************************
  subroutine calc_Lstar2(Lstar,Lstar_max,rc)
!
! Routine calcuates L* with K (second adiabatic invariant) dependence.
!**************************************************************************

  use ModNumConst, only:pi=>cPi
  use ModCimiPlanet, only:re_m,xme=>dipmom
  use ModFieldTrace,only: bm,ro,iba,rb
  use ModCimiGrid, only:xlati=>xlatr,ir=>np,ip=>nt,xmlt,ik=>nk
  implicit none

  real BE,Bflux(ir,ip),dBdlat(ir),Lstar(ir,ip,ik),Lstar_max(ik),rc,&
       dmlt(ip),mlt1,mlt2,dp,dmlt1,xlat_i,dBdlat_i,b1,b2,mlat2,logBm,&
       logbm_i(ir)
  integer i,j,m,j0,ib0,ii,ii1,jdawn,jdusk,jnoon,j1,j2,B_island(ir,ip)
 
  jdawn=nint(float(ip)*0.75)+1     ! MLT index at dawn
  jdusk=nint(float(ip)*0.25)+1     ! MLT index at dusk
  jnoon=1      ! MLT index at noon
  dp=24./ip                        ! MLT grid size
  ib0=11          ! ro(ib0,*)=2.1, inner boundary of L* calculation

  BE=xme/re_m**3                    ! mangetic field at L* = 1 in

! Calulates dmlt at MLT grid 
  dmlt(1:ip)=24./ip

  do m=1,ik
! Find magnetic island
     B_island(1:ir,1:ip)=0    ! 0=normal, 1=island
     do j=1,ip
        ii=1000
        do i=iba(j)-1,1,-1
           if (ii.lt.i) go to 20 
           if (bm(i,j,m).lt.bm(i+1,j,m)) then
              B_island(i,j)=1          
              ii=i
10            if(bm(ii-1,j,m).lt.bm(i+1,j,m)) then
                ii=ii-1
                B_island(ii,j)=1
              go to 10
              endif
              write(*,'(" bm(i,j,m) has a island at ,m,j = ",2i3,"  i = ",i3," -",i3 )') m,j,ii,i
           endif
20         continue
        enddo
     enddo
   
! Caculates L* at each grid
     Bflux=0.  ! Initialize magnetic fluxes
     Lstar(1:ib0-1,1:ip,m)=ro(1:ib0-1,1:ip)  ! set L*=ro at i<ib0
     do j0=1,ip                     ! MLT grid
        do i=ib0,iba(j0)
           Bflux(i,j0)=0.             ! initialization of magnetic flux  
           if (B_island(i,j0).eq.1) then 
              lstar(i,j0,m)=-1.*rb      
              go to 40                ! skipping magnetic island
           endif
           logBm=log10(bm(i,j0,m))
           do j=1,ip
              logbm_i(1:iba(j))=log10(bm(1:iba(j),j,m))            ! Drift path is estimated with bm 
              if (j.ne.j0.and.logBm.lt.logbm_i(1).and.logBm.ge.logbm_i(iba(j))) &
              call locateB(logbm_i(1:iba(j)),iba(j),logBm,ii) ! Find the MLAT index for the same bo
              if (logBm.ge.logbm_i(1)) ii=1
              if (logBm.lt.logbm_i(iba(j))) ii=iba(j)
              if (j.eq.j0) ii=i
              if (ii.lt.iba(j)) then
                 ii1=ii
30               ii1=ii1+1
                 if(B_island(ii1,j).eq.1) go to 30
                 b2=logbm_i(ii1)-logBm
                 b1=logBm-logbm_i(ii)
                 xlat_i=(xlati(ii)*b2+xlati(ii1)*b1)/(b2+b1) ! interpolates xlati
              else
                 xlat_i=xlati(iba(j))                   
              endif
              dBdlat_i=dBdMLAT(BE,rc,xlat_i)
              Bflux(i,j0)=Bflux(i,j0)+dBdlat_i*dmlt(j)*pi/12. ! analytic intergration of magnetic flux
           enddo ! end of j
           lstar(i,j0,m)=2.*pi*BE/Bflux(i,j0)
40         continue
        enddo    ! end of i
        if (iba(j0).lt.ir) lstar(iba(j0)+1:ir,j0,m)=rb     ! arbitrary big L* 
     enddo       ! end of j0

  ! Determines L* max
     Lstar_max(m)=Lstar(iba(jnoon),jnoon,m) ! L* at noon magnetopause
     do j=jdusk,jdawn                  ! find L* max from dawn to dusk
        if (Lstar_max(m).gt.Lstar(iba(j),j,m)) Lstar_max(m)=Lstar(iba(j),j,m)
     enddo
  enddo          ! end of k
write(*,'("L*max =",8f8.2)') (Lstar_max(m),m=1,ik)
write(*,'("L* =",10f8.2)') (Lstar(i,1,1),i=1,ir)

  end subroutine calc_Lstar2

!**************************************************************************
      subroutine locateB(xx,n,x,j)
!  
!  Routine modified from "locate1" in the CIMI model
!
!  Routine is return a value of j such that x is between xx(j) and xx(j+1),
!   excluding small humps and dips.
!  Routine is still valid if xx is a  monotonically increasing 
!   or decreasing function including small humps or dips by passing those.
!
!**************************************************************************

      implicit none
 
      integer n,i,ii,iii,j,jl,ju,jm,j1,j2
      integer js(n)
      real xx(n),x

!  Make sure xx is increasing or decreasing monotonically
      js(1:n)=1000     ! arbitrary big value 
      ii=0        ! dummy index to count magnetic island
      i=n
10    i=i-1
      if (xx(i).le.xx(i+1)) then
         ii=ii+1
         js(ii)=i
20       if(xx(i-1).lt.xx(i+1)) then
           ii=ii+1
           i=i-1
           js(ii)=i
           go to 20
         endif
      endif
      if (i.gt.1) go to 10
       

   if (js(1).gt.999) then
      jl=0
      ju=n+1
40    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      go to 40
      endif
      j=jl
   endif
   if (js(1).le.n) then
      jl=0
      ju=js(ii)
50    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      go to 50
      endif
      j1=jl
      jl=js(1)
      ju=n+1
60    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      go to 60
      endif
      j2=jl
      ! avoid a hump or dip
      j=j1
      if (x.le.xx(js(ii)-1).and.x.gt.xx(js(1)+1)) j=js(ii)-1
      if (x.le.xx(js(1)+1)) j=j2
   endif

      end subroutine locateB


!**************************************************************************
!                             dBdMLAT
!  Function calculates interation of dipole Br * dMLAT * r**2 at ionosphere
!  form mlat1 to MLAT = 90 degree
!**************************************************************************
     function dBdMLAT(BE,rc,mlat)

     implicit none
     real BE,rc,mlat,dBdMLAT

     dBdMLAT=BE/rc*(1.-sin(mlat)**2)

     end function dBdMLAT
  end Module ModLstar
