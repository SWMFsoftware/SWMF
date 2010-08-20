Module ModFieldTrace
  use ModCrcmGrid,ONLY: ir=>np, ip=>nt, iw=>nm , ik=>nk, neng, energy
  use ModCrcmPlanet,ONLY: nspec, amu_I
  implicit none
  real, allocatable :: &
       bo(:,:), ro(:,:), xmlto(:,:), sinA(:,:,:), &
       Have(:,:,:), pp(:,:,:,:,:), vel(:,:,:,:,:),&
       ekev(:,:,:,:), rmir(:,:,:), alscone(:,:,:,:,:),&
       tanA2(:,:,:), volume(:,:), bm(:,:,:), gamma(:,:,:,:),&
       xo(:,:), yo(:,:), tya(:,:,:), gridoc(:,:)

  real:: parmod(10), rb

  integer :: irm(ip),irm0(ip),iba(ip)

  integer :: iw2(ik)

  logical :: UseEllipse = .true.
contains
  subroutine init_mod_field_trace
    
    if(allocated(bo)) RETURN
    allocate( bo(ir,ip),ro(ir,ip),xmlto(ir,ip),sinA(ir,ip,ik),&
         Have(ir,ip,ik),pp(nspec,ir,ip,iw,ik),vel(nspec,ir,ip,iw,ik),&
         ekev(ir,ip,iw,ik),rmir(ir,ip,ik),alscone(nspec,ir,ip,iw,ik),&
         tanA2(ir,ip,0:ik+1),&
         volume(ir,ip),bm(ir,ip,ik),gamma(ir,ip,iw,ik),&
         xo(ir,ip),yo(ir,ip),tya(ir,ip,0:ik+1),gridoc(ir,ip) )

  end subroutine init_mod_field_trace

  !***********************************************************************
  !                             fieldpara
  ! Routine calculates kinetic energy, velocity, y, latitude and altitude
  ! at mirror point, etc, for given magnetic moment, K and position for a
  ! given magnetic field configuration.
  ! Output: iba,irm,iw2,vel,ekev,pp,sinA,Have,alscone             
  !***********************************************************************
  subroutine fieldpara(t,dt,c,q,rc,re,xlati,xmlt,phi,si,&
       xme)! uncomment when T04 tracing fixed ! ,xmp)
    use ModTsyInput, ONLY: tf,dsth,tdst,byw,bzw,timf,xnswa,vswa,tsw,&
         ndst,nimf,nsw,js,iyear,iday,imod
    use ModNumConst, ONLY: pi => cPi, cDegToRad
    use ModCrcmInitialize, ONLY: w=>xmm

    ! uncomment when T04 Tracing fixed
    ! common/geopack/aa(10),sps,cps,bb(3),ps,cc(11),kk(2),dd(8)
    ! external tsyndipoleSM,MHD_B
    
    integer, parameter :: np=2545,nd=3
    real :: rc, re,xme,dt,t,c
    real xlati(ir),phi(ip),si(ik),&
         si3(np),bm1(np),rm(np),rs(np),dss(np),&
         h3(np),bs(np),bba(np),&
         x1(ir),xmlt(ip),xli(0:ir),&
         ra(np),dssa(np),tya3(np)
    ! coeff and integrals in Taylor expansion
    integer,parameter ::nTaylor=10
    real a_I(0:nTaylor),b_I(0:nTaylor),sumBn(0:nTaylor),sumhBn(0:nTaylor), &
         BnI(0:nTaylor,np),hBnI(0:nTaylor,np)
    integer :: n,i,j,k,m,mir,npf,npf1,im,im1,im2,igood,ii,iTaylor
    integer :: iopt, n5, iout,n8,n7,m0,n70,ib
    real    :: rlim,xmltlim,dre,xlati1,phi1,xmlt1,ro1,volume1,bo1,dss2
    real    :: dssm,rm1,rme,rl,cost,bsn,bsndss,dssp,bsp,ss,ss1,ss2,bm_n
    real    :: bnibm_n,sim,bmmx,rmm, bs_n,tya33,h33,xmass,c2mo,c4mo2,ro2,pp1
    real    :: pijkm,pc,c2m,e,q,tcone1,tcone2,x
    integer :: iLatTest = -1, iLonTest=-1
!    integer :: iLatTest = 51, iLonTest=27

    integer :: imax
    real :: R_12,R_24,xmltr,xBoundary(ip),MajorAxis,MinorAxis,&
            MajorAxis2,MinorAxis2, sin2, Req2, xo1,xc, xCenter,rell2
    real, parameter :: LengthMax = 50.0
    !--------------------------------------------------------------------------

    a_I(0)=1.
    b_I(0)=1.
    do iTaylor=1,nTaylor
       a_I(iTaylor)=a_I(iTaylor-1)*(2.*iTaylor-3.)/(2.*iTaylor)
       b_I(iTaylor)=b_I(iTaylor-1)*(2.*iTaylor-1.)/(2.*iTaylor)
    enddo

    rb=10.               ! nightside outer boundary in RE
    iopt=1               ! dummy parameter in t96_01 and t04_s 
    rlim=2.*rb
    xmltlim=2.           ! limit of field line warping in hour
    dre=0.06             ! r interval below the surface of the Earth
    n5=16                ! no. of point below the surface of the Earth

    ! Save irm0
    irm0(1:ip)=irm(1:ip)

! uncomment when T04 Tracing fixed
!    if (imod <= 2) then
!       !  Determine parmod
!       call TsyParmod(t,tf,tsw,xnswa,vswa,nsw,tdst,dsth,ndst,timf,byw,bzw,&
!            nimf,xmp,imod,parmod)
!
!       !  Call recalc to calculate the dipole tilt
!       isec=mod(ifix(t),60)
!       min1=ifix(t)/60
!       minu=mod(min1,60)
!       ihour1=ifix(t)/3600
!       ihour=mod(ihour1,24)
!       iday1=iday+ifix(t)/86400
!       if (imod.le.2) then     
!          ps=0.                ! force ps = 0 when using Tsy models
!          cps=cos(ps)
!          sps=sin(ps)
!       else
!          call recalc(iyear,iday1,ihour,min,isec)
!       endif
!    endif
    !  Start field line tracing.  
    call timing_start('rbe_trace')
    LONGITUDE: do j=1,ip
       irm(j)=ir
       LATITUDE: do i=1,ir
          iout=0
          xlati1=xlati(i)*cDegToRad
          xli(i)=rc/cos(xlati1)/cos(xlati1)
          phi1=phi(j)                  ! +x corresponing to noon

          ! uncomment when T04 Tracing fixed
          !if (imod.le.2) call tsy_trace(i,rlim,re,rc,xlati1,phi1,t,ps,parmod,&
          !     imod,np,npf1,dssa,bba,volume1,ro1,xmlt1,bo1,ra)
          if (imod.eq.3) call Mhd_trace_IM(xlati1,phi(j),re,i,j,np, &
               npf1,dssa,bba,volume1,ro1,xmlt1,bo1,ra)


          if (i==iLatTest .and. j==iLonTest) then
             write(*,*) npf1,xlati1*180.0/3.14,xmlt1,ro1
             call IM_plot_fieldline(npf1,i,j,dssa,ra,bba) 
          endif
          volume(i,j)=volume1
          ro(i,j)=ro1
          if (i.gt.1) then
             if (ro(i,j).lt.ro(i-1,j)) ro(i,j)=ro(i-1,j)
          endif
          xmlto(i,j)=xmlt1
          bo(i,j)=bo1
          phi1=xmlt1*pi/12.+pi         ! phi1=0 corresponing to noon

          xo(i,j)=ro1*cos(phi1)
          yo(i,j)=ro1*sin(phi1)
          gridoc(i,j)=1.
          if (npf1.eq.0) gridoc(i,j)=0.

          !xmlto(i,j)=xmlt_1
          !if (j.gt.1.and.xmlto(i,j).lt.0.) xmlto(i,j)=xmlto(i,j)+24.
          !ro(i,j)=ra(ieq)
          !bo(i,j)=bba(ieq)
          !xo(i,j)=xa(ieq)
          !yo(i,j)=ya(ieq)
          !if (iout.ge.1) gridoc(i,j)=0.
          !if (iout.eq.0) gridoc(i,j)=1.

          if (npf1.eq.0) then             ! open field line
             irm(j)=i-1
             exit LATITUDE                   ! do next j                 
          endif

          if (xmlt1.gt.xmltlim.and.xmlt1.lt.(24.-xmltlim).and.&
               abs(xmlt1-xmlt(j)).gt.xmltlim) then   ! big warping in mlt
             irm(j)=i-1
             exit LATITUDE
          endif
          
          ! Excessively long lines are considered open 
          if (dssa(npf1) > LengthMax) then
             irm(j)=i-1
             exit LATITUDE
          endif

          dss2=dssa(npf1)/2.      ! find the middle point
          !write(*,*) '!!! start1, iLat,iLon,npf1',i,j,npf1
          call locate1IM(dssa,npf1,dss2,im)
          !write(*,*) '!!! end1'
          im1=im
          if ((dssa(im+1)-dss2).lt.(dss2-dssa(im))) im1=im+1

          npf=n5           ! make sure B decreases to bba(im1) and rises
          dssm=0.
          do m=1,npf1
             if (m.lt.npf1) dssm=dssm+(dssa(m+1)-dssa(m))
             igood=-1
             if (m.eq.1.or.m.eq.im1) igood=1
             if (m.gt.1.and.m.lt.im1.and.bba(m).gt.bba(im1).and.&
                  bba(m).lt.bm1(npf)) igood=1     ! B should be decreasing
             if (m.gt.im1.and.bba(m).gt.bba(im1).and.&
                  bba(m).gt.bm1(npf)) igood=1     ! B should be increasing
             if (igood.eq.1) then
                npf=npf+1
                bm1(npf)=bba(m)
                rm(npf)=ra(m)
                if (m.lt.npf1) dss(npf)=dssm
                dssm=0.                      ! reset dssm
                if (m.eq.im1) im2=npf        ! new im1
             endif
          enddo

          do m=1,n5               ! Add n5 points below rc
             rm1=rc-m*dre
             rme=rm1*re  
             rm(n5-m+1)=rm1
             rm(npf+m)=rm1
             bm1(n5-m+1)=sqrt(4.-3.*rm1/xli(i))*xme/rme**3 !assume dipole
             bm1(npf+m)=bm1(n5-m+1)
          enddo

          npf=npf+n5
          n=npf-1  ! new no. of intervals from N to S hemisphere
          do m=1,n
             rs(m)=0.5*(rm(m+1)+rm(m))
             bs(m)=0.5*(bm1(m+1)+bm1(m))
          enddo
          do m=1,n
             if (m.le.n5.or.m.ge.(npf-n5)) then
                rl=rs(m)/xli(i)
                cost=sqrt(1.-rl*rl)
                dss(m)=dre*sqrt(3.*cost*cost+1.)/2./cost
             endif
          enddo

          !                                <--dss(m)-->
          !      rs(1)                         rs(m)                     rs(n)
          !      bs(1)                         bs(m)                     bs(n)
          !    |-------|-----|......|------|----------|----|.......|---|-------|
          !  rm(1)                       rm(m)                               rm(n+1)
          !  bm(1)                       bm(m)                               bm(n+1)


          ! Field line integration using Taylor expansion method
          call timing_start('rbe_taylor')
          sumBn(0:nTaylor)=0.
          sumhBn(0:nTaylor)=0.
          do m=im2-1,1,-1 !im2 = middle of field line
             ! search for the southern conjugate point
             n8=npf
             SEARCH: do ii=m,n
                if (bm1(ii+1).ge.bm1(m)) then
                   n8=ii+1
                   EXIT SEARCH
                endif
             enddo SEARCH
             n7=n8-1


             ! field line integration at the northern hemisphere
             bs_n=1.
             do iTaylor=0,nTaylor
                bsndss=bs_n*dss(m)
                ! Sum_m=1^im2-1 (ds*B^n)
                sumBn(iTaylor)=sumBn(iTaylor)+bsndss
                sumhBn(iTaylor)=sumhBn(iTaylor)+hden(rs(m))*bsndss
                bs_n=bs_n*bs(m)
             enddo

             ! field line integration at the southern hemisphere
             m0=m+1
             if (m.lt.(im2-1)) m0=n70
             do mir=m0,n7-1
                bs_n=1.
                do iTaylor=0,nTaylor
                   bsndss=bs_n*dss(mir)
                   ! Sum_m=im2^n7-1 (ds*B^n)
                   sumBn(iTaylor)=sumBn(iTaylor)+bsndss
                   sumhBn(iTaylor)=sumhBn(iTaylor)+hden(rs(mir))*bsndss
                   bs_n=bs_n*bs(mir)
                enddo
             enddo

             ! add the partial segment near the southern conjugate point
             dssp=dss(n7)*(bm1(m)-bm1(n7))/(bm1(n8)-bm1(n7)) ! partial ds
             bsp =0.5*(bm1(m)+bm1(n7))
             bs_n=1.
             do iTaylor=0,nTaylor
                bsndss=bs_n*dssp
                ! Sum_m=1^n7-1 (ds*B^n) + correction
                BnI(iTaylor,m)=sumBn(iTaylor)+bsndss
                hBnI(iTaylor,m)=sumhBn(iTaylor)+hden(rs(n7))*bsndss
                bs_n=bs_n*bsp
             enddo

             n70=n7
          enddo

          ! Set up arrarys: si3, tya3, h3
          do m=1,im2-1
             ss=0.
             ss1=0.
             ss2=0.
             bm_n=1.
             do iTaylor=0,nTaylor
                BnIbm_n=BnI(iTaylor,m)/bm_n

                ! ss = int_s(1)^s(m) sqrt(B(m)-B(s)) ds
                ss=ss+a_I(iTaylor)*BnIbm_n
                ! ss1 = int_s(1)^s(m) 1/sqrt(B(m)-B(s)) ds              
                ss1=ss1+b_I(iTaylor)*BnIbm_n
                ! ss2 = int_s(1)^s(m) n_H(s)/sqrt(B(m)-B(s)) ds              
                ss2=ss2+b_I(iTaylor)*hBnI(iTaylor,m)/bm_n
                bm_n=bm_n*bm1(m)
             enddo
             si3(m)=re*sqrt(bm1(m))*ss
             tya3(m)=ss1/ro1/2.
             h3(m)=ss2/ss1
          enddo

          call timing_stop('rbe_taylor')

          si3(im2)=0.               ! equatorially mirroring
          tya3(im2)=tya3(im2-1)
          h3(im2)=hden(rm(im2))

          ! Calculate y, rmir (dist. of mirror point), T(y), bounced average [H]
          do m=1,ik
             !write(*,*) '!!! iLat,iLon',i,j
             sim=si(m)                 ! get Bm @ given K & location
             call lintpIM(si3,bm1,im2,sim,bmmx)
             bm(i,j,m)=bmmx
             sinA(i,j,m)=sqrt(bo(i,j)/bmmx)
             if (sinA(i,j,m).gt.1.) sinA(i,j,m)=1.
             call lintpIM(si3,rm,im2,sim,rmm)
             rmir(i,j,m)=rmm
             call lintpIM(si3,tya3,im2,sim,tya33)
             tya(i,j,m)=tya33
             call lintpIM(si3,h3,im2,sim,h33)
             Have(i,j,m)=h33  ! bounce-ave [H]
          enddo

       enddo LATITUDE                              ! end of i loop
    enddo LONGITUDE                              ! end of j loop

    ! Fill in volumes and coordinates for open (?) field lines
    do j = 1, ip
       do i=irm(j)+1,ir
          volume(i,j)=volume(irm(j),j)
          xo(i,j)=xo(irm(j),j)
          yo(i,j)=yo(irm(j),j)
       enddo
    end do

    call timing_stop('rbe_trace')

    ! Peridic boundary condition
    !do i=1,ir        
    !   volume(i,ip+1)=volume(i,1)
    !   xo(i,ip+1)=xo(i,1)
    !   yo(i,ip+1)=yo(i,1)
    !   gridoc(i,ip+1)=gridoc(i,1)
    !enddo

    !.....Calculate p, E, v, et (mc^2) at grid boundaries, and lifetime for
    !     loss cone particles
    do n=1,nspec
       xmass=1.673e-27*amu_I(n)
       c2mo=c*c*xmass
       c4mo2=c2mo*c2mo
       do j=1,ip
          do i=1,irm(j)
             ro2=2.*ro(i,j)*re
             do m=1,ik
                pp1=sqrt(2.*xmass*bm(i,j,m))
                tcone1=ro2*tya(i,j,m) 
                do k=1,iw
                   pijkm=pp1*sqrt(w(k))
                   pc=pijkm*c
                   c2m=sqrt(pc*pc+c4mo2)
                   e=c2m-c2mo                 ! E in J
                   ekev(i,j,k,m)=e/1000./q    ! E in keV
                   gamma(i,j,k,m)=c2m/c2mo
                   pp(n,i,j,k,m)=pijkm  
                   vel(n,i,j,k,m)=pc*c/c2m
                   alscone(n,i,j,k,m)=1.
                   if (rmir(i,j,m).le.rc) then
                      tcone2=tcone1/vel(n,i,j,k,m)      ! Tbounce/2
                      x=dt/tcone2
                      alscone(n,i,j,k,m)=0.
                      if (x.le.80.) alscone(n,i,j,k,m)=exp(-x)
                   endif
                enddo
                
             enddo
          enddo
       enddo
    enddo

    ! Reduce irm by 1 for convenience in calculating vl at i=irm+0.5
    do j=1,ip
       irm(j)=irm(j)-1
    enddo
    
    ! Find iba
    if (UseEllipse) then
       R_24=rb                 ! boundary distance at midnight
       do j=1,ip
          imax=irm(j)
          xmltr=xmlto(imax,j)*pi/12.
          xBoundary(j)=-ro(imax,j)*cos(xmltr)
       enddo
       R_12=0.95*maxval(xBoundary)    ! boundary distance at noon
       MajorAxis=0.5*(R_12+R_24)      ! major axis
       MinorAxis=min(R_12,R_24)       ! minor axis
       xCenter=0.5*(R_12-R_24)        ! center on x-axis
       MajorAxis2=MajorAxis*MajorAxis
       MinorAxis2=MinorAxis*MinorAxis
       do j=1,ip
          find_ib: do i=irm(j),1,-1
             xmltr=xmlto(i,j)*pi/12.
             sin2=sin(xmltr)*sin(xmltr)
             Req2=ro(i,j)*ro(i,j)
             xo1=-ro(i,j)*cos(xmltr)
             xc=xo1-xCenter
             rell2=MinorAxis2*(1.-xc*xc/MajorAxis2)/sin2 ! r^2 onellipse at x=xc
             if (Req2.le.rell2) then
                iba(j)=i
                exit find_ib
             endif
          enddo find_ib
       enddo
    else
       !use circle
       do j=1,ip
          do i=1,irm(j)
             x1(i)=ro(i,j)
          enddo
          call locate1IM(x1,irm(j),rb,ib)
          iba(j)=ib
       enddo
    endif

    ! Find iw2(m)
    do m=1,ik
       iw2(m)=iw
       find_iw2: do k=1,iw
          if (ekev(irm(1),1,k,m).gt.energy(neng)) then
             iw2(m)=k
             exit find_iw2
          endif
       enddo find_iw2
    enddo

  end subroutine fieldpara

! uncomment when T04 Tracing fixed
!
!  !*******************************************************************************
!  !                             TsyParmod
!  !  Rountine calculates the parmod in Tsyganenko model.
!  !*******************************************************************************
!  subroutine TsyParmod(t,tf,tsw,xnswa,vswa,nsw,tdst,dsth,ndst,timf,&
!       byw,bzw,nimf,xmp,imod,parmod)
!    use rbe_cread2, ONLY: ismo
!    real tsw(nsw),xnswa(nsw),vswa(nsw),tdst(ndst),dsth(ndst),timf(nimf),&
!         byw(nimf),bzw(nimf),parmod(10),w04(6),rr(6),xlamb(6),beta(6),gamm(6)
!
!    ! Parameters for T04_S model
!    data rr/0.39,0.7,0.031,0.58,1.15,0.88/     ! relaxation rate in hour^-1
!    data xlamb/0.39,0.46,0.39,0.42,0.41,1.29/
!    data beta/0.8,0.18,2.32,1.25,1.6,2.4/
!    data gamm/0.87,0.67,1.32,1.29,0.69,0.53/
!
!    parmod(1:10)=0.             ! initial values
!
!    tsmo=1.7                    ! running window time normalized by tf
!    if (ismo.eq.0) tsmo=0.
!    tti=t-tsmo*tf
!    ttf=t+(1.+tsmo)*tf
!
!    !  parmod(1): solar wind pressure in nPa
!    call locate1IM(tsw,nsw,tti,j1)
!    call locate1IM(tsw,nsw,ttf,j_2)
!    j2=j_2+1
!    if (j1.eq.0) j1=1
!    if (j2.gt.nsw) j2=nsw
!    xnsw=0.
!    vsw=0.
!    do j=j1,j2
!       xnsw=xnsw+xnswa(j)/(j2-j1+1)        ! find average SW ./. tti and ttf
!       vsw=vsw+vswa(j)/(j2-j1+1)
!    enddo
!    parmod(1)=xmp*xnsw*1.e6*vsw*vsw*1.e6/1.e-9    ! Pdyn in nPa
!    if (parmod(1).lt.2.0) parmod(1)=2.0  ! set min parmod(1) to 2.0
!
!    !  parmod(2): Dst
!    call locate1IM(tdst,ndst,tti,j1)
!    call locate1IM(tdst,ndst,ttf,j_2)
!    j2=j_2+1
!    if (j1.eq.0) j1=1
!    if (j2.gt.ndst) j2=ndst
!    dst=0.
!    do j=j1,j2
!       dst=dst+dsth(j)/(j2-j1+1)        ! find average dst ./. tti and ttf
!    enddo
!    parmod(2)=dst
!
!    !  parmod(3:4): IMF By, Bz in nT
!    call locate1IM(timf,nimf,tti,j1)
!    call locate1IM(timf,nimf,ttf,j_2)
!    j2=j_2+1
!    if (j1.eq.0) j1=1
!    if (j2.gt.nimf) j2=nimf
!    byimf=0.
!    bzimf=0.
!    do j=j1,j2
!       byimf=byimf+byw(j)/(j2-j1+1)    ! find average IMF ./. tti and ttf
!       bzimf=bzimf+bzw(j)/(j2-j1+1)
!    enddo
!    parmod(3)=byimf
!    parmod(4)=bzimf
!
!    !  Limit the values of parmod(1:4) in t96 model (imod=1)
!    if (imod.eq.1) then
!       if (parmod(1).gt.10.) parmod(1)=10.              
!       if (parmod(2).lt.-100.) parmod(2)=-100.
!       if (parmod(2).gt.20.) parmod(2)=20.
!       if (parmod(3).lt.-10.) parmod(3)=-10.
!       if (parmod(3).gt.10.) parmod(3)=10.
!       if (parmod(4).lt.-10.) parmod(4)=-10.
!       if (parmod(4).gt.10.) parmod(4)=10.
!    endif
!
!    !  parmod(5:10) for t04_s: w04(1:6) defined in Tsyganenko and Sitnov, 2005
!    if (imod.eq.2) then
!       tti=t-100.*3600.              ! 100 hours before t
!       call locate1IM(tsw,nsw,tti,j1)    
!       if (j1.eq.0) j1=1
!       call locate1IM(tsw,nsw,t,j2)
!       if (j2.eq.0) j2=1
!       w04(1:6)=0.
!       do j=j1,j2      ! average over preceding hours
!          tk=tsw(j)
!          tdiff=(tk-t)/3600.   ! time difference in hour
!          if (tk.lt.timf(1)) tk=timf(1)
!          if (tk.gt.timf(nimf)) tk=timf(nimf)        
!          call lintp(timf,bzw,nimf,tk,bz1)
!          if (bz1.lt.0.) Bs1=-bz1
!          if (bz1.ge.0.) goto 1       ! +ve Bz, no contribution to w04
!          xnsw_n=xnswa(j)/5.          ! normalized sw density
!          vsw_n=vswa(j)/400.          ! normalized sw velocity
!          Bs_n=Bs1/5.                 ! normalized Bs
!          do m=1,6
!             ert=exp(rr(m)*tdiff)
!             Sk=xnsw_n**xlamb(m)*vsw_n**beta(m)*Bs_n**gamm(m)
!             w04(m)=w04(m)+Sk*ert
!          enddo
!1         continue
!       enddo
!       del_t=(tsw(j2)-tsw(j1))/(j2-j1+1)/3600.      ! delta t in hour
!       if (del_t.le.0.) del_t=1./12.                
!       do m=1,6
!          w04(m)=w04(m)*rr(m)*del_t
!          parmod(m+4)=w04(m)
!       enddo
!    endif        ! end of if (imod.eq.2) 
!
!    ! Set limit to parmod for t04_s
!    if (imod.eq.2) then
!       if (parmod(2).lt.-300.) parmod(2)=-300.     ! Dst
!       if (parmod(8).gt.25.) parmod(8)=25.         ! partial ring current
!       if (parmod(10).gt.100) parmod(10)=100.      ! region 2 current
!    endif
!  end subroutine TsyParmod
!  
!  !=============================================================================
!  subroutine tsy_trace(iLat,rlim,re,rc,xlati1,phi1,t,ps,parmod,imod,np, &
!       npf1,dssa,bba,volume1,ro1,xmlt1,bo1,ra)
!    ! Routine does field line tracing in Tsyganenko field.For a given xlati1 and
!    ! phi1, output distance from the ionosphere,magnetic field, flux tube volume
!    ! per unit flux and equatorial crossing point.
!    ! npf1 is the number of point along that field line.npf1=0 if this isan open
!    ! field line.
!    !
!    ! Input: re,rc,xlati1,phi1,t,ps,parmod,imod,np
!    ! Output: npf1,dssa,bba,volume1,ro1,xmlt1,bo1    ! bba, bo1 in Tesla 
!
!    use ModSort, ONLY: sort_quick
!    use ModNumConst, ONLY: cPi
!    implicit none
!    external tsyndipoleSM
!
!    integer, parameter :: nd=3
!    integer imod,np,npf1,i,j,n,ii,iopt,ind(np)
!    integer, intent(in) :: iLat  
!    real, intent(out)   :: ra(np)
!    integer:: ieq
!    real rlim,re,rc,xlati1,phi1,t,ps,parmod(10),dssa(np),bba(np),volume1,ro1,&
!         xmlt1,bo1
!    real xa(np),ya(np),za(np),x0(3),xend(3),f(3),t0,tend,h,h1,aza(np)
!    real dir,pas,xwrk(4,nd),b_mid,dss(np),ss,yint(np)
!
!    iopt=1               ! dummy variable for tsy models
!    !  rlim=20.
!    dir=-1.              ! start fieldline tracing from Northern hemisphere
!    pas=0.1              ! fieldline tracing step in RE
!    h=pas*dir
!
!    ! initial
!    x0(1)=rc*cos(xlati1)*cos(phi1)  ! sm x
!    x0(2)=rc*cos(xlati1)*sin(phi1)  ! sm y
!    x0(3)=rc*sin(xlati1)            ! sm z
!    t0=0.
!    npf1=1
!    call tsyndipoleSM(imod,iopt,parmod,ps,t,x0(1),x0(2),x0(3),f(1),f(2),f(3))
!    bba(1)=sqrt(f(1)*f(1)+f(2)*f(2)+f(3)*f(3))*1.e-9   ! B in T
!    xa(1)=x0(1)
!    ya(1)=x0(2)
!    za(1)=x0(3)
!    ra(1)=sqrt(xa(1)*xa(1)+ya(1)*ya(1)+za(1)*za(1))
!    dssa(1)=0.
!
!    ! start tracing
!    trace: do
!       call rk4(tsyndipoleSM,imod,iopt,parmod,ps,t,t0, &
!            h,x0,xend,xwrk,nd,f,tend)
!       npf1=npf1+1
!       ra(npf1)=sqrt(xend(1)*xend(1)+xend(2)*xend(2)+xend(3)*xend(3))
!       xa(npf1)=xend(1)
!       ya(npf1)=xend(2)
!       za(npf1)=xend(3)
!       bba(npf1)=sqrt(f(1)*f(1)+f(2)*f(2)+f(3)*f(3))*1.e-9    ! B in T
!       dssa(npf1)=dssa(npf1-1)+abs(h)
!
!       if (ra(npf1).gt.rlim.or.npf1.gt.np) then
!          npf1=0              ! open field line
!          !        write(*,*) 'iLat,Lat,mlt',iLat,xlati1*180/3.14, xmlt1
!          !        exit trace              
!          return
!       endif
!
!       if (ra(npf1).le.rc) then               ! at south hemisphere
!          ! reduce step size such that ra(npf1) is at rc
!          h1=(ra(npf1-1)-rc)/(ra(npf1-1)-ra(npf1))
!          ra(npf1)=rc
!          xa(npf1)=xa(npf1-1)+(xend(1)-xa(npf1-1))*h1
!          ya(npf1)=ya(npf1-1)+(xend(2)-ya(npf1-1))*h1
!          za(npf1)=za(npf1-1)+(xend(3)-za(npf1-1))*h1
!          call tsyndipoleSM(imod,iopt,parmod,ps,t, &
!               xa(npf1),ya(npf1),za(npf1),f(1),f(2),f(3))
!          bba(npf1)=sqrt(f(1)*f(1)+f(2)*f(2)+f(3)*f(3))*1.e-9   ! B in T
!          dssa(npf1)=dssa(npf1-1)+abs(h)*h1
!
!          ! Calculate the flux tube volume per magnetic flux (volume1)
!          n=npf1-1 ! n = no. of intervals from N(rc) to S(rc) hemisphere
!          do ii=1,n
!             b_mid=0.5*(bba(ii)+bba(ii+1))
!             dss(ii)=dssa(ii+1)-dssa(ii)
!             yint(ii)=1./b_mid
!          enddo
!          call closed(n,yint,dss,ss)  ! use closed form
!          if (iLat.ge.1.and.iLat.le.ir) volume1=ss*re   ! volume / flux
!          exit trace      ! finish tracing this field line
!       endif
!
!       x0(1:nd)=xend(1:nd)
!       t0=tend
!    enddo trace           ! continue tracing along this field line
!
!    ! find the equatorial crossing point
!    aza(1:npf1)=abs(za(1:npf1))
!
!    call sort_quick(npf1,aza,ind)    ! find the equatorial crossing point
!    ieq=ind(1)
!    ro1=ra(ieq)
!    xmlt1=atan2(-ya(ieq),-xa(ieq))*12./cPi   ! mlt in hr
!    if (xmlt1.lt.0.) xmlt1=xmlt1+24.
!    bo1=bba(ieq)
!
!  end subroutine tsy_trace

  !=============================================================================
  subroutine mhd_trace_IM (Lat,Lon,re,iLat,iLon,np, &
       nAlt,FieldLength_I,Bfield_I,volume1,ro1,xmlt1,bo1,RadialDist_I)

    use ModGmCrcm
    use ModNumConst, ONLY: cPi

    integer,intent(in)  :: iLat,iLon,np
    real,   intent(in)  :: Lat,Lon,re
    ! bba, bo1 in Tesla 
    real,   intent(out) :: RadialDist_I(np),FieldLength_I(np),Bfield_I(np),&
         volume1,ro1,xmlt1,bo1
    integer,intent(out) :: nAlt

    ! Number of points covering the gap from 1 Re to rBody 
    integer, parameter :: MinAlt = 25

    integer, parameter :: nd=3
    integer ::i,j,n,ii,iopt,ind(np),iPoint,iAlt
    integer, parameter :: I_=1,S_=2,R_=3,B_=4
    real,    parameter :: LatMin = .886  ! 50.7degrees, below this, fieldline 
    ! extends below 2.5Re

    real xa(np),ya(np),za(np),x0(3),xend(3),f(3),t0,tend,h,h1,aza(np)
    real dir,pas,xwrk(4,nd),rlim,Bmid,dss(np),ss,yint(np)
    character(len=*), parameter :: NameSub='mhd_trace_IM'
    Logical IsFoundLine,UseDipole
    !---------------------------------------------------------------------------

    ! Put BufferLine_VI indexed by line number into StateLine_CIIV

    ! Start after MinAlt points inside rBody
    iAlt = MinAlt
    IsFoundLine=.false.
    UseDipole  =.false.
    FieldTrace: do iPoint = 1,nPoint

!       if(iLon==5 .and. iLat==22) &
!            write(*,*)'iLineIndex_II(iLon,iLat),int(StateLine_VI(I_,iPoint))'&
!            ,iLineIndex_II(iLon,iLat),int(StateLine_VI(I_,iPoint))

       if (iLineIndex_II(iLon,iLat) == int(StateLine_VI(I_,iPoint)))then
          !if(iLon==5)write(*,*)'!!! Found line for iLat,iLon',iLat,iLon
          !when line index found, fill in output arrays
          Bfield_I(iAlt)     = StateLine_VI(B_,iPoint)
          FieldLength_I(iAlt)= StateLine_VI(S_,iPoint)
          RadialDist_I(iAlt) = StateLine_VI(R_,iPoint)

          iAlt = iAlt+1
          IsFoundLine=.true.
       elseif (iLineIndex_II(iLon,iLat) /= int(StateLine_VI(I_,iPoint)) &
            .and. IsFoundLine) then
          exit FieldTrace 
       endif
    end do FieldTrace

    ! Add points for the other end inside rBody
    nAlt = iAlt-1 + MinAlt

    ! Fill in points below rBody
    if (IsFoundLine) &
         call trace_dipoleIM(Re,Lat,nAlt,MinAlt,FieldLength_I,&
         Bfield_I,RadialDist_I,Ro1)

    ! Field lines fully inside rBody
    if (Lat <= Latmin) then
       nAlt=2*MinAlt
       call trace_dipoleIM(Re,Lat,nAlt,MinAlt,FieldLength_I,&
            Bfield_I,RadialDist_I,Ro1)
       xmlt1= mod((Lon)*12./cPi+12.0,24.0)   ! mlt in hr           
       if (xmlt1.lt.0.) xmlt1=xmlt1+24.
       bo1=Bfield_I(nAlt/2)
       IsFoundLine = .true.
       UseDipole   = .true.
    end if

    !Check if FieldLine is open
    if (.not. IsFoundLine) then
       nAlt=0
       return
    endif

    if (.not. UseDipole) then
       ro1=sqrt(sum(StateIntegral_IIV(iLat,iLon,1:2)**2.0))
       !if (iLat==23 .and. iLon==3 )
       !write(*,*) 'Lat,Lon,iLat,iLon,iLineIndex_II(iLon,iLat)',Lat*180.0/cPi,Lon*180.0/cPi,&
       !     iLat,iLon,iLineIndex_II(iLon,iLat)
       !write(*,*) 'ro1,StateIntegral_IIV(iLat,iLon,1:2),maxval(RadialDist_I(1:nAlt))',&
       !     ro1,StateIntegral_IIV(iLat,iLon,1:2),maxval(RadialDist_I(1:nAlt))

       xmlt1=&
            (atan2(-StateIntegral_IIV(iLat,iLon,2),-StateIntegral_IIV(iLat,iLon,1))&
            )&
            *12./cPi   ! mlt in hr
       if (xmlt1 < 0.) xmlt1=xmlt1+24.
       bo1=StateIntegral_IIV(iLat,iLon,3)
    endif

    !Check that nAlt < np
    if (nAlt > np) then 
       !write(*,*) 'nAlt,np',nAlt,np
       !call CON_STOP('IM error: nAlt > np in mhd_trace_IM. Increase np and recompile.')
       !Treat line as open
       nAlt=0
       return
    endif
    
    ! Calculate the flux tube volume per magnetic flux (volume1)

    !write(*,*) '!!! iLat,iLon,nAlt',iLat,iLon,nAlt
    do ii=1,nAlt-1
       Bmid=0.5*(Bfield_I(ii)+Bfield_I(ii+1))
       dss(ii)=FieldLength_I(ii+1)-FieldLength_I(ii)
       yint(ii)=1./Bmid
    enddo
    call closed(nAlt-1,yint,dss,ss)  ! use closed form
    if (iLat >= 1 .and. iLat <= ir) volume1=ss*re   ! volume / flux
  end subroutine mhd_trace_IM
  
  !============================================================================
  subroutine lintpIM(xx,yy,n,x,y)
    !-----------------------------------------------------------------------
    !  Routine does 1-D interpolation.  xx must be increasing or decreasing
    !  monotonically.  x is between xx(j) and xx(j+1)
    integer :: n
    real xx(n),yy(n)
    integer :: ier = 0, i, j, jl, ju, jm 
    real    :: x, d, y
    !  Make sure xx is increasing or decreasing monotonically
    do i=2,n
       if (xx(n).gt.xx(1).and.xx(i).lt.xx(i-1)) then
          write(*,*) ' lintpIM: xx is not increasing monotonically '
          write(*,*) n,(xx(j),j=1,n)
          call CON_stop('IM ERROR')
       endif
       if (xx(n).lt.xx(1).and.xx(i).gt.xx(i-1)) then
          write(*,*) ' lintpIM: xx is not decreasing monotonically '
          write(*,*) 'i,xx(i),xx(i-1) ',i,xx(i-1),xx(i)
          write(*,*) n,(xx(j),j=1,n)
          call CON_stop('IM ERROR')
       endif
    enddo
    
    !  Set ier=1 if out of range
    if (xx(n).gt.xx(1)) then
       if (x.lt.xx(1).or.x.gt.xx(n)) ier=1
    else
       if (x.gt.xx(1).or.x.lt.xx(n)) ier=1
    endif
    if (ier.eq.1) then
       write(*,*) ' Error: ier.eq.1'
       write(*,*) ' x  ',x
       write(*,*) ' xx  ',xx
       call CON_stop('IM ERROR')
    endif
    !
    !    initialize lower and upper values
    !
    jl=1
    ju=n
    !
    !    if not dne compute a midpoint
    !
10  if(ju-jl.gt.1)then
       jm=(ju+jl)/2
       !
       !    now replace lower or upper limit
       !
       if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
       else
          ju=jm
       endif
       !
       !    try again
       !
       go to 10
    endif
    !
    !    this is j
    !
    j=jl      ! if x.le.xx(1) then j=1
    !                 if x.gt.x(j).and.x.le.x(j+1) then j=j
    !                 if x.gt.x(n) then j=n-1
    d=xx(j+1)-xx(j)
    y=(yy(j)*(xx(j+1)-x)+yy(j+1)*(x-xx(j)))/d
    
  end subroutine lintpIM
  

  !============================================================================
  
  function Hden(x)
    !--------------------------------------------------------------------------
    ! Chamberlain model of [H] fitted by an exponential function.
    ! The fit matches Rairden et al. [1986] for radial distance
    ! from 1.08(rexob) to 12 Re.
    
    implicit none
    
    real x,rexob,rr,ar,Hden
    
    rexob=1.08       ! 1.08 = exobase in Rairden et al. [1986]
    rr=x
    if (rr.lt.rexob) rr=rexob
    ar=log(rr/rexob)       ! rexob = exobase in Rairden et al. [1986]
    Hden=1.e6*exp(10.692-4.4431*ar**0.715831)      ! H density in m^-3
    
  end function Hden
  
  !***********************************************************************
  !                          closed
  ! Routine performs numerical integration using closed form 
  !
  !  S = y1dx1 + y2dx2 + .. + yidxi + .. + yn-1dxn-1 + yndxn
  !                 
  !  where yi is the value at the middle of dxi
  !***********************************************************************
  subroutine closed(n,y,dx,s)
    integer, intent(in) :: n
    real,    intent(in) :: y(n),dx(n)
    real,    intent(out):: s
    integer :: i
    !--------------------------------------------------------------------------
    s=0.           
    do i=1,n     
       s=s+y(i)*dx(i) 
    enddo
    
  end subroutine closed
  

end Module ModFieldTrace
