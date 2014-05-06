module ModUserSetup

  implicit none

  private

  public :: grid
  public :: field_init

contains

  subroutine grid
    use ModPar
    use ModSolar
    use ModPhysunits
    use ModDomain
    use ModGrid
    use ModBack
    use ModCoeff
    use ModBoundary
    use ModInterp,     ONLY:lint

    implicit none

    real :: rp,yp,re_p,pe_p,gacc_p
    real :: delta_domain(inmax), sr(inmax)
    real :: sigma_nd,solflux,fluxup,fluxdn,bvfq2

    integer :: i, j, k
    real x1min,x2min,x3min,x1max,x2max,x3max, &
         delx1,delx2,delx3
    real, dimension(:), allocatable :: vol1a, vol1b, &
         vol2a, vol2b, vol3a, vol3b
    real :: deltas,deltath, svarmean,areath
    !-------------------------------------------------------
    allocate(vol1a(in), vol1b(in), vol2a(jn), vol2b(jn), &
         vol3a(kn), vol3b(kn))
    !
    x1min=rmin/unit_l
    x1max=rmax/unit_l
    x2min=thmin
    x2max=thmax
    x3min=phmin
    x3max=phmax
    delx1=(x1max-x1min)/dble(inmax-5)
    delx3=(x3max-x3min)/dble(knmax-5)
    delx2=(x2max-x2min)/dble(jnmax-5)

    do i = 1, inmax
       xxa(i)=dble(i-is)*delx1+x1min
       dxxa(i)=delx1
       xxb(i)=xxa(i)+0.5D0*delx1
       dxxb(i)=delx1
       g2xxa(i)=xxa(i)
       g2xxb(i)=xxb(i)
       g31xxa(i)=xxa(i)
       g31xxb(i)=xxb(i)
    enddo

    do i = 1, in
       x1a(i)=xxa(i+myid1*(in-5))
       dx1a(i)=dxxa(i+myid1*(in-5))
       x1b(i)=xxb(i+myid1*(in-5))
       dx1b(i)=dxxb(i+myid1*(in-5))
    enddo

    do j = 1,jnmax
       yya(j)=dble(j-js)*delx2+x2min
       dyya(j)=delx2
       yyb(j)=yya(j)+0.5D0*delx2
       dyyb(j)=delx2
       g32yya   (j) = abs(sin ( yya(j) ))
       g32yyb   (j) = abs(sin ( yyb(j) ))
    enddo
    if(yya(js) .eq. 0.D0) then
       g32yya(js)=0.D0
       g32yya(jsm1)=g32yya(jsp1)
       g32yya(jsm2)=g32yya(jsp2)
       g32yyb(jsm1)=g32yyb(js)
       g32yyb(jsm2)=g32yyb(jsp1)
    endif
    if(yya(jnmax-2) .eq. pi) then
       g32yya(jnmax-2)=0.D0
       g32yya(jnmax-1)=g32yya(jnmax-3)
       g32yya(jnmax)=g32yya(jnmax-4)
       g32yyb(jnmax-2)=g32yyb(jnmax-3)
       g32yyb(jnmax-1)=g32yyb(jnmax-4)
    endif

    do j = 1, jn
       x2a(j)=yya(j+myid2*(jn-5))
       dx2a(j)=dyya(j+myid2*(jn-5))
       x2b(j)=yyb(j+myid2*(jn-5))
       dx2b(j)=dyyb(j+myid2*(jn-5))
    enddo

    do k = 1, kn
       x3a(k)=dble(k-ks)*delx3+x3min
       dx3a(k)=delx3
       x3b(k)=x3a(k)+0.5D0*delx3
       dx3b(k)=delx3
    enddo

    do i = 1, in
       g2a(i)=x1a(i)
       g31a(i)=x1a(i)
       dg2bd1(i)=1.D0
       dg31bd1(i)=1.D0
       vol1a(i)=x1a(i)**3/3.D0
       g2b(i)=x1b(i)
       g31b(i)=x1b(i)
       dg2ad1(i)=1.D0
       dg31ad1(i)=1.D0
       vol1b(i)=x1b(i)**3/3.D0
    enddo

    do j = 1, jn
       g32a   (j) = abs(sin ( x2a(j) ))
       vol2a  (j) =-cos ( x2a(j) )
       g32b   (j) = abs(sin ( x2b(j) ))
       vol2b  (j) =-cos ( x2b(j) )
    enddo

    if(myid2 .eq. 0) then

       if(x2a(js) .eq. 0.D0) then
          g32a(js)=0.D0
          g32a(jsm1)=g32a(jsp1)
          g32a(jsm2)=g32a(jsp2)
          g32b(jsm1)=g32b(js)
          g32b(jsm2)=g32b(jsp1)
       endif

    endif

    if(myid2 .eq. nproc2-1) then

       if(x2a(jep1) .eq. pi) then
          g32a(jep1)=0.D0
          g32a(jep2)=g32a(je)
          g32a(jep3)=g32a(jem1)
          g32b(jep1)=g32b(je)
          g32b(jep2)=g32b(jem1)
       endif

    endif

    do k = 1, kn
       vol3a(k) = x3a(k)
       vol3b(k) = x3b(k)
    enddo

    do j = 1, jn
       dg32bd2(j) = cos(x2a(j))
       dg32ad2(j) = cos(x2b(j))
    enddo

    if(myid2 .eq. 0) then

       if(x2a(js) .eq. 0.D0) then
          dg32bd2(js) = 0.D0
          dg32bd2(jsm1) = - dg32bd2(jsp1)
          dg32bd2(jsm2) = - dg32bd2(jsp2)
          dg32ad2(jsm1) = - dg32ad2(js)
          dg32ad2(jsm2) = - dg32ad2(jsp1)
       endif

    endif

    if(myid2 .eq. nproc2-1) then

       if(x2a(jep1) .eq. pi) then
          dg32bd2(jep1) = 0.D0
          dg32bd2(jep2) = - dg32bd2(je)
          dg32bd2(jep3) = - dg32bd2(jem1)
          dg32ad2(jep1) = - dg32ad2(je)
          dg32ad2(jep2) = - dg32ad2(jem1)
       endif

    endif


    do i = ism2, iep2
       dvl1a(i)=x1b(i)*x1b(i)*dx1a(i)
    enddo
    do i = ism2, iep3
       dvl1b(i)=x1a(i)*x1a(i)*dx1b(i)
    enddo
    !     dvl1b not used here

    do j = jsm2, jep2
       dvl2a(j) = sin(x2b(j))*dx2a(j)
    enddo
    do j = jsm2, jep3
       dvl2b(j) = sin(x2a(j))*dx2b(j)
    enddo
    !     dvl2b not used here

    if(myid2 .eq. 0) then

       if(x2a(js) .eq. 0.D0) then
          dvl2a(jsm1) = dvl2a(js)
          dvl2a(jsm2) = dvl2a(jsp1)
          dvl2b(js) = 0.D0
          dvl2b(jsm1) = dvl2b(jsp1)
          dvl2b(jsm2) = dvl2b(jsp2)
       endif

    endif

    if(myid2 .eq. nproc2-1) then

       if(x2a(jep1) .eq. pi) then
          dvl2a(jep1) = dvl2a(je)
          dvl2a(jep2) = dvl2a(jem1)
          dvl2b(jep1) = 0.D0
          dvl2b(jep2) = dvl2b(je)
          dvl2b(jep3) = dvl2b(jem1)
       endif

    endif


    do k=ksm2,kep2
       dvl3a(k) = dx3a(k)
    enddo
    do k=ksm2,kep3
       dvl3b(k) = dx3b(k)
    enddo

    do i = ism2, iep3
       dx1bi (i) = 1.D0 / ( dx1b (i) + tiny )
       g2ai  (i) = 1.D0 / ( g2a  (i) + tiny )
       g31ai (i) = 1.D0 / ( g31a (i) + tiny )
       dvl1bi(i) = 1.D0 / ( dvl1b(i) + tiny )
    enddo
    do i = ism2, iep2
       dx1ai (i) = 1.D0 / ( dx1a (i) + tiny )
       g2bi  (i) = 1.D0 / ( g2b  (i) + tiny )
       g31bi (i) = 1.D0 / ( g31b (i) + tiny )
       dvl1ai(i) = 1.D0 / ( dvl1a(i) + tiny )
    enddo

    do j=jsm2,jep3
       dx2bi (j) = 1.D0 / ( dx2b (j) + tiny )
       g32ai (j) = 1.D0 / ( g32a (j) + tiny )
       dvl2bi(j) = 1.D0 / ( dvl2b(j) + tiny )
    enddo
    do j=jsm2,jep2
       dx2ai (j) = 1.D0 / ( dx2a (j) + tiny )
       g32bi (j) = 1.D0 / ( g32b (j) + tiny )
       dvl2ai(j) = 1.D0 / ( dvl2a(j) + tiny )
    enddo
    !     SAFEGUARD

    if(myid2 .eq. 0) then

       if(x2a(js) .eq. 0.D0) then
          g32ai(js) = 0.D0
          dvl2bi(js) = 0.D0
       endif

    endif

    if(myid2 .eq. nproc2-1) then

       if(x2a(jep1) .eq. pi) then
          g32ai(jep1)=0.D0
          dvl2bi(jep1)=0.D0
       endif

    endif


    do k=ksm2,kep3
       dx3bi (k) = 1.D0 / ( dx3b (k) + tiny )
       dvl3bi(k) = 1.D0 / ( dvl3b(k) + tiny )
    enddo
    do k=ksm2,kep2
       dx3ai (k) = 1.D0 / ( dx3a (k) + tiny )
       dvl3ai(k) = 1.D0 / ( dvl3a(k) + tiny )
    enddo
    !
    !------------------------------------------------------
    !     background stratification
    !------------------------------------------------------
    !
    !     specify d0, gacc_str=gacc0/gacc_czb,
    !             gamma, temp0 (only lower rboundary required),
    !             kappa, delta_domain
    !
    do i=1,inmax
       rp=xxb(i)*unit_l
       call lint(rtab,re,nptjcd,rp,re_p)
       d0(i)=re_p/unit_d
       call lint(rtab,gacc,nptjcd,rp,gacc_p)
       gacc_str(i)=gacc_p/gacc_czb
       call lint(rtab,gamaad,nptjcd,rp,yp)
       gamma(i)=yp
       call lint(rtab,te,nptjcd,rp,yp)
       temp0(i)=yp/unit_temp
       call lint(rtab,ross_kappa,nptjcd,rp,yp)
       kappa(i)=yp*unit_d*unit_l
       call lint(rtab,delta,nptjcd,rp,yp)
       delta_domain(i)=0.D0
    enddo
    !
    !     compute temp0 using temp0(is) and dtemp0/dr = -(1-1/gamma)*gacc_str
    !     in the convection zone
    !
    temp0(is)=1.D0-(1.D0-1.D0/gamma(is))*gacc_str(is) &
         *dxxa(is)*0.5D0
    do i=is+1,inmax
       temp0(i)=temp0(i-1) &
            -(1.D0-1.D0/gamma(i-1))*gacc_str(i-1)*dxxa(i-1)*0.5D0 &
            -(1.D0-1.D0/gamma(i))*gacc_str(i)*dxxa(i)*0.5D0
    enddo
    !
    !     compute ovdxxb, hp, ovhp and fact
    !
    do i=1,inmax
       ovdxxb(i)=1.D0/dxxb(i)
       hp(i)=temp0(i)/gacc_str(i)
       ovhp(i)=1.D0/hp(i)
       fact(i)=1.D0/d0(i)
    enddo
    !
    !     compute sr from delta_domain, sr at lower boundary set to 0
    !
    sr(is)=-1.D0/(1.D0-1.D0/gamma(is))*(hp_czb*gacc_czb/unit_v**2) &
         *delta_domain(is)*ovhp(is)*dxxa(is)*0.5D0
    do i=is+1,inmax
       sr(i)=sr(i-1) &
            -1.D0/(1.D0-1.D0/gamma(i-1))*(hp_czb*gacc_czb/unit_v**2) &
            *delta_domain(i-1)*ovhp(i-1)*dxxa(i-1)*0.5D0 &
            -1.D0/(1.D0-1.D0/gamma(i))*(hp_czb*gacc_czb/unit_v**2) &
            *delta_domain(i)*ovhp(i)*dxxa(i)*0.5D0
    enddo
    do i=is-1,1,-1
       sr(i)=sr(i+1) &
            +1.D0/(1.D0-1.D0/gamma(i+1))*(hp_czb*gacc_czb/unit_v**2) &
            *delta_domain(i+1)*ovhp(i+1)*dxxa(i+1)*0.5D0 &
            +1.D0/(1.D0-1.D0/gamma(i))*(hp_czb*gacc_czb/unit_v**2) &
            *delta_domain(i)*ovhp(i)*dxxa(i)*0.5D0
    enddo
    !
    !     specify enforced entropy variation siib1 at the lower boundary
    deltas=0.2D0
    deltath=(thmax-thmin)*0.5D0
    svarmean=0.D0
    areath=0.D0
    do j=js,jnmax-3
       svarmean=svarmean &
            - deltas*cos((yyb(j)-pi*0.5D0)/deltath*pi) &
            *g32yyb(j)*dyya(j)
       areath=areath+g32yyb(j)*dyya(j)
    enddo
    svarmean=svarmean/areath
    do k=ks,ke
       do j=js,je
          siib1(j,k)=-svarmean &
               -deltas*cos((x2b(j)-pi*0.5D0)/deltath*pi)
       enddo
    enddo
    !
    !     compute radiative diffusive heat flux: radflux
    !     specify radiative heating: heatrad
    sigma_nd=sigma*unit_temp**4/unit_d/unit_v**3
    solflux=lsol/4.D0/pi/(unit_d*unit_v**3*unit_l**2)
    do i=3,inmax-3
       fluxdn=-16.D0*sigma_nd/(3.D0*0.5D0*(kappa(i-1)+kappa(i))) &
            *(0.5D0*(temp0(i-1)+temp0(i)))**3 &
            /(0.5D0*(d0(i-1)+d0(i))) &
            *(temp0(i)-temp0(i-1))/dxxb(i) &
            *xxa(i)*xxa(i)
       if(i .eq. 3) then
          fluxdn=solflux
       endif
       fluxup=-16.D0*sigma_nd/(3.D0*0.5D0*(kappa(i)+kappa(i+1))) &
            *(0.5D0*(temp0(i)+temp0(i+1)))**3 &
            /(0.5D0*(d0(i)+d0(i+1))) &
            *(temp0(i+1)-temp0(i))/dxxb(i+1) &
            *xxa(i+1)*xxa(i+1)
       heatrad(i)=-1.D0/xxb(i)/xxb(i) &
            *(fluxup-fluxdn)/dxxa(i)
       radflux(i)=fluxdn
    enddo
    heatrad(1)=0.D0
    heatrad(2)=0.D0
    heatrad(inmax-2)=0.D0
    heatrad(inmax-1)=0.D0
    heatrad(inmax)=0.D0
    radflux(1)=solflux
    radflux(2)=solflux
    radflux(inmax-2)=-16.D0*sigma_nd &
         /(3.D0*0.5D0*(kappa(inmax-3)+kappa(inmax-2))) &
         *(0.5D0*(temp0(inmax-3)+temp0(inmax-2)))**3 &
         /(0.5D0*(d0(inmax-3)+d0(inmax-2))) &
         *(temp0(inmax-2)-temp0(inmax-3))/dxxb(inmax-2) &
         *xxa(inmax-2)*xxa(inmax-2)
    radflux(inmax-1)=radflux(inmax-2)
    radflux(inmax)=radflux(inmax-2)
    !
    !     specify ovrth, ovrevar, and ovrmvar
    !     determine entropy profile s0 which is the initial <s>:
    !     here it is such that conduction with the specified
    !     ovrth balances the radiative heating
    !
    s0(1)=0.D0
    ovrth(1)=ovrt*sqrt(d0(inmax-3)/d0(1))
    ovrevar(1)=ovre*sqrt(d0(inmax-3)/d0(1))
    ovrmvar(1)=ovrm*sqrt(d0(inmax-3)/d0(1))
    do i=2,inmax
       ovrth(i)=ovrt*sqrt(d0(inmax-3)/d0(i))
       ovrevar(i)=ovre*sqrt(d0(inmax-3)/d0(i))
       ovrmvar(i)=ovrm*sqrt(d0(inmax-3)/d0(i))
       s0(i)=s0(i-1)-(solflux-radflux(i))/((xxa(i))**2 &
            *0.5D0*(d0(i-1)+d0(i)) &
            *0.5D0*(temp0(i-1)+temp0(i)) &
            *0.5D0*(ovrth(i-1)+ovrth(i)))*dxxb(i)
    enddo
    !
    !     compute ovbvfq given delta_domain
    !
    bvfq2=tiny
    do i=is,inmax-3
       bvfq2=max(bvfq2, -(hp_czb*gacc_czb/unit_v**2) &
            *delta_domain(i)*gacc_str(i)*ovhp(i))
    enddo
    ovbvfq=1.D0/sqrt(bvfq2)
    !     checking ovbvfq
    if(myid .eq. 0) then
       write(6,*) 'ovbvfq=',ovbvfq
    endif
    !
    !-----checking background
    if(myid .eq. 0) then
       open(unit=18,file='backstr.dat', &
            form='unformatted',access='stream')
       write(18) (xxb(i),i=1,inmax)
       write(18) (hp(i),i=1,inmax)
       write(18) (ovhp(i),i=1,inmax)
       write(18) (gamma(i),i=1,inmax)
       write(18) (d0(i),i=1,inmax)
       write(18) (gacc_str(i),i=1,inmax)
       write(18) (temp0(i),i=1,inmax)
       write(18) (s0(i),i=1,inmax)
       write(18) (fact(i),i=1,inmax)
       write(18) (ovrth(i),i=1,inmax)
       write(18) (kappa(i),i=1,inmax)
       write(18) (heatrad(i),i=1,inmax)
       write(18) (radflux(i),i=1,inmax)
       write(18) (ovrevar(i),i=1,inmax)
       write(18) (ovrmvar(i),i=1,inmax)
       write(18) (delta_domain(i),i=1,inmax)
       write(18) (sr(i),i=1,inmax)
       close(18)
    endif

    !     compute prbt, prtp, srbt, srtp for boundary conditions of p
    prbt=1.D0
    srbt=-2.D0/(1.D0/fact(is-1)+1.D0/fact(is)) &
         *dxxb(is)
    prtp=1.D0
    srtp=2.D0/(1.D0/fact(inmax-3)+1.D0/fact(inmax-2)) &
         *dxxb(inmax-2)

    deallocate(vol1a, vol1b, vol2a, vol2b, &
         vol3a, vol3b)
    deallocate(rtab, te, pe, re, &
         gacc, bvfsq, gamaad, gradad, &
         delta, gradrad, vc, ross_kappa)

  end subroutine grid

  subroutine field_init
    use ModPar
    use ModSolar
    use ModPhysunits
    use ModDomain
    use ModGrid
    use ModBack,    ONLY: d0, gacc_str, ovhp, gamma, bquench_k
    use ModField
    use ModBoundary
    use ModBval
    use ModIoFSAM

    implicit none

    integer :: i,j,k
    real amp
    real delthp,delthm,delphp,delphm,dwdrp,dwdrm
    real thc,thw,phc,phw,rc,rw,rr2,omg0
    real, allocatable :: afunc(:,:)
    real, allocatable :: wfunc(:,:,:)
    !------------------------------------------------------------
    allocate(afunc(in,jn))
    allocate(wfunc(in,jn,kn))
    !
    !-----------------------------------------------------------------
    ! specify boundary condition in the top and bottom boundarys
    !-----------------------------------------------------------------
    !
    niib=1
    noib=5
    nijb=1
    nojb=1
    !
    !-----------------------------------------------------------------
    ! initialize fields
    !-----------------------------------------------------------------
    !
    !     initialize magnetic field
    !
    bquench_k=0.D0
    amp=10.D0/unit_b*(pi-2.D0*thmin)/pi*rmax/unit_l &
         *sin(thmin)*2.D0/(1.D0-rmin/rmax)
    do j=jsm2,jep3
       do i=ism2,iep3
          afunc(i,j)=amp*(x1a(i)-rmin/unit_l) &
               *(1.D0-0.5D0*(x1a(i)-rmin/unit_l) &
               /(rmax/unit_l-rmin/unit_l)) &
               *sin((x2a(j)-thmin)/(pi-2.D0*thmin)*pi)
       enddo
    enddo
    do k=ksm2,kep2
       do j=jsm2,jep2
          do i=ism2,iep2
             b3(i,j,k)=0.D0
          enddo
       enddo
    enddo

    do k=ksm2,kep2
       do j=jsm2,jep2
          do i=ism2,iep2
             b1(i,j,k)=(afunc(i,j+1)-afunc(i,j))*dx2ai(j) &
                  *g2ai(i)**2/sin(x2b(j))
          enddo
       enddo
    enddo
    !     boundary of b1
    if(myid1 .eq. 0) then
       do k=ksm2,kep2
          do j=jsm2,jep2
             b1(is,j,k)=0.D0
          enddo
       enddo
    endif
    if(myid1 .eq. nproc1-1) then
       do k=ksm2,kep2
          do j=jsm2,jep2
             b1(iep2,j,k)=b1(ie,j,k)
          enddo
       enddo
    endif
    !
    do k=ksm2,kep2
       do j=jsm2,jep2
          do i=ism2,iep2
             b2(i,j,k)=-(afunc(i+1,j)-afunc(i,j))*dx1ai(i) &
                  *g2bi(i)*g32ai(j)
          enddo
       enddo
    enddo
    !     boundary of b2
    if(myid2 .eq. 0) then
       do k=ksm2,kep2
          do i=ism2,iep2
             b2(i,js,k)=0.D0
          enddo
       enddo
    endif
    if(myid2 .eq. nproc2-1) then
       do k=ksm2,kep2
          do i=ism2,iep2
             b2(i,jep1,k)=0.D0
          enddo
       enddo
    endif
    if(myid1 .eq. nproc1-1) then
       do k=ksm2,kep2
          do j=jsm2,jep2
             b2(iep1,j,k)=-b2(ie,j,k)
             b2(iep2,j,k)=-b2(iem1,j,k)
          enddo
       enddo
    endif
    !
    !     initialize velocity field
    !
    omg0=0.D0
    thc=0.5D0*(pi/6.D0+pi/2.D0)
    thw=(pi/2.D0-pi/6.D0)/6.D0
    phc=pi/2.D0/2.D0
    phw=pi/2.D0/6.D0
    rc=(r_sun-2.D9+r_czb)/2.D0/unit_l
    rw=(r_sun-2.D9-r_czb)/6.D0/unit_l
    do k=ksm2,kep2
       do j=jsm2,jep2
          do i=ism2,iep2
             rr2=(x1a(i)-rc)**2/rw**2 &
                  +(x2b(j)-thc)**2/thw**2 &
                  +(x3b(k)-phc)**2/phw**2
             if(rr2 .lt. 2.75D0**2) then
                wfunc(i,j,k)=1.D-1*exp(-rr2)
             else
                wfunc(i,j,k)=0.D0
             endif
          enddo
       enddo
    enddo

    do k=ks,ke
       do j=js,je
          do i=is,ie
             delthp=(wfunc(i,j+1,k)-wfunc(i,j,k))*dx2bi(j+1)
             delthm=(wfunc(i,j,k)-wfunc(i,j-1,k))*dx2bi(j)
             delphp=(wfunc(i,j,k+1)-wfunc(i,j,k))*dx3bi(k+1)
             delphm=(wfunc(i,j,k)-wfunc(i,j,k-1))*dx3bi(k)
             v1(i,j,k)=-((g32a(j+1)*delthp-g32a(j)*delthm) &
                  *g32bi(j)*dx2ai(j) &
                  +g32bi(j)*g32bi(j)*(delphp-delphm) &
                  *dx3ai(k)) &
                  /x1a(i)/x1a(i) &
                  /(0.5D0*(d0(i-1+myid1*(in-5))+d0(i+myid1*(in-5))))
          enddo
       enddo
    enddo

    do k=ks,ke
       do j=js,je
          do i=is,ie
             dwdrp=(wfunc(i+1,j,k)-wfunc(i,j,k))*dx1ai(i)
             dwdrm=(wfunc(i+1,j-1,k)-wfunc(i,j-1,k))*dx1ai(i)
             v2(i,j,k)=(dwdrp-dwdrm)*dx2bi(j)/x1b(i) &
                  /d0(i+myid1*(in-5))
             dwdrm=(wfunc(i+1,j,k-1)-wfunc(i,j,k-1))*dx1ai(i)
             v3(i,j,k)=(dwdrp-dwdrm)*dx3bi(k)/x1b(i)*g32bi(j) &
                  /d0(i+myid1*(in-5))
          enddo
       enddo
    enddo

    call bvalv

    !
    !     initialize entropy field
    !
    do k=ks,ke
       do j=js,je
          do i=is,ie
             s(i,j,k) = siib1(j,k)
          enddo
       enddo
    enddo
    call bvals
    !
    !     initialize pressure
    !
    do k=ksm2,kep2
       do j=jsm2,jep2
          do i=ism2,iep2
             p(i,j,k)=0.D0
          enddo
       enddo
    enddo

    deallocate(wfunc)
    deallocate(afunc)

  end subroutine field_init

end module ModUserSetup
