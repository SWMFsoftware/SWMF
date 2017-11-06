
*******************************************
*******************************************

!            transprt

*******************************************
*******************************************

      subroutine transprt (nfl,nll,p_crit)

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc'

      real prod(nz,nion),loss(nz,nion),lossr,
     .     phprodr(nz,nion),chrate(nz,nchem),
     .     chloss(nz,nion),chprod(nz,nion),relossr(nz,nion)
      real deni_old(nz,nion),te_old(nz),ti_old(nz,nion),vsi_old(nz,nion)
      real tvn(nz,nl)
      real nuin(nz,nion,nneut),
     .     nuij(nz,nion,nion),sumnuj(nz,nion)
      real vsin(nz,nion),vsidn(nz,nion),denin(nz,nion),cs(nz,nion)
      real ten(nz),tin(nz,nion)
      real p_crit(nnx-1)

! calculation of production and loss
!   phprodr: photo production rates
!   chrate:  chemical rates (ichem)
!   chloss:  chemical loss term
!   chprod:  chemical production term
!   relossr: recombination loss rates

! initialize tvn and gs and cfs (centrifugal force)

      do i = 1,nz
        tvn(i,nll) = 0.
        gs(i,nll)  = 0.
        cfs(i,nll)  = 0.
      enddo

      do i = 1,nz
        ne(i,nfl,nll)   = 1.
        te_old(i)       = te(i,nfl,nll)
        do j = nion1,nion2
          deni_old(i,j) = deni(i,nfl,nll,j)
          ne(i,nfl,nll) = ne(i,nfl,nll) + deni(i,nfl,nll,j)  
          ti_old(i,j)   = ti(i,nfl,nll,j)
          vsi_old(i,j)  = vsi(i,nfl,nll,j)
        enddo
      enddo

      call photprod ( phprodr,nfl,nll               ) ! calculates phprodr
      call chemrate ( chrate,nfl,nll                ) ! calculates chrate
      call chempl   ( chrate,chloss,chprod,nfl,nll  ) ! calcualtes chloss,chprod
      call recorate ( relossr,nfl,nll               ) ! calculates relossr

      do i = 1,nz 
        do j = nion1,nion2
          prod  (i,j) =  phprodr(i,j) * denn(i,nfl,nll,j)        
     .                   + chprod(i,j)
          lossr       =  relossr(i,j) * deni(i,nfl,nll,j) * 
     .                   ne(i,nfl,nll) 
     .                   + chloss(i,j)
          loss (i,j)  =  lossr / deni(i,nfl,nll,j)
        enddo

!     loss term for hydrogen and helium

        if ( alts(i,nfl,nll) .gt. pcrit*re ) then
          loss(i,pthp)  = loss(i,pthp)  + 1./decay_time 
          loss(i,pthep) = loss(i,pthep) + 1./decay_time 
! add loss term for O+ (JK)
!          loss(i,ptop) = loss(i,ptop) + 1./decay_time 
        endif

        gs(i,nll)   =  gzero * xrg(i,nfl,nll)
     .                 * ( re / (re + alts(i,nfl,nll)) ) ** 2

!        if (nll.eq.nl/2 .and. nfl.eq.44)
!     .   print *,'old s',i,gs(i,nll)

!        gs(i,nll)   = -gzero 
!     .                 * ( re / (re + alts(i,nfl,nll)) ) ** 2
!     .                 * ( gsrx(i,nfl,nll)*bdirsx(i,nfl,nll) +
!     .                     gsry(i,nfl,nll)*bdirsy(i,nfl,nll) +
!     .                     gsrz(i,nfl,nll)*bdirsz(i,nfl,nll)  )


!JK     centrifugal force (see notes 2012/01/04)
        fzero = 3.369
        clat = cos(pie*glats(i,nfl,nll)/180.0)
        slat = sin(pie*glats(i,nfl,nll)/180.0)
        cfs(i,nll)   =  -fzero * 
     .                 (clat*xrg(i,nfl,nll) + slat*xthg(i,nfl,nll))
     .                 * (re + alts(i,nfl,nll)) * clat / re 

!        if (nll.eq.nl/2 .and. nfl.eq.44)
!     .   print *,'new s',i,gs(i,nll)

!       approximation: not good for offset dipole
!       note: sign of gp expicitly accounted for
!       in derivation, i.e., g = -gp phat
!       so gp is positive here

!       gp(i,nfl,nll) = sqrt( gzero**2 - gs(i,nll)**2 ) 

!        if (nll.eq.nl/2 .and. nfl.eq.144)
!     .   print *,'old p',i,gp(i,nfl,nll)

!       should be good for offset dipole
!       note: sign of gp expicitly accounted for
!       in derivation, i.e., g = -gp phat
!       so gp is positive here (JH)

        gp(i,nfl,nll)   = gzero 
     .                 * ( re / (re + alts(i,nfl,nll)) ) ** 2
     .                 * ( gsrx(i,nfl,nll)*vpsnx(i,nfl,nll) +
     .                     gsry(i,nfl,nll)*vpsny(i,nfl,nll) +
     .                     gsrz(i,nfl,nll)*vpsnz(i,nfl,nll)  )


!        if (nll.eq.nl/2 .and. nfl.eq.144)
!     .   print *,'new p',i,gp(i,nfl,nll)

        vnq(i,nfl,nll) = v(i,nfl,nll) *
     .                   ( gsthetax(i,nfl,nll) * bdirsx(i,nfl,nll) +
     .                     gsthetay(i,nfl,nll) * bdirsy(i,nfl,nll) +
     .                     gsthetaz(i,nfl,nll) * bdirsz(i,nfl,nll)   ) +
     .                   u(i,nfl,nll) *
     .                   ( gsphix(i,nfl,nll) * bdirsx(i,nfl,nll) +
     .                     gsphiy(i,nfl,nll) * bdirsy(i,nfl,nll) +
     .                     gsphiz(i,nfl,nll) * bdirsz(i,nfl,nll)   )   +
     .                   w(i,nfl,nll) *
     .                   ( gsrx(i,nfl,nll) * bdirsx(i,nfl,nll) +
     .                     gsry(i,nfl,nll) * bdirsy(i,nfl,nll) +
     .                     gsrz(i,nfl,nll) * bdirsz(i,nfl,nll)   ) 

        vnp(i,nfl,nll) = v(i,nfl,nll) *
     .                   ( gsthetax(i,nfl,nll) * vpsnx(i,nfl,nll) +
     .                     gsthetay(i,nfl,nll) * vpsny(i,nfl,nll) +
     .                     gsthetaz(i,nfl,nll) * vpsnz(i,nfl,nll)   ) +
     .                   u(i,nfl,nll) *
     .                   ( gsphix(i,nfl,nll) * vpsnx(i,nfl,nll) +
     .                     gsphiy(i,nfl,nll) * vpsny(i,nfl,nll) +
     .                     gsphiz(i,nfl,nll) * vpsnz(i,nfl,nll)   )   +
     .                   w(i,nfl,nll) *
     .                   ( gsrx(i,nfl,nll) * vpsnx(i,nfl,nll) +
     .                     gsry(i,nfl,nll) * vpsny(i,nfl,nll) +
     .                     gsrz(i,nfl,nll) * vpsnz(i,nfl,nll)   ) 

        vnphi(i,nfl,nll) = v(i,nfl,nll) *
     .                   ( gsthetax(i,nfl,nll) * vhsnx(i,nfl,nll) +
     .                     gsthetay(i,nfl,nll) * vhsny(i,nfl,nll) +
     .                     gsthetaz(i,nfl,nll) * vhsnz(i,nfl,nll)   ) +
     .                   u(i,nfl,nll) *
     .                   ( gsphix(i,nfl,nll) * vhsnx(i,nfl,nll) +
     .                     gsphiy(i,nfl,nll) * vhsny(i,nfl,nll) +
     .                     gsphiz(i,nfl,nll) * vhsnz(i,nfl,nll)   )   +
     .                   w(i,nfl,nll) *
     .                   ( gsrx(i,nfl,nll) * vhsnx(i,nfl,nll) +
     .                     gsry(i,nfl,nll) * vhsny(i,nfl,nll) +
     .                     gsrz(i,nfl,nll) * vhsnz(i,nfl,nll)   ) 

!        u3(i,nfl,nll) = vnq(i,nfl,nll)
!        u4(i,nfl,nll) = vnp(i,nfl,nll)

        tvn(i,nll)    = vnq(i,nfl,nll) 
          
      enddo

      call update ( tvn,nuin,sumnuj,nuij,nfl,nll )

      do i = 1,nz
        do nni = nion1,nion2
          sumvsi(i,nfl,nll,nni) = 0.
          do nj = nion1,nion2
          sumvsi(i,nfl,nll,nni) =   sumvsi(i,nfl,nll,nni) +
     .                              nuij(i,nni,nj)*vsi(i,nfl,nll,nj)
          enddo
        enddo
      enddo

! define new arrays for velocity and density

      do ni = nion1,nion2
        do i = 1,nz
          vsin (i,ni) = vsi(i,nfl,nll,ni)
          vsidn(i,ni) = vsid(i,nfl,nll,ni)
          denin(i,ni) = deni(i,nfl,nll,ni)
        enddo
      enddo

! define sound velocity used in vsisolv

      do ni = nion1,nion2
        do i = 1,nz
          cfac     = 1.6667 * 8.6174e-5 * te(i,nfl,nll) / ami(ni)
          cs(i,ni) = 9.79e5 * sqrt(cfac)
        enddo
      enddo

! update variables

      do ni = nion1,nion2

        call vsisolv ( vsin(1,ni),vsidn(1,ni),vsi_old(1,ni),
     .                 sumnuj(1,ni),nfl,nll,cs(1,ni) )

! compensating filter

       call smoothz ( vsin(1,ni), 1 )

! put stuff back into velocity array

        do i = 1,nz
          vsi(i,nfl,nll,ni)  = vsin(i,ni)
          vsid(i,nfl,nll,ni) = vsidn(i,ni)
!          if ( alts(i,nfl,nll) .gt. 6.*re .and. ni .eq. pthp ) then
!            vsi(i,nfl,nll,ni)  = 0.
!            vsid(i,nfl,nll,ni) = 0.
!          endif
        enddo

        call densolv2 ( ni,denin(1,ni),
     .       prod(1,ni),loss(1,ni),deni_old(1,ni),nfl,nll )

! put stuff back into density array

        do i = 1,nz
          deni(i,nfl,nll,ni) = denin(i,ni)
        enddo

! put floor on density

        do i = 1,nz
          deni(i,nfl,nll,ni) = amax1 ( deni(i,nfl,nll,ni), denmin )
! below commented out (JK)
          if ( alts(i,nfl,nll) .gt. pcrit*re .and. ni .eq. pthp ) 
     .         deni(i,nfl,nll,ni) = amax1 ( deni(i,nfl,nll,ni), .1 )
          if ( alts(i,nfl,nll) .gt. pcrit*re .and. ni .eq. pthep ) 
     .         deni(i,nfl,nll,ni) = amax1 ( deni(i,nfl,nll,ni), .01 )

        enddo


      enddo

!JK   Basic plasmapause:  for L>__ and nz/2-_ < i < nz/2+_, reduce density
!     L = 4 is height 3*re = 19,110 km
!     L = 3 is height 2*re = 12,740 km
!      altloss = 4.0*re
!      framp = 0.99
!      nz_1 = nz/2 - 3
!      nz_2 = nz/2 + 3
!      do i=nz_1,nz_2
!        if (alts(i,nfl,nll) .gt. altloss) then
!          do ni=1,4
!            if (ni .eq. 1) ip = pthp
!            if (ni .eq. 2) ip = pthep
!            if (ni .eq. 3) ip = ptop
!            if (ni .eq. 4) ip = ptnp
!            deni(i,nfl,nll,ip)=framp*deni(i,nfl,nll,ip)
!            deni(i,nfl,nll,ip) = amax1 ( deni(i,nfl,nll,ip), denmin ) 
!          enddo
!        endif
!      enddo

! JH kill high altitude density 

!      framp = 0.999
!      nz_1 = nz/2 - 7
!      nz_2 = nz/2 + 7

!     epsden is decay rate of density (i.e., loss rate)

!!      epsden = 1.e-2
!!      framp  = (1. - epsden)
!!      pcrit = 6.

!      pcrit = 6.
!      nz_1 = nz/2 - 20
!      nz_2 = nz/2 + 20
!       kk = (taskid - 1) * (nl - 2) + (nll - 1)
!       if (alts(nz/2,nfl,nll)/re .ge. p_crit(kk) ) then

!!        do i=1,nz
!!          if (alts(i,nfl,nll)/re .ge. pcrit ) then
!!            do ni=1,4
!!              if (ni .eq. 1) ip  = pthp
!!              if (ni .eq. 2) ip  = pthep
!!              if (ni .eq. 3) ip  = ptop
!!              if (ni .eq. 4) ip  = ptnp
!!              deni(i,nfl,nll,ip) = framp*deni(i,nfl,nll,ip)
!!              deni(i,nfl,nll,ip) = amax1 ( deni(i,nfl,nll,ip), 10. )
!!            enddo
!!          endif
!!        enddo

! define new arrays for temperature

      do ni = nion1,nion2
        do i = 1,nz
          tin(i,ni)  = ti(i,nfl,nll,ni)
        enddo
      enddo

      do i = 1,nz
        ten(i)  = te(i,nfl,nll)
      enddo

! temperatures (with floors and warnings)

      tn0 = 200. ! floor on temperature

      call etemp  (ten,te_old,phprodr,nfl,nll)
      do i = 1,nz
        te(i,nfl,nll)  = amax1(tn(i,nfl,nll),ten(i))
!        te(i,nfl,nll)  = amax1(tn0,ten(i))
        te(i,nfl,nll)  = amin1(te(i,nfl,nll),1.e4)
        if ( te(i,nfl,nll) .lt. 0 ) then
          print *,' T(e) negative: i,nfl,nll taskid',i,nfl,nll,taskid
          stop
        endif
      enddo



      call htemp  (tin(1,pthp) ,ti_old(1,pthp) ,tvn,nuin,nfl,nll)
      do i = 1,nz
        ti(i,nfl,nll,pthp)  = amax1(tn(i,nfl,nll),tin(i,pthp))
!        ti(i,nfl,nll,pthp)  = amax1(tn0,tin(i,pthp))
        ti(i,nfl,nll,pthp)  = amin1(ti(i,nfl,nll,pthp),1.e4)
        if ( ti(i,nfl,nll,pthp) .lt. 0 ) then
          print *,' T(H) negative: i,nfl,nll',i,nfl,nll 
          stop
        endif
      enddo

      call hetemp (tin(1,pthep),ti_old(1,pthep),tvn,nuin,nfl,nll)
      do i = 1,nz
        ti(i,nfl,nll,pthep)  = amax1(tn(i,nfl,nll),tin(i,pthep))
!        ti(i,nfl,nll,pthep)  = amax1(tn0,tin(i,pthep))
        ti(i,nfl,nll,pthep)  = amin1(ti(i,nfl,nll,pthep),1.e4)
        if ( ti(i,nfl,nll,pthep) .lt. 0 ) then
          print *,' T(He) negative: i,nfl,nll',i,nfl,nll 
          stop
        endif
      enddo

      call otemp  (tin(1,ptop) ,ti_old(1,ptop) ,tvn,nuin,nfl,nll)
      do i = 1,nz
        ti(i,nfl,nll,ptop)  = amax1(tn(i,nfl,nll),tin(i,ptop))
!        ti(i,nfl,nll,ptop)  = amax1(tn0,tin(i,ptop))
        ti(i,nfl,nll,ptop)  = amin1(ti(i,nfl,nll,ptop),1.e4)
        if ( ti(i,nfl,nll,ptop) .lt. 0 ) then
          print *,' T(O) negative: i,nfl,nll',i,nfl,nll 
          stop
        endif
      enddo

      do i = 1,nz
        ti(i,nfl,nll,ptnp )    = ti(i,nfl,nll,ptop)
        ti(i,nfl,nll,ptn2p)    = ti(i,nfl,nll,ptop)
        ti(i,nfl,nll,ptnop)    = ti(i,nfl,nll,ptop)
        ti(i,nfl,nll,pto2p)    = ti(i,nfl,nll,ptop)
      enddo

      return
      end

