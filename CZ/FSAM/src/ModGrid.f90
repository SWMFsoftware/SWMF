module ModGrid
  use ModPar,   ONLY: inmax, jnmax, in, jn, kn
  implicit none
  
  integer, parameter :: is = 3, js = 3, ks = 3, ie = in-3, je = jn-3, ke = kn-3
  integer, parameter :: ione = 1, jone = 1, kone = 1
  integer, parameter :: ism1 = is - ione, jsm1 = js - jone, ksm1 = ks - kone
  integer, parameter :: ism2 = is - ione - ione, jsm2 = js - jone - jone, &
       ksm2 = ks - kone - kone
  integer, parameter :: isp1 = is + ione, jsp1 = js + jone, ksp1 = ks + kone
  integer, parameter :: isp2 = is + ione + ione, jsp2 = js + jone + jone, &
       ksp2 = ks + kone + kone
  integer, parameter :: iem1 = ie - ione, jem1 = je - jone, kem1 = ke - kone
  integer, parameter :: iem2 = ie - ione - ione, jem2 = je - jone - jone, &
       kem2 = ke - kone - kone
  integer, parameter :: iep1 = ie + ione, jep1 = je + jone, kep1 = ke + kone
  integer, parameter :: iep2 = ie + ione + ione, jep2 = je + jone + jone, &
       kep2 = ke + kone + kone
  integer, parameter :: iep3 = ie + ione + ione + ione, &
       jep3 = je + jone + jone + jone, kep3 = ke + kone + kone + kone
  
  real    :: dx1a(in), dx2a(jn), dx3a(kn)
  real    :: dx1b(in), dx2b(jn), dx3b(kn)
  real    :: dx1ai(in), dx2ai(jn), dx3ai(kn)
  real    :: dx1bi(in), dx2bi(jn), dx3bi(kn)
  real    :: x1a(in), x2a(jn), x3a(kn)
  real    :: x1b(in), x2b(jn), x3b(kn)

  real    :: g2a(in), dg2ad1(in), g2ai(in)
  real    :: g2b(in), dg2bd1(in), g2bi(in)
  real    :: g31a(in), dg31ad1(in), g31ai(in)
  real    :: g31b(in), dg31bd1(in), g31bi(in)
  real    :: g32a(jn), dg32ad2(jn), g32ai(jn)
  real    :: g32b(jn), dg32bd2(jn), g32bi(jn)
  
  real    :: dvl1a(in), dvl1ai(in), dvl1b(in), dvl1bi(in)
  real    :: dvl2a(jn), dvl2ai(jn), dvl2b(jn), dvl2bi(jn)
  real    :: dvl3a(kn), dvl3ai(kn), dvl3b(kn), dvl3bi(kn)

  real    :: xxa(inmax), xxb(inmax), dxxa(inmax), dxxb(inmax), &
             g2xxa(inmax), g2xxb(inmax), g31xxa(inmax), g31xxb(inmax)

  real    :: yya(jnmax), yyb(jnmax), dyya(jnmax), dyyb(jnmax), &
             g32yya(jnmax),g32yyb(jnmax)

end module ModGrid
