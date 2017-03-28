; read output from test
filename='test-spectrum.out'
read_data

; set margin on horizontal sides
plot_spacex = 10

; set plotting environment
!p.charsize=2
!p.color=0  
!p.background=255
!p.thick = 4
bottomline=0

plottitle = ' '
!x.title = 'wavelength [A] ' 
!y.title = 'intensity [erg cm!U-2 !Nsr!U-1 !Ns!U-1 !NA!U-1]'

; plot test output
func='flux'

window,3
cut=grid(30000:*)
plot_data

; analytic result
; calculate references                                                              
kB = 1.38e-23
T  = 10.^6.05
mpc2 = 1.6726219e-27*(3e8)^2

cut = grid(30000:*)
x0 = x[30000:*]

; reference lines 949.354, 949.745, 950.157

m0 = 949.354
s0 = m0 * sqrt(kB * T / mpc2/28.0850000000000009)
p0 =  10.^(-25.1540009173210670)*1e16/(4*3.1415926)/(sqrt(2.*3.1415926)*s0)

y0 = exp(-(x0-m0)^2/(2*s0^2))*p0

m1 = 949.745
s1 = m1 * sqrt(kB * T / mpc2/1.0080000000000000)
p1 =  10.^(-27.0884370735838544)*1e16/(4*3.1415926)/(sqrt(2.*3.1415926)*s1)

y1 = exp(-(x0-m1)^2/(2*s1^2))*p1

m4 = 950.157
s4 = m4 * sqrt(kB * T / mpc2/28.0850000000000009 )
p4 =  10.^(-25.3722235261173061)*1e16/(4*3.1415926)/(sqrt(2.*3.1415926)*s4)

y4 = exp(-(x0-m4)^2/(2*s4^2))*p4

; reference lines 972.538, 972.539, 973.330 

m2 = 972.538
s2 = m2 * sqrt(kB * T / mpc2/1.0080000000000000)
p2 = 10.^(-26.4422997950043985)*1e16/(4*3.1415926)/(sqrt(2.*3.1415926)*s2)

y2 = exp(-(x0-m2)^2/(2*s2^2))*p2

m3 = 972.539
s3 = m3 * sqrt(kB * T / mpc2/1.0080000000000000 )
p3 = 10.^(-26.7442997950043981)*1e16/(4*3.1415926)/(sqrt(2.*3.1415926)*s3)

y3 = exp(-(x0-m3)^2/(2*s3^2))*p3

m5 = 973.330
s5 = m5 * sqrt(kB * T / mpc2/20.1799999999999997)
p5 = 10.^(-26.6680009173210664)*1e16/(4*3.1415926)/(sqrt(2.*3.1415926)*s5)

y5 = exp(-(x0-m5)^2/(2*s5^2))*p5

; reference lines 962.425, 967.540, 966.599 too faint

; plot reference lines                                                              
y = y0 + y1 + y2 + y3 + y4 + y5

oplot,x0,y,color=250,linestyle=3, thick = 4

