; read output from test
filename='spectrum.out'
.r getpict

; set margin on horizontal sides
plot_spacex = 10

; set plotting environment
!p.charsize=1.5
!p.color=0  
!p.background=255

; open plotting device
set_device,'test-spectrum.eps',/eps,/land
!x.margin=[10,10] 
!y.range=[0,8]

; cut a 1D line and find peaks
n = size(grid(*,0,0))
hn = n[1]/2
cut = grid(*,0,0)
x0 = x[*,0,0,0]

; plot test output
func='flux'
.r plotfunc

; calculate references
kB = 1.38e-23
T  = 6e6
mc2 = 1.6726219e-27*(3e8)^2 

; reference line 949.745
m1 = 949.745
p1 = max(w[0:hn,0,0])
s1 = sqrt( (949.745*1e-10)^2.0 * kB * T / mc2 )*1e10 

oplot,x0,exp(-(x0-m1)^2/(2*s1^2))*p1,color=150,linestyle=2

; reference lines 972.538 and 972.539
m2 = 972.5385
p2 = max(w[hn+1:n[1]-1,0,0])
s2 = sqrt( (972.5385*1e-10)^2.0 * kB * T / mc2 )*1e10

oplot,x0,exp(-(x0-m2)^2/(2*s2^2))*p2,color=250,linestyle=3

; close device and idl
close_device
exit
