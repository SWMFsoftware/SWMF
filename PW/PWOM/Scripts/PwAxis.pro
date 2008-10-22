r80=sin(10.*!pi/180)
r70=sin(20.*!pi/180)
r60=sin(30.*!pi/180)
r50=sin(40.*!pi/180)

xTheta=findgen(360)*!pi/180.


oplot,r80*cos(xTheta),r80*sin(xTheta),thick=2
oplot,r70*cos(xTheta),r70*sin(xTheta),thick=2
oplot,r60*cos(xTheta),r60*sin(xTheta),thick=2
oplot,r50*cos(xTheta),r50*sin(xTheta),thick=2

xyouts,(r80+.05)*cos(135.*!pi/180.),(r80+.05)*sin(135.*!pi/180.),'80',charsize=1.

xyouts,(r70+.05)*cos(135.*!pi/180.),(r70+.05)*sin(135.*!pi/180.),'70',charsize=1.

xyouts,(r60+.05)*cos(135.*!pi/180.),(r60+.05)*sin(135.*!pi/180.),'60',charsize=1.

xyouts,(r50+.05)*cos(135.*!pi/180.),(r50+.05)*sin(135.*!pi/180.),'50',charsize=1.

xyouts,(r50+.05)*cos(90.*!pi/180.),(r50+.05)*sin(90.*!pi/180.),'12',charsize=1.5,charthick=2
xyouts,(r50+.05)*cos(0.*!pi/180.),(r50+.05)*sin(0.*!pi/180.),'6',charsize=1.5,charthick=2
xyouts,(r50+.07)*cos(180.*!pi/180.),(r50+.07)*sin(180.*!pi/180.),'18',charsize=1.5,charthick=2
xyouts,(r50+.07)*cos(270.*!pi/180.),(r50+.07)*sin(270.*!pi/180.),'0',charsize=1.5,charthick=2

end
