import numpy as np
from PIL import Image
im = Image.open("europa2_out.jpg") #Can be many different formats.
pix = im.load()
print im.size #Get the width and hight of the image for iterating over
fo = open("jpgmapdata.dat","w+")
foo = open("jpgmapdata_Xyz.dat","w+")
R=1560E3 #m
#fo.write("longitude latitude     R     G     B"+"\n")

Step=5

fo.write('TITLE = "JPG data "'+'\n')
fo.write('VARIABLES = "Longtitude (degree) "'+'\n')
fo.write('"Latitude (degree) "'+'\n')
fo.write('"Altitude (km) "'+'\n')
fo.write('"R "'+'\n')
fo.write('"G "'+'\n')
fo.write('"B "'+'\n')
fo.write('ZONE T="map with cloud"'+'\n')
fo.write('STRANDID=0, SOLUTIONTIME=0'+'\n')
fo.write('I=%04d'%im.size[0]+', '+'J=%03d'%im.size[1]+', '+'K=1, ZONETYPE=Ordered'+'\n')
fo.write('DATAPACKING=POINT'+'\n')
fo.write('DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE )'+'\n')

for j in range(im.size[1]):
    for i in range(im.size[0]):
        lon=-180.+360./im.size[0]*i
        lat=90.-180./im.size[1]*j
        #print lon, lat
        fo.write("%6.2f"%lon+'  '+"%6.2f"%lat+'  3396.00'"%7d%7d%7d\n" % pix[i,j])
        #fo.write("%10d%10d%10d\n" % pix[i,j])
        #print pix[i,j] #Get the RGBA Value of the a pixel of an image
        #pix[x,y] = value # Set the RGBA Value of the image (tuple)

foo.write('TITLE = "JPG data "'+'\n')
foo.write('VARIABLES = "X"'+'\n')
foo.write('"Y"'+'\n')
foo.write('"Z"'+'\n')
foo.write('"R "'+'\n')
foo.write('"G "'+'\n')
foo.write('"B "'+'\n')
foo.write('ZONE T="map with cloud"'+'\n')
foo.write('STRANDID=0, SOLUTIONTIME=0'+'\n')
#foo.write('I=%04d'%im.size[0]+', '+'J=%03d'%im.size[1]+', '+'K=1, ZONETYPE=Ordered'+'\n')

iNewPoints=1+int((im.size[0]-1)/Step)
jNewPoints=1+int((im.size[1]-1)/Step)

foo.write('N=%04d'%(iNewPoints*jNewPoints)+', E=%04d'%((iNewPoints-1)*(jNewPoints-1))+', ZONETYPE=FEQUADRILATERAL\n')
foo.write('DATAPACKING=POINT'+'\n')
foo.write('DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE )'+'\n')

for j in range(jNewPoints):
    for i in range(iNewPoints):
        phi=-180.+360./(iNewPoints-1)*i
        phi=phi*np.pi/180.
        theta=180./(jNewPoints-1)*j
        theta=theta*np.pi/180.

        #print phi, theta
        x = R*np.sin(theta)*np.cos(phi)
        y = R*np.sin(theta)*np.sin(phi)
        z = R*np.cos(theta)

        foo.write("%f"%x+'  '+"%f"%y+'  '+"%f"%z+\
        #'  '+"%6.2f"%theta+'  '+"%6.2f"%phi+'  '+\
                "%7d%7d%7d\n" % pix[Step*i,Step*j])

#print the connectivity list
for iZenith in range(jNewPoints-1):
    for iAzimuthal in range(iNewPoints-1):
      nd0=1+iAzimuthal+iZenith*iNewPoints  
      nd1=2+iAzimuthal+iZenith*iNewPoints  
      nd2=2+iAzimuthal+(iZenith+1)*iNewPoints  
      nd3=1+iAzimuthal+(iZenith+1)*iNewPoints  
                  
      foo.write('%04d'%nd0+' '+'%04d'%nd1+' '+'%04d'%nd2+' '+'%04d'%nd3+'\n')



print theta, np.cos(theta)
