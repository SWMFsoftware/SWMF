import numpy as np
from PIL import Image
im = Image.open("map1.jpg") #Can be many different formats.
pix = im.load()
print im.size #Get the width and hight of the image for iterating over
fo = open("jpgmapdata.dat","w+")
foo = open("jpgmapdata_Xyz.dat","w+")
R=3396 #km
#fo.write("longitude latitude     R     G     B"+"\n")

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
foo.write('VARIABLES = "X [R] "'+'\n')
foo.write('"Y [R] "'+'\n')
foo.write('"Z [R] "'+'\n')
foo.write('"R "'+'\n')
foo.write('"G "'+'\n')
foo.write('"B "'+'\n')
foo.write('ZONE T="map with cloud"'+'\n')
foo.write('STRANDID=0, SOLUTIONTIME=0'+'\n')
foo.write('I=%04d'%im.size[0]+', '+'J=%03d'%im.size[1]+', '+'K=1, ZONETYPE=Ordered'+'\n')
foo.write('DATAPACKING=POINT'+'\n')
foo.write('DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE )'+'\n')

for j in range(im.size[1]):
    for i in range(im.size[0]):
        phi=-180.+360./im.size[0]*i
        phi=phi*np.pi/180.
        theta=180./im.size[1]*j
        theta=theta*np.pi/180.
        #print phi, theta
        x = np.sin(theta)*np.cos(phi)
        y = np.sin(theta)*np.sin(phi)
        z = np.cos(theta)
        foo.write("%6.2f"%x+'  '+"%6.2f"%y+'  '+"%6.2f"%z+\
        #'  '+"%6.2f"%theta+'  '+"%6.2f"%phi+'  '+\
                "%7d%7d%7d\n" % pix[i,j])
print theta, np.cos(theta)
