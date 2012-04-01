#!/usr/bin/env python

'''
Open a GITM 3D file adn create a plot similar to the example given by Aaron.
Note that as pybats.gitm is more developed, a plot like this should be made
using syntax like,
>>>a=gitm.GitmBin('filename')
>>>a.add_alt_slice(0, 'Rho', add_cbar=True)
That's how most pybats stuff works right now.
'''

def quickplot(infile, outfile, min, max):
    '''
    Takes file name, infile, and creates a PNG plot thing.  Yeah.
    '''

    # Import shit.  I needed a lot of shit this time.  
    import numpy as np
    from spacepy.pybats import gitm
    import matplotlib.pyplot as plt
    from matplotlib.cm import get_cmap
    from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

    # Open file.
    a=gitm.GitmBin(infile)
    base = gitm.GitmBin('3DALL_t020321_120100.bin')

    # Make contour of rho at lowest altitude (index 0).
    # Convert lat lon from rad to degrees.
    p=180.0/np.pi
    
    f=plt.figure(figsize=(6,7)) #make a fig.
    ax=f.add_subplot(111) #make an ax.

    # Create the contour for an altitude slice and call it 'cnt'.
    # The '61' is the number of contours; you could use a vector of values to set
    # levels manually if you wish.  get_cmap accepts any of the color map names
    # from the colormap demo pic from the Matplotlib gallery; adding '_r' 
    # reverses the colormap.
    t = a['Temperature'][:,:,:]
    sz = t.shape
    c = sz[0]/2.0
    diff = 100*(a['Rho'][c,:,:] - base['Rho'][c,:,:])/base['Rho'][c,:,:]

    diff[diff>max] = max
    diff[diff<min] = min

#    diff = 100*(a['Temperature'][c,:,:] - base['Temperature'][c,:,:])/base['Temperature'][c,:,:]
#    lev = (np.arange(60)-30.0)

#    diff = np.log10(a['Rho'][c,:,:])
#    diff = a['e-'][c,:,:]

    u = a['V!Dn!N (north)'][c,:,:]
    v = a['V!Dn!N (up)'][c,:,:]

    x = a['Latitude'][c,:,:]*p * 111.0
    y = a['Altitude'][c,:,:]/1000.0

    if (min > -1.0e31 and max < 1.0e32):
        lev = (np.arange(61))/60.0 * (max-min) + min
        cnt=ax.contourf(x,y,
                        diff, 61, cmap=get_cmap('RdBu_r'),
                        levels=lev)
    else:
        cnt=ax.contourf(x,y,
                        diff, 61, cmap=get_cmap('RdBu_r')) #,
#        cnt=ax.contourf(a['Latitude'][c,:,:]*p,
#                        a['Altitude'][c,:,:]/1000.0,
#                        diff, 61, cmap=get_cmap('RdBu_r')) #,
#                    levels=lev)
#                    diff, 61, cmap=get_cmap('BrBG'),
#                    np.log10(a['Rho'][5,:,:]), 61, cmap=get_cmap('Spectral_r'))

    # Configure axis.
    ax.set_xlabel('Latitude (in km)')
    ax.set_ylabel('Altitude (km)')
    ax.set_title(r'Pecent Diff of Rho at longitude=%5.2f$^{\circ}$' % (a['Longitude'][c,0,0]*p))
    f.suptitle('File=%s'%(a.attrs['file']))
    
    q = plt.quiver(x,y,u,v,angles='xy',scale=5000,color='b')

    # Add a colorbar and set the tick format to exponential notation.
    cb=plt.colorbar(cnt)
    cb.formatter=FormatStrFormatter('%7.2E')
    cb.update_ticks()

    # Add the quivers.
#ax.quiver(a['Longitude'][:,:,0]*p, p*a['Latitude'][:,:,0],
#          a['V!Dn!N (east)'][:,:,0],a['V!Dn!N (north)'][:,:,0])

#f.savefig('filename.png')
    if outfile[-4]!='.png':
        outfile+='.png'
    f.savefig(outfile)

# Draw to screen.
    if plt.isinteractive():
        plt.draw() #In interactive mode, you just "draw".
    #else:
        # W/o interactive mode, "show" stops the user from typing more 
        # at the terminal until plots are drawn.
        #plt.show()
    plt.close(f)

from glob import glob
import sys
import re

min = -1.0e32
max =  1.0e32

outfile = 'image'

for arg in sys.argv:
    matchObj = re.match(r'-min=(.*)',arg)
    if matchObj:
        min=float(matchObj.group(1))
        print 'min set to ',repr(min)
    matchObj = re.match(r'-max=(.*)',arg)
    if matchObj:
        max=float(matchObj.group(1))
        print 'max set to ',repr(max)
    matchObj = re.match(r'-out=(.*)',arg)
    if matchObj:
        outfile=matchObj.group(1)
        print 'output file set to ',outfile

i = 0
for f in glob('3DALL*.bin'):
    s = repr(i).zfill(4)
    o = outfile+'.{}'.format(s)
    print f
    quickplot(f, o, min, max)
    i=i+1
