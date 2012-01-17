#!/usr/bin/env python

'''
Open a GITM 3D file adn create a plot similar to the example given by Aaron.
Note that as pybats.gitm is more developed, a plot like this should be made
using syntax like,
>>>a=gitm.GitmBin('filename')
>>>a.add_alt_slice(0, 'Rho', add_cbar=True)
That's how most pybats stuff works right now.
'''

def quickplot(infile, outfile):
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
    base = gitm.GitmBin('3DALL_t020321_120010.bin')

    # Make contour of rho at lowest altitude (index 0).
    # Convert lat lon from rad to degrees.
    p=180.0/np.pi
    
    f=plt.figure(figsize=(5,7)) #make a fig.
    ax=f.add_subplot(111) #make an ax.

    # Create the contour for an altitude slice and call it 'cnt' (no jokes, please.)
    # The '61' is the number of contours; you could use a vector of values to set
    # levels manually if you wish.  get_cmap accepts any of the color map names
    # from the colormap demo pic from the Matplotlib gallery; adding '_r' 
    # reverses the colormap.
    diff = (a['Rho'][5,:,:] - base['Rho'][5,:,:])/base['Rho'][5,:,:]

    cnt=ax.contourf(a['Latitude'][5,:,:]*p,
                    a['Altitude'][5,:,:],
                    diff, 61, cmap=get_cmap('Spectral_r'))
#                    np.log10(a['Rho'][5,:,:]), 61, cmap=get_cmap('Spectral_r'))

    # Configure axis.
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(r'$\rho$ at Altitude=%5.2f$km$' % (a['Altitude'][0,0,0]/1000.0))
    f.suptitle('File=%s'%(a.attrs['file']))
    
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
i = 0
for f in glob('3DALL*.bin'):
    s = repr(i).zfill(4)
    o = 'image.{}'.format(s)
    quickplot(f, o)
    i=i+1
