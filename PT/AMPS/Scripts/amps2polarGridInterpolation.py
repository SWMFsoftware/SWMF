#!/usr/bin/env python
'''
call with: python amps2polarGridInterpolation [filename, variable1 variable2 variable3 ...]

filename: path+filename of the amps datafile
variables: name of the variables to interpolate on the polar grid (must be the same variable name as in the amps datafile)

reads amps output file specified by filename and interpolates  selected variables
onto a stretched log grid.

number of grid points in radial direction = N_r
number of grid points in azimuthal direction = N_phi
'''



from __future__ import division
import sys
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def reduceWhiteSpaces(line):
    # reduces number of consecutive white spaces in line to 1
    
    newLine = ' '
    for character in line:
        if character != " ":
            newLine += character
        else:
            if newLine[-1] == ' ':
                pass
            else:
                newLine += character
                
    return newLine[1:]



if len(sys.argv) > 1:
    filename = sys.argv[1]
    for variable in sys.argv[2:-1]:
        variableListUser.append(variable)
    
else:
    filename = 'F21.H2O.dsmc.dat'
    variableListUser = ['n','Vx','Vy','T']
    
nrOfVars = len(variableListUser)

#########################################################
# READ HEADER from amps output file to get variable names
#########################################################
file = open(filename, 'r')

for line in file:
    if "VARIABLES" in line:
        variableList = line.replace(" ","").split('"')[1:-1:2] # remove all white spaces of line, then split into list by '"' and use only every second element, starting with 1
        break
file.close()


variableIndices = []

for var in variableListUser:
    variableIndices.append(variableList.index(var))

for var,varI in zip(variableListUser,variableIndices):
    print var, varI
    
#########################################################
# re open amps output file and LOAD DATA into x, y and z
#########################################################

x = []
y = []
z = [[] for i in range(nrOfVars)]               # holds values of variables to interpolate from

file = open(filename,'r')
for line in file:
    data = reduceWhiteSpaces(line).split(" ")
    try:
        if len(data) > 6:
            x.append(np.float(data[0]))
            y.append(np.float(data[1]))
            j = 0
            for varI in variableIndices:
                z[j].append(np.float(data[varI]))
                j+=1
    except:
        pass
file.close()


x = np.array(x)
y = np.array(y)
for i in range(nrOfVars):
    z[i] = np.array([value+1 for value in z[i]])

print 'data loaded'

triangles = mpl.tri.Triangulation(x,y)
print 'tirangulation done'

Interpolator = [mpl.tri.LinearTriInterpolator(triangles,z[i]) for i in range(nrOfVars)]

N_phi = 181
phiMin = 0
phiMax = np.pi
dphi = (phiMax - phiMin) / (N_phi-1)

N_r = 800
rMin = np.log10(2000)
rMax = np.log10(1e9)
dr = (rMax - rMin) / (N_r-1)

xInterpolate = []
yInterpolate = []
zInterpolate = [[] for i in range(nrOfVars)]

#####################################################################
# build coordinates to interpolate to (polar coordinates)
######################################################################
for phi in np.arange(phiMin, phiMax+dphi, dphi):
    for r in np.arange(rMin,rMax+dr, dr):
        xInterpolate.append(10**r * np.cos(phi))
        yInterpolate.append(10**r * np.sin(phi))
        
for i in range(nrOfVars):
    zInterpolate[i] = Interpolator[i].__call__(xInterpolate,yInterpolate)
print "interpolation done"


###########################################
# write results in file
###########################################
file = open(filename + ".interpol", 'w')

file.write("Data interpolated from %s with %s (%sT%s)\n"% (filename, sys.argv[0],time.strftime('%x'), time.strftime('%X')))
file.write("0 0.0 2 1 %i\n" % (nrOfVars+2))   # (timeSteps Time nDim nParameter nVariables+2) --> +2 because x and y are also written to file
file.write("%i %i\n" %(N_r, N_phi))
file.write("2000\n")                                       # Parameter
file.write("logR phi ")
for var in variableListUser:
    file.write(var + " ")
    
file.write("x y rBody\n")

i=0
for phi in np.arange(phiMin, phiMax+dphi, dphi):
    for r in np.arange(rMin,rMax+dr, dr):
        
        file.write("%f %f " %(r, phi/np.pi*180))
        
        
        for j in range(len(zInterpolate)):
            file.write("%e " %(zInterpolate[j][i]))
            
        file.write("%e %e" %(xInterpolate[i], yInterpolate[i]))
        file.write("\n")
        i += 1
            
file.close()


##################################
# make plot of results
##################################
try:
    xPlotMax = 5e4
    xPlotMin = -xPlotMax
    yPlotMax = xPlotMax
    yPlotMin = 0
    
    rBody = 2000
    tPlt = np.arange(0,np.pi,0.01)
    xPlt = np.cos(tPlt) * rBody
    yPlt = np.sin(tPlt) * rBody
    
    i = 0
    for variable in variableListUser:
        
        if variable == 'n':
            z[i] = np.log10(np.array(z[i]+1))                   # cheat for plotting (number densities only): add 1 to all the datapoints --> log10 does not complain about zeros
            zInterpolate[i] = np.log10(zInterpolate[i]+1)
            
        i +=1
        
    j = 0
    
    for i in range(1,2*len(z)+1,2):
        
        subplotKey1 = 100 * len(z) + 10 *2 + i
        subplotKey2 = subplotKey1+1
        
        
        plt.subplot(subplotKey1)
        plt.tricontourf(xInterpolate,yInterpolate,zInterpolate[j],100)
        plt.fill(xPlt,yPlt,'k')
        plt.xlim((xPlotMin,xPlotMax))
        plt.ylim((yPlotMin,yPlotMax))
        plt.colorbar(label=variableListUser[j])
        
        plt.subplot(subplotKey2)
        plt.tricontourf(x,y,z[j],100)
        plt.fill(xPlt,yPlt,'k')
        plt.xlim((xPlotMin,xPlotMax))
        plt.ylim((yPlotMin,yPlotMax))
        plt.colorbar(label=variableListUser[j])
        j +=1
    plt.show()
except Exception,e:
    print "*******************************"
    print "Error when plotting results:"
    print e
    print "*******************************"
