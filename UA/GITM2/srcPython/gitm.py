#!/usr/bin/env python
#-----------------------------------------------------------------------------
# $Id$
#
# GITM.py, Dan Welling, UMich
#
# Comments: Defines a class for GITM binary output files and modifies the input
#           to improve the data analysis and plotting experience
#
# Contains: class GitmBin    - The class for the GITM binary, which will read
#                              a single GITM output binary
#           def calc_magvel  - Reads a single GITM ion output binary and
#                              uses data from the standard GITM output (3DAll)
#                              to compute the ion characteristics in magnetic
#                              coordinates
#           def calc_deg     - Computes and appends latitude and longitude
#                              in degrees from values in radians
#           def calc_lt      - Computes and appends local time in hours from
#                              the universal time for the file and longitude
#           def append_units - Appends unit, descriptive name, and scale
#                              attributes to each known data type
#           def append_data  - Appends a list of data variables to a GitmBin
#                              data structure, where only a limited number of
#                              data variables from that file have been read
#                              in before
#           def calc_2dion   - Calculate the 2D ionospheric parameters (VTEC,
#                              hmF2, NmF2)
#
# Updates: Angeline Burrell (AGB) - 1/7/13: Added calc_lt, append_units, and
#                                           calc_magvel
#------------------------------------------------------------------------------

'''
PyBats submodule for handling input/output for the Global 
Ionosphere-Thermosphere Model (GITM), one of the choices for the UA module
in the SWMF.
'''

# Global imports:
import numpy as np
import datetime as dt
from spacepy.pybats import PbData
from spacepy.datamodel import dmarray

class GitmBin(PbData):
    '''
    Object to open, manipulate and visualize 1-3 dimensional GITM output
    stored in binary format.  Object inherits from spacepy.pybats.PbData; see
    that documentation for general information on how these objects work.

    GITM index ordering is [lon, lat, altitude]; data arrays read from file
    will always be of the same shape and size.

    kwargs may be specified for:
    magfile = 3DION or 3DMAG file, allows computation of velocities in magnetic
              coordinates
    varlist = list of variable keys.  Will limit the variables read in to those
              listed.  Time and position will always be read in.  If the list
              is empty, all variables will be read in.
    '''

    def __init__(self, filename, *args, **kwargs):
        # Remove any known kwargs, as we don't want them to be included in
        # the GITM data keys!

        if not kwargs.has_key('varlist'):
            varlist = list()
        else:
            varlist = kwargs.pop('varlist')

        if not kwargs.has_key('magfile'):
            magfile = None
        else:
            magfile = kwargs.pop('magfile')

        # Initialize the GITM data structure
        super(GitmBin, self).__init__(*args, **kwargs) # Init as PbData.
        self.attrs['file']=filename

        # Load the GITM data
        self._read(varlist)
        self.calc_deg()
        self.calc_lt()
        self.append_units()

        if magfile:
            self.attrs['magfile']=magfile
            self.calc_magvel()

    def __repr__(self):
        return 'GITM binary output file %s' % (self.attrs['file'])

    def _read(self, varlist, newfile=True):
        '''
        Read binary file.
        '''

        from re import sub
        from struct import unpack
        
        # Read data and header info
        f=open(self.attrs['file'], 'rb')

        # Using the first FORTRAN header, determine endian.
        # Default is little.
        self.attrs['endian']='little'
        endChar='>'
        rawRecLen=f.read(4)
        recLen=(unpack(endChar+'l',rawRecLen))[0]
        if (recLen>10000)or(recLen<0):
            # Ridiculous record length implies wrong endian.
            self.attrs['endian']='big'
            endChar='<'
            recLen=(unpack(endChar+'l',rawRecLen))[0]

        # Read version; read fortran footer+header.
        self.attrs['version']=unpack(endChar+'d',f.read(recLen))[0]
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Read grid size information.
        (self.attrs['nLon'],self.attrs['nLat'],self.attrs['nAlt'])=\
            unpack(endChar+'lll',f.read(recLen))
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Read number of variables.
        self.attrs['nVars']=unpack(endChar+'l',f.read(recLen))[0]
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Collect variable names.
        var=[]
        for i in range(self.attrs['nVars']):
            var.append(unpack(endChar+'%is'%(recLen),f.read(recLen))[0])
            (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Extract time. 
        (yy,mm,dd,hh,mn,ss,ms)=unpack(endChar+'lllllll',f.read(recLen))
        self['time']=dt.datetime(yy,mm,dd,hh,mn,ss,ms/1000)
        (oldLen)=unpack(endChar+'l',f.read(4))


        # Read the rest of the data.
        nTotal=self.attrs['nLon']*self.attrs['nLat']*self.attrs['nAlt']
        for val in var:
            # Trim variable names.
            v=sub('\[|\]', '', val).strip()
            s=unpack(endChar+'l',f.read(4))[0]
            # Test to see if this variable is desired
            gvar=True
            if len(varlist) > 0:
                try:
                    varlist.index(v)
                except ValueError:
                    if((v.find('Altitude') < 0 and v.find('Longitude') < 0
                        and v.find('Latitude') < 0) or not newfile):
                        gvar=False
            # Unpack the data and save, if desired
            temp=unpack(endChar+'%id'%(nTotal),f.read(s))
            if gvar:
                self[v]=dmarray(np.array(temp))
                # Reshape arrays, note that ordering in file is Fortran-like.
                self[v]=self[v].reshape( 
                    (self.attrs['nLon'],self.attrs['nLat'],self.attrs['nAlt']),
                    order='fortran')
                
            f.read(4)


    def calc_deg(self):
        '''
        Gitm defaults to radians for lat and lon, which is sometimes difficult
        to use.  This routine leaves the existing latitude and longitude
        arrays intact and creates *dLat* and *dLon*, which contain the lat and
        lon in degrees.
        '''
        from numpy import pi
        import string

        self['dLat'] = dmarray(self['Latitude']*180.0/pi, 
                               attrs={'units':'degrees', 'scale':'linear',
                                      'name':'Latitude'})
        self['dLon'] = dmarray(self['Longitude']*180.0/pi, 
                               attrs={'units':'degrees', 'scale':'linear',
                                      'name':'Longitude'})

        # Do not correct for over/under limits because these are duplicates
        # that allow contour plots (and more) to display correctly
        #
        #for i in range(self.attrs['nLon']):
        #    for j in range(self.attrs['nLat']):
        #        if self['dLon'][i][j][0] < 0.0:
        #            self['dLon'][i][j] += 360.0
        #        elif self['dLon'][i][j][0] >= 360.0:
        #            self['dLon'][i][j] -= 360.0

    def calc_lt(self):
        '''
        Gitm defaults to universal time.  Compute local time from date and
        longitude.
        '''

        from numpy import pi
        import math

        ut = (self['time'].hour * 3600 + self['time'].minute * 60
              + self['time'].second + self['time'].microsecond * 1e-6) / 3600.0
        self['LT'] = dmarray(ut + self['Longitude']*12.0/pi,
                             attrs={'units':'hours', 'scale':'linear',
                                    'name':'Local Time'})

        # Because datetime won't do lists or arrays
        if dmarray.max(self['LT']) >= 24.0:
            for i in range(self.attrs['nLon']):
                # All local times will be the same for each lat/alt
                # in the same longitude index
                ltmax = dmarray.max(self['LT'][i])
                if ltmax >= 24.0:
                    self['LT'][i] -= 24.0 * math.floor(ltmax / 24.0)

    def calc_magvel(self):
        '''
        Gitm 3DIon files contain the magnetic coordinates that allow the
        field-aligned and field-perpendicular velocities to be computed.
        The 3DIon file must be from the same run as the 3DAll file so that
        the geographic grid is the same.  If a 3D Ion file was not produced
        in the original run, don't fret!  You can get one by using the same
        UAM.in file.  Unless the magnetic field is varying secularly, you
        can get away with producing just one 3DIon file.  It is better to have
        a matching 3D Ion file for every 3D All file, however.
        '''
        import math
        import string
        import sys

        if not self.attrs.has_key('magfile'):
            print "No 3D MAG/ION file associated with this GITM Binary"
        elif(self.attrs['magfile'].find("ION") <= 0 and
             self.attrs['magfile'].find("MAG") <= 0):
            print "No 3D MAG/ION file associated with this GITM Binary"
        else:
            ion = GitmBin(self.attrs['magfile'])

            # Compute the field-aligned unit vector in East, North,
            # and Vertical coordinates

            bhat_e = ion['B.F. East'] / ion['B.F. Magnitude']
            bhat_n = ion['B.F. North'] / ion['B.F. Magnitude']
            bhat_v = ion['B.F. Vertical'] / ion['B.F. Magnitude']

            # Compute the zonal unit vector in East, North, Vertical coord

            mag = np.sqrt(np.square(ion['B.F. East'])
                          + np.square(ion['B.F. North']))

            zhat_e = -ion['B.F. North'] / mag
            zhat_n = ion['B.F. East'] / mag
            # zhat_v is identically zero

            # Compute the meridional unit vector in East, North, Vertical coord

            mhat_e = (-ion['B.F. East']*ion['B.F. Vertical']
                       / (mag * ion['B.F. Magnitude']))
            mhat_n = (-ion['B.F. North']*ion['B.F. Vertical']
                       / (mag * ion['B.F. Magnitude']))
            mhat_v = mag / ion['B.F. Magnitude']

            # Compute the magnetic coordinate velocities for each overlapping
            # latitude, longitude, and altitude.  Also include the mag coord.

            eps  = 1.0e-3
            skey = self.keys()
            vkey = dict()
            nkey = dict()

            for item in skey:
                if string.find(item, "V!") >= 0:
                    sp = string.split(item)
                    if not sp[0] in vkey.keys():
                        east  = string.join([sp[0], "(east)"], " ")
                        north = string.join([sp[0], "(north)"], " ")
                        up    = string.join([sp[0], "(up)"], " ")
                        par   = string.join([sp[0], "(par)"], " ")
                        zon   = string.join([sp[0], "(zon)"], " ")
                        mer   = string.join([sp[0], "(mer)"], " ")

                        vp=bhat_e*self[east]+bhat_n*self[north]+bhat_v*self[up]
                        vz=zhat_e*self[east]+zhat_n*self[north]
                        vm=mhat_e*self[east]+mhat_n*self[north]+mhat_v*self[up]
                        self[par] = dmarray.copy(vp)
                        self[par].attrs = {'units':'m s^{-1}{\mathdefault{,\,positive\,mag\,north}}',
                                           'scale':'linear',
                                           'name':'v$_\parallel$'}
                        self[zon] = dmarray.copy(vz)
                        self[zon].attrs = {'units':'m s^{-1}{\mathdefault{,\,positive\,east}}',
                                           'scale':'linear', 'name':'v$_{zon}$'}
                        self[mer] = dmarray.copy(vm)
                        self[mer].attrs = {'units':'m s^{-1}{\mathdefault{,\,positive\,up}}',
                                           'scale':'linear', 'name':'v$_{mer}$'}
        
            self['B.F. East']          = dmarray.copy(ion['B.F. East'])
            self['B.F. North']         = dmarray.copy(ion['B.F. North'])
            self['B.F. Vertical']      = dmarray.copy(ion['B.F. Vertical'])
            self['B.F. Magnitude']     = dmarray.copy(ion['B.F. Magnitude'])
            self['Magnetic Latitude']  = dmarray.copy(ion['Magnetic Latitude'])
            self['Magnetic Longitude'] = dmarray.copy(ion['Magnetic Longitude'])

    def append_units(self):
        '''
        Append units, descriptive names, and plot scaling (e.g. linear,
        exponential) to the attributes of known data types
        '''
        unit_dict = {"Altitude":"m", "Ar Mixing Ratio":"", "Ar":"m^{-3}",
                     "CH4 Mixing Ratio":"", "Conduction":"W m$^{-1}$ K$^{-1}$",
                     "EuvHeating":"", "H":"m^{-3}", "H!U+!N":"m^{-3}",
                     "H2 Mixing Ratio":"", "HCN Mixing Ratio":"", "He":"m^{-3}",
                     "He!U+!N":"m^{-3}", "Heating Efficiency":"",
                     "Heat Balance Total":"", "Latitude":"radians",
                     "Longitude":"radians", "N!D2!N":"m^{-3}",
                     "N!D2!U+!N":"m^{-3}", "N!U+!N":"m^{-3}",
                     "N(!U2!ND)":"m^{-3}", "N(!U2!NP)":"m^{-3}",
                     "N(!U4!NS)":"m^{-3}", "N2 Mixing Ratio":"", "NO":"m^{-3}",
                     "NO!U+!N":"m^{-3}", "O!D2!N":"m^{-3}",
                     "O(!U1!ND)":"m^{-3}", "O!D2!U+!N":"m^{-3}",
                     "O(!U2!ND)!":"m^{-3}", "O(!U2!ND)!U+!N":"m^{-3}",
                     "O(!U2!NP)!U+!N":"m^{-3}", "O(!U2!NP)!U+!N":"m^{-3}",
                     "O(!U3!NP)":"m^{-3}", "O_4SP_!U+!N":"m^{-3}",
                     "RadCooling":"", "Rho":"m^{-3}", "Temperature":"K",
                     "V!Di!N (east)":"m s^{-1}", "V!Di!N (north)":"m s^{-1}",
                     "V!Di!N (up)":"m s^{-1}", "V!Dn!N (east)":"m s^{-1}",
                     "V!Dn!N (north)":"m s^{-1}", "V!Dn!N (up)":"m s^{-1}",
                     "V!Dn!N (up,N!D2!N              )":"m s^{-1}",
                     "V!Dn!N (up,N(!U4!NS)           )":"m s^{-1}",
                     "V!Dn!N (up,NO                  )":"m s^{-1}",
                     "V!Dn!N (up,O!D2!N              )":"m s^{-1}",
                     "V!Dn!N (up,O(!U3!NP)           )":"m s^{-1}",
                     "e-":"m^{-3}", "Electron_Average_Energy":"J",
                     "eTemperature":"K", "iTemperature":"K", "LT":"h",
                     "Solar Zenith Angle":"radians", "Vertical TEC":"TECU",
                     "CO!D2!N":"m^{-3}",  "DivJu FL":"", "DivJuAlt":"",
                     "Electron_Energy_Flux":"J m$^{-2}$", "FL Length":"m",
                     "Pedersen FL Conductance":"S m^{-1}", "dLon":"degrees",
                     "Pedersen Conductance":"S m^{-1}", "dLat":"degrees",
                     "Hall FL Conductance":"S m^{-1}", "Potential":"V",
                     "Hall Conductance":"S m^{-1}", "Je2":"A m^{-2}",
                     "Je1":"A m^{-2}", "Magnetic Longitude":"degrees",
                     "E.F. Vertical":"V m^{-1}", "E.F. East":"V m^{-1}",
                     "E.F. North":"V m^{-1}", "E.F. Magnitude":"V m^{-1}",
                     "B.F. Vertical":"nT", "B.F. East":"nT", "B.F. North":"nT",
                     "B.F. Magnitude":"nT", "Magnetic Latitude":"degrees",
                     "Ed1":"", "Ed2":""}

        scale_dict = {"Altitude":"linear", "Ar Mixing Ratio":"linear",
                      "Ar":"exponential", "CH4 Mixing Ratio":"linear",
                      "Conduction":"linear", "EuvHeating":"linear",
                      "H":"exponential", "H!U+!N":"exponential",
                      "H2 Mixing Ratio":"linear", "HCN Mixing Ratio":"linear",
                      "He":"exponential", "He!U+!N":"exponential",
                      "Heating Efficiency":"linear", "DivJuAlt":"linear",
                      "Heat Balance Total":"linear", "Latitude":"linear",
                      "Longitude":"linear", "N!D2!N":"exponential",
                      "N!D2!U+!N":"exponential", "N!U+!N":"exponential",
                      "N(!U2!ND)":"exponential", "N(!U2!NP)":"exponential",
                      "N(!U4!NS)":"exponential", "N2 Mixing Ratio":"linear",
                      "NO":"exponential", "NO!U+!N":"exponential",
                      "O!D2!N":"exponential", "O!D2!U+!N":"exponential",
                      "O(!U2!ND)!":"exponential", "O(!U1!ND)":"exponential",
                      "O(!U2!ND)!U+!N":"exponential", "CO!D2!N":"exponential",
                      "O(!U2!NP)!U+!N":"exponential", "DivJu FL":"",
                      "O(!U2!NP)!U+!N":"exponential", "O(!U3!NP)":"exponential",
                      "O_4SP_!U+!N":"exponential", "RadCooling":"linear",
                      "Rho":"exponential", "Temperature":"linear",
                      "V!Di!N (east)":"linear", "V!Di!N (north)":"linear",
                      "V!Di!N (up)":"linear", "V!Dn!N (east)":"linear",
                      "V!Dn!N (north)":"linear", "V!Dn!N (up)":"linear",
                      "V!Dn!N (up,N!D2!N              )":"linear",
                      "V!Dn!N (up,N(!U4!NS)           )":"linear",
                      "V!Dn!N (up,NO                  )":"linear",
                      "V!Dn!N (up,O!D2!N              )":"linear",
                      "V!Dn!N (up,O(!U3!NP)           )":"linear",
                      "e-":"exponential", "Electron_Average_Energy":"linear",
                      "eTemperature":"linear", "iTemperature":"linear",
                      "Solar Zenith Angle":"linear", "Vertical TEC":"linear",
                      "Electron_Energy_Flux":"exponential",
                      "FL Length":"linear", "Pedersen FL Conductance":"linear",
                      "Hall Conductance":"linear", "Potential":"linear",
                      "Hall FL Conductance":"linear", "dLon":"linear",
                      "Pedersen Conductance":"linear", "Je2":"linear",
                      "Je1":"linear", "Ed1":"linear", "Ed2":"linear",
                      "E.F. Vertical":"linear", "E.F. East":"linear",
                      "E.F. North":"linear", "E.F. Magnitude":"linear",
                      "B.F. Vertical":"linear", "B.F. East":"linear",
                      "B.F. North":"linear", "B.F. Magnitude":"linear",
                      "Magnetic Latitude":"linear", "LT":"linear",
                      "Magnetic Longitude":"linear", "dLat":"linear"}

        name_dict = {"Altitude":"Altitude",
                     "Ar Mixing Ratio":"Argon Mixing Ratio", "Ar":"[Ar]",
                     "CH4 Mixing Ratio":"Methane Mixing Ratio",
                     "Conduction":"Conduction", "EuvHeating":"EUV Heating",
                     "H":"[H]", "H!U+!N":"[H$^+$]",
                     "H2 Mixing Ratio":"H$_2$ Mixing Ratio",
                     "HCN Mixing Ratio":"Hydrogen Cyanide Mixing Ratio",
                     "He":"[He]", "He!U+!N":"[He$^+$]",
                     "Heating Efficiency":"Heating Efficiency",
                     "Heat Balance Total":"Heat Balance Total",
                     "Latitude":"Latitude", "Longitude":"Longitude",
                     "N!D2!N":"[N$_2$]", "N!D2!U+!N":"[N$_2$$^+$]",
                     "N!U+!N":"[N$^+$]", "N(!U2!ND)":"[N($^2$D)]",
                     "N(!U2!NP)":"[N($^2$P)]", "N(!U4!NS)":"[N($^4$S)]",
                     "N2 Mixing Ratio":"N$_2$ Mixing Ratio",
                     "NO":"[NO]", "NO!U+!N":"[NO$^+$]", "O!D2!N":"[O$_2$]",
                     "O(!U1!ND)":"[O($^1$D)]", "O!D2!U+!N":"[O$_2$$^+$]",
                     "O(!U2!ND)!":"[O($^2$D)]", "O(!U2!ND)!U+!N":"[O($^2$D)]",
                     "O(!U2!NP)!U+!N":"[O($^2$P)$^+$]",
                     "O(!U2!NP)!U+!N":"[O($^2$P)]", "O(!U3!NP)":"[O($^3$P)]",
                     "O_4SP_!U+!N":"[O($^4$SP)$^+$]",
                     "RadCooling":"Radiative Cooling", "Rho":"Neutral Density",
                     "Temperature":"T$_n$", "V!Di!N (east)":"v$_{East}$",
                     "V!Di!N (north)":"v$_{North}$", "V!Di!N (up)":"v$_{Up}$",
                     "V!Dn!N (east)":"u$_{East}$",
                     "V!Dn!N (north)":"u$_{North}$", "V!Dn!N (up)":"u$_{Up}$",
                     "V!Dn!N (up,N!D2!N              )":"u$_{Up, N_2}$",
                     "V!Dn!N (up,N(!U4!NS)           )":"u$_{Up, N(^4S)}$",
                     "V!Dn!N (up,NO                  )":"u$_{Up, NO}$",
                     "V!Dn!N (up,O!D2!N              )":"u$_{Up, O_2}$",
                     "V!Dn!N (up,O(!U3!NP)           )":"u$_{Up, O(^3P)}$",
                     "e-":"[e-]",
                     "Electron_Average_Energy":"Electron Average Energy",
                     "eTemperature":"T$_e$", "iTemperature":"T$_i$",
                     "Solar Zenith Angle":"Solar Zenith Angle",
                     "Vertical TEC":"VTEC", "CO!D2!N":"[CO$_2$]",
                     "DivJu FL":"DivJu FL", "DivJuAlt":"DivJuAlt",
                     "Electron_Energy_Flux":"Electron Energy Flux",
                     "FL Length":"Field Line Length",
                     "Pedersen FL Conductance":"$\sigma_P$",
                     "Pedersen Conductance":"$\Sigma_P$",
                     "Hall FL Conductance":"$\sigma_H$",
                     "Potential":"Potential", "Hall Conductance":"$\Sigma_H$",
                     "Je2":"Region 2 Current", "Je1":"Region 1 Current",
                     "Ed1":"Ed1", "Ed2":"Ed2", "LT":"Solar Local Time",
                     "E.F. Vertical":"Vertical Electric Field",
                     "E.F. East":"Eastward Electric Field",
                     "E.F. North":"Northward Electric Field",
                     "E.F. Magnitude":"Electric Field Magnitude",
                     "B.F. Vertical":"Vertical Magnetic Field",
                     "B.F. East":"Eastward Magnetic Field",
                     "B.F. North":"Northward Magnetic Field",
                     "B.F. Magnitude":"Magnetic Field Magnitude",
                     "Magnetic Latitude":"Magnetic Latitude",
                     "Magnetic Longitude":"Magnetic Longitude",
                     "dLat":"Latitude", "dLon":"Longitude"}

        for k in self.keys():
            if type(self[k]) is dmarray:
                nk = k
                # Temporary fix for misspelled key (9/30/13)
                if nk.find("Heaing Efficiency") >= 0:
                    nk = "Heating Efficiency"

                # Different versions of GITM differ in header capitalization
                if not name_dict.has_key(nk):
                    # Try to capitalize or lowercase the key
                    if nk == nk.capitalize():
                        nk = k.lower()
                    else:
                        nk = k.capitalize()

                if not self[k].attrs.has_key('units'):
                    self[k].attrs['units'] = unit_dict[nk]
                if not self[k].attrs.has_key('scale'):
                    self[k].attrs['scale'] = scale_dict[nk]
                if not self[k].attrs.has_key('name'):
                    self[k].attrs['name'] = name_dict[nk]

    def append_data(self, varlist):
        '''
        A routine to append more variables to an existing GITM data structure.
        New variables can only be added from the same file (specified in
        the 'file' attribute).
        '''

        temp = GitmBin(self.attrs['file'], varlist, False)

        for nkey in temp.keys():
            if not self.has_key(nkey):
                self[nkey] = dmarray.copy(temp[nkey])

    def calc_2dion(self):
        '''
        A routine to calculate the 2D ionospheric parameters: VTEC, hmF2, NmF2.
        To perform these calculations, electron density ("e-") must be one of
        the available data types.
        '''

        import scipy.integrate as integ
        from scipy.interpolate import interp1d

        if self.has_key('e-'):
            self['VTEC'] = dmarray(self['e-'] * 1.0e-16,
                                   attrs={'units':'TECU', 'scale':'linear',
                                          'name':'Vertical TEC'})
            self['NmF2'] = dmarray.copy(self['e-'])
            self['NmF2'].attrs['name'] = 'N$_m$F$_2$'
            self['hmF2'] = dmarray(self['Altitude'] / 1000.0,
                                   attrs={'units':'km', 'scale':'linear',
                                          'name':'h$_m$F$_2$'})

            for ilon in range(self.attrs['nLon']):
                for ilat in range(self.attrs['nLat']):
                    # Integrate electron density over altitude
                    vtec = integ.simps(self['VTEC'][ilon,ilat,:],
                                       self['Altitude'][ilon,ilat,:], "avg")
                    self['VTEC'][ilon,ilat,:] = vtec

                    # Interpolate over the electron density altitude profile
                    eprof = interp1d(self['hmF2'][ilon,ilat,:],
                                     self['NmF2'][ilon,ilat,:], kind="cubic")
                    alt = np.linspace(200.0,max(self['hmF2'][ilon,ilat,:]),1000)
                    edens = list(eprof(alt))
                    self['NmF2'][ilon,ilat,:] = max(edens)
                    self['hmF2'][ilon,ilat,:] = alt[edens.index(max(edens))]
