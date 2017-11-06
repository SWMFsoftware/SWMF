Update, May 2015: This is an updated user-friendly weimer folder. 

SAMI3/Weimer05: 

This is set up so we first run the Weimer05 model to get the potential
on the SAMI3 grid for a multi-day series of runs. 
The resulting file is phi_weimer.inp, which is read by SAMI3. 
Upon restart, SAMI3 also reads nweimer.rst, which tells SAMI3 how many
records in phi_weimer.inp to skip over before using the data.

1. Run SAMI3 for one time step to get grid files blatpu.dat and blonpu.dat
   Copy those files into this directory also copy 'param3_mpi-1.90w_p.inc'
   and 'com3_mpi-1.90w_p.inc' (exact names may vary)
   Needed files: blatpu.dat, blonpu.dat, param3_mpi..., com3_mpi...

2. Modify weimer_grid.f with the correct param3... and com3... names.

3. Generate the weimer grid from the SAMI3 grid using weimer_grid.f
   >ifort -s weimer_grid.f -o grid.x
   ...
   >grid.x

   This creates file weimer_grid.dat

4. Create a solar wind input file for weimer, 'weimer_input.inp'.

   The input file looks like the lines below, it is smoothed from the 
   measured/processed data (i.e., the OMNI SW data) with a 20-minute 
   window, output every 15 minutes, and covers the entire length of the
   planned run series. 
   Decimal times are monotonically increasing and are relative to the 
   specified DOY (e.g. DOY 30, hUT 25.0 = 0100 UT on DOY 31).

   The input file looks like this (the first 4 lines are info text):
    
     Data comes from file(s): OMNI_HRO_1MIN_day30-35.txt      
     Cadence, window (min):    15.00000       20.00000    
     Npts, Year, DOY, <VSW>[km/s], spacecraft-to-1AU time [h]
     hUT[decimal], By[nt], Bz[nt], -Vx[km/s], n[cm^-3], tilt[degrees]
       671    2001      30      375.45419        0.00000
         0.16667       -1.20000        4.02800      456.76321        4.74526      -20.54445
         0.41667       -1.64810        3.80714      457.85718        4.41381      -21.19313
         ...

5. Make sure the input file looks OK by running plot_weimer_input.pro
   >idl_7.0
   IDL>.r plot_weimer_input (creates file fig_sw.ps)

6. Run the makefile and run the code
   Make sure test_w05sc_mod2.f90 reads in "weimer_grid.dat"
   Make sure, in test_w05sc_mod2.f90, nlat is set to nf+1 and 
   nlon is set to nl+1:
     integer,parameter :: nlat=125, nlon=97

   Needed are: Makefile read_data.f90 w05sc.f90 test_w05sc_mod2.f90
               W05scEpot.dat SCHAtable.dat W05scBndy.dat W05scBpot.dat

   > make
   >...
   > w05sc.x

   This produces phi_weimer.inp

7. Check the output.
   Need files: read-basic.pro read-weimer.pro map8-phi.pro
   
   >IDL_7.0
   IDL> .r read-basic  (modify to read from directory where
                       blatpu.dat and blonpu.dat were created)
   IDL> .r read-weimer  (this reads the phi_weimer.inp that was just
                        generated; nt should be set <= the number 
                        of lines in weimer_input.inp)
   IDL> phi=phiw
   IDL> time=timew
   IDL> ntm=0
   IDL> iwin=0
   IDL> .r map8-phi

   Results should vary with ntm and look reasonable in terms of input 
   Bz and strength of the convection potential. That is, when Bz is large
   and negative, there should be lots of contour lines.

8. Use phi_weimer.inp as a SAMI3 input file. Set logical variables in 
   the SAMI3 namelist file to specify a SAMI3/Weimer05 run.
   It is recommended that this entire Weimer directory, particularly 
   the weimer_grid.dat and weimer_input.inp files, be copied and 
   saved for future reference.


