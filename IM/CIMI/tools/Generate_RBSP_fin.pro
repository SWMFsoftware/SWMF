pro Generate_RBSP_fin, $
;  This idl program opens ASCII files generated from CDAWeb
;  (http://cdaweb.gsfc.nasa.gov/istp_public/) for the Van Allen Probes
;  (RBSP) A/B Spacecraft and converts the measure fluxes to .fin files
;  for initializing CIMI.  The variables needed for generating ASCII
;  files to be processed are as follows: 
;   
;  1) Select Van Allen Probes (RBSP)
;
;  2) Select the following data products from either RBSPA/B for the
;     desired time frame
;        -  RBSPA_REL02_ECT-HOPE-SCI-L2SA
;        -  RBSPA_REL03_ECT-MAGEIS-L2
;        -  RBSPA_REL03_ECT-REPT-SCI-L2
;
;  3) Select "List Data (ASCII)..." and select "Output listing
;     times as year and seconds of year" and select the
;     following variable parameters for each intruement:
;
;     NOTE:  EACH FILE MUST BE GENERATED SEPARATELY WITH NO OTHER
;            INTRUMENT OR VARIABLES SPECIFIED!   
;
;     HOPE O+
;        - Calculated McIlwains L parameter (ion timebase)
;        - HOPE Spin averaged differential oxygen flux (as spectrogram)
;
;     HOPE H+
;        - Calculated McIlwains L parameter (ion timebase)
;        - HOPE Spin averaged differential proton flux 
;
;     HOPE e-
;        - Calculated McIlwains L parameter (electron timebase)
;        - HOPE Spin averaged differential electron flux 
;
;     MagEIS H+
;        - FPSA: Spin-Averaged Differential Proton Flux ~50-1500 keV 
;
;     MagEIS e-
;        - FESA: Spin-Averaged Differential Electron Flux ~20-4000 keV
;        - [Time Series] ---> Calculated McIlwains L parameter (Earths
;                             radii)
;
;     REPT e-
;        - FESA: Spin-Averaged Differential Electron Flux 2-20 MeV, in
;                12 energy bands (spectrogram) 
;        - [Time Series] ---> Calculated McIlwains L parameter (Earths radii) 
;
;  4) Save files to CIMI rundir/IM/  Default file names are:
;        - RBSP_HOPE_Oplus.txt
;        - RBSP_HOPE_Hplus.txt
;        - RBSP_HOPE_e.txt   
;        - RBSP_MagEIS_Hplus.txt
;        - RBSP_MagEIS_e.txt
;        - RBSP_REPT_e.txt
;
;  5) Run this program through IDL in rundir/IM to output CIMI
;     initialization files:
;        - RBSP_o.fin
;        - RBSP_h.fin
;        - RBSP_e.fin
;
;  Version History:
;   09/09/2015 - CMK added comments and full functionality to output
;                .fin files.


   
; VARIABLE DEFINITIONS:
;   
; INPUT: N/A

; OUTPUT:   
   HOPE_OPLUS_DATA=data_HOPE_Oplus, $
   HOPE_HPLUS_DATA=data_HOPE_Hplus, $
   HOPE_E_DATA=data_HOPE_e, $
   MAGEIS_HPLUS_DATA=data_MagEIS_Hplus, $
   MAGEIS_E_DATA=data_MagEIS_e, $
   REPT_E_DATA=data_REPT_e, $

; FLAGS:   
   DEBUG=debug_check

; Defines the default filenames to be read  
  filename_HOPE_Hplus="RBSP_HOPE_Hplus.txt"
  filename_HOPE_Oplus="RBSP_HOPE_Oplus.txt"
  filename_HOPE_e="RBSP_HOPE_e.txt"

  filename_MagEIS_Hplus="RBSP_MagEIS_Hplus.txt"
  filename_MagEIS_e="RBSP_MagEIS_e.txt"

  filename_REPT_e="RBSP_REPT_e.txt"

; Prompts the user if they want to change the default filename.
  filename_check=""
  print,"Default filename for HOPE H+: ",filename_HOPE_Hplus
  read,"Do you want to change it (y/n)? ",filename_check
  if (filename_check eq "y") then $
     read,"Enter filename for HOPE H+: ",filename_HOPE_Hplus

  print,"Default filename for HOPE O+: ",filename_HOPE_Oplus
  read,"Do you want to change it (y/n)? ",filename_check
  if (filename_check eq "y") then $
     read,"Enter filename for HOPE O+: ",filename_HOPE_Oplus

  print,"Default filename for HOPE e-: ",filename_HOPE_e
  read,"Do you want to change it (y/n)? ",filename_check
  if (filename_check eq "y") then $
     read,"Enter filename for HOPE e-: ",filename_HOPE_e
  
  print,"Default filename for MagEIS H+: ",filename_MagEIS_Hplus
  read,"Do you want to change it (y/n)? ",filename_check
  if (filename_check eq "y") then $
     read,"Enter filename for MagEIS H+: ",filename_MagEIS_Hplus

  print,"Default filename for MagEIS e-: ",filename_MagEIS_e
  read,"Do you want to change it (y/n)? ",filename_check
  if (filename_check eq "y") then $
     read,"Enter filename for MagEIS e-: ",filename_MagEIS_e
  
  print,"Default filename for REPT e-: ",filename_REPT_e
  read,"Do you want to change it (y/n)? ",filename_check
  if (filename_check eq "y") then $
     read,"Enter filename for REPT e-: ",filename_REPT_e

; Defines the length of the header and the data structure in
; preparation for reading the HOPE data files.
  
  HOPE_Oplus_header_length=60
  HOPE_Oplus_header=STRARR(HOPE_Oplus_header_length)

  HOPE_Hplus_header_length=60
  HOPE_Hplus_header=STRARR(HOPE_Hplus_header_length)

  HOPE_e_header_length=60
  HOPE_e_header=STRARR(HOPE_e_header_length)

  num_HOPE_Energy_Channels=72
  dataStruct_HOPE={year:0S,seconds:0.0D,L:0.0,delta_T:0.0, $
                   FLUXES:fltarr(num_HOPE_Energy_Channels), $
                   ENERGY:fltarr(num_HOPE_Energy_Channels)}

; Defines the length of the header and the data structure in
; preparation for reading the MagEIS H+ data files.
  
  MagEIS_Hplus_header_length=50
  MagEIS_Hplus_header=STRARR(MagEIS_Hplus_header_length)

  num_MagEIS_Hplus_Energy_channels=31
  dataStruct_MagEIS_Hplus={year:0S,seconds:0.0D, $
                           FLUXES:fltarr(num_MagEIS_Hplus_Energy_channels), $
                           ENERGY:fltarr(num_MagEIS_Hplus_Energy_channels)}

; Defines the length of the header and the data structure in
; preparation for reading the MagEIS e- data files.
  
  MagEIS_e_header_length=52
  MagEIS_e_header=STRARR(MagEIS_e_header_length)
  
  num_MagEIS_e_Energy_channels=25
  dataStruct_MagEIS_e={year:0S,seconds:0.0D,L:0.0, $
                       FLUXES:fltarr(num_MagEIS_e_Energy_channels), $
                       ENERGY:fltarr(num_MagEIS_e_Energy_channels)}

; REPT Energy channels converted from  MeV to keV
  num_REPT_e_Energy_channels=12
  REPT_e_Energy_arr= $
     [1.80E+00, 2.10E+00, 2.60E+00, 3.40E+00, 4.20E+00, $
      5.20E+00, 6.30E+00, 7.70E+00, 9.90E+00, 1.23E+01, $
      1.52E+01, 2.00E+01]*1000.
  
; Defines the length of the header and the data structure in
; preparation for reading the REPT e- data files.
  
  REPT_e_header_length=52
  REPT_e_header=STRARR(REPT_e_header_length)
  dataStruct_REPT_e={year:0S,seconds:0.0D,L:0.0, $
                     FLUXES:fltarr(num_REPT_e_Energy_channels)}

; CDAWeb data files have three lines as the footer of each file, which
; need to be substracted from the total length of the file for the
; appropriate length of the data files.
  
  num_lines_footer=3

; Closes all open files in preparation of reading the files.  
  close,/all
  file_index=50

; Finds the total number of lines in each ASCII data file and makes
; data structures of the appropriate size for each species and
; instrument.
  
  nrows_HOPE_Oplus=file_lines(filename_HOPE_Oplus)- $
                   HOPE_Oplus_header_length-num_lines_footer
  data_HOPE_Oplus=Replicate(dataStruct_HOPE,nrows_HOPE_Oplus)
               
  nrows_HOPE_Hplus=file_lines(filename_HOPE_Hplus)- $
                   HOPE_Hplus_header_length-num_lines_footer
  data_HOPE_Hplus=Replicate(dataStruct_HOPE,nrows_HOPE_Hplus)
               
  nrows_HOPE_e=file_lines(filename_HOPE_e)- $
                   HOPE_e_header_length-num_lines_footer
  data_HOPE_E=Replicate(dataStruct_HOPE,nrows_HOPE_e)
               
  nrows_MagEIS_Hplus=file_lines(filename_MagEIS_Hplus)- $
                   MagEIS_Hplus_header_length-num_lines_footer
  data_MagEIS_Hplus=Replicate(dataStruct_MagEIS_Hplus,nrows_MagEIS_Hplus)
               
  nrows_MagEIS_e=file_lines(filename_MagEIS_e)- $
                   MagEIS_e_header_length-num_lines_footer
  data_MagEIS_E=Replicate(dataStruct_MagEIS_e,nrows_MagEIS_e)
               
  nrows_REPT_e=file_lines(filename_REPT_e)- $
                   REPT_e_header_length-num_lines_footer
  data_REPT_E=Replicate(dataStruct_REPT_e,nrows_REPT_e)
               
; Opens each instrument and species file.
  
  openr,file_index,filename_HOPE_Oplus
  readf,file_index,HOPE_Oplus_header
  readf,file_index,data_HOPE_Oplus
  file_index=file_index+1

  openr,file_index,filename_HOPE_Hplus
  readf,file_index,HOPE_Hplus_header
  readf,file_index,data_HOPE_Hplus
  file_index=file_index+1

  openr,file_index,filename_HOPE_e
  readf,file_index,HOPE_e_header
  readf,file_index,data_HOPE_e
  file_index=file_index+1

  openr,file_index,filename_MagEIS_Hplus
  readf,file_index,MagEIS_Hplus_header
  readf,file_index,data_MagEIS_Hplus
  file_index=file_index+1

  openr,file_index,filename_MagEIS_e
  readf,file_index,MagEIS_e_header
  readf,file_index,data_MagEIS_e
  file_index=file_index+1

  openr,file_index,filename_REPT_e
  readf,file_index,REPT_e_header
  readf,file_index,data_REPT_e
  file_index=file_index+1

; Closes all open files.
  
  close,/all

; Extracts the time variable and number of time steps from each HOPE
; species.
  
  time_HOPE_Oplus=data_HOPE_Oplus.seconds
  num_time_HOPE_Oplus=size(time_HOPE_Oplus,/n_elements)

  time_HOPE_Hplus=data_HOPE_Hplus.seconds
  num_time_HOPE_Hplus=size(time_HOPE_Hplus,/n_elements)

  time_HOPE_e=data_HOPE_e.seconds
  num_time_HOPE_e=size(time_HOPE_e,/n_elements)

; Extracts the measured fluxes from the HOPE data structures.
  
  HOPE_Oplus_fluxes=data_HOPE_Oplus.fluxes
  HOPE_Hplus_fluxes=data_HOPE_Hplus.fluxes
  HOPE_e_fluxes=data_HOPE_e.fluxes

; Extracts the energy channels from the HOPE data structures and
; converts them from eV to keV.
  
  HOPE_Oplus_Energy=data_HOPE_Oplus.energy/1000.
  HOPE_Hplus_Energy=data_HOPE_Hplus.energy/1000.
  HOPE_e_Energy=data_HOPE_e.energy/1000.

; Instantiates the low energy flux arrays for each HOPE species.
  
  lowE_Oplus_fluxes=fltarr(num_time_HOPE_Oplus,num_HOPE_Energy_channels)
  lowE_Hplus_fluxes=fltarr(num_time_HOPE_Hplus,num_HOPE_Energy_channels)
  lowE_e_fluxes=fltarr(num_time_HOPE_e,num_HOPE_Energy_channels)

; Cycles through each measured flux for the HOPE instrument to throw
; out those measurements that are undefined, i.e., -1e-31. 
  
  for $
     n_time_HOPE_Oplus=0,num_time_HOPE_Oplus-1 $
  do $
     for $
     	n_HOPE_Energy_channel=0,num_HOPE_Energy_channels-1 $
     do $
        if $
     	   (HOPE_Oplus_fluxes(n_HOPE_Energy_Channel, $
                              n_time_HOPE_Oplus) gt -1e28) $
        then $
           lowE_Oplus_fluxes(n_time_HOPE_Oplus, $
                             n_HOPE_Energy_Channel)= $
     		HOPE_Oplus_fluxes(n_HOPE_Energy_Channel, $
                                  n_time_HOPE_Oplus) $
        else $
           lowE_Oplus_fluxes(n_time_HOPE_Oplus, $
                             n_HOPE_Energy_Channel)=0.

  for $
     n_time_HOPE_Hplus=0,num_time_HOPE_Hplus-1 $
  do $
     for $
     	n_HOPE_Energy_channel=0,num_HOPE_Energy_channels-1 $
     do $
        if $
     	   (HOPE_Hplus_fluxes(n_HOPE_Energy_Channel, $
                              n_time_HOPE_Hplus) gt -1e28) $
        then $
           lowE_Hplus_fluxes(n_time_HOPE_Hplus, $
                             n_HOPE_Energy_Channel)= $
     		HOPE_Hplus_fluxes(n_HOPE_Energy_Channel, $
                                  n_time_HOPE_Hplus) $
        else $
           lowE_Hplus_fluxes(n_time_HOPE_Hplus, $
                             n_HOPE_Energy_Channel)=0.

  for $
     n_time_HOPE_e=0,num_time_HOPE_e-1 $
  do $
     for $
     	n_HOPE_Energy_channel=0,num_HOPE_Energy_channels-1 $
     do $
        if $
     	   (HOPE_e_fluxes(n_HOPE_Energy_Channel, $
                              n_time_HOPE_e) gt -1e28) $
        then $
           lowE_e_fluxes(n_time_HOPE_e, $
                             n_HOPE_Energy_Channel)= $
     		HOPE_e_fluxes(n_HOPE_Energy_Channel, $
                                  n_time_HOPE_e) $
        else $
           lowE_e_fluxes(n_time_HOPE_e, $
                             n_HOPE_Energy_Channel)=0.

; Extracts the time variable and number of time steps from each MagEIS
; species.
  
  time_MagEIS_Hplus=data_MagEIS_Hplus.seconds
  num_time_MagEIS_Hplus=size(time_MagEIS_Hplus,/n_elements)

  time_MagEIS_e=data_MagEIS_Hplus.seconds
  num_time_MagEIS_e=size(time_MagEIS_e,/n_elements)

; Extracts the measured fluxes from the MagEIS data structures.

  MagEIS_Hplus_fluxes=data_MagEIS_Hplus.fluxes
  MagEIS_e_fluxes=data_MagEIS_e.fluxes

; Extracts the energy channels from the MagEIS data structures.
  
  MagEIS_Hplus_Energy=data_MagEIS_Hplus.energy
  MagEIS_e_Energy=data_MagEIS_e.energy
  
; Instantiates the mid energy flux arrays for each MagEIS species.
  
  midE_Hplus_fluxes= $
     fltarr(num_time_MagEIS_Hplus,num_MagEIS_Hplus_Energy_channels)
  midE_e_fluxes= $
     fltarr(num_time_MagEIS_e,num_MagEIS_e_Energy_channels)

; Cycles through each measured flux for the MagEIS instrument to throw
; out those measurements that are undefined, i.e., -1e-31. 
  
  for $
     n_time_MagEIS_Hplus=0,num_time_MagEIS_Hplus-1 $
  do $
     for $
     	n_MagEIS_Energy_channel=0,num_MagEIS_Hplus_Energy_channels-1 $
     do $
        if $
     	   (MagEIS_Hplus_fluxes(n_MagEIS_Energy_Channel, $
                                n_time_MagEIS_Hplus) gt -1e28) $
        then $
           midE_Hplus_fluxes(n_time_MagEIS_Hplus, $
                             n_MagEIS_Energy_Channel)= $
     		MagEIS_Hplus_fluxes(n_MagEIS_Energy_Channel, $
                                    n_time_MagEIS_Hplus) $
        else $
           midE_Hplus_fluxes(n_time_MagEIS_Hplus, $
                             n_MagEIS_Energy_Channel)=0.
  
  for $
     n_time_MagEIS_e=0,num_time_MagEIS_e-1 $
  do $
     for $
     	n_MagEIS_Energy_channel=0,num_MagEIS_e_Energy_channels-1 $
     do $
        if $
     	   (MagEIS_e_fluxes(n_MagEIS_Energy_Channel, $
                            n_time_MagEIS_e) gt -1e28) $
        then $
           midE_e_fluxes(n_time_MagEIS_e, $
                         n_MagEIS_Energy_Channel)= $
     		MagEIS_e_fluxes(n_MagEIS_Energy_Channel, $
                                n_time_MagEIS_e) $
        else $
           midE_e_fluxes(n_time_MagEIS_e, $
                         n_MagEIS_Energy_Channel)=0.

; Extracts the time variable and number of time steps for the REPT
; electrons.
    
  time_REPT_e=data_REPT_e.seconds
  num_time_REPT_e=size(time_REPT_e,/n_elements)

; Extracts the measured fluxes from the REPT data structure and
; converts from (cm^2 s sr MeV)^-1 to (cm^2 s sr keV)^-1.

  REPT_e_fluxes=data_REPT_e.fluxes/1000.

; Instantiates the relativistic energy flux arrays for REPT electrons.
  
  highE_e_fluxes= $
     fltarr(num_time_REPT_e,num_REPT_e_Energy_channels)

; Cycles through the measured fluxes for the REPT electrons to throw
; out those measurements that are undefined, i.e., -1e-31. 
  
  for $
     n_time_REPT_e=0L,num_time_REPT_e-1 $
  do $
     for $
     	n_REPT_Energy_channel=0,num_REPT_e_Energy_Channels-1 $
     do $
        if $
     	   (REPT_e_fluxes(n_REPT_Energy_Channel, $
                            n_time_REPT_e) gt -1e28) $
        then $
           highE_e_fluxes(n_time_REPT_e, $
                          n_REPT_Energy_channel)= $
     		REPT_e_fluxes(n_REPT_Energy_Channel, $
                              n_time_REPT_e) $
        else $
           highE_e_fluxes(n_time_REPT_e, $
                          n_REPT_Energy_Channel)=0.

; CIMI's Default Latitude grid.
  
  lat_CIMI_grid= $
     [11.812,13.777,15.742,17.705,19.665,21.622,23.576,25.527,27.473, $
      29.414,31.350,33.279,35.200,37.112,39.012,40.897,42.763,44.604, $
      46.409,48.163,49.837,51.382,52.725,53.823,54.720,55.488,56.175, $
      56.812,57.413,57.990,58.547,59.090,59.622,60.144,60.659,61.168, $
      61.671,62.170,62.666,63.159,63.649,64.137,64.624,65.109,65.593, $
      66.077,66.560,67.043,67.526,68.009,68.492,68.975,69.458]

; Converts magnetic latitude to L in a dipolar field
  
  L_CIMI_grid=cos(lat_CIMI_grid*!dtor)^(-2)
  
; Get indices for each species in each instrument's dataset to determine
; the array indices over which the flux will be output to the .fin
; file.
;  
; HOPE First  

  Lmin_HOPE_Oplus= $
     min(data_HOPE_Oplus[where(data_HOPE_Oplus[*].L gt -1e28)].L)
  imin_HOPE_Oplus= $
     min(where(data_HOPE_Oplus[*].L eq Lmin_HOPE_Oplus))
  Lmax_HOPE_Oplus= $
     max(data_HOPE_Oplus[where(data_HOPE_Oplus[*].L gt -1e28)].L)
  imax_HOPE_Oplus= $
     max(where(data_HOPE_Oplus[*].L eq Lmax_HOPE_Oplus))

  Lmin_HOPE_Hplus= $
     min(data_HOPE_Hplus[where(data_HOPE_Hplus[*].L gt -1e28)].L)
  imin_HOPE_Hplus= $
     min(where(data_HOPE_Hplus[*].L eq Lmin_HOPE_Hplus))
  Lmax_HOPE_Hplus= $
     max(data_HOPE_Hplus[where(data_HOPE_Hplus[*].L gt -1e28)].L)
  imax_HOPE_Hplus= $
     max(where(data_HOPE_Hplus[*].L eq Lmax_HOPE_Hplus))

  Lmin_HOPE_e= $
     min(data_HOPE_e[where(data_HOPE_e[*].L gt -1e28)].L)
  imin_HOPE_e= $
     min(where(data_HOPE_e[*].L eq Lmin_HOPE_e))
  Lmax_HOPE_e= $
     max(data_HOPE_e[where(data_HOPE_e[*].L gt -1e28)].L)
  imax_HOPE_e= $
     max(where(data_HOPE_e[*].L eq Lmax_HOPE_e))

;  MagEIS' turn.
;  
;  NOTE:
;     L is only defined in the MagEIS electron flux file. Therefore,
;     Lmin and Lmax are defined from the electron dataset, and the
;     proton flux will be collected from those times that are
;     relatively close to Lmin and Lmax. 

  Lmin_MagEIS_e= $
     min(data_MagEIS_e[where(data_MagEIS_e[*].L gt -1e28)].L)
  imin_MagEIS_e= $
     min(where(data_MagEIS_e[*].L eq Lmin_MagEIS_e))
  tmin_MagEIS_e= $
     data_MagEIS_e[imin_MagEIS_e].seconds
     
  Lmax_MagEIS_e= $
     max(data_MagEIS_e[where(data_MagEIS_e[*].L gt -1e28)].L)
  imax_MagEIS_e= $
     max(where(data_MagEIS_e[*].L eq Lmax_MagEIS_e))
  tmax_MagEIS_e= $
     data_MagEIS_e[imax_MagEIS_e].seconds

  imin_MagEIS_Hplus= $
     min(where(data_MagEIS_Hplus[*].seconds ge tmin_MagEIS_e))
  imax_MagEIS_Hplus= $
     min(where(data_MagEIS_Hplus[*].seconds ge tmax_MagEIS_e))

; And finally REPT
  
  Lmin_REPT_e= $
     min(data_REPT_e[where(data_REPT_e[*].L gt -1e28)].L)
  imin_REPT_e= $
     min(where(data_REPT_e[*].L eq Lmin_REPT_e))
     
  Lmax_REPT_e= $
     max(data_REPT_e[where(data_REPT_e[*].L gt -1e28)].L)
  imax_REPT_e= $
     max(where(data_REPT_e[*].L eq Lmax_REPT_e))

; Sorts the indices to account for RBSP being an inbound or outbound
; orbit.
  
  istart_HOPE_Oplus=min([imin_HOPE_Oplus,imax_HOPE_Oplus])
  istop_HOPE_Oplus=max([imin_HOPE_Oplus,imax_HOPE_Oplus])
  
  istart_HOPE_Hplus=min([imin_HOPE_Hplus,imax_HOPE_Hplus])
  istop_HOPE_Hplus=max([imin_HOPE_Hplus,imax_HOPE_Hplus])
  
  istart_HOPE_e=min([imin_HOPE_e,imax_HOPE_e])
  istop_HOPE_e=max([imin_HOPE_e,imax_HOPE_e])
  
  istart_MagEIS_Hplus=min([imin_MagEIS_Hplus,imax_MagEIS_Hplus])
  istop_MagEIS_Hplus=max([imin_MagEIS_Hplus,imax_MagEIS_Hplus])
  
  istart_MagEIS_e=min([imin_MagEIS_e,imax_MagEIS_e])
  istop_MagEIS_e=max([imin_MagEIS_e,imax_MagEIS_e])
  
  istart_REPT_e=min([imin_REPT_e,imax_REPT_e])
  istop_REPT_e=max([imin_REPT_e,imax_REPT_e])
  
; Finds the locations in the CIMI L grid that bounding the RBSP
; half-orbit.
  
  i_CIMI_Grid_start= $
     max(where(Lmin_HOPE_Oplus ge L_CIMI_grid))
  i_CIMI_Grid_end= $
     max(where(Lmax_HOPE_Oplus ge L_CIMI_grid))

; Instantiates the L_output grid from the minimum and maximum L of the
; RBSP half-orbit, bounding the L_CIMI_grid.
  
  num_L_fin_output=(i_Cimi_grid_end-i_CIMI_grid_start)+2
  L_fin_output=fltarr(num_L_fin_output)
  L_fin_output(0)= $
     Lmin_HOPE_Oplus
  L_fin_output(1:num_L_fin_output-2)= $
     L_CIMI_Grid(i_Cimi_grid_start+1:i_CIMI_grid_end)
  L_fin_output(num_L_fin_output-1)= $
     Lmax_HOPE_Oplus

; Defines the average flux for a given L location and the
; instrument's measured energies.
  
  average_fluxes_HOPE_Oplus= $
     fltarr(num_L_fin_output,num_HOPE_Energy_channels)

  average_fluxes_HOPE_Hplus= $
     fltarr(num_L_fin_output,num_HOPE_Energy_channels)

  average_fluxes_HOPE_e= $
     fltarr(num_L_fin_output,num_HOPE_Energy_channels)

  average_fluxes_MagEIS_Hplus= $
     fltarr(num_L_fin_output,num_MagEIS_Hplus_Energy_channels)

  average_fluxes_MagEIS_e= $
     fltarr(num_L_fin_output,num_MagEIS_e_Energy_channels)

  average_fluxes_REPT_e= $
     fltarr(num_L_fin_output,num_REPT_e_Energy_channels)

  
; Defines the bounds in L to cycle through the RBSP data.
  
  L_start=L_fin_output(0)
  L_stop= $
     (L_fin_output(0)+L_fin_output(1))/2.
  
; Obtains the indices to cycle over for the HOPE O+ data set.
  
  HOPE_Oplus_summation_indices= $
     where((data_HOPE_Oplus[istart_HOPE_Oplus:istop_HOPE_Oplus].L $
            ge L_start) and $
           (data_HOPE_Oplus[istart_HOPE_Oplus:istop_HOPE_Oplus].L $
            le L_stop))+istart_HOPE_Oplus

; Finds the average HOPE O+ flux for the given radial range and sets
; the minimum flux at 1e-10.
  
  for $
     num_HOPE_Energy=0,num_HOPE_Energy_channels-1 $
  do $
     average_fluxes_HOPE_Oplus(0,num_HOPE_Energy)= $
	mean(lowE_Oplus_Fluxes(HOPE_Oplus_summation_indices, $
                                num_HOPE_Energy))+1e-10

; Obtains the indices to cycle over for the HOPE H+ data set.

  HOPE_Hplus_summation_indices= $
     where((data_HOPE_Hplus[istart_HOPE_Hplus:istop_HOPE_Hplus].L $
            ge L_start) and $
           (data_HOPE_Hplus[istart_HOPE_Hplus:istop_HOPE_Hplus].L $
            le L_stop))+istart_HOPE_Hplus

; Finds the average HOPE H+ flux for the given radial range and sets
; the minimum flux at 1e-10. 
  
  for $
     num_HOPE_Energy=0,num_HOPE_Energy_channels-1 $
  do $
     average_fluxes_HOPE_Hplus(0,num_HOPE_Energy)= $
	mean(lowE_Hplus_Fluxes(HOPE_Hplus_summation_indices, $
                               num_HOPE_Energy))+1e-10

; Obtains the indices to cycle over for the HOPE e- data set.
  
  HOPE_e_summation_indices= $
     where((data_HOPE_e[istart_HOPE_e:istop_HOPE_e].L $
            ge L_start) and $
           (data_HOPE_e[istart_HOPE_e:istop_HOPE_e].L $
            le L_stop))+istart_HOPE_e

; Finds the average HOPE e- flux for the given radial range and sets
; the minimum flux at 1e-10.  

  for $
     num_HOPE_Energy=0,num_HOPE_Energy_channels-1 $
  do $
     average_fluxes_HOPE_e(0,num_HOPE_Energy)= $
	mean(lowE_e_Fluxes(HOPE_e_summation_indices, $
                           num_HOPE_Energy))+1e-10

; Obtains the indices to cycle over for the MagEIS e- data set.
  
  MagEIS_e_summation_indices= $
     where((data_MagEIS_e[istart_MagEIS_e:istop_MagEIS_e].L $
            ge L_start) and $
           (data_MagEIS_e[istart_MagEIS_e:istop_MagEIS_e].L $
            le L_stop))+istart_MagEIS_e

; Finds the average MagEIS e- flux for the given radial range and sets
; the minimum flux at 1e-10.  

  for $
     num_MagEIS_Energy=0,num_MagEIS_e_Energy_channels-1 $
  do $
     average_fluxes_MagEIS_e(0,num_MagEIS_Energy)= $
	mean(midE_e_Fluxes(MagEIS_e_summation_indices, $
                           num_MagEIS_Energy))+1e-10

; Finds the indices in the MagEIS H+ dataset corresponding to the
; minimum and maximum stop time.
  
  t_MagEIS_Hplus_start= $
     min(data_MagEIS_e[MagEIS_e_summation_indices].seconds)
  t_MagEIS_Hplus_stop= $
     max(data_MagEIS_e[MagEIS_e_summation_indices].seconds)
  
; Obtains the indices to cycle over for the MagEIS H+ data set.
  
  MagEIS_Hplus_summation_indices= $
     where((data_MagEIS_Hplus[istart_MagEIS_Hplus:istop_MagEIS_Hplus].seconds $
            ge t_MagEIS_Hplus_start) and $
           (data_MagEIS_Hplus[istart_MagEIS_Hplus:istop_MagEIS_Hplus].seconds $
            le t_MagEIS_Hplus_stop))+istart_MagEIS_Hplus

; Finds the average MagEIS H+ flux for the given radial range and sets
; the minimum flux at 1e-10.
  
  for $
     num_MagEIS_Energy=0,num_MagEIS_Hplus_Energy_channels-1 $
  do $
     average_fluxes_MagEIS_Hplus(0,num_MagEIS_Energy)= $
	mean(midE_Hplus_Fluxes(MagEIS_Hplus_summation_indices, $
                               num_MagEIS_Energy))+1e-10

; Obtains the indices to cycle over for the REPT e- data set.

  REPT_e_summation_indices= $
     where((data_REPT_e[istart_REPT_e:istop_REPT_e].L $
            ge L_start) and $
           (data_REPT_e[istart_REPT_e:istop_REPT_e].L $
            le L_stop))+istart_REPT_e

; Finds the average REPT e- flux for the given radial range and sets
; the minimum flux at 1e-10.
  
  for $
     num_REPT_Energy=0,num_REPT_e_Energy_channels-1 $
  do $
     average_fluxes_REPT_e(0,num_REPT_Energy)= $
	mean(highE_e_Fluxes(REPT_e_summation_indices, $
                           num_REPT_Energy))+1e-10

; For loop that cycles over the L_fin_output grid.

  for $
     num_L_loc=1,num_L_fin_output-2 $
  do begin $

; Defines the bounds in L to cycle through the RBSP data.
     
     L_start=(L_fin_output(num_L_loc-1)+L_fin_output(num_L_loc))/2.
     L_stop=(L_fin_output(num_L_loc+1)+L_fin_output(num_L_loc))/2.
     
; Obtains the indices to cycle over for the HOPE O+ data set.
  
     HOPE_Oplus_summation_indices= $
        where((data_HOPE_Oplus[istart_HOPE_Oplus:istop_HOPE_Oplus].L $
               ge L_start) and $
              (data_HOPE_Oplus[istart_HOPE_Oplus:istop_HOPE_Oplus].L $
               le L_stop))+istart_HOPE_Oplus

; Finds the average HOPE O+ flux for the given radial range and sets
; the minimum flux at 1e-10.  

     for $
        num_HOPE_Energy=0,num_HOPE_Energy_channels-1 $
     do $
        average_fluxes_HOPE_Oplus(num_L_loc,num_HOPE_Energy)= $
           mean(lowE_Oplus_Fluxes(HOPE_Oplus_summation_indices, $
                                   num_HOPE_Energy))+1e-10

; Obtains the indices to cycle over for the HOPE H+ data set.
  
     HOPE_Hplus_summation_indices= $
        where((data_HOPE_Hplus[istart_HOPE_Hplus:istop_HOPE_Hplus].L $
               ge L_start) and $
              (data_HOPE_Hplus[istart_HOPE_Hplus:istop_HOPE_Hplus].L $
               le L_stop))+istart_HOPE_Hplus

; Finds the average HOPE H+ flux for the given radial range and sets
; the minimum flux at 1e-10.  

     for $
        num_HOPE_Energy=0,num_HOPE_Energy_channels-1 $
     do $
        average_fluxes_HOPE_Hplus(num_L_loc,num_HOPE_Energy)= $
           mean(lowE_Hplus_Fluxes(HOPE_Hplus_summation_indices, $
                                  num_HOPE_Energy))+1e-10

; Obtains the indices to cycle over for the HOPE e- data set.
  
     HOPE_e_summation_indices= $
        where((data_HOPE_e[istart_HOPE_e:istop_HOPE_e].L $
               ge L_start) and $
              (data_HOPE_e[istart_HOPE_e:istop_HOPE_e].L $
               le L_stop))+istart_HOPE_e

; Finds the average HOPE e- flux for the given radial range and sets
; the minimum flux at 1e-10.  

     for $
        num_HOPE_Energy=0,num_HOPE_Energy_channels-1 $
     do $
        average_fluxes_HOPE_e(num_L_loc,num_HOPE_Energy)= $
           mean(lowE_e_Fluxes(HOPE_e_summation_indices, $
                              num_HOPE_Energy))+1e-10

; Obtains the indices to cycle over for the MagEIS H+ data set.
  
     MagEIS_e_summation_indices= $
        where((data_MagEIS_e[istart_MagEIS_e:istop_MagEIS_e].L $
               ge L_start) and $
              (data_MagEIS_e[istart_MagEIS_e:istop_MagEIS_e].L $
               le L_stop))+istart_MagEIS_e

; Finds the average MagEIS e- flux for the given radial range and sets
; the minimum flux at 1e-10.  

     for $
        num_MagEIS_Energy=0,num_MagEIS_e_Energy_channels-1 $
     do $
        average_fluxes_MagEIS_e(num_L_loc,num_MagEIS_Energy)= $
           mean(midE_e_Fluxes(MagEIS_e_summation_indices, $
                              num_MagEIS_Energy))+1e-10

; Finds the indices in the MagEIS H+ dataset corresponding to the
; minimum and maximum stop time.
  
     t_MagEIS_Hplus_start= $
        min(data_MagEIS_e[MagEIS_e_summation_indices].seconds)
     t_MagEIS_Hplus_stop= $
        max(data_MagEIS_e[MagEIS_e_summation_indices].seconds)
     
; Obtains the indices to cycle over for the MagEIS H+ data set.
  
     MagEIS_Hplus_summation_indices= $
        where((data_MagEIS_Hplus[istart_MagEIS_Hplus:istop_MagEIS_Hplus].seconds $
               ge t_MagEIS_Hplus_start) and $
              (data_MagEIS_Hplus[istart_MagEIS_Hplus:istop_MagEIS_Hplus].seconds $
               le t_MagEIS_Hplus_stop))+istart_MagEIS_Hplus

; Finds the average MagEIS H+ flux for the given radial range and sets
; the minimum flux at 1e-10.  

     for $
        num_MagEIS_Energy=0,num_MagEIS_Hplus_Energy_channels-1 $
     do $
        average_fluxes_MagEIS_Hplus(num_L_loc,num_MagEIS_Energy)= $
           mean(midE_Hplus_Fluxes(MagEIS_Hplus_summation_indices, $
                                  num_MagEIS_Energy))+1e-10

; Obtains the indices to cycle over for the REPT e- data set.
  
     REPT_e_summation_indices= $
        where((data_REPT_e[istart_REPT_e:istop_REPT_e].L $
               ge L_start) and $
              (data_REPT_e[istart_REPT_e:istop_REPT_e].L $
               le L_stop))+istart_REPT_e

; Finds the average REPT e- flux for the given radial range and sets
; the minimum flux at 1e-10.  

     for $
        num_REPT_Energy=0,num_REPT_e_Energy_channels-1 $
     do $
        average_fluxes_REPT_e(num_L_loc,num_REPT_Energy)= $
           mean(highE_e_Fluxes(REPT_e_summation_indices, $
                              num_REPT_Energy))+1e-10

  endfor

  L_start=L_stop
  L_stop= $
     L_fin_output(num_L_fin_output-1)

; Obtains the indices to cycle over for the HOPE O+ data set.
  
  HOPE_Oplus_summation_indices= $
     where((data_HOPE_Oplus[istart_HOPE_Oplus:istop_HOPE_Oplus].L $
            ge L_start) and $
           (data_HOPE_Oplus[istart_HOPE_Oplus:istop_HOPE_Oplus].L $
            le L_stop))+istart_HOPE_Oplus

; Finds the average HOPE O+ flux for the given radial range and sets
; the minimum flux at 1e-10.  

  for $
     num_HOPE_Energy=0,num_HOPE_Energy_channels-1 $
  do $
     average_fluxes_HOPE_Oplus(num_L_fin_output-1,num_HOPE_Energy)= $
     	mean(lowE_Oplus_Fluxes(HOPE_Oplus_summation_indices, $
                                num_HOPE_Energy))+1e-10

; Obtains the indices to cycle over for the HOPE H+ data set.
  
  HOPE_Hplus_summation_indices= $
     where((data_HOPE_Hplus[istart_HOPE_Hplus:istop_HOPE_Hplus].L $
            ge L_start) and $
           (data_HOPE_Hplus[istart_HOPE_Hplus:istop_HOPE_Hplus].L $
            le L_stop))+istart_HOPE_Hplus

; Finds the average HOPE H+ flux for the given radial range and sets
; the minimum flux at 1e-10.  

  for $
     num_HOPE_Energy=0,num_HOPE_Energy_channels-1 $
  do $
     average_fluxes_HOPE_Hplus(num_L_fin_output-1,num_HOPE_Energy)= $
     	mean(lowE_Hplus_Fluxes(HOPE_Hplus_summation_indices, $
                                num_HOPE_Energy))+1e-10

; Obtains the indices to cycle over for the HOPE e- data set.
  
  HOPE_e_summation_indices= $
     where((data_HOPE_e[istart_HOPE_e:istop_HOPE_e].L $
            ge L_start) and $
           (data_HOPE_e[istart_HOPE_e:istop_HOPE_e].L $
            le L_stop))+istart_HOPE_e

; Finds the average HOPE e- flux for the given radial range and sets
; the minimum flux at 1e-10.  

  for $
     num_HOPE_Energy=0,num_HOPE_Energy_channels-1 $
  do $
     average_fluxes_HOPE_e(num_L_fin_output-1,num_HOPE_Energy)= $
     	mean(lowE_e_Fluxes(HOPE_e_summation_indices, $
                                num_HOPE_Energy))+1e-10

; Obtains the indices to cycle over for the MagEIS e- data set.
  
  MagEIS_e_summation_indices= $
     where((data_MagEIS_e[istart_MagEIS_e:istop_MagEIS_e].L $
            ge L_start) and $
           (data_MagEIS_e[istart_MagEIS_e:istop_MagEIS_e].L $
            le L_stop))+istart_MagEIS_e

; Finds the average MagEIS e- flux for the given radial range and sets
; the minimum flux at 1e-10.  

  for $
     num_MagEIS_Energy=0,num_MagEIS_e_Energy_channels-1 $
  do $
     average_fluxes_MagEIS_e(num_L_fin_output-1,num_MagEIS_Energy)= $
     	mean(midE_e_Fluxes(MagEIS_e_summation_indices, $
                           num_MagEIS_Energy))+1e-10
  
; Finds the indices in the MagEIS H+ dataset corresponding to the
; minimum and maximum stop time.
  
  t_MagEIS_Hplus_start= $
     min(data_MagEIS_e[MagEIS_e_summation_indices].seconds)
  t_MagEIS_Hplus_stop= $
     max(data_MagEIS_e[MagEIS_e_summation_indices].seconds)
     
; Obtains the indices to cycle over for the MagEIS H+ data set.
  
  MagEIS_Hplus_summation_indices= $
     where((data_MagEIS_Hplus[istart_MagEIS_Hplus:istop_MagEIS_Hplus].seconds $
            ge t_MagEIS_Hplus_start) and $
           (data_MagEIS_Hplus[istart_MagEIS_Hplus:istop_MagEIS_Hplus].seconds $
            le t_MagEIS_Hplus_stop))+istart_MagEIS_Hplus

; Finds the average MagEIS H+ flux for the given radial range and sets
; the minimum flux at 1e-10.  

  for $
     num_MagEIS_Energy=0,num_MagEIS_Hplus_Energy_channels-1 $
  do $
     average_fluxes_MagEIS_Hplus(num_L_fin_output-1,num_MagEIS_Energy)= $
     	mean(midE_Hplus_Fluxes(MagEIS_Hplus_summation_indices, $
                           num_MagEIS_Energy))+1e-10

; Obtains the indices to cycle over for the REPT e- data set.

  REPT_e_summation_indices= $
     where((data_REPT_e[istart_REPT_e:istop_REPT_e].L $
            ge L_start) and $
           (data_REPT_e[istart_REPT_e:istop_REPT_e].L $
            le L_stop))+istart_REPT_e
  
; Finds the average REPT e- flux for the given radial range and sets
; the minimum flux at 1e-10.  

  for $
     num_REPT_Energy=0,num_REPT_e_Energy_channels-1 $
  do $
     average_fluxes_REPT_e(num_L_fin_output-1,num_REPT_Energy)= $
     	mean(highE_e_Fluxes(REPT_e_summation_indices, $
                           num_REPT_Energy))+1e-10

; Sets up the energy output arrays for each species.  The first four
; channels of MagEIS H+ data are ignored because they are orders of
; magnitude higher than HOPE's measured fluxes in the highest
; energy channels.  The last few channels are undefined.  Likewise,
; MagEIS e- fluxes for the first 6 channels are undefined and overlap
; with HOPE and the last several overlap with REPT.  Channel 17 is
; ignored as it is undefined.
  
  Oplus_Energy_Output=HOPE_Oplus_Energy(*,0)
  Hplus_Energy_Output=[HOPE_Hplus_Energy(*,0), $
                       MagEIS_Hplus_Energy(4:20,0)]
  e_Energy_Output=[HOPE_e_Energy(*,0), $
                   MagEIS_e_Energy(6:15,0), $
                   MagEIS_e_Energy(18:19,0), $
                   REPT_e_energy_arr]

  print,"O+ output array: ",Oplus_Energy_Output
  print,"H+ output array: ",Hplus_Energy_Output
  print,"e- output array: ",e_Energy_Output
  
; The number of points in the L array.
  
  IL=num_L_fin_output

; The number of points in the each species' energy output array.
  
  IE_Oplus=size(Oplus_Energy_Output,/n_elements)
  IE_Hplus=size(Hplus_Energy_Output,/n_elements)
  IE_e=size(e_Energy_Output,/n_elements)

; iunit is always 1 since the average RBSP fluxes are measured in
; (cm^2 s sr keV)^{-1}.
  
  iunit=1

; Closes all open files.
  
  close,/all

; Outputs the O+ average flux in the standard .fin format.
  
  openw,file_index,"RBSP_o.fin"
  printf,file_index, $
         ' '+STRTRIM(STRING(IL),2)+' '+STRTRIM(STRING(IE_Oplus),2)+ $
         '    ! IL, IE, dimension of RBSP flux'
  printf,file_index,'  '+STRTRIM(STRING(iunit),2)+ $
         '       ! iunit: 1=CCE flux in (cm^2 s sr keV)^-1, ' + $
         '2=flux in (cm^2 s MeV)^-1'
  printf,file_index,L_fin_output
  printf,file_index,Oplus_Energy_Output
  for num_energy_HOPE=0,num_HOPE_Energy_Channels-1 do $
     printf,file_index,average_fluxes_HOPE_Oplus(*,num_Energy_HOPE)
  file_index=file_index+1

; Outputs the H+ average flux in the standard .fin format.

  openw,file_index,"RBSP_h.fin"
  printf,file_index, $
         ' '+STRTRIM(STRING(IL),2)+' '+STRTRIM(STRING(IE_Hplus),2)+ $
         '    ! IL, IE, dimension of RBSP flux'
  printf,file_index,'  '+STRTRIM(STRING(iunit),2)+ $
         '       ! iunit: 1=CCE flux in (cm^2 s sr keV)^-1, ' + $
         '2=flux in (cm^2 s MeV)^-1'
  printf,file_index,L_fin_output
  printf,file_index,Hplus_Energy_Output
  for num_energy_HOPE=0,num_HOPE_Energy_Channels-1 do $
     printf,file_index,average_fluxes_HOPE_Hplus(*,num_Energy_HOPE)
  for num_energy_MagEIS_Hplus=4,20 do $
     printf, $
        file_index,average_fluxes_MagEIS_Hplus(*,num_Energy_MagEIS_Hplus)
  file_index=file_index+1

; Outputs the e- average flux in the standard .fin format.

  openw,file_index,"RBSP_e.fin"
  printf,file_index, $
         ' '+STRTRIM(STRING(IL),2)+' '+STRTRIM(STRING(IE_e),2)+ $
         '    ! IL, IE, dimension of RBSP flux'
  printf,file_index,'  '+STRTRIM(STRING(iunit),2)+ $
         '       ! iunit: 1=CCE flux in (cm^2 s sr keV)^-1, ' + $
         '2=flux in (cm^2 s MeV)^-1'
  printf,file_index,L_fin_output
  printf,file_index,e_Energy_Output
  for num_energy_HOPE=0,num_HOPE_Energy_Channels-1 do $
     printf,file_index,average_fluxes_HOPE_e(*,num_Energy_HOPE)
  for num_energy_MagEIS_e=6,15 do $
     printf, $
        file_index,average_fluxes_MagEIS_e(*,num_Energy_MagEIS_e)
  printf,file_index,average_fluxes_MagEIS_e(*,18)
  printf,file_index,average_fluxes_MagEIS_e(*,19)
  for num_energy_REPT_e=0,num_REPT_e_Energy_channels-1 do $
     printf, $
        file_index,average_fluxes_REPT_e(*,num_Energy_REPT_e)
  file_index=file_index+1

; Closes all open files.
  
  close,/all
  
; The debug flag routine which plots several diagnostic information of
; the ASCII files to make sure they were correctly read in.  
  
  if $
     keyword_set(debug_check) $
  then begin

; Defines the debugging string to be read
     
     debug_string=''

     wset,0
     w,1,1
     plot,lat_CIMI_grid,L_CIMI_grid, $
          xtitle='Latitude (Deg.)',xr=[0,90], $
          ytitle=TEXTOIDL('L (R_E)'),yr=[0,10], $
          psym=-6
     makelinex,Lmin_HOPE_Oplus,linesytle=2
     makelinex,Lmax_HOPE_Oplus,linesytle=2
     makelinex,Lmin_HOPE_Hplus,linesytle=2
     makelinex,Lmax_HOPE_Hplus,linesytle=2
     makelinex,Lmin_HOPE_e,linesytle=2
     makelinex,Lmax_HOPE_e,linesytle=2
     
     print,''
     print,'Lmin_HOPE_Oplus (R_E): ',Lmin_HOPE_Oplus
     print,'imin_HOPE_Oplus: ',imin_HOPE_Oplus
     print,'Lmax_HOPE_Oplus (R_E): ',Lmax_HOPE_Oplus
     print,'imax_HOPE_Oplus: ',imax_HOPE_Oplus
     print,'istart_HOPE_Oplus: ',istart_HOPE_Oplus
     print,'istop_HOPE_Oplus: ',istop_HOPE_Oplus
     
     print,''
     print,'Lmin_HOPE_Hplus (R_E): ',Lmin_HOPE_Hplus
     print,'imin_HOPE_Hplus: ',imin_HOPE_Hplus
     print,'Lmax_HOPE_Hplus (R_E): ',Lmax_HOPE_Hplus
     print,'imax_HOPE_Hplus: ',imax_HOPE_Hplus
     print,'istart_HOPE_Hplus: ',istart_HOPE_Hplus
     print,'istop_HOPE_Hplus: ',istop_HOPE_Hplus
     
     print,''
     print,'Lmin_HOPE_e (R_E): ',Lmin_HOPE_e
     print,'imin_HOPE_e: ',imin_HOPE_e
     print,'Lmax_HOPE_e (R_E): ',Lmax_HOPE_e
     print,'imax_HOPE_e: ',imax_HOPE_e
     print,'istart_HOPE_e: ',istart_HOPE_e
     print,'istop_HOPE_e: ',istop_HOPE_e
     
     print,''
     print,'Lmin_MagEIS_e (R_E): ',Lmin_MagEIS_e
     print,'imin_MagEIS_e: ',imin_MagEIS_e
     print,'Lmax_MagEIS_e (R_E): ',Lmax_MagEIS_e
     print,'imax_MagEIS_e: ',imax_MagEIS_e
     print,'istart_MagEIS_e: ',istart_MagEIS_e
     print,'istop_MagEIS_e: ',istop_MagEIS_e
     
     print,''
     print,'tmin_MagEIS_e (s): ',tmin_MagEIS_e
     print,'imin_MagEIS_Hplus: ',imin_MagEIS_Hplus
     print,'tmin_MagEIS_Hplus (s): ', $
           data_MagEIS_Hplus[imin_MagEIS_Hplus].seconds
     
     print,''
     print,'tmax_MagEIS_e (s): ',tmax_MagEIS_e
     print,'imax_MagEIS_Hplus: ',imax_MagEIS_Hplus
     print,'tmax_MagEIS_Hplus (s): ', $
           data_MagEIS_Hplus[imax_MagEIS_Hplus].seconds

     print,'istart_MagEIS_Hplus: ',istart_MagEIS_Hplus
     print,'istop_MagEIS_Hplus: ',istop_MagEIS_Hplus

     print,''
     print,'Lmin_REPT_e (R_E): ',Lmin_REPT_e
     print,'imin_REPT_e: ',imin_REPT_e
     print,'Lmax_REPT_e (R_E): ',Lmax_REPT_e
     print,'imax_REPT_e: ',imax_REPT_e

     print,'istart_REPT_e: ',istart_REPT_e
     print,'istop_REPT_e: ',istop_REPT_e

     print,''
     print,'L_CIMI_min(n), L_RBSP, L_CIMI_min(n+1): '
     print,L_CIMI_grid(i_CIMI_grid_start),Lmin_HOPE_Oplus, $
           L_CIMI_grid(i_CIMI_grid_start+1)
     print,''
     print,'L_CIMI_max(n), L_RBSP, L_CIMI_max(n+1): '
     print,L_CIMI_grid(i_CIMI_grid_end),Lmax_HOPE_Oplus, $
           L_CIMI_grid(i_CIMI_grid_end+1)
     print,'L_fin_output: '
     print,L_fin_output

     char_size=2.5

     window,0,title='HOPE Data'
     loadct,13,/silent
     w,1,5
     plot,data_HOPE_Oplus.seconds,data_HOPE_Oplus.L,yr=[1,7]
     plot,data_HOPE_e.seconds,data_HOPE_e.L,yr=[1,7]
     image_cont,alog10(lowE_Oplus_fluxes), $
                data_HOPE_Oplus.seconds,HOPE_Oplus_Energy(*,0), $
                /ylog,minval=0, $
                charsize=char_size
     image_cont,alog10(lowE_Hplus_fluxes), $
                data_HOPE_Hplus.seconds,HOPE_Hplus_Energy(*,0), $
                /ylog,minval=0, $
                charsize=char_size
     image_cont,alog10(lowE_e_fluxes), $
                data_HOPE_e.seconds,HOPE_e_Energy(*,0), $
                /ylog,minval=0, $
                charsize=char_size
     
     window,1,title='MagEIS Data'
     loadct,13,/silent
     w,1,3
     image_cont,alog10(midE_e_fluxes), $
                data_MagEIS_e.seconds,MagEIS_e_Energy(*,0), $
                /ylog,yr=[20.,2500],minval=0, $
                charsize=char_size
     image_cont,alog10(midE_Hplus_fluxes), $
                data_MagEIS_Hplus.seconds,MagEIS_Hplus_Energy(*,0), $
                /ylog,minval=0, $
                charsize=char_size
     plot,data_MagEIS_e.seconds,data_MagEIS_e.L,yr=[1,7], $
          charsize=char_size
     
     window,2,title='REPT Data'
     w,1,2
     loadct,13,/silent
     image_cont,alog10(highE_e_fluxes*1000.), $
                time_REPT_e,REPT_e_Energy_arr,minval=0
                
     plot,time_REPT_e,data_REPT_e[*].L,yr=[1,6]

     n_time_max=2
     legend_pos_arr=[1e6,7.]
     window,4,title='Species Flux v. Energy'
     loadct,0,/silent
     tvlct,000,000,255,102
     tvlct,255,000,000,101
     tvlct,255,255,255,100
     w,1,3
     plot,HOPE_Oplus_Energy(*,0),alog10(lowE_Oplus_fluxes(0,*)), $
          xr=[.001,60000000],/xlog,xtitle='Energy (keV)', $
          yr=[-2,10],ytitle=TEXTOIDL('log_{10}(O^{+} Flux) ' + $
                                     '(cm^{2} s sr keV)^{-1}'), $
          psym=6,charsize=char_size,/nodata
     legend,['HOPE','MagEIS','REPT', $
             'L='+STRTRIM(STRING(data_HOPE_Oplus[istart_HOPE_Oplus].L, $
                                 FORMAT='(F4.2)'),2), $
             'L='+STRTRIM(STRING(data_HOPE_Oplus[istart_HOPE_Oplus+1].L, $
                                 FORMAT='(F4.2)'),2), $
             'L='+STRTRIM(STRING(data_HOPE_Oplus[istart_HOPE_Oplus+2].L, $
                                 FORMAT='(F4.2)'),2)], $
            psym=[6,5,4,0,0,0],color=[100,100,100,100,101,102], $
            pos=legend_pos_arr
     for $
        n_time=0,n_time_max $
     do $
        oplot,HOPE_Oplus_Energy(*,n_time), $
              alog10(lowE_Oplus_fluxes(n_time,*)), $
              psym=6,color=100+n_time
     
     plot,HOPE_Hplus_Energy(*,0),alog10(lowE_Hplus_fluxes(0,*)), $
          xr=[.001,60000000],/xlog,xtitle='Energy (keV)', $
          yr=[-2,10],ytitle=TEXTOIDL('log_{10}(H^{+} Flux) ' + $
                                     '(cm^{2} s sr keV)^{-1}'), $
          psym=6,charsize=char_size,/nodata
     legend,['HOPE','MagEIS','REPT', $
             'L='+STRTRIM(STRING(data_HOPE_Oplus[0].L,FORMAT='(F4.2)'),2), $
             'L='+STRTRIM(STRING(data_HOPE_Oplus[0].L,FORMAT='(F4.2)'),2), $
             'L='+STRTRIM(STRING(data_HOPE_Oplus[0].L,FORMAT='(F4.2)'),2)], $
            psym=[6,5,4,0,0,0],color=[100,100,100,100,101,102], $
            pos=legend_pos_arr
     for $
        n_time=0,n_time_max $
     do begin
        
        oplot,HOPE_Hplus_Energy(*,n_time), $
              alog10(lowE_Hplus_fluxes(n_time,*)), $
              psym=6,color=100+n_time
        oplot,MagEIS_Hplus_Energy(*,n_time), $
              alog10(midE_Hplus_fluxes(n_time,*)), $
              psym=5,color=100+n_time
        
     endfor
     
     plot,HOPE_e_Energy(*,0), $
          alog10(lowE_e_fluxes(0,*)), $
          xr=[.001,60000000],/xlog,xtitle='Energy (keV)', $
          yr=[-2,10],ytitle=TEXTOIDL('log_{10}(e^{-} Flux) ' + $
                                    '(cm^{2} s sr keV)^{-1}'), $
          psym=6,charsize=char_size,/nodata
     legend,['HOPE','MagEIS','REPT', $
             'L='+STRTRIM(STRING(data_HOPE_Oplus[0].L,FORMAT='(F4.2)'),2), $
             'L='+STRTRIM(STRING(data_HOPE_Oplus[0].L,FORMAT='(F4.2)'),2), $
             'L='+STRTRIM(STRING(data_HOPE_Oplus[0].L,FORMAT='(F4.2)'),2)], $
            psym=[6,5,4,0,0,0],color=[100,100,100,100,101,102], $
            pos=legend_pos_arr
     for $
        n_time=0,n_time_max $
     do begin
        
        oplot,HOPE_e_Energy(*,n_time), $
              alog10(lowE_e_fluxes(n_time,*)), $
              psym=6,color=100+n_time
        oplot,MagEIS_e_Energy(*,n_time), $
              alog10(midE_e_fluxes(n_time,*)), $
              psym=5,color=100+n_time
        oplot,REPT_e_Energy_arr, $
              alog10(highE_e_fluxes(n_time,*)), $
              psym=4,color=100+n_time
        
     endfor
     
     window,7,xsize=1280,ysize=1600,title='.fin Map'
     loadct,13,/silent
     w,1,6
     char_size=3
     image_cont,alog10(average_fluxes_HOPE_Oplus), $
                L_fin_output,HOPE_Oplus_Energy(*,0),minval=-2, $
                xtitle=TEXTOIDL('L (R_E)'), $
                ytitle=TEXTOIDL('Energy (keV)'),/ylog, $
                charsize=char_size
     image_cont,alog10(average_fluxes_HOPE_Hplus), $
                L_fin_output,HOPE_Hplus_Energy(*,0),minval=-2, $
                xtitle=TEXTOIDL('L (R_E)'), $
                ytitle=TEXTOIDL('Energy (keV)'),/ylog, $
                charsize=char_size
     image_cont,alog10(average_fluxes_HOPE_e), $
                L_fin_output,HOPE_e_Energy(*,0),minval=-2, $
                xtitle=TEXTOIDL('L (R_E)'), $
                ytitle=TEXTOIDL('Energy (keV)'),/ylog, $
                charsize=char_size
     
     image_cont,alog10(average_fluxes_MagEIS_Hplus), $
                L_fin_output,MagEIS_Hplus_Energy(*,0),minval=-2, $
                xtitle=TEXTOIDL('L (R_E)'), $
                ytitle=TEXTOIDL('Energy (keV)'),/ylog, $
                charsize=char_size
     image_cont,alog10(average_fluxes_MagEIS_e), $
                L_fin_output,MagEIS_e_Energy(*,0),minval=-2, $
                xtitle=TEXTOIDL('L (R_E)'), $
                ytitle=TEXTOIDL('Energy (keV)'),/ylog, $
                charsize=char_size
     
     image_cont,alog10(average_fluxes_REPT_e), $
                L_fin_output,REPT_e_Energy_arr,minval=-2, $
                xtitle=TEXTOIDL('L (R_E)'), $
                ytitle=TEXTOIDL('Energy (keV)'),/ylog, $
                charsize=char_size
     
     set_plot,'ps'
     device,filename='RBSP_Fluxes.eps', $
            xsize=8.50,xoffset=1., $
            ysize=11.0,yoffset=1., $
            /inches,/encapsulated,/portrait,/color,bits_per_pixel=8
     
     char_size=2.5
     legend_pos_arr=[1e5,7.]
     loadct,0,/silent
     tvlct,000,000,255,102
     tvlct,255,000,000,101
     tvlct,000,000,000,100
     w,1,3
     plot,HOPE_Oplus_Energy(*,0),alog10(lowE_Oplus_fluxes(0,*)), $
          xr=[.001,60000000],/xlog,xtitle='Energy (keV)', $
          yr=[0,10],ytitle=TEXTOIDL('log_{10}(O^{+} Flux) ' + $
                                     '(cm^{2} s sr keV)^{-1}'), $
          psym=6,charsize=char_size,/nodata
     legend,['HOPE','MagEIS','REPT', $
             'L='+STRTRIM(STRING(data_HOPE_Oplus[0].L,FORMAT='(F4.2)'),2), $
             'L='+STRTRIM(STRING(data_HOPE_Oplus[0].L,FORMAT='(F4.2)'),2), $
             'L='+STRTRIM(STRING(data_HOPE_Oplus[0].L,FORMAT='(F4.2)'),2)], $
            psym=[6,5,4,0,0,0],color=[100,100,100,100,101,102], $
            pos=legend_pos_arr
     for $
        n_time=0,n_time_max $
     do $
        oplot,HOPE_Oplus_Energy(*,n_time), $
              alog10(lowE_Oplus_fluxes(n_time,*)), $
              psym=6,color=100+n_time
     
     plot,HOPE_Hplus_Energy(*,0),alog10(lowE_Hplus_fluxes(0,*)), $
          xr=[.001,60000000],/xlog,xtitle='Energy (keV)', $
          yr=[0,10],ytitle=TEXTOIDL('log_{10}(H^{+} Flux) ' + $
                                     '(cm^{2} s sr keV)^{-1}'), $
          psym=6,charsize=char_size,/nodata
     legend,['HOPE','MagEIS','REPT', $
             'L='+STRTRIM(STRING(data_HOPE_Oplus[0].L,FORMAT='(F4.2)'),2), $
             'L='+STRTRIM(STRING(data_HOPE_Oplus[0].L,FORMAT='(F4.2)'),2), $
             'L='+STRTRIM(STRING(data_HOPE_Oplus[0].L,FORMAT='(F4.2)'),2)], $
            psym=[6,5,4,0,0,0],color=[100,100,100,100,101,102], $
            pos=legend_pos_arr
     for $
        n_time=0,n_time_max $
     do begin
        
        oplot,HOPE_Hplus_Energy(*,n_time), $
              alog10(lowE_Hplus_fluxes(n_time,*)), $
              psym=6,color=100+n_time
        oplot,MagEIS_Hplus_Energy(*,n_time), $
              alog10(midE_Hplus_fluxes(n_time,*)), $
              psym=5,color=100+n_time
        
     endfor
     
     plot,HOPE_e_Energy(*,0), $
          alog10(lowE_e_fluxes(0,*)), $
          xr=[.001,60000000],/xlog,xtitle='Energy (keV)', $
          yr=[0,10],ytitle=TEXTOIDL('log_{10}(e^{-} Flux) ' + $
                                    '(cm^{2} s sr keV)^{-1}'), $
          psym=6,charsize=char_size,/nodata
     legend,['HOPE','MagEIS','REPT', $
             'L='+STRTRIM(STRING(data_HOPE_Oplus[0].L,FORMAT='(F4.2)'),2), $
             'L='+STRTRIM(STRING(data_HOPE_Oplus[0].L,FORMAT='(F4.2)'),2), $
             'L='+STRTRIM(STRING(data_HOPE_Oplus[0].L,FORMAT='(F4.2)'),2)], $
            psym=[6,5,4,0,0,0],color=[100,100,100,100,101,102], $
            pos=legend_pos_arr
     for $
        n_time=0,n_time_max $
     do begin
        
        oplot,HOPE_e_Energy(*,n_time), $
              alog10(lowE_e_fluxes(n_time,*)), $
              psym=6,color=100+n_time
        oplot,MagEIS_e_Energy(*,n_time), $
              alog10(midE_e_fluxes(n_time,*)), $
              psym=5,color=100+n_time
        oplot,REPT_e_Energy_arr,alog10(highE_e_fluxes(n_time,*)), $
              psym=4,color=100+n_time
        
     endfor
     
     device,/close

     device,filename='RBSP_Oplus_fin_Map.eps', $
            xsize=6.8,xoffset=1., $
            ysize=6.8,yoffset=1., $
            /inches,/encapsulated,/portrait,/color,bits_per_pixel=8
     
     loadct,3,/silent
     w,1,1
     image_cont,alog10(average_fluxes_HOPE_Oplus), $
                L_fin_output,HOPE_Oplus_Energy,minval=-11, $
                xtitle=TEXTOIDL('L (R_E)'), $
                ytitle=TEXTOIDL('Energy (keV)')

     device,/close
     set_plot,'x'
     
  endif

  return

end
