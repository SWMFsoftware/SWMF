
;function read_density, file

  device,true=24
  device,de=0
  device,retain=2
  loadct,41

;file='02doy107.nc'
;file = '/home/jack/huba/ionosphere/sami3tgcm/day91_stand_alone/091a_new.nc'

;file = '/home/jack/huba/ionosphere/sami3tgcm/v1.82/ben_foster/timegcm1.42.scntr_eqnx_smin_001.nc'

file = '/home/jack/huba/ionosphere/sami3tgcm/v1.82/ben_foster/two/sami3_r694_export/timegcm1.42.scntr_eqnx_smin_hilatoff_001.nc'

ncid = ncdf_open(file)

for i=0,15 do begin
   case i of
      0: name = 'lon'
      1: name = 'lat'
      2: name = 'lev'
      3: name = 'year'
      4: name = 'day'
      5: name = 'ut'
      6: name = 'TN'
      7: name = 'O1'
      8: name = 'O2'
      9: name = 'Z'
     10: name = 'NO'
     11: name = 'H'
     12: name = 'N4S'
     13: name = 'UN'
     14: name = 'VN'
     15: name = 'POTEN'
    endcase
   
   vid = ncdf_varid(ncid, name)
   ncdf_varget, ncid, vid, value
   
   case	i of
      0: lons = value
      1: lats = value
      2: levs = value
      3: year = value
      4: day  = value
      5: ut   = value
      6: tn_data = value
      7: o1_data = value
      8: o2_data = value
      9: z_data  = value/1d5
     10: no_data = value
     11: h_data = value
     12: n4s_data = value
     13: un_data  = value/1d2
     14: vn_data  = value/1d2
     15: phi_tgcm = value
   endcase
endfor

n_lons = n_elements(lons)
n_lats = n_elements(lats)
n_levs = n_elements(levs)
n_hours = n_elements(ut)

;tn_data = correct_data(tn_data, n_hours, n_levs)

;o1_data = correct_data(o1_data, n_hours, n_levs, /NON_ION_SPECIES)

;o2_data = correct_data(o2_data, n_hours, n_levs, /O_2, /NON_ION_SPECIES)

;no_data = correct_data(no_data, n_hours, n_levs, /NON_ION_SPECIES)

;h_data = correct_data(h_data, n_hours, n_levs, /NON_ION_SPECIES)

;n4s_data = correct_data(n4s_data, n_hours, n_levs, /NON_ION_SPECIES)

;; calc_n2
n2_data = 1.0 - o2_data - o1_data

;; calc_pressure
p_naught = 5.0e-4
pressure_data = transpose(cmreplicate(p_naught * exp(-1.0 * levs) / 10.0,[n_lons,n_lats,n_hours]),[1,2,0,3])


;; calc_tnd
boltz = 1.3807e-23
tnd_data = (pressure_data / boltz / tn_data) / 1.0e6


;; calc_mmm
o1_amu  = 16.0
o2_amu  = 32.0
n2_amu  = 28.0
no_amu  = 30.0
h_amu   =  1.0
n4s_amu = 14.0
mmm_data = 1.0 / ( (o2_data / o2_amu) + (o1_data / o1_amu) + (n2_data / n2_amu) )


;; calc_rho
rho_data = tnd_data * mmm_data / 6.02e23

;; calc o1_density
o1_density = o1_data * tnd_data * mmm_data / o1_amu

;; calc o2_density
o2_density = o2_data * tnd_data * mmm_data / o2_amu

;; calc n2_density
n2_density = n2_data * tnd_data * mmm_data / n2_amu

;; calc no_density
no_density = no_data * tnd_data * mmm_data / no_amu

;; calc h_density
h_density = h_data * tnd_data * mmm_data / h_amu

;; calc n_density
n_density = n4s_data * tnd_data * mmm_data / n4s_amu

thesize = size(o1_density)
JLON    = thesize(1)
JLAT    = thesize(2)
JALT    = thesize(3)
JTIME   = thesize(4)

print,JLON,JLAT,JALT,JTIME

end
