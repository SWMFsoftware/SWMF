

; READ_INDICES_DIR.PRO

; set the color table 
  
  device,true=24
  device,de=0
  device,retain=2
  loadct,41

; sami3tgcm

  nz    = 160
  nion  = 7
  nf    = 50
  nl    = 96
  nt    = 11

; kate zawdie

;  nz    = 201
;  nion  = 7
;  nf    = 200
;  nl    = 96
;  nt    = 14

; esf

;  nz = 101
;  nf = 202
;  nion = 7
;  nl = 96
;  nt = 101

; volland-stern

;  nz  = 160
;  nf  = 200
;  nl  =  90
;  nt  = 187


  nx = 100  ; latitude
  ny = 100  ; altitude

  nnx = nl+1
  nny = nf-1

;  dir='/media/disk-5/huba/ionosphere/sami3_mpi-1.80cr_theta/volland_stern_5cr_mjm/'
;  dir='/media/disk-5/huba/ionosphere/sami3_mpi-1.80_theta/volland_stern_2_mjm/'
;  dir = '/home/jack/huba/ionosphere/sami3/sami3_mpi-1.84_p/euvac_hwm07/run0/'

;  dir = '/home/jack/huba/ionosphere/sami3/sami3_mpi-1.84_p/euvac_hwm07/primo_run4_no_ntphotoion/'
;  dir = '/home/jack/huba/ionosphere/sami3/sami3_mpi-1.84_p/euvac_hwm07/primo_run4_no_ntphotoion_tvn1.5/'

;  dir='/home/jack/huba/tmp/zawdie/run_sami3/'

dir='/home/jack/huba/ionosphere/sami3/sami3_mpi-1.84_p_fac/euvac_hwm07/code/'

  end






