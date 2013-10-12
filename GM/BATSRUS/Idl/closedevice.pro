;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
;^CFG COPYRIGHT UM

pro closedevice

  if !d.name ne 'PS' then return
  device, /close
  set_plot, 'X'
  !p.font=-1

end
