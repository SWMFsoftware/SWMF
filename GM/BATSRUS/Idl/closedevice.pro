;^CFG COPYRIGHT UM

pro closedevice

  if !d.name ne 'PS' then return
  device, /close
  set_plot, 'X'
  !p.font=-1

end
