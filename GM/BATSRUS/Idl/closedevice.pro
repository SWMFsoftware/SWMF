;^CFG COPYRIGHT UM

pro closedevice

  if !d.name eq 'PS' then begin
    device, /close
    set_plot, 'X'
    !p.font=-1
  endif

  return

end
