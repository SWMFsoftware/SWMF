pro closedevice

  if !d.name eq 'PS' then begin
    device, /close
    set_plot, 'X'
  endif

  return

end
