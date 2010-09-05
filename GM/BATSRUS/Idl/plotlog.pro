askstr,'logfunc(s)     ',logfunc,doask

if !d.name eq 'X' then begin
    if !d.window lt 0 then window
    wshow
endif

plot_log, logfilenames, logfunc, wlog,  wlognames,         $
  wlog1, wlognames1, wlog2, wlognames2, wlog3, wlognames3, $
  wlog4, wlognames4, wlog5, wlognames5, wlog6, wlognames6, $
  wlog7, wlognames7, wlog8, wlognames8, wlog9, wlognames9, $
  xrange=xrange, yranges=yranges, timeshifts=timeshifts, $
  smooths=smooths,                                       $
  colors=colors, symbols=symbols, linestyles=linestyles, $
  title=title, xtitle=xtitle, ytitles=ytitles, timeunit=timeunit

end
