askstr,'logfilename(s) ',logfilename,doask
askstr,'logfunc(s)     ',logfunc,doask

nlog=0
str2arr,logfilename,logfilenames,nlog

plot_log, logfilename, logfunc,                          $
  wlog,wlognames, wlog1, wlognames1, wlog2, wlognames2,  $
  xrange=xrange, yranges=yranges, timeshifts=timeshifts, $
  smooths=smooths,                                       $
  colors=colors, symbols=symbols, linestyles=linestyles, $
  title=title, xtitle=xtitle, ytitles=ytitles

end
