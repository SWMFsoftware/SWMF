;  Copyright (C) 2002 Regents of the University of Michigan, 
;  portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf

common ask_param
common getlog_param

askstr,'logfunc(s)     ',logfunc,doask

if !d.name eq 'X' then begin
    if !d.window lt 0 then window
    wshow
endif

plot_log

end
