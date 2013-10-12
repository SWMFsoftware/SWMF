#!/bin/sh
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

IDL_PATH=Idl:${IDL_PATH}
IDL_STARTUP=Idl/idlrc_gui
export IDL_PATH IDL_STARTUP

echo show_head,"'"$1"'" | idl
