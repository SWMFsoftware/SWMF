#!/bin/sh

IDL_PATH=Idl:${IDL_PATH}
IDL_STARTUP=Idl/idlrc_gui
export IDL_PATH IDL_STARTUP

echo show_head,"'"$1"'" | idl
