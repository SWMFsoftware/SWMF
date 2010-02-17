#!/bin/csh -f
setenv ssw_dir /csem1/ssw
setenv SSW ${ssw_dir}
setenv SSW_INSTR "eit sxt chianti"
source $SSW/gen/setup/setup.ssw
${SSW}/gen/setup/ssw_idl call_euv.pro
exit
