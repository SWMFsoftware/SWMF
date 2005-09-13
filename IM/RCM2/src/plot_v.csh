#!/usr/bin/csh
set rec_num = 1
while ($rec_num <= 25)
      ./plot_v.out 2 $rec_num gif 1 2 'RCM-BATSRUS TEST RUN 1'
    @ rec_num ++
end
exit
