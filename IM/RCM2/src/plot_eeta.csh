#!/usr/bin/csh
set kc = $argv[1]
set rec_num = 1
while ($rec_num <= 25)
      ./plot_eeta.out $rec_num $kc $kc .false. gif 1 1 1
    @ rec_num ++
end
exit
