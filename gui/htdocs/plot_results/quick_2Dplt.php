<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php
  Exec("rsync -av $imagedir/style.sty $imagedir/script.mcr $imagedir/tecplot.mcr $batchdir/$tmpdir/;
        rsync -av $imagedir/batch-${number}.mcr $batchdir/$tmpdir/batch.mcr");
  Exec("cd $batchdir/$tmpdir;
        echo '#!/bin/sh' > batchscript.sh;
        echo '' >> batchscript.sh;
        echo '${tecplot} -p batch.mcr -b ../../../${cmp}/${plotfile}' >> batchscript.sh;
        echo 'cp -f batch.log ../.' >> batchscript.sh;
        echo 'mv print.cps ${file1}' >> batchscript.sh;
        echo '' >> batchscript.sh;
        chmod 755 batchscript.sh");
?>
