<?php
  $runidlexists = "0";
  $dir2 = opendir($imagedir);
  while( $file = readdir( $dir2 ) ) {
    if (eregi("runidl.sh", $file)) {
      $runidlexists = "1";
    }
  }
  if(! $runidlexists) {
    Exec("cd $imagedir;
          echo '#!/bin/sh' > runidl.sh;
          echo '' >> runidl.sh;
          echo 'IDL_PATH=../../../Idl:\${IDL_PATH}' >> runidl.sh;
          echo 'IDL_STARTUP=../../../Idl/idlrc_gui' >> runidl.sh;
          echo 'export IDL_PATH IDL_STARTUP' >> runidl.sh;
          echo '' >> runidl.sh;
          echo 'idl batch$macroextension' >> runidl.sh;
          chmod 755 runidl.sh");
  }
  Exec("rsync -av $imagedir/runidl.sh $batchdir/$tmpdir/;
        cd $batchdir/$tmpdir;
        ln -s $cwd/$filedir/$plotfile file.out;
        rsync -av $imagedir/batch-001.pro $batchdir/$tmpdir/batch.pro");
  Exec("cd $batchdir/$tmpdir;
        echo '#!/bin/sh' > batchscript.sh;
        echo '' >> batchscript.sh;
        echo './runidl.sh >& batch.log' >> batchscript.sh;
        echo 'cp -f batch.log ../.' >> batchscript.sh;
        echo 'mv *ps $file1' >> batchscript.sh;
        chmod 755 batchscript.sh");
?>
