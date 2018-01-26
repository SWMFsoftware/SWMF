#!/bin/csh

source /home3/vtenishe/.cshrc

set WorkDir = /nobackup/`whoami`
cd $WorkDir/Tmp_AMPS_test

rm -rf qsub.log

foreach job (test_amps.pleiades.all.*job) #
  rm -rf AmpsTestComplete
  @ iTest=0

  while ($iTest<2)
    echo Submitting $job on `date`  >> qsub.log
    /PBS/bin/qsub $job >>& qsub.log

    echo $job is submitted
    set JobCompletedFlag='false'

    while ($JobCompletedFlag == 'false')
      sleep 60

      k=`/PBS/bin/qstat | grep AMPS_pfe`

      echo k=$k on `date`  >> qsub.log

      if ("$k" == "") then
        set JobCompletedFlag='true'

        echo $k >> qsub.log
        echo `/PBS/bin/qstat | grep vtenishe` >> qsub.log
        echo $job is completed on `date`  >> qsub.log
      endif
    end

    echo $job is completed
    sleep 180

    if (-e AmpsTestComplete) then
      break
    endif

    @ iTest++
  end
end

echo The test routine is completed on `date` 

