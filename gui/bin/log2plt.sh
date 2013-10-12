#!/bin/sh
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

if [ -n "$1" ]
then
  echo ""
  echo "This script converts the standard GM log output into tecplot format"
  echo "  and preplots the file if possible."
  echo "It looks for files matching */log*log or */sat*sat and takes no arguements."
  echo ""
  echo "  Usage: `basename $0`"
  echo ""
  exit 0
fi

for FILEIN in `ls */log*log`
do
  # Set filenames
  FILEDAT=`dirname $FILEIN`/`basename $FILEIN .log`.dat
  FILEPLT=`dirname $FILEIN`/`basename $FILEIN .log`.plt
  FILELOG=`dirname $FILEIN`/process.log

  # Echo message
  MSG="Processing file $FILEIN into file $FILEDAT ..."
  echo ""
  echo $MSG
  echo $MSG > $FILELOG

  # Parse header
  STR=`head -1 $FILEIN`
  echo "TITLE = \"$STR\"" > $FILEDAT
  STR=`head -2 $FILEIN | tail -1`
  echo "VARIABLES = " >> $FILEDAT
  for VAR in ${STR[*]}
  do
    echo " \"$VAR\"" >> $FILEDAT
  done

  # Get number of lines in file
  WC=`wc $FILEIN`
  for VALUES in ${WC[*]}
  do
    LINES=$VALUES
    break
  done
  let "LINES-=2"

  # Write zone and data
  echo "ZONE T=\"logfile\"" >> $FILEDAT
  echo " I=$LINES, J=1, K=1, ZONETYPE=Ordered, DATAPACKING=POINT" >> $FILEDAT
  tail -$LINES $FILEIN >> $FILEDAT

  # Preplot
  PP=`which preplot | grep -v 'not found'`
  if [ $PP ]
  then
    preplot $FILEDAT >> $FILELOG
    if [ -e $FILEPLT ]
    then
      rm -f $FILEDAT $FILELOG
    else
      mv -f $FILELOG ${FILELOG}_keep
      echo "  WARNING: preplot failed, not removing .dat file."
    fi
  else
    echo "  WARNING: preplot is not in your path, skipping."
  fi
done

for FILEIN in `ls */sat*sat`
do
  # Set filenames
  FILEDAT=`dirname $FILEIN`/`basename $FILEIN .sat`.dat
  FILEPLT=`dirname $FILEIN`/`basename $FILEIN .sat`.plt
  FILELOG=`dirname $FILEIN`/process.log

  # Echo message
  MSG="Processing file $FILEIN into file $FILEDAT ..."
  echo ""
  echo $MSG
  echo $MSG > $FILELOG

  # Parse header
  STR=`head -1 $FILEIN`
  echo "TITLE = \"$STR\"" > $FILEDAT
  STR=`head -2 $FILEIN | tail -1`
  echo "VARIABLES = " >> $FILEDAT
  for VAR in ${STR[*]}
  do
    echo " \"$VAR\"" >> $FILEDAT
  done

  # Get number of lines in file
  WC=`wc $FILEIN`
  for VALUES in ${WC[*]}
  do
    LINES=$VALUES
    break
  done
  let "LINES-=2"

  if [ "$LINES" -gt "1" ]
  then
    # Write zone and data
    echo "ZONE T=\"satfile\"" >> $FILEDAT
    echo " I=$LINES, J=1, K=1, ZONETYPE=Ordered, DATAPACKING=POINT" >> $FILEDAT
    tail -$LINES $FILEIN >> $FILEDAT
  
    # Preplot
    PP=`which preplot | grep -v 'not found'`
    if [ $PP ]
    then
      preplot $FILEDAT >> $FILELOG
      if [ -e $FILEPLT ]
      then
        rm -f $FILEDAT $FILELOG
      else
        mv -f $FILELOG ${FILELOG}_keep
        echo "  WARNING: preplot failed, not removing .dat file."
      fi
    else
      echo "  WARNING: preplot is not in your path, skipping."
    fi
  else
    echo "ERROR, empty file: $FILEIN"
    rm -f $FILEDAT
  fi
done

echo "... done."
echo ""

exit 0
