#!/bin/bash
#
# FixCVSRoot.sh: This script will replace all occurrences of /p/volov/DEVTREE
#   with /FRAMEWORK in all files found below the current directory.

echo Replacing /p/volov/DEVTREE with /FRAMEWORK in $PWD

FILES=`find . -type f -name Root`

for file in $FILES
do
  if grep -q '/p/volov/DEVTREE' $file
  then
    echo Fixing $file
    cat $file | sed s/\\/p\\/volov\\/DEVTREE/\\/FRAMEWORK/ > fixedfile
    mv -f fixedfile $file
  fi
done

FILES=`find . -type f -name Repository`

for file in $FILES
do
  if grep -q '/p/volov/DEVTREE' $file
  then
    echo Fixing $file
    cat $file | sed s/\\/p\\/volov\\/DEVTREE/\\/FRAMEWORK/ > fixedfile
    mv -f fixedfile $file
  fi
done

exit 0
