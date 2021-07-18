#!/bin/sh

# this script will list the modules not appearing in
# any testscripts in $AAHOME/developer/testscripts

# we assume the script lives in $AAHOME/developer
# (we need this to find the aa_modules directory!) 

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

AAHOME=`dirname $SCRIPTPATH`

MCOUNT=0
TCOUNT=0

for FNAME in $AAHOME/aa_modules/*.xml
do
   let MCOUNT=MCOUNT+1
   FILENAME=`basename $FNAME`
   MODULENAME="${FILENAME%%.*}"
   if ! grep -q $MODULENAME $SCRIPTPATH/testscripts/*.xml; then
     echo $MODULENAME not tested
   else
     let TCOUNT=TCOUNT+1
   fi
done

echo $TCOUNT out of $MCOUNT modules tested
