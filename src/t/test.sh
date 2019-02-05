#!/bin/bash
datafile=../../data/human.dat
#set -o errexit
set -o nounset
set -o pipefail

rm -f ./test.out

../hsubgroup ./test.pir > test.out

diff -w test.out.compare test.out

if [ $? -ne 0 ]; then
   echo "hsubgroup (no data file): unexpected output!";
   exit 1
else
   echo "hsubgroup (no data file): test passed";
fi

rm -f ./test.out

../hsubgroup -d $datafile ./test.pir > test.out

diff -w test.out.compare test.out

if [ $? -ne 0 ]; then
   echo "hsubgroup (with datafile): unexpected output!";
   exit 1
else
   echo "hsubgroup (with datafile): test passed";
fi
