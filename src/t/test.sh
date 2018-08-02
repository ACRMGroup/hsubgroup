#!/bin/bash

#set -o errexit
set -o nounset
set -o pipefail

rm -f ./test.out

../hsubgroup ./test.pir > test.out

diff -w test.out.compare test.out

if [ $? -ne 0 ]; then
   echo "hsubgroup: unexpected output!";
   exit 1
else
   echo "hsubgroup: test passed";
fi
