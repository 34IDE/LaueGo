#!/bin/bash
clear
echo " "
echo " "
echo " "
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "****"
# echo DIR = $DIR
cd "$DIR"
echo "pwd"
pwd
echo "****"

find . -name "*.xtal" -type f -exec xmllint -noout {} \;
find . -name "*.xml" -type f -exec xmllint -noout {} \;
echo " "
diff -x *.ipf -x *.pxp -x checkAll -x ".DS_Store" --brief -r . /Users/tischler/Documents/materials/standard\ xtals
