#!/bin/bash
if [[ -z $1 ]]; then
	printf '\e[8;25;120t'
	clear
fi

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
echo "**********  Not doing a diff, this is the main directory for xtal files.  **********"
echo " "
echo " "
diff -x *.ipf -x *.pxp -x "checkAll.command" -x ".DS_Store" --brief -r . /Users/tischler/Documents/materials/standard\ xtals
