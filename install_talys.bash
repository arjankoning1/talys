#!/bin/bash
unalias -a
echo
echo "       TALYS installation (Version December 29 2024) (C) Copyright 2024 Arjan Koning All Rights Reserved"
echo
echo " Two ways to use this script:"
echo
echo " install_talys.bash 'Arjan Koning'"
echo "     Replace  Arjan Koning by your own name"
echo
echo " or"
echo
echo " install_talys.bash"
echo "     after which you will be prompted to input your name"
echo
if [ $# -eq 1 ] ; then
  yourname=$1
else
  echo 'Enter your name (which will appear in the output files): '
  read yourname
fi
echo ${yourname}
pfile=path_change.bash
if [ -e $pfile ] ; then
  sed "s/user=.*/user='${yourname}'/" $pfile > tmp
  mv -f tmp $pfile
fi
chmod a+x $pfile
code_build.bash talys
