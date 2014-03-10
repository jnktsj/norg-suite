#!/bin/bash

# Copyright 2012, 2013, 2014 Junko Tsuji

#============================================================
# DESCRIPTION:
#   This script downloads software packages required for
#   norg-seq and norg-age runs. Installed tools are:
#
#     * norg-seq:
#       - LAST  : http://last.cbrc.jp/
#       - tantan: http://www.cbrc.jp/tantan
#
#     * norg-age:
#       - MAFFT : http://mafft.cbrc.jp/alignment/software/
#       - RAxML : http://sco.h-its.org/exelixis/software.html
#
# USAGE: fetch.sh [options]
#
#     -h: show options
#     -a: fetch all software packages
#     -l: fetch LAST
#     -t: fetch tantan
#     -m: fetch MAFFT
#     -r: fetch RAxML
#============================================================


######## START: parameter ########

# flags
tantan="FALSE"
last="FALSE"
mafft="FALSE"
raxml="FALSE"

# software directories
tantan_dir=""
last_dir=""
mafft_dir=""
raxml_dir=""

########  END: parameter  ########


######## START: functions ########

function help {
  echo ""
  echo "This script downloads software packages required for"
  echo "norg-seq and norg-age runs. Installed tools are:"
  echo ""
  echo "  * norg-seq:"
  echo "    - LAST  : http://last.cbrc.jp/"
  echo "    - tantan: http://www.cbrc.jp/tantan/"
  echo ""
  echo "  * norg-age:"
  echo "    - MAFFT : http://mafft.cbrc.jp/alignment/software/"
  echo "    - RAxML : http://sco.h-its.org/exelixis/software.html"
  usage
}

function usage {
  echo ""
  echo "Usage: fetch.sh [options]"
  echo ""
  echo "    -h: show options"
  echo "    -a: fetch all softwares"
  echo "    -l: fetch LAST"
  echo "    -t: fetch tantan"
  echo "    -m: fetch MAFFT"
  echo "    -r: fetch RAxML"
  echo ""
}

########  END:  functions ########


if [ $# -eq 0 ]
then
  usage
  exit 1;
fi

# Check parameters
while getopts "haltmr" OPT
do
  case $OPT in
    "a") FLG_A="TRUE" ;;
    "t") FLG_T="TRUE" ;;
    "l") FLG_L="TRUE" ;;
    "m") FLG_M="TRUE" ;;
    "r") FLG_R="TRUE" ;;
    "h") FLG_H="TRUE"
         help
         exit 1;;
     * ) usage
         exit 1;;
  esac
done

echo ":: Set parameters ..."
if [ $FLG_A = "TRUE" ]
then
  tantan="TRUE"
  last="TRUE"
  mafft="TRUE"
  raxml="TRUE"
fi
if [ $FLG_T = "TRUE" ]
then
  tantan="TRUE"
fi
if [ $FLG_L = "TRUE" ]
then
  last="TRUE"
fi
if [ $FLG_M = "TRUE" ]
then
  mafft="TRUE"
fi
if [ $FLG_R = "TRUE" ]
then
  raxml="TRUE"
fi


#### Download software packages (from "zip" to "bin") ####
echo ":: Download software packages"
if [ $tantan = "TRUE" ]
then
  wget -N http://www.cbrc.jp/tantan/
  latestVersionNum=$(grep '<a href="tantan-[0-9]*.zip"' index.html |\
    cut -d "-" -f2 | cut -d "." -f1)
  tantan_dir=tantan-$latestVersionNum
  wget http://www.cbrc.jp/tantan/$tantan_dir.zip
  rm index.html
fi

if [ $last = "TRUE" ]
then
  wget -N http://last.cbrc.jp/archive/
  latestVersionNum=$(grep "last-[0-9]*.zip" index.html |\
    cut -d "-" -f3 | cut -d "." -f1 | sort -k1,1n | tail -n1)
  last_dir=last-$latestVersionNum
  wget http://last.cbrc.jp/archive/$last_dir.zip
  rm index.html
fi

if [ $mafft = "TRUE" ]
then
  wget -N http://mafft.cbrc.jp/alignment/software/source.html
  latestVersionNum=$(grep -i "without-extension" source.html |\
    cut -d '"' -f2 | awk 'BEGIN{ FS="-src" } { print $1 }')
  mafft_dir=$latestVersionNum
  wget http://mafft.cbrc.jp/alignment/software/$mafft_dir-src.tgz
  rm source.html
fi

if [ $raxml = "TRUE" ]
then
  wget -N http://sco.h-its.org/exelixis/software.html
  latestVersionNum=$(grep -i "countSource" software.html |\
    head -n1 | cut -d '"' -f 2)
  wget http://sco.h-its.org/exelixis/$latestVersionNum
  rm software.html
fi

echo ":: Completed"

