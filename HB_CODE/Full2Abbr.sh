#! bin/bash 
# for all compounds in the database,they all have both full name and abbravation
# full name is for accurate recongnize the specise while abbrevation is for results summary
# in this context, a script for intertransting both is demanded urgently
# This script reads the file named Index.txt under this dir to build a connection between full name and abbrevation
# then replace the full name or abbrevation ever seen in the file metioned by the  argument with its counterpart
# Furthermore, this script offers a method for replace text by batch operation

# 2020/05/26 11:57 by zxwu

helpout() {
cat <<EOF
			Batch_Replace

	GENERAL SHELL-SCRIPT FOR REPLACE VARITY TEXT 

[ -h ] or [ -help ] view this help

The script asks for one argument for normal excute
Index.txt containing the relationship between the text to be dealt and counterpart after treating should be prepared.
The second one is the file to be operated.
EOF
}

###########################################################
# Script STARTS HERE
###########################################################

if [ -n $1 ]; then
	case $1 in 
		"-h" |  "-?" | "-help" | "--help") helpout ; exit 0 ;;
	esac
fi
if [ $? -ne 1 ]; then
	helpout
	exit 0
fi

if [ ! -f Index.txt ]; then
cat <<EOF
	Index.txt should be prepared !!!
EOF
	exit 0
fi

TXT_FLAG=`cat Index.txt | awk 'NF!=2'`
if [ -n "$TXT_FLAG" ]; then
	cat <<EOF
	Index.txt should have only two item each lines
EOF
	exit 0
fi

if [ ! -f $1 ]; then
cat <<EOF
	Please be sure that file $1 exist
EOF
	exit 0
fi

for a in $(cat $1)
do
	c=`cat Index.txt | awk -v b=$a '$1==b {print $2}'`
	
