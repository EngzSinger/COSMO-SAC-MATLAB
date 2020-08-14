#! /bin/bash 

a=`cat Index.txt | awk 'NF!=2'`
if [ -n $a ]; then echo ringt; else echo no; fi
echo $a
if [ $a ]
then
	echo right
else
	echo no
fi

b="abc"
c="abc"
if [ $b = $c ]
then
	echo 1
fi
if [ "$b" = "$c" ]
then 
	echo 2
fi
echo $b $c
echo "$b" "$c"
if [ -z "$b" ]
then 
	echo non
else 
	echo o
fi

