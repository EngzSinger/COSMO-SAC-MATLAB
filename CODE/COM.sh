#! /bin/bash
# combine all cation and anion from database 
# 2020/05/21 20:10 by zxwu

if [ -f CombList.txt ]
then
	rm -f CombList.txt
fi

t=(ANION CATION SOLUTE)
for a in ${t[@]}
do
	cat /home/zxwu/WORK_SPACE/COSMO_SAC/DATABASE/COSMO/$a/Index.txt | awk '{print $2}' >$a
	dos2unix $a
done

for b in $(cat CATION)
do
	for c in $(cat ANION)
	do
		for d in $(cat SOLUTE)
		do
			echo  $d $b $c  298.15 >> CombList.txt
		done
	done
done

