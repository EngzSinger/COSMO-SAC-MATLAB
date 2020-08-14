#! /bin/sh
# set s serise of chb and alpha initial value 
# sub matlab job
# 2020/05/19 21:34 by zxwu

alpha=(5000 7000 9000 11000)
chb=(2000 3000 4000)

for a in ${alpha[@]}
do
	for b in ${chb[@]}
	do
		if [ -d ${a}_${b} ]
		then 
			rm -r ${a}_${b}
		fi
		mkdir ${a}_${b}
		cp *m ${a}_${b}
		cp *txt ${a}_${b}
		cp matlab.pbs ${a}_${b}
		cd ${a}_${b}
		sed -i "9c x0=[$a,$b];" CosmoParaOptim.m
		sed -i "3c\#PBS -l nodes=1:ppn=24" matlab.pbs
		qsub matlab.pbs
		cd ..
	done
done
