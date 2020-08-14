#! /bin/bash
# 2020/05/20 10:47 by zxwu

IFS=$'\n'
grep -r 'Error using parpool' > .tmp
sed -i '/reRun.sh/d' .tmp
for a in $(cat .tmp)
do
	echo $a > .tmpl
	O_FILE=`cat .tmpl | awk -F: '{print $1}'`
	FILE_PATH=`dirname $O_FILE`
	rm $O_FILE
	echo $FILE_PATH
	cd $FILE_PATH 
	qsub matlab.pbs
	cd ..
	rm .tmpl
done
rm .tmp
