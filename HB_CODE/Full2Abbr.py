#! /home/zxwu/SOFTWARE/anaconda3/bin/python3
# batch replace of text is needed 
# while shell can not satisfy this requirement
# a package named pandas is used here
# 2020/05/27 8:53

import pandas as pd
import os 
import sys

def helpout():
	print("no help for you")

def read_row(file_name,row_num):
	if not os.path.exists(file_name):
		return 0
	list1 = []
	with open(file_name) as f:
		line = f.readline()
		while line:
			a = line.split()
			b = a[row_num-1]
			list1.append(b)
			line = f.readline()
	return list1

if len(sys.argv) != 2:
	print("only 1 argument")
	sys.exit(2) 

if not os.path.exists("Index.txt"):
	print("Index.txt should be prepared")
	sys.exit(2)

if not os.path.exists(sys.argv[1]):
	print(sys.argv[1] + "should exist")
	sys.exit(2)

FullName = read_row('Index.txt',1)
AbbrName = read_row('Index.txt',2)
print(FullName)
print(AbbrName)
with open(sys.argv[1],'r') as f:
	TreatData = pd.read_csv(f,sep='\t',header=None)
	TreatData.replace(FullName,AbbrName,inplace=True)
with open(sys.argv[1],'w') as f:
	TreatData.to_csv(f,sep=' ',header=False,index=False)
