#!/usr/bin/python
import sys
import os

os.system("rm stats.txt")
kset=["24","28","32"]
for k in kset:
	os.system("rm distfile.txt")
	old=0
	TP=0
	TN=0
	FP=0
	FN=0
	borderhit=0
	bordermiss=0
	t=.15
	b=.1
	for i in os.listdir(os.getcwd()):
		if i[-4:]!=".fas":
			continue
		wr="mash sketch -k "+k+" "+i
		os.system(wr)
	for i in os.listdir(os.getcwd()):
		if i[-4:]==".msh":
			if old==2:
				writestring="mash dist "+i+" olpasty.msh >> distfile.txt"
				os.system(writestring)
				merges="mash paste PastedMash "+i+" olpasty.msh"
				os.system(merges)
				mover="mv PastedMash.msh olpasty.msh"
				os.system(mover)
			elif old==1:
		                old=2
				writestring="mash dist "+i+" "+first+" >> distfile.txt"
				os.system(writestring)
				merges="mash paste olpasty "+i+" "+first
				os.system(merges)
			elif old==0:
				first=i
				old=1
	os.system("rm *.msh")
	f=open("distfile.txt","r")
	lines=f.readlines()
	f.close()
	lineno=0
	for line in lines:
		splitline=line.split("\t")
		d=float(splitline[2])
		if splitline[0][:2]!=splitline[1][:2]:
			if d<b:
				FP+=1
			elif d>t:
				TN+=1
			else:
				bordermiss+=1
		else:
			if d<b:
				TP+=1
			elif d>t:
				FN+=1
			else:
					borderhit+=1
	f=open("stats.txt","a")
	headerline="TP\tTN\tFP\tFN\tiffyhit\tiffymiss   for k="+k+", upperT="+str(t)+", lowerT="+str(b)+"\n"
	dataline=str(TP)+"\t"+str(TN)+"\t"+str(FP)+"\t"+str(FN)+"\t"+str(borderhit)+"\t"+str(bordermiss)+"\n"
	f.write(headerline)
	f.write(dataline)
f.close()
os.system("cat stats.txt")
