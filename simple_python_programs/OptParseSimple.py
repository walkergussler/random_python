#!/usr/bin/python
import sys
import os
import argparse

parser = argparse.ArgumentParser(description='Get majority sequence from raw reads files.')#, epilog="Now we can do options~!~!~!")#,add_help=False)
parser.add_argument('-m','--midlen', type=int,required=True, help="The length of the MID to be trimmed")
parser.add_argument('-p1','--R1primer', type=int,required=True, help="R1 file primer length to be trimmed")
parser.add_argument('-p2','--R2primer', type=int,required=True, help="R2 file primer length to be trimmed")
parser.add_argument('-c','--cleanupFiles', action='store_true',default=False, help="Add this flag to remove intermediate files. Default is for intermediate files to not be deleted")
parser.add_argument('-f1','--file1',type=str,required=True, help="The name of the R1 file to be analyzed")
parser.add_argument('-f2','--file2',type=str,required=True, help="The name of the R2 file to be analyzed")
parser.add_argument('-ec','--errorcorrector',type=str,default="lighter",choices=['lighter','bayeshammer','karect','coral'],help='Type of error corrector to use')

args = parser.parse_args()

mid=args.midlen
clean=args.cleanupFiles
p1=args.R1primer
p2=args.R2primer
f1=args.file1
f2=args.file2
ec=args.errorcorrector
print("midlen: ",mid)
print("cleanupfiles: ",clean)
print("r1primerlen: ",p1)
print("r2primerlen: ",p2)
print("file1: ",f1)
print("file2: ",f2)
print("ec: ",ec)

#lines=f2.readlines()
#for line in lines:
#    print(line)
