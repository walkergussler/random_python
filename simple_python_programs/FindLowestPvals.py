f=open("asdf.txt","r")
lines=f.readlines()
f.close()
l=[]
min=1
for line in lines:
    splitline=line.split("\t")
    if float(splitline[3])<min and float(splitline[3])!=0:
        print(splitline[3])
        min=float(splitline[3])
