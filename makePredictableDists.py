strLen=264
acount=-1
f=open("predictable.fas","w")
for seq in range(strLen):
    acount+=1
    tcount=0
    working=""
    for a in range(strLen):
        tcount+=1
        if tcount<acount:
            working+="A"
        else:
            working+="T"
    f.write(">SEQ_1\n")
    f.write(working+"\n")
