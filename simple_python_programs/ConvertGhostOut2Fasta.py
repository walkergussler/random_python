

#files=["ghost_01_fix.fasta","ghost_02_1a_fix.fasta","ghost_02_1b_fix.fasta","ghost_03_fix.fasta","ghost_04_fix.fasta","ghost_05_1a_fix.fasta","ghost_05_1b_fix.fasta"]
files=["ghost_07.fasta","ghost_08.fasta"]

for file in files:
    num_lines = sum(1 for line in open(file))
    print(file,num_lines)
    #for my in mine:
    f=open(file,"r")
    newfile=file[:len(file)-6]+"_fix"+file[-6:]
##    #newfile=file[:-10]+file[-6:]
    g=open(newfile,"w")
    lines=f.readlines()
    f.close()
    for line in lines:
        if line[0]==">":
            g.write(line)
            l=[]
        else:
            if len(line)==61:
                asdf=line[:-1]
                l.append(asdf)
            else:
                l.append(line)
                seq="".join(l)
                g.write(seq)
##
###prog()
