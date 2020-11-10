import sys
with open(sys.argv[1],"r") as f:
    lines=f.readlines()
with open("asdf.py","w") as f:
    for line in lines:  
        numSpaces=0
        for char in line:
            if char!=" ":
                break    
            numSpaces+=2
        if len(line)>1:
            l=" "*numSpaces+line.lstrip()
            f.write(l)
        else:
            f.write(line)


