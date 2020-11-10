import sys
def main(input,col,delim):
    col=int(col)
    with open(input) as f:
        lines=f.readlines()
    l=[]
    for line in lines:
        s=line.split(',')
        index=col-1
        l.append(float(s[index]))
    return sum(l)/float(len(l))
    
if __name__=='__main__':
    if len(sys.argv)==1:
        print('python printAverage.py input col delim')
    if len(sys.argv)==2:
        print(main(sys.argv[1],1,','))
    if len(sys.argv)==3:
        print(main(sys.argv[1],sys.argv[2],','))
    if len(sys.argv)==4:
        print(main(sys.argv[1],sys.argv[2],sys.argv[3]))