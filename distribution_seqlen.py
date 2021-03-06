from Bio import SeqIO
import os,sys
from collections import defaultdict

def see(file):
    print(file)
    lens=defaultdict(int)
    with open(file) as f:
        for record in SeqIO.parse(f,'fasta'):
            lens[len(record.seq)]+=1
            # if len(record.seq)!=264:
                # print(record.id)
                # print(record.seq)
                # input()
    for item in lens:
        print(item,lens[item])
    # print('')


    
if __name__=="__main__":
    s=len(sys.argv)
    if s==1:
        for file in os.listdir(os.getcwd()):
            if file.endswith('fas'):
            # if 'NH' in file and '264' in file:
                see(file)
    elif s==2:
        print(see(sys.argv[1]))
        


#################TO ALIGN IF THEY ARENT##########################
#        if not see(file):
#            print('mafft --quiet --auto --thread 20 --preservecase '+file+' > aligned/'+file)
#            os.system('mafft --quiet --auto --thread 20 --preservecase '+file+' > aligned/'+file)
#        else:
#            print('cp '+file+' aligned/'+file)
#            os.system('cp '+file+' aligned/'+file)
