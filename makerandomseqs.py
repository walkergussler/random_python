import random
seqnum=100
seqlen=100
nucs=['A','T','C','G']

for seqid in range(seqnum):
    seq=''
    for letid in range(seqlen):
        chosen=random.choice(nucs)
        seq+=chosen
    print('>seq_21\n'+seq)
