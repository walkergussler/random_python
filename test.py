def prealign(f1,f2):
    verbose=True
    uid1=re.findall('([^_]*)_.*',f1)[0]
    uid2=re.findall('([^_]*)_.*',f2)[0]
    with NamedTemporaryFile(delete=False) as tmpcat:
        catname=tmpcat.name
        try:
            subprocess.check_call(['cat',f1,f2],stdout=tmpcat)
        except:
            print("Invalid file. Check your -s argument if in csv mode")
            return 0
    with NamedTemporaryFile(delete=False) as aligned:
        alignname=aligned.name
        subprocess.check_call(['mafft', '--quiet', '--auto', '--thread', '20', '--preservecase', catname], stdout=aligned) 
    os.unlink(catname)
    seqs=SeqIO.parse(alignname,'fasta')
    current=uid1
    fnameIt=0
    seenSecond=False
    seqs1=[]
    while not seenSecond:
        record=next(seqs)
        thisseq=re.findall('([^_]*)_.*',record.id)[0]
        if thisseq!=uid1:
            seenSecond=True
        else:
            seqs1.append(record)
    with NamedTemporaryFile(delete=False) as align1:
        SeqIO.write(seqs1,align1,'fasta')
        f1name=align1.name
    seqs2=[]
    seqs2.append(record)#write record already accessed and ignored because it doesn't belong in first file
    co2=1
    done=False
    while not done:     
        try:
            record=next(seqs)
            seqs2.append(record)
        except StopIteration:
            done=True
    with NamedTemporaryFile(delete=False) as align2:
        SeqIO.write(seqs2,align2,'fasta')
        f2name=align2.name
    os.unlink(alignname)    
    if verbose:
        print("Running on: %s" %", ".join([f1,f2]))
    return(f1name,f2name)

from Bio import SeqIO
import os
import re
import subprocess
from tempfile import NamedTemporaryFile

f1name,f2name=prealign('AA20_unique_1b_43.fas','AA45_unique_1b_161.fas')
# f1name="align1.fas"
# f2name="align2.fas"
seqs1=SeqIO.parse(f1name,'fasta')
seqs2=SeqIO.parse(f2name,'fasta')
for seq in seqs1:
    print(seq.seq)
for seq in seqs2:
    print(seq.seq)