#!/usr/bin/python
import sys
import os
import re
import operator

#add things to make mid, primer length parameters, argparse, remove temp files
print("Usage: python full_assembly_pipeline_test.py <R1_file.fastq> <R2_file.fastq>")
os.system("echo  ")
os.system("echo remember to load your modules!")
os.system("echo you will need cutadapt/1.8, Python, pandaseq/2.7, Lighter, Kanalyze to run this program")
os.system("echo I only work well with paired end reads files with deep coverage of a ~300bp sequence")

def finder(file1,file2):
    r1=file1
    r2=file2
    l=[r1,r2]
    looptrack=0
    for fle in l:
        pasdf="echo Now processing "+fle
        os.system(pasdf)
        looptrack+=1
        f=open(fle,"r")
        lines=f.readlines()
        f.close()
        counter=0
        types={}
        ptypes={}
        for line in lines:
            counter+=1
            if counter%4==1:
                seqid=line
            elif counter%4==2:
                mid=line[:10]
                if looptrack==1:
                    primer=line[10:32]
                elif looptrack==2:
                    primer=line[10:30] 
                if primer not in ptypes:
                    ptypes[primer]=1
                else:
                    ptypes[primer]+=1
                if mid not in types:
                    types[mid]=1
                else:
                    types[mid]+=1
        MIDS = sorted(types.items(), key=operator.itemgetter(1))
        PRIM =sorted(ptypes.items(), key=operator.itemgetter(1))
        #print(MIDS[-1],PRIM[-1],MIDS[-2],PRIM[-2])
        MainMID=MIDS[-1][0]
        MainPRI=PRIM[-1][0]
        if PRIM[-2][1]*10>PRIM[-1][1]:
            os.system("echo WARNING: The second most common primer appears in more than 10% of the reads with your most common primer, trimming may remove a significant portion of your sequences")
        if MIDS[-2][1]*10>MIDS[-1][1]:
            os.system("echo WARNING: The second most common MID appears in more than 10% of the reads with your most common MID, trimming may remove a significant portion of your sequences")
        if looptrack==1:
            r1seq=MainMID+MainPRI
        elif looptrack==2:
            r2seq=MainMID+MainPRI
    s1="echo Sequence to be trimmed from "+r1+" : "+r1seq
    s2="echo Sequence to be trimmed from "+r2+" : "+r2seq
    os.system(s1)
    os.system(s2)
    return(r1seq, r2seq)

f1=sys.argv[1]
f2=sys.argv[2]
(pr1,pr2)=finder(f1,f2)

ca="cutadapt -g "+pr1+" -G "+pr2+" --discard-untrimmed -m 20 -o R1_cutadapt_temp.fq -p R2_cutadapt_temp.fq "+f1+" "+f2
os.system(ca)
os.system("echo ---------------------------------------------------------------------------------")
os.system("echo Now, I will correct your reads with Lighter")
os.system("echo ---------------------------------------------------------------------------------")
os.system("lighter -r R1_cutadapt_temp.fq -r R2_cutadapt_temp.fq -K 15 300 -t 16") #maybe make 'expected genome size a parameter, remove cutadapt files?
os.system("echo ---------------------------------------------------------------------------------")
os.system("echo Now, I will assemble your reads together with Pandaseq")
os.system("echo ---------------------------------------------------------------------------------")
os.system("pandaseq -F -f R1_cutadapt_temp.cor.fq -r R2_cutadapt_temp.cor.fq -u unmerged_pandaseq.fq 2> pandaseq_stats.txt 1> pandaseq_output.fastq")
os.system("java -jar /scicomp/home/mnz0/bin/kanalyze-2.0.0/kanalyze.jar count -d 16 -l 16 -t 16 -g 1000000 -c kmercount:1000 -v -k 215 -f fastq -o Kanalyze.txt pandaseq_output.fastq")
ks="sort -k 2 -n -r Kanalyze.txt > Kanalyze_sorted_"+startLetter+".txt"
os.system(ks)

sfstr="Kanalyze_sorted_"+startLetter+".txt"
sf(sfstr)

#find the longest common substring between any 2 input strings
def lcs(S,T):
    m = len(S)
    n = len(T)
    counter = [[0]*(n+1) for x in range(m+1)]
    longest = 0
    lcs_set = set()
    for i in range(m):
        for j in range(n):
            if S[i] == T[j]:
                c = counter[i][j] + 1
                counter[i+1][j+1] = c
                if c > longest:
                    lcs_set = set()
                    longest = c
                    lcs_set.add(S[i-c+1:i+1])
                elif c == longest:
                    lcs_set.add(S[i-c+1:i+1])
    return list(lcs_set)[0]

#Iterate through sorted k-mer file: determine which k-mers are the 'major' ones, and build a list of them
def sf(file):
    f=open(file,"r")
    lines=f.readlines()
    f.close()
    seqlist=[]
    numlist=[]
    lineno=0
    for line in lines:
        lineno+=1
        splitline=line.split("\t")
        num=int(splitline[1])
        if lineno!=1:
            dec=num/prevnum
            numlist.append(dec)
        prevnum=num
    val=min(numlist)
    Crit=numlist.index(min(numlist))
    cstr="echo Only looking at the first "+str(Crit)+" lines of your sorted Kmer file, there was a "+str(val)+" drop in frequency at that line"
    os.system(cstr)
    lineno=0
    for line in lines:
        lineno+=1
        splitline=line.split("\t")
        seqlist.append(splitline[0])
        if lineno+1==Crit:
            break
    seqlistcleaner(seqlist,1)

#determine distinct groups of kmers which have large portions in common
def seqlistcleaner(seqlist,it):
    badseqs=[]
    for seq1 in seqlist:
        for seq2 in seqlist:
            if len(lcs(seq1,seq2))*2<len(seq1):
                if len(lcs(seq1,seqlist[0]))*2<len(seq1):
                    seqlist.remove(seq1)
                    badseqs.append(seq1)
                elif len(lcs(seq2,seqlist[0]))*2<len(seq2):
                    seqlist.remove(seq2)
                    badseqs.append(seq2)
    contigFinder(seqlist)
    if len(badseqs)!=0:
        seqlistcleaner(badseqs,it+1)
    else:
        if it==1:
            os.system("echo I found one genotype in your file")
        else:
            itsr="echo I found "+str(it)+" different genotypes in your file"
            os.system(itsr)

#looks at set of frequent kmers which have large portions of their sequence in common, finds sequence
def contigFinder(seqs):
    seqlist=seqs
    tracker=0
    for seq in seqlist:
        if tracker==0:
            tracker+=1
            cons=seq
            continue
        temp=lcs(seq,cons)
        if len(temp)*2<len(cons):
            os.system("echo Error: new sequence found - this shouldn't have happened")
        cons=temp
    poslist=[]
    for seq in seqlist:
        pos=re.search(cons,seq).start()
        poslist.append(pos)
    first=poslist.index(max(poslist))
    last=poslist.index(min(poslist))
    num=max(poslist)*-1
    se=seqlist[first]+seqlist[last][num:]
    print(se)
