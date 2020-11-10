def prealign2(files): #this program depends on every sequence starting with the same string as the name of the file. By starting with,I mean preceeding the first underscore
    print("prealigning")
    # get names for files
    starts=[]
    id=0
    fnames=[]
    for file in files:
        id+=1
        starts.append(re.findall('([^_]*)_.*',file)[0])
        fnames.append("align"+str(id)+".fas")
    # align with os.system calls
    jstr='cat '+" ".join(files)+" > tmp.fas"
    os.system(jstr)
    mstr='mafft --quiet --auto --thread 20 --preservecase tmp.fas > align.fas'
    os.system(mstr)
    # split file into smaller files
    current=starts[0]
    fnameIt=0
    g=open(fnames[fnameIt],"w")
    with open("align.fas","rU") as f:
        for record in SeqIO.parse(f,"fasta"):
            thisseq=re.findall('([^_]*)_.*',record.id)[0]
            if thisseq!=current:
                print(thisseq)
                print(current)
                g.close()
                fnameIt+=1
                current=thisseq
                g=open(fnames[fnameIt],"w")
            SeqIO.write(record,g,"fasta")
    g.close()
    os.remove("tmp.fas")
    os.remove("align.fas")
    print("Running on: %s" %",".join(files))
    return fnames
    
def align(seqs,output):
    print("aligning")
    with NamedTemporaryFile() as f,NamedTemporaryFile() as ali_fd:
        for seq in seqs:
            f.write(">"+"_".join(seqs[seq])+"\n")
            f.write(str(seq)+"\n")
        f.flush()
        subprocess.check_call(['mafft','--quiet','--auto','--thread','20','--preservecase',f.name],stdout=ali_fd)
        ali_fd.seek(0,os.SEEK_SET)
        outs={}
        for record in SeqIO.parse(ali_fd,'fasta'):
            if len(record.seq) > 0 and len(record.seq) < 50000:
                if collections.Counter(record.seq).most_common(1)[0][0]!="-":
                    outs[record.seq]=record.id.split("_")
    return outs 