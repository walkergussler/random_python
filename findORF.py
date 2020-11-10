from Bio import SeqIO

files=['AA20_unique_1b_43.fas','AA8_unique_1b_88.fas','AA45_unique_1b_161.fas','AC122_unique_1a_48.fas','AC123_unique_1a_48.fas','AC124_unique_1a_47.fas','AC121_unique_1a_48.fas','AJ199_unique_1b_29.fas','AJ6_unique_1b_68.fas','AJ86_unique_1b_32.fas','AW28_unique_1b_48.fas','AW12_unique_1b_48.fas','AW21_unique_1b_48.fas','AW18_unique_1b_48.fas','AW10_unique_1b_48.fas','AW35_unique_1b_48.fas','AW34_unique_1b_48.fas','AW7_unique_1b_47.fas','AW22_unique_1b_47.fas','AW20_unique_1b_47.fas','AW16_unique_1b_52.fas','AW1_unique_1b_38.fas','AW9_unique_1b_48.fas','AW15_unique_1b_6.fas','AW27_unique_1b_47.fas','AW2_unique_1b_49.fas','AW4_unique_1b_47.fas','AW8_unique_1b_48.fas','AW13_unique_1b_48.fas','BA2_unique_1b_10.fas','BA3_unique_1b_52.fas','BA4_unique_1b_23.fas','BA8_unique_1b_33.fas','BA16_unique_1b_15.fas','BA1_unique_1b_9.fas','BB29_unique_1a_59.fas','BB42_unique_1a_84.fas','BB45_unique_1a_177.fas','BB41_unique_1a_30.fas','BB1_unique_1a_99.fas','BB44_unique_1a_12.fas','BB31_unique_1a_94.fas','BC46_unique_1a_69.fas','BC30_unique_1a_62.fas','BJ23_unique_1a_7.fas','BJ25_unique_1a_51.fas','BJ30_unique_1a_6.fas','BJ28_unique_1a_66.fas','AQ28_unique_1a_5.fas','AQ6_unique_1a_10.fas','AQ4_unique_1a_39.fas','AQ90_unique_1a_16.fas','AQ89_unique_1a_15.fas','AQ25_unique_1a_12.fas','AQ13_unique_1a_9.fas','AQ16_unique_1a_11.fas','AQ92_unique_1a_23.fas']
stops=['TAA','TGA','TAG']



for file in files:
    stops1=[]
    stops2=[]
    stops3=[]

    seqs=[]
    with open(files[0]) as f:
        for record in SeqIO.parse(f,'fasta'):
            seqs.append(record.seq)
            
    starts1=[]
    starts2=[]
    starts3=[]
    HAPLOSIZE=len(seqs[0])
    for seq in seqs:
        s1=0
        s2=0
        s3=0
        n=0
        while n<HAPLOSIZE-4:
            starts1.append(str(seq[n:n+3]))
            starts2.append(str(seq[n+1:n+4]))
            starts3.append(str(seq[n+2:n+5]))
            n+=3
        if n+5==HAPLOSIZE:
            starts1.append(str(seq[n:n+3]))
            starts2.append(str(seq[n+1:n+4]))
            starts3.append(str(seq[n+2:n+5]))
        elif n+4==HAPLOSIZE:
            starts1.append(str(seq[n:n+3]))
            starts2.append(str(seq[n+1:n+4]))
        elif n+3==HAPLOSIZE:
            starts1.append(str(seq[n:n+3]))
        for item in starts1:
            if item in stops:
                s1+=1
        for item in starts2:
            if item in stops:
                s2+=1
        for item in starts3:
            if item in stops:
                s3+=1
        stops1.append(s1)
        stops2.append(s2)
        stops3.append(s3)
    m1=min(stops1)
    m2=min(stops2)
    m3=min(stops3)
    if m3<m1:
        if m3<m2:
            print(file,'+2')
        else:
            print(file,'+1')
    else:
        if m1<m2:
            print(file,'+0')
        else:
            print(file,'+1')
    print(m1,m2,m3)
    
# seqs={}
# with open('aatest') as f:
    # for record in SeqIO.parse(f,'fasta'):
        # seqs[record.seq]=record.id
# with open('AA20_unique_1b_43.fas') as f:
    # for record in SeqIO.parse(f,'fasta'):
        # if record.seq in seqs:
            # print("in")
        # else:
            # print("out")
            # print(seq)