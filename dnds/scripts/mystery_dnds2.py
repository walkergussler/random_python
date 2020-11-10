import numpy as np 
import pickle
from  itertools import permutations
import copy # @todo: remove the deepcopying?
import pdb
import math
import sys
import warnings
import time 
start_time = time.time()  # @time

def align_query_vs_reference( qry_seq, ref_seq, aln_gap_open = -10, aln_gap_extend = -0.5 ):
   from Bio import pairwise2
   alns = pairwise2.align.globalxs(qry_seq, ref_seq, aln_gap_open, aln_gap_extend)  # @todo: make sure open and extend are in right order

   top_aln = alns[0] # sorte list of alternative alignments, highest-scoring alignent is alns[0]

   qry_seq_aligned = top_aln[0]
   ref_seq_aligned = top_aln[1]
   return qry_seq_aligned, ref_seq_aligned


def trim_gaps_from_aligned_seqs( qry_seq_aln, ref_seq_aln ):
   qry_R=np.zeros(len(qry_seq_aln))

   j=0
   for i in range(len(qry_seq_aln)):
      if("-"!=qry_seq_aln[i]): 
         qry_R[i]=j
         j=j+1
   
   qry_and_ref_arr = [(j,ref_seq_aln[i],qry_R[i]) for i,j in enumerate(qry_seq_aln) if (("-" != ref_seq_aln[i]) and ("-" != j))]
   
   qry_and_ref_arr = np.array(qry_and_ref_arr)

   qry_trimmed= "".join(list(qry_and_ref_arr[:,0]))
   ref_trimmed= "".join(list(qry_and_ref_arr[:,1]))
   indeces_pro= list(qry_and_ref_arr[:,2])
   indeces_pro1= indeces_pro[0::3]
   qry_indices= ",".join(indeces_pro1)
   return qry_trimmed, ref_trimmed, qry_indices
      
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    """ Check if two values are almost equal """
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def geneticCode(name):

    """ Dictionary that maps codons to amino acids """ 

    if name == 'standard':
        gc = {  'AAA':'K', 'AAC':'N', 'AAG':'K', 'AAT':'N', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'AGA':'R', 'AGC':'S', 'AGG':'R', \
                'AGT':'S','ATA':'I','ATC':'I','ATG':'M','ATT':'I','CAA':'Q','CAC':'H','CAG':'Q','CAT':'H','CCA':'P','CCC':'P','CCG':'P', \
                'CCT':'P','CGA':'R','CGC':'R','CGG':'R','CGT':'R','CTA':'L','CTC':'L','CTG':'L','CTT':'L','GAA':'E','GAC':'D','GAG':'E', \
                'GAT':'D','GCA':'A','GCC':'A','GCG':'A','GCT':'A','GGA':'G','GGC':'G','GGG':'G','GGT':'G','GTA':'V','GTC':'V','GTG':'V', \
                'GTT':'V','TAA':'*','TAC':'Y','TAG':'*','TAT':'Y','TCA':'S','TCC':'S','TCG':'S','TCT':'S','TGA':'*','TGC':'C','TGG':'W', \
                'TGT':'C','TTA':'L','TTC':'F','TTG':'L','TTT':'F'  }
    return gc

def potential_changes_dict(nt_to_aa):

    potential_changes = {   'S': {  'AAA':0.0,'AAC':0.0,'AAG':0.0,'AAT':0.0, 'ACA':0.0, 'ACC':0.0, 'ACG':0.0, 'ACT':0.0, 'AGA':0.0, 'AGC':0.0, \
                                    'AGG':0.0, 'AGT':0.0, 'ATA':0.0, 'ATC':0.0, 'ATG':0.0, 'ATT':0.0, 'CAA':0.0, 'CAC':0.0, 'CAG':0.0, 'CAT':0.0, \
                                    'CCA':0.0,'CCC':0.0,'CCG':0.0,'CCT':0.0,'CGA':0.0,'CGC':0.0,'CGG':0.0,'CGT':0.0,'CTA':0.0,'CTC':0.0,'CTG':0.0, \
                                    'CTT':0.0,'GAA':0.0,'GAC':0.0,'GAG':0.0,'GAT':0.0,'GCA':0.0,'GCC':0.0,'GCG':0.0,'GCT':0.0,'GGA':0.0,'GGC':0.0, \
                                    'GGG':0.0,'GGT':0.0,'GTA':0.0,'GTC':0.0,'GTG':0.0,'GTT':0.0,'TAA':0.0,'TAC':0.0,'TAG':0.0,'TAT':0.0,'TCA':0.0, \
                                    'TCC':0.0,'TCG':0.0,'TCT':0.0,'TGA':0.0,'TGC':0.0,'TGG':0.0,'TGT':0.0,'TTA':0.0,'TTC':0.0,'TTG':0.0,'TTT':0.0},

                            'N': {  'AAA':0.0, 'AAC':0.0, 'AAG':0.0, 'AAT':0.0, 'ACA':0.0, 'ACC':0.0, 'ACG':0.0, 'ACT':0.0, 'AGA':0.0, 'AGC':0.0, 'AGG':0.0, \
                                    'AGT':0.0,'ATA':0.0,'ATC':0.0,'ATG':0.0,'ATT':0.0,'CAA':0.0,'CAC':0.0,'CAG':0.0,'CAT':0.0,'CCA':0.0,'CCC':0.0,'CCG':0.0, \
                                    'CCT':0.0,'CGA':0.0,'CGC':0.0,'CGG':0.0,'CGT':0.0,'CTA':0.0,'CTC':0.0,'CTG':0.0,'CTT':0.0,'GAA':0.0,'GAC':0.0,'GAG':0.0, \
                                    'GAT':0.0,'GCA':0.0,'GCC':0.0,'GCG':0.0,'GCT':0.0,'GGA':0.0,'GGC':0.0,'GGG':0.0,'GGT':0.0,'GTA':0.0,'GTC':0.0,'GTG':0.0, \
                                    'GTT':0.0,'TAA':0.0,'TAC':0.0,'TAG':0.0,'TAT':0.0,'TCA':0.0,'TCC':0.0,'TCG':0.0,'TCT':0.0,'TGA':0.0,'TGC':0.0,'TGG':0.0, \
                                    'TGT':0.0,'TTA':0.0,'TTC':0.0,'TTG':0.0,'TTT':0.0}}   


    # Mutate (substitutions) all possible codons in the given genetic code, and count proportions of mutations that are synonymous and non-synonmyous
    for codon in nt_to_aa.keys():
        for codon_p in range(0,2+1):

            nts = ['A','G','T','C']  # @DONE: refactor away, A: we can't, since the next line

            nts.remove(codon[codon_p]) # we do not consider self substitutions, e.g. A->A

            # ...and for each nucleotide that the bp can change 
            # into... 
            for nt in nts:

                codon_mutated = list(copy.deepcopy(codon))
                #codon_mutated = codon
                codon_mutated[codon_p] = nt  # mutate the basepair
                codon_mutated = ''.join(codon_mutated)
                
                # ...count how many of them are synonymous.
                if nt_to_aa[codon]==nt_to_aa[codon_mutated]:
                    potential_changes['S'][codon]+=1.0/3.0 #@DONE: Q: but why 1/3? to me it should be 1/(3*4), A: since it represents 1/3 of a "site"
                else:
                    potential_changes['N'][codon]+=1.0/3.0 #@DONE: Q: but why 1/3? to me it should be 1/(3*4), A: since it represents 1/3 of a "site"

    codons      = nt_to_aa.keys()
    codonPairs  = list(permutations(codons,2))
    selfies     = [(i,i) for i in codons]
    codonPairs  = codonPairs + selfies 
    
    codonPair_to_potential = {}

    for pair in codonPairs:

        codon1 = pair[0]
        codon2 = pair[1]
        pn1 = potential_changes['N'][codon1]
        pn2 = potential_changes['N'][codon2]
        ps1 = potential_changes['S'][codon1]
        ps2 = potential_changes['S'][codon2]
        codonPair_to_potential[pair] = {'N':(pn1+pn2)/2.,'S':(ps1+ps2)/2.}

    f = open('potential_changes_dict.p','w')
    pickle.dump(codonPair_to_potential,f)
    f.close()

def observed_changes_dict(nt_to_aa):

    codons      = nt_to_aa.keys()
    codonPairs  = list(permutations(codons,2))
    selfies     = [(i,i) for i in codons]
    codonPairs  = codonPairs + selfies 
    
    codonPair_to_observed = {}

    for pair in codonPairs:
        codon1 = pair[0]
        codon2 = pair[1]
        indices_to_permute = []

        # Collect the position of the letters (1, 2, 3) where the two sequences differ... 
        for letter_i in range(0,3):
            if not codon1[letter_i] == codon2[letter_i]:
                indices_to_permute.append(letter_i)

        # We now have all the possible mutational pathways, represented as indices 
        permuted_indices = list(permutations(indices_to_permute))
        syn = []
        non = []

        for i,path in enumerate(permuted_indices):
            syn.append(int()) 
            non.append(int()) 

            codon1_path1 = list(codon1) # copies of seqs for 'mutating'
            codon2_path1 = list(codon2)

            for site in path:

                codon1_past         = ''.join(codon1_path1)
                codon1_path1[site]  = codon2_path1[site]        # s1 = 'TTT' , s2 = 'ATA'  ==> 'TTT' --> 'ATT' 
                codon1_path1        = ''.join(codon1_path1)
                if nt_to_aa[codon1_path1] == nt_to_aa[codon1_past]:  # 'TTT --> 'ATT'
                    syn[i] = syn[i] + 1 
                    non[i] = non[i] + 0
                else:
                    syn[i] = syn[i] + 0
                    non[i] = non[i] + 1
                codon1_path1 = list(codon1_path1)

        try:
            assert isclose(np.mean(syn)+np.mean(non),float(len(path)))
        except AssertionError:
            raise ValueError("Calculations are incorrect, mutation pathways calculation failed...")


        codonPair_to_observed[pair] = {'S':np.mean(syn),'N':np.mean(non)}
    f2 = open('observed_changes_dict.p','wb')
    pickle.dump(codonPair_to_observed,f2)
    f2.close()    

def dnds( seq1, seq2, changes_potential, changes_observed, msCorrect='approximate', sliding=False, windowLength=3, stepLength=1):
    """ Perform dN/dS analysis, using the 'NG' algoritm, includes both whole sequence or sliding window, and either an approximate or exact multiple-substiution correction method. (@todo: make sure it actually is exact... it could be
             something else)
    ARGS:
        seq1,  a DNA sequence as string of letters, AGTC. Seq1 must be equal in length 
            to, and aligned with, seq2, with gaps trimmed away. @todo: how on earth can 
            we reliably make this work on the web service?

        seq2,  a DNA sequence similar to seq1 but with differences (substitutions), 
            representing a CDS orthologue of seq1 from a different species. Read 
            description of seq1 for other required similarities to avoid errors.

        changes_potential, a dict, with key=pair of codons tuple, e.g. ('ATG','ATG'), and value=['S':<S>,'N':<N>]. Where <S> is the number of potential 
        synonmyous sites for each codon (averaged between the two codons), and <N> is the same but for non-synonymous sites
            e.g. changes.potential_changes_dict(...)  (see: ./changes.py)
        changes_observed, @todo
        msCorrect, a string to toggle between multiple-substitution correction methods:
            "approximate", "exact" (@todo: make sure it actually is exact... it could be
             something else)
            e.g. msCorrect = 'approximate'
        sliding, a boolean to toggle between sliding window analysis (vector of dN/dS values at successive chunks of sequence) or whole sequence analysis (a single 
            dN/dS value for the given pair of input sequences), either: True, False
            e.g. sliding = False
        windowLength, an integer specifying the width of the sliding window, measured in no. of codons in the window to measure dN/dS over, from 1-to-length(seq1)
            e.g. windowLength = 50
        stepLength, an integer specifying no. of codons to shift the sliding window with each iteration. If stepLength < windowLength then windows will overlap, overlapping is dealt with prior to plotting (acts to smooth values, averages along the overlaps are taken as a dN/dS value for any codon).
            e.g. stepLength = 1
    """

    def chunks(l, n):
        """ Yield successive n-sized chunks from l. """
        for i in xrange(0, len(l), n):
            yield l[i:i+n]

    # todo: stop codons to deal with, reject
    # todo: ambiguous bases to deal with: 
        # gaps,  @done 
        # Ns,    @todo
        # Xs     @todo

    warning_count = 0

    # STATS per CODON-PAIR:
    codons_seq1   = [codon for codon in chunks(seq1,3)]  #splits
    codons_seq2   = [codon for codon in chunks(seq2,3)]
    print(codons_seq1)
    print(len(codons_seq1))
    print(codons_seq2)
    print(len(codons_seq2))
    codons_paired = [pair for pair in zip(codons_seq1,codons_seq2) if (len(pair[0])+len(pair[1]))==6] # aligned codons are paired into tuples, excess codons are truncated, @todo: in main example, we lose 5 bps of data
    print(codons_paired)
    print(len(codons_paired))
    changes_all = {'observed':{'S':[],'N':[]},'potential':{'S':[],'N':[]}}

    for pair in codons_paired:
        changes_all['potential']['S'].append(changes_potential[pair]['S'])
        changes_all['potential']['N'].append(changes_potential[pair]['N'])
        changes_all['observed']['S'].append(changes_observed[pair]['S'])
        changes_all['observed']['N'].append(changes_observed[pair]['N'])

    list_S  = changes_all['potential']['S']
    list_Sd = changes_all['observed']['S']

    list_N  = changes_all['potential']['N']
    list_Nd = changes_all['observed']['N']

    if sliding:
        # STATS for each WINDOW seq
        intervals    = range(0,len(codons_paired)-windowLength+1,stepLength)
        print(intervals)
        windows      = zip(intervals,[i + windowLength - 1 for i in intervals]) 
        print(windows)
        window_stats = {}

        #window_stats_list = []

        # @done: test against matlab's sliding window, also @todo: find out what stepLength does, @todo: try to plot the sliding window version

        for window_i,window in enumerate(windows):
            print(window_i,window)
            start = window[0]
            end   = window[1]+1
            print(start,end)

            window_stats[window] = {    'S':sum(list_S[start:end]),
                                        'Sd':sum(list_Sd[start:end]),
                                        'N': sum(list_N[start:end]),
                                        'Nd':sum(list_Nd[start:end])    }


            pS = window_stats[window]['Sd']/window_stats[window]['S']
            pN = window_stats[window]['Nd']/window_stats[window]['N']

            # try:
            if msCorrect=='approximate':
                dN = -(3./4.)*math.log(1.-(4./3.)*pN)
                dS = -(3./4.)*math.log(1.-(4./3.)*pS)
                print(pS,dS)

            else: 
                dN = pN
                dS = pS
            window_stats[window]['dNdS'] = dN/dS
            # except ZeroDivisionError:
                # warning_count += 1
                # window_stats[window]['dNdS'] = float('Inf')
            # except ValueError:
                # warning_count += 1
                # window_stats[window]['dNdS'] = float('nan')

        return window_stats,warning_count  # list of dnds per window interval // dict of dnds, key=(<from #base pair>,<to #base pair>), value=<dN/dS of the window specified in the key>
    else:
        # STATS for WHOLE SEQ
        S   = sum(list_S)
        Sd  = sum(list_Sd)
        pS  = Sd/S 
        N   = sum(list_N)
        Nd  = sum(list_Nd)
        pN  = Nd/N

        # try:
        if msCorrect=='approximate':

            if (pS>=3./4.):
                raise ValueError("Query and reference sequences are too divergent. Approximate multiple-substitutions correction cannot be achieved, SYNONYMOUS changes per synonymous site, pS>=3/4, log() operation will yeild return undefined. Try alternative value for argument: msCorrect (e.g. 'exact')...") 

            if (pN>=3./4.):
                raise ValueError("Query and reference sequences are too divergent. Approximate multiple-substitutions correction cannot be achieved, NON-SYNONYMOUS changes per synonymous site, pN>=3/4, log() operation will yeild return undefined. Try alternative value for argument: msCorrect (e.g. 'exact')...") 

            dS  = -(3./4.)*math.log(1.-((4./3.)*pS))
            dN  = -(3./4.)*math.log(1.-((4./3.)*pN))
            dN_dS = dN/dS

        else: # @todo: is this the exact one? Or something else? 
            
            # @DONE: one day the following three lines of code will error, giving a ZeroDivisionError, this needs to be handled with try
            dS = pS  # i.e. dS = Sd/S
            dN = pN
            dN_dS = dN/dS
        # except ValueError:
            # warning_count += 1
            # warnings.warn("Query and reference sequencea are too divergent. ValueError: Approximate multiple-substitutions correction cannot be achieved: UNKNOWN reason, probably due to illegal numbers in a log() function...\n") 
            # dN_dS = float("nan")
        # except ZeroDivisionError:
            # warning_count += 1
            # warnings.warn("Query and reference sequences are too divergent. ZeroDiviSionError: Approximate multiple-substitutions correction cannot be achieved: UNKNOWN reason, probably due to illegal numbers in a log() function...\n") 
            # dN_dS = float('Inf')

        return dN_dS, warning_count  # i.e. omega = dN/dS = (Nd/N)/(Sd/S)


def plot_dnds_sliding(dnds_slide_dict):
  
    window_intervals = dnds_slide_dict.keys() # @todo: I dont think sorting the windows will make a difference to final result, but it will make it slower, @todo: test this just in case

    # print(window_intervals)
    max_window_position = np.amax(window_intervals) # e.g. 243 @done: a whole order of magnitude faster than: max_window_interval = max(window_intervals, key=lambda x: x[1])
    overlap_matrix      = np.empty((len(window_intervals),max_window_position+1))  # @todo: are you sure it's +1? initialize empty matrix, note: entries are not actually NaN, just near-zero
    overlap_matrix[:]   = np.NAN # initiate empty np array with NaN, so later we can mask

    for window_i,window in enumerate(window_intervals):

        start = window[0] # e.g. 0
        end   = window[1] # e.g. 49 

        # in the i-th row, fill all elements from the window[0]-th to window[1]-th with the dN/dS value for this window
        overlap_matrix[window_i,start:end+1] = dnds_slide_dict[window]['dNdS'] # @todo: are you sure it's +1? test, keep in mind for these indices it does -1 for the "to" part

    nan_masker              = ~np.isfinite(overlap_matrix) # boolean matrix, True if element is finite, False if element is Inf or NaN
    overlap_matrix_masked   = np.ma.masked_array(overlap_matrix,mask=nan_masker)
    overlap_matrix_avg      = overlap_matrix_masked.mean(axis=0)
    return list(overlap_matrix_avg), overlap_matrix_avg.mean()

def dnds_pipeline(qry_seq_in, ref_seq_in):
    try:       
        f1 = open('observed_changes_dict.p','rb')
        observed_changes  = pickle.load(f1)
        f1.close()

        f2 = open('potential_changes_dict.p','rb')
        potential_changes = pickle.load(f2)
        f2.close()
        q=[len(observed_changes),len(potential_changes),type(observed_changes),type(potential_changes)]
        raw_input(' '.join(map(str,q)))
        for item in observed_changes:
            print(item,observed_changes[item])
        for item in potential_changes:
            print(item,potential_changes[item])
    except IOError:
        sys.exit("cannot find potential_changes_dict or observed_changes_dict. Please palce in CWD") # @todo:REMOVE
    qry_seq_raw = qry_seq_in
    ref_seq_raw = ref_seq_in
    aln_gap_open = -10      
    aln_gap_extend = -0.5 # @todo: make cmd args later?

    qry_seq_aln, ref_seq_aln                          = align_query_vs_reference(qry_seq_raw, ref_seq_raw, aln_gap_open, aln_gap_extend)
    qry_seq_trimmed, ref_seq_trimmed, qry_seq_indices = trim_gaps_from_aligned_seqs( qry_seq_aln, ref_seq_aln )

    dnds_slide_dict, warning_count = dnds( qry_seq_trimmed, ref_seq_trimmed, potential_changes, observed_changes, msCorrect='approximate', sliding=True, windowLength=25, stepLength=1 )
    #
    # Plot the sliding window values
    #
    dnds_sliding_vec, dnds_sliding_mean = plot_dnds_sliding(dnds_slide_dict)

    print("\tElapsed time: "+str(time.time() - start_time))   # @time
    print("\tTotal warnings (missing values): "+str(warning_count))
    print "\tAvg. dN/dS over all windows: "+str(dnds_sliding_mean)
    print("")

    return dnds_sliding_mean

if __name__ == "__main__":
    try:       
        f1 = open('observed_changes_dict.p','rb')
        observed_changes  = pickle.load(f1)
        f1.close()

        f2 = open('potential_changes_dict.p','rb')
        potential_changes = pickle.load(f2)
        f2.close()
        print "LOADED!!" # @todo:REMOVE

    except IOError:
        nt_to_aa_dict     = geneticCode("standard")

        potential_changes_dict(nt_to_aa_dict)
        observed_changes_dict(nt_to_aa_dict)

        f1 = open('observed_changes_dict.p','rb') # @todo:coderedundancey = bad, wrap these lines into a function and call the function
        observed_changes  = pickle.load(f1)
        f1.close()

        f2 = open('potential_changes_dict.p','rb')
        potential_changes = pickle.load(f2)
        f2.close()

        print "CREATED!!" # @todo:REMOVE

    # qry_seq_raw = sys.argv[1]
    qry_seq_raw = "ATGTGTGAATATTCGAATACGCGCAACAAGATGAGCAACCTGGTCGTCGTCCTCGTCCTGCTGACGATGTACATTGTGCTTTCGGCCCCATTCGAAATACCGGACCGGTACAAAAAGCCGGCTAAAATGTTGCACGAAATTTGTATCGCCGAGTCGGGCGCCTCGGAGGAGCAGCTGCGCACCTGTCTCGATGGAACCGTACCGACAGCTCCGGCCGCCAAGTGCTACATCCACTGCCTGTTCGACAAGATCGACGTGGTGGACGAGGCGACTGGGCGCATCCTGCTCGACCGACTGCTTTACATCATCCCGGACGACGTGAAGGCAGCGGTGGACCATTTAACGCGCGAATGTAGCCACATCGTAACGCCGGATAAGTGCGAAACCGCCTACGAGACGGTCAAATGTTATTTCAATGCGCACGACGAGGTGATCAAATTCTGCCACCTACTAGTGCTGGAGTGA"
    # ref_seq_raw = sys.argv[2]
    ref_seq_raw = "ATGATGGAACAGCTTATGCTGGCAGTTTTGCTGGCGGTTTTTCTCGGGCTCGTAGCAGATGTTACGATGGCCGCTCAAATCAAGGACAATTTGGAGCTACCCGAATATTACAAACGTCCGGCCAAAATTCTGCACAACATCTGTCTGGCAGAATCCGGTGCCATGGAGAGCAAACTAAAGCAGTGCATGGACGGAGTGCTTCATGACGACCGGGAAGTCAAGTGCTACATCCATTGTCTATTCGACAAGGTGGACGTAATCGACGAAGCAACCGGGCAGATCCTGTTGGACCGATTGGCACCACTGGCACCGGACAACGATGTGAAGGATGTGTTCAATCATTTGACCAAAGAGTGTGGTCATATCAAACTACAAGATTCCTGCGATACGGCGTACGAAGTGGCCAAATGTTACTTCGCGGCACACGATCAGGTCGTCAAATTCTGTCACCTGTTGATGGCTGATGTTACCAGCTAG"
    dnds_pipeline(qry_seq_raw,ref_seq_raw)
