import numpy as np 
import sys
import pickle
import math
import warnings
import changes 
import time 
start_time = time.time()  # @time

def dnds( seq1, seq2, changes_potential, changes_observed, msCorrect='approximate', sliding=False, windowLength=3, stepLength=1):
    def chunks(l, n):
        for i in xrange(0, len(l), n):
            yield l[i:i+n]
    warning_count = 0
    codons_seq1   = [codon for codon in chunks(seq1,3)]  #splits
    codons_seq2   = [codon for codon in chunks(seq2,3)]  
    codons_paired = [pair for pair in zip(codons_seq1,codons_seq2) if (len(pair[0])+len(pair[1]))==6] # aligned codons are paired into tuples, excess codons are truncated, @todo: in main example, we lose 5 bps of data
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
        intervals    = range(0,len(codons_paired)-windowLength+1,stepLength)
        windows      = zip(intervals,[i + windowLength - 1 for i in intervals]) 
        window_stats = {}
        for window_i,window in enumerate(windows):
            start = window[0]
            end   = window[1]+1
            window_stats[window] = {    'S':sum(list_S[start:end]),
                                        'Sd':sum(list_Sd[start:end]),
                                        'N': sum(list_N[start:end]),
                                        'Nd':sum(list_Nd[start:end])    }
            pS = window_stats[window]['Sd']/window_stats[window]['S']
            pN = window_stats[window]['Nd']/window_stats[window]['N']
            try:
                if msCorrect=='approximate':
                    dN = -(3./4.)*math.log(1.-(4./3.)*pN)
                    dS = -(3./4.)*math.log(1.-(4./3.)*pS)
                else:
                    dN = pN
                    dS = pS
                window_stats[window]['dNdS'] = dN/dS
            except ZeroDivisionError:
                warning_count += 1
                window_stats[window]['dNdS'] = float('Inf')
            except ValueError:
                warning_count += 1
                window_stats[window]['dNdS'] = float('nan')
        return window_stats,warning_count  # list of dnds per window interval // dict of dnds, key=(<from #base pair>,<to #base pair>), value=<dN/dS of the window specified in the key>
    else:
        S   = sum(list_S)
        Sd  = sum(list_Sd)
        pS  = Sd/S 
        N   = sum(list_N)
        Nd  = sum(list_Nd)
        pN  = Nd/N
        try:
            if msCorrect=='approximate':
                if (pS>=3./4.):
                    raise ValueError("Query and reference sequences are too divergent. Approximate multiple-substitutions correction cannot be achieved, SYNONYMOUS changes per synonymous site, pS>=3/4, log() operation will yeild return undefined. Try alternative value for argument: msCorrect (e.g. 'exact')...") 
                if (pN>=3./4.):
                    raise ValueError("Query and reference sequences are too divergent. Approximate multiple-substitutions correction cannot be achieved, NON-SYNONYMOUS changes per synonymous site, pN>=3/4, log() operation will yeild return undefined. Try alternative value for argument: msCorrect (e.g. 'exact')...") 
                dS  = -(3./4.)*math.log(1.-((4./3.)*pS))
                dN  = -(3./4.)*math.log(1.-((4./3.)*pN))
                dN_dS = dN/dS
            else:
                dS = pS  # i.e. dS = Sd/S
                dN = pN
                dN_dS = dN/dS
        except ValueError:
            warning_count += 1
            warnings.warn("Query and reference sequencea are too divergent. ValueError: Approximate multiple-substitutions correction cannot be achieved: UNKNOWN reason, probably due to illegal numbers in a log() function...\n") 
            dN_dS = float("nan")
        except ZeroDivisionError:
            warning_count += 1
            warnings.warn("Query and reference sequences are too divergent. ZeroDiviSionError: Approximate multiple-substitutions correction cannot be achieved: UNKNOWN reason, probably due to illegal numbers in a log() function...\n") 
            dN_dS = float('Inf')
        return dN_dS, warning_count  # i.e. omega = dN/dS = (Nd/N)/(Sd/S)

def plot_dnds_sliding(dnds_slide_dict):
    window_intervals = dnds_slide_dict.keys() # @todo: I dont think sorting the windows will make a difference to final result, but it will make it slower, @todo: test this just in case
    max_window_position = np.amax(window_intervals) # e.g. 243 @done: a whole order of magnitude faster than: max_window_interval = max(window_intervals, key=lambda x: x[1])
    overlap_matrix      = np.empty((len(window_intervals),max_window_position+1))  # @todo: are you sure it's +1? initialize empty matrix, note: entries are not actually NaN, just near-zero
    overlap_matrix[:]   = np.NAN # initiate empty np array with NaN, so later we can mask
    for window_i,window in enumerate(window_intervals):
        start = window[0] # e.g. 0
        end   = window[1] # e.g. 49 
        overlap_matrix[window_i,start:end+1] = dnds_slide_dict[window]['dNdS'] # @todo: are you sure it's +1? test, keep in mind for these indices it does -1 for the "to" part
    nan_masker              = ~np.isfinite(overlap_matrix) # boolean matrix, True if element is finite, False if element is Inf or NaN
    overlap_matrix_masked   = np.ma.masked_array(overlap_matrix,mask=nan_masker)
    overlap_matrix_avg      = overlap_matrix_masked.mean(axis=0)
    return list(overlap_matrix_avg), overlap_matrix_avg.mean()

def dnds_pipeline(qry_seq_in, ref_seq_in):
    print("\tdN/dS sliding analysis: pre-cached statistics for all possible codon pairs...")
    try:       
        f1 = open('../data/observed_changes_dict.p','rb')
        observed_changes  = pickle.load(f1)
        f1.close()
        f2 = open('../data/potential_changes_dict.p','rb')
        potential_changes = pickle.load(f2)
        f2.close()
        print("\t\tLOADED!") # @todo:REMOVE
    except IOError:
        nt_to_aa_dict     = changes.geneticCode("standard")
        changes.potential_changes_dict(nt_to_aa_dict)
        changes.observed_changes_dict(nt_to_aa_dict)
        f1 = open('../data/observed_changes_dict.p','r') # @todo:coderedundancey = bad, wrap these lines into a function and call the function
        observed_changes  = pickle.load(f1)
        f1.close()
        f2 = open('../data/potential_changes_dict.p','r')
        potential_changes = pickle.load(f2)
        f2.close()
    qry_seq_raw = qry_seq_in
    ref_seq_raw = ref_seq_in
    print("\tProcessing heuristic codon alignments....")
    aln_gap_open = -10      
    aln_gap_extend = -0.5 # @todo: make cmd args later?
    qry_seq_aln, ref_seq_aln                          = align_query_vs_reference(qry_seq_raw, ref_seq_raw, aln_gap_open, aln_gap_extend)
    qry_seq_trimmed, ref_seq_trimmed, qry_seq_indices = trim_gaps_from_aligned_seqs( qry_seq_aln, ref_seq_aln )
    print("\tSliding window analysis....")
    dnds_slide_dict, warning_count = dnds( qry_seq_trimmed, ref_seq_trimmed, potential_changes, observed_changes, msCorrect='approximate', sliding=True, windowLength=25, stepLength=1 )
    for item in dnds_slide_dict:
        print(item,dnds_slide_dict[item])
    raw_input()
    dnds_sliding_vec, dnds_sliding_mean = plot_dnds_sliding(dnds_slide_dict)
    print("\tElapsed time: "+str(time.time() - start_time))   # @time
    print("\tTotal warnings (missing values): "+str(warning_count))
    print("\tAvg. dN/dS over all windows: "+str(dnds_sliding_mean))
    return dnds_sliding_vec, qry_seq_indices

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
   print(type(qry_trimmed))
   print(len(qry_trimmed))
   print(qry_trimmed)
   print(type(ref_trimmed))
   
   print(len(ref_trimmed))
   print(ref_trimmed)
   return (qry_trimmed, ref_trimmed)
   # return 1

    
if __name__ == "__main__":
    try:       
        f1 = open('./py/data/observed_changes_dict.p','rb')
        observed_changes  = pickle.load(f1)
        f1.close()
        f2 = open('./py/data/potential_changes_dict.p','rb')
        potential_changes = pickle.load(f2)
        f2.close()
        print("LOADED!!")# @todo:REMOVE
    except IOError:
        nt_to_aa_dict     = changes.geneticCode("standard")
        changes.potential_changes_dict(nt_to_aa_dict)
        changes.observed_changes_dict(nt_to_aa_dict)
        f1 = open('../data/observed_changes_dict.p','rb') # @todo:coderedundancey = bad, wrap these lines into a function and call the function
        observed_changes  = pickle.load(f1)
        f1.close()

        f2 = open('../data/potential_changes_dict.p','rb')
        potential_changes = pickle.load(f2)
        f2.close()
        print("CREATED!!")# @todo:REMOVE
    # qry_seq_raw = sys.argv[1]
    qry_seq_raw = "ATGTGTGAATATTCGAATACGCGCAACAAGATGAGCAACCTGGTCGTCGTCCTCGTCCTGCTGACGATGTACATTGTGCTTTCGGCCCCATTCGAAATACCGGACCGGTACAAAAAGCCGGCTAAAATGTTGCACGAAATTTGTATCGCCGAGTCGGGCGCCTCGGAGGAGCAGCTGCGCACCTGTCTCGATGGAACCGTACCGACAGCTCCGGCCGCCAAGTGCTACATCCACTGCCTGTTCGACAAGATCGACGTGGTGGACGAGGCGACTGGGCGCATCCTGCTCGACCGACTGCTTTACATCATCCCGGACGACGTGAAGGCAGCGGTGGACCATTTAACGCGCGAATGTAGCCACATCGTAACGCCGGATAAGTGCGAAACCGCCTACGAGACGGTCAAATGTTATTTCAATGCGCACGACGAGGTGATCAAATTCTGCCACCTACTAGTGCTGGAGTGA"
    # ref_seq_raw = sys.argv[2]
    ref_seq_raw = "ATGATGGAACAGCTTATGCTGGCAGTTTTGCTGGCGGTTTTTCTCGGGCTCGTAGCAGATGTTACGATGGCCGCTCAAATCAAGGACAATTTGGAGCTACCCGAATATTACAAACGTCCGGCCAAAATTCTGCACAACATCTGTCTGGCAGAATCCGGTGCCATGGAGAGCAAACTAAAGCAGTGCATGGACGGAGTGCTTCATGACGACCGGGAAGTCAAGTGCTACATCCATTGTCTATTCGACAAGGTGGACGTAATCGACGAAGCAACCGGGCAGATCCTGTTGGACCGATTGGCACCACTGGCACCGGACAACGATGTGAAGGATGTGTTCAATCATTTGACCAAAGAGTGTGGTCATATCAAACTACAAGATTCCTGCGATACGGCGTACGAAGTGGCCAAATGTTACTTCGCGGCACACGATCAGGTCGTCAAATTCTGTCACCTGTTGATGGCTGATGTTACCAGCTAG"

    aln_gap_open = -10      
    aln_gap_extend = -0.5 # @todo: make cmd args later?

    # qry_seq_aln, ref_seq_aln            = align_query_vs_reference(qry_seq_raw, ref_seq_raw, aln_gap_open, aln_gap_extend)
    # (query, reference)= trim_gaps_from_aligned_seqs( qry_seq_aln, ref_seq_aln )
    # a = trim_gaps_from_aligned_seqs( qry_seq_aln, ref_seq_aln )
    a,b=dnds_pipeline('ATGTGTGAATATTCGAATACGCGCAACAAGATGAGCAACCTGGTCGTCGTCCTCGTCCTGCTGACGATGTACATTGTGCTTTCGGCCCCATTCGAAATACCGGACCGGTACAAAAAGCCGGCTAAAATGTTGCACGAAATTTGTATCGCCGAGTCGGGCGCCTCGGAGGAGCAGCTGCGCACCTGTCTCGATGGAACCGTACCGACAGCTCCGGCCGCCAAGTGCTACATCCACTGCCTGTTCGACAAGATCGACGTGGTGGACGAGGCGACTGGGCGCATCCTGCTCGACCGACTGCTTTACATCATCCCGGACGACGTGAAGGCAGCGGTGGACCATTTAACGCGCGAATGTAGCCACATCGTAACGCCGGATAAGTGCGAAACCGCCTACGAGACGGTCAAATGTTATTTCAATGCGCACGACGAGGTGATCAAATTCTGCCACCTACTAGTGCTGGAGTGA','ATGATGGAACAGCTTATGCTGGCAGTTTTGCTGGCGGTTTTTCTCGGGCTCGTAGCAGATGTTACGATGGCCGCTCAAATCAAGGACAATTTGGAGCTACCCGAATATTACAAACGTCCGGCCAAAATTCTGCACAACATCTGTCTGGCAGAATCCGGTGCCATGGAGAGCAAACTAAAGCAGTGCATGGACGGAGTGCTTCATGACGACCGGGAAGTCAAGTGCTACATCCATTGTCTATTCGACAAGGTGGACGTAATCGACGAAGCAACCGGGCAGATCCTGTTGGACCGATTGGCACCACTGGCACCGGACGATGTGAAGGATGTGTTCAATCATTTGACCAAAGAGTGTGGTCATATCAAACTACAAGATTCCTGCGATACGGCGTACGAAGTGGCCAAATGTTACTTCGCGGCACACGATCAGGTCGTCAAATTCTGTCACCTGTTGATGGCTGACTAG')
    # dnds_slide_dict, warning_count      = dnds( qry_seq_trimmed, ref_seq_trimmed, potential_changes, observed_changes, msCorrect='approximate', sliding=True, windowLength=50, stepLength=1 )
    # dnds_sliding_vec, dnds_sliding_mean = plot_dnds_sliding(dnds_slide_dict)
    print('-')
    print(a)
    print('-')
    print(b)
    print('-')
    # print "Mean dN/dS over all windows: "+str(dnds_sliding_mean)
