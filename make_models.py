from numpy import array

kmeans_models={'5':['meanConsensus','degree_assortativity','inscape_prot_kmer_3','phacelia_score','corrPageRankfreq'],
'10_1':['atchley_st1','max_degree','VonEntropy','KSEntropy','corrPageRankfreq','spectralRadiusAdj','inscape_nuc_kmer_7','inscape_prot_kmer_3','atchley_wa0','freq_corr'],
'13':['protein_nuc_entropy','pelin_haplo_freq','corrPageRankfreq','max_degree','corr_path_hamm_dist','inscape_prot_kmer_3','inscape_nuc_kmer_7','KSEntropy','VonEntropy','spectralRadiusAdj','degree_assortativity','atchley_wa0','num_1step_components'],
'16':['atchley_st1','meanConsensus','corrPageRankfreq','max_degree','corr_path_hamm_dist','KSEntropy','degree_assortativity','spectralRadiusAdj','pelin_haplo_freq','VonEntropy','freq_corr','atchley_wa0','inscape_prot_kmer_3','atchley_wa3','atchley_st2','RMSE_path_hamm_dist'],
'17':['trans_mut','atchley_st3','atchley_wa0','meanConsensus','corrPageRankfreq','max_degree','genetic_load','inscape_prot_kmer_3','degree_assortativity','VonEntropy','spectralRadiusAdj','protein_nuc_entropy','KSEntropy','corr_path_hamm_dist','atchley_wa3','pelin_haplo_freq','atchley_st1'],
'12':['VonEntropy','KSEntropy','degreeEntropy','std_dev','meanConsensus','degree_assortativity','inscape_nuc_kmer_7','inscape_prot_kmer_3','phacelia_score','atchley_wa0','corr_bet_cent_freq','corrPageRankfreq'],
'10_2':['VonEntropy','std_dev','meanConsensus','degree_assortativity','inscape_nuc_kmer_7','inscape_prot_kmer_3','phacelia_score','atchley_st3','atchley_wa0','corrPageRankfreq']}

def ml_data_parser_for_printing_models(file):
    # print(file)
    data_end=-3
    groups_col=-2
    c=0
    X=[]
    y=[]
    g1=[]
    g2=[]
    samples=[]
    with open(file) as f:
        for line in f.readlines():
            c+=1
            if c!=1:
                s=line.strip().split(',') 
                xx=s[1:data_end]
                samples.append(s[0])
                X.append(list(map(float,xx)))
                g1.append(int(s[-1]))
                g2.append(int(s[-2]))
                y.append(int(s[data_end]))
            else:
                s=line.strip().split(',')
                names=s[1:data_end]
                # groupname=s[groups_col]
    return array(X),array(y),names,samples,g1,g2

def get_indices(names,candidate):
    out=[]
    for i in range(len(names)):
        name=names[i]
        if name in candidate:
            out.append(i)
    return out
    
def list2str_special(x):
    return ','.join(map(str,x))+'\n'
    
X,y,names,samples,g1,g2=ml_data_parser_for_printing_models('74_data.csv')
for model in kmeans_models:
    out_data=[]
    indices=get_indices(names,kmeans_models[model])
    our_x=X[:,indices]
    file=model+'_model.csv'
    with open(file,'w') as f:
        f.write(list2str_special(kmeans_models[model]))
        for i in range(len(X)):
            row=list(our_x[i])
            row.insert(0,samples[i])
            row.append(y[i])
            row.append(g2[i])
            row.append(g1[i])
            f.write(list2str_special(row))
            