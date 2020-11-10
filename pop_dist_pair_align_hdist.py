#!/usr/bin/env python

__author__ = 'David S. Campo, Ph.D.'
'''
pop_dist_pair_align. Program for calculating hamming distances between unaligned files.
1. It requires as input the location of the fasta files. 
3. It saves a .csv file, whre for every pair of files, it calculates the average hamming distance, the minimal distance, the maximum distance and the haplotype overlap.

Activate ghost environment
source /scicomp/groups/OID/NCHHSTP/DVH/AMD/ghost3_stable4/bin/activate



Usage with default parameters:
python pop_dist_pair_align.py -i input_folder_location -o output_folder_location

Usage,if you want to align each file pair:
python pop_dist_pair_align.py -i input_folder_location -o output_folder_location -a 1

'''

import time
import os, re
import numpy as np
#import math
import glob
import optparse
from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import generic_dna
#from collections import Counter, OrderedDict
#
#
#from scipy import spatial
#from Bio.Align.Applications import MafftCommandline
#import operator
#import scipy.sparse as ss
#from scipy.sparse.csgraph import connected_components
#from scipy.sparse.csgraph import csgraph_from_dense
from ghost.util.distance import hamming
#from pyseqdist import hamming
from scipy.stats.mstats import gmean
from scipy.stats import spearmanr


def calcDistanceMatrix(seqs1, seqs2): #calculate distance matrix from the 1-step list
    arr=[]
    hdist=hamming(seqs1, seqs2,ignore_gaps=False)
    for i in hdist:
        for j in i:
            arr.append(j[0])
    return arr

if __name__ == '__main__':
    try:
        argument_parser = optparse.OptionParser()            
        argument_parser.add_option('-i', metavar='IDIR', action='store', type=str, dest='input_directory', default='.', help='Input file')
        argument_parser.add_option('-o', metavar='ODIR', action='store', type=str, dest='output_directory', default='.', help='Output directory')  
        argument_parser.add_option('-a', action='store', type=int, dest='alignFlag', default= 0, help='align?')  
        argument_parser.add_option('-w', action='store', type=int, dest='withinFlag', default= 0, help='do corrected?')                                  
        argument_parser.add_option('-c', action='store', type=int, dest='compFlag', default= 0, help='do only components within a patient?')                                  

        options, args = argument_parser.parse_args()    
        pathInput = options.input_directory
        pathOutput = options.output_directory    
        alignFlag = options.alignFlag
        withinFlag = options.withinFlag
        compFlag = options.compFlag
        
        if not os.path.exists(pathOutput):
            os.makedirs(pathOutput)

        fileList = sorted(glob.glob(os.path.join(pathInput, '*.fas')))
        fileNum = len(fileList)
        nComp = ((fileNum*fileNum)- fileNum) / 2
        print('Number of pairwise comparisons =', nComp)
        #allFileNames1 = [['array%i' %i] for i in range(nComp)]
        #allFileNames2 = [['array%i' %i] for i in range(nComp)]
  
        
        ## open output files
        outputName1 = pathOutput + '/allDistances.csv'                
        outputHandle1 = open(os.path.join(outputName1), "w")
        outputHandle1.write('fileName1,fileName2,patientFlag,haploSize,majorDist,avgDist,minDist,maxDist,corrDist,geoDist,weiDist,hapOverlap,hapOverlapWei,thresholdPair,thresholdPairWei,thresholdFile,sharedCorrP,sharedCorrS' + '\n')


        outputName2 = pathOutput + '/chosenDist.csv'                
        outputHandle2 = open(os.path.join(outputName2), "w")
        outputHandle2.write('fileName1,fileName2,dist' + '\n')

        # within file 1
        if withinFlag == 1:
            #print('with within')
            f = 0
            withinHamming = []
            for pathFileInput1 in fileList:
                dirInput1,fileInput1 = os.path.split(pathFileInput1)
                fileName1, fileExtension1 = os.path.splitext(fileInput1)
                ##parse input             
                reads0 = []
                counts0 = []
                inputHandle = open(pathFileInput1) 
                for record in SeqIO.parse(inputHandle, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
                        reads0.append(str(record.seq))
                        #counts.append(1) 
                        #takes into account the frequencies
                        counts0.append(int(re.split('_', record.id)[-1])) 
                inputHandle.close()
                haploSize0 = len(reads0[0])

                withinFileList = calcDistanceMatrix(reads0, reads0)
#                for haplotype1 in reads0:    
#                    for haplotype2 in reads0:
#                        rawDist = sum(1 for a, b in zip(haplotype1, haplotype2) if a != b)
#                        withinFileList.append(rawDist)
                withinHamming.append(np.divide(np.mean(withinFileList), haploSize0, dtype = float))
        else:
            withinHamming = np.ones(len(fileList))
        
        startTime = time.time()
        
        f1 = 0
        #p = 0
        for pathFileInput1 in fileList[0:-1]:
            dirInput1,fileInput1 = os.path.split(pathFileInput1)
            fileName1, fileExtension1 = os.path.splitext(fileInput1)
            fileName1Cut = fileName1[0:-17]
            
            if compFlag == 1:
                patientNameOnly1 = str(re.split('@', fileName1)[-2])
            else:     
                patientNameOnly1 = fileName1         

            ##parse input             
            reads1 = []
            names1 = []
            readCounts1 = []
            source1 = []
            inputHandle = open(pathFileInput1) 
            for record in SeqIO.parse(inputHandle, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
                    reads1.append(str(record.seq))
                    names1.append(record.id)
                    #counts.append(1) 
                    #takes into account the frequencies
                    readCounts1.append(int(re.split('_', record.id)[-1]))
                    source1.append(1) 
            inputHandle.close()
            haploNum1 = len(reads1)
            #print(haploNum1)
            ##mark the start time
            


    
            ## file 2
            f2= f1 + 1
            for pathFileInput2 in fileList[f2:]:
                dirInput2,fileInput2 = os.path.split(pathFileInput2)
                fileName2, fileExtension2 = os.path.splitext(fileInput2)
                fileName2Cut = fileName2[0:-17]
                #print('file 1:', fileName1, 'file 2:', fileName2)
                outputName4 = pathOutput + '/' + fileName1 + '_' + fileName2 + '.fas'            
                outputName5 = pathOutput + '/' + fileName1 + '_' + fileName2 + '_aligned.fas'    
                if compFlag == 1:
                    patientNameOnly2 = str(re.split('@', fileName2)[-2])
                else:     
                    patientNameOnly2 = fileName2     


                    
                if patientNameOnly1 == patientNameOnly2:
                    patientFlag = 1
                else:
                    patientFlag = 0
                
                if patientFlag == 1 or (compFlag == 0 and patientFlag == 0):
                    print('file 1:', fileName1, 'file 2:', fileName2)
                    ##parse input             
                    inputHandle1 = open(pathFileInput2) 
                    reads2 = []
                    readCounts2 = []
                    source2 = []
                    names2 = []
                    for record in SeqIO.parse(inputHandle1, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
                            names1.append(record.id)
                            reads2.append(str(record.seq))
                            #counts.append(1) 
                            #takes into account the frequencies
                            readCounts2.append(int(re.split('_', record.id)[-1])) 
                            source2.append(2) 
                    inputHandle1.close()
                    haploNum2 = len(reads2)



                    haploNum = len(reads1)
                    haploSize = len(reads1[0])
                    ## align pairfile
                    if alignFlag == 1:
                        print('aligning')
                        seqs_aligned = []
                        mafftstr='mafft --auto --thread 20 --quiet --preservecase '+outputName4+' > '+outputName5
                        os.system(mafftstr)
                        seqs_aligned = []
                        inputHandle = open(outputName5)
                        for record in SeqIO.parse(inputHandle, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
                            seqs_aligned.append(record.seq)
                        inputHandle.close()
                        seqs = seqs_aligned
                        seqs1 = seqs[0:haploNum1]    
                        seqs2 = seqs[haploNum1:len(seqs)]    
                    else:
                        seqs1 = reads1    
                        seqs2 = reads2
                    # Calculate distances
                    #finalDist = []
                    #for r1 in range(haploNum-1):
                        #haplotype1 = seqs[r1]
                        #for r2 in range(r1+1, haploNum):
                            #haplotype2 = seqs[r2]
                            #tempDist = 0
                            #for a, b in izip(haplotype1, haplotype2):
                                #if a != b:
                                    #tempDist = tempDist + 1
                            #finalDist.append(tempDist)
                    readNum1 = sum(readCounts1)
                    readNum2 = sum(readCounts2)
                    freq1 = np.divide(readCounts1, readNum1, dtype = float)
                    freq2 = np.divide(readCounts2, readNum2, dtype = float)
                    freqPair = []
                    freqPair1 = []
                    freqPair2 = []
                    h1 = []
                    h2 = []
                    
                   
                    
                    for i in range(haploNum1):
                        for j in range(haploNum2):
                            hf1 = freq1[i]
                            hf2 = freq2[j]
                            h1.append(i)
                            h2.append(j)
                            freqPair.append(hf1*hf2)
                            freqPair1.append(hf1)
                            freqPair2.append(hf2)
                    finalDist = calcDistanceMatrix(seqs2, seqs1)    

                    majorDist = 0
                    major1 = [seqs1[0]]
                    major2 = [seqs2[0]]
                    majorDist = np.divide(sum(1 for a, b in zip(major1, major2) if a != b), haploSize, dtype = float)

                    noZeroFinalDist = []
                    shared = 0
                    sharedFreqList1 = []
                    sharedFreqList2 = []
                    thresholdNum = 0
                    thresholdNumWei = 0
                    h1T = []
                    h2T = []
                    for i in range(len(finalDist)):
                        d = finalDist[i]
                        if d == 0:
                            shared = shared + 1                    
                            noZeroFinalDist.append(0.5)
                            sharedFreqList1.append(freqPair1[i])
                            sharedFreqList2.append(freqPair2[i])
                        if np.divide(d, haploSize, dtype = float) <= np.divide(9, 264, dtype = float):
                            thresholdNum = thresholdNum + 1
                            thresholdNumWei = thresholdNumWei + freqPair[i]
                            h1T.append(h1[i])
                            h2T.append(h2[i])
                        if d > 0:
                            noZeroFinalDist.append(d) 
                    
                    tf1 = 0
                    tf2 = 0
                    for i in range(haploNum1):
                        if i in h1T:
                            tf1 =  tf1 + 1
                    for i in range(haploNum2):
                        if i in h2T:
                            tf2 =  tf2 + 1                            
                    

                    #print (haploNum1,haploNum2,haploNum1*haploNum2,shared,sum(freqPair),len(finalDist))
                    #print(freqPair)
                    
                    #print(finalDist)
                    # Statistics
                    
                    weiDist = np.divide(np.average(finalDist, weights=freqPair), haploSize, dtype = float)
                    avgDist = np.divide(np.mean(finalDist), haploSize, dtype = float)
                    geoDist = np.divide(gmean(noZeroFinalDist), haploSize, dtype = float)
                    minDist = np.divide(min(finalDist), haploSize, dtype = float)
                    maxDist = np.divide(max(finalDist), haploSize, dtype = float)
                    corrFactor = np.divide((withinHamming[f1] + withinHamming[f2]), 2, dtype = float)
                    corrDist = avgDist - corrFactor
                    hapOverlap1 = np.divide(shared, haploNum1, dtype = float)
                    hapOverlap2 = np.divide(shared, haploNum2, dtype = float)
                    hapOverlap = np.mean([hapOverlap1,hapOverlap2])
                    hapOverlapWei1 = sum(sharedFreqList1)
                    hapOverlapWei2 = sum(sharedFreqList2)
                    hapOverlapWei = np.mean([hapOverlapWei1,hapOverlapWei2])
                    thresholdPair = np.divide(thresholdNum, len(finalDist), dtype = float)
                    thresholdFile1 = np.divide(tf1, haploNum1, dtype = float)
                    thresholdFile2 = np.divide(tf2, haploNum2, dtype = float)
                    thresholdFile = max([thresholdFile1,thresholdFile2])
                    #print(thresholdFile1,thresholdFile2)
                    if shared > 1:
                        sharedCorrP = np.corrcoef(sharedFreqList1, sharedFreqList2)[0][1]
                        sharedCorrS = spearmanr(sharedFreqList1, sharedFreqList2)[0]
                    else:
                        sharedCorrP = 0
                        sharedCorrS = 0
                
                    #save .csv fle
                    outputHandle1.write(str(fileName1) + ',' + str(fileName2) + ',' + str(patientFlag) + ',' + str(haploSize) + ',' + str(majorDist) + ',' + str(avgDist) + ',' + str(minDist) + ',' + str(maxDist) + ',' + str(corrDist) + ',' + str(geoDist) + ',' + str(weiDist) + ',' +  str(hapOverlap) + ',' + str(hapOverlapWei) + ','  + str(thresholdPair) + ','  + str(thresholdNumWei) + ','  + str(thresholdFile) + ',' + str(sharedCorrP) + ',' + str(sharedCorrS) + '\n')
                    chosenDist = geoDist
                    if chosenDist == 0:
                        chosenDist = np.divide(0.5, haploSize, dtype = float)
                    outputHandle2.write(str(fileName1) + ',' + str(fileName2) + ',' + str(chosenDist) + '\n')
                

            f1 = f1 + 1
        outputHandle1.close()
        outputHandle2.close()        
                ##mark the end time
        endTime = time.time()
        workTime =  endTime - startTime
        print(workTime)

    except KeyboardInterrupt:
        exit(-1)

