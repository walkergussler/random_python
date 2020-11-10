import numpy as np
import sys, os, re, time, math, glob, optparse, itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from scipy.sparse.csgraph import connected_components, csgraph_from_dense, shortest_path

if __name__ == '__main__':
    try:
        argument_parser = optparse.OptionParser()            
        argument_parser.add_option('-i', metavar='IDIR', action='store', type=str, dest='input_directory', default=os.getcwd(), help='Input file')
        argument_parser.add_option('-o', metavar='ODIR', action='store', type=str, dest='output_directory', default='9step', help='Output directory')
        argument_parser.add_option('-d', action='store', dest='distanceLimit', type=int, default= 9, help='distanceLimit')
        argument_parser.add_option('-t', action='store', dest='allowedCompFraction', type=float, default= 0.05, help='allowedCompFraction')

        options, args = argument_parser.parse_args()    
        pathInput = options.input_directory
        pathOutput = options.output_directory   
        allowedCompFraction = options.allowedCompFraction   
        distanceLimit = options.distanceLimit

        if not os.path.exists(pathOutput):
            os.makedirs(pathOutput)
        
        # save stats
        outputHandle = open(os.path.join(pathOutput, 'big_comp_stats.csv'), "w")
        outputHandle.write('file_name,haploNum,reads,Total_number_components,number_big_components,frequency_threshold,total_big_frequency,workTime \n')

        #sorted file list
        fileList = sorted(glob.glob(os.path.join(pathInput, '*.fas')))

        for pathFileInput in fileList:
        
            seqs = []
            counts = []
            seqNames = []
            timeList = []

            # Mark the start time
            startTime = time.time() 
                    
            # Load sequences    
            dirInput,fileInput = os.path.split(pathFileInput)
            fileName, fileExtension = os.path.splitext(fileInput)
            # raw_input(fileName)

            
            ##parse input & calculate frequences

            input_handle = open(pathFileInput) 
            for record in SeqIO.parse(input_handle, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
                if len(record.seq) > 0 and len(record.seq) < 50000:
                    seqs.append(record.seq)
                    #counts.append(1) 
                    counts.append(int(re.split('_', record.id)[-1]))
                    seqNames.append(record.id) 
                    # takes into account the frequencies

            input_handle.close()

            haploNum = len(seqs)
            haploSize = len(seqs[0])
            #print 'file name =', fileName, '; n of haplotypes =', haploNum, '; size of haplotypes =', haploSize
            nReads = sum(counts)

            ## calculate the ordered frequencies
            
            freqCount = np.zeros((haploSize, 5))
            productVector = np.zeros((haploSize, 5))
            entropyVector = np.zeros(haploSize)
            print fileName
            for read in seqs:
                for pos in range(haploSize):
                    if read[pos] == 'A':
                        freqCount[pos, 0] = freqCount[pos, 0] + 1
                    elif read[pos] == 'C':
                        freqCount[pos, 1] = freqCount[pos, 1] + 1
                    elif read[pos] == 'G':
                        freqCount[pos, 2] = freqCount[pos, 2] + 1
                    elif read[pos] == 'T':
                        freqCount[pos, 3] = freqCount[pos, 3] + 1
                    elif read[pos] == '-':
                        freqCount[pos, 4] = freqCount[pos, 4] + 1
            freqRel = np.divide(freqCount, haploNum, dtype = float)
            for pos in range(haploSize):
                for i in range(5):
                    freqPos = freqRel[pos, i]
                    if freqPos > 0:
                        logFreqRel = math.log(freqPos, 2)
                        productVector[pos, i] = -1*(np.multiply(freqPos, logFreqRel, dtype = float))
                        
            hVector = np.sum(productVector, axis = 1)
            # order positions
            invH =  np.multiply(-1, hVector, dtype = float)
            ordH = np.argsort(invH)
        
            # reorder the sequences by entropy
            ordSeqs = []
            for i in seqs:
                newOne = ''
                for p in range(haploSize):
                    newOne = ''.join([newOne, i[ordH[p]]])
                ordSeqs.append(newOne)
    
            # Calculate distances
            compNum = haploSize
            compList = range(haploNum)
            adjMatrix = np.zeros((haploNum, haploNum))

            # check sequence pair below distLimit
            for r1 in range(haploNum-1):
                haplotype1 = ordSeqs[r1]
                for r2 in range(r1+1, haploNum):
                    haplotype2 = ordSeqs[r2]
                    tempDist = 0
                    for a, b in itertools.izip(haplotype1, haplotype2):
                        if a != b:
                            tempDist = tempDist + 1
                            if tempDist > distanceLimit:
                                break
                    if tempDist <= distanceLimit: 
                        adjMatrix[r1][r2] = 1
                        adjMatrix[r2][r1] = 1
            # calculate components
            sparseMatrix = csgraph_from_dense(adjMatrix)
            connected = connected_components(sparseMatrix, directed=False, connection='weak', return_labels=True)
            compNum = connected[0]
            compList = connected[1]
            haploNumList = range(haploNum)
            bigCompNum = 0
            bigFreq = 0
            for i in range(compNum):
                freqSub = []
                indices = []
                idx = 0
                for comp in compList:
                    if comp == i:
                        indices.append(haploNumList[idx])
                        freqSub.append(counts[idx])
                    idx = idx+1
                tempFreqCount = np.sum(freqSub)
                tempFreqFraction = np.divide(tempFreqCount, nReads, dtype = float)
                if tempFreqFraction >= allowedCompFraction:
                    bigFreq = bigFreq + tempFreqFraction
                    bigCompNum = bigCompNum + 1
                    # outputName1 = pathOutput + '/' + fileName + '_' + str(i) + '_links.csv'
                    outputName2 = pathOutput + '/' + fileName + '_' + str(i) + '.fas'
                    compSize = len(indices)
                    
                    # adjMatrixComp = np.zeros((compSize, compSize))
                    # calculate distance between every pair of sequences    
                    for s1 in range(compSize-1):
                        haplotype1 = seqs[indices[s1]]
                        for s2 in range(s1+1, compSize):
                            haplotype2 = seqs[indices[s2]]
                            tempDist = 0
                            for a, b in itertools.izip(haplotype1, haplotype2):
                                if a != b:
                                    tempDist = tempDist + 1
                            # if tempDist <= distanceLimit:
                                # adjMatrixComp[s1][s2] = 1
                                # adjMatrixComp[s2][s1] = 1
                    # sparseMatrixComp = csgraph_from_dense(adjMatrixComp)
                    # pathDistMat = shortest_path(sparseMatrixComp, method='auto', directed=False, return_predecessors=False, unweighted=True, overwrite=False)   
                    
                    #save .fas file
                    with open(outputName2, "w") as outputHandle2:
                        n = 0
                        for k in indices:
                            haplotype = str(seqs[k])
                            n += 1
                            seq_id = str(n) + '_' + str(counts[k])
                            seq_dna = Seq(haplotype,generic_dna)
                            seqFASTA = SeqRecord(seq_dna, id = seq_id, name = "", description = "")     
                            SeqIO.write(seqFASTA, outputHandle2, "fasta")   
                        
            print(fileName)
            outputHandle.write(fileName + ',' + str(haploNum) + ',' + str(nReads) + ',' + str(compNum) + ',' + str(bigCompNum) + ',' + str(allowedCompFraction) + ',' + str(bigFreq) + '\n')
        outputHandle.close()



    except KeyboardInterrupt:
        exit(-1)

