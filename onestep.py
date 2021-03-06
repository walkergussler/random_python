import os, re, glob, sys
from optparse import OptionParser
from math import log
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from scipy.sparse.csgraph import connected_components, csgraph_from_dense, shortest_path
from tempfile import NamedTemporaryFile
from subprocess import check_call
from collections import defaultdict

def align(seqs):
    with NamedTemporaryFile(delete=False, mode='w') as seqdump:
        catname=seqdump.name
        for seq in seqs:
            seqdump.write('>seq_'+str(seqs[seq])+'\n'+str(seq)+'\n')
    with NamedTemporaryFile(delete=False) as aligned:
        alignname=aligned.name
        check_call(['mafft', '--quiet', '--auto', '--thread', '20', '--preservecase', catname], stdout=aligned)
    os.unlink(catname)
    seqs,_=parse_input(alignname)
    if type(seqs)==bool:
        print('error parsing aligned seqs!')
        sys.exit()
    return seqs

def parse_input(input): #get sequences from a file
    seqs=defaultdict(int)
    with open(input,'r') as f:
        for record in SeqIO.parse(f,'fasta'):
            freq = int(record.id.split('_')[-1])
            seq=record.seq.upper()
            seqs[seq]+=freq
    return seqs

if __name__ == '__main__':
    try:
        argument_parser = OptionParser()            
        argument_parser.add_option('-i', metavar='IDIR', action='store', type=str, dest='input_directory', default=os.getcwd(), help='Input file')
        argument_parser.add_option('-o', metavar='ODIR', action='store', type=str, dest='output_directory', default='one_step', help='Output directory')
        argument_parser.add_option('-d', action='store', dest='distanceLimit', type=int, default= 1, help='distanceLimit')
        argument_parser.add_option('-t', action='store', dest='allowedCompFraction', type=float, default= 0.05, help='allowedCompFraction')
        argument_parser.add_option('-s', action='store', dest='saveBigOnly', type=int, default= 1, help='save only biggest')
        options, args = argument_parser.parse_args()    
        pathInput = options.input_directory
        pathOutput = options.output_directory    
        allowedCompFraction = options.allowedCompFraction    
        distanceLimit = options.distanceLimit
        saveBigOnly = options.saveBigOnly
        if not os.path.exists(pathOutput):
            os.makedirs(pathOutput)
        #sorted file list
        fileList = sorted(glob.glob(os.path.join(pathInput, '*.fa*')))
        finished=[]
        with open('onestep_trash.txt') as f:
            for line in f.readlines():
                finished.append(line.strip())
        files=[]
        for file in fileList:
            if file not in finished:
                files.append(file)
            else:
                print(file)
        # sys.exit()
        for pathFileInput in fileList:
            print("running"+pathFileInput)
            seqs = []
            counts = []
            seqNames = []
            # Load sequences    
            _,fileInput = os.path.split(pathFileInput)
            fileName, _ = os.path.splitext(fileInput)
            ##parse input & calculate frequences
            seqsDict=align(pathFileInput)
            seqs=list(seqsDict.keys())
            counts=list(seqsDict.values())
            haploNum = len(seqs)
            haploSize = len(seqs[0])
            nReads = sum(counts)
            ## calculate the ordered frequencies
            freqCount = np.zeros((haploSize, 5))
            productVector = np.zeros((haploSize, 5))
            entropyVector = np.zeros(haploSize)
            print(fileName)
            readCounter = 0
            stateList = ['A', 'C', 'G', 'T', '-']
            for read in seqs:
                for pos in range(haploSize):
                    print(pos,haploSize)
                    if read[pos] == stateList[0]:
                        freqCount[pos, 0] = freqCount[pos, 0] + counts[readCounter]
                    elif read[pos] == stateList[1]:
                        freqCount[pos, 1] = freqCount[pos, 1] + counts[readCounter]
                    elif read[pos] == stateList[2]:
                        freqCount[pos, 2] = freqCount[pos, 2] + counts[readCounter]
                    elif read[pos] == stateList[3]:
                        freqCount[pos, 3] = freqCount[pos, 3] + counts[readCounter]
                    elif read[pos] == stateList[4]:
                        freqCount[pos, 4] = freqCount[pos, 4] + counts[readCounter]
            freqRel = np.divide(freqCount, haploNum, dtype = float)
            for pos in range(haploSize):
                for i in range(5):
                    freqPos = freqRel[pos, i]
                    if freqPos > 0:
                        logFreqRel = log(freqPos, 2)
                        productVector[pos, i] = -1*(np.multiply(freqPos, logFreqRel, dtype = float))
            hVector = np.sum(productVector, axis = 1)
            # order positions
            invH =  np.multiply(-1, hVector, dtype = float)
            ordH = np.argsort(invH)
            # reorder the sequences by entropy
            ordSeqs = []
            for i in seqs:
                newOne = []
                for p in range(haploSize):
                    newOne.append(str(i[ordH[p]]))
                ordSeqs.append(''.join(newOne))
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
                    for a, b in zip(haplotype1, haplotype2):
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
            ## component sizes
            finalCompFreqList = []
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
                finalCompFreqList.append(np.divide(tempFreqCount, nReads, dtype = float))
            ## find biggest component
            bigIdxFinal = np.argmax(finalCompFreqList)
            if saveBigOnly == 0:
                bigFlagList = np.ones(compNum)
            else:
                bigFlagList = np.zeros(compNum)
                for i in range(compNum):
                    if i == bigIdxFinal:
                        bigFlagList[i] = 1
            ## results for each component
            for i in range(compNum):
                tempFreqFraction = finalCompFreqList[i]
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
                if tempFreqFraction >= allowedCompFraction and bigFlagList[i] == 1:
                    bigFreq = bigFreq + tempFreqFraction
                    bigCompNum = bigCompNum + 1
                    outputName2 = pathOutput + '/' + fileName + '_onestep.fas'
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
    except KeyboardInterrupt:
        exit(-1)