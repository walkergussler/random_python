#!/usr/bin/env python

__author__ = 'David S. Campo, Ph.D.'
'''
entropy_dist. Program for calculating the pairs of sequences between two files that are below the relatednes threshold.
1. The entropy of each position in file 1 is calcualted
2. The sequences are re-ordered according to this value, so higher entropy positions go first.
3. For every pair of sequences, hamming distances are calculated, but now it stops when the distance passes the relatedness threshold.


Usage:
module load Python/2.7.3

With default parameters:
python entropy_dist.py -i database_folder_location -q query_folder_location -o output_folder_location

Changing the threshold:
python entropy_dist.py -i database_folder_location -q query_folder_location -o output_folder_location -t 9

input Folders must contain Fasta files, aligned to each other.
'''

import time
import sys, os, re
import numpy as np
import math
import glob
import optparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import itertools
from string import maketrans
from multiprocessing import Pool
from collections import Counter, OrderedDict
import operator
import functools
from scipy.stats import binom
from random import randrange
from collections import defaultdict
from itertools import izip


if __name__ == '__main__':
	try:
        	argument_parser = optparse.OptionParser()            
        	argument_parser.add_option('-i',
            		metavar='IDIR', action='store', type=str, dest='input_directory_database', default='.',
            		help='Input directory')
        	argument_parser.add_option('-q',
            		metavar='IDIR', action='store', type=str, dest='input_directory_query', default='.',
            		help='Input directory')
        	argument_parser.add_option('-o',
            		metavar='ODIR', action='store', type=str, dest='output_directory', default='.',
            		help='Output directory')  
        	argument_parser.add_option('-t',
            		action='store', dest='t', type= int, default= 9, help='threshold')  
        	options, args = argument_parser.parse_args()
		
			
		pathInput = options.input_directory_database
		pathQuery = options.input_directory_query
		pathOutput = options.output_directory	
		t = options.t

		
		if not os.path.exists(pathOutput):
			os.makedirs(pathOutput)
		
		seqs = []
		counts = []
		seqNames = []
		timeList = []
		#print 'reading database'
		
			
		for pathFileInput in glob.glob(os.path.join(pathInput, '*.fas')):
			
			dirInput,fileInput = os.path.split(pathFileInput)
			fileName, fileExtension = os.path.splitext(fileInput)
			
			##parse input & calculate frequences

			input_handle = open(pathFileInput) 
			for record in SeqIO.parse(input_handle, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
				seqs.append(record.seq)
				#counts.append(1) 
				counts.append(int(re.split('_', record.id)[-1]))
				seqNames.append(record.id) 
				# takes into account the frequencies
			
			input_handle.close()
		haploNum = len(seqs)
		haploSize = len(seqs[0])
		
			
		## calculate the ordered frequencies
		
		freqCount = np.zeros((haploSize, 5))
		productVector = np.zeros((haploSize, 5))
		entropyVector = np.zeros(haploSize)
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
		
		finalDist = np.zeros(haploNum) + (t + 1)
		
		for i in range(10):
			## bring the queries			
			#print 'reading query'
			for pathFileQuery in glob.glob(os.path.join(pathQuery, '*.fas')):
				dirQuery,fileQuery = os.path.split(pathFileQuery)
				fileNameQuery, fileExtensionQuery = os.path.splitext(fileQuery)
				#print 'query name', fileNameQuery
					
				outputHandle1 = open(os.path.join(pathOutput, fileNameQuery + '_multi_hash.txt'), "w")
				outputHandle2 = open(os.path.join(pathOutput, fileNameQuery + '_time.txt'), "w")
				
				##parse input & calculate frequences
				input_handle = open(pathFileQuery) 
				seqsQ = []
				for record in SeqIO.parse(input_handle, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
					seqsQ.append(record.seq)
					#counts.append(1) 
					counts.append(int(re.split('_', record.id)[-1])) 
					# takes into account the frequencies
				input_handle.close()
				haploNumQuery = len(seqsQ)
				haploSizeQuery = len(seqs[0])
	
				# reorder the sequences by entropy
				ordSeqsQ = []
				for i in seqsQ:
					newOneQ = ''
					for p in range(haploSizeQuery):
						newOneQ = ''.join([newOneQ, i[ordH[p]]])
					ordSeqsQ.append(newOne)
	
				# mark the start time
				startTime = time.time()
	
				# Check each query sequence	
				candidList = np.zeros(haploNum)
				for r1 in range(haploNum):
					if candidList[r1] == 0:
						haplotype1 = ordSeqs[r1]
						for r2 in range(haploNumQuery):
							haplotype2 = ordSeqsQ[r2]
							tempDist = 0
							for a, b in izip(haplotype1, haplotype2):
								if a != b:
									tempDist = tempDist + 1
									if tempDist > t:
										candidList[r1] = 1
										break
							finalDist[r1] = tempDist
				##mark the end time
				endTime = time.time()
				workTime =  endTime - startTime
				
				#print 'saving ouput files'	
				#save names of sequences below threshold
				l = 0
				passed = 0
				for fd in finalDist:
					if candidList[l] == 1:  
						passed = passed + 1
						#outputHandle1.write(str(seqNames[l]) + '\t' + str(fd) + '\n')
					l = l + 1	
				outputHandle1.close()
				
			
				##mark the end time
				#endTime = time.time()
				#workTime =  endTime - startTime
				timeList.append(workTime)
				#print 'seconds =', workTime
				#print 'sequences passed =', passed
				
				outputHandle2.write(str(fileNameQuery) + '\t' + str(workTime) + '\n')
				outputHandle2.close()
		print np.mean(timeList)	
	except KeyboardInterrupt:
		exit(-1)

