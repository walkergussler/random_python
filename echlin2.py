#!/usr/bin/env python

__author__ = 'David S. Campo, Ph.D.'
'''
one_step_stats. Program for calculating network statistics of one-step networks. 
1.Input is an edge list in csv format. The list was produced by one_step_definition.py and includes: node1, node2, hamming distance, one-step logical, path distances, frequency of node 1 and frequency of node 2.  
2. Output is .csv file with the summary of several topological properties of each network

usage:
python one_step_stats.py -i input_folder_location -o output_folder_location


'''
#import pandas as pd
import os, optparse, csv
import numpy as np
import math, glob
import networkx as nx
#from collections import defaultdict, Counter, OrderedDict
from itertools import combinations
#from itertools import izip, combinations, permutations
from scipy.sparse.csgraph import csgraph_from_dense, shortest_path
from scipy.stats import linregress
from scipy.stats import variation
from sklearn.metrics import mean_squared_error
                                       

def motifDictCalc():
    
    global motifDict
    motifDict = {}
    ##generate all possible graph4 motifs
    #tempAdj = np.zeros((4,4))
    for e1 in range(2):
        for e2 in range(2):
            for e3 in range(2):
                for e4 in range(2):
                    for e5 in range(2):
                        for e6 in range(2):
                            #tempAdj[0][1] = e1
                            #tempAdj[0][2] = e2
                            #tempAdj[0][3] = e3
                            #tempAdj[1][2] = e4
                            #tempAdj[1][3] = e5
                            #tempAdj[2][3] = e6
                            sumEdges = e1+e2+e3+e4+e5+e6    
                            degreeList = [e1+e2+e3,e1+e4+e5,e2+e4+e6,e3+e5+e6]
                            maxSubgraphK = max(degreeList)
                            minSubgraphK = min(degreeList)
                            Q = (e1,e2,e3,e4,e5,e6)
                            if sumEdges < 3 or minSubgraphK < 1:
                                motifDict[Q] = 7
                            else:
                                if sumEdges == 3:
                                    if maxSubgraphK == 3:# and minSubgraphK == 1:# and subgraph_0_1 == 1  and subgraph_0_2 == 1 and subgraph_0_3 == 1:
                                        motifDict[Q] = 1
                                    elif maxSubgraphK == 2:# and minSubgraphK == 1:# and minSubgraphK == 1:# and subgraph_0_1 == 1 and subgraph_0_2 == 1 and subgraph_1_3 == 1:
                                        motifDict[Q] = 2
                                elif sumEdges == 4:
                                    if maxSubgraphK == 2:# and minSubgraphK == 2:# and subgraph_0_1 == 1  and subgraph_1_2 == 1 and subgraph_0_3 == 1 and subgraph_2_3 == 1:
                                        motifDict[Q] = 3
                                    elif maxSubgraphK == 3:# and minSubgraphK == 1:# and subgraph_0_1 == 1  and subgraph_1_2 == 1 and subgraph_0_3 == 1 and subgraph_0_2 == 1:
                                        motifDict[Q] = 4
                                elif sumEdges == 5:# and maxSubgraphK == 3:# and minSubgraphK == 2:# and minSubgraphK == 2:# and subgraph_0_1 == 1  and subgraph_1_2 == 1 and subgraph_0_3 == 1 and subgraph_0_2 == 1 and subgraph_1_3 == 1:
                                    motifDict[Q] = 5
                                elif sumEdges == 6:# and maxSubgraphK == 3:# and minSubgraphK == 3:# and subgraph_0_1 == 1  and subgraph_1_2 == 1 and subgraph_2_3 == 1 and subgraph_0_3 == 1 and subgraph_1_3 == 1 and subgraph_0_2 == 1:
                                    motifDict[Q] = 6
    return motifDict

def graphEntropyCalc(colors): 
    colorList = list(set(list(colors.values())))
    nodeList = list(set(list(colors.keys())))
    nodeNum = len(nodeList)
    p = []
    equalFreq = np.divide(1, nodeNum, dtype= float)
    for i in range(nodeNum):
        p.append(equalFreq)
    colorFreq = []
    for j in colorList:
        colorTemp = 0
        for i in nodeList:
            if colors[i] == j:
                colorTemp = colorTemp + 1
        colorFreq.append(np.divide(colorTemp, nodeNum, dtype = float))
    colorEntropy = []
    for j in colorFreq:
        hTemp = []
        for i in p:
            hTemp.append(i*math.log(np.divide(1, j, dtype = float),2))
            
        colorEntropy.append(sum(hTemp))
    
    graphEntropy = min(colorEntropy)
    return graphEntropy

def motifCalc4(G):
    motif1Star, motif2Path, motif3Cycle, motif4TailedTriangle, motif5Envelope, motif6Clique = 0,0,0,0,0,0
    motifTempList = []
    nodeNames = list(GM.nodes())
    nodeNum = len(nodeNames)
    for n0 in range(nodeNum-3):
        for n1 in range(n0+1, nodeNum-2):
            if G4.has_edge(nodeNames[n0], nodeNames[n1]):
                e1 = GM.has_edge(nodeNames[n0], nodeNames[n1])
                for n2 in range(n1+1, nodeNum-1):
                    if G4.has_edge(nodeNames[n0], nodeNames[n2]) and G4.has_edge(nodeNames[n1], nodeNames[n2]):
                        e2 = GM.has_edge(nodeNames[n0], nodeNames[n2])
                        for n3 in range(n2 + 1, nodeNum):
                            if G4.has_edge(nodeNames[n0], nodeNames[n3]):# and G4.has_edge(n1, n3) and G4.has_edge(n2, n3):
                                e3 = GM.has_edge(nodeNames[n0], nodeNames[n3])
                                if (e1+e2+e3) > 0:
                                    Q = (e1,e2,e3,GM.has_edge(nodeNames[n1], nodeNames[n2]),GM.has_edge(nodeNames[n1], nodeNames[n3]),GM.has_edge(nodeNames[n2], nodeNames[n3]))
                                    motifTempList.append(motifDict[Q])    

    for i in motifTempList:
        if i == 1:
            motif1Star = motif1Star + 1
        elif i == 2:
            motif2Path = motif2Path + 1
        elif i == 3:
            motif3Cycle = motif3Cycle + 1
        elif i == 4:
            motif4TailedTriangle = motif4TailedTriangle + 1
        elif i == 5:
            motif5Envelope = motif5Envelope + 1
        elif i == 6:
            motif6Clique = motif6Clique + 1

    nMotifs = motif1Star + motif2Path + motif3Cycle + motif4TailedTriangle + motif5Envelope  + motif6Clique
    #print motif1Star, motif2Path, motif3Cycle, motif4TailedTriangle, motif5Envelope ,motif6Clique
    motif1StarFrac=np.divide(motif1Star, nMotifs, dtype = float)
    motif2PathFrac=np.divide(motif2Path, nMotifs, dtype = float)
    motif3CycleFrac=np.divide(motif3Cycle, nMotifs, dtype = float)
    motif4TailedTriangleFrac=np.divide(motif4TailedTriangle, nMotifs, dtype = float)
    motif5EnvelopeFrac=np.divide(motif5Envelope, nMotifs, dtype = float)
    motif6CliqueFrac=np.divide(motif6Clique, nMotifs, dtype = float)
    freqM = [motif1StarFrac, motif2PathFrac, motif3CycleFrac, motif4TailedTriangleFrac, motif5EnvelopeFrac, motif6CliqueFrac]
    return freqM

def entropyCalc(freqM):#different results than scipy.stats.entropy - how to resolve?
    productVectorM = 0
    for i in freqM:
        if i > 0:
            productVectorM = productVectorM + (i*math.log(i, 2))
    entropy = -1*productVectorM
    return entropy

def localEfficiencyCalc(GL):
    nodeNum = len(GL)
    localEfficiency = 0
    if nodeNum > 1:
        local1 = []
        for boxName in nx.nodes_iter(GL):
            radiusNodeList = GL.neighbors(boxName)
            boxNet = nx.Graph(GL.subgraph(radiusNodeList))
            boxNodes = len(boxNet)
            boxMat = nx.to_numpy_matrix(boxNet)
            boxSparse = csgraph_from_dense(boxMat)
            boxMatPath = shortest_path(boxSparse, method='auto', directed=False, return_predecessors=False, unweighted=True, overwrite=False)    
            boxPathList = []
            for i in range(boxNodes-1):
                for j in range(i+1, boxNodes):
                    tempDist = boxMatPath[i][j]
                    if np.isfinite(tempDist):
                        boxPathList.append(np.divide(1, tempDist, dtype = float))
            if len(boxPathList) > 0:
                local1.append(np.mean(boxPathList))
            else:
                local1.append(0)    
        localEfficiency = np.mean(local1)    
    
    return localEfficiency                    

def boxStats(boxNet): #fordavid other three calculated here?
    ## matrices    
    boxNodes = len(boxNet)
    boxMat = nx.to_numpy_matrix(boxNet)
    boxSparse = csgraph_from_dense(boxMat)
    boxMatPath = shortest_path(boxSparse, method='auto', directed=False, return_predecessors=False, unweighted=True, overwrite=False)    
    boxPathList = []
    pairsNumBox = len(list(combinations(range(boxNodes), 2)))
    for i in range(boxNodes-1):
        for j in range(i+1, boxNodes):
            tempDist = boxMatPath[i][j]
            if tempDist > 0 and np.isfinite(tempDist):
                boxPathList.append(tempDist)
    
    ##boxNet characteristics
    degreeRaw = list(boxNet.degree())
    degreeBox = []
    for i in degreeRaw:
        degreeBox.append(i)
    degreeNormBox = np.divide(degreeBox, np.sum(degreeBox), dtype = float)
    
    diameterPathBox = np.max(boxPathList)
    avgPathDistBox = np.mean(boxPathList)
    nEdgesBox = np.divide(np.sum(degreeBox), 2, dtype = float)
    edgePBox = nx.density(boxNet)
    globalEfficiencyBox = np.divide(sum(np.divide(1, boxPathList, dtype = float)),pairsNumBox , dtype = float)
    radiusBox = nx.radius(boxNet)
    kCoreBox = max(list(nx.core_number(boxNet).values()))
    degreeAssortBox = nx.degree_assortativity_coefficient(boxNet)
    avgDegreeBox = np.mean(degreeBox)
    maxDegreeBox = max(degreeBox)
    eValsBox = np.linalg.eigvals(boxMat)
    spectralRadiusAdjBox = max(abs(eValsBox))
    eigenCentDictBox = nx.eigenvector_centrality_numpy(boxNet, weight=None)
    eigenCentRawBox = list(eigenCentDictBox.values())
    eigenCentBox = np.divide(eigenCentRawBox, sum(eigenCentRawBox), dtype = float)
    colorsBox = nx.coloring.greedy_color(boxNet, strategy=nx.coloring.strategy_connected_sequential_bfs)
    colorNumBox = len(list(set(list(colorsBox.values()))))
    avgClustCoeffBox = nx.average_clustering(boxNet)                        
    scaledSpectralRadiusBox = np.divide(spectralRadiusAdjBox, avgDegreeBox, dtype = float)
    if motifChoice == 1:
        freqMBox = motifCalc4(boxNet)
    else:
        freqMBox =  [0.166666667, 0.166666667, 0.166666667, 0.166666667, 0.166666667, 0.166666667]
    # network entropy
    lapMatBox= np.asarray(nx.to_numpy_matrix(nx.from_scipy_sparse_matrix(nx.laplacian_matrix(boxNet))))
    eValsLapBox = np.linalg.eigvals(lapMatBox)
    eValsLapBoxSorted = sorted(np.real(eValsLapBox))
    spectralGapBox = eValsLapBoxSorted[1]
    degreeSumBox = np.sum(degreeBox)
    lapMatBoxNorm =  np.divide(lapMatBox, degreeSumBox, dtype = float)
    eValsLapBoxNorm = np.linalg.eigvals(lapMatBoxNorm)
    eValsLapNonZeroBoxNorm = []
    for i in eValsLapBoxNorm:
        j = abs(i)
        if j > 0:
            eValsLapNonZeroBoxNorm.append(j)
    vonEntropyBox = np.divide(entropyCalc(eValsLapNonZeroBoxNorm), math.log(boxNodes,2), dtype = float)
    degreeEntropyBox = np.divide(entropyCalc(degreeNormBox), math.log(boxNodes,2), dtype = float)
    KSEntropyBox = np.divide(math.log(spectralRadiusAdjBox, 2), math.log(boxNodes-1,2), dtype = float)
    motifEntropyBox = np.divide(entropyCalc(freqMBox), math.log(len(freqMBox),2), dtype = float)
    popEntropyBox = np.divide(entropyCalc(eigenCentBox), math.log(boxNodes,2), dtype = float)
    graphEntropyBox = np.divide(graphEntropyCalc(colorsBox), math.log(boxNodes,2), dtype = float)
    
    return edgePBox, radiusBox, kCoreBox, degreeAssortBox, diameterPathBox, avgPathDistBox, nEdgesBox, globalEfficiencyBox, avgDegreeBox, maxDegreeBox, spectralRadiusAdjBox, spectralGapBox, scaledSpectralRadiusBox, colorNumBox, avgClustCoeffBox, freqMBox, motifEntropyBox, vonEntropyBox, graphEntropyBox, popEntropyBox, KSEntropyBox, degreeEntropyBox

def comparisonList(edgePBox, radiusBox,  kCoreBox,  degreeAssortBox, diameterPathBox, avgPathDistBox, nEdgesBox, globalEfficiencyBox, avgDegreeBox, maxDegreeBox, spectralRadiusAdjBox, spectralGapBox, scaledSpectralRadiusBox, colorNumBox, avgClustCoeffBox, freqMBox, motifEntropyBox, vonEntropyBox, graphEntropyBox, popEntropyBox, KSEntropyBox, degreeEntropyBox):
    # list of comparisons
    diameterPathBoxList.append(np.absolute(diameterPath - diameterPathBox)) 
    avgPathDistBoxList.append(np.absolute(avgPathDist - avgPathDistBox)) 
    nEdgesBoxList.append(np.absolute(nEdges - nEdgesBox)) 
    edgePBoxList.append(np.absolute(edgeP - edgePBox))
    radiusBoxList.append(np.absolute(radius - radiusBox))  
    kCoreBoxList.append(np.absolute(kCore - kCoreBox))  
    degreeAssortBoxList.append(np.absolute(degreeAssort - degreeAssortBox))  
    globalEfficiencyBoxList.append(np.absolute(globalEfficiency - globalEfficiencyBox)) 
    avgDegreeBoxList.append(np.absolute(avgDegree - avgDegreeBox)) 
    maxDegreeBoxList.append(np.absolute(maxDegree - maxDegreeBox)) 
    spectralRadiusAdjBoxList.append(np.absolute(spectralRadiusAdj - spectralRadiusAdjBox)) 
    spectralGapBoxList.append(np.absolute(spectralGap - spectralGapBox))  
    scaledSpectralRadiusBoxList.append(np.absolute(scaledSpectralRadius - scaledSpectralRadiusBox)) 
    colorNumBoxList.append(np.absolute(colorNum - colorNumBox)) 
    avgClustCoeffBoxList.append(np.absolute(avgClustCoeff - avgClustCoeffBox))                     
    freqMBoxListEach = []
    for i in range(6):
        freqMBoxListEach.append(np.absolute(freqM[i] - freqMBox[i])) 
    freqMBoxList.append(np.mean(freqMBoxListEach))
    motifEntropyBoxList.append(np.absolute(motifEntropy - motifEntropyBox))                                 
    vonEntropyBoxList.append(np.absolute(vonEntropy - vonEntropyBox)) 
    graphEntropyBoxList.append(np.absolute(graphEntropy - graphEntropyBox))                
    popEntropyBoxList.append(np.absolute(popEntropy - popEntropyBox)) 
    KSEntropyBoxList.append(np.absolute(KSEntropy - KSEntropyBox))
    degreeEntropyBoxList.append(np.absolute(degreeEntropy - degreeEntropyBox))

def fractalCalc(dist, nodeNum, methodChoice):
    pairsNum = len(dist)
    diameter = max(dist)
    ## if only one sequence
    if nodeNum < 2:
        sumBoxes, hammDb, hammRSquare, hammDbConstant, hammModularity, hammModularityConstant, hammModularityRSquare, diameterPathSelf, avgPathDistSelf, nEdgesSelf, globalEfficiencySelf, avgDegreeSelf, maxDegreeSelf, spectralRadiusAdjSelf, spectralGapSelf, popEntropySelf, scaledSpectralRadiusSelf, colorNumSelf, avgClustCoeffSelf, freqMBoxSelf, motifEntropySelf, graphEntropySelf = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    else:
        lb = 1
        nBoxesAll = []
        modLbAll = []
        # initial box size 0
        nBoxes = nodeNum
        nBoxesAll.append(nBoxes)
        modLb = 0
        boxWeightList = []
        
        ## self similarity lists
        global radiusBoxList, kCoreBoxList, degreeAssortBoxList, diameterPathBoxList, avgPathDistBoxList, nEdgesBoxList, globalEfficiencyBoxList, avgDegreeBoxList, maxDegreeBoxList, spectralRadiusAdjBoxList, spectralGapBoxList, colorNumBoxList, avgClustCoeffBoxList, freqMBoxList, motifEntropyBoxList, graphEntropyBoxList, scaledSpectralRadiusBoxList, edgePBoxList, popEntropyBoxList, vonEntropyBoxList, KSEntropyBoxList, degreeEntropyBoxList
        radiusBoxList, kCoreBoxList, degreeAssortBoxList, diameterPathBoxList, avgPathDistBoxList, nEdgesBoxList, globalEfficiencyBoxList, avgDegreeBoxList, maxDegreeBoxList, spectralRadiusAdjBoxList, spectralGapBoxList, colorNumBoxList, avgClustCoeffBoxList,     freqMBoxList, motifEntropyBoxList, graphEntropyBoxList, scaledSpectralRadiusBoxList, edgePBoxList, popEntropyBoxList, vonEntropyBoxList, KSEntropyBoxList, degreeEntropyBoxList = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []


        # random nodeList for growth method (it has to be outside the loop because it is done just once)
        numChosen = 100
        nodeListRandom = []
        for i in range(numChosen):
            nodeListRandom.append(np.random.choice(nodeList))
        
        while lb < diameter:
            lb = lb +1
            if lb not in dist:
                nBoxesAll.append(nBoxes)
                modLbAll.append(modLb)
            else:
                
                # make new M  graph
                edgeListStep = []
                for i in range(pairsNum):
                    if dist[i] >= lb:
                        edgeListStep.append([node1[i], node2[i]])
                M = nx.Graph()
                M.add_nodes_from(nodeList)
                M.add_edges_from(edgeListStep)
                
                # coloring
                boxes = nx.coloring.greedy_color(M, strategy=nx.coloring.strategy_saturation_largest_first)
                boxesList = list(set(list(boxes.values())))
                nBoxes = len(boxesList)
                nBoxesAll.append(nBoxes)
                withinBox = 1
                betweenBox = 1

                if methodChoice == 0:
                    # box network and box Modularity
                    allBoxesDict = {}                
                    for boxName in boxesList:
                        allBoxesDict[boxName] = nx.Graph()
                    for i in range(pairsNum):
                        if dist[i] == 1:
                            if boxes[node1[i]] == boxes[node2[i]]:
                                withinBox = withinBox + 1
                                allBoxesDict[boxes[node1[i]]].add_edge(node1[i], node2[i])
                            else:
                                betweenBox = betweenBox + 1    
                    modLb = np.divide(np.divide(withinBox, betweenBox, dtype = float), nBoxes, dtype = float)
                    modLbAll.append(modLb)
                    ## stats
                    for boxName in boxesList:
                        boxNet = allBoxesDict[boxName]
                        boxNodes = len(boxNet)
                        if  boxNodes > 4 and nx.is_connected(boxNet):
                            boxWeight = np.divide(boxNodes, nodeNum, dtype = float)
                            boxWeightList.append(boxWeight)
                            edgePBox, radiusBox, kCoreBox, degreeAssortBox, diameterPathBox, avgPathDistBox, nEdgesBox, globalEfficiencyBox, avgDegreeBox, maxDegreeBox, spectralRadiusAdjBox, spectralGapBox, scaledSpectralRadiusBox, colorNumBox, avgClustCoeffBox, freqMBox, motifEntropyBox, vonEntropyBox, graphEntropyBox, popEntropyBox, KSEntropyBox, degreeEntropyBox = boxStats(boxNet)
                            #print freqMBox
                            
                            # write to box file
                            outputHandleBox.write(str(fileName1) + ',' + str(lb) + ',' + str(nBoxes) + ',' + str(boxName) + ',' + str(boxNodes) + ',' + str(diameterPathBox) + ',' + str(avgPathDistBox) + ',' + str(nEdgesBox) + ',' + str(edgePBox) + ',' + str(radiusBox) + ',' + str(kCoreBox) + ',' + str(degreeAssortBox) + ',' + str(globalEfficiencyBox) + ',' + str(avgDegreeBox) + ',' + str(maxDegreeBox) + ',' + str(spectralRadiusAdjBox) + ',' + str(spectralGapBox) + ',' + str(popEntropyBox) + ',' + str(colorNumBox) + ',' + str(avgClustCoeffBox) + ',' + str(scaledSpectralRadiusBox) + ',' + str(vonEntropyBox) + ',' + str(KSEntropyBox) + ',' + str(degreeEntropyBox) + ',' + str(graphEntropyBox) + ',' + str(motifEntropyBox) + ',' + str(100*freqMBox[0]) + ',' + str(100*freqMBox[1]) + ',' + str(100*freqMBox[2]) + ',' + str(100*freqMBox[3]) + ',' + str(100*freqMBox[4]) + ',' + str(100*freqMBox[5])     + '\n')

                            # list of comparisons
                            comparisonList(edgePBox, radiusBox,  kCoreBox,  degreeAssortBox, diameterPathBox, avgPathDistBox, nEdgesBox, globalEfficiencyBox, avgDegreeBox, maxDegreeBox, spectralRadiusAdjBox, spectralGapBox, scaledSpectralRadiusBox, colorNumBox, avgClustCoeffBox, freqMBox, motifEntropyBox, vonEntropyBox, graphEntropyBox, popEntropyBox, KSEntropyBox, degreeEntropyBox)
                            
                elif methodChoice == 1:
                    # box network and box Modularity
                    boxNet = nx.Graph()
                    boxNet.add_nodes_from(boxesList)
                    boxNodes = len(boxNet)
                    boxName = str(lb)
                    for i in range(pairsNum):
                        if dist[i] > 0:
                            if boxes[node1[i]] == boxes[node2[i]]:
                                withinBox = withinBox + 1
                            else:
                                betweenBox = betweenBox + 1    
                                boxNet.add_edge(boxes[node1[i]], boxes[node2[i]])
                    modLb = np.divide(np.divide(withinBox, betweenBox, dtype = float), nBoxes, dtype = float)
                    modLbAll.append(modLb)
                    ## stats
                    if  boxNodes > 4 and nx.is_connected(boxNet):
                        boxWeight = np.divide(boxNodes, nodeNum, dtype = float)
                        boxWeightList.append(boxWeight)
                        edgePBox, radiusBox, kCoreBox, degreeAssortBox, diameterPathBox, avgPathDistBox, nEdgesBox, globalEfficiencyBox, avgDegreeBox, maxDegreeBox, spectralRadiusAdjBox, spectralGapBox, scaledSpectralRadiusBox, colorNumBox, avgClustCoeffBox, freqMBox, motifEntropyBox, vonEntropyBox, graphEntropyBox, popEntropyBox, KSEntropyBox, degreeEntropyBox = boxStats(boxNet)
            
                        # write to box file
                        outputHandleBox.write(str(fileName1) + ',' + str(lb) + ',' + str(nBoxes) + ',' + str(boxName) + ',' + str(boxNodes) + ',' + str(diameterPathBox) + ',' + str(avgPathDistBox) + ',' + str(nEdgesBox) + ',' + str(edgePBox) + ',' + str(radiusBox) + ',' + str(kCoreBox) + ',' + str(degreeAssortBox) + ',' + str(globalEfficiencyBox)  + ',' + str(avgDegreeBox) + ',' + str(maxDegreeBox) + ',' + str(spectralRadiusAdjBox) + ',' + str(spectralGapBox) + ',' + str(popEntropyBox) + ',' + str(colorNumBox) + ',' + str(avgClustCoeffBox) + ',' + str(scaledSpectralRadiusBox) + ',' + str(vonEntropyBox) + ',' + str(KSEntropyBox) + ',' + str(degreeEntropyBox) + ',' + str(graphEntropyBox) + ',' + str(motifEntropyBox) + ',' + str(100*freqMBox[0]) + ',' + str(100*freqMBox[1]) + ',' + str(100*freqMBox[2]) + ',' + str(100*freqMBox[3]) + ',' + str(100*freqMBox[4]) + ',' + str(100*freqMBox[5])     + '\n')
                    
                        # list of comparisons
                        comparisonList(edgePBox, radiusBox,  kCoreBox,  degreeAssortBox, diameterPathBox, avgPathDistBox, nEdgesBox, globalEfficiencyBox, avgDegreeBox, maxDegreeBox, spectralRadiusAdjBox, spectralGapBox, scaledSpectralRadiusBox, colorNumBox, avgClustCoeffBox, freqMBox, motifEntropyBox, vonEntropyBox, graphEntropyBox, popEntropyBox, KSEntropyBox, degreeEntropyBox)

                
                elif methodChoice == 2:
                    # box network and box Modularity
                    for i in range(pairsNum):
                        if dist[i] == 1:
                            if boxes[node1[i]] == boxes[node2[i]]:
                                withinBox = withinBox + 1
                            else:
                                betweenBox = betweenBox + 1    
                    modLb = np.divide(np.divide(withinBox, betweenBox, dtype = float), nBoxes, dtype = float)
                    modLbAll.append(modLb)
                    ## make new graph of this radius
                    edgeListStep = []
                    for i in range(pairsNum):
                        if dist[i] < lb:
                            edgeListStep.append([node1[i], node2[i]])
                    R = nx.Graph()
                    R.add_nodes_from(nodeList)
                    R.add_edges_from(edgeListStep)    
                    ## stats
                    
                    for boxName in nodeListRandom:
                        radiusNodeList = R.neighbors(boxName)
                        radiusNodeList.append(boxName)
                        boxNet = nx.Graph(G.subgraph(radiusNodeList))
                        boxNodes = len(boxNet)
                        if  boxNodes > 4 and nx.is_connected(boxNet):
                            boxWeight = np.divide(boxNodes, nodeNum, dtype = float)
                            boxWeightList.append(boxWeight)
                            edgePBox, radiusBox, kCoreBox, degreeAssortBox, diameterPathBox, avgPathDistBox, nEdgesBox, globalEfficiencyBox, avgDegreeBox, maxDegreeBox, spectralRadiusAdjBox, spectralGapBox, scaledSpectralRadiusBox, colorNumBox, avgClustCoeffBox, freqMBox, motifEntropyBox, vonEntropyBox, graphEntropyBox, popEntropyBox, KSEntropyBox, degreeEntropyBox = boxStats(boxNet)
                
                            # write to box file
                            outputHandleBox.write(str(fileName1) + ',' + str(lb) + ',' + str(nBoxes) + ',' + str(boxName) + ',' + str(boxNodes) + ',' + str(diameterPathBox) + ',' + str(avgPathDistBox) + ',' + str(nEdgesBox) + ',' + str(edgePBox) + ',' + str(radiusBox) + ',' + str(kCoreBox) + ',' + str(degreeAssortBox) + ',' + str(globalEfficiencyBox) + ',' + str(avgDegreeBox) + ',' + str(maxDegreeBox) + ',' + str(spectralRadiusAdjBox) + ',' + str(spectralGapBox) + ',' + str(popEntropyBox) + ',' + str(colorNumBox) + ',' + str(avgClustCoeffBox) + ',' + str(scaledSpectralRadiusBox) + ',' + str(vonEntropyBox) + ',' + str(KSEntropyBox) + ',' + str(degreeEntropyBox) + ',' + str(graphEntropyBox) + ',' + str(motifEntropyBox) + ',' + str(100*freqMBox[0]) + ',' + str(100*freqMBox[1]) + ',' + str(100*freqMBox[2]) + ',' + str(100*freqMBox[3]) + ',' + str(100*freqMBox[4]) + ',' + str(100*freqMBox[5])     + '\n')

                            # list of comparisons
                            comparisonList(edgePBox, radiusBox,  kCoreBox,  degreeAssortBox,  diameterPathBox, avgPathDistBox, nEdgesBox, globalEfficiencyBox, avgDegreeBox, maxDegreeBox, spectralRadiusAdjBox, spectralGapBox, scaledSpectralRadiusBox, colorNumBox, avgClustCoeffBox, freqMBox, motifEntropyBox, vonEntropyBox, graphEntropyBox, popEntropyBox, KSEntropyBox, degreeEntropyBox)

                if methodChoice == 3:
                    # box network and box Modularity
                    allBoxesDict = {}                
                    for boxName in boxesList:
                        allBoxesDict[boxName] = nx.Graph()
                    for i in range(pairsNum):
                        if dist[i] == 1:
                            if boxes[node1[i]] == boxes[node2[i]]:
                                withinBox = withinBox + 1
                                allBoxesDict[boxes[node1[i]]].add_edge(node1[i], node2[i])
                            else:
                                betweenBox = betweenBox + 1    
                    modLb = np.divide(np.divide(withinBox, betweenBox, dtype = float), nBoxes, dtype = float)
                    modLbAll.append(modLb)

        if methodChoice == 3:
            degreeEntropySelf,diameterPathSelf,avgPathDistSelf,nEdgesSelf,edgePSelf,radiusSelf,kCoreSelf ,degreeAssortSelf,globalEfficiencySelf,avgDegreeSelf,maxDegreeSelf,spectralRadiusAdjSelf,spectralGapSelf,popEntropySelf,scaledSpectralRadiusSelf,colorNumSelf,avgClustCoeffSelf,freqMBoxSelf,graphEntropySelf,motifEntropySelf,vonEntropySelf,KSEntropySelf =  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] 
        else:    
            ## final self similarity
            degreeEntropySelf = 100*np.divide(np.average(degreeEntropyBoxList, weights=boxWeightList), degreeEntropy, dtype = float)
            diameterPathSelf = 100*np.divide(np.average(diameterPathBoxList, weights=boxWeightList), diameterPath, dtype = float)
            avgPathDistSelf = 100*np.divide(np.average(avgPathDistBoxList, weights=boxWeightList), avgPathDist, dtype = float)
            nEdgesSelf = 100*np.divide(np.average(nEdgesBoxList, weights=boxWeightList), nEdges, dtype = float)
            edgePSelf = 100*np.divide(np.average(edgePBoxList, weights=boxWeightList), edgeP, dtype = float)
            radiusSelf = 100*np.divide(np.average(radiusBoxList, weights=boxWeightList), radius, dtype = float)
            kCoreSelf = 100*np.divide(np.average(kCoreBoxList, weights=boxWeightList), kCore, dtype = float)
            degreeAssortSelf = 100*np.divide(np.average(degreeAssortBoxList, weights=boxWeightList), degreeAssort, dtype = float) 
            globalEfficiencySelf= 100*np.divide(np.average(globalEfficiencyBoxList, weights=boxWeightList), globalEfficiency, dtype = float)
            avgDegreeSelf = 100*np.divide(np.average(avgDegreeBoxList, weights=boxWeightList),avgDegree , dtype = float)
            maxDegreeSelf = 100*np.divide(np.average(maxDegreeBoxList, weights=boxWeightList), maxDegree, dtype = float)
            spectralRadiusAdjSelf = 100*np.divide(np.average(spectralRadiusAdjBoxList, weights=boxWeightList), spectralRadiusAdj, dtype = float)
            spectralGapSelf = 100*np.divide(np.average(spectralGapBoxList, weights=boxWeightList), spectralGap, dtype = float)
            popEntropySelf = 100*np.divide(np.average(popEntropyBoxList, weights=boxWeightList), popEntropy, dtype = float)
            scaledSpectralRadiusSelf = 100*np.divide(np.average(scaledSpectralRadiusBoxList, weights=boxWeightList), scaledSpectralRadius, dtype = float)
            colorNumSelf = 100*np.divide(np.average(colorNumBoxList, weights=boxWeightList), colorNum, dtype = float)
            avgClustCoeffSelf = 100*np.divide(np.average(avgClustCoeffBoxList, weights=boxWeightList), avgClustCoeff,  dtype = float)            
            freqMBoxSelf = 100*np.average(freqMBoxList, weights=boxWeightList)
            graphEntropySelf = 100*np.divide(np.average(graphEntropyBoxList, weights=boxWeightList),graphEntropy , dtype = float) 
            motifEntropySelf = 100*np.divide(np.average(motifEntropyBoxList, weights=boxWeightList),motifEntropy , dtype = float)                 
            vonEntropySelf = 100*np.divide(np.mean(vonEntropyBoxList), vonEntropy, dtype = float)
            KSEntropySelf = 100*np.divide(np.mean(KSEntropyBoxList), KSEntropy, dtype = float)    
            
        #final db        
        nBoxesAll.append(1)
        steps = range(1, int(diameter)+2, 1)
        x = np.log(steps)
        y = np.log(nBoxesAll)
        slope, intercept, corr, p_value, stdErr = linregress(x,y)
        fractalDb = -1*slope
        fractalDbConstant = intercept
        fractalRSquare = np.power(corr, 2, dtype = float)
        #final modularity
        steps = range(2, int(diameter)+1, 1)
        x = np.log(steps)
        y = np.log(modLbAll)
        slope, intercept, corr, p_value, stdErr = linregress(x,y)
        fractalModularity = slope
        fractalModularityConstant = intercept
        fractalModularityRSquare = np.power(corr, 2, dtype = float)
        return fractalDb, fractalRSquare, fractalDbConstant, fractalModularity, fractalModularityConstant, fractalModularityRSquare, diameterPathSelf, avgPathDistSelf, nEdgesSelf, edgePSelf, radiusSelf, kCoreSelf, degreeAssortSelf, globalEfficiencySelf, avgDegreeSelf, maxDegreeSelf, spectralRadiusAdjSelf, spectralGapSelf, popEntropySelf,scaledSpectralRadiusSelf, colorNumSelf, avgClustCoeffSelf, vonEntropySelf, KSEntropySelf, degreeEntropySelf, graphEntropySelf, motifEntropySelf, freqMBoxSelf
                        

if __name__ == '__main__':
    try:
        argument_parser = optparse.OptionParser()            
        argument_parser.add_option('-i', metavar='IDIR', action='store', type=str, dest='input_file', default='.', help='Input file')
        argument_parser.add_option('-o', metavar='ODIR', action='store', type=str, dest='output_directory', default='output', help='Output directory')
        argument_parser.add_option('-m', action='store', dest='methodChoice', type=int, default= 3, help='self method Choice, 3 no boxing done')
        argument_parser.add_option('-q', action='store', dest='motifChoice', type=int, default= 0, help='motifChoice')

        vars2add=['VonEntropy','KSEntropy','degreeEntropy','corrPageRankfreq']
        options, args = argument_parser.parse_args()    
        pathInput = options.input_file
        pathOutput = options.output_directory    
        methodChoice = options.methodChoice
        motifChoice = options.motifChoice
        methodList = ['box', 'renormalization', 'growth','noSelf']
        #startTime = time.time()
        
        if not os.path.exists(pathOutput):
            os.makedirs(pathOutput)

        fileList = sorted(glob.glob(os.path.join(pathInput, '*links.csv')))
        # fileList = [pathInput]
        fileNum = len(fileList)

        #create motif dictionary
        if motifChoice == 1:
            motifDictCalc()

        for pathFileInput1 in fileList:
            
            dirInput1,fileInput1 = os.path.split(pathFileInput1)
            fileName1, fileExtension1 = os.path.splitext(fileInput1)

            outputName1 = pathOutput + '/' + fileName1  + '_stats_summary.csv'
            outputHandle1 = open(os.path.join(outputName1), "w")
            outputNameBox = pathOutput + '/' + fileName1  +  '_box_summary.csv'
            outputHandle1.write('file_name,comp_size,n_edges,edgeP,radius,kCore,degreeAssort,corr_path_hamm_dist,RMSE_path_hamm_dist,avg_degree,max_degree,mean_hamm_dist,mean_path_dist,diameter_Path,diameter_Hamm,corr_degree_freq,corr_eigen_cent_freq,corr_close_cent_freq,corr_bet_cent_freq,corrPageRankfreq,genetic_load,CV_freq,localOptFrac,viable_fraction,scaledSpectralRadius,spectralRadiusAdj,spectralGap,spectralRadiusHamm,popEntropy,VonEntropy,KSEntropy,degreeEntropy,average_clustering_coeficcient,global_efficiency,motif_1_star,motif_2_path,motif_3_cycle,motif_4_tailed_triangle,motif_5_envelope,motif_6_clique,graphEntropy,motif_entropy,colorNum,pathDb,pathDbConstant,pathRSquare,pathModularity,pathModularityConstant,pathModularityRSquare,diameterPathSelf,avgPathDistSelf,nEdgesSelf,edgePSelf,globalEfficiencySelf,avgDegreeSelf,maxDegreeSelf,spectralRadiusAdjSelf,spectralGapSelf,popEntropySelf,scaledSpectralRadiusSelf,colorNumSelf,avgClustCoeffSelf,vonEntropySelf,KSEntropySelf,degreeEntropySelf,graphEntropySelf,motifEntropySelf,freqMBoxSelf\n')
            if methodChoice != 3:
                outputHandleBox = open(os.path.join(outputNameBox), "w")
                outputHandleBox.write('file_name,lb,nBoxes,boxName,boxNodes,diameterPathBox,avgPathDistBox,nEdgesBox,edgePBox,radiusBox,kCoreBox,degreeAssortBox,globalEfficiencyBox,avgDegreeBox,maxDegreeBox,spectralRadiusAdjBox,spectralGapBox,popEntropyBox,colorNumBox,avgClustCoeffBox,scaledSpectralRadiusBox,vonEntropyBox,KSEntropyBox,degreeEntropyBox,graphEntropyBox,motifEntropyBox,motif_1_star,motif_2_path,motif_3_cycle,motif_4_tailed_triangle,motif_5_envelope,motif_6_clique\n')


            print('loading:', fileName1, '; Fractal_method:', methodList[methodChoice])
            ## get data
            G = nx.Graph() #fordavid what is G supposed to be
            G4 = nx.Graph()
            
            with open(pathFileInput1, 'r') as csvfile1:
                original = csv.reader(csvfile1, delimiter=',')
                #Params = next(original, None)
                node1 = []
                node2 = []
                hammDist  = []
                pathDist = []
                freqDict = {}
                nodeDict = {}
                for row in original:
                    node1Temp = int(row[0])
                    node2Temp = int(row[1])
                    hammDistTemp = int(float(row[2]))
                    pathDistTemp = int(float(row[4]))
                    freqDict[node1Temp] = float(row[5])
                    freqDict[node2Temp] = float(row[6])
                    nodeDict[node1Temp] = 1
                    nodeDict[node2Temp] = 1
                    node1.append(node1Temp)
                    node2.append(node2Temp)
                    hammDist.append(hammDistTemp)
                    pathDist.append(pathDistTemp)
                    if pathDistTemp == 1:
                        G.add_edge(node1Temp, node2Temp)
                        G4.add_edge(node1Temp, node2Temp)
                    elif pathDistTemp < 4:
                        G4.add_edge(node1Temp, node2Temp)
                # print(G.edges())
                # exit()
               
            ## get number of nodes
            nodes = len(nodeDict.keys())
            del nodeDict
            pairsNum = len(node1)
            nodeList = range(nodes)
            ## node frequencies
            freqCount = np.zeros(nodes)
            for k in nodeList:
                freqCount[k] = freqDict[k]
            del freqDict


            #constants
            posNum = 15
            letters = 2
            u = 0.000115

            if nodes > 3:


                #print 'correlations'
                ## distance properties
                diameterHamm = max(hammDist)
                avgHammDist = np.mean(hammDist)
                corrPathHammDist = np.corrcoef(hammDist, pathDist)
                RMSEPathHammDist = mean_squared_error(hammDist, pathDist)
                del hammDist, 
                degreeRaw = list(G.degree())
                degree = []
                for i in degreeRaw:
                    degree.append(i)
                degreeNorm = np.divide(degree, np.sum(degree), dtype = float)
                eigenCentRaw = list(nx.eigenvector_centrality_numpy(G, weight=None).values())
                eigenCent = np.divide(eigenCentRaw, sum(eigenCentRaw), dtype = float)
                closeCent = list(nx.closeness_centrality(G).values())
                betCent = list(nx.betweenness_centrality(G).values())
                pageRank = list(nx.pagerank_numpy(G).values())
                # correlations
                corrDegreeFreq = np.corrcoef(freqCount, degree)
                corrEigenCentFreq = np.corrcoef(freqCount, eigenCent)
                corrCloseCentFreq = np.corrcoef(freqCount, closeCent)
                corrBetCentFreq = np.corrcoef(freqCount, betCent)
                corrPageRankfreq = np.corrcoef(freqCount, pageRank)                
                # del degreeNorm, eigenCentRaw,  closeCent, betCent, pageRank
                print(freqCount)
                print(pageRank)
                print(corrPageRankfreq[0][1])
                # exit()

            
                #print 'net stats'
                edgeP, radius, kCore, degreeAssort, diameterPath, avgPathDist, nEdges, globalEfficiency, avgDegree, maxDegree, spectralRadiusAdj, spectralGap, scaledSpectralRadius, colorNum, avgClustCoeff, freqM, motifEntropy, vonEntropy, graphEntropy, popEntropy, KSEntropy, degreeEntropy = boxStats(G)
                spectralRadiusHamm = 1#max(abs(eValsH))

                #print 'genetic stats'
                #Nimwegen stats, change the sequence length (posNum) if needed
                nReads = sum(freqCount)
                freqCountRel = np.divide(freqCount, nReads, dtype = float)
                d1 = sum(degree*freqCountRel)
                neighbors = posNum*(letters-1)
                uPosNum = posNum * u
                geneticLoad = uPosNum * (1 - np.divide(scaledSpectralRadius, neighbors, dtype = float)) 
                # coefficient of variation of the frequencies
                CV = variation(freqCountRel)
                #local optima
                localOptNum = 0
                for n0 in G.nodes():#range(76,77):
                    flagMax = 0
                    flagMin = 0
                    for n1 in G.neighbors(n0):
                        #print freqCount[n0], freqCount[n1]
                        if freqCountRel[n0] > freqCountRel[n1]:
                           flagMax = 1
                        if freqCountRel[n0] < freqCountRel[n1]:
                           flagMin = 1                          
                           
                    if flagMax == 1 and flagMin == 0:
                        #rint n0, n1, freqCount[n0], freqCount[n1]
                        localOptNum = localOptNum + 1 
                localOptFrac = np.divide(localOptNum, nodes, dtype = float)          
                
                
                vs = 1 - uPosNum * (1 - (np.divide(degree, neighbors, dtype = float)))
                viableFraction = sum(vs*eigenCent)
                del freqCount, eigenCent, freqCountRel, d1, vs

                # Fractal dimension
                #print 'fractal'
                if methodChoice == 3:
                    fractalDb, fractalRSquare, fractalDbConstant, fractalModularity, fractalModularityConstant, fractalModularityRSquare, diameterPathSelf, avgPathDistSelf, nEdgesSelf, edgePSelf, radiusSelf, kCoreSelf, degreeAssortSelf, globalEfficiencySelf, avgDegreeSelf, maxDegreeSelf, spectralRadiusAdjSelf, spectralGapSelf, popEntropySelf,scaledSpectralRadiusSelf, colorNumSelf, avgClustCoeffSelf, vonEntropySelf, KSEntropySelf, degreeEntropySelf, graphEntropySelf, motifEntropySelf, freqMBoxSelf  = fractalCalc(pathDist, nodes, methodChoice)    
             
                
                if methodChoice == 0 or methodChoice == 2:
                    # placing the original network as the last box
                    fractalDb, fractalRSquare, fractalDbConstant, fractalModularity, fractalModularityConstant, fractalModularityRSquare, diameterPathSelf, avgPathDistSelf, nEdgesSelf, edgePSelf, radiusSelf, kCoreSelf, degreeAssortSelf, globalEfficiencySelf, avgDegreeSelf, maxDegreeSelf, spectralRadiusAdjSelf, spectralGapSelf, popEntropySelf,scaledSpectralRadiusSelf, colorNumSelf, avgClustCoeffSelf, vonEntropySelf, KSEntropySelf, degreeEntropySelf, graphEntropySelf, motifEntropySelf, freqMBoxSelf  = fractalCalc(pathDist, nodes, methodChoice)    
                    outputHandleBox.write(str(fileName1)  +  ',' +  str(diameterPath + 1) + ',' + str(1) + ',' + str('original') + ',' + str(nodes) + ',' + str(diameterPath) + ',' + str(avgPathDist) + ',' + str(nEdges) + ',' + str(edgeP) + ',' + str(radius) + ',' + str(kCore) + ',' + str(degreeAssort) + ',' + str(globalEfficiency) + ',' + str(avgDegree) + ',' + str(maxDegree) + ',' + str(spectralRadiusAdj) + ',' + str(spectralGap) + ',' + str(popEntropy) + ',' + str(colorNum) + ',' + str(avgClustCoeff) + ',' + str(scaledSpectralRadius) + ',' + str(vonEntropy) + ',' + str(KSEntropy) + ',' + str(degreeEntropy) + ',' + str(graphEntropy) + ',' + str(motifEntropy) + ',' +  str(100*freqM[0]) + ',' +  str(100*freqM[1]) + ',' +  str(100*freqM[2]) + ',' +  str(100*freqM[3]) + ',' +  str(100*freqM[4]) + ',' +  str(100*freqM[5]) + '\n')
                elif methodChoice == 1:    
                    # placing the original network as the first box
                    outputHandleBox.write(str(fileName1)  +  ',' +  str(1) + ',' + str(nodes) + ',' + str('original') + ',' + str(nodes) + ',' + str(diameterPath) + ',' + str(avgPathDist) + ',' + str(nEdges) + ',' + str(edgeP) + ',' + str(radius) + ',' + str(kCore) + ',' + str(degreeAssort) + ',' + str(globalEfficiency) + ',' + str(avgDegree) + ',' + str(maxDegree) + ',' + str(spectralRadiusAdj) + ',' + str(spectralGap) + ',' + str(popEntropy) + ',' + str(colorNum) + ',' + str(avgClustCoeff) + ',' + str(scaledSpectralRadius) + ',' + str(vonEntropy) + ',' + str(KSEntropy) + ',' + str(degreeEntropy) + ',' + str(graphEntropy) + ',' + str(motifEntropy) + ',' +  str(100*freqM[0]) + ',' +  str(100*freqM[1]) + ',' +  str(100*freqM[2]) + ',' +  str(100*freqM[3]) + ',' +  str(100*freqM[4]) + ',' +  str(100*freqM[5]) + '\n')
                    fractalDb, fractalRSquare, fractalDbConstant, fractalModularity, fractalModularityConstant, fractalModularityRSquare, diameterPathSelf, avgPathDistSelf, nEdgesSelf, edgePSelf, radiusSelf, kCoreSelf, degreeAssortSelf, globalEfficiencySelf, avgDegreeSelf, maxDegreeSelf, spectralRadiusAdjSelf, spectralGapSelf, popEntropySelf,scaledSpectralRadiusSelf, colorNumSelf, avgClustCoeffSelf, vonEntropySelf, KSEntropySelf, degreeEntropySelf, graphEntropySelf, motifEntropySelf, freqMBoxSelf  = fractalCalc(pathDist, nodes, methodChoice)    

                #save summary file
                print('we did it')
                outputHandle1.write(str(fileName1) + ',' +  str(nodes) + ',' +  str(nEdges) + ',' +  str(edgeP) + ',' + str(radius) + ',' + str(kCore) + ',' + str(degreeAssort) + ',' +  str(corrPathHammDist[0][1]) +  ',' +  str(RMSEPathHammDist) +  ',' +  str(avgDegree) + ',' +  str(maxDegree) +  ',' +  str(avgHammDist) + ',' +  str(avgPathDist) + ',' +  str(diameterPath) + ',' +  str(diameterHamm) + ',' +  str(corrDegreeFreq[0][1]) + ',' +  str(corrEigenCentFreq[0][1]) + ',' +  str(corrCloseCentFreq[0][1]) + ',' +  str(corrBetCentFreq[0][1]) + ',' +  str(corrPageRankfreq[0][1]) + ',' +  str(geneticLoad) + ',' +  str(CV) + ',' +  str(localOptFrac) + ',' +  str(viableFraction) +  ',' + str(scaledSpectralRadius) + ',' +  str(spectralRadiusAdj) + ',' + str(spectralGap) + ',' +  str(spectralRadiusHamm) + ',' +  str(popEntropy) + ',' +  str(vonEntropy) + ',' +  str(KSEntropy)  + ',' +  str(degreeEntropy) + ',' +  str(avgClustCoeff) + ',' +  str(globalEfficiency) + ',' +  str(100*freqM[0]) + ',' +  str(100*freqM[1]) + ',' +  str(100*freqM[2]) + ',' +  str(100*freqM[3]) + ',' +  str(100*freqM[4]) + ',' +  str(100*freqM[5]) + ',' +  str(graphEntropy) + ',' +  str(motifEntropy) +  ',' +  str(colorNum) + ',' +  str(fractalDb) + ',' +  str(fractalDbConstant) + ',' +  str(fractalRSquare) + ',' +  str(fractalModularity) + ',' +  str(fractalModularityConstant) + ',' +  str(fractalModularityRSquare) + ',' +  str(diameterPathSelf) + ',' +  str(avgPathDistSelf) + ',' +  str(nEdgesSelf) + ',' +  str(edgePSelf) + ',' +  str(globalEfficiencySelf) + ',' +  str(avgDegreeSelf) + ',' +  str(maxDegreeSelf) + ',' +  str(spectralRadiusAdjSelf) + ',' + str(spectralGapSelf) + ',' +  str(popEntropySelf)+ ',' +  str(scaledSpectralRadiusSelf) + ',' +  str(colorNumSelf) + ',' +  str(avgClustCoeffSelf) + ',' +  str(vonEntropySelf)+ ',' +  str(KSEntropySelf)+','  +  str(degreeEntropySelf) + ',' +  str(graphEntropySelf) + ',' +  str(motifEntropySelf) + ',' +  str(freqMBoxSelf)     + '\n')

        if methodChoice != 3:
            outputHandleBox.close()
        outputHandle1.close()
        ###mark the end time
        #endTime = time.time()
        #workTime =  endTime - startTime
        #print 'time', workTime

    except KeyboardInterrupt:
        exit(-1)

