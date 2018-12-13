#!/usr/bin/python

import sys
import numpy as np
import pandas as pd 
from collections import defaultdict


# Generate the the matrix file from the network
# relationship file

file = sys.argv[1]

def generateMatrix(file):
	f = open(file, "r")
	header = f.readline()
	allNodes = defaultdict(int)

	# Loop through lines
	networkData = list()
	for line in f:
		(start, end, desc, direction) = line.strip().split("\t")
		networkData.append([start, end, desc, direction])
		allNodes[start] += 1
		allNodes[end] += 1
	f.close()

	# Count number of genes
	geneCount = len(allNodes)
	allGenes = allNodes.keys()

	# Build the matrix of the pathway
	pathwayMat = np.zeros((geneCount,geneCount), dtype=int)
	#print allNodes.keys().index("")
	for data in networkData:
		(start, end, desc, direction) = data
		startIndex = allGenes.index(start)
		endIndex = allGenes.index(end)

		if direction == "Enhance":
			pathwayMat[startIndex, endIndex] = 1
		elif direction == "Repress":
			pathwayMat[startIndex, endIndex] = -1
			



	df = pd.DataFrame(pathwayMat, columns=allGenes, index=allGenes)
	df.to_csv("PDL1_pathway.matrix.csv")
	#df.to_csv("PDL1_pathway.matrix.csv", columns=allGenes, index=allGenes)
	print df

	# pathwayMat


	# return pathwayMat


print file
generateMatrix(file)