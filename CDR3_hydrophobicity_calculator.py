#!/usr/bin/env python3
# coding = "utf-8
# Estabished date: 1/11/2017	Last modification: 2/11/2017
# project description: 
'''
Aim: to calculate hydrophobicity of amino acid sequence of CDR3 using the 
Hydropathy scale from Kyte, J. and Doolittle, R.F., J. Mol. Biol., 157, 105-132 
(1982). (also available on 
http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/abbreviation.html)
'''
## submitted as part of the MSc internship project at Utrecht University ##
################################################################################
import pandas as pd
import numpy as np
import sys
# initialise the hydropathy scale
hydropathy_scale = {'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 
					'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'S': -0.8,
					'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5,
					'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5}
					
def calculateHydrophobicity(cdr3_seq, hydropathy_scale):
	"""
	calculate hydropathy score by sum the values of hydropathy scales of all
	amino acids in a certain cdr3 sequence then, devideded by the seq length 
	"""
	cdr3_ls = list(cdr3_seq)
	sum_hdp = sum([hydropathy_scale[aa] for aa in cdr3_ls])
	result = sum_hdp/len(cdr3_ls)
	return(result)
	
def importDataframe(input_file):
	# each row of an input file contains rowname\tseq
	"""
	import the cdr3 column from cdr3_df dataframe
	"""
	cdr3_dict = dict()
	with open(input_file,'r') as f:
		lines = [line.rstrip() for line in f]
		for line in lines:
			k,v = line.split("\t")
			cdr3_dict[k] = v
	return(cdr3_dict)
	
def iterateCalculation(cdr3_seq, hydropathy_scale):
	"""
	calculate hydrophobicity while iterate through a list of CDR3 sequeces 
	"""
	cdr3_ar = np.array(cdr3_seq)
	vectorize_fun = np.vectorize(calculateHydrophobicity)
	cdr3_hydropathy_score = vectorize_fun(cdr3_ar, hydropathy_scale)
	return(cdr3_hydropathy_score)
	
def exportHydropathy(score, output_file):
	"""
	write the hydropathy score in a txt file
	"""
	with open(output_file, "w") as of:
		# convert an array to list
		score_num = score.tolist()
		# convert numeric score to character
		score_char = ['{:.2f}'.format(score) for score in score_num]
		# write in the score in the separate line
		of.write("\n".join(score_char))
	
	
################################################################################
if __name__ == "__main__":
	cdr3_dict = importDataframe(sys.argv[1]) 
	# extract only the cdr3 column
	cdr3_seq = list(cdr3_dict.values())
	# calculate hydropathy score then coverst an array of score into a list
	hydropathy_calculation = iterateCalculation(cdr3_seq, hydropathy_scale)
	# export the score
	exportHydropathy(hydropathy_calculation, sys.argv[2])
	
	
	
