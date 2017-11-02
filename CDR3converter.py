#!/usr/bin/env python3
# coding = "utf-8
# Estabished date: 4/10/2017	Last modification: 1/11/2017
# project description: 
'''
Aim: to convert the amino acid sequence of CDR3 to nucleotide sequences, which
will be used as the input for IGoR algorithm
'''
## submitted as part of the MSc internship project at Utrecht University ##
################################################################################
import sys
from sys import argv
import itertools
import numpy

RNA_codon_table = {
    'F': ['TTT', 'TTC'], 'L': ['CTT', 'CTC', 'TTA', 'TTG', 'CTA', 'CTG'],	
	'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'Y': ['TAT', 'TAC'], 
	'C': ['TGT', 'TGC'], '*': ['TAA', 'TAG', 'TGA'], 'H': ['CAT', 'CAC'],
	'I': ['ATT', 'ATC', 'ATA'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'W': ['TGG'],
	'Q': ['CAA', 'CAG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG'], 'E': ['GAA', 'GAG'],
	'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'],
	'R': ['AGA', 'AGG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'D': ['GAT', 'GAC'],
	'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'], 
	'M': ['ATG']}
# '*' means stop codon

def readFasta(file):
	"""
	read CDR3 amino acid sequence from fasta file, return dict of sequences
	"""
	fasta = dict()
	with open(file,'r') as f:
		for line in f:
			if line.startswith('>'):
				header = line.rstrip()
				seq = ''
			else:
				seq += line.rstrip()
			fasta[header]=seq
	return fasta

def convertAminoacid(cdr3, RNA_codon_table):
	"""
	convert amino acid to the corresponding mRNA, 
	return {cdr3:[[mrna1, mrna1.1],[mrna2,mrna2.1, mrna2.2], ...]}
	the number (nr.) of sulists = nr. of amino acids in certain cdr3
	"""
	cdr3mrna_dict = dict()
	# assign the CDR3 seq as the key and a list of its amino acid as the values
	aa_list = [[e] for e in list(cdr3)]
	# look up in the codon tabke the corresponding mRNA while looping
	i = 0
	while i < len(aa_list):
		aa = aa_list[i][0]
		aa_list[i].extend(RNA_codon_table[aa])
		i += 1			
	# remove the first element, that is amino acid
	mrna = [e[1:] for e in aa_list]
	# add a list of corresponding mRNA to cdr3mrna_dict
	cdr3mrna_dict[cdr3] = mrna
	return(cdr3mrna_dict)
	
def combinemRNA(cdr3mrna_dict):
	"""
	generate all possible mRNA seqs by combining the mRNA sequences correspoding
	to each amino acid in certain CDR3
	"""
	mrna_sublist = list(cdr3mrna_dict.values())[0]
	premrna = list(itertools.product(*mrna_sublist))
	mrna = list(''.join(mrna) for mrna in premrna)
	return(mrna)
	
################################################################################
if __name__ == "__main__":
	x = readFasta(argv[1])
	for v in list(x.values()):
		y = convertAminoacid(v, RNA_codon_table)	
		z = combinemRNA(y)
		print(len(z))
		# get the header of each mRNA seq (identicial to original CDR3 header)
		k = [header for header, seq in x.items() if seq == v][0]
		# print the sequences in fasta format using identical fasta header as
		# the CDR3 but added the number of corresponding sequences
		with open(argv[2], 'w') as out:
			for j in range(len(z)):
				out.write(k+str(j)+'\n'+z[j]+'\n')