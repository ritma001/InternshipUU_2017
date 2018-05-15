#!/usr/bin/env python3
# coding = "utf-8
# Estabished date: 20/12/2017	Last modification: /2017
# project description: 
'''
Aim: to convert inserted nucleotide sequences of naive TCR 
to the amino acid sequence
'''
## submitted as part of the MSc internship project at Utrecht University ##
################################################################################
import sys
from sys import argv
import itertools
import numpy as np

codon_table = {
    'F': ['TTT', 'TTC'], 'L': ['CTT', 'CTC', 'TTA', 'TTG', 'CTA', 'CTG'],	
	'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'Y': ['TAT', 'TAC'], 
	'C': ['TGT', 'TGC'], '*': ['TAA', 'TAG', 'TGA'], 'H': ['CAT', 'CAC'],
	'I': ['ATT', 'ATC', 'ATA'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'W': ['TGG'],
	'Q': ['CAA', 'CAG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 
	'E': ['GAA', 'GAG'],  'G': ['GGT', 'GGC', 'GGA', 'GGG'], 'M': ['ATG'], 
	'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'],
	'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'D': ['GAT', 'GAC'],
	'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'NA':[]}
# '*' means stop codon

def getAminoacid(codon, codon_table):
	"""
	extract amino acid from the corresponding codon
	"""
	for k,v in codon_table.items():
		if codon in v:
			return(k)

def getCodon(nucleotide):
	"""
	divide a nucleotide sequence into the codon 
	"""
	codon_list = []
	i = 0
	while i < len(nucleotide)-2:
		codon = nucleotide[i:i+3]
		codon_list.append(codon)
		i += 3
	return(codon_list)
		
def convertCodon(nucleotide, codon_table):
	"""
	convert amino acid a list of codons to peptide sequences
	"""
	peptide = ''.join([getAminoacid(cd, codon_table) for cd in getCodon(nucleotide)])
	return(peptide)
	
################################################################################
if __name__ == "__main__":
	# create a list of cdr3 of input file
	output = open(argv[2], "w")
	with open(argv[1], "r") as cdr3_list:
		for cdr3 in cdr3_list:
			x = convertCodon(cdr3.rstrip(), codon_table)
			output.write(x+"\n")
	output.close()
