#!/home/rs2474/anaconda3/envs/imktData/bin/python


import sys
import os
import time
import numpy as np
import pandas as pd
import pyfaidx as px


"""
	Convert fasta file to iMKT input files

    These are the requirements of the multi-FASTA files:

    1. Sequences need to be protein-coding sequences (CDS) and aligned.
    2. The first sequence in the file is the reference sequence, which will be used to compute the degeneracy of each nucleotide position: 0 as 0-fold site, 2 as 2-fold site, and 4 as 4-fold site.
    3. The last sequence in the file is the outgroup sequence, which will be used to infer the ancestral alleles and compute divergence data.
    4. Any other sequences between the first and the last sequence are considered polymorphic sequences, which will be used to compute the Derived Allele Frequency of polymorphic sites. At least two sequences are needed.
"""


def degenerancy(data, codonDict='standard'):

	# DEGENERANCY DICTIONARIES
	standardDict = {
		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202',
		'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004',
		'TAT': '002', 'TAC': '002', 'TAA': '022', 'TAG': '002',
		'TGT': '002', 'TGC': '002', 'TGA': '020', 'TGG': '000',
		'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204',
		'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004',
		'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002',
		'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204',
		'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000',
		'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004',
		'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002',
		'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202',
		'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004',
		'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004',
		'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002',
		'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}

	if(codonDict == 'standard'):
		degenerateCodonTable = standardDict

	degenerancy = ''
	for i in range(0, len(data),3):
		codon = data[i:i+3]
		if('N' in codon or '-' in codon):
			degenerancy += codon
		else:
			degenerancy += degenerateCodonTable[codon]

	return(degenerancy)


def aa(data):

	#AA DICTIONARIES
	
	# DEGENERANCY DICTIONARIES	
	degenerateCodonTable = {
		"TTC":'PHE',"TTT":'PHE',#PHE 002
		"TTA":'LEU',"TTG":'LEU',#LEU 002
		"CTT":'LEU',"CTC":'LEU',"CTA":'LEU',"CTG":'LEU',#LEU 204
		"ATG":'MET',#MET 000
		"ATT":'ILE',"ATC":'ILE',"ATA":'ILE',#ILE 003
		"GTC":'VAL',"GTT":'VAL',"GTA":'VAL',"GTG":'VAL',#VAL 004
		"TCT":'SER',"TCA":'SER',"TCC":'SER',"TCG":'SER',#SER 004
		"CCT":'PRO',"CCA":'PRO',"CCC":'PRO',"CCG":'PRO',#PRO 004
		"ACT":'THR',"ACA":'THR',"ACC":'THR',"ACG":'THR',#THR 004
		"GCT":'ALA',"GCA":'ALA',"GCC":'ALA',"GCG":'ALA',#ALA 004
		"TAT":'TRY',"TAC":'TRY',#TYR 002
		"TAA":'STP',"TAG":'STP',"TGA":'STP',#STOP 000
		"CAT":'HIS',"CAC":'HIS',#HIS 002
		"CAA":'GLN',"CAG":'GLN',#GLN 002
		"AAT":'ASN',"AAC":'ASN',#ASN 002
		"AAG":'LYS',"AAA":'LYS',#LYS 002
		"GAT":'ASP',"GAC":'ASP',#ASP 002
		"GAA":'GLU',"GAG":'GLU',#GLU 002
		"TGT":'CYS',"TGC":'CYS',#CYS 002
		"TGG":'TRP',#TRP 000
		"CGT":'ARG',"CGG":'ARG',"CGA":'ARG',"CGC":'ARG',#ARG 204
		"AGA":'ARG',"AGG":'ARG',#ARG 002
		"AGT":'SER',"AGC":'SER',#SER 002
		"GGT":'GLY',"GGA":'GLY',"GGC":'GLY',"GGG":'GLY',#GLY 004
	}

	degenerancy = ''
	for i in range(0, len(data),3):
		codon = data[i:i+3]
		if('N' in codon or '-' in codon or 'M' in codon):
			degenerancy += codon
		else:
			degenerancy += degenerateCodonTable[codon]

	return(degenerancy)


def sequencesToMatrix(multiFasta):
	"""
		Convert the multi-FASTA file to a matrix of arrays with degenerancy codes

		multiFasta: the multiFasta file read from px.Fasta

		Output: a matrix of dim NxM, where N is the sample size and M is the number of sites per sequence
	"""
	# Extract samples from fastas
	samples = list(multiFasta.keys())
	# get the length of the sequence (should be the same across the population) 
	seqLen = len(multiFasta[samples[0]][:].seq)

	if((seqLen % 3) != 0):
		print('cdsLength')
		sys.exit('cdsLength')

	# Create empty array with ndimesions equal to multi-Fasta lines and length
	matrix = np.empty([len(samples),len(multiFasta[samples[0]][:].seq)],dtype='str')

	# List to append indexes if whole sequence at any population is len(seq) * 'N'
	deleteIndex = list()

	# Iter fasta to add sequence to matrix
	for i in range(1,len(samples)):
		# Extract each sample sequence
		tmp = multiFasta[samples[i]][:].seq
		if(len(tmp) != seqLen):
			print('errorAlign')
			sys.exit('errorAlign')

		if(tmp == ('N' * len(tmp))):
			deleteIndex.append(i)
		else:
			matrix[i] = list(tmp)

	# Delete lines
	matrix = np.delete(matrix,deleteIndex,0)

	degenCode = degenerancy(multiFasta[samples[0]][:].seq)
	# Put degenerancy in first ndarray element
	matrix[0] = list(degenCode)
	
	return(matrix)


def uSfsFromFastaDegen(sequenceMatrix):
	output = list()
	it = -1
	cpos=np.array([[i,i+1,i+2] for i in range(0,sequenceMatrix.shape[1]-1,3)])
	remove_sites = {'0': 0, '2': 0, '3': 0, '4': 0}

	for x in np.nditer(sequenceMatrix, order='F',flags=['external_loop']):
		degen = x[0]
		AA = x[-1]
		it += 1
		# Undefined Ancestral Allele. Try to clean out of the loop
		if(AA == 'N' or AA == '-' or degen == '-'):
			next
		# Monomorphic sites. Try to clean out of the loop
		elif(np.unique(x[1:][np.where(x[1:]!='N')]).shape[0] == 1):
			next
		else:
			idx = cpos[np.where(cpos==it)[0],:].flatten()
			codons = sequenceMatrix[:,idx]
			
			## filter out complex codons (a codon with >2 site changes)
			unique_change_pos = np.unique(codons[1:],axis=1)
			unique_change_pos = np.sort(unique_change_pos,axis=0)
			number_unique_changes = (unique_change_pos[1:,:]!=unique_change_pos[:-1,:]).sum(axis=0)
			if (len(number_unique_changes[number_unique_changes >= 1]) > 1):
				remove_sites[degen] += 1
				next

			# check the sampled alleles for polymorphisms
			pol = x[x!='N'][1:-1]

			if(degen == '0'):
				functionalClass = '0fold'
			elif(degen == '4'):
				functionalClass = '4fold'
			else:
				## filter out complex codons (one codon site with >2 changes)
				uniqueCodons = np.unique(codons[1:],axis=0)
				aaList = np.unique([aa("".join(item)) for item in uniqueCodons])
				if(aaList.shape[0] <= 2):
					if(aaList.shape[0] == 1):
						functionalClass = '4fold'
					else:
						functionalClass = '0fold'
				else:
					remove_sites[degen] += 1
					next

			# Check if pol != AA and monomorphic
			if((np.unique(pol).shape[0] == 1) and (np.unique(pol)[0] != AA)):
				div = 1; AF = 0
				tmp = [AF,div,functionalClass]
				output.append(tmp)
			else:
				AN = x[1:-1].shape[0]
				AC = pd.DataFrame(data=np.unique(x[1:-1], return_counts=True)[1],index=np.unique(x[1:-1], return_counts=True)[0])
				div = 0
				if(AA not in AC.index):
					next
				else:
					AC = AC[AC.index!=AA]
					if(len(AC) == 0):
						next
					else:
						AF = AC.iloc[0]/AN
						AF = AF.iloc[0]
						tmp = [AF,div,functionalClass]
						output.append(tmp)

	return (output, remove_sites)


def formatSfs(sequenceMatrix,rawSfsOutput,remove_sites,dafFile,divFile,append=False):

	df = pd.DataFrame(rawSfsOutput)
	df['id'] = 'uploaded'
	df.columns = ['derivedAlleleFrequency','d','functionalClass','id']

	# Extract divergence data
	div = df[['id','functionalClass','d']]
	div = div[div['d']!=0]
	div = div.groupby(['id','functionalClass'])['d'].count().reset_index()
	div = div.pivot_table(index=['id'],columns=['functionalClass'],values='d').reset_index()
	try:
		div = div[['0fold','4fold']]
	except:
		if('4fold' in div.columns):
			div = div[['4fold']]
			div['0fold'] = 0
		elif('0fold' in div.columns):
			div = div[['0fold']]
			div['4fold'] = 0

	## calculate the total number of analyzed sites in the MKT analysis
	## mi: analyzed non-synonymous sites = 0fold site + 2/3 * 2fold sites
	## m0: analyzed synonymous sites = 4fold site + 1/3 * 2fold sites
	## di, d0: divergent non-synonymous sites vs synonymous sites
	div['mi'] =  sequenceMatrix[0][np.where(sequenceMatrix[0]=='0')].shape[0]-remove_sites['0'] + (sequenceMatrix[0][np.where(sequenceMatrix[0]=='2')].shape[0]-remove_sites['2'])//3*2
	div['m0'] =  sequenceMatrix[0][np.where(sequenceMatrix[0]=='4')].shape[0]-remove_sites['4'] + (sequenceMatrix[0][np.where(sequenceMatrix[0]=='2')].shape[0]-remove_sites['2'])//3
	div = div.rename(columns={'0fold':'Di','4fold':'D0','mi':'mi','m0':'m0'})
	div = div[['Di', 'D0', 'mi', 'm0']]

	# Create SFS pd.DataFrame by functionClass and 20 frequency bin
	daf = df[df['d']!=1][['derivedAlleleFrequency','functionalClass','id']]
	bins = np.arange(0.,1.05,0.05)
	labels = np.arange(0.025,1.0,0.05).tolist()
	daf['categories'] = pd.cut(daf['derivedAlleleFrequency'],bins=bins,labels=labels)
	daf = daf.groupby(['functionalClass','id','categories']).count().reset_index()
	sfs = pd.DataFrame({'daf':daf['categories'].unique(),'P0':daf[daf['functionalClass']=='4fold']['derivedAlleleFrequency'].reset_index(drop=True),'Pi':daf[daf['functionalClass']=='0fold']['derivedAlleleFrequency'].reset_index(drop=True)})

	sfs = sfs[['daf','P0','Pi']]
	sfs['P0'] = sfs['P0'].fillna(0)
	sfs['Pi'] = sfs['Pi'].fillna(0)
	sfs['daf'] = sfs['daf'].apply(lambda x: round(x,3))

	if(append is True):
		sfs.to_csv(dafFile,sep='\t',header=True,index=False,mode='a')
		div.to_csv(divFile,sep='\t',header=True,index=False,mode='a')
	else:
		sfs.to_csv(dafFile,sep='\t',header=True,index=False)
		div.to_csv(divFile,sep='\t',header=True,index=False)

        
def main():
    multiFasta = sys.argv[1]
    dafFile = sys.argv[2]
    divFile = sys.argv[3]
#     dafFile = multiFasta.split('.')[0]+'.daf'
#     divFile = multiFasta.split('.')[0]+'.div'
    
    # Open multi-Fasta
    fasta_file = px.Fasta(multiFasta,duplicate_action='first',sequence_always_upper=True,read_long_names=True)

    # Create ndarray with sequences
    multiFastaMatrix = sequencesToMatrix(fasta_file)
    if(multiFastaMatrix.shape[0] < 4):
        print('numberOfLines')
        sys.exit('numberOfLines')
        
    # Estimating SFS
    rawSfs, remove_sites = uSfsFromFastaDegen(multiFastaMatrix)
    formatSfs(multiFastaMatrix,rawSfs,remove_sites,dafFile,divFile)
    
    print('Convert multiFASTA file {} to {} and {}'.format(multiFasta, dafFile, divFile))
    
if __name__ == "__main__":
    main()
