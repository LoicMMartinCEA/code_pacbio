#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''

'''

### MODULES IMPORT ###
import locale
import subprocess
import argparse
import gzip
import pickle
import re
import scipy.stats as stats
import statistics
import warnings

warnings.filterwarnings("ignore")

# American global settings
locale.setlocale(locale.LC_ALL, 'C')


### FUNCTIONs ###
# execute terminal commands and save the output to a logfile
def run(cmd, logfile):
        p = subprocess.Popen(cmd, shell=True, universal_newlines=True, stdout=logfile)
        ret_code = p.wait()
        logfile.flush()
        return ret_code
        
#retrieve list of lists based on 2 separators
def double_split(Text, sep1, sep2):
	# further parse options
	WholeL = []
	MainL = Text.split(sep1)
	for text in MainL:
		WholeL.append(text.strip().split(sep2))

	return WholeL
          
# get the reverse complement of a DNA sequence
def rc_dna(sequence):
	
	rc_sequence = ''
	reverse = sequence[::-1]
	
	# complementary bases dictionary
	complement_nt = { 'A':'T', 'C':'G', 'G':'C', 'T':'A', 'R':'Y', 'Y':'R', 'S':'S', 'W':'W', 'K':'M', 'M':'K', 'B':'V', 'V':'B', 'D':'H', 'H':'D', 'N':'N' }

	for n in range(0, len(reverse)):
		if reverse[n] in complement_nt.keys():
			rc_sequence += complement_nt.get(reverse[n], "N")

	return rc_sequence

# isolate the region of interest
def find_between( text, left, right ):
	try:
		start = text.index( left ) + len( left )
		end = text.index( right, start )
		return text[start:end]
	except ValueError:
		return ""
		
# convert IUPAC nucleotide sequences to regular expressions
def Nt2RE(seq):
	Nt2REDic = {
	'A':'A', 'C':'C', 'G':'G', 'T':'T', 'C':'C', 'W':'[AT]', 'S':'[CG]', 'M':'[AC]', 'K':'[GT]', 'R':'[AG]', 'Y':'[CT]', 'B':'[CGT]', 'D':'[AGT]', 'H':'[ACT]', 'V':'[ACG]', 'N':'[ACGT]' 
	}
	for nt in Nt2REDic.keys():
		seq = seq.replace(nt, Nt2REDic[nt])
	
	return seq

# translate a DNA sequence to protein sequence
def translate_dna(sequence):
	proteinsequence = ''

	codontable = {
	    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
	    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
	    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
	    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
	    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
	    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
	    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
	    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
	    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
	    }

	for n in range(0,len(sequence),3):
		if sequence[n:n+3] in codontable.keys():
		    proteinsequence += codontable[sequence[n:n+3]]
		else:
		    proteinsequence += "#" # non-expected codon (maybe codon composed of non standard DNA bases A, C, G, T or incomplete)
	sequence = ''
	return proteinsequence 

# get fasta file information to a dictionary
def get_fasta(FastaFN):
	FastaDic = {}
	Fasta = get_file_content(FastaFN)
		
	for line in Fasta:
		splitted_line = line.split(":")
		if line[0] == ">":
			desc = line.strip(">")
			FastaDic[desc] = ""
		else:
			line = line.strip("\n").strip(" ")
			if len(line) > 1 and re.match(r'[ACGTacgt]', line, re.M|re.I):
				FastaDic[desc] = FastaDic[desc] + line
	
	return FastaDic

# transfer the content of a file to the RAM  (handle gzip compressed or uncompressed files by extension)
def get_file_content(filename):		
	extension = filename.split(".")[-1]	
	if  extension == "gz" or extension == "gzip": # if compressed
		try:
			with gzip.open(filename) as f:
				      data = f.read().decode('UTF-8')
				      my_splitlines = data.splitlines()
				      f.close()
				      return my_splitlines
		except FileNotFoundError:
			print(f"The file {filename} was not found !")
			
	else: # if considered not compressed
		try:
			with open(filename) as f:
				      data = f.read()
				      my_splitlines = data.splitlines()
				      f.close()
				      return my_splitlines
		except FileNotFoundError:
			print(f"The file {filename} was not found !")

# get run parameters from options file		
def get_Options(OptionsDic, AllowedPars):
	# Parse of command line options and provides help to user
	parser = argparse.ArgumentParser(description='Process Pacbio CCSs ', prefix_chars='-+')
	parser.add_argument("-op", help="the path to the options file.")

	#actually parse command line options
	args = parser.parse_args()

	if args.op is None: # for command line options
		print('Please provide an option file')
		quit()
	else: # for options file
		try:
			OptFN = args.op
			optionsfile=open(OptFN, "r")
			for line in optionsfile :
				if line[:1] != "#" and len(line.strip("\n".strip(""))) > 1 :
					splitted_line = line.split(":")
					key=splitted_line[0].split(" ")[0]
					value=splitted_line[1].strip("\n").split(" ")[0]
					if key in AllowedPars :
						OptionsDic[key] = value
		except:
			print('Error reading %s file !' % (OptFN))
			exit(1)
			
		return OptionsDic

# get run parameters from options file (This function was adapted from "get_Options" to be used from jupyter notebooks)
def get_OptionsIPyNB(OptFN, AllowedPars):
	OptionsDic = {}
	
	try:
		optionsfile=open(OptFN, "r")
		for line in optionsfile :
			if line[:1] != "#" and len(line.strip("\n".strip(""))) > 1 :
				splitted_line = line.split(":")
				key=splitted_line[0].split(" ")[0]
				value=splitted_line[1].strip("\n").split(" ")[0]
				if key in AllowedPars :
					OptionsDic[key] = value
	except:
		print('Error reading %s file !' % (OptFN))
		exit(1)
		
	return OptionsDic


# get reference sequence and related annotations
def get_Refs(RefsFN):
	RefsDic = {}
	Refs = get_file_content(RefsFN)
	Ref_name = None
	
	for line in Refs:
		splitted_line = line.split(":")
		if len(splitted_line) > 1 and line[0] == ">" :
			Ref_name = splitted_line[0].strip(">") # key is the reference name
			CDS = splitted_line[1].split("-") 
			Anchors = splitted_line[2].split("-") 
			BC = splitted_line[3].split("-") 
			Idx = splitted_line[4].split("-") 
			''''' RefsDic Structure: key is the reference name and value is a list of 8 variables
			RefsDic[Ref_name] = [CDS positions, Anchors sequences, BC positions, Idx positions, Full ref sequene, CDS sequence, BC sequence, Index sequence]
			field 0 (CDS) comprise the start and end position of the fusion partner
			field 1 (Anchros) are the 5' and 3' anchor sequences flanking the BC 
			field 2 (BC) comprise the start and end position of the BC
			field 3 (Index) comprise the start and end position of the Index (for sample demultiplexing)
			field 4 is reserved for full reference sequence
			field 5 is for CDS
			field 6 is for BC sequence
			field 7 is for index sequence
			'''''
			RefsDic[Ref_name] = [CDS,Anchors,BC,Idx,"", "","",""]
		else:
			line = line.strip("\n").strip(" ")
			if len(line) > 1 and re.match(r'[ACGTWSMKRYBDHVNacgtwsmkrybdhvnacgtwsmkrybdhvn]', line, re.M|re.I): # get full reference sequence
				RefsDic[Ref_name][4] = RefsDic[Ref_name][4] + line
			if len(RefsDic[Ref_name][4]) >= int(RefsDic[Ref_name][0][1]): # get the CDS sequence
				CDSseq = RefsDic[Ref_name][4][int( RefsDic[Ref_name][0][0]) -1 : int(RefsDic[Ref_name][0][1])]
				RefsDic[Ref_name][5] = CDSseq
			if len(RefsDic[Ref_name][4]) >= int(RefsDic[Ref_name][2][1]): # get the BC sequence
				BC_seq = RefsDic[Ref_name][4][int( RefsDic[Ref_name][2][0]) -1 : int(RefsDic[Ref_name][2][1])]
				RefsDic[Ref_name][6] = BC_seq
			if len(RefsDic[Ref_name][4]) >= int(RefsDic[Ref_name][3][1]): # get the index sequence
				Idx_seq = RefsDic[Ref_name][4][int( RefsDic[Ref_name][3][0]) -1 : int(RefsDic[Ref_name][3][1])]
				RefsDic[Ref_name][7] = Idx_seq
				
	return RefsDic
	
# store data as pickle
def data2pickle(data, Ccomplex, FN):
	data_slice = {}
	data_slice[Ccomplex] = data[Ccomplex]
	with open(FN, 'wb') as handle:
		pickle.dump(data_slice, handle, protocol=pickle.HIGHEST_PROTOCOL)
	handle.close()

# store data as pickle
def data2pickleLR(data, FN):
	data_slice = data
	with open(FN, 'wb') as handle:
		pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)
	handle.close()

# clean intermediary files given in a list
def clean(IntFNList, LogFile):
		FileList = "* ".join(IntFNList)
		# run Nucmer using CDS of each reference and the CCS that were validated for the same reference 
		cmd = (f"rm {FileList}*  2>/dev/null") # Redirect errors to /dev/null; Nucmer version 4rc1
		run(cmd, LogFile)
		print("\nIntermediary files deleted.")
		
# remove outliers
def remove_outiers(MyList,ZScoreTreshold):
	MyList2 = []
	ZScoresL = stats.zscore(MyList)
	if sum(i < abs(ZScoreTreshold) for i in ZScoresL) >=3 :
		# remove outliers
		for n in range(len(ZScoresL)):
			if abs(ZScoresL[n]) > ZScoreTreshold:
				continue
			else:
				MyList2.append(MyList[n])
	else:
		MyList2 = MyList 
	
	return MyList2
