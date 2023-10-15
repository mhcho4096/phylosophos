#!/usr/bin/env python

################################################################
#### 
#### PROJECT PHYLOSOPHOS: FORMAL VERSION
#### 
#### INITIALIZE MODULE - DATA IMPORT SUBMODULE
#### 
#### ORIGINAL SCRIPT WRITTEN BY MIN HYUNG CHO, PH.D.
#### 
#### BIOINFORMATICS AND MOLECULAR DESIGN RESEARCH CENTER
#### 
################################################################

#### Library import

import datetime
import gc
import itertools
import numpy
import os
import sys

#### Preset parameters

alpha_numeric = {"a":0, "b":1, "c":2, "d":3, "e":4, "f":5, "g":6, "h":7, "i":8, "j":9, 
"k":10, "l":11, "m":12, "n":13, "o":14, "p":15, "q":16, "r":17, "s":18, "t":19, 
"u":20, "v":21, "w":22, "x":23, "y":24, "z":25}

#### Function definition

def phylosophos_initialization(arg_list):

	# STEP 0: USER GUIDE PRINT

	if "-h" in arg_list or "-help" in arg_list:
		phylosophos_help()
		exit()

	# STEP 1: REFERENCE DATASET SELECTION

	print("#### PhyloSophos initialization started ####")

	ref_dir = os.getcwd()+"\\pp_ref\\"
	ref_pool = [j.split("_")[0] for j in os.listdir(ref_dir) if j[-4:] == ".txt"]
	default_reference = "ncbi" # Change default setting as appropriate
	ref_type = ""
	ref_stat = 0
	if "-ref" in arg_list or "-r" in arg_list:
		ref_det = 0
		for pi_1 in range(len(arg_list)-1):
			if arg_list[pi_1].lower() in ["-ref", "-r"]:
				ref_det += 1
				if arg_list[pi_1+1].lower() in ref_pool:
					ref_type += arg_list[pi_1+1]	
				else:
					ref_type += default_reference
					ref_stat += 1
				break
		if ref_det < 1:
			ref_type += default_reference
	else:
		ref_type += default_reference
	if ref_stat == 0:
		print("## Reference type:", ref_type)
	elif ref_stat == 1:
		print("## Reference type: ncbi (reference not found: default option applied)")

	# STEP 2: INPUT TYPE SELECTION

	input_type = 0 # 0: all files in the input directory, 1: individual file
	input_file_name = ""
	if "-input" in arg_list or "-i" in arg_list:
		for pi_2 in range(len(arg_list)-1):
			if arg_list[pi_2].lower() in ["-input", "-i"]:
				input_type = 1
				input_file_name += arg_list[pi_2+1]
	if input_type == 0:
		print("## Input type: default directory")
	elif input_type == 1:
		print("## Input type: single file //", input_file_name)

	# STEP 3: LEVENSHTEIN DISTANCE CUTOFF

	default_cutoff = 3
	if "-lev" in arg_list or "-l" in arg_list or "-cutoff" in arg_list:
		for pi_3 in range(len(arg_list)-1):
			if arg_list[pi_3].lower() in ["-lev", "-l", "-cutoff"]:
				default_cutoff = int(arg_list[pi_3+1])
	print("## Levenshtein cutoff:", default_cutoff)

	# STEP 4: MANUAL CURATION USAGE

	mc_stat = False
	if "-m" in arg_list or "-manual" in arg_list or "-curation" in arg_list:
		for pi_3 in range(len(arg_list)-1):
			if arg_list[pi_3].lower() in ["-m", "-manual", "-curation"]:
				mc_stat = bool(int(arg_list[pi_3+1]))
				break
	print("## Manual curation data usage:", mc_stat)

	print("#### PhyloSophos analysis started ####")

	# RETURN PARAMETERS

	return ref_type, input_type, input_file_name, default_cutoff, mc_stat

def phylosophos_help():

	ref_dir = os.getcwd()+"\\pp_ref\\"
	ref_pool = list(numpy.unique([j.split("_")[0] for j in os.listdir(ref_dir) if j[-4:] == ".txt"]))

	#

	print("#### PHYLOSOPHOS TAXONOMIC MAPPING ALGORITHM ####")
	print("##")
	print("## DESIGNED AND IMPLEMENTED BY MIN HYUNG CHO, PH.D.")
	print("")
	print("[1] IN A NUTSHELL")
	print("")
	print("> phylosophos_initialize_update.py [OPTIONAL_UPDATE_PARAMETER]")
	print("> phylosophos_core.py [[OPTIONAL_PARAMETER_TYPE] [OPTIONAL_PARAMETER_VALUE]]")
	print("")
	print("[2] OPTIONAL CORE PARAMETER TYPES")
	print("")
	print("* HELP (-h, -help, -guide): PRINTS SHORT GUIDE (AS YOU CAN SEE HERE)")
	print("")
	print("* REFERENCE (-r, -ref): REFERENCE TYPE CHANGE (DEFAULT = 'ncbi')")
	print("    CURRENTLY AVAILABLE REFERENCES: "+str(ref_pool))
	print("")
	print("* INPUT (-i, -input): INPUT TYPE MODIFICATION")
	print("    (DEFAULT SETTING: IMPORT EVERY FILES WITHIN INPUT DIRECTORY AS SCIENTIFIC NAME INPUTS)")
	print("    (IF FILE PATH IS GIVEN: IMPORT SINGLE DATA FILE AS A SCIENTIFIC NAME INPUT)")
	print("")
	print("* DL DISTANCE CUTOFF (-l, -lev, -cutoff): DAMERAU-LEVENSHTEIN DISTANCE CUTOFF (DEFAULT = 3)")
	print("")
	print("* MANUAL CURATION (-m, -manual, -curation): APPLICATION OF MANUALLY CURATED DATA (DEFAULT = FALSE)")
	print("")
	print("[3] OPTIONAL UPDATE PARAMETER")
	print("")
	print("* REFERENCE RAW DATA UPDATE: FORCED TAXONOMIC METADATA DOWNLOAD & UPDATE")
	print("    (DEFAULT SETTING = 0: DOWNLOAD STEP IS SKIPPED IF CORRESPONDING FILES EXIST IN /EXTERNAL_FILES DIRECTORY)")
	print("    (IF INTEGER VALUE OTHER THAN 0 IS GIVEN: DOWNLOAD REFERENCE METADATA AND OVERRIDE PRE-EXISTING LOCAL FILES)")
	print("")
	print("## FOR DETAILED INFORMATION, PLEASE REFER TO PHYLOSOPHOS_GUIDE.PDF FILE")

##

def phylosophos_input_import(i_type, i_fname):

	input_dir = os.getcwd()+"\\input\\"
	pii_input_list = []

	print("#### Input file import started ####", end='\r')

	# IMPORT FROM INPUT DIRECTORY

	if i_type == 0:		
		for pii_1 in os.listdir(input_dir):
			pii_sub = []
			with open(input_dir+pii_1, encoding = 'utf-8') as inp_f:
				for line in inp_f:
					pii_sub.append(line.rstrip('\n'))
			pii_input_list.append([pii_1, pii_sub])

	# IMPORT FROM SPECIFIC FILE

	elif i_type == 1:
		if i_fname in os.listdir(input_dir):
			pii_sub = []
			with open(input_dir+i_fname, encoding = 'utf-8') as inp_f:
				for line in inp_f:
					pii_sub.append(line.rstrip('\n'))
			pii_input_list.append([i_fname, pii_sub])
		else:
			try:
				pii_sub = []
				with open(i_fname, encoding = 'utf-8') as inp_f:
					for line in inp_f:
						pii_sub.append(line.rstrip('\n'))
				pii_input_list.append([i_fname, pii_sub])
			except:
				print("- ERROR OCCURRED: INPUT FILE NOT FOUND")
				print("- ANALYSIS TERMINATED")
				exit()

	print("#### Input file import completed ####")

	return pii_input_list

##

def alpha_count(string):
	lst = string.lower()
	tar_code = [0 for j in range(27)]
	for ac_1 in range(len(lst)):
		if lst[ac_1] in alpha_numeric:
			tar_code[alpha_numeric[lst[ac_1]]] += 1
		else:
			tar_code[26] += 1
	return tar_code

def ref_partial_correction(string):
	#
	ignore_list = ['sp.', 'ssp.', 'genomosp.', 'genosp.', 'subsp.', 'var.', 'str.', 'f.', 'pv.', 'bv.', 's.', 's.l.', 
	'al.', 'sect.', 'subgen.', 'nom.', 'no.', "species"]
	dropout_list = ['cf.', 'aff.', 'nr.', 'n.', 's.n.', 'nov.', 'gen.', 'inval.']
	hybrid_list = ['×', 'x', 'x']
	#
	stat_code = [0, 0, 0]
	tar_temp = []
	#
	ssl = string.lower().split()
	for ic_1 in ssl:
		if ic_1 in ignore_list:
			stat_code[0] = 1
		elif ic_1 in dropout_list:
			stat_code[1] = 1
		elif ic_1 in hybrid_list:
			stat_code[2] = 1
		elif '.' in ic_1 or ',' in ic_1:
			if len(ssl) >= 2:
				continue
		elif '<' in ic_1 or '>' in ic_1:
			tar_temp.append(ic_1.split('>')[0].split('<')[0])
		else:
			tar_temp.append(ic_1)
	#
	return ' '.join(tar_temp), stat_code

def phylosophos_ref_import(ref_path):

	print("#### Reference file import started ####", end='\r')
	time_start = datetime.datetime.now()

	tar_ref = sorted(os.listdir(ref_path))
	ref_types = []
	for ri_1 in tar_ref:
		if ri_1[-4:] == ".txt":
			ref_types.append(ri_1.split("_")[0])
	ref_types = list(sorted(numpy.unique(ref_types)))

	raw_names_list = {}
	whole_names_list = {}
	whole_genus_list = {}

	for ri_2 in range(len(ref_types)):
		n_ref_file = ref_types[ri_2]+"_node_dict.txt"
		g_ref_file = ref_types[ri_2]+"_genus_dict.txt"
		#
		raw_names_list[ref_types[ri_2]] = {}
		whole_names_list[ref_types[ri_2]] = {}
		whole_genus_list[ref_types[ri_2]] = {}
		canon_name_list = {}
		#
		with open(ref_path+n_ref_file) as inp_f:
			for line in inp_f:
				ssl = line.rstrip('\n').split('\t')
				whole_names_list[ref_types[ri_2]][ssl[0]] = ssl[0:6]
				#
				canon_name_list[ssl[1].lower()] = ssl[0]
				if ssl[1].lower() not in raw_names_list[ref_types[ri_2]]:
					raw_names_list[ref_types[ri_2]][ssl[1].lower()] = []
				raw_names_list[ref_types[ri_2]][ssl[1].lower()].append(ssl[0])
		#
		with open(ref_path+n_ref_file) as inp_f:
			for line in inp_f:
				ssl = line.rstrip('\n').split('\t')
				if len(ssl[2]) >= 1:
					ssl_2 = [j for j in ssl[2].lower().split("|") if len(j) >= 1]
					for ri_3 in ssl_2:
						if ri_3 not in canon_name_list:
							if ri_3 not in raw_names_list[ref_types[ri_2]]:
								raw_names_list[ref_types[ri_2]][ri_3] = []
							raw_names_list[ref_types[ri_2]][ri_3].append(ssl[0])
				if len(ssl[1].split()) >= 4:
					corr_string, corr_stat = ref_partial_correction(ssl[1].lower())
					if corr_stat[0] == 1 and sum(corr_stat) < 2 and corr_string not in raw_names_list and ssl[1].split()[1].lower() != 'sp.':
						raw_names_list[ref_types[ri_2]][corr_string] = [ssl[0]]
		#
		with open(ref_path+g_ref_file) as inp_f:
			for line in inp_f:
				ssl = line.rstrip('\n').split('\t')
				whole_genus_list[ref_types[ri_2]][ssl[1].lower()] = [alpha_count(ssl[1]), ssl[2].split("|")]
		#
		print(ref_types[ri_2], len(whole_genus_list[ref_types[ri_2]]), len(whole_names_list[ref_types[ri_2]]), len(raw_names_list[ref_types[ri_2]]), datetime.datetime.now()-time_start, "                ", end='\r')

	print("                                                                                ", end='\r')
	print("#### Reference file import completed ####")
	return ref_types, whole_names_list, whole_genus_list, raw_names_list

#

