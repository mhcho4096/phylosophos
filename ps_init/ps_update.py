#!/usr/bin/env python

################################################################
#### 
#### PROJECT PHYLOSOPHOS: FORMAL VERSION
#### 
#### INITIALIZE MODULE - REFERENCE DATASET UPDATE SUBMODULE
#### 
#### ORIGINAL SCRIPT WRITTEN BY MIN HYUNG CHO, PH.D.
#### 
#### BIOINFORMATICS AND MOLECULAR DESIGN RESEARCH CENTER
#### 
################################################################

#### Library import

import numpy
import os
import ssl
import sys
import tarfile
import zipfile

from urllib import request

#### Function definition

def ncbi_tax_update_new(ustat, fdir, rdir):

	# Download step

	context = ssl._create_unverified_context()

	tar_ncbi_url = "https://ftp.ncbi.nih.gov/pub/taxonomy/"
	tar_ncbi_file = "taxdump.tar.gz"
	fdlist = os.listdir(fdir)
	if ustat >= 1:
		request.urlretrieve(tar_ncbi_url+tar_ncbi_file, fdir+tar_ncbi_file, context=context)
		print("## NCBI taxonomy updated file download completed")
	else:
		if tar_ncbi_file not in fdlist:
			request.urlretrieve(tar_ncbi_url+tar_ncbi_file, fdir+tar_ncbi_file, context=context)
			print("## NCBI taxonomy file download completed")
		else:
			print("## NCBI taxonomy file found in local system")

	# Unzip & processing step

	try:
		os.makedirs(fdir+"ncbi_taxdump")
	except:
		pass
	ncbi_tar = tarfile.open(fdir+tar_ncbi_file, "r:gz")
	ncbi_tar.extractall(path=fdir+"ncbi_taxdump\\")
	ncbi_tar.close()
	print("## NCBI taxonomy file unzip successful")

	# Extraction step

	taxid_names_dict = {}
	with open(fdir+"ncbi_taxdump\\names.dmp") as inp_f:
		for line in inp_f:
			ssl = [j.replace('\t', '') for j in line.rstrip('\n').split('|')]
			if ssl[0] not in taxid_names_dict:
				taxid_names_dict[ssl[0]] = ["", []]
			#
			if ssl[3] == "scientific name":
				if taxid_names_dict[ssl[0]][0] == "":
					taxid_names_dict[ssl[0]][0] = ssl[1]
			elif ssl[3] in ["acronym", "blast name", "common name", "equivalent name", "genbank common name", "genbank synonym", "synonym"]:
					taxid_names_dict[ssl[0]][1].append(ssl[1])
	print("## NCBI taxonomy name extraction completed", len(taxid_names_dict))
	taxid_parent_dict = {}
	taxid_level_code = {"superkingdom":1, "kingdom":2, "phylum":3, "class":4, "order":5, "family":6, "genus":7, "species":8}
	with open(fdir+"ncbi_taxdump\\nodes.dmp") as inp_f:
		for line in inp_f:
			ssl = [j.replace('\t', '') for j in line.rstrip('\n').split('|')]
			if ssl[2] in taxid_level_code:
				taxid_parent_dict[ssl[0]] = [ssl[1], ssl[2], taxid_level_code[ssl[2]]]
			else:
				taxid_parent_dict[ssl[0]] = [ssl[1], ssl[2], 0]
	print("## NCBI taxonomy node extraction completed", len(taxid_parent_dict))

	# Export step #1: address mapping

	taxid_list = list(taxid_parent_dict.keys())
	taxid_res = []
	genus_dict = {}

	for i_1 in range(len(taxid_list)):
		taxid_code = [[taxid_list[i_1]], [taxid_parent_dict[taxid_list[i_1]][2]]]
		t_det = 0
		while t_det == 0:
			tar_parent = taxid_parent_dict[taxid_code[0][-1]]
			taxid_code[0].append(tar_parent[0])
			taxid_code[1].append(tar_parent[2])
			if tar_parent[0] == "1":
				t_det += 1
		taxid_code[1].append(0)
		taxid_res.append([taxid_list[i_1], taxid_names_dict[taxid_list[i_1]][0], '|'.join(taxid_names_dict[taxid_list[i_1]][1]), taxid_parent_dict[taxid_list[i_1]][1], '|'.join(taxid_code[0]), '|'.join([str(j) for j in taxid_code[1][1:]])])
		if max(taxid_code[1]) >= 7:
			taxid_parse = [j.lower() for j in taxid_names_dict[taxid_list[i_1]][0].split()]
			if taxid_parse[0] not in genus_dict:
				genus_dict[taxid_parse[0]] = []
			genus_dict[taxid_parse[0]].append(taxid_list[i_1])
			if len(taxid_names_dict[taxid_list[i_1]][1]) >= 1:
				for i_2 in range(len(taxid_names_dict[taxid_list[i_1]][1])):
					taxid_px = [j.lower() for j in taxid_names_dict[taxid_list[i_1]][1][i_2].split()]
					if taxid_px[0] not in genus_dict:
						genus_dict[taxid_px[0]] = []
					genus_dict[taxid_px[0]].append(taxid_list[i_1])				
		print(i_1, taxid_list[i_1], taxid_names_dict[taxid_list[i_1]][0][0:min(20, len(taxid_names_dict[taxid_list[i_1]][0]))], "                        ", end='\r')

	print("## NCBI taxonomy parent mapping completed // file exported")

	with open(rdir+"ncbi_node_dict.txt", 'w') as res_f:
		for i_1 in range(len(taxid_res)):
			res_f.write('\t'.join(taxid_res[i_1]) + '\n')

	genus_list = list(genus_dict.keys())
	with open(rdir+"ncbi_genus_dict.txt", 'w') as res_f:
		for i_1 in range(len(genus_list)):
			res_f.write(str(i_1) + '\t')
			res_f.write(genus_list[i_1] + '\t')
			res_f.write('|'.join(list(numpy.unique(genus_dict[genus_list[i_1]]))) + '\n')

	print("## NCBI taxonomy genus mapping completed // file exported", len(genus_list))

	#

	print("#### NCBI taxonomy file preprocessing completed ####")

#

def col_tax_update_new(ustat, fdir, rdir, f_vname):

	# Download step

	context = ssl._create_unverified_context()

	tar_col_url = "http://download.catalogueoflife.org/col/"
	tar_col_file = f_vname
	fdlist = os.listdir(fdir)
	if ustat >= 1:
		request.urlretrieve(tar_col_url+tar_col_file, fdir+tar_col_file, context=context)
		print("## CoL taxonomy updated file download completed")
	else:
		if tar_col_file not in fdlist:
			request.urlretrieve(tar_col_url+tar_col_file, fdir+tar_col_file, context=context)
			print("## CoL taxonomy file download completed")
		else:
			print("## CoL taxonomy file found in local system")

	# Unzip & processing step

	try:
		os.makedirs(fdir+"col_taxdump")
	except:
		pass
	col_tar = zipfile.ZipFile(fdir+tar_col_file)
	col_tar.extractall(path=fdir+"col_taxdump\\")
	col_tar.close()
	print("## CoL taxonomy file unzip successful")

	# Extraction step

	taxid_dict = {}
	taxid_level_code = {"domain":1, "kingdom":2, "phylum":3, "class":4, "order":5, "family":6, "genus":7, "species":8}
	tax_types = []

	with open(fdir+"col_taxdump\\Taxon.tsv", errors="ignore") as inp_f:
		inp_f.readline()
		for line in inp_f:
			ssl = line.rstrip('\n').split('\t')
			clade_code = 0
			tax_types.append(ssl[6])
			if ssl[7] in taxid_level_code:
				clade_code += taxid_level_code[ssl[7]]
			#
			if ssl[0] not in taxid_dict and ssl[6] in ["accepted" , "provisionally accepted"]:
				taxid_dict[ssl[0]] = [ssl[0], ssl[1], clade_code, "", [], ssl[7]]
			#
			if ssl[6] in ["accepted" , "provisionally accepted"]:
				if len(ssl[11]) == 0:
					taxid_dict[ssl[0]][3] += ssl[8]
				else:
					t_temp = [ssl[11]]
					t_temp.extend(ssl[13:16])
					taxid_dict[ssl[0]][3] += ' '.join([j for j in t_temp if len(j) >= 1])
	print("## CoL taxonomic accepted data extraction completed", len(taxid_dict))
	with open(fdir+"col_taxdump\\Taxon.tsv", errors="ignore") as inp_f:
		inp_f.readline()
		for line in inp_f:
			ssl = line.rstrip('\n').split('\t')
			if ssl[6] == "synonym":
				if ssl[2] in taxid_dict:
					if len(ssl[11]) == 0:
						taxid_dict[ssl[2]][4].append(ssl[8])
					else:
						t_temp = [ssl[11]]
						t_temp.extend(ssl[13:16])
						taxid_dict[ssl[2]][4].append(' '.join([j for j in t_temp if len(j) >= 1]))
	print("## CoL taxonomic information extraction completed", len(taxid_dict))

	# Export step
	
	taxid_list = list(taxid_dict.keys())
	taxid_res = []
	genus_dict = {}

	for i_1 in range(len(taxid_list)):
		taxid_code = [[taxid_list[i_1]], [taxid_dict[taxid_list[i_1]][2]]]
		t_det = 0
		while t_det < 1:
			if len(taxid_dict[taxid_code[0][-1]][1]) < 1:
				t_det += 1
				break
			target_tax = taxid_dict[taxid_dict[taxid_code[0][-1]][1]]
			taxid_code[0].append(target_tax[0])
			taxid_code[1].append(target_tax[2])
			if len(target_tax[1]) < 1:
				t_det += 1
		taxid_code[1][-1] = 1
		#
		if max(taxid_code[1]) <= 7:
			taxid_res.append([taxid_list[i_1], taxid_dict[taxid_list[i_1]][3].split()[0], '|'.join(taxid_dict[taxid_list[i_1]][4]), taxid_dict[taxid_list[i_1]][5], '|'.join(taxid_code[0]), '|'.join([str(j) for j in taxid_code[1]])])
		elif taxid_dict[taxid_list[i_1]][5] == "species" and " " in taxid_dict[taxid_list[i_1]][3] and taxid_dict[taxid_list[i_1]][3][0].isupper() == True:
			taxid_res.append([taxid_list[i_1], taxid_dict[taxid_list[i_1]][3].split()[0] + " " + taxid_dict[taxid_list[i_1]][3].split()[1], '|'.join(taxid_dict[taxid_list[i_1]][4]), taxid_dict[taxid_list[i_1]][5], '|'.join(taxid_code[0]), '|'.join([str(j) for j in taxid_code[1]])])
		else:
			taxid_res.append([taxid_list[i_1], taxid_dict[taxid_list[i_1]][3], '|'.join(taxid_dict[taxid_list[i_1]][4]), taxid_dict[taxid_list[i_1]][5], '|'.join(taxid_code[0]), '|'.join([str(j) for j in taxid_code[1]])])
		#
		if max(taxid_code[1]) >= 7:
			taxid_parse = [j.lower() for j in taxid_dict[taxid_list[i_1]][3].split()]
			if taxid_parse[0] not in genus_dict:
				genus_dict[taxid_parse[0]] = []
			genus_dict[taxid_parse[0]].append(taxid_list[i_1])
			if len(taxid_dict[taxid_list[i_1]][4]) >= 1:
				for i_2 in range(len(taxid_dict[taxid_list[i_1]][4])):
					taxid_px = [j.lower() for j in taxid_dict[taxid_list[i_1]][4][i_2].split()]
					if taxid_px[0] not in genus_dict:
						genus_dict[taxid_px[0]] = []
					genus_dict[taxid_px[0]].append(taxid_list[i_1])
		print(i_1, taxid_list[i_1], taxid_dict[taxid_list[i_1]][3][0:min(20, len(taxid_dict[taxid_list[i_1]][3]))], "                        ", end='\r')

	print("## CoL taxonomy parent mapping completed // file exported")

	with open(rdir+"col_node_dict.txt", 'w') as res_f:
		for i_1 in range(len(taxid_res)):
			res_f.write('\t'.join(taxid_res[i_1]) + '\n')

	genus_list = list(genus_dict.keys())
	with open(rdir+"col_genus_dict.txt", 'w') as res_f:
		for i_1 in range(len(genus_list)):
			res_f.write(str(i_1) + '\t')
			res_f.write(genus_list[i_1] + '\t')
			res_f.write('|'.join(genus_dict[genus_list[i_1]]) + '\n')

	print("## CoL taxonomy genus mapping completed // file exported", len(genus_list))

	#

	print("#### CoL taxonomy file preprocessing completed ####")

#

def eol_tax_update_new(ustat, fdir, rdir, f_vname):

	# Download step

	context = ssl._create_unverified_context()

	tar_eol_url = "https://opendata.eol.org/dataset/0a023d9a-f8c3-4c80-a8d1-1702475cda18/resource/00adb47b-57ed-4f6b-8f66-83bfdb5120e8/download/"
	tar_eol_file = f_vname
	fdlist = os.listdir(fdir)
	if ustat >= 1:
		request.urlretrieve(tar_eol_url+tar_eol_file, fdir+tar_eol_file, context=context)
		print("## EoL taxonomy updated file download completed")
	else:
		if tar_eol_file not in fdlist:
			request.urlretrieve(tar_eol_url+tar_eol_file, fdir+tar_eol_file, context=context)
			print("## EoL taxonomy file download completed")
		else:
			print("## EoL taxonomy file found in local system")

	# Unzip & processing step

	try:
		os.makedirs(fdir+"eol_taxdump")
	except:
		pass
	
	if tar_eol_file[-3:] == ".gz":
		eol_tar = tarfile.open(fdir+tar_eol_file, "r:gz")
		eol_tar.extractall(path=fdir+"eol_taxdump\\")
		eol_tar.close()
	elif tar_eol_file[-3:] == "zip":
		with zipfile.ZipFile(fdir+tar_eol_file) as zip_ref:
			zip_ref.extractall(path=fdir+"eol_taxdump\\")

	print("## EoL taxonomy file unzip successful")

	# Extraction step

	taxid_dict = {}
	taxid_level_code = {"domain":1, "kingdom":2, "phylum":3, "class":4, "order":5, "family":6, "genus":7, "species":8}

	with open(fdir+"eol_taxdump\\taxon.tab", errors="ignore") as inp_f:
		inp_f.readline()
		for line in inp_f:
			ssl = line.rstrip('\n').split('\t')
			tar_code = ssl[0]
			parent_code = ssl[4]
			clade_code = 0
			if ssl[6] in taxid_level_code:
				clade_code += taxid_level_code[ssl[6]]
			#
			if tar_code not in taxid_dict and ssl[7] in ["valid", "accepted"]:
				taxid_dict[tar_code] = [tar_code, parent_code, clade_code, "", [], ssl[6]]
			#
			if ssl[7] in ["valid", "accepted"]:
				taxid_dict[tar_code][3] = ssl[5]
				if ssl[9] != ssl[5] and len(ssl[9]) >= 1:
					taxid_dict[tar_code][4].append(ssl[9])
			elif ssl[7] in ["synonym", "ambiguous synonym"]:
				taxid_dict[ssl[3]][4].append(ssl[5])
				if ssl[9] != ssl[5] and ssl[9] != taxid_dict[ssl[3]][3] and len(ssl[9]) >= 1:
					taxid_dict[ssl[3]][4].append(ssl[9])
	print("## EoL taxonomic information extraction completed", len(taxid_dict))

	# Export step
	
	taxid_list = list(taxid_dict.keys())
	taxid_res = []
	genus_dict = {}

	for i_1 in range(len(taxid_list)):
		taxid_code = [[taxid_list[i_1]], [taxid_dict[taxid_list[i_1]][2]]]
		t_det = 0
		while t_det < 1:
			if len(taxid_dict[taxid_code[0][-1]][1]) < 1:
				t_det += 1
				break
			target_tax = taxid_dict[taxid_dict[taxid_code[0][-1]][1]]
			taxid_code[0].append(target_tax[0])
			taxid_code[1].append(target_tax[2])
			if len(target_tax[1]) < 1:
				t_det += 1
		taxid_code[1].append(0)
		if max(taxid_code[1]) <= 7:
			taxid_res.append([taxid_list[i_1], taxid_dict[taxid_list[i_1]][3].split()[0], '|'.join(taxid_dict[taxid_list[i_1]][4]), taxid_dict[taxid_list[i_1]][5], '|'.join(taxid_code[0]), '|'.join([str(j) for j in taxid_code[1][:-1]])])
		elif taxid_dict[taxid_list[i_1]][5] == "species" and " " in taxid_dict[taxid_list[i_1]][3] and taxid_dict[taxid_list[i_1]][3][0].isupper() == True:
			taxid_res.append([taxid_list[i_1], taxid_dict[taxid_list[i_1]][3].split()[0] + " " + taxid_dict[taxid_list[i_1]][3].split()[1], '|'.join(taxid_dict[taxid_list[i_1]][4]), taxid_dict[taxid_list[i_1]][5], '|'.join(taxid_code[0]), '|'.join([str(j) for j in taxid_code[1][:-1]])])
		else:
			taxid_res.append([taxid_list[i_1], taxid_dict[taxid_list[i_1]][3], '|'.join(taxid_dict[taxid_list[i_1]][4]), taxid_dict[taxid_list[i_1]][5], '|'.join(taxid_code[0]), '|'.join([str(j) for j in taxid_code[1][:-1]])])
		if max(taxid_code[1]) >= 7:
			taxid_parse = [j.lower() for j in taxid_dict[taxid_list[i_1]][3].split()]
			if taxid_parse[0] not in genus_dict:
				genus_dict[taxid_parse[0]] = []
			genus_dict[taxid_parse[0]].append(taxid_list[i_1])
			if len(taxid_dict[taxid_list[i_1]][4]) >= 1:
				for i_2 in range(len(taxid_dict[taxid_list[i_1]][4])):
					taxid_px = [j.lower() for j in taxid_dict[taxid_list[i_1]][4][i_2].split()]
					if taxid_px[0] not in genus_dict:
						genus_dict[taxid_px[0]] = []
					genus_dict[taxid_px[0]].append(taxid_list[i_1])
		print(i_1, taxid_list[i_1], taxid_dict[taxid_list[i_1]][3][0:min(20, len(taxid_dict[taxid_list[i_1]][3]))], "                        ", end='\r')

	print("## EoL taxonomy parent mapping completed // file exported")

	with open(rdir+"eol_node_dict.txt", 'w') as res_f:
		for i_1 in range(len(taxid_res)):
			res_f.write('\t'.join(taxid_res[i_1]) + '\n')

	genus_list = list(genus_dict.keys())
	with open(rdir+"eol_genus_dict.txt", 'w') as res_f:
		for i_1 in range(len(genus_list)):
			res_f.write(str(i_1) + '\t')
			res_f.write(genus_list[i_1] + '\t')
			res_f.write('|'.join(genus_dict[genus_list[i_1]]) + '\n')

	print("## EoL taxonomy genus mapping completed // file exported", len(genus_list))

	#

	print("#### EoL taxonomy file preprocessing completed ####")

####
