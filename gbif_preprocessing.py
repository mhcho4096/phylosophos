#!/usr/bin/env python

################################################################
#### 
#### PROJECT PHYLOSOPHOS: ACCESSORY SCRIPT
#### 
#### EXAMPLE ADDITIONAL DB PREPARATION SCRIPT // GBIF-BASED
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

####

base_path = os.path.realpath(__file__)[:-21]
ef_path = "external_files\\"
pp_path = "pp_ref\\"

ssl._create_default_https_context = ssl._create_unverified_context

####

def gbif_tax_update(ustat, fdir, rdir):

	# Download step	

	tar_gbif_url = "https://hosted-datasets.gbif.org/datasets/backbone/current/"
	tar_gbif_file = "backbone.zip"
	fdlist = os.listdir(fdir)
	if ustat >= 1:
		request.urlretrieve(tar_gbif_url+tar_gbif_file, fdir+tar_gbif_file)
		print("## GBIF taxonomy updated file download completed")
	else:
		if tar_gbif_file not in fdlist:
			request.urlretrieve(tar_gbif_url+tar_gbif_file, fdir+tar_gbif_file)
			print("## GBIF taxonomy file download completed")
		else:
			print("## GBIF taxonomy file found in local system")

	# Unzip & processing step

	try:
		os.makedirs(fdir+"gbif_taxdump")
	except:
		pass
	gbif_tar = zipfile.ZipFile(fdir+tar_gbif_file)
	gbif_tar.extractall(path=fdir+"gbif_taxdump\\")
	gbif_tar.close()
	print("## GBIF taxonomy file unzip successful")

	# Extraction step

	taxid_raw = {}
	accepted_names = {}

	with open(fdir+"gbif_taxdump\\backbone\\"+"Taxon.tsv", encoding="UTF-8", errors='ignore') as inp_f:
		inp_f.readline()
		for line in inp_f:
			ssl = line.rstrip('\n').split('\t')
			if ssl[0] == "0":
				taxid_raw["-"] = ["-", ssl[2], ssl[3], "Root", "accepted", "root"]
			#
			if len(ssl[7]) >= 1:
				if ssl[11] in ['species', 'subspecies', 'form', 'variety']:
					if len(ssl[9]) < 1:
						taxid_raw[ssl[0]] = [ssl[0], ssl[2], ssl[3], ssl[5], ssl[14], ssl[11]]
					else:
						taxid_raw[ssl[0]] = [ssl[0], ssl[2], ssl[3], ssl[7], ssl[14], ssl[11]]
				else:
					taxid_raw[ssl[0]] = [ssl[0], ssl[2], ssl[3], ssl[7], ssl[14], ssl[11]]
			else:
				taxid_raw[ssl[0]] = [ssl[0], ssl[2], ssl[3], ssl[5], ssl[14], ssl[11]]
			#
			if ssl[14] == "accepted":
				accepted_names[taxid_raw[ssl[0]][3]] = ssl[0]
			#
			if len(taxid_raw) % 1000 == 0:
				print(ssl[0], ssl[7][0:min(20, len(ssl[7]))], "        ", end='\r')

	print("## GBIF taxonomy node extraction completed", len(taxid_raw), len(accepted_names),  "        ")
	taxid_key = list(taxid_raw.keys())

	#

	taxid_names_dict = {}
	taxid_parent_dict = {}
	taxid_genus_dict = {}
	taxid_temp = {}

	for i_1 in range(len(taxid_key)):
		if len(taxid_raw[taxid_key[i_1]][3]) >= 1:
			taxid_temp[taxid_raw[taxid_key[i_1]][0]] = taxid_raw[taxid_key[i_1]][3]
			if taxid_raw[taxid_key[i_1]][4] == "accepted":
				if taxid_raw[taxid_key[i_1]][0] not in taxid_names_dict:
					taxid_names_dict[taxid_raw[taxid_key[i_1]][0]] = ["", []]
				taxid_names_dict[taxid_raw[taxid_key[i_1]][0]][0] = taxid_raw[taxid_key[i_1]][3]

			elif taxid_raw[taxid_key[i_1]][4] == "doubtful" and taxid_raw[taxid_key[i_1]][3] not in accepted_names:
				if taxid_raw[taxid_key[i_1]][0] not in taxid_names_dict:
					taxid_names_dict[taxid_raw[taxid_key[i_1]][0]] = ["", []]
				taxid_names_dict[taxid_raw[taxid_key[i_1]][0]][0] = taxid_raw[taxid_key[i_1]][3]

			elif len(taxid_raw[taxid_key[i_1]][2]) >= 1:
				if taxid_raw[taxid_key[i_1]][2] not in taxid_names_dict:
					taxid_names_dict[taxid_raw[taxid_key[i_1]][2]] = ["", []]
				taxid_names_dict[taxid_raw[taxid_key[i_1]][2]][1].append(taxid_raw[taxid_key[i_1]][3])
				print(i_1, taxid_key[i_1], "        ", end='\r')

	print("## GBIF taxonomy node classification completed", len(taxid_names_dict), "                                ")
	names_key = list(taxid_names_dict.keys())

	for i_2 in range(len(names_key)):
		#
		if taxid_names_dict[names_key[i_2]][0] == "":
			taxid_names_dict[names_key[i_2]][0] += taxid_temp[names_key[i_2]]
		#
		taxid_level_code = {"root":1, "kingdom":2, "phylum":3, "class":4, "order":5, "family":6, "genus":7, "species":8}
		par_sub = [[taxid_raw[names_key[i_2]][0]], [0]]
		if taxid_raw[names_key[i_2]][5] in taxid_level_code:
			par_sub[1][-1] += taxid_level_code[taxid_raw[names_key[i_2]][5]]
		#
		if names_key[i_2] == "-":
			taxid_parent_dict[names_key[i_2]] = par_sub
			taxid_parent_dict[names_key[i_2]][1].append(1)
			print(i_2, names_key[i_2], taxid_raw[names_key[i_2]][3][:min(40, len(taxid_raw[names_key[i_2]][3]))], "                ", end='\r')
			pass
		else:
			p_det = 0
			while p_det < 1:
				if taxid_raw[par_sub[0][-1]][1] == "":
					par_sub[0].append("-")
					par_sub[1].append(2)
					p_det += 1
				else:
					if taxid_raw[par_sub[0][-1]][5] in taxid_level_code:
						par_sub[1].append(taxid_level_code[taxid_raw[par_sub[0][-1]][5]])
					else:
						par_sub[1].append(0)
					par_sub[0].append(taxid_raw[par_sub[0][-1]][1])
			par_sub[1].append(1)
			taxid_parent_dict[names_key[i_2]] = par_sub
			print(i_2, names_key[i_2], taxid_raw[names_key[i_2]][3][:min(40, len(taxid_raw[names_key[i_2]][3]))], "                ", end='\r')
		#
		if max(par_sub[1]) >= 2:
			ssl = []
			if len(taxid_names_dict[names_key[i_2]][0]) >= 1:
				ssl.append(taxid_names_dict[names_key[i_2]][0].split()[0])
			if len(taxid_names_dict[names_key[i_2]][1]) >= 1:
				ssl.extend([j.split()[0] for j in taxid_names_dict[names_key[i_2]][1]])
			#
			for i_4 in range(len(ssl)):
				if ssl[i_4] not in taxid_genus_dict:
					taxid_genus_dict[ssl[i_4]] = []
				taxid_genus_dict[ssl[i_4]].append(names_key[i_2])

	print("## GBIF taxonomy parent mapping completed",  "        ")

	#

	with open(rdir+"gbif_node_dict.txt", 'w', errors='ignore') as res_f:
		for i_3 in range(len(names_key)):
			#
			taxid_names_dict[names_key[i_3]][1] = list(sorted(numpy.unique(taxid_names_dict[names_key[i_3]][1])))
			#
			res_f.write(names_key[i_3] + '\t')
			res_f.write(taxid_names_dict[names_key[i_3]][0] + '\t')
			res_f.write('|'.join(taxid_names_dict[names_key[i_3]][1]) + '\t')
			res_f.write(taxid_raw[names_key[i_3]][5] + '\t')
			res_f.write('|'.join(taxid_parent_dict[names_key[i_3]][0]) + '\t')
			res_f.write('|'.join([str(j) for j in taxid_parent_dict[names_key[i_3]][1]][1:]) + '\n')

	genus_key = list(taxid_genus_dict.keys())

	with open(rdir+"gbif_genus_dict.txt", 'w', errors='ignore') as res_f:
		for i_3 in range(len(genus_key)):
			res_f.write(str(i_3) + '\t')
			res_f.write(genus_key[i_3] + '\t')
			res_f.write('|'.join(taxid_genus_dict[genus_key[i_3]]) + '\n')

	print("## GBIF taxonomy mapping file export completed",  "        ")

	##

	print("## GBIF taxonomy mapping completed",  "        ")

#### Main part

def main():
	update_stat = int(sys.argv[1])
	gbif_tax_update(update_stat, ef_path, pp_path)

if __name__ == "__main__":
	main()



