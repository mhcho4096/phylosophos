#!/usr/bin/env python

################################################################
#### 
#### PROJECT PHYLOSOPHOS: FORMAL VERSION
#### 
#### ANALYSIS/PROCESSING MODULE - DATA IMPORT SUBMODULE
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

## String correction functions

def lev_dist(string_1, string_2, cut_dist): # Damerau-Levenshtein distance
	# Parameter setup
	dl_xlen = len(string_1)+1
	dl_ylen = len(string_2)+1
	dl_det = 0 # Status determinant for early termination
	# Matrix setup
	dl_matrix = numpy.zeros((dl_xlen, dl_ylen))
	for dl_1 in range(dl_xlen):
		dl_matrix[dl_1, 0] = dl_1
	for dl_2 in range(dl_ylen):
		dl_matrix[0, dl_2] = dl_2
	# Score calculation
	for dl_3 in range(1, dl_xlen):
		for dl_4 in range(1, dl_ylen):
			if string_1[dl_3-1] == string_2[dl_4-1]:
				dl_matrix[dl_3, dl_4] = min(dl_matrix[dl_3-1, dl_4-1], dl_matrix[dl_3-1, dl_4]+1, dl_matrix[dl_3, dl_4-1]+1)
			elif string_1[dl_3-2] == string_2[dl_4-1] and string_1[dl_3-1] == string_2[dl_4-2]:
				dl_matrix[dl_3, dl_4] = min(dl_matrix[dl_3-2, dl_4-2]+1, dl_matrix[dl_3-1, dl_4]+1, dl_matrix[dl_3, dl_4-1]+1)
			else:
				dl_matrix[dl_3, dl_4] = min(dl_matrix[dl_3-1, dl_4-1]+1, dl_matrix[dl_3-1, dl_4]+1, dl_matrix[dl_3, dl_4-1]+1)
			#
			if dl_3 == dl_4 and min(dl_matrix[dl_3-1]) > cut_dist:
				dl_det += 1
				break
		if dl_det >= 1:
			break
	# return value
	if dl_det >= 1:
		return max(dl_xlen, dl_ylen)
	else:
		return int(dl_matrix[dl_xlen-1, dl_ylen-1])

def alpha_count(string):
	lst = string.lower()
	tar_code = [0 for j in range(27)]
	for ac_1 in range(len(lst)):
		if lst[ac_1] in alpha_numeric:
			tar_code[alpha_numeric[lst[ac_1]]] += 1
		else:
			tar_code[26] += 1
	return tar_code

def alpha_compare(list_1, list_2, cutoff):
	count_score = 0
	for i_1 in range(len(list_1)):
		if list_1[i_1] != list_2[i_1]:
			count_score += 1
			if count_score > cutoff:
				break
	return count_score

## Input correction functions

def string_split_correction(string):
	sst = ""
	if len(string) >= 1:
		for ssc_1 in range(len(string)):
			if ssc_1 >= 1:
				if string[ssc_1-1] in [",", ".", ")"] and string[ssc_1] not in [" ", ")"]:
					sst += " "
				elif string[ssc_1] in ["("] and string[ssc_1-1] not in [" "]:
					sst += " "
			#
			if string[ssc_1] == "," and ssc_1 == len(string)-1:
				continue
			if string[ssc_1] == " " and ssc_1 >= 1:
				if string[ssc_1-1] == "-":
					continue
			if string[ssc_1] == " " and ssc_1 < len(string)-1:
				if string[ssc_1+1] in [".", ",", "-"]:
					continue
			#
			sst += string[ssc_1]
	return sst

def input_correction(string):
	#
	ignore_list = ['sp.', 'ssp.', 'genomosp.', 'genosp.', 'subsp.', 'var.', 'str.', 'f.', 'pv.', 'bv.', 's.', 's.l.', 
	'al.', 'sect.', 'subgen.', 'nom.', 'no.', "species", "var."]
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
	tar_str = ' '.join(tar_temp)
	if "[syn." in tar_str:
		tar_str = tar_str.split("[syn.")[0].rstrip(' ')
	#
	return tar_str, stat_code

def latin_correction(tar_psd):
	#
	removed_list = ["bulbus", "carapax", "caulis", "concha", "cortex", "embryo", "exocarpium", ]
	removed_list.extend(["flos", "folium", "fructus", "herba", "lignum", "ligum", "oleum", "ootheca", "pedicellus"])
	removed_list.extend(["pericarpium", "pollen", "radix", "rhizoma", "sanguis", "semen", "strobilus"])
	removed_list.extend(["cum", "et", "radicis", "cornu", "praeparata", "preparata", "succus", "squama"])
	removed_list.extend(["stigma", "oviductus", "acidum", "fel", "stamen", "pulveratum", "resina", "excrementum"])
	removed_list.extend(["periostracum", "plumula", "corium", "dens", "extractum", "liquidum", "fossilia", "gummi"])
	removed_list.extend(["plumula", "massa", "fermentata", "medulla", "nidus", "penis", "corticis", "folii"])
	removed_list.extend(["nodus", "rhizomatis", "plastrum", "ramulus", "testa", "testis", "exodermis", "germinatus"])
	removed_list.extend(["arillus", "alburnum", "cacumen", "os", "petiolus", "pix", "pulvis", "spina", "receptaculum"])
	removed_list.extend(["leaf", "of", "root", "twig", "the", "family", "et", "from"])
	#
	drop_list = ['spp.', 'sp.', "ssp", "spp", "sp", 'var.', 'f.', "L", "L.", "(L.)", "Roxb", "Roxb.", "Fr.", "Fr", "DC", "DC."]
	#
	corr_det = 0
	corr_temp = []
	for i_1 in range(len(tar_psd)):
		if tar_psd[i_1] in removed_list:
			corr_det += 1
		elif i_1 == len(tar_psd)-1 and tar_psd[i_1] in drop_list:
			pass
		elif "(" in tar_psd[i_1] or ")" in tar_psd[i_1]:
			if tar_psd[i_1][0] == "(":
				pass
			else:
				t_1 = tar_psd[i_1].split("(")
				for i_2 in range(len(t_1)):
					t_2 = t_1[i_2].split(")")
					corr_temp.extend([str(j) for j in t_2 if len(j) >= 1])
		elif "." in tar_psd[i_1] and i_1 < 2 and len(tar_psd[i_1]) <= 2:
			pass
		elif len(tar_psd[i_1]) == 0:
			pass
		else:
			corr_temp.append(tar_psd[i_1])
	if corr_det < 1:
		return [[j] for j in corr_temp]
	else:
		corr_list = [[] for j in corr_temp]
		for i_2 in range(len(corr_temp)):
			if corr_temp[i_2][-2:] == "ae":
				corr_list[i_2].append(corr_temp[i_2][:-1])
				corr_list[i_2].append(corr_temp[i_2])
			elif corr_temp[i_2][-2:] == "is":
				corr_list[i_2].append(corr_temp[i_2])
				if corr_temp[i_2][-3:] == "cis":
					if corr_temp[i_2][-4:] == "icis":
						corr_list[i_2].append(corr_temp[i_2][:-4] + "ex")
					corr_list[i_2].append(corr_temp[i_2][:-3] + "x")
				elif corr_temp[i_2][-3:] == "dis":
					corr_list[i_2].append(corr_temp[i_2][:-3] + "s")
					corr_list[i_2].append(corr_temp[i_2][:-3] + "des")
				elif corr_temp[i_2][-4:] == "inis" or corr_temp[i_2][-4:] == "onis":
					corr_list[i_2].append(corr_temp[i_2][:-4] + "o")
				elif corr_temp[i_2][-4:] == "ntis":
					corr_list[i_2].append(corr_temp[i_2][:-4] + "s")
				elif corr_temp[i_2][-4:] == "itis":
					corr_list[i_2].append(corr_temp[i_2][:-3] + "s")
				elif corr_temp[i_2][-4:] == "this":
					corr_list[i_2].append(corr_temp[i_2][:-4] + "thes")
				else:					
					corr_list[i_2].append(corr_temp[i_2][:-2])
			elif corr_temp[i_2][-1] == "i":
				corr_list[i_2].append(corr_temp[i_2])
				corr_list[i_2].append(corr_temp[i_2][:-1]+"os")
				corr_list[i_2].append(corr_temp[i_2][:-1]+"on")
				corr_list[i_2].append(corr_temp[i_2][:-1]+"um")
				corr_list[i_2].append(corr_temp[i_2][:-1]+"us")
			else:
				corr_list[i_2].append(corr_temp[i_2])
		return corr_list

## Mapping subprocesses

def first_bare_match(string, tr_list, rn_dict, rg_dict, rr_dict):
	#
	tar_str = string.lower()
	tar_gen = tar_str.split()[0]
	tar_corr, tar_cstat = input_correction(tar_str)
	fb_match = [[] for j in tr_list]
	fb_stat = [1000 for j in tr_list]

	# Exact matching

	for fb_1 in range(len(tr_list)):
		if tar_str in rr_dict[tr_list[fb_1]]: # 
			fb_match[fb_1].extend(rr_dict[tr_list[fb_1]][tar_str])
			fb_stat[fb_1] = 0
		elif tar_cstat[1] == 0 and tar_corr in rr_dict[tr_list[fb_1]]:
			fb_match[fb_1].extend(rr_dict[tr_list[fb_1]][tar_corr])
			fb_stat[fb_1] = 3
		elif tar_gen in rr_dict[tr_list[fb_1]]:
			fb_match[fb_1].extend(rr_dict[tr_list[fb_1]][tar_gen])
			fb_stat[fb_1] = 100

	# Recursive mapping

	for fb_2 in range(len(tr_list)):
		if min(fb_stat) < 10 and fb_stat[fb_2] >= 10:
			syn_sub = []
			for fb_3 in range(len(tr_list)):
				if fb_stat[fb_3] == 0:
					for fb_4 in rr_dict[tr_list[fb_3]][tar_str]:
						syn_sub.append(rn_dict[tr_list[fb_3]][fb_4][1])
						ssl_2 = rn_dict[tr_list[fb_3]][fb_4][2].split("|")
						for fb_5 in ssl_2:
							if len(fb_5) >= 1 and fb_5 in rn_dict[tr_list[fb_3]][fb_4][1]:
								syn_sub.append(fb_5)
				elif fb_stat[fb_3] == 3:
					for fb_4 in rr_dict[tr_list[fb_3]][tar_corr]:
						syn_sub.append(rn_dict[tr_list[fb_3]][fb_4][1])
						ssl_2 = rn_dict[tr_list[fb_3]][fb_4][2].split("|")
						for fb_5 in ssl_2:
							if len(fb_5) >= 1 and fb_5 in rn_dict[tr_list[fb_3]][fb_4][1]:
								syn_sub.append(fb_5)
			syn_sub = [j.lower() for j in numpy.unique(syn_sub) if len(j) >= 1]
			syn_tar = []
			for fb_5 in syn_sub:
				if fb_5 in rr_dict[tr_list[fb_2]]:
					syn_tar.extend(rr_dict[tr_list[fb_2]][fb_5])
					fb_stat[fb_2] = 6
			if len(syn_tar) >= 1:
				fb_match[fb_2][:] = syn_tar
	
	# Clearance & return

	for fb_7 in range(len(tr_list)):
		if len(fb_match[fb_7]) >= 2:
			fb_match[fb_7][:] = list(numpy.unique(fb_match[fb_7]))
		#
		if len(fb_match[fb_7]) == 1:
			if fb_stat[fb_7] == 0:
				if rn_dict[tr_list[fb_7]][fb_match[fb_7][0]][1].lower() == tar_str:
					pass
				elif tar_str in rn_dict[tr_list[fb_7]][fb_match[fb_7][0]][2].lower():
					fb_stat[fb_7] += 1
			elif fb_stat[fb_7] == 3:
				if rn_dict[tr_list[fb_7]][fb_match[fb_7][0]][1].lower() == tar_corr:
					pass
				elif tar_corr in rn_dict[tr_list[fb_7]][fb_match[fb_7][0]][2].lower():
					fb_stat[fb_7] += 1
		elif len(fb_match[fb_7]) >= 2:
			if fb_stat[fb_7] < 10:
				fb_stat[fb_7] += 2
		#
	return fb_match, fb_stat
	
	#

def lowest_taxon_match(prev_map, prev_stat, tr_list, rn_dict, rr_dict):
	#
	broad_type = {
	"|R|": "Archaea", "|B6LM6|": "Bacteria", "|F|": "Fungi", "|N|": "Metazoa", "|P|": "Plant", 
	"|EOL-000000024748|": "Archaea", "|EOL-000000000003|": "Bacteria", "|EOL-000002172573|": "Fungi", "|EOL-000000541397|": "Metazoa", "|EOL-000000097815|": "Plant", 
	"|2|-": "Archaea", "|3|": "Bacteria", "|5|": "Fungi", "|1|": "Metazoa", "|6|": "Plant", 
	"|2157|": "Archaea", "|2|131567|": "Bacteria", "|4751|": "Fungi", "|33208|": "Metazoa", "|33090|": "Plant"
	}
	#
	lt_match = [[] for j in tr_list]
	lt_stat = [j for j in prev_stat]
	#
	for ltm_1 in range(len(tr_list)):
		if prev_stat[ltm_1] >= 10:
			match_type_r = "(others)"
			match_r_list = []
			#
			for ltm_2 in range(len(tr_list)):
				if prev_stat[ltm_2] < 5: # Accept exact match only
					#
					for ltm_3 in prev_map[ltm_2]:
						phylo_tree = rn_dict[tr_list[ltm_2]][ltm_3][4].split("|")
						phylo_rank = [int(j) for j in rn_dict[tr_list[ltm_2]][ltm_3][5].split("|")]
						for ltm_11 in broad_type:
							if ltm_11 in rn_dict[tr_list[ltm_2]][ltm_3][4] and match_type_r == "(others)":
								match_type_r = broad_type[ltm_11]
								for ltm_12 in phylo_tree:
									if ltm_12 in rn_dict[tr_list[ltm_2]]:
										match_r_list.append(rn_dict[tr_list[ltm_2]][ltm_12][1].lower())
								break
			#
			ltm_sub = {}
			ltm_depth = {}
			for ltm_2 in prev_map[ltm_1]:
				phylo_tree = rn_dict[tr_list[ltm_1]][ltm_2][4].split("|")
				phylo_rank = [int(j) for j in rn_dict[tr_list[ltm_1]][ltm_2][5].split("|")]
				match_type_c = "(others)"
				for ltm_12 in broad_type:
					if ltm_12 in rn_dict[tr_list[ltm_1]][ltm_2][4]:
						match_type_c = broad_type[ltm_12]
						break
				if match_type_c == match_type_r:
					ltm_sub[ltm_2] = max(phylo_rank)
					phylo_tree_names = []
					phylo_score = 0
					for ltm_13 in phylo_tree:
						if ltm_13 in rn_dict[tr_list[ltm_1]]:
							phylo_tree_names.append(rn_dict[tr_list[ltm_1]][ltm_13][1].lower())
					for ltm_14 in phylo_tree_names:
						if ltm_14 in match_r_list:
							phylo_score += 1
					ltm_depth[ltm_2] = phylo_score
			#
			if len(ltm_sub) == 0:
				for ltm_10 in prev_map[ltm_1]:
					ltm_sub[ltm_10] = max([int(j) for j in rn_dict[tr_list[ltm_1]][ltm_10][5].split("|")])
					ltm_depth[ltm_10] = 0
			elif max(ltm_sub.values()) < 7 and len(prev_map[ltm_1]) >= 1:
				for ltm_10 in prev_map[ltm_1]:
					ltm_sub[ltm_10] = max([int(j) for j in rn_dict[tr_list[ltm_1]][ltm_10][5].split("|")])
					ltm_depth[ltm_10] = 0
			#
			if len(ltm_sub) >= 1:
				ltm_sub = dict(sorted(ltm_sub.items(), key=lambda item: item[1], reverse = True))
				lt_match[ltm_1][:] = []
				ltm_max_depth = max(ltm_sub.values())
				ltm_max_score = max(ltm_depth.values())
				for ltm_6 in ltm_sub:
					if ltm_sub[ltm_6] == ltm_max_depth and ltm_depth[ltm_6] == ltm_max_score:
						lt_match[ltm_1].append(ltm_6)
						lt_stat[ltm_1] = 18-ltm_sub[ltm_6]
					else:
						break
			#
	for ltm_7 in range(len(lt_match)):
		if len(lt_match[ltm_7]) == 0:
			lt_match[ltm_7].extend(prev_map[ltm_7])
	#
	return lt_match, lt_stat

def tax_rule_screening(string):
	#
	string_list = string.split()

	# Non-organism flag
	remove_key = ['collection', 'construct', 'library', 'plasmid', 'sequence', 'transposon', 'vector', "<<"]
	for trs_1 in remove_key:
		if trs_1 in string_list:
			return 90

	# Unclassified/Uncultured/Unidentified flag
	unclear_key = ['unclassified', 'uncultured', 'unidentified']
	for trs_2 in unclear_key:
		if trs_2 in string_list:
			return 91

	# Environmental sample/Enrichment culture flag
	env_key = ['environmental', 'sample', 'enrichment', 'culture']
	for trs_3 in env_key:
		if trs_3 in string_list:
			return 92

	# Virus flag
	virus_key = ['virus', 'phage']
	for trs_4 in virus_key:
		if trs_4 in string_list:
			return 93

	# Phytoplasma flag
	if 'phytoplasma' in string_list:
		return 94

	# Endosymbiont flag
	symbiont_key = ["symbiont", "endosymbiont", "symbiotic"]
	for trs_5 in symbiont_key:
		if trs_5 in string_list:
			return 95

	# Hybrid flag
	hybrid_key = ["x", "×"]
	for trs_6 in hybrid_key:
		if trs_6 in string_list:
			return 96

	# Materia medica // complicated input flag
	multiple_key = ['lac', 'massa fermentata', 'sal', 'seu']
	for trs_7 in multiple_key:
		if trs_7 in string_list:
			return 97

	# (all other else)
	return 0

def intra_generic_edit_dist(string, cut_dist, tr_list, rn_dict, rg_dict, rr_dict):
	
	# Initialization

	string_block = string.lower().split()
	first_block = string_block[0]
	second_block = ""
	if len(string_block) >= 2:
		second_block += string_block[1]
	igd_match = [[] for j in tr_list]
	igd_stat = [1000 for j in tr_list]

	# Candidate 1st mapping

	lowest_case = []

	for igd_1 in numpy.arange(len(string_block), 0, -1):
		sst = ' '.join(string_block[0:igd_1])
		for igd_2 in range(len(tr_list)):
			if sst in rr_dict[tr_list[igd_2]]:
				lowest_case.extend(sst.split())
				break
		if len(lowest_case) >= 1:
			break

	# 

	if len(lowest_case) == 0:
		return igd_match, igd_stat

	if any(char.isdigit() for char in string) == True:
		for igd_3 in range(len(tr_list)):
			for igd_4 in numpy.arange(len(lowest_case), 0, -1):
				if ' '.join(lowest_case[0:igd_4]) in rr_dict[tr_list[igd_3]]:
					igd_match[igd_3].extend(rr_dict[tr_list[igd_3]][' '.join(lowest_case[0:igd_4])])
					igd_stat[igd_3] = 25
					break
		return igd_match, igd_stat

	# Candidate 2nd mapping
	
	ig_temp_names = []
	for igd_5 in range(len(tr_list)):
		i_sub = []
		if lowest_case[0] in rg_dict[tr_list[igd_5]]:
			if len(lowest_case) >= 2:
				sst = ' '.join(lowest_case)
				for igd_6 in rg_dict[tr_list[igd_5]][lowest_case[0]][1]:
					if sst in rn_dict[tr_list[igd_5]][igd_6][1].lower() or sst in rn_dict[tr_list[igd_5]][igd_6][2].lower():
						i_sub.append(igd_6)
			elif len(lowest_case) == 1:
				i_sub.extend(rg_dict[tr_list[igd_5]][lowest_case[0]][1])
		#
		if len(i_sub) >= 1:
			for igd_7 in i_sub:
				ssl = rn_dict[tr_list[igd_5]][igd_7][1].lower().split()
				if ssl[0] == lowest_case[0] and len(ssl) >= len(lowest_case):
					if ' '.join(lowest_case) in rn_dict[tr_list[igd_5]][igd_7][1].lower():
						if len(second_block) >= 1 and len(ssl) >= 2:
							if ssl[1][0] == second_block[0] or ssl[1][1:min(5, len(ssl[1]))] in second_block:
								ig_temp_names.append(rn_dict[tr_list[igd_5]][igd_7][1].lower())
						else:
							ig_temp_names.append(rn_dict[tr_list[igd_5]][igd_7][1].lower())
				if len(rn_dict[tr_list[igd_5]][igd_7][2]) >= 1:
					ssl_2 = rn_dict[tr_list[igd_5]][igd_7][2].lower().split("|")
					for igd_8 in ssl_2:
						ssl_3 = igd_8.split()
						if ssl_3[0] == lowest_case[0] and len(ssl_3) == len(lowest_case):
							if ' '.join(lowest_case) in igd_8:
								if len(second_block) >= 1 and len(ssl_3) >= 2:
									if ssl_3[1][0] == second_block[0] or ssl_3[1][1:min(5, len(ssl_3[1]))] in second_block:
										ig_temp_names.append(igd_8)
								else:
									ig_temp_names.append(igd_8)

	ig_temp_names[:] = list(numpy.unique(ig_temp_names))

	#

	sorted_list = {}
	for igd_5 in ig_temp_names:
		if numpy.fabs(len(igd_5) - len(string)) <= cut_dist:
			ed_value = lev_dist(string.lower(), igd_5, cut_dist)
			if ed_value <= cut_dist and ed_value < float(len(string)-len(first_block))*(0.3333):
				if string.lower() in igd_5 or igd_5 in string.lower():
					sorted_list[igd_5] = min(ed_value-1, 0)
				else:
					sorted_list[igd_5] = ed_value
	sorted_list = dict(sorted(sorted_list.items(), key=lambda item: item[1]))

	#

	if len(sorted_list) >= 1:
		min_ed = min(sorted_list.values())
		igd_match_temp = [[], []]
		for igd_6 in sorted_list:
			if sorted_list[igd_6] == min_ed:
				pc_map_igd, pc_stat_igd = first_bare_match(igd_6, tr_list, rn_dict, rg_dict, rr_dict)
				igd_match_temp[0].append(pc_map_igd)
				igd_match_temp[1].append(pc_stat_igd)
		for igd_7 in range(len(tr_list)):
			min_map = min([j[igd_7] for j in igd_match_temp[1]])
			for igd_8 in range(len(igd_match_temp[1])):
				if igd_match_temp[1][igd_8][igd_7] == min_map:
					igd_match[igd_7].extend(igd_match_temp[0][igd_8][igd_7])
					igd_stat[igd_7] = min_map

	#

	if min(igd_stat) < 10:
		igd_lt_map, igd_lt_stat = lowest_taxon_match(igd_match, igd_stat, tr_list, rn_dict, rr_dict)

		#

		for igd_9 in range(len(igd_stat)):
			if igd_stat[igd_9] in [0, 1, 3, 4]:
				igd_stat[igd_9] = 20
			elif igd_stat[igd_9] in [2, 5]:
				igd_stat[igd_9] = 21
			elif igd_stat[igd_9] == 6:
				igd_stat[igd_9] = 22
			elif igd_stat[igd_9] == 8:
				igd_stat[igd_9] = 23
			elif igd_lt_stat[igd_9] < 20:
				igd_stat[igd_9] = 24
				igd_match[igd_9][:] = igd_lt_map[igd_9][:]
			else:
				igd_stat[igd_9] = 1000
				igd_match[igd_9][:] = []
	else:
		for igd_9 in range(len(igd_stat)):
			igd_stat[igd_9] = 1000
			igd_match[igd_9][:] = []
	#
	#print(string, len(ig_temp_names), len(sorted_list), sorted_list, igd_match, igd_stat)
	return igd_match, igd_stat

def in_depth_edit_dist(string, cut_dist, tr_list, rn_dict, rg_dict, rr_dict):
	#
	idd_match = [[] for j in tr_list]
	idd_stat = [1000 for j in tr_list]
	# Latin correction raw match
	latin_corr = [' '.join(j) for j in list(itertools.product(*latin_correction(string.split())))]
	latin_first_block = ""
	if len(string) >= 1:
		if len(latin_corr) >= 1:
			if len(latin_corr[0]) >= 1:
				latin_first_block += latin_corr[0].split()[0]
			else:
				latin_first_block += string.split()[0]
		else:
			latin_first_block += string.split()[0]
	#
	for idd_2 in range(len(tr_list)):
		if idd_stat[idd_2] > 100 and len(latin_corr) >= 1:
			for idd_3 in latin_corr:
				if idd_3 in rr_dict[tr_list[idd_2]]:
					idd_match[idd_2].extend(rr_dict[tr_list[idd_2]][idd_3])
					idd_stat[idd_2] = 0
					break
		#
		if idd_stat[idd_2] > 100 and len(latin_corr) >= 1:		
			for idd_3 in latin_corr:
				ssl = idd_3.split()
				ssl_len = len(ssl)
				if len(ssl) >= 3:
					for idd_4 in numpy.arange(len(ssl)-1, 1, -1):
						if ' '.join(ssl[0:idd_4]) in rr_dict[tr_list[idd_2]]:
							idd_match[idd_2].extend(rr_dict[tr_list[idd_2]][' '.join(ssl[0:idd_4])])
							if idd_4 >= 3:
								idd_stat[idd_2] = 0
							else:
								idd_stat[idd_2] = 100
							break
		#
		if idd_stat[idd_2] > 100 and len(latin_corr) >= 1:		
			for idd_3 in latin_corr:
				ssl = idd_3.split()
				if len(ssl) >= 1:
					if ssl[0] in rr_dict[tr_list[idd_2]]:
						idd_match[idd_2].extend(rr_dict[tr_list[idd_2]][ssl[0]])
						idd_stat[idd_2] = 100
						latin_first_block = ssl[0]
						break
	#
	if any(char.isdigit() for char in string) == True:
		pass
	else:
		if idd_stat[idd_2] > 10 and len(latin_corr) >= 1:
			if len(latin_corr[0].split()) > 1:
				rg_min = min([len(rg_dict[j]) for j in rg_dict])
				first_alpha = alpha_count(latin_first_block)
				first_len = len(latin_first_block)
				first_raw = {}
				for idd_4 in range(len(tr_list)):
					if len(rg_dict[tr_list[idd_4]]) == rg_min:
						for idd_5 in rg_dict[tr_list[idd_4]]:
							if numpy.fabs(len(idd_5)-first_len) < cut_dist:
								if alpha_compare(rg_dict[tr_list[idd_4]][idd_5][0], first_alpha, cut_dist) <=  cut_dist:
									ed_score = lev_dist(idd_5, latin_first_block, cut_dist)
									if ed_score <= min(cut_dist, int(float(first_len)/3.0)):
										first_raw[idd_5] = ed_score
				first_raw = dict(sorted(first_raw.items(), key=lambda item: item[1]))
				#
				second_raw = []
				for idd_6 in range(len(tr_list)):
					for idd_7 in first_raw:
						if idd_7 in rg_dict[tr_list[idd_6]]:
							for idd_8 in rg_dict[tr_list[idd_6]][idd_7][1]:
								sst = rn_dict[tr_list[idd_6]][idd_8][1].lower()
								if idd_7 in sst:
									second_raw.append(sst)
								if len(rn_dict[tr_list[idd_6]][idd_8][2]) >= 1:
									ssl = rn_dict[tr_list[idd_6]][idd_8][2].lower().split("|")
									for idd_9 in ssl:
										if idd_7 in idd_9:
											second_raw.append(idd_9)
				second_raw[:] = list(numpy.unique(second_raw))
				#
				sorted_list = {}
				for idd_10 in second_raw:
					if numpy.fabs(len(string) - len(idd_10)) <= cut_dist:
						ed_score = lev_dist(string, idd_10, cut_dist)
						if ed_score <= cut_dist:
							if idd_10 in sorted_list:
								if ed_score < sorted_list[idd_10]:
									sorted_list[idd_10] = ed_score
							else:
								sorted_list[idd_10] = ed_score
					#
					if len(latin_corr) >= 1:
						if numpy.fabs(len(latin_corr[0]) - len(idd_10)) <= cut_dist:
							ed_score_2 = lev_dist(latin_corr[0], idd_10, cut_dist)
							if ed_score_2 <= cut_dist:
								if idd_10 in sorted_list:
									if ed_score_2 < sorted_list[idd_10]:
										sorted_list[idd_10] = ed_score_2
								else:
									sorted_list[idd_10] = ed_score_2
				sorted_list = dict(sorted(sorted_list.items(), key=lambda item: item[1]))
				if len(sorted_list) >= 1:
					sorted_min = min(sorted_list.values())
					for idd_11 in sorted_list:
						if sorted_list[idd_11] == sorted_min:
							for idd_12 in range(len(tr_list)):
								if idd_11 in rr_dict[tr_list[idd_12]]:
									idd_match[idd_12][:] = rr_dict[tr_list[idd_12]][idd_11]
									idd_stat[idd_12] = 3
						else:
							break
		#
		for idd_13 in range(len(tr_list)):
			if min(idd_stat) < 10 and idd_stat[idd_13] >= 10:
				syn_sub = []
				for idd_14 in range(len(tr_list)):
					if idd_stat[idd_14] < 10:
						for idd_15 in idd_match[idd_14]:
							syn_sub.append(rn_dict[tr_list[idd_14]][idd_15][1])
							ssl_2 = rn_dict[tr_list[idd_14]][idd_15][2].split("|")
							for idd_16 in ssl_2:
								if len(idd_16) >= 1 and idd_16 in rn_dict[tr_list[idd_14]][idd_15][1]:
									syn_sub.append(idd_16)
				syn_sub = [j.lower() for j in numpy.unique(syn_sub) if len(j) >= 1]
				syn_tar = []
				for idd_16 in syn_sub:
					if idd_16 in rr_dict[tr_list[idd_13]]:
						syn_tar.extend(rr_dict[tr_list[idd_13]][idd_16])
						idd_stat[idd_13] = 6
				if len(syn_tar) >= 1:
					idd_match[idd_13][:] = list(numpy.unique(syn_tar))
	#
	if min(idd_stat) < 10:
		idd_lt_map, idd_lt_stat = lowest_taxon_match(idd_match, idd_stat, tr_list, rn_dict, rr_dict)
		#
		for idd_17 in range(len(idd_stat)):
			if idd_stat[idd_17] == 0:
				idd_stat[idd_17] = 30
			elif idd_stat[idd_17] == 1:
				idd_stat[idd_17] = 31
			elif idd_stat[idd_17] == 3:
				if len(idd_match[idd_17]) == 1:
					idd_stat[idd_17] = 32
				else:
					idd_stat[idd_17] = 33
			elif idd_stat[idd_17] == 6:
				if len(idd_match[idd_17]) == 1:
					idd_stat[idd_17] = 34
				else:
					idd_stat[idd_17] = 35
			elif idd_stat[idd_17] < 20:
				idd_stat[idd_17] = 36
			elif idd_stat[idd_17] < 1000:
				idd_stat[idd_17] = 100
			else:
				idd_stat[idd_17] = 1000
				idd_match[idd_17][:] = []
	else:
		for idd_17 in range(len(idd_stat)):
			idd_stat[idd_17] = 1000
			idd_match[idd_17][:] = []
	#
	#print(string, idd_match, idd_stat)
	#
	return idd_match, idd_stat

## Partial mapping sequence

def partial_mapping(string, tr_list, rn_dict, rg_dict, rr_dict):
	#
	tar_str = string.lower()
	pm_match = [[] for j in tr_list]
	pm_stat = [1000 for j in tr_list]
	#
	tar_split_list = tar_str.split()
	for pm_1 in range(len(tar_split_list)):
		tar_partial = ' '.join(tar_split_list[:(len(tar_split_list)-pm_1)])
		tar_part_corr, tar_part_cstat = input_correction(tar_str)
		for pm_2 in range(len(tr_list)):
			if pm_stat[pm_2] > 100:
				tar_det = 0
				if tar_partial in rr_dict[tr_list[pm_2]]:
					pm_match[pm_2].extend(rr_dict[tr_list[pm_2]][tar_partial])
					pm_stat[pm_2] = 100
					tar_det += 1
				elif tar_part_cstat[1] == 0 and tar_part_corr in rr_dict[tr_list[pm_2]]:
					pm_match[pm_2].extend(rr_dict[tr_list[pm_2]][tar_part_corr])
					pm_stat[pm_2] = 100
					tar_det += 1
				#
				if tar_det >= 1:
					for pm_3 in range(len(tr_list)):
						if pm_stat[pm_3] > 100:
							syn_sub = []
							for pm_4 in pm_match[pm_2]:
								syn_sub.append(rn_dict[tr_list[pm_2]][pm_4][1])
								if len(rn_dict[tr_list[pm_2]][pm_4][2]) >= 1:
									ssl_2 = rn_dict[tr_list[pm_2]][pm_4][2].split("|")
									for pm_5 in ssl_2:
										if len(pm_5) >= 1 and pm_5 in rn_dict[tr_list[pm_2]][pm_4][1]:
											syn_sub.append(pm_5)
							syn_sub = [j.lower() for j in numpy.unique(syn_sub) if len(j) >= 1]
							syn_tar = []
							for pm_6 in syn_sub:
								if pm_6 in rr_dict[tr_list[pm_3]]:
									syn_tar.extend(rr_dict[tr_list[pm_3]][pm_6])
							if len(syn_tar) >= 1:
								pm_match[pm_3][:] = syn_tar
								pm_stat[pm_3] = 100
	return pm_match, pm_stat
	#


## Core mapping sequence

def phylosophos_sequential_mapping(string, ref_sel, cut_dist, tr_list, rn_dict, rg_dict, rr_dict):

	# Initialization

	ref_ord = 0
	for pcm_1 in range(len(tr_list)):
		if ref_sel == tr_list[pcm_1]:
			ref_ord += pcm_1
			break
	pc_map_map = [[] for j in tr_list]
	pc_map_stat = [1000 for j in tr_list]

	# Step 1. Exact matching

	pc_map_1, pc_stat_1 = first_bare_match(string, tr_list, rn_dict, rg_dict, rr_dict)
	for pc_1 in range(len(tr_list)):
		if pc_stat_1[pc_1] < pc_map_stat[pc_1]:
			pc_map_map[pc_1][:] = pc_map_1[pc_1][:]
			pc_map_stat[pc_1] = pc_stat_1[pc_1]

	# Step 2. Closest taxon mapping

	if min(pc_map_stat) < 10:
		pc_map_2, pc_stat_2 = lowest_taxon_match(pc_map_1, pc_stat_1, tr_list, rn_dict, rr_dict)
		for pc_1 in range(len(tr_list)):
			if pc_stat_2[pc_1] < pc_map_stat[pc_1]:
				pc_map_map[pc_1][:] = pc_map_2[pc_1][:]
				pc_map_stat[pc_1] = pc_stat_2[pc_1]
		if pc_map_stat[ref_ord] < 20:
			return pc_map_map, pc_map_stat

	# Step 3. Rule-based screening

	pc_rule_stat = tax_rule_screening(string)

	if pc_rule_stat >= 90:
		pc_map_stat[:] = [pc_rule_stat for j in tr_list]
		return [[], [], [], []], pc_map_stat

	if cut_dist == 0: # Ignore edit-distance based mapping process if cutoff distance is zero
		return pc_map_map, pc_map_stat

	# Step 4. (intrageneric) edit distance-based mapping

	corr_string = input_correction(string)[0]

	if sum([len(j) for j in pc_map_1]) >= 1:
		pc_map_3, pc_stat_3 = intra_generic_edit_dist(corr_string, cut_dist, tr_list, rn_dict, rg_dict, rr_dict)
		for pc_1 in range(len(tr_list)):
			if pc_stat_3[pc_1] < pc_map_stat[pc_1]:
				pc_map_map[pc_1][:] = pc_map_3[pc_1][:]
				pc_map_stat[pc_1] = pc_stat_3[pc_1]
		if min(pc_map_stat) < 30:
			return pc_map_map, pc_map_stat

	# Step 5. (intergeneric) edit distance-based mapping

	pc_map_4, pc_stat_4 = in_depth_edit_dist(corr_string, cut_dist, tr_list, rn_dict, rg_dict, rr_dict)
	for pc_1 in range(len(tr_list)):
		if pc_stat_4[pc_1] < pc_map_stat[pc_1]:
			pc_map_map[pc_1][:] = pc_map_4[pc_1][:]
			pc_map_stat[pc_1] = pc_stat_4[pc_1]

	# Step 6. Mapping status return

	pc_map_5, pc_stat_5 = partial_mapping(string, tr_list, rn_dict, rg_dict, rr_dict)
	for pc_1 in range(len(tr_list)):
		if pc_stat_5[pc_1] <= pc_map_stat[pc_1]:
			pc_map_map[pc_1][:] = pc_map_5[pc_1][:]
			pc_map_stat[pc_1] = pc_stat_5[pc_1]

	return pc_map_map, pc_map_stat

## Export function

def phylosophos_result_export(input_name, raw_list, precalc_list, mapping_results, ref_type, tr_list, rn_dict):

	# Preset information

	mapping_status_dict = {
	"0":"Raw / Exact DB / Canonical match", 
	"1":"Raw / Exact DB / Synonym match", 
	"2":"Raw / Exact DB / Multiple match", 
	"3":"Simple corrected / Exact DB / Canonical match", 
	"4":"Simple corrected / Exact DB / Synonym match", 
	"5":"Simple corrected / Exact DB / Multiple match", 
	"6":"Recursive / Single match", 
	"8":"Recursive / Multiple match", 
	"10":"Recursive / Nearest match / Species level", 
	"11":"Recursive / Nearest match / Genus level", 
	"12":"Recursive / Nearest match / Family level", 
	"13":"Recursive / Nearest match / Class level", 
	"14":"Recursive / Nearest match / Order level", 
	"15":"Recursive / Nearest match / Phylum level", 
	"16":"Recursive / Nearest match / Kingdom level", 
	"17":"Recursive / Nearest match / Domain level", 
	"20":"Specific epithet corrected / Exact DB / Single match", 
	"21":"Specific epithet corrected / Exact DB / Multiple match", 
	"22":"Specific epithet corrected / Recursive / Single match", 
	"23":"Specific epithet corrected / Recursive / Multiple match", 
	"24":"Specific epithet corrected / Recursive / Nearest match", 
	"30":"Latin inflection corrected / Exact DB / Single match", 
	"31":"Latin inflection corrected / Exact DB / Nearest match", 
	"32":"Generic epithet corrected / Exact DB / Single match",
	"33":"Generic epithet corrected / Exact DB / Multiple match",
	"34":"Generic epithet corrected / Recursive / Single match",
	"35":"Generic epithet corrected / Recursive / Multiple match",
	"36":"Generic epithet corrected / Recursive / Nearest match",
	"40":"Correction denied / Strain name involved / Nearest match",
	"41":"Correction denied / Similarity-related abbreviation identified / Nearest match",
	"90":"Rule-based screening / Non-organism", 
	"91":"Rule-based screening / Unclassified-Uncultured-Unidentified", 
	"92":"Rule-based screening / Environmental sample", 
	"93":"Rule-based screening / Virus or phage - manual check required", 
	"94":"Rule-based screening / Phytoplasma - manual check required", 
	"95":"Rule-based screening / (endo)symbiont - manual check required", 
	"96":"Rule-based screening / Unresolvable hybrid - manual check required", 
	"97":"Rule-based screening / Multiple materia medica - manual check required", 
	"100":"Unmapped / partial match (most likely genus level)",  
	"1000":"Unmapped"
	}

	export_header = ["Input_file_name", "Input_original_order", "Raw_name_input", "Pre_corrected_input", 
	"Chosen_reference", "Chosen_reference_mapped_ID", "Chosen_reference_scientific_name", 
	"Chosen_reference_mapping_status_code", "Chosen_reference_mapping_status_description"]

	# Initialization

	export_time = str(datetime.datetime.now())[2:19].replace("-", "").replace(":", "").replace(" ", "_")
	export_file_name = "phylosophos_result_"+export_time+"_"+input_name.split("\\")[-1]

	for pre_1 in tr_list:
		export_header.append(pre_1+"_mapped_ID")
		export_header.append(pre_1+"_scientific_name")
		export_header.append(pre_1+"_mapping_status_code")

	export_header.append("Manual_curation_recommended")

	ref_ord = 0
	for pcm_1 in range(len(tr_list)):
		if ref_type == tr_list[pcm_1]:
			ref_ord += pcm_1
			break

	# Similar species mapping correction

	for i_1 in range(len(mapping_results)):
		id_corr, id_stat = input_correction(str(raw_list[i_1]))
		if id_stat[1] == 1:
			for i_2 in range(4):
				if mapping_results[i_1][1][i_2] >= 20 and mapping_results[i_1][1][i_2] < 100 and len(mapping_results[i_1][0][i_2]) == 1:
					higher_id = []
					phylo_code = rn_dict[tr_list[i_2]][mapping_results[i_1][0][i_2][0]][4].split("|")
					phylo_level = [int(j) for j in rn_dict[tr_list[i_2]][mapping_results[i_1][0][i_2][0]][5].split("|")]
					if len(str(raw_list[i_1]).split()) <= 2:
						for i_3 in numpy.arange(1, len(phylo_code), 1):
							if phylo_level[i_3] >= 1:
								higher_id.append(phylo_code[i_3])
								break
					else:
						if max(phylo_level) <= 7:
							higher_id.append(phylo_code[0])
							break
						else:
							for i_3 in numpy.arange(1, len(phylo_code), 1):
								if phylo_level[i_3] >= 1:
									higher_id.append(phylo_code[i_3])
									break
					mapping_results[i_1][0][i_2][:] = higher_id[:]
					mapping_results[i_1][1][i_2] = 41
		for i_4 in range(4):
			if mapping_results[i_1][1][i_4] == 25:
				mapping_results[i_1][1][i_4] = 40

	# Mapping data export

	b_path = os.getcwd()+"\\"
	with open(b_path+"result\\"+export_file_name, 'w', encoding = 'UTF-8', errors = 'ignore') as res_f:
		res_f.write('\t'.join(export_header) + '\n')
		for i_1 in range(len(mapping_results)):
			res_f.write(input_name.split("\\")[-1] + '\t' + str(i_1+1) + '\t')
			res_f.write(str(raw_list[i_1]) + '\t' + str(precalc_list[i_1]) + '\t' + ref_type + '\t')
			res_f.write('|'.join(mapping_results[i_1][0][ref_ord]) + '\t')
			res_f.write('|'.join([rn_dict[ref_type][j][1] for j in mapping_results[i_1][0][ref_ord]]) + '\t')
			res_f.write(str(mapping_results[i_1][1][ref_ord]) + '\t')
			res_f.write(mapping_status_dict[str(mapping_results[i_1][1][ref_ord])])
			for i_2 in range(len(tr_list)):
				res_f.write('\t' + '|'.join(mapping_results[i_1][0][i_2]))
				res_f.write('\t' + '|'.join([rn_dict[tr_list[i_2]][j][1] for j in mapping_results[i_1][0][i_2]]))
				res_f.write('\t' + str(mapping_results[i_1][1][i_2]))
			#
			if str(mapping_results[i_1][1][ref_ord]) in ["0", "1", "3", "4", "6", "10", "20", "22", "30", "31", "32", "34"]:
				res_f.write('\t' + "NO")
			elif str(mapping_results[i_1][1][ref_ord]) in ["11", "12", "13", "14", "15", "16", "17", "24", "36"]:
				res_f.write('\t' + "MAYBE")
			else:
				res_f.write('\t' + "YES")
			res_f.write('\n')

## Analysis function

def phylosophos_core_analysis(input_list, tr_list, rn_dict, rg_dict, rr_dict, ref_type, lev_cutoff, manual_stat):

	# Initialization

	ref_ord = 0
	for pcm_1 in range(len(tr_list)):
		if ref_type == tr_list[pcm_1]:
			ref_ord += pcm_1
			break
	
	# Manual curation file import

	manual_dict = {}

	if manual_stat == True:
		b_path = os.getcwd()+"\\"
		with open(b_path+"pp_learning\\manual_curation_list.tsv", encoding = "UTF-8") as inp_f:
			inp_f.readline()
			for line in inp_f:
				ssl = line.rstrip('\n').split('\t')
				manual_dict[ssl[0].lower()] = ssl[1]

	# Input names simple correction

	precalc_list = []

	for pca_1 in range(len(input_list[1])):		
		if input_list[1][pca_1].lower() in manual_dict:
			precalc_list.append(manual_dict[input_list[1][pca_1].lower()])
		else:
			precalc_list.append(string_split_correction(input_list[1][pca_1]))

	# Core analysis

	mapping_results = []

	for pca_2 in precalc_list:
		map_label, map_stat = phylosophos_sequential_mapping(pca_2, ref_type, lev_cutoff, tr_list, rn_dict, rg_dict, rr_dict)
		mapping_results.append([map_label, map_stat])
		print(input_list[0], len(mapping_results), "/", len(precalc_list), end='\r')

	# Result export

	print("-", input_list[0], len(precalc_list), "analysis completed")
	phylosophos_result_export(input_list[0], input_list[1], precalc_list, mapping_results, ref_type, tr_list, rn_dict)

#

