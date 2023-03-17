#!/usr/bin/env python

################################################################
#### 
#### PROJECT PHYLOSOPHOS: FORMAL VERSION
#### 
#### INITIALIZATION SCRIPT: REFERENCE UPDATE MODULE
#### 
#### ORIGINAL SCRIPT WRITTEN BY MIN HYUNG CHO, PH.D.
#### 
#### BIOINFORMATICS AND MOLECULAR DESIGN RESEARCH CENTER
#### 
################################################################

#### Library import

## Base library

import datetime
import numpy
import os
import ssl
import sys
import tarfile
import zipfile

from urllib import request

## Custom library

from ps_init import ps_update

#### Function definition

#### Global path parameters

base_path = os.path.realpath(__file__)[:-32]
ef_path = "external_files\\"
pp_path = "pp_ref\\"

#### Main part

def main():
	# Parameter setup
	update_stat = 0
	if len(sys.argv) >= 2:
		if int(sys.argv[1]) != 0:
			update_stat += 1

	# Path setup
	if "external_files" not in os.listdir(base_path):
		os.mkdir(base_path+ef_path)
	if "pp_ref" not in os.listdir(base_path):
		os.mkdir(base_path+pp_path)
		os.mkdir(base_path+pp_path+"previous\\")
	if "input" not in os.listdir(base_path):
		os.mkdir(base_path+"input\\")
	if "result" not in os.listdir(base_path):
		os.mkdir(base_path+"result\\")

	# Archiving previous dictionary
	pdir_flist = os.listdir(base_path + pp_path)
	if len(pdir_flist) > 1:
		prev_path = ""
		if str(datetime.datetime.now())[2:10].replace("-", "") not in os.listdir(base_path + pp_path + "previous\\"):
			os.mkdir(base_path + pp_path + "previous\\" + str(datetime.datetime.now())[2:10].replace("-", ""))
		prev_path = base_path + pp_path + "previous\\" + str(datetime.datetime.now())[2:10].replace("-", "") + "\\"
		for i_1 in range(len(pdir_flist)):
			if "_dict.txt" in pdir_flist[i_1]:
				os.rename(base_path + pp_path + pdir_flist[i_1], prev_path + pdir_flist[i_1])
		print("#### Previous reference files: archived ####")
	else:
		print("#### No previous reference files found: initialize ####")

	# Update protocol
	ps_update.ncbi_tax_update_new(update_stat, base_path+ef_path, base_path+pp_path)
	ps_update.col_tax_update_new(update_stat, base_path+ef_path, base_path+pp_path, "latest_dwca.zip") # Change it as appropriate
	ps_update.eol_tax_update_new(update_stat, base_path+ef_path, base_path+pp_path, "dhv21.zip") # Change it as appropriate

if __name__ == "__main__":
	main()

#### Script ended
