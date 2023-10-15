#!/usr/bin/env python

################################################################
#### 
#### PROJECT PHYLOSOPHOS: MODULAR VERSION 231013
#### 
#### PREDICTOR SCRIPT
#### 
#### ORIGINAL SCRIPT WRITTEN BY MIN HYUNG CHO, PH.D.
#### 
#### BIOINFORMATICS AND MOLECULAR DESIGN RESEARCH CENTER
#### 
################################################################

#### LIBRARY IMPORT

## GENERIC LIBRARIES

import datetime
import gc
import itertools
import numpy
import os
import sys

## CUSTOM LIBRARIES

from phylosophos import ps_analysis, ps_initialize

#### GLOBAL PARAMETERS

alpha_numeric = {"a":0, "b":1, "c":2, "d":3, "e":4, "f":5, "g":6, "h":7, "i":8, "j":9, 
"k":10, "l":11, "m":12, "n":13, "o":14, "p":15, "q":16, "r":17, "s":18, "t":19, 
"u":20, "v":21, "w":22, "x":23, "y":24, "z":25}

#### MAIN FUNCTION

def phylosophos():
	base_path = os.getcwd()+"\\"
	ref_path = base_path + "pp_ref\\"

	# Initialization

	ref_type, input_type, input_file_name, default_cutoff, mc_stat = ps_initialize.phylosophos_initialization(sys.argv)

	# Input files import

	input_list = ps_initialize.phylosophos_input_import(input_type, input_file_name)

	# Reference files import

	tax_ref_list, ref_names_dict, ref_genus_dict, ref_raw_names = ps_initialize.phylosophos_ref_import(ref_path)

	# Core analysis & export

	for i_1 in input_list:
		ps_analysis.phylosophos_core_analysis(i_1, tax_ref_list, ref_names_dict, ref_genus_dict, ref_raw_names, ref_type, default_cutoff, mc_stat)

	#

	print("#### PhyloSophos analysis completed ####")

#### END OF SCRIPT
