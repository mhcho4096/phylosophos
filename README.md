# PhyloSophos

A python-based high-throughput scientific name mapping algorithm

## Table of contents

- [In a nutshell](#in-a-nutshell)
- [Summary](#summary)
- [Installation guide](#installation-guide)
- [Update guide](#update-guide)
- [Usage guide](#usage-guide)
- [Result format guide](#result-format-guide)
- [Citing PhyloSophos](#citing-phylosophos)
- [License](#license)
- [Contact](#contact)

## In a nutshell

```
python phylosophos_initialize_update.py [optional_update_parameter]

python phylosophos_core.py [[optional_parameter_type] [optional_parameter_value]]
```

## Summary

PhyloSophos is a high-throughput scientific name processor which achieves greater mapping performance by referencing multiple taxonomic references and recognizing the semantic structure of scientific names. It also corrects common Latin variants and vernacular names, which often appear in various biological databases and resources. 

Please refer to "phylosophos_guide.pdf" for detailed information on the package and instructions on how to use and modify it.

## Installation guide

PhyloSophos requires Numpy, along with other basic libraries of Python 3. Before initializing PhyloSophos to fit with your environment, please set up an environment with a python version >=3.7 (preferably with Conda), and install numpy (https://numpy.org/install/) as appropriate.

A script for initialization and update is phylosophos_initialize_update.py. This script will download raw taxonomic metadata from Catalogue of Life (CoL), Encyclopedia of Life (EoL), and NCBI Taxonomy, and then process it into reference files. You may execute this script as:

<pre><code>python phylosophos_initialize_update.py [optional_update_parameter]</code></pre>

For initialization purposes, the optional parameter will not affect downstream processes. It is possible for Catalogue of Life and Encyclopedia of Life to change the name of the metadata file in their FTP server. If this happens, please change file names in the phylosophos_initialize_update.py (see lines 76 & 77) as appropriate.

In addition, there was an instance when the taxonomic databases changed the format of their metadata, making the initialization script non-functional. If this happens, please notify me immediately (see contact information).

## Update guide

Again, taxonomic reference update uses phylosophos_initialize_update.py script. You may execute this script as:

<pre><code>python phylosophos_initialize_update.py 1</code></pre>

In this case, optional update parameter affects the downstream process: by default, if a database metadata file is found in the /external files directory, the metadata download step is skipped to reduce processing time (allowing manual download of reference metadata). If an integer value other than 0 is given as an optional argument, the update script will start downloading the reference metadata file, overriding any pre-existing data.

## Usage guide

PhyloSophos analysis could be performed by executing phylosophos_core.py script. You may execute this script as:

<pre><code>python phylosophos_core.py [[optional_parameter_type] [optional_parameter_value]]</code></pre>

PhyloSophos currently recognizes five types of optional parameter types.

* Help (-h, -help, -guide): if one of these arguments is given, a hard-coded guide to PhyloSophos will appear in the console. This will provide simple instructions on how to customize PhyloSophos mapping parameters. No following parameter value is required.
* Reference type change (-r, -ref): if one of these arguments is given, PhyloSophos will change the database of choice to the one specified by the following argument. The default setting is 'ncbi', while 'col' and 'eol' are also available in basic PhyloSophos system. You may change the default setting by modifying /ps_init/ps_initialize.py (see lines 56, 60 & 62). If you want to include other types of references into PhyloSophos system, please read chapter 6.
* Input type change (-i, -input): if one of these arguments is given, along with the name of the input file, PhyloSophos will specifically import the given file as an input. If not (as a default setting), PhyloSophos will consider all files within the /input directory to be scientific name input files.
* Levenshtein distance cutoff (-l, -lev, -cutoff): if one of these arguments is given, along with an integer value, PhyloSophos will change the edit distance cutoff (default setting = 3) to the specified value. 
* Manual curation status (-m, -manual, -curation): if one of these arguments is given, along with a value 1, PhyloSophos will import /pp_learning/manual_curation_list.tsv and utilize this information to pre-process inputs. If not (as a default setting), PhyloSophos will not import extra information other than reference data files within /pp_ref directory.

## Result format guide

The results of the PhyloSophos analysis will be deposited in the /result directory. The name of the result file will be "phylosophos_result_[export_date]\_[export_time]\_[input_file_name]". [export_date] and [export_time] will be six-digit numbers. 

PhyloSophos result file has the following format (delimited by tabs). Each column shows the information about:

* Input_file_name: Name of input file which contains given scientific name input
* Input_original_order: The order in which the scientific name is contained in the file
* Raw_name_input: Raw scientific name input, as it is appeared in the input file
* Pre_corrected_input: Pre-corrected scientific name input
* Chosen_reference: Reference of choice (e.g. 'ncbi', 'col', 'eol')
* Chosen_reference_mapped_ID: Taxonomic entry ID(s) mapped to given scientific name input (based on reference of choice)
* Chosen_reference_scientific_name: Canonical scientific name(s) associated with mapped taxonomic ID(s) (based on reference of choice)
* Chosen_reference_mapping_status_code: PhyloSophos mapping status code (based on reference of choice)
* Chosen_reference_mapping_status_description: Short description of mapping status code (based on reference of choice)
* (specific_reference)\_mapped\_id: Taxonomic entry ID(s) mapped to given scientific name input (based on specific reference)
* (specific_reference)\_scientific\_name: Canonical scientific name(s) associated with mapped taxonomic ID(s) (based on specific reference)
* (specific_reference)\_mapping\_status_code: PhyloSophos mapping status code (based on specific reference)
* Manual_curation_recommended: Mapping status code-based opinion on whether manual curation is needed for this mapping result

Each taxonomic reference found within /pp_ref directory provides [(specific\_reference)\_mapped\_id] - [(specific\_reference)\_scientific\_name] - [(specific\_reference)\_mapping_status\_code] column triplet. Base PhyloSophos provides 3 column triplets for CoL/EoL/NCBI taxonomy respectively.

## Citing PhyloSophos

We recommend that those wishing to cite PhyloSophos use the following citation:
Cho MH, No KT. PhyloSophos: a high-throughput scientific name mapping algorithm augmented with explicit consideration of taxonomic science, and its application on natural product (NP) occurrence database processing. BMC Bioinformatics. (under review: to be updated soon)

## License

PhyloSophos is released with MIT license. This is open to all and can be used for any purpose. We look forward to an increase in the development of this field with the help of contributors who share their derived work publicly.

## Contact

mhcho@bmdrc.org
