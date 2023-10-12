# PhyloSophos

A python-based high-throughput scientific name mapping algorithm

## In a nutshell

```
python phylosophos_initialize_update.py [optional_update_parameter]
python phylosophos_core.py [[optional_parameter_type] [optional_parameter_value]]
```

## Summary

PhyloSophos is a high-throughput scientific name processor which achieves greater mapping performance by referencing multiple taxonomic references and recognizing the semantic structure of scientific names. It also corrects common Latin variants and vernacular names, which often appear in various biological databases and resources. 

Please refer to phylosophos_guide.pdf for detailed information on the package and instructions on how to install it.

## Installation guide

PhyloSophos requires Numpy, along with other basic libraries of Python 3. Before initializing PhyloSophos to fit with your environment, please set up an environment with a python version >=3.7 (preferably with Conda), and install numpy (https://numpy.org/install/) as appropriate.

A script for initialization and update is phylosophos_initialize_update.py. This script will download raw taxonomic metadata from Catalogue of Life (CoL), Encyclopedia of Life (EoL), and NCBI Taxonomy, and then process it into reference files. You may execute this script as:

<pre><code>python phylosophos_initialize_update.py [optional_update_parameter]</code></pre>

For initialization purposes, the optional parameter will not affect downstream processes. It is possible for Catalogue of Life and Encyclopedia of Life to change the name of the metadata file in their FTP server. If this happens, please change file names in the script (see lines 76 & 77) as appropriate.

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

The results of the PhyloSophos analysis will be deposited in the /result directory. The name of the result file will be 'phylosophos_result_[export_date]_[export_time]_[input_file_name]'. [export_date] and [export_time] will be six-digit numbers. 

## Citing PhyloSophos

We recommend that those wishing to cite PhyloSophos use the following citation:
Cho MH, No KT. PhyloSophos: a high-throughput scientific name mapping algorithm augmented with explicit consideration of taxonomic science, and its application on natural product (NP) occurrence database processing. BMC Bioinformatics. (under review: to be updated soon)

## License

PhyloSophos is released with MIT license. This is open to all and can be used for any purpose. We look forward to an increase in the development of this field with the help of contributors who share their derived work publicly.

## Contact

mhcho@bmdrc.org
