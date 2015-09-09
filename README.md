multi MLST
==============
##Introduction

This pipeline is designed to perform multi-locus sequence typing of prokaryotic genomes. Allele databases are updated with novel alleles as necessary. The
criteria for a novel allele is that it must be the same length and have at least 98% sequence identity to the closest reference allele. Sequence type profiles
can be optionally updated. The updated profiles are not curated - if you use a typing scheme that is not designed for your query organism, then you could be
adding nonsense to your profile database.

##Files
In order for the pipeline to function, the following files must be provided

1. One or more query genomes in fasta format. The file extension must be .fa*. This means that .fa, .fas, .fasta are all acceptable, while .ffn is not.
2. Allele files. These must be in fasta format with a file extension of .*fa*. The addition of the second wildcard allows for the extension used in 
rMLST profiles, .tfa, to be used. The default path the pipeline uses to find allele files is "[path]/Organism/[organismName]/[analysisName]/alleles/*.*fa*". 
3. Profile database. The default path the pipeline uses to find profile database is "[path]/Organism/[organismName]/[analysisName]/profile/*.txt". __Please note that you will not just be able to download profiles from pubmlst.org and immediately run analyses. These default profiles often come with a
column of data with a header called "clonal complex". This column (and any other column other than columns containing the sequence type, and the allele data)
must be removed. Additionally, the first column must be renamed "ST". Otherwise, the pipeline will not be able to parse the sheet correctly. Finally, in order
to ensure that the profile database file has been manipulated, the extension must be changed to .txt.__

##Use
Run mMLST.py from the console.
The following arguments are valid:

* -p', --path, required=True, Specify path for custom folder locations. If you don't supply additional paths e.g. sequencePath, allelePath, 
or organismPath, then the program will look for MLST files in [path]/Organism, and the query sequences in [path]/sequences
* -c, --cutoff, required=False, default=98, The percent identity cutoff value for BLAST matches. Default is 98%
* -s, --sequencePath, required=False, The location of the query sequence files
* -a, --alleleProfilePath, required=False, The path of the folder containing the two folders containing the allele files, and the profile 
file e.g. [path]/Organism/Salmonella/cgMLST Please note the requirements for the profile database above
* -O, --organismPath, required=False, The path of the folder containing the organism folders e.g. [path]/Organism
* -o, --organism, required=False, The name of the organism you wish to type. Must match the folder name containing the schemes e.g. Salmonella
* -S, --scheme, required=False, The scheme you wish to use. Must match the folder name containing the scheme e.g. cgMLST. Furthermore, this folder 
must contain two folders: "alleles" and "profile". The alleles folder contains the allele files in fasta format, and the profile folder contains
the profile in .txt format. Please note the requirements for the profile database above')
* -u, --updateProfileFalse, required=False, default=True, By default, the program automatically creates new sequence profiles and appends these
profiles to the profile file. If, instead, you wish to wish to see the closest match of a query genome to known reference profiles, set this to False
* -U', '--updateAlleleFalse', required=False, default=True, By default, the program automatically creates new allels and appends these alleles to the 
appropriate file. If, instead, you wish to wish to see the closest match of a query genome to known reference alleles, set this to False.')

The pipeline has two modes: interactive and non-interactive. The interactive mode occurs only when insufficient arguments are provided during the running of 
the pipeline. For instance, if the -o and/or -S flags are not provided, then the pipeline will scan the paths provided to try and find the folder containing
the data. For instance, if neither flag is provided, then the pipeline lists all folders in [path]/Organism (or the -O path) and prompts the user to indicate
which organism to use. Once that input is received, the pipeline lists the folders in that directory and prompts the user to select which typing scheme to use.
If only one flag is provided, the pipeline should only prompt the user once, as well as use the provided argument.
                         
                         
Example commands would be:
* mMLST.py -p [path] # Use this if you are using the default folder structure and want an interactive experience 
* mMLST.py -p [path] -o Vibrio -S MLST # Use this if you are using the default folder structure and don't want the program to prompt you for the organism and
the scheme
* mMLST.py -p [path] -s [path1]/sequences -O [path2]/Organism  -o Vibrio -S MLST # Use this if you are not using the default folder structure
* mMLST.py -p [path] -o Vibrio -S MLST -c 95 -u False -U False # Use this if you want to check to see if there are close matches that just fall below 
the 98% cutoff threshold for writing novel alleles. The update profile and allele options are set to false, as you might not want to update the alleles and 
profile database with these lower identity alleles 

##Requirements
* Somewhat powerful Linux workstation - no recommendations on specs, but due to the multi-threaded nature of the script, performance should scale with the system
* Python
    * Biopython
* Blast+

##Outputs
This pipeline generates multiple outputs:

1. The MLST report
2. (Optional) updated allele file(s)
3. (Optional) updated profile database


