# FastKraken
Extra fast Kraken2 implementation based on modified scripts from lexinwei and jenniferlu717. Useful for quick and dirty classification of a large number of bins when doing MAG assembly on a dataset where you might expect a variety of prokaryotic or eukaryotic bins, or for when you have a very large database that takes a long time to load. 

# Installation and usage
This script requires you to have kraken2 already installed and in your path.

Download FastKraken from Github. Inside the FastKraken directory is the bin folder, which contains all the scripts you need. The easiest way to get going is to copy the contents to a folder that is in your path. You can also navigate to the FastKraken folder and run the script from inside the folder, specifying the bin as the location of the scripts, like so:

./FastKraken -i folder_of_paired_fastq -d /path/tp/KRAKEN_DB -s R1.fq,R2.fq -o out_test -t 60 -l bin

Options:
```
FastKraken

This script requires that kraken2 be installed in your path.

-i  Directory of sequencing files as input for Kraken2 (REQUIRED)
-d  Path to Kraken2 database (REQUIRED)
-s  Suffix of sequencing files, e.g. 'R1.fastq,R2.fastq' for paired-end files, '.fq' for single-end files. (REQUIRED)
-o  Directory in which to store results (REQUIRED)
-t  Number of threads (OPTIONAL)
-g  Is file gzip compressed? options: true or false (defaults to false). Example: -g true (OPTIONAL)
-l  Specify location of FastKraken scripts if not placed in path (OPTIONAL)
-r  Use read lengths rather than read counts when generating kreport options: true or false (defaults to true). Better option when inputs are contigs or long reads. Example: -g true (OPTIONAL)
```

# About

Kraken2 and its various cousins provide very fast classification for a single file, but the database has to be loaded each time for each file. If you have a very large database, a large number of files, or both, this means that Kraken2 can take quite awhile.

Lexi Wei (lexinwei) created a python script that smushes all sequencing files of interest into a single file or single set or forward/reverse files; the temporary file is run through Kraken2 to create an enormous results file and the results are divided back into separate files afterwards. This allows for your Kraken2 database to only be loaded into memory once. The script makes use of Jennifer Lu's (jenniferlu717) KrakenTools to create a taxonomy file (usually already present in your database) and convert your results into a Kraken report; however, when you have a large number of files that you've classified, this again takes quite awhile because the taxonomy db must be loaded into memory each time as a python dictionary. My modifications to both scripts allows for loading the taxonomy db into memory only once for creating the kreports, even if you have hundreds of files.


