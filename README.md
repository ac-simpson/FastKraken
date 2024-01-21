# FastKraken
Extra fast Kraken2 implementation based on modified scripts from lexinwei and jenniferlu717. Especially useful for quick and dirty classification of a large number of bins when doing MAG assembly on a dataset where you might expect eukaryotic bins.

Kraken2 and its various cousins provide very fast classification for a single file, but the database has to be loaded each time for each file. If you have a very large database, a large number of files, or both, this means that Kraken2 can take quite awhile.

Lexi Wei (lexinwei) created a python script that smushes all sequencing files of interest into a single file or single set or forward/reverse files; the temporary file is run through Kraken2 to create an enormous results file and the results are divided back into separate files afterwards. This allows for your Kraken2 database to only be loaded into memory once. The script makes use of Jennifer Lu's (jenniferlu717) KrakenTools to create a taxonomy file (usually already present in your database) and convert your results into a Kraken report; however, when you have a large number of files that you've classified, this again takes quite awhile because the taxonomy db must be loaded into memory each time as a python dictionary. My modifications to both scripts allows for loading the taxonomy db into memory only once for creating the kreports, even if you have hundreds of files.

This script requires you to have kraken2 already installed and in your path.
