This version was developed for the following citation:
P Lavretsky, et al. (2020) Assessing changes in genomic divergence following a century of human‐mediated secondary contact among wild and captive‐bred ducks. Molecular ecology 29(3):578-595

usage: process_sequences.py [-h] -i DIRNAME -b FILENAME -o DIRNAME
                            [-a FILENAME] [-t INT]

Processing Bait Sequencing Files

optional arguments:
  -h, --help            show this help message and exit
  -i DIRNAME, --input_directory DIRNAME
                        Required: Input directory with raw fastq files.
  -b FILENAME, --bait_reference FILENAME
                        Required: Fasta sequence containing bait sequences.
  -o DIRNAME, --output_dir DIRNAME
                        Required: Name of output directory to create.
  -a FILENAME, --ancient_list FILENAME
                        This file contains the names of the files that may
                        need the MapDamage correction
  -t INT, --threads INT
