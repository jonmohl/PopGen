usage: process_sequences.py [-h] -i DIRNAME -b FILENAME -o DIRNAME
                            [-a FILENAME] [-t INT] [-s STRING] [-PE]

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
  -s STRING, --start_at STRING
                        This option allows the user to start at a particular
                        portion of the analysis pipeline. (1) Full pipeline
                        starting at quality trimming, (2) BWA alignment, (3)
                        Mapdamage, (4) Samtools rmdup, (5) Samtools index and
                        sort, (6) Samtools mpileup. Example: '-s' 3 would
                        start the process at the Mapdamage step.
  -PE, --PairedEnd      Is used for PE fastq files

Please cite: P Lavretsky, et al. (2020) Assessing changes in genomic divergence following a century of human‐mediated secondary contact among wild and captive‐bred ducks. Molecular ecology 29(3):578-595
