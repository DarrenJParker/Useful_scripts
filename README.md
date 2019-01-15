# Useful_scripts

Some potentially useful scripts.

If you use any of the scripts in this repository please cite it as follows:

Parker, D.J. 2018. Useful_scripts. _GitHub repository_: https://github.com/DarrenJParker/Useful_scripts/

NOTE: help messages can be displayed for all scripts by specifying -h (i.e. python [name of script] -h) .
All scripts written in python3.4 or higher.

## Scripts

* **Kalli_to_edgeR.py** | Takes a directory of Kallisto produced directories a single csv file for use in EdgeR (or similar), along with a stat file.

* **HTSeq accessory scripts**
  * **Maker_gff_to_HTseq_gff.py** | Takes a gff produced by Maker2 and converts it for use in HTseq.
  * **HTSeq_to_edgeR.py** | Takes a directory of read count files from HTseq and produces a single csv file for use in EdgeR (or similar), along with a stat file.

* **gff_feature_lengths.py** | Gets lengths of features in a gff file, summed by a parent ID.

* **fasta_file_tidier.py** | Takes a fasta file, orders sequences by size, and renames them sequentially. Also filters out small contigs if required.
