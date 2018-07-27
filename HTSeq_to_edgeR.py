### HTseq_to_edgeR.py

import sys
import os
import getopt
import decimal
from decimal import *

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:o:e:t:Zh')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

in_dir_name = "NOTHINGSET"
out_base_name = "NOTHINGSET"
do_orthos = "NOTHINGSET"
orth_split = "NOTHINGSET"
ending = ".counts"


#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** HTseq_to_edgeR.py | Written by DJP, 24/07/17 in Python 3.5 in Lausanne ****\n")
		print("This program takes a directory count files from HTseq and produces a count file for EdgeR (csv), plus a stats file")
		print("Note I assume the extension of the count files is .counts ***If this is different please use the -e option***")
		print("NOTE: Sample names in the count file will be the same as the filename they are from.")
		
		print("\n**** USAGE **** \n")
		
		print("python3 HTseq_to_edgeR.py -i [input directory] -o [output_file_base] [options] \n")
		print("\n**** USAGE OPTIONS ****\n")
		print("\n-e\tcount file extension. set to whatever the file extension of the count files is. Default is .counts")
		print("\n-Z\tortholog option. Default OFF. When OFF all count files sould have genes with the exactly the same names, e.g.\n")
		print("\tFile_1:")
		print("\tOG-1000_Tbi_TRINITY_DN59892_c0_g1_i1	1793")
		print("\tOG-1001_Tbi_TRINITY_DN23248_c0_g1_i1	2519")
		print("\tOG-1002_Tbi_TRINITY_DN50274_c0_g1_i2	2848")
		print("\tFile_2:")
		print("\tOG-1000_Tbi_TRINITY_DN59892_c0_g1_i1	700")
		print("\tOG-1001_Tbi_TRINITY_DN23248_c0_g1_i1	251")
		print("\tOG-1002_Tbi_TRINITY_DN50274_c0_g1_i2	288")
		print("\n\tWhen -Z is specified (ON) only the first part of the genename has to be the same between samples (everything before the first underscore)e.g.\n")
		print("\tFile_1:")
		print("\tOG-1000_Tbi_TRINITY_DN59892_c0_g1_i1	1793")
		print("\tOG-1001_Tbi_TRINITY_DN23248_c0_g1_i1	2519")
		print("\tOG-1002_Tbi_TRINITY_DN50274_c0_g1_i2	2848")
		print("\tFile_2:")
		print("\tOG-1000_Tte_TRINITY_DN57940_c0_g1_i1	700")
		print("\tOG-1001_Tte_TRINITY_DN838_c0_g1_i1	251")
		print("\tOG-1002_Tte_TRINITY_DN23660_c0_g1_i2	280")


		print("\n-t\tortholog splitter option. Everything before this will be used as the gene name. Default is a underscore.\n\n")
		

		sys.exit(2)
	elif opt in ('-i'):
		in_dir_name = arg
	elif opt in ('-o'):
		out_base_name = arg
	elif opt in ('-Z'):
		do_orthos = 'YES'
	elif opt in ('-t'):
		orth_split = arg
	elif opt in ('-e'):
		ending = arg
	else:
		print("i dont know")
		sys.exit(2)

if in_dir_name == "NOTHINGSET":
	print("Please specify a input dir! For more info see help with option -h")
	sys.exit(2)

if out_base_name == "NOTHINGSET":
	print("Please specify an output base name! For more info see help with option -h")
	sys.exit(2)

if do_orthos == "NOTHINGSET":
	print("\nOrtholog mode is OFF!.")
else:
	print("\nOrtholog mode is ON!.")

if orth_split == "NOTHINGSET":
	orth_split = "_"


file_names_set = set()
path = in_dir_name
for path, subdirs, files in os.walk(path):
	for name in files:
		file_path = os.path.join(path, name)
		
		#print(file_path)
		
		if file_path.endswith(ending):
			file_names_set.add(file_path)

if len(file_names_set) == 0:	
	print("found 0 files ending with " + ending + " in " + in_dir_name + "\n\nNote expected file extention can be specified with -e\n\n\n")
	sys.exit(2)

##### read files and add counts to dict 

count_dict = {}
gene_list = []
first_file = ""
gene_seen = set()
file_N = 0
gene_N_read = []
sample_list = []


stat_file_name = out_base_name + "_H2E_stat_counts.txt"

stat_file = open(stat_file_name, "w")
stat_file.write("sample_name\ttotal_reads\tno_feature\tambiguous\ttoo_low_aQual\tnot_aligned\talignment_not_unique\n")

path = in_dir_name
for path, subdirs, files in os.walk(path):
	for name in files:
		file_path = os.path.join(path, name)
		
		#print(file_path)
		
		if file_path.endswith(ending):
			curr_file = open(file_path)
			curr_file_gene_N = 0
			file_N = file_N + 1
			total_reads = 0
			no_feature = 0
			ambiguous = 0
			too_low_aQual = 0
			not_aligned = 0
			alignment_not_unique = 0		
			sample_name = file_path.rstrip("/").split("/")[-1].split(ending)[0]
			sample_list.append(sample_name)
			if file_N == 1:
				first_file = file_path
				# print(file_path)
				line_N = 0
				for line in curr_file:
					line_N = line_N + 1
					if line_N != 0:
						line = line.rstrip("\n").split("\t")
						gene_name = line[0]
						
						## get non-gene lines
						
						if line[0] == "__no_feature":
							no_feature = line[1]
							total_reads = total_reads + int(line[1])
						elif line[0] == "__ambiguous":
							ambiguous = line[1]
							total_reads = total_reads + int(line[1])
						elif line[0] == "__too_low_aQual":
							too_low_aQual = line[1]
							total_reads = total_reads + int(line[1])
						elif line[0] == "__not_aligned":
							not_aligned = line[1]
							total_reads = total_reads + int(line[1])							
						elif line[0] == "__alignment_not_unique":
							alignment_not_unique = line[1]
							total_reads = total_reads + int(line[1])
						
						else:
							
							## orth option
							if do_orthos == "YES":
								gene_name = gene_name.split(orth_split)[0]
							
							cnt_val = int(line[1])
							if gene_name not in gene_seen:
								gene_seen.add(gene_name)
								count_dict[gene_name] = [cnt_val]
								gene_list.append(gene_name)
								curr_file_gene_N = curr_file_gene_N + 1
								total_reads = total_reads + cnt_val
							else:
								print("Gene names (or orth names if using -Z option) are NOT unique in the first file read in (" + file_path + "), Please fix this before continuing, Exiting!" )
								sys.exit(2)
							
			else:
				#print(file_path)
				line_N = 0
				for line in curr_file:
					line_N = line_N + 1
					if line_N != 0:
						line = line.rstrip("\n").split("\t")
						gene_name = line[0]
						
						## get non-gene lines
						
						if line[0] == "__no_feature":
							no_feature = line[1]
							total_reads = total_reads + int(line[1])
						elif line[0] == "__ambiguous":
							ambiguous = line[1]
							total_reads = total_reads + int(line[1])
						elif line[0] == "__too_low_aQual":
							too_low_aQual = line[1]
							total_reads = total_reads + int(line[1])
						elif line[0] == "__not_aligned":
							not_aligned = line[1]
							total_reads = total_reads + int(line[1])							
						elif line[0] == "__alignment_not_unique":
							alignment_not_unique = line[1]
							total_reads = total_reads + int(line[1])
						
						else:
												
							## orth option
							if do_orthos == "YES":
								gene_name = gene_name.split(orth_split)[0]
	
							cnt_val = int(line[1])
							rec = count_dict.get(gene_name)
							if rec == None:
								print("There are additional Gene names (or orth names if using -Z option) in " + file_path + " that were not in the first file read in (" +  first_file + "), Please fix this before continuing, Exiting!" )
								sys.exit(2)
							else:
								rec.append(cnt_val)
								count_dict[gene_name] = rec
								curr_file_gene_N = curr_file_gene_N + 1
								total_reads = total_reads + cnt_val
								
			gene_N_read.append(curr_file_gene_N)
			# print(curr_file_gene_N)
			# print(total_reads_mapped)
			
			stat_file.write(sample_name + "\t" + str(total_reads)  + "\t" + str(no_feature) + "\t" + str(ambiguous) + "\t" + str(too_low_aQual) + "\t" + str(not_aligned) + "\t" + str(alignment_not_unique) + "\n")
			
			
			
	
#### check read correct number of genes from each file


gene_N_read_set = set(gene_N_read)
if len(gene_N_read_set) != 1:
	print("Some files do not have the correct number of genes in. Please fix this before continuing, Exiting!" )
	sys.exit(2)
	

print("\nFound " + str(len(sample_list)) + " samples each with " +  str(len(gene_list)) + " genes.\n")

#### output counts

count_file_name = out_base_name + "_H2E.counts.csv"

count_file = open(count_file_name, "w")

head_cnt = "Gene_name"
for el in sample_list:
	head_cnt = head_cnt + "," + el 

count_file.write(head_cnt + "\n")

for el in gene_list:
	counts_cs = el
	rec = count_dict.get(el)
	for c in rec:
		counts_cs = counts_cs + "," + str(c)
	
	count_file.write(counts_cs + "\n")
	
	#print(counts_cs)

print("Outputted counts to: " + count_file_name)
print("Outputted total counts to: " + stat_file_name)


print("\n\nFinished, O'Neill.\n\n")


