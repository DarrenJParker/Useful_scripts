### Kalli_to_edgeR.py

import sys
import os
import getopt
import decimal
from decimal import *

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:o:t:s:Zh')
																						
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
sum_iso_f_name = "NOTHINGSET"

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** Kalli_to_edgeR.py | Written by DJP, 17/07/17 in Python 3.5 in Lausanne ****\n")
		print("This program takes a directory of Kallisto produced directories and produces a count file for EgdeR")
		print("It returns a csv file of counts for use in edgeR, and total counts per sample")
		print("NOTE: Sample names in the count file will be the same as the directory they are in.")
		print("NOTE: Estimates of counts are rounded to integers.")
		print("\n**** USAGE **** \n")
		print("python  Kalli_to_edgeR.py -i input directory -o output_file_base [options] \n")
		print("\n**** USAGE OPTIONS ****\n")
		print("\n-s\tisoform option. Specify gene-to-isoform file here to have the program sum isoform counts to genes. Default OFF\n")
		print("\tgene-to-isoform file:")
		print("\tgene1\tiso12\n\tgene1\tiso122\n\tgene2\tiso1288\n\tgene2\tiso1222")
		
		print("\n-Z\tortholog option. Default OFF. When OFF all count files sould have genes with the exactly the same names, e.g.\n")
		print("\tFile_1:")
		print("\ttarget_id	length	eff_length	est_counts	tpm")
		print("\tOG-1000_Tbi_TRINITY_DN59892_c0_g1_i1	1793	1584.00	70	60.53")
		print("\tOG-1001_Tbi_TRINITY_DN23248_c0_g1_i1	2519	2310.00	52	45.25")
		print("\tOG-1002_Tbi_TRINITY_DN50274_c0_g1_i2	2848	2639.00	74	64.54")
		print("\tFile_2:")
		print("\ttarget_id	length	eff_length	est_counts	tpm")
		print("\tOG-1000_Tbi_TRINITY_DN59892_c0_g1_i1	1793	1584.00	71	64.27")
		print("\tOG-1001_Tbi_TRINITY_DN23248_c0_g1_i1	2519	2310.00	67	60.78")
		print("\tOG-1002_Tbi_TRINITY_DN50274_c0_g1_i2	2848	2639.00	61	55.45")
		print("\n\tWhen -Z is specified (ON) only the first part of the genename has to be the same between samples (everything before the first underscore)e.g.\n")
		print("\tFile_1:")
		print("\ttarget_id	length	eff_length	est_counts	tpm")
		print("\tOG-1000_Tbi_TRINITY_DN59892_c0_g1_i1	1793	1584.00	70	60.53")
		print("\tOG-1001_Tbi_TRINITY_DN23248_c0_g1_i1	2519	2310.00	52	45.25")
		print("\tOG-1002_Tbi_TRINITY_DN50274_c0_g1_i2	2848	2639.00	74	64.54")
		print("\tFile_2:")
		print("\ttarget_id	length	eff_length	est_counts	tpm")
		print("\tOG-1000_Tte_TRINITY_DN57940_c0_g1_i1	1794	1585.00	58	49.89")
		print("\tOG-1001_Tte_TRINITY_DN838_c0_g1_i1	2519	2310.00	64	55.11")
		print("\tOG-1002_Tte_TRINITY_DN23660_c0_g1_i2	2859	2650.00	76	64.74")

		print("\n-t\tortholog splitter option. Everything before this will be used as the gene name. Default is a underscore.\n\n")
		

		sys.exit(2)
	elif opt in ('-i'):
		in_dir_name = arg
	elif opt in ('-o'):
		out_base_name = arg
	elif opt in ('-Z'):
		do_orthos = 'YES'
	elif opt in ('-s'):
		sum_iso_f_name = arg
	elif opt in ('-t'):
		orth_split = arg
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
	
	
if sum_iso_f_name  == "NOTHINGSET":
	print("Not summing isoform counts into gene counts (Default). Use -s option to change this\n")
else:
	print("Summing isoform counts into gene counts using info from: " + sum_iso_f_name + "\n")

	

### if need to sum isoforms up to genes do it BEFORE the next bit (run file through - sum - export - pass to next bit)


if sum_iso_f_name  != "NOTHINGSET":
	iso_to_gene_dict = {}
	iso_set = set()
	gene_set = set()
	sum_iso_f = open(sum_iso_f_name)
	for line in sum_iso_f:
		line = line.rstrip("\n").split("\t")
		gene_id = line[0]
		gene_set.add(gene_id)
		tran_id = line[1]
		if tran_id not in iso_set:
			iso_set.add(tran_id)
		else:
			print("There are duplicate transcript IDs in " + sum_iso_f_name + "! Sort this out. Exiting!")
			sys.exit(2)
			
		iso_to_gene_dict[tran_id] = gene_id
		
		#print(line)
		
	print("Found " + str(len(iso_set)) + " transcript ids from "  + str(len(gene_set)) + " gene ids")
	
	all_sample_summed_by_gene = {}
	samp_l = []
	path = in_dir_name
	for path, subdirs, files in os.walk(path):
		for name in files:
			file_path = os.path.join(path, name)
			if file_path.endswith("abundance.tsv"):	
				#print(file_path)	
				curr_file = open(file_path)
				sample_name = file_path.rstrip("/").split("/")[-2]
				for_out_p = ""
				a2 = file_path.rstrip("/").split("/")
				for i in range(0,len(a2) -2):
					for_out_p = for_out_p + "/" + a2[i]
				for_out_p = for_out_p.lstrip("/")
				
				samp_l.append(sample_name)
				line_N = 0
				
				seen_gene = set()
				for line in curr_file:
					line_N = line_N + 1
					if line_N != 1:
						line = line.rstrip("\n").split("\t")
						trans_name = line[0]
						gene_name = iso_to_gene_dict.get(trans_name)
						gene_name_and_sample = gene_name + "IAMASPLITTER" + sample_name 
						## not round here as round in next part
						cnt_val = decimal.Decimal(line[3])
						
						if gene_name_and_sample not in seen_gene:
							seen_gene.add(gene_name_and_sample)
							all_sample_summed_by_gene[gene_name_and_sample] = cnt_val
						else:
							old_cnt = all_sample_summed_by_gene.get(gene_name_and_sample)
							new_cnt = old_cnt + cnt_val
							all_sample_summed_by_gene[gene_name_and_sample] = new_cnt

	#### open files
	
	for el in samp_l:
		out_file_n = os.path.join(for_out_p,el,"abundance.gene.txt")
		outfile = open(out_file_n, "w")
		outfile.write("target_id\tlength\teff_length\test_counts\ttpm\n")
		#print(out_file_n)
		outfile.close()
		
	for el in all_sample_summed_by_gene:
		g_count = all_sample_summed_by_gene.get(el)
		sample_name = el.split("IAMASPLITTER")[1]
		gene_name =  el.split("IAMASPLITTER")[0]
		out_file_n = os.path.join(for_out_p,sample_name,"abundance.gene.txt")
		outfile = open(out_file_n, "a")
		outfile.write(gene_name + "\tNA\tNA\t"+ str(g_count) + "\tNA\n")
		outfile.close()
		
		
		#print(el)
	
	#print(len(all_sample_summed_by_gene))					
	#print(for_out_p)
	


##### read files and add counts to dict 

count_dict = {}
gene_list = []
first_file = ""
gene_seen = set()
file_N = 0
gene_N_read = []
sample_list = []

if sum_iso_f_name == "NOTHINGSET":
	stat_file_name = out_base_name + "_K2edge.total_counts.txt"
else:
	stat_file_name = out_base_name + "_jgenes_K2edge.total_counts.txt"

stat_file = open(stat_file_name, "w")
stat_file.write("sample_name\ttotal_reads_mapped\n")

path = in_dir_name
for path, subdirs, files in os.walk(path):
	for name in files:
		file_path = os.path.join(path, name)
		
		if sum_iso_f_name == "NOTHINGSET":
			ending = "abundance.tsv"
		else:
			ending = "abundance.gene.txt"
		
		if file_path.endswith(ending):
			curr_file = open(file_path)
			curr_file_gene_N = 0
			file_N = file_N + 1
			total_reads_mapped = 0
			sample_name = file_path.rstrip("/").split("/")[-2]
			sample_list.append(sample_name)
			if file_N == 1:
				first_file = file_path
				# print(file_path)
				line_N = 0
				for line in curr_file:
					line_N = line_N + 1
					if line_N != 1:
						line = line.rstrip("\n").split("\t")
						gene_name = line[0]
						
						## orth option
						if do_orthos == "YES":
							gene_name = gene_name.split(orth_split)[0]
						
						cnt_val = decimal.Decimal(line[3]).to_integral_value()
						if gene_name not in gene_seen:
							gene_seen.add(gene_name)
							count_dict[gene_name] = [cnt_val]
							gene_list.append(gene_name)
							curr_file_gene_N = curr_file_gene_N + 1
							total_reads_mapped = total_reads_mapped + cnt_val
						else:
							print("Gene names (or orth names if using -Z option) are NOT unique in the first file read in (" + file_path + "), Please fix this before continuing, Exiting!" )
							sys.exit(2)
							
			else:
				#print(file_path)
				line_N = 0
				for line in curr_file:
					line_N = line_N + 1
					if line_N != 1:
						line = line.rstrip("\n").split("\t")
						gene_name = line[0]
						
						## orth option
						if do_orthos == "YES":
							gene_name = gene_name.split(orth_split)[0]

						cnt_val = decimal.Decimal(line[3]).to_integral_value()
						rec = count_dict.get(gene_name)
						if rec == None:
							print("There are additional Gene names (or orth names if using -Z option) in " + file_path + " that were not in the first file read in (" +  first_file + "), Please fix this before continuing, Exiting!" )
							sys.exit(2)
						else:
							rec.append(cnt_val)
							count_dict[gene_name] = rec
							curr_file_gene_N = curr_file_gene_N + 1
							total_reads_mapped = total_reads_mapped + cnt_val
							
			gene_N_read.append(curr_file_gene_N)
			# print(curr_file_gene_N)
			# print(total_reads_mapped)
			
			stat_file.write(sample_name + "\t" + str(total_reads_mapped) + "\n")
			
#### check read correct number of genes from each file

gene_N_read_set = set(gene_N_read)
if len(gene_N_read_set) != 1:
	print("Some files do not have the correct number of genes in. Please fix this before continuing, Exiting!" )
	sys.exit(2)
	

print("\nFound " + str(len(sample_list)) + " samples each with " +  str(len(gene_list)) + " genes.\n")

#### output counts

if sum_iso_f_name == "NOTHINGSET":
	count_file_name = out_base_name + "_K2edge.counts"
else:
	count_file_name = out_base_name + "_jgenes_K2edge.counts"

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










