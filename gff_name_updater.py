## gff_name_updater.py

import sys
import os
import getopt
import decimal
import numpy as np
import re


try:
	opts, args = getopt.getopt(sys.argv[1:], 'f:n:o:g:SFh')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

old_fasta_name = None
new_fasta_name = None
old_gff_name   = None
new_gff_name   = "New.gff"
remove_after_space = "NO"
add_fasta_at_the_end = "NO"

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** gff_name_updater.py | Written by DJP, 16/02/18 in Python 3.5 in Lausanne ****\n")
		print("Updates the contig names of gff file made on one fasta file with those in a new fasta file (e.g. a renamed subset of the original file)")
		
		print("\n***** USAGE *****\n")		
		print("\npython3 gff_name_updater.py -f [original fasta file] -n [new fasta file] -g [original gff] -o [output gff name] [options]\n\n")

		print("\n***** Options *****\n")		
		print("\t-S\tremove_after_space: Discarded information in fasta header following a space. [Default: OFF]")
		print("\t-F\tadd fasta sequence at the end of the new gff. [Default: OFF]\n\n\n")	
		
		sys.exit(2)
		
	elif opt in ('-f'):
		old_fasta_name = arg
	elif opt in ('-n'):
		new_fasta_name = arg
	elif opt in ('-o'):
		new_gff_name  = arg
	elif opt in ('-g'):
		old_gff_name = arg
	elif opt in ('-S'):
		remove_after_space = "YES"
	elif opt in ('-F'):
		add_fasta_at_the_end = "YES"
	else:
		print("i dont know")
		sys.exit(2)

if remove_after_space == "YES":
	print("\n\nremove_after_space option (-S) specified. Information in fasta header following a space will be discarded.")


### read old fasta file in
##### FIRST unwrap fasta - precautionary will be necessary for some files 
### note making a temp unwrapped fasta file  - removed at end

output_fasta_name = old_fasta_name + ".TEMP_extract_fasta_file" 

output_file = open(output_fasta_name, "w")
print("\nUnwrapping fasta file")
count = 0
in_file = open(old_fasta_name)
for line in in_file:
	count = count + 1
	line = line.rstrip("\n")
	if line.startswith(">") and count == 1:
		output_file.write(line + "\n")
	elif line.startswith(">") and count > 1:
		output_file.write("\n" + line + "\n")
	else: 
		output_file.write(line)	

output_file.close()

### add seqs to dictionary
name_list = []
seq_list = []
len_list = []
seq_dict = {}

done = 0
seq_file_1 = open(output_fasta_name)
for line in seq_file_1:
	lineA = line.rstrip("\n")
	if lineA.startswith(">"):
		lineB = lineA.replace(">", "")
		
		if remove_after_space == "YES":
			name_list.append(lineB.split(" ")[0])			
		else:
			name_list.append(lineB)
	else:
		seq_list.append(lineA)
		done = done + 1
		seq_len = len(lineA)
		len_list.append(seq_len)
			

for element in range(0,len(name_list)):
	name1 = name_list[element]
	seq1 = seq_list[element].replace(" ", "").upper() ## remove gaps if seq comes from gblocks 
	seq_dict[name1] = seq1

#print(seq_dict)
## tidyup
seq_file_1.close()
os.remove(output_fasta_name)

print("Read " + str(done) + " sequences from " + old_fasta_name)



## check none of the seqs names are duplicates in the old fasta (this is a code assumption)

if len(name_list) != len(seq_dict):
	print("Error! Some sequences in the old fasta file are duplicates! Maybe you are using the -S option and shouldn't be?\n\nExiting\n\n")
	sys.exit(2)

### read new fasta file in
##### FIRST unwrap fasta - precautionary will be necessary for some files 
### note making a temp unwrapped fasta file  - removed at end


output_fasta_name = new_fasta_name + ".TEMP_extract_fasta_file" 

output_file = open(output_fasta_name, "w")
print("\nUnwrapping fasta file")
count = 0
in_file = open(new_fasta_name)
for line in in_file:
	count = count + 1
	line = line.rstrip("\n")
	if line.startswith(">") and count == 1:
		output_file.write(line + "\n")
	elif line.startswith(">") and count > 1:
		output_file.write("\n" + line + "\n")
	else: 
		output_file.write(line)	

output_file.close()


### add seqs to dictionary
name_list = []
seq_list = []
len_list = []
seq_dict_new = {}
seq_dict_new2 = {}

done = 0
seq_file_1 = open(output_fasta_name)
for line in seq_file_1:
	lineA = line.rstrip("\n")
	if lineA.startswith(">"):
		lineB = lineA.replace(">", "")
		if remove_after_space == "YES":
			name_list.append(lineB.split(" ")[0])			
		else:
			name_list.append(lineB)
	else:
		seq_list.append(lineA)
		done = done + 1
		seq_len = len(lineA)
		len_list.append(seq_len)
			

for element in range(0,len(name_list)):
	name1 = name_list[element]
	seq1 = seq_list[element].replace(" ", "").upper() ## remove gaps if seq comes from gblocks 
	seq_dict_new[seq1] = name1 ## seq is the key 
	seq_dict_new2[name1] = seq1

## tidyup
seq_file_1.close()
os.remove(output_fasta_name)

print("Read " + str(done) + " sequences from " + new_fasta_name)


## check none of the seqs are duplicates in the new fasta (this is a code assumption)

if len(name_list) != len(seq_dict_new):
	print("Error! Some sequences in the new fasta file are duplicates! \n\nExiting\n\n")
	sys.exit(2)


##################################################################
## update gff

genes_in_old_gff = 0
genes_in_new_gff = 0

contigs_with_annot = []
seen_contigs = set()

new_gff = open(new_gff_name, "w")
old_gff = open(old_gff_name)
for line in old_gff:
	line = line.strip("\n").split("\t")
	contig_name  = line[0]
	contig_seq = seq_dict.get(contig_name)
	
	if len(line) > 3:
		seq_type = line[2]
		if seq_type == "gene":
			genes_in_old_gff = genes_in_old_gff + 1
	
	if contig_seq != None:
		
		
		if seq_type == "gene":
			genes_in_new_gff = genes_in_new_gff + 1
		new_contig_name = seq_dict_new.get(contig_seq)
		
		new_line = ""
		for i in range(1,len(line)):
			new_line = new_line + "\t" + line[i]
			
		new_line = new_contig_name + new_line
		new_gff.write(new_line + "\n")
		
		if new_contig_name not in seen_contigs:
			seen_contigs.add(new_contig_name)
			contigs_with_annot.append(new_contig_name)
		
		

print("\nNumber of genes in " + old_gff_name + " = " + str(genes_in_old_gff))
print("Number of genes in " + new_gff_name + " = " + str(genes_in_new_gff))

#############################################################################
## output wrapped fasta at the end


	
if add_fasta_at_the_end == "YES":
	print("\nAdding fasta sequence to end of " + new_gff_name)
	new_gff.write("##FASTA\n")
		
def insert_newlines(string, every=60):
    return '\n'.join(string[i:i+every] for i in range(0, len(string), every))

if add_fasta_at_the_end == "YES":
	for el in contigs_with_annot:
		seq = insert_newlines(seq_dict_new2.get(new_contig_name))
		new_gff.write(">" + el + "\n" + seq + "\n")
		



#### 


print("\n\nFinished Wallfacer Luo\n\n")







