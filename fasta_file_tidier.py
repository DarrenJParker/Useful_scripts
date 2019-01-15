## fasta_file_tidier.py

import sys
import os
import getopt
import decimal
import numpy as np
import re


try:
	opts, args = getopt.getopt(sys.argv[1:], 'f:m:o:s:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

input_fasta = None
min_seq_len_thresh = 1
output_fasta_filename = "fasta_tidy_out.fa"
seq_prefix = "contig_"

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** fasta_file_tidier.py | Written by DJP, 15/01/18 in Python 3.5 in Lausanne ****\n")
		print("Takes a fasta file, sorts it by length, and renames sequences sequentially")
		
		print("\n***** USAGE *****\n")		
		print("\npython3 fasta_file_tidier.py -f [fasta file name] -m [min seq length, Default: 1] -s [seq name prefix, Default: contig_] -o [output filename]\n\n")
		
		sys.exit(2)
		
	elif opt in ('-f'):
		input_fasta = arg
	elif opt in ('-m'):
		min_seq_len_thresh = arg
	elif opt in ('-o'):
		output_fasta_filename = arg
	elif opt in ('-s'):
		seq_prefix = arg
	else:
		print("i dont know")
		sys.exit(2)


##### FIRST unwrap fasta - precautionary will be necessary for some files 
### note making a temp unwrapped fasta file  - removed at end

output_fasta_name = input_fasta + ".TEMP_extract_fasta_file" 

output_file = open(output_fasta_name, "w")
print("\nUnwrapping fasta file")
count = 0
in_file = open(input_fasta)
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
		name_list.append(lineB)
	else:
		seq_list.append(lineA)
		done = done + 1
		seq_len = len(lineA)
		len_list.append(seq_len)
			

for element in range(0,len(name_list)):
	name1 = name_list[element]
	seq1 = seq_list[element].replace(" ", "") ## remove gaps if seq comes from gblocks 
	seq_dict[name1] = seq1

#print(seq_dict)
## tidyup
seq_file_1.close()
os.remove(output_fasta_name)

print("Read " + str(done) + " sequences from " + input_fasta)

# print(seq_dict)

###########################################################################################
### sort by len

seq_name_len_list = []

for i in range(0,len(name_list)):
	seq_name = name_list[i]
	seq_len  = len_list[i]
	seq_name_len_list.append((seq_name,seq_len))
	
seq_name_len_list_sorted = sorted(seq_name_len_list, key=lambda x: x[1], reverse=True) 
	

###########################################################################################
### output in length order, renamed, with min seq len (if selected)

output_fasta_file = open(output_fasta_filename, "w")

seq_N = 0
excluded_N = 0

for el in seq_name_len_list_sorted:
	seq_len = el[1]
	if seq_len >= int(min_seq_len_thresh):
		seq_N = seq_N + 1
		new_name = seq_prefix + str(seq_N)
		seq = seq_dict.get(el[0])
		output_fasta_file.write(">" + new_name + "\n" + seq + "\n")
	else:
		excluded_N = excluded_N + 1
		
print("Number of seqs excluded as length less than " + str(min_seq_len_thresh) + ": " +  str(excluded_N))
print("Kept seqs: " + str(seq_N))


print("\n\nFinished Helward\n\n")











