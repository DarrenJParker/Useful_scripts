# gff_feature_lengths.py

import sys
import os
import getopt
import decimal


try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:o:f:p:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


in_file_name = "NOTHINGSET"
out_prefix = "NOTHINGSET"
feature_want = "exon"
parent_ID = "gene_id"


#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** gff_feature_lengths.py | Written by DJP, 30/07/18 in Python 3.5 in Lausanne, Swiss ****\n")
		print("Gets lengths of features in a gff file, summed by a group ID (e.g. exons lengths by gene)")
		
		print("\n**** USAGE **** \n")
		print("gff_feature_lengths.py -i [gff file] -o [out prefix] [options] \n")
		
		print("\n**** OPTIONS ****\n")
		print("-f\tfeature type wanted - Default: exon")
		print("-p\tparent ID (what ID should be used to group lengths by (in the description column) - Default: gene_ID")
		
		print("\n**** common options ****\n")
		print("For exons in genes (from a gff processed by Maker_gff_to_HTseq_gff.py):\ngff_feature_lengths.py -i [gff file] -o [out prefix] -f exon -p gene_id\n")
		print("For exons in mRNAs:\ngff_feature_lengths.py -i [gff file] -o [out prefix] -f exon -p Parent\n")		
		print("For gene lengths: \ngff_feature_lengths.py -i [gff file] -o [out prefix] -f gene -p ID\n\n")		

		sys.exit(2)
		
	elif opt in ('-i'):
		in_file_name = arg
	elif opt in ('-o'):
		out_prefix = arg
	elif opt in ('-f'):
		feature_want = arg
	elif opt in ('-p'):
		parent_ID = arg
	else:
		print("i dont know")
		sys.exit(2)

feature_len_dict = {}

grouping_feature_names = set()
N_features = 0

in_gff = open(in_file_name)
for line in in_gff:
	line = line.rstrip("\n")
	feature = line.split("\t")[2]
	if feature == feature_want:
		descrip_line = line.split("\t")[8]
		
		if len(descrip_line.split(parent_ID + "=")) == 1 :		
			print("\n" + parent_ID + " not found in all " + feature_want + " lines. Please sort this.\n\nExiting\n\n\n")
			sys.exit(2)
		grouping_feature = descrip_line.split(parent_ID + "=")[1].split(";")[0]
		
		
		N_features = N_features + 1
		
		start_coord = int(line.split("\t")[3])
		end_coord   = int(line.split("\t")[4])

		
		if start_coord > end_coord:
			print("Start coord is larger than end coord. I do not know how to deal with this... \n\nExiting.\n\n\n")
			sys.exit(2)
			
		else:
			
			feat_len = (end_coord - start_coord) + 1
			# print(line)
			# print(feat_len)
			
			## add to dict
			
			if grouping_feature not in grouping_feature_names:
				feature_len_dict[grouping_feature] = [feat_len, 1]
				grouping_feature_names.add(grouping_feature)
			else:
				rec_len = feature_len_dict.get(grouping_feature)[0]
				new_len = rec_len + feat_len
				N_feature = feature_len_dict.get(grouping_feature)[1] + 1
				feature_len_dict[grouping_feature] = [new_len,N_feature]
				#print(new_len)
				
in_gff.close()		

print("\nNumber of unique " + parent_ID + "s: " + str(len(grouping_feature_names)))
print("Total number of " + feature_want + "s: " + str(N_features))


#####################################################
### output

out_file = open(out_prefix + "_" + feature_want + "_by_" + parent_ID + "_lengths.csv", "w")
out_file.write("Id,total_" + feature_want + "_length,N_" + feature_want + "_length\n")

for el in feature_len_dict:
	f_len = feature_len_dict.get(el)[0]
	n_feat = feature_len_dict.get(el)[1]
	out_file.write(el + "," + str(f_len) + "," + str(n_feat) +  "\n")


print("\n\nFinished, Elise.\n\n")










