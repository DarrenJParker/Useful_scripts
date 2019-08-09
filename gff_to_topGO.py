import sys
import os
import getopt


try:
	opts, args = getopt.getopt(sys.argv[1:], 'g:s:t:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

gff_file_name = None
GO_split_value = None
gene_name_tag = "Name="

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** gff_to_topGO.py | Written by DJP, 23/10/18 in Python 3.5 in Lausanne, Swiss ****\n")
		print("Takes a gff file with GO terms and adds makes a file for TopGO")
		print("\n**** Usage ****\n")
		print("python3 gff_to_topGO.py -g [gff file] -s [GO term tag] [options]")
		print("\n**** OPTIONS ****\n")
		print("-t\tGene tag: tag for the gene name in the gff file. Default = Name")
		print("-s\tGO term tag: tag for the GO terms in the gff file.")
		print("-g\tgff file name\n\n\n")
								

		sys.exit(2)
		

	elif opt in ('-g'):
		gff_file_name = arg
	elif opt in ('-s'):
		GO_split_value = arg
	elif opt in ('-t'):
		gene_name_tag = arg
	else:
		print("i dont know")
		sys.exit(2)

gene_name_tag = gene_name_tag.rstrip("=") + "="
print("\ngene_name_tag: " + gene_name_tag)

if GO_split_value == None:
	print("\n\nERROR: No GO tag provided, please provide using -s\n\n")
	sys.exit(2)

GO_split_value = GO_split_value.rstrip("=") + "="
print("\nGO_tag used: " + GO_split_value)



outfile_name = gff_file_name + "_" + GO_split_value.rstrip("=") + "_fortopgo.txt"
outfile      = open(outfile_name, "w")

gff_file = open(gff_file_name)


N_with_GO = 0
N_without_GO = 0

for line in gff_file:
	line = line.rstrip("\n")
	line = line.split("\t")
	feature = line[2]
	desc_line = line[8]
	if feature == "gene":
		gene_name = desc_line.split(gene_name_tag)[1].split(";")[0]
		GO_line = desc_line.split(GO_split_value)
		
		if len(GO_line) == 1:
			GO_line = ""
			N_without_GO = N_without_GO + 1
		else:
			N_with_GO = N_with_GO + 1
			GO_line = GO_line[1].split(";")[0]
		
		GO_line = GO_line.replace(",", "\t")

		outfile.write(gene_name + "\t" + GO_line + "\n")


print("Number of genes with GO annot: " + str(N_with_GO))
print("Number of genes without GO annot: " + str(N_without_GO))

print("\n\nFinished Jeli\n\n\n")











