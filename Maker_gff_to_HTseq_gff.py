### Maker_gff_to_HTseq_gff.py

import sys
import os
import getopt
import decimal


try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:o:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


in_file_name = "NOTHINGSET"
out_prefix = "NOTHINGSET"

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** Maker_gff_to_HTseq_gff.py | Written by DJP, 23/07/18 in Python 3.5 in Lausanne, Swiss ****\n")
		print("Takes gffs from the Maker pipeline and converts them so they can be used properly with HTseq")
		print("Maker gffs do not work properly with HTseq as parent IDs for features below mRNA (exon, CDS, UTR) point to the mRNA rather than to the gene. \ni.e. there is no gene ID for these features in the same line")		
		print("This program adds a gene ID to these features.")
		
		print("\n**** USAGE **** \n")
		print("python3 Maker_gff_to_HTseq_gff.py -i [Maker gff file] -o [out prefix] \n")
		
		print("\n**** Use in HTseq ****\n")
		print("\n**** !!!NOTE the way to access HTseq may differ on your cluster! ****\n")		
		print("After conversion:\n\tfor gene level counts for reads mapping to exons (not introns!) (HTseq's default run mode when given a gtf):")
		print("\tpython2.7 -m HTSeq.scripts.count --order=name --stranded=reverse --format=bam --type=exon --idattr=gene_id my.sorted.bam converted_gff_forHTSeq.gff > my.counts")		
		print("\n\tfor gene level counts for reads mapping within genes (incl. introns!):")
		print("\tpython2.7 -m HTSeq.scripts.count --order=name --stranded=reverse --format=bam --type=gene --idattr=ID my.sorted.bam converted_gff_forHTSeq.gff > my.counts\n\n\n")		

		sys.exit(2)
		
	elif opt in ('-i'):
		in_file_name = arg
	elif opt in ('-o'):
		out_prefix = arg
	else:
		print("i dont know")
		sys.exit(2)



### read gff | get link between transcript name and gene name

trans_to_gene_dict = {}

gene_names = set()
trans_names = set()

in_gff = open(in_file_name)
for line in in_gff:
	line = line.rstrip("\n")
	feature = line.split("\t")[2]
	if feature == "mRNA":
		descrip_line = line.split("\t")[8]
		gene_name = descrip_line.split("Parent=")[1].split(";")[0]
		trans_name = descrip_line.split("ID=")[1].split(";")[0]
		gene_names.add(gene_name)
		trans_names.add(trans_name)
		trans_to_gene_dict[trans_name] = gene_name
in_gff.close()		


print("\nNumber of genes: " + str(len(gene_names)))
print("Number of mRNAs: " + str(len(trans_names)))


### go through gff and add gene_name to all features that are not genes

out_gff_file =  open(out_prefix + "_forHTSeq.gff", "w")
skipped_features = set()


in_gff = open(in_file_name)
for line in in_gff:
	line = line.rstrip("\n")
	feature = line.split("\t")[2]
	descrip_line = line.split("\t")[8]
	if feature == "gene":
		out_gff_file.write(line.rstrip(";")  + ";" + "\n")
	elif feature == "mRNA":
		trans_ID = descrip_line.split("ID=")[1].split(";")[0]
		gene_name =  trans_to_gene_dict.get(trans_ID)
		if gene_name == None:
			print("Something has gone wrong, Exiting!")
			sys.exit(2)
		
		out_gff_file.write(line.rstrip(";") + ";" + "gene_id=" + gene_name + ";" + "\n")
			
			
	else:
		try:
			parents = descrip_line.split("Parent=")[1].split(";")[0].split(",")
			
			#print(parents)
			gene_name = ""
			
			if len(parents) == 1:
				parent = parents[0]
				gene_name =  trans_to_gene_dict.get(parent)
			else:
				### check features come from the same gene when multiple mRNAs
				
				parent_set = set()
				for el in parents:
					gene_name_t =  trans_to_gene_dict.get(el)	
					parent_set.add(gene_name_t)
				if len(parent_set) != 1:
					print("Features have different gene parents..., Exiting!")
					sys.exit(2)
					
				parent = parents[0]	
				gene_name =  trans_to_gene_dict.get(parent)
			
			if gene_name == None:
				print("Something has gone wrong, Exiting!")
				sys.exit(2)
	
			if gene_name == "":
				print("Something has gone wrong, Exiting!")
				sys.exit(2)
				
			out_gff_file.write(line.rstrip(";") + ";" + "gene_id=" + gene_name + ";" + "\n")
		
		## SKIp features that do not have Parent IDs (these should be masking annotations)
		except:
			out_gff_file.write(line.rstrip(";") + ";" + "\n")
			skipped_features.add(feature)

print("\nThe following features were skipped as the have no Parent ID. (Note - only really worry if any of these are exon. having 'match' listed here is expected)")
print(skipped_features)			
	
print("\n\n\nFinished, Shadow\n\n\n")
