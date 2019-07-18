#!/usr/bin/env python3

#Our goal with this assignment is to generate some measures of accuracy for whole genome 
#assemblies and then use these measures to assay our success in several Velvet assemblies.

#use argparse to call out your arguments  
import argparse
def get_arguments():
    parser = argparse.ArgumentParser(description="K-mer Nomalization")
    parser.add_argument("-f", "--filename", help="Filename should be a string", required=True, type=str)
    parser.add_argument("-nf", "--new_file", help="Type new output file name", required=False, type=str)
    return parser.parse_args()

#global variables                         
args = get_arguments() #this will print the argument given to -k or --kmer_size

f = args.filename 
nf = args.new_file

#1. Parse the contigs.fa file that is output by velvetg. Extract the FASTA ID lines as 
#you parse the file (remember: these strings will begin with the “>” character).

#2. You can use the sample data contigs.fa from Talapas to test your code.

#3. Using Python regular expressions, extract k-mer length of each contig (in red below). In
#addition, extract the k-mer coverage for the contig (in blue). Assume a k-mer length of 49.
#>NODE_11_length_3717_cov_19.315845

#example     >NODE_1_length_1154_cov_536.720947

kmer_length = []
kmer_dist = []
contigs_length_list = []
cov_blue_list = []
contigs_N50_list = []
import re

#to count number of contigs within the file
contig_counter = 0

with open(f, "r") as fh:
	i=1
	for line in fh:
		#pattern = '(\A>[A-Z]+_[0-9]+_[a-zA-Z]+_[0-9]+_[a-zA-Z]+_[0-9]+.[0-9]+)'
		#pattern = r'NODE_2_length_822_cov_81.681267'
		#use r'' for REs, note >NODE_, _length_, and _cov_ are same in each line
		#() groups what you want for output ex: (\d+) gets digits
		pattern = r'^>NODE_\d+_length_(\d+)_cov_(\d+\.\d+)$'
		pattern = re.findall(pattern, line)
		#pattern = re.match(pattern, line) #creates tuple, want only the tuples with numbers	
		i+=1
		#print(pattern)
		j=0
		for x in pattern:
			#print(x)
			kmer_length = x[0]
			kmer_length = int(kmer_length)
			#print(kmer_length)
			cov = x[1]
			cov = float(cov)
			#Adjust the k-mer length to represent the physical length
			#physical length (nucleotides) from kmer length ex: adapt (kcnt = length - kmer-size + 1)
			physical_length = kmer_length + 48 
			#print(physical_length)
			#print(contigs)
			#print(type(contigs))
			contigs_length_list.append(physical_length)
			contigs_length_list.sort()
			cov_blue_list.append(cov)
			contigs_N50_list.append(physical_length)
			contigs_N50_list.sort(reverse=True)
			



			#print(cov_blue_list)

			j+=1
			contig_counter+=1
			
#to calculate n50			
total = 0
for contigs in contigs_N50_list:
	total+=contigs
	#print(total)
	#print(contigs)
	if total >= sum(contigs_N50_list)/2:
		N50=contigs
		break
print("n50:", contigs)


	




genome = sum(contigs_length_list)


tot_coverage = (sum(contigs_length_list)*contig_counter)/genome



print("Maximum contig length:", contigs_length_list[-1]) 
print("Mean length of contigs:", sum(contigs_length_list)/contig_counter)
print("Length of genome across contigs:", genome)
print("Number of contigs:", contig_counter)
print("Mean depth of coverage for contigs:", sum(cov_blue_list)/contig_counter)

#Calculate the N50 value of your assembly.









#grep "^>" tester.fa in command line to see what the header lines look like

# Program to extract numbers from a string
#import re
#string = 'hello 12 hi 89. Howdy 34'
#pattern = '\d+'
#result = re.findall(pattern, string) 
#print(result)
# Output: ['12', '89', '34']




	