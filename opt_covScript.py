#!/usr/bin/env python3

#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=SAMtoBaM      ### Job Name
#SBATCH --output=SAMtoBaM.out         ### File in which to store job output
#SBATCH --error=SAMtoBaM.err          ### File in which to store job error messages
#SBATCH --time=0-00:60:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --account=bgmp      ### Account used for job submission
#SBATCH --cpus-per-task=5

#conda deactivate
#conda deactivate 
#conda deactivate 
#conda deactivate 
#conda deactivate  
#conda activate bgmp_py3

file1="/projects/bgmp/msconce/PS6/800_3_PE5_interleaved.fq_1"
file2="/projects/bgmp/msconce/PS6/800_3_PE5_interleaved.fq_2"
file3="/projects/bgmp/msconce/PS6/800_3_PE5_interleaved.fq.unmatched"


#a. Our fosmid library is comprised of 50 fosmids, each (approximately) 40 Kb long. 
#How many total nt should be in this fosmid library?
fosmid_library=50
fosmid_size=40000
expected_genome_ln=fosmid_library*fosmid_size
coverage=()
#/usr/bin/time -v 

#b. Coverage is a metric applicable to total fosmid library size just as it is applicable to 
#genome (estimated) size. Think of all of these fosmids combined as the genome we’re trying 
#to assemble. Calculate the expected coverage.
#	i. NOTE: Reads are varying lengths due to adapter clipping and quality trimming. You must calculate the total number of nt in the input FASTQ files in order to calculate the expected coverage.
nucleotide_length_file1=[]
with open(file1, "r") as fh:
	i = 0
	for line in fh:
		i+=1
		line = line.strip('\n')
		if i%4 ==2:
			#print(len(line))
			nucleotide_length_file1.append(len(line))
#print(sum(nucleotide_length_file1))

nucleotide_length_file2=[]
with open(file2, "r") as fh:
	i = 0
	for line2 in fh:
		i+=1
		line2 = line2.strip('\n')
		if i%4 ==2:
			#print(len(line))
			nucleotide_length_file2.append(len(line2))

nucleotide_length_file3=[]
with open(file3, "r") as fh:
	i = 0
	for line3 in fh:
		i+=1
		line3 = line3.strip('\n')
		if i%4 ==2:
			#print(len(line))
			nucleotide_length_file3.append(len(line3))

file1_nuc_tot=sum(nucleotide_length_file1)
print("file 1 nucleotide total is:", file1_nuc_tot)
file2_nuc_tot=sum(nucleotide_length_file2)
print("file 2 nucleotide total is:", file2_nuc_tot)
file3_nuc_tot=sum(nucleotide_length_file3)
print("file 3 nucleotide total is:", file3_nuc_tot)

all_files_tot=(file1_nuc_tot)+(file2_nuc_tot)+(file3_nuc_tot)
print("Total nucleotides:", all_files_tot)

#calculating expected k-mer coverage
#c. Given the calculated coverage from 3.b (all_files_tot) and total fosmid library size (expected_genome_ln), 
#calculate the k-mer coverage.

#Ck =C*(L-K+1)/L
#Ck: k-mer coverage
#C: coverage
#L: length of reads (bp) 
#K: k-mer length

#Ck1 =C*(L1 -K+1)/L1 
#Ck2 =C*(L2 -K+1)/L2 
#Ck3 =C*(L3 -K+1)/L3
#Ck = W1Ck1 + W2Ck2 + W3Ck3


 #W: percent of dataset composed
     #of reads of length k1

k=40

k_cov_1=float(all_files_tot*(len(nucleotide_length_file1)-k+1)/len(nucleotide_length_file1))
k_cov_2=float(all_files_tot*(len(nucleotide_length_file2)-k+1)/len(nucleotide_length_file2))
k_cov_3=float(all_files_tot*(len(nucleotide_length_file3)-k+1)/len(nucleotide_length_file3))

CK = w1*(k_cov_1)+w2*(k_cov_2)

kmer_coverage=float(all_files_tot/expected_genome_ln)
print("k-mer coverage from calc. coverage and fosmid library size:", kmer_coverage)

#calculating average coverage/file
ave_length_file1=all_files_tot/len(nucleotide_length_file1)
ave_length_file2=all_files_tot/len(nucleotide_length_file2)
ave_length_file3=all_files_tot/len(nucleotide_length_file3)

average=(ave_length_file1+ave_length_file2+ave_length_file3)/3
print("average coverage per file:", average)


#e. Run velveth/velvetg with k-mer sizes of 31, 41, and 49.

#f. Use your code from Part 1 to collect assembly statistics on each result.

#make ’MAXKMERLENGTH=57’



