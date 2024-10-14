# Goal: simulate RE digest of A. thaliana genome
# output: fasta file of fragments of the correct size, bed file of fragments of the correct size, some sort of .csv file with per chromosome coverage 

import re
import copy 
import numpy as np 
import pandas as pd 
import os
import glob as glob
import string

#for the test, start with one enzyme and one chromosome 
parsed={}

#requried function from rapyds 
def parse_enzymedb(enzyme_db_file):
	"""
	function to parse the enzyme database file.
	Returns a dictionary of restriction enzyme and its site.
	input format: each enzyme with its restriction site in separate lines
		ex.	SbfI,CCTGCA|GG
			ApeKI,G|CWGC
	"""
	list_enzymes = {}
	input_db = open(enzyme_db_file, "r+")
	line_no = 1
	for line in input_db:
		line = line.strip().rstrip().split(",")		## strip strip and split

		## catch in RE DB
		if(len(line)!=2):
			print("Error in restriction enzyme database file line no "+str(line_no))
			raise SystemExit
		else:
			search=re.compile(r'[GCATNMRWYSKHBVD]+[|]+[GCATNMRWYSKHBVD]+').search
			if(bool(search(line[1])) == False):
				print("Error in restriction enzyme database file line no "+str(line_no)+". Invalid letters in sequence.")
				raise SystemExit
			else:
				line_no+=1
				list_enzymes[line[0]] = line[1]

	# if there enzymes parsed in the database file, raise an error
	if(len(list_enzymes) < 0):
		print("No restriction enzymes found in "+enzyme_db_file)
		raise SystemExit

	return list_enzymes

def parse_input(input_name):
	"""
		function that parses the input sequence/genome file (fasta format).
		Returns a list of [0] sequence name and its corresponding [1] genome/dna sequence
	"""
	input_file  = open(input_name, "r+")

	new_genome_name = ""
	new_genome_seq = ""
	list_genome = []
	for line in input_file:
		line = line.strip().rstrip()
		if(len(line) == 0):
			continue
		if(">" in line):	## append to list after seeing ">"
			list_genome.append([new_genome_name,new_genome_seq])
			new_genome_seq =""
			new_genome_name = line
		else:
			new_genome_seq += line.upper().strip().rstrip()	## strip whitespaces
	list_genome.append([new_genome_name,new_genome_seq])
	return list_genome[1:]


#load enzyme db 
parsed['db']=parse_enzymedb('/global/scratch/users/chandlersutherland/e19/re_digest/abascal_enzymes.txt')
#load genome 
genome=parse_input('/global/scratch/users/chandlersutherland/e19/re_digest/input_genome/Athaliana_clean.fa')

#adapt raypds
def digest(enzyme, chromosome):
	enzyme_regex=parsed['db'][enzyme]
	p5, p3 = enzyme_regex.split("|")
	fragments=re.split(p5+p3, chromosome[1]) #chr 1 fragments 

	#add index information 
	index = 0
	start = 0
	curr_len = start				## temporary holder for current base in genome
	new_fragments = []

	p5_orig = copy.copy(p5)
	p3_orig = copy.copy(p3)

	for i in range(0,len(fragments)-1):
		temp_frag = []			## temporary list holder for current fragment
		temp_start = curr_len	## take note of curr_len, it will be the start location of the fragment
		index += len(fragments[i])
		curr_len += len(fragments[i])

		temp_p5 = copy.copy(p5_orig)
		temp_p3 = copy.copy(p3_orig)

		#at least for Alu, with no degenerate sequences, not required
		## Add the 5' part
		p5_replace = [m.start() for m in re.finditer('[NMRWYSKHBVD]', p5_orig)]

		for replace_index in p5_replace:
			temp_p5 = temp_p5[:replace_index] + chromosome[index+replace_index] + temp_p5[replace_index+1:]

		fragments[i] = fragments[i]+temp_p5
		index += len(temp_p5)
		curr_len += len(temp_p5)

		## Add the 3' part
		p3_replace = [m.start() for m in re.finditer('[NMRWYSKHBVD]', p3_orig)]

		for replace_index in p3_replace:
			temp_p3 = temp_p3[:replace_index] + chromosome[index+replace_index] + temp_p3[replace_index+1:]

		fragments[i+1] = temp_p3+fragments[i+1]
		temp_frag.append(fragments[i])
		temp_frag.append(temp_start)
		temp_frag.append(curr_len-1)
		new_fragments.append(temp_frag)

	## append last fragment to list
	len_genome = start + len(chromosome[1])
	temp_frag = []
	temp_frag.append(fragments[-1])
	temp_frag.append(len_genome-len(fragments[-1]))
	temp_frag.append(len_genome-1)
	new_fragments.append(temp_frag)
	return new_fragments 

def size_selection(new_fragments, min_size, max_size):
	selected_fragments=[]
	for frag in new_fragments:
		if((len(frag[0]) > min_size) & (len(frag[0]) < max_size)):
			selected_fragments.append(frag)

	return selected_fragments

def perc_coverage(selected_fragments, len_chrom):
	#get % of covered chromosome
	ratio_frag=[]
	for frag in selected_fragments:
		ratio_frag.append(len(frag[0]))
		
	cov_sel=float(sum(ratio_frag))/len_chrom
	return(cov_sel)

def write_fasta(selected_fragments, enzyme, chrom_name):
#write fasta file 
	infile=open('/global/scratch/users/chandlersutherland/e19/re_digest/'+enzyme+'/'+chrom_name+'_fragment.fasta', 'w')
	i=0
	for frag in selected_fragments:
		#header includes arbitrary fragment #, Chr1, and then assigned index 
		header='> '+'frag'+ str(i) + ' '+chrom_name+' '+str(frag[1])+'-'+str(frag[2])+'\n'
		seq=frag[0]+'\n'
		infile.write(header)
		infile.write(seq)
		i += 1

	infile.close()

def write_bed(selected_fragments, enzyme, chrom_name, entry_type):
	#selected_fragments is a list of lists, where each entry is a [DNA_sequence, start index, stop index]
	#entry type is fragment of sequencing (150 PE)
	#write bed file 
	if chrom_name == 'ChrC':
		chrom_num = 'chloroplast'
	elif chrom_name == 'ChrM':
		chrom_num = 'mitochondria'
	else:
		chrom_num=chrom_name.strip('Chr')
	infile=open('/global/scratch/users/chandlersutherland/e19/re_digest/'+enzyme+'/'+chrom_name+'_'+entry_type+'.bed', 'w')
	#header='chrom'+'\t'+'chromStart'+'\t'+'chromEnd'+'\t'+'frag'+'\n'
	#infile.write(header)
	i=0
	for frag in selected_fragments:
		seq_name='frag'+str(i)
		left=str(frag[1])
		right=str(frag[2]+1) #bed coordinates are non-inclusive 
		line=chrom_num+"\t"+left+"\t"+right+'\t'+seq_name+"\n"
		infile.write(line)
		i += 1

	infile.close()


#simulate 150bp PE reads
def write_illumina_fasta_forward(selected_fragments, enzyme, chrom_name):
#write fasta file 
	infile=open('/global/scratch/users/chandlersutherland/e19/re_digest/'+enzyme+'/'+chrom_name+'_seq_1.fasta', 'w')
	i=0
	forward_read=[]
	for frag in selected_fragments:
		#header includes arbitrary fragment #, Chr1, and then assigned index 
		header='> '+'frag'+ str(i) + ' '+chrom_name+'\n'
		seq=frag[0][0:150]+'\n'
		start=frag[1]
		stop=frag[1]+149
		infile.write(header)
		infile.write(seq)
		forward_read.append([seq.strip('\n'), start, stop])
		i += 1
	print('finished forward')
	infile.close()

	return forward_read

def write_illumina_fasta_reverse(selected_fragments, enzyme, chrom_name):
	infile=open('/global/scratch/users/chandlersutherland/e19/re_digest/'+enzyme+'/'+chrom_name+'_seq_2.fasta', 'w')
	
	table = str.maketrans({
		"A": "T",
		"T": "A",
		"G": "C",
		"C": "G"})

	i=0
	reverse_read=[]
	for frag in selected_fragments:
		#header includes arbitrary fragment #, Chr1, and then assigned index 
		header='> '+'frag'+ str(i) + ' '+chrom_name+'\n'
		seq=frag[0][-150:].translate(table)[::-1]+'\n'  #rev comp 
		stop=frag[2]
		start=frag[2]-149
		infile.write(header)
		infile.write(seq)
		reverse_read.append([seq.strip('\n'), start, stop])
		i += 1
	infile.close()
	print('finished reverse')
	return reverse_read 


#cov=perc_coverage(selected_fragments)

#enzyme='AluI'
#cov_list=[]
#frag_to_write = [] 
#for chrom in genome:
#	chrom_name=chrom[0].strip('>')
#	temp_frag = []
#	temp_select_frag = []
#	print('starting '+chrom_name)
#	temp_frag=digest(enzyme, chrom)
#	temp_select_frag=size_selection(temp_frag, 250, 500)
#	len_chrom=len(chrom[1])
#	cov=perc_coverage(temp_select_frag, len_chrom)
#	cov_output=[enzyme, chrom_name, cov]
#	cov_list.append(cov_output)
#	print(cov_output)
#	frag_to_write.append([chrom[0], temp_select_frag])

#for chrom in frag_to_write:
#	chrom_name=chrom[0].strip('>')
#	write_fasta(chrom[1], enzyme, chrom_name)
#	write_bed(chrom[1], enzyme, chrom_name)


# output coverage file
# columns: enzyme, chromosome, coverage 
#output_coverage='/global/scratch/users/chandlersutherland/e19/re_digest/'+enzyme+'/coverage.csv'
#cov_df=pd.DataFrame(cov_list)
#cov_df.columns = ['enzyme', 'chrom', 'cov']
#cov_df.to_csv(output_coverage, sep='\t')

def wrapper(enzyme):
	os.makedirs('/global/scratch/users/chandlersutherland/e19/re_digest/'+enzyme, exist_ok=True)
	cov_list=[]
	frag_to_write = [] 
	for chrom in genome:
		chrom_name=chrom[0].strip('>')
		temp_frag = []
		temp_select_frag = []
		print('starting '+chrom_name+ ' '+enzyme)
		temp_frag=digest(enzyme, chrom)
		temp_select_frag=size_selection(temp_frag, 250, 500)
		len_chrom=len(chrom[1])
		cov=perc_coverage(temp_select_frag, len_chrom)
		
		cov_output=[enzyme, chrom_name, cov]
		cov_list.append(cov_output)
		frag_to_write.append([chrom[0], temp_select_frag])

	print('Fragments identified. Writing fasta and bed file')
	for chrom in frag_to_write:
		chrom_name=chrom[0].strip('>')
		write_fasta(chrom[1], enzyme, chrom_name)
		write_bed(chrom[1], enzyme, chrom_name, 'fragment')
		forward=write_illumina_fasta_forward(chrom[1], enzyme, chrom_name)
		reverse=write_illumina_fasta_reverse(chrom[1], enzyme, chrom_name)
		write_bed(forward, enzyme, chrom_name, 'seq_forward')
		write_bed(reverse, enzyme, chrom_name, 'seq_reverse')
	
	# output coverage file
	# columns: enzyme, chromosome, coverage 
	output_coverage='/global/scratch/users/chandlersutherland/e19/re_digest/'+enzyme+'/coverage.csv'
	cov_df=pd.DataFrame(cov_list)
	cov_df.columns = ['enzyme', 'chrom', 'cov']
	cov_df.to_csv(output_coverage, sep='\t')
	#return(frag_to_write) #just in case needed again


for enzyme in parsed['db']:
	wrapper(enzyme)

#glob.glob('')