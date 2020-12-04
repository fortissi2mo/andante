#python


import sys
import re
import glob


input_file = '/home/yohan1633/HLAexp/hla.dat'

output_file = '/home/yohan1633/HLAexp/hla.pos.idx.txt'


INPUT = open(input_file , 'r')


#OUTPUT
prot_fasta
gen_fasta
exon_fasta

digit4_fasta
digit6_fasta
digit8_fasta


output_header = '#HLA_ALLELE>trs_info:pos/trs_info2:pos2/'
OUTPUT.write(output_header + '\n')
total_line = INPUT.read()


allele_list = total_line.split('//')
hla_list=[]
for allele in allele_list : 
	seq_info = allele.split('SQ')[1]
	header=allele.split('\n') 
	
	w_pos_line = '' # UTR:1..300/exon:301..373
	
	for info in header : 

		if info[0:2] == 'FT' :
			i = info.split(' ')[-1]

	for idx, line in enumerate(seq_info.split('\n')) : 
    if idx == 0 : continue
		seq_group = line.lstrip(' ').split(' ')
		seq_line = seq_group[0] + seq_group[1] + seq_group[2] + seq_group[3] + seq_group[4] + seq_group[5]
		seq_line = seq_line.upper()



			hla=info.lstrip('DE  ').split(',')[0] # HLA-A*01:01:01:01
		if info[0:2] != 'FT' : continue
		trs_idx=info.split(' ')[3]
		if trs_idx == 'UTR' or trs_idx == 'exon' or trs_idx == 'intron' : 
			pos = info.split(' ')[-1]
			w_pos_line = w_pos_line + trs_idx + ':' + pos + '/'

	w_line = hla + '>' + w_pos_line
	if hla not in hla_list : 
		if len(w_pos_line) != 0 : 
			OUTPUT.write(w_line + '\n') 
			hla_list.append(hla)



OUTPUT.close()
INPUT.close()
