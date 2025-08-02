
import gzip
import numpy as np
import random

def signal_fileread(filename,using_start,using_end):
	"""
	kmer    kmer_contig_readindex_tranpos   P1_mean P1_std  P1_length       P0_mean P0_std  P0_length       N1_mean N1_std      N1_length       baseflank
	AGGGACACT       AGGGACACT_ENST00000373812.8_9505cb45-2be2-4951-9c4b-65f4a4fc284f_1759   -0.478046       0.0762121  0.00325  1.23739 0.236583        0.002   1.20801 0.145114        0.00575 TAGGGACACTC
	AAAGACCGT       AAAGACCGT_ENST00000371281.4_749d9694-a44a-4b07-bd4d-44dc7092c122_1763   0.165908        0.264325   0.0075   2.02827 0.675255        0.0095  1.96534 0.149658        0.00125 TAAAGACCGTA
	sitename = list_x[index_dict['kmer_contig_readindex_tranpos']]
	flank_sequecne = list_x[index_dict['baseflank']]
	kmer_sequecne = list_x[index_dict['kmer']]
	####################################################
	chrom   genome_position genome_motif    strand  transcript_id   transcript_position     readname        up_and_down_sequecne        mean    stdv    length  skew    kurtosis        median  quantile_25     quantile_75     max     minbase_calling_sequecne
	DRACH_oligo_0   82      CTGCAGACTGATA   +       DRACH_oligo_0   78      7bc298ab-face-4966-945d-b73b04c5bc31    CTGCAGACTGATA       0.6269|0.6269|0.6269|-0.3477|0.2284     0.1104|0.1104|0.1104|0.09477|0.1005     0.0085|0.0085|0.0085|0.004|0.00425  -0.8512|-0.8512|-0.8512|-0.5198|-0.6454 1.196|1.196|1.196|-1.279|-0.5155        0.6155|0.6155|0.6155|-0.3092|0.2537 0.5552|0.5552|0.5552|-0.4583|0.18       0.7076|0.7076|0.7076|-0.269|0.3073      0.8098|0.8098|0.8098|-0.2154|0.381  0.2738|0.2738|0.2738|-0.5102|0.005744   CTGCAGACTGATA
	DRACH_oligo_0   64      GACATAACTGTGC   +       DRACH_oligo_0   60      7bc298ab-face-4966-945d-b73b04c5bc31    GACATAACTGTGC       -0.03254|0.1217|0.3918|-0.2096|-0.1552  0.0424|0.05065|0.1942|0.04975|0.1052    0.021|0.00325|0.05125|0.0035|0.0115 -0.04476|-0.457|0.2986|0.7642|-0.3807   0.7279|0.03287|-0.8621|-0.7703|-0.9533  -0.02776|0.1331|0.3877|-0.2321|-0.1249      -0.06126|0.08615|0.2001|-0.2489|-0.2589 -0.007657|0.1465|0.5552|-0.1785|-0.08136        0.09285|0.2068|0.9237|-0.1082|0.02585       -0.1618|0.005744|0.05265|-0.269|-0.3628 GACATAACTGTGC

	"""
	f = gzip.open(filename,'rb')
	motif2count_dict = {}
	motif2sitename_dict = {}
	motiflist = []
	for str_x in f:
		str_x = str_x.decode()
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "chrom":
			index_dict = dict(zip(list_x,list(range(len(list_x)))))
			continue
		###############################
		chrom = list_x[index_dict['chrom']]
		genome_position = list_x[index_dict['genome_position']]
		strand = list_x[index_dict['strand']]
		readname = list_x[index_dict['readname']]
		###############################
		flank_sequecne = list_x[index_dict['genome_motif']]
		up_and_down_sequecne = list_x[index_dict['up_and_down_sequecne']]
		sitename = ";".join([chrom,genome_position,strand,readname,flank_sequecne])
		
		###############################
		###############################
		try:
			motif2count_dict[flank_sequecne] += 1
			motif2sitename_dict[flank_sequecne].append(sitename)
		except:
			motif2count_dict[flank_sequecne] = 1
			motiflist.append(flank_sequecne)
			motif2sitename_dict[flank_sequecne] = [sitename]
	return motif2count_dict,motif2sitename_dict,motiflist

def motif_list_merge(pos_motif2count_dict,pos_motiflist,neg_motif2count_dict,neg_motiflist,min_count):
	merged_motif_set = set(pos_motiflist + neg_motiflist)
	merged_motif_list = list(merged_motif_set)
	out_motif_list = []
	for motif in merged_motif_list:
		try:
			pos_count = pos_motif2count_dict[motif]
		except:
			pos_count = 0
		try:
			neg_count = neg_motif2count_dict[motif]
		except:
			neg_count = 0
		if pos_count >= min_count and neg_count >= min_count:
			out_motif_list.append(motif)
		else:
			pass
	return out_motif_list

def balance_motif_based_on_count(pos_motif2count_dict,neg_motif2count_dict,pos_motif2sitename_dict,neg_motif2sitename_dict,merged_motif_list,label_filename):
	label_file = gzip.open(label_filename,'ab')
	for motif in merged_motif_list:
		pos_count = pos_motif2count_dict[motif]
		neg_count = neg_motif2count_dict[motif]
		print(motif,pos_count,neg_count)
		pos_sitename_list = pos_motif2sitename_dict[motif]
		neg_sitename_list = neg_motif2sitename_dict[motif]
		#############################################
		if neg_count >= pos_count:
			random.shuffle(pos_sitename_list)
			random.shuffle(neg_sitename_list)
			tmp_pos_sitename_list = pos_sitename_list
			tmp_neg_sitename_list = neg_sitename_list[:pos_count]
			resultwrite(
				label_file=label_file,
				pos_sitename_list=tmp_pos_sitename_list,
				neg_sitename_list=tmp_neg_sitename_list,
				motif=motif
				)
		#############################################
		else:
			multivalue = int(pos_count/neg_count) + 1
			multi_neg_sitename_list = neg_sitename_list * multivalue
			random.shuffle(pos_sitename_list)
			random.shuffle(multi_neg_sitename_list)
			tmp_pos_sitename_list = pos_sitename_list
			tmp_neg_sitename_list = multi_neg_sitename_list[:pos_count]
			resultwrite(
				label_file=label_file,
				pos_sitename_list=tmp_pos_sitename_list,
				neg_sitename_list=tmp_neg_sitename_list,
				motif=motif
				)
		#############################################

def resultwrite(label_file,pos_sitename_list,neg_sitename_list,motif):
	for sitename in pos_sitename_list:
		str2write = "\t".join([sitename,'1.0','0.0']) + "\n"
		label_file.write(str2write.encode())
	for sitename in neg_sitename_list:
		str2write = "\t".join([sitename,'0.0','1.0']) + "\n"
		label_file.write(str2write.encode())
	print(motif,len(pos_sitename_list),len(neg_sitename_list))

def args_make():
	import argparse
	parser = argparse.ArgumentParser(description='label file make based on motif')
	parser.add_argument('--pos_filename', required=True, help="pos_filename")
	parser.add_argument('--neg_filename', required=True, help="neg_filename")
	parser.add_argument('--min_motif_count', required=True, help="min motif count")
	parser.add_argument('--label_filename', required=True, help="label_filename")
	args = parser.parse_args()
	return args


if __name__=="__main__":
	import sys
	args = args_make()
	pos_filename = args.pos_filename
	neg_filename = args.neg_filename
	min_motif_count = int(args.min_motif_count)
	label_filename = args.label_filename
	##########################
	pos_motif2count_dict,pos_motif2sitename_dict,pos_motiflist = signal_fileread(
		filename=pos_filename,
		using_start=0,
		using_end=5
		)
	neg_motif2count_dict,neg_motif2sitename_dict,neg_motiflist = signal_fileread(
		filename=neg_filename,
		using_start=0,
		using_end=5
		)
	merged_motif_list = motif_list_merge(
		pos_motif2count_dict=pos_motif2count_dict,
		pos_motiflist=pos_motiflist,
		neg_motif2count_dict=neg_motif2count_dict,
		neg_motiflist=neg_motiflist,
		min_count=min_motif_count
		)
	balance_motif_based_on_count(
		pos_motif2count_dict=pos_motif2count_dict,
		neg_motif2count_dict=neg_motif2count_dict,
		pos_motif2sitename_dict=pos_motif2sitename_dict,
		neg_motif2sitename_dict=neg_motif2sitename_dict,
		merged_motif_list=merged_motif_list,
		label_filename=label_filename
		)







	