import random

def datasplit_by_chrom(r_filename,d1_filename,d2_filename,valid_chrom_list,bal_motif_sitename_dict):
	"""
	chrom   genepos kmer    kmer_contig_readindex_tranpos   N2_mean N2_std  N2_lengthN2_quality       N1_mean N1_std  N1_length       N1_quality      P0_mean P0_std  P0_length P0_quality      P1_mean P1_std  P1_length       P1_quality      P2_mean P2_std    P2_length       P2_quality
	chrX    129906214       TGACT   TGACT_ENST00000394422_2937577_31        85.1    2.1205208333333334        0.0067643749999999996   15      94.6    5.596   0.00267 7131.8    11.697000000000001      0.00433 26      101.6   4.875684210526315       0.008456842105263157      31      86.6    5.685   0.010329999999999999    28
	chrX    129906214       TGACT   TGACT_ENST00000394422_2937578_31        85.4    1.11      0.00367 13      87.3    2.793096774193548       0.0067196774193548385   7123.9    8.542   0.006670000000000001    17      96.7    3.501957575757576       0.011879636363636365      25      87.3    4.66147619047619        0.003508095238095238      22
	"""
	f = open(r_filename,'r')
	d1 = open(d1_filename,'a')
	d2 = open(d2_filename,'a')
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0]=="chrom" or list_x[0]=="kmer":
			d1.write(str_x+"\n")
			d2.write(str_x+"\n")

			motif_index = list_x.index("kmer")
			readname_index = list_x.index("kmer_contig_readindex_tranpos")
			continue
		motif = list_x[motif_index]
		sitename = list_x[0] + ";" + list_x[1]
		readname = list_x[readname_index]
		R = random.random()
		#########################
		if bal_motif_sitename_dict == None:
			if R < 0.1:
				d1.write(str_x+"\n")
			else:
				d2.write(str_x+"\n")
		else:
			try:
				bal_readname_dict = bal_motif_sitename_dict[motif]
				try:
					v = bal_readname_dict[readname]
					if R < 0.1:
						d1.write(str_x+"\n")
					else:
						d2.write(str_x+"\n")
				except:
					pass
			except:
				pass

def file_split_by_motif(filename):
	f = open(filename,'r')
	outdict = {}
	outlist = []
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0]=="chrom" or list_x[0] == "kmer":
			continue
		motif = list_x[2]
		sitename = list_x[0] + ";" + list_x[1]
		readname = list_x[3]
		try:
			outdict[motif].append(readname)
		except:
			outdict[motif] = [readname]
			outlist.append(motif)
	return outdict,outlist

def motif_balance(pos_motif_to_sitename_dict,neg_motif_to_sitename_dict,motif_list,multi_number):
	for motif in motif_list:
		pos_sitename_list = pos_motif_to_sitename_dict[motif]
		neg_sitename_list = neg_motif_to_sitename_dict[motif]
		random.seed(1234)
		random.shuffle(neg_sitename_list)
		###########
		balance_number = multi_number * len(pos_sitename_list)
		###########
		tmp_v_list = [1 for x in neg_sitename_list[:balance_number]]
		tmp_k_list = neg_sitename_list[:balance_number]
		neg_motif_to_sitename_dict[motif] = dict(zip(tmp_k_list,tmp_v_list))
	return neg_motif_to_sitename_dict,pos_motif_to_sitename_dict


def main():
	apobe_neg_filename = "../data2/apobec.training.negative.genepos.sorted.baseflank.tsv"
	apobe_pos_filename = "../data2/apobec.training.positive.merge.tsv"

	neg_motif_to_sitename_dict,neg_motif_list = file_split_by_motif(filename=apobe_neg_filename)
	pos_motif_to_sitename_dict,pos_motif_list = file_split_by_motif(filename=apobe_pos_filename)

	print(neg_motif_list)
	print(pos_motif_list)

	bal_neg_motif_to_sitename_dict,bal_pos_motif_to_sitename_dict = motif_balance(
		pos_motif_to_sitename_dict=pos_motif_to_sitename_dict,
		neg_motif_to_sitename_dict=neg_motif_to_sitename_dict,
		motif_list=pos_motif_list,
		multi_number=30
		)
	apobe_neg_filename_train = "../data2/train_apobec.training.negative.genepos.sorted.tsv"
	apobe_pos_filename_train = "../data2/train_apobec.training.positive.genepos.sorted.tsv"
	##########################
	apobe_neg_filename_valid = "../data2/valid_apobec.training.negative.genepos.sorted.tsv"
	apobe_pos_filename_valid = "../data2/valid_apobec.training.positive.genepos.sorted.tsv"
	##########################
	datasplit_by_chrom(r_filename=apobe_neg_filename,d1_filename=apobe_neg_filename_valid,d2_filename=apobe_neg_filename_train,valid_chrom_list=["chr20","chr21","chr22","chrX"],bal_motif_sitename_dict=None)
	datasplit_by_chrom(r_filename=apobe_pos_filename,d1_filename=apobe_pos_filename_valid,d2_filename=apobe_pos_filename_train,valid_chrom_list=["chr20","chr21","chr22","chrX"],bal_motif_sitename_dict=None)

if __name__ == "__main__":
	main()