import random

def positive_data_motif_dict_make(pos_filename):
	"""
	kmer    kmer_contig_readindex_tranpos   P1_mean P1_std  P1_length       P0_mean P0_std  P0_length       N1_mean    N1_std  N1_length       baseflank
	GGACC   GGACC_ENST00000394422_2937577_2264      79.2    2.3865462184873953      0.016093361344537815    123.9      3.3989999999999996      0.00533 121.6   2.7858185053380775      0.007602491103202847    CGGACCT
	GGACC   GGACC_ENST00000394422_2937578_2264      74.0    1.527   0.008   108.8   10.017999999999999      0.002      119.1   8.494   0.01733 CGGACCT
	"""
	f = open(pos_filename,'r')
	outdict = {}
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0]=="chrom" or list_x[0]=="kmer":
			motif_index = list_x.index("kmer")
			baseflank_index = list_x.index("baseflank")
			continue
		motif = list_x[motif_index]
		baseflank = list_x[baseflank_index]
		try:
			outdict[baseflank] +=1
		except:
			outdict[baseflank] = 1
	return outdict

def negative_data_motif_dict_make(neg_filename):
	"""
	kmer    kmer_contig_readindex_tranpos   P1_mean P1_std  P1_length       P0_mean P0_std  P0_length       N1_mean    N1_std  N1_length       baseflank
	GGACC   GGACC_ENST00000394422_2937577_2264      79.2    2.3865462184873953      0.016093361344537815    123.9      3.3989999999999996      0.00533 121.6   2.7858185053380775      0.007602491103202847    CGGACCT
	GGACC   GGACC_ENST00000394422_2937578_2264      74.0    1.527   0.008   108.8   10.017999999999999      0.002      119.1   8.494   0.01733 CGGACCT
	"""
	f = open(neg_filename,'r')
	neg_motif2readname_dict = {}
	out_readname2value_dict = {}
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0]=="chrom" or list_x[0]=="kmer":
			motif_index = list_x.index("kmer")
			baseflank_index = list_x.index("baseflank")
			readname_index = list_x.index("kmer_contig_readindex_tranpos")
			continue
		baseflank = list_x[baseflank_index]
		readname = list_x[readname_index]
		out_readname2value_dict[readname] = list_x
		######################
		######################
		try:
			neg_motif2readname_dict[baseflank].append(readname)
		except:
			neg_motif2readname_dict[baseflank] = [readname]
	return neg_motif2readname_dict,out_readname2value_dict

def motif_balance(pos_motif2count_dict,neg_motif2readname_dict,min_negative_sample_number=10000,negative_over_sample_number=25):
	outdict = {}
	for motif in neg_motif2readname_dict.keys():
		try:
			pos_count = pos_motif2count_dict[motif] * negative_over_sample_number
			if pos_count >= min_negative_sample_number:
				pass
			else:
				pos_count = min_negative_sample_number
		except:
			pos_count = min_negative_sample_number
		###################################################
		neg_readname_list = neg_motif2readname_dict[motif]
		###################################################
		if len(neg_readname_list) >= pos_count:
			random.shuffle(neg_readname_list)
			random.shuffle(neg_readname_list)
			outdict[motif] = neg_readname_list[:pos_count]
		else:
			multi_number = (pos_count//len(neg_readname_list))+1
			multi_neg_readname_list = multi_number*neg_readname_list
			random.shuffle(multi_neg_readname_list)
			random.shuffle(multi_neg_readname_list)
			outdict[motif] = multi_neg_readname_list[:pos_count]
	return outdict

def balance_neg_result_write(filename,balanced_neg_motif2readname_dict,neg_readname2value_dict):
	d = open(filename,'a')
	###############################
	str2write = "\t".join(['kmer','kmer_contig_readindex_tranpos','P1_mean','P1_std','P1_length','P0_mean','P0_std','P0_length','N1_mean','N1_std','N1_length','baseflank']) + "\n"
	d.write(str2write)
	###############################
	repcount_dict = {}
	for motif in balanced_neg_motif2readname_dict.keys():
		tmp_neg_readname_list = balanced_neg_motif2readname_dict[motif]
		for readname_i in tmp_neg_readname_list:
			list_x = neg_readname2value_dict[readname_i]
			try:
				repcount_dict[readname_i] += 1
				readname_rep = list_x[1]+";"+str(repcount_dict[readname_i])
			except:
				repcount_dict[readname_i] = 0
				readname_rep = list_x[1]
			tmplist = list_x.copy()
			tmplist[1] = readname_rep
			str2write = "\t".join(tmplist) + "\n"
			d.write(str2write)

def motif_balance_neg_data_make(apobe_neg_filename,apobe_pos_filename,balan_neg_filename,min_negative_sample_number=1000,negative_over_sample_number=25):
	pos_motif2count_dict = positive_data_motif_dict_make(pos_filename=apobe_pos_filename)
	neg_motif2readname_dict,neg_readname2value_dict = negative_data_motif_dict_make(neg_filename=apobe_neg_filename)
	balanced_neg_motif2readname_dict = motif_balance(
		pos_motif2count_dict=pos_motif2count_dict,
		neg_motif2readname_dict=neg_motif2readname_dict,
		min_negative_sample_number=min_negative_sample_number,
		negative_over_sample_number=negative_over_sample_number
		)
	balance_neg_result_write(
		filename=balan_neg_filename,
		balanced_neg_motif2readname_dict=balanced_neg_motif2readname_dict,
		neg_readname2value_dict=neg_readname2value_dict,
		)
