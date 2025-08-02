import gzip
import numpy as np

def signal_data_read(filename,using_start,using_end):
	"""
	kmer    kmer_contig_readindex_tranpos   P1_mean P1_std  P1_length       P0_mean P0_std  P0_length       N1_mean N1_std      N1_length       baseflank
	AGGGACACT       AGGGACACT_ENST00000373812.8_9505cb45-2be2-4951-9c4b-65f4a4fc284f_1759   -0.478046       0.0762121  0.00325  1.23739 0.236583        0.002   1.20801 0.145114        0.00575 TAGGGACACTC
	AAAGACCGT       AAAGACCGT_ENST00000371281.4_749d9694-a44a-4b07-bd4d-44dc7092c122_1763   0.165908        0.264325   0.0075   2.02827 0.675255        0.0095  1.96534 0.149658        0.00125 TAAAGACCGTA
	"""
	f = gzip.open(filename,'rb')
	name2sequence_dict = {}
	name2signal_dict = {}
	sitename_list = []
	for str_x in f:
		str_x = str_x.decode()
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "kmer":
			index_dict = dict(zip(list_x,list(range(len(list_x)))))
			continue
		###############################
		sitename = list_x[index_dict['kmer_contig_readindex_tranpos']]
		###############################
		P1_mean = list_x[index_dict['P1_mean']]
		P0_mean = list_x[index_dict['P0_mean']]
		N1_mean = list_x[index_dict['N1_mean']]
		P1_std = list_x[index_dict['P1_std']]
		P0_std = list_x[index_dict['P0_std']]
		N1_std = list_x[index_dict['N1_std']]
		P1_length = list_x[index_dict['P1_length']]
		P0_length = list_x[index_dict['P0_length']]
		N1_length = list_x[index_dict['N1_length']]
		################################
		flank_sequecne = list_x[index_dict['baseflank']]
		################################
		mean = [P1_mean,P0_mean,N1_mean]
		stdv = [P1_std,P0_std,N1_std]
		length = [P1_length,P0_length,N1_length]
		################################
		flank_sequecne = sequence_abstract(sequence=flank_sequecne,using_start=using_start,using_end=using_end)
		flank_sequecne_coding = oneHot_encoding_sequence(sequence=flank_sequecne)
		###############################
		mean = signal_abstract(signal=mean,using_start=using_start,using_end=using_end)
		stdv = signal_abstract(signal=stdv,using_start=using_start,using_end=using_end)
		length = signal_abstract(signal=length,using_start=using_start,using_end=using_end)
		################################
		################################
		signal_matrix = signal_covert_to_matrix(
			mean_list=mean,
			stdv_list=stdv,
			dtime_list=length
			)
		#########################################
		sequence_matrix = flank_sequecne_coding
		#########################################
		name2sequence_dict[sitename] = sequence_matrix
		name2signal_dict[sitename] = signal_matrix
		sitename_list.append(sitename)
		
	return name2sequence_dict,name2signal_dict,sitename_list

def sitename_to_label_dictmake(filename):
	f = gzip.open(filename,'rb')
	sitename2label_dict = {}
	sitename_list = []
	for str_x in f:
		str_x = str_x.decode()
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		sitename = list_x[0]
		value = [float(x) for x in list_x[1:]]
		sitename2label_dict[sitename] = value
		sitename_list.append(sitename)
		#print(sitename,value)
	return sitename2label_dict,sitename_list

def sitename_id_make(sitename_list):
	out_id2sitename_dict = {}
	out_id_list = []
	for i in range(len(sitename_list)):
		out_id2sitename_dict[i] = sitename_list[i]
		out_id_list.append(i)
	return out_id2sitename_dict,out_id_list


def signal_covert_to_matrix(mean_list,stdv_list,dtime_list):
	outlist = []
	for i in range(len(mean_list)):
		x = [mean_list[i],stdv_list[i],dtime_list[i]]
		outlist.append(x)
	return outlist

def sequence_convert_to_matrix(coding_sequence_1,coding_sequence_2):
	outlist = []
	for i in range(len(coding_sequence_1)):
		merge_list_i = coding_sequence_1[i] + coding_sequence_2[i]
		outlist.append(merge_list_i)
	return outlist

def signal_abstract(signal,using_start,using_end):
	signal_list = signal
	signal_list = [float(x) for x in signal_list]
	####################################################
	v = [float(x) for x in signal_list[using_start:using_end]]
	return v

def sequence_abstract(sequence,using_start,using_end):
	#print(sequence)
	tmp_start = using_start
	tmp_end = len(sequence) - ((len(sequence) - 9 + 1) - using_end)
	out_sequence = sequence[tmp_start:tmp_end]
	#print(out_sequence)
	return out_sequence

def oneHot_encoding_sequence(sequence):
	sequence = sequence
	outlist = []
	for i in range(len(sequence)):
		base_i = sequence[i]
		if base_i == "A":
			outlist.append([1.0,0.0,0.0,0.0,0.0])
		elif base_i == "T":
			outlist.append([0.0,1.0,0.0,0.0,0.0])
		elif base_i == "C":
			outlist.append([0.0,0.0,1.0,0.0,0.0])
		elif base_i == "G":
			outlist.append([0.0,0.0,0.0,1.0,0.0])
		elif base_i == ".":
			outlist.append([0.0,0.0,0.0,0.0,1.0])
		elif base_i == "N":
			outlist.append([0.25,0.25,0.25,0.25,0.25])
		else:
			raise Exception("oneHot_encoding_sequence error")
	return outlist

if __name__=="__main__":
	import sys
	#signal_data_read(filename=sys.argv[1],using_start=0,using_end=3)
	sitename_to_label_dictmake(filename=sys.argv[1])