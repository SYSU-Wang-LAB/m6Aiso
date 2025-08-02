
import gzip
import numpy as np

def signal_data_read(filename,using_start,using_end):
	"""
	chrom   genome_position genome_motif    strand  transcript_id   transcript_position     readname        up_and_down_sequecne        mean    stdv    length  skew    kurtosis        median  quantile_25     quantile_75     max     min base_calling_sequecne
	chr4    672443  CTGTGTTGTGTGA   -       ENST00000304312.5       290     c42ece04-b898-41c7-b12a-f9ddc604f473    CTGTGTTGTGTGA       0.1099|0.6998|-0.4487|-0.2136|-0.01117  0.1776|0.4867|0.1276|0.5915|0.2472      0.009|0.0075|0.007|0.0095|0.003     -0.2128|-0.04567|0.2801|1.113|-1.15     1.786|-1.11|2.104|1.038|0.4595  0.09483|0.6475|-0.4666|-0.2999|0.04219      0.004909|0.3141|-0.5105|-0.7035|-0.09378        0.1913|1.143|-0.403|0.09263|0.1891      0.5334|1.533|-0.1069|1.569|0.244    -0.4578|-0.2122|-0.8087|-0.914|-0.6157  CTGTGTTGTGTGA
	chr4    672444  CCTGTGTTGTGTG   -       ENST00000304312.5       289     c42ece04-b898-41c7-b12a-f9ddc604f473    CCTGTGTTGTGTG       -0.09507|0.1099|0.6998|-0.4487|-0.2136  0.496|0.1776|0.4867|0.1276|0.5915       0.00425|0.009|0.0075|0.007|0.0095   2.034|-0.2128|-0.04567|0.2801|1.113     4.226|1.786|-1.11|2.104|1.038   -0.2297|0.09483|0.6475|-0.4666|-0.2999      -0.4227|0.004909|0.3141|-0.5105|-0.7035 0.05096|0.1913|1.143|-0.403|0.09263     1.542|0.5334|1.533|-0.1069|1.569    -0.5631|-0.4578|-0.2122|-0.8087|-0.914  CCTGTGTTGTGTG
	"""
	f = gzip.open(filename,'rb')
	name2sequence_dict = {}
	name2signal_dict = {}
	sitename_list = []
	for str_x in f:
		str_x = str_x.decode()
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "chrom":
			index_dict = dict(zip(list_x,list(range(len(list_x)))))
			continue
		###############################
		chrom = list_x[index_dict['chrom']]
		genomepos = list_x[index_dict['genome_position']]
		strand = list_x[index_dict['strand']]
		readname = list_x[index_dict['readname']]

		###############################
		genome_motif = list_x[index_dict['genome_motif']]
		up_and_down_sequecne = list_x[index_dict['up_and_down_sequecne']]
		base_calling_sequecne = list_x[index_dict['base_calling_sequecne']]

		sitename = ";".join([chrom,genomepos,strand,readname,genome_motif])

		#print(sitename)

		mean = list_x[index_dict['mean']]
		stdv = list_x[index_dict['stdv']]
		length = list_x[index_dict['length']]
		skew = list_x[index_dict['skew']]
		kurtosis = list_x[index_dict['kurtosis']]
		median = list_x[index_dict['median']]
		quantile_25 = list_x[index_dict['quantile_25']]
		quantile_75 = list_x[index_dict['quantile_75']]
		value_max = list_x[index_dict['max']]
		value_min = list_x[index_dict['min']]
		###############################
		#print(base_calling_sequecne)
		###############################
		genome_motif = sequence_abstract(sequence=genome_motif,using_start=using_start,using_end=using_end)
		up_and_down_sequecne = sequence_abstract(sequence=up_and_down_sequecne,using_start=using_start,using_end=using_end)
		base_calling_sequecne = sequence_abstract(sequence=base_calling_sequecne,using_start=using_start,using_end=using_end)
		genome_motif_coding = oneHot_encoding_sequence(sequence=genome_motif)
		up_and_down_sequecne_coding = oneHot_encoding_sequence(sequence=up_and_down_sequecne)
		base_calling_sequecne_coding = oneHot_encoding_sequence(sequence=base_calling_sequecne)
		###############################
		mean = signal_abstract(signal=mean,using_start=using_start,using_end=using_end)
		stdv = signal_abstract(signal=stdv,using_start=using_start,using_end=using_end)
		length = signal_abstract(signal=length,using_start=using_start,using_end=using_end)
		skew = signal_abstract(signal=skew,using_start=using_start,using_end=using_end)
		kurtosis = signal_abstract(signal=kurtosis,using_start=using_start,using_end=using_end)
		median = signal_abstract(signal=median,using_start=using_start,using_end=using_end)
		quantile_25 = signal_abstract(signal=quantile_25,using_start=using_start,using_end=using_end)
		quantile_75 = signal_abstract(signal=quantile_75,using_start=using_start,using_end=using_end)
		value_max = signal_abstract(signal=value_max,using_start=using_start,using_end=using_end)
		value_min = signal_abstract(signal=value_min,using_start=using_start,using_end=using_end)
		###############################
		###############################
		signal_matrix = signal_covert_to_matrix(
			mean_list=mean,
			stdv_list=stdv,
			dtime_list=length,
			skew_list=skew,
			kurtosis_list=kurtosis,
			median_list=median,
			quantile_25_list=quantile_25,
			quantile_75_list=quantile_75,
			value_max_list=value_max,
			value_min_list=value_min
			)
		sequence_matrix = sequence_convert_to_matrix(
			coding_sequence_1=genome_motif_coding,
			coding_sequence_2=base_calling_sequecne_coding
			)
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


def signal_covert_to_matrix(mean_list,stdv_list,dtime_list,skew_list,kurtosis_list,median_list,quantile_25_list,quantile_75_list,value_max_list,value_min_list):
	outlist = []
	for i in range(len(mean_list)):
		x = [mean_list[i],stdv_list[i],dtime_list[i],skew_list[i],kurtosis_list[i],value_min_list[i],quantile_25_list[i],median_list[i],quantile_75_list[i],value_max_list[i]]
		outlist.append(x)
	return outlist

def sequence_convert_to_matrix(coding_sequence_1,coding_sequence_2):
	outlist = []
	for i in range(len(coding_sequence_1)):
		merge_list_i = coding_sequence_1[i] + coding_sequence_2[i]
		outlist.append(merge_list_i)
	return outlist

def signal_abstract(signal,using_start,using_end):
	signal_list = signal.split("|")
	signal_list = [float(x) for x in signal_list]
	####################################################
	v = [float(x) for x in signal_list[using_start:using_end]]
	return v

def sequence_abstract(sequence,using_start,using_end):
	#print(sequence)
	#sequence = sequence[3:-3]
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
		elif base_i == "*":
			outlist.append([0.25,0.25,0.25,0.25,0.25])
		else:
			raise Exception("oneHot_encoding_sequence error")
	return outlist

if __name__=="__main__":
	import sys
	signal_data_read(filename=sys.argv[1],using_start=1,using_end=4)