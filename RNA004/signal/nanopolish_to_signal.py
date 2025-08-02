import gzip
import random
import math
import numpy as np
import scipy.stats as st

def nanopolish_fileread(nanopolish_filename,signal_filename,chrom2read_dict,transcript_to_chrom_dict,transcript_to_exon_dict,transcript_to_strand_dict,sample_select_length,signal_extend_length):
	"""
	contig  position        reference_kmer  read_index      strand  event_index     event_level_mean        event_stdv event_length    model_kmer      model_mean      model_stdv      standardized_level      start_idx end_idx  samples
	ENST00000390447.3       180     CAGAA   0       t       436     107.09  6.521   0.00500 CAGAA   108.60  5.87       -0.23   44472   44487   125.118,110.627,98.1655,104.252,110.917,109.178,101.064,105.411,106.86,109.178,117.873,108.599,94.6877,105.701,123.379
	ENST00000390447.3       181     AGAAG   0       t       437     129.19  4.614   0.00600 AGAAG   123.66  5.56       0.87    44454   44472   136.421,129.465,134.972,128.016,125.988,125.408,122.8,125.408,129.465,135.551,137.58,117.003,132.943,127.147,132.943,132.653,130.335,132.653
	contig  position        reference_kmer  read_name       strand  event_index     event_level_mean        event_stdv event_length    model_kmer      model_mean      model_stdv      standardized_level      start_idx end_idx  samples
	ENST00000440196.3  1       GGGAA   f1584f7c-675d-4efe-82da-9ccc9fca7195    t       2122    117.67     1.553   0.00232 GGGAA   120.77  3.17    -0.82   30326   30333   117.254,115.209,120.759,119.299,117.692,116.815,116.669
	ENST00000440196.3  2       GGAAG   f1584f7c-675d-4efe-82da-9ccc9fca7195    t       2123    114.16     2.417   0.00631 GGAAG   115.76  5.56    -0.24   30307   30326   115.355,116.231,114.771,115.501,113.748,116.523,114.478,112.872,114.917,118.86,111.849,118.568,114.04,114.186,114.478,111.265,108.636,114.04,108.636
	"""
	f = gzip.open(nanopolish_filename,'rb')
	writefile = gzip.open(signal_filename,'ab')
	same_readname_signal_list = []
	readname_to_lenght_dict = {}
	for str_x in f:
		str_x = str_x.decode()
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		list_x = list_x[1:]
		if list_x[0] == "chrom":
			try:
				read_id_index = list_x.index("read_index")
				reference_motif_index = list_x.index("reference_kmer")
				model_motif_index = list_x.index("model_kmer")
			except:
				read_id_index = list_x.index("read_name")
				reference_motif_index = list_x.index("reference_kmer")
				model_motif_index = list_x.index("model_kmer")
			str2write = "\t".join(['chrom','genome_position','genome_motif','strand','transcript_id','transcript_position','readname','up_and_down_sequecne','mean','stdv','length','skew','kurtosis','median','quantile_25','quantile_75','max','min']) + "\n"
			writefile.write(str2write.encode())
			continue
		########################
		readname = list_x[read_id_index]
		reference_kmer = list_x[reference_motif_index]
		model_kmer = list_x[model_motif_index]
		########################
		if not reference_kmer == model_kmer:
			continue
		########################
		try:
			readname_to_lenght_dict[readname] += 1
			same_readname_signal_list.append(list_x)
		except:
			readname_to_lenght_dict[readname] = 1
			################
			if len(same_readname_signal_list) == 0:
				same_readname_signal_list = []
				same_readname_signal_list.append(list_x)
			else:
				#print("####################################")
				#print(len(same_readname_signal_list))
				#print(same_readname_signal_list[0])
				genome_annotated_signal_list = repilicate_signal_merge(
					same_readname_signal_list=same_readname_signal_list,
					chrom2read_dict=chrom2read_dict,
					transcript_to_chrom_dict=transcript_to_chrom_dict,
					transcript_to_exon_dict=transcript_to_exon_dict,
					transcript_to_strand_dict=transcript_to_strand_dict,
					select_value_length=sample_select_length,
					)
				if len(genome_annotated_signal_list) == 0:
					pass
				else:
					signal_merge_and_result_write(tmp_signal_list=genome_annotated_signal_list,length=signal_extend_length,writefile=writefile)
				same_readname_signal_list = []
				same_readname_signal_list.append(list_x)


def multi_value_make(input_signal_list):
	samples = input_signal_list[-1]
	samples_value_list = samples.split(",")
	samples_value_list = [float(x) for x in samples_value_list]
	v1 = st.skew(samples_value_list)
	v2 = st.kurtosis(samples_value_list)
	v3 = np.median(samples_value_list)
	v4 = np.quantile(samples_value_list,0.25)
	v5 = np.quantile(samples_value_list,0.75)
	v6 = max(samples_value_list)
	v7 = min(samples_value_list)
	out_value_list = input_signal_list[:-1] + [v1,v2,v3,v4,v5,v6,v7]
	return out_value_list

def repilicate_signal_merge(same_readname_signal_list,chrom2read_dict,transcript_to_chrom_dict,transcript_to_exon_dict,transcript_to_strand_dict,select_value_length=20):
	"""
	contig  position        reference_kmer  read_index      strand  event_index     event_level_mean        event_stdv event_length    model_kmer      model_mean      model_stdv      standardized_level      start_idx end_idx  samples
	ENST00000390447.3       193     AAGCG   0       t       465     106.57  2.488   0.00267 AAGCG   109.95  3.02       -0.99   44041   44049   109.758,108.599,106.57,109.468,102.223,105.411,110.627,103.672
	ENST00000390447.3       193     AAGCG   0       t       466     112.01  11.292  0.00867 AAGCG   109.95  3.02       0.60    44015   44041   94.6877,117.873,117.873,111.207,123.089,120.771,103.382,108.019,120.771,120.191,121.061,107.44,123.089,101.933,76.1396,75.8498,113.236,122.22,124.539,104.252,116.424,116.134,119.322,114.975,114.975,107.729
	"""
	position_dict = {}
	position_list = []
	for list_x in same_readname_signal_list:
		pos = list_x[1]
		try:
			position_dict[pos].append(list_x)
		except:
			position_dict[pos] = [list_x]
			position_list.append(pos)
	####################
	out_genome_annotated_signal_list = []
	####################
	for pos in position_list:
		same_pos_list = position_dict[pos]
		if len(same_pos_list) > 1:
			##################################
			samples_merge = ",".join([x[13] for x in same_pos_list])
			samples_merge_list = [float(x) for x in samples_merge.split(",")]
			event_mean_merge = np.mean(samples_merge_list)
			event_std_merge = np.std(samples_merge_list)
			event_length_list = [float(x[8]) for x in same_pos_list]
			event_length_merge = sum(event_length_list)
			##################################
			merge_contig = same_pos_list[0][0]
			merge_position = same_pos_list[0][1]
			merge_reference_kmer = same_pos_list[0][2]
			merge_read_index = same_pos_list[0][3]
			merge_strand = same_pos_list[0][4]
			merge_event_index = same_pos_list[0][5]
			##################################
			merge_event_level_mean = str(event_mean_merge)
			merge_event_level_std = str(event_std_merge)
			merge_event_level_length = str(event_length_merge)
			##################################
			merge_model_kmer =  same_pos_list[0][9]
			merge_model_mean = same_pos_list[0][10]
			merge_model_stdv = same_pos_list[0][11]
			merge_standardized_level = same_pos_list[0][12]
			##################################
			#merge_start_idx = str(min([int(x[13]) for x in same_pos_list]))
			#merge_end_idx = str(max([int(x[14]) for x in same_pos_list]))
			merge_samples = samples_merge
			###################################
			tmp_merged_list = [
			merge_contig,merge_position,merge_reference_kmer,merge_read_index,
			merge_strand,merge_event_index,merge_event_level_mean,merge_event_level_std,
			merge_event_level_length,merge_model_kmer,merge_model_mean,merge_model_stdv,
			merge_standardized_level,merge_samples
			]
		else:
			tmp_merged_list = same_pos_list[0]

		##############################################################
		sample_selected_signal_list = multi_value_make(input_signal_list=tmp_merged_list)
		genome_annotation_list = genome_position_and_motif_annotation(
			input_signal_list=sample_selected_signal_list,
			transcript_to_chrom_dict=transcript_to_chrom_dict,
			transcript_to_exon_dict=transcript_to_exon_dict,
			transcript_to_strand_dict=transcript_to_strand_dict,
			chrom2read_dict=chrom2read_dict
			)
		if genome_annotation_list == []:
			pass
		else:
			out_genome_annotated_signal_list.append(genome_annotation_list + sample_selected_signal_list)
	return out_genome_annotated_signal_list

def gtf_fileread(filename):
	"""
	ENST00000456328.2       chr1    +       11868   14409   14409   14409   3       11868,12612,13220,      12227,12721,14409,
	ENST00000450305.2       chr1    +       12009   13670   13670   13670   6       12009,12178,12612,12974,13220,13452,       12057,12227,12697,13052,13374,13670,
	ENST00000488147.1       chr1    -       14403   29570   29570   29570   11      14403,15004,15795,16606,16857,17232,17605,17914,18267,24737,29533, 14501,15038,15947,16765,17055,17368,17742,18061,18366,24891,29570,
	"""
	f = open(filename,'r')
	out_transcript_to_exon_dict = {}
	out_transcript_to_strand_dict = {}
	out_transcript_to_chrom_dict = {}
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		#############################
		transid = list_x[0].split(".")[0]
		chrom = list_x[1]
		strand = list_x[2]
		exon_start_list = [int(x) for x in list_x[8].split(",")[:-1]]
		exon_end_list = [int(x) for x in list_x[9].split(",")[:-1]]
		exon_list = []
		for i in range(len(exon_start_list)):
			exon_list.append([exon_start_list[i],exon_end_list[i]])
		#############################
		out_transcript_to_exon_dict[transid] = exon_list
		out_transcript_to_strand_dict[transid] = strand
		out_transcript_to_chrom_dict[transid] = chrom
		#############################
	return out_transcript_to_exon_dict,out_transcript_to_strand_dict,out_transcript_to_chrom_dict

def fasta_fileread(filename):
	f = open(filename,'r')
	out_read_dict = {}
	out_chrom_list = []
	for str_x in f:
		str_x = str_x.strip("\n")
		############
		if str_x[0] == ">":
			tmp_chrom = str_x[1:].split(" ")[0]
			out_read_dict[tmp_chrom] = []
			out_chrom_list.append(tmp_chrom)
		else:
			out_read_dict[tmp_chrom].append(str_x)
		############
	out_chrom_to_read_dict = {}
	for chrom in out_chrom_list:
		tmpread = "".join(out_read_dict[chrom])
		out_chrom_to_read_dict[chrom] = tmpread
	return out_chrom_to_read_dict

def genomePosition_to_transcriptPosition(genome_position,exon_list):
	transcript_position = 0
	for exon in exon_list:
		exon_start = exon[0]
		exon_end = exon[1]
		if genome_position >= exon_end and genome_position > exon_start:
			transcript_position += (exon_end - exon_start)
		elif genome_position < exon_end and genome_position >= exon_start:
			transcript_position += (genome_position - exon_start)
		elif genome_position < exon_end and genome_position < exon_start:
			transcript_position += 0
		else:
			raise Exception("genomePosition_to_transcriptPosition error")
	return transcript_position

def exon_length_calculate(exon_list):
	length = 0
	for i in range(len(exon_list)):
		exon_start_i = exon_list[i][0]
		exon_end_i = exon_list[i][1]
		length += (exon_end_i - exon_start_i)
	return length

def transcriptPosition_to_genomePosition(transcript_position,exon_list):
	genome_position = exon_list[0][0]
	tmp_transcript_position = transcript_position
	for i in range(len(exon_list)):
		exon_start = exon_list[i][0]
		exon_end = exon_list[i][1]
		tmp_genome_position = genome_position + tmp_transcript_position
		if tmp_genome_position > exon_end and tmp_genome_position > exon_start:
			genome_position = exon_list[i+1][0]
			tmp_transcript_position = tmp_transcript_position - (exon_end - exon_start)
		elif tmp_genome_position < exon_end and tmp_genome_position > exon_start:
			genome_position = tmp_genome_position
			tmp_transcript_position = tmp_transcript_position - (tmp_genome_position - exon_start)
		elif tmp_genome_position == exon_end and tmp_genome_position > exon_start:
			if i == (len(exon_list)-1):
				genome_position = tmp_genome_position
				return genome_position
			else:
				genome_position = exon_list[i+1][0]
				tmp_transcript_position = tmp_transcript_position - (exon_end - exon_start)
		elif tmp_genome_position < exon_end and tmp_genome_position == exon_start:
			genome_position = tmp_genome_position
			tmp_transcript_position = tmp_transcript_position - (tmp_genome_position - exon_start)
		elif tmp_genome_position < exon_end and tmp_genome_position < exon_start:
			pass
		else:
			raise Exception("transcriptPosition_to_genomePosition error")
	return genome_position


def genome_position_and_motif_annotation(input_signal_list,transcript_to_chrom_dict,transcript_to_exon_dict,transcript_to_strand_dict,chrom2read_dict):
	transid = input_signal_list[0].split(".")[0]
	try:
		chrom = transcript_to_chrom_dict[transid]
	except:
		return []
	tmp_read = chrom2read_dict[chrom]
	tmp_transcript_pos = int(input_signal_list[1])
	tmp_exon_list = transcript_to_exon_dict[transid]
	tmp_strand = transcript_to_strand_dict[transid]
	####################
	transcript_length = exon_length_calculate(exon_list=tmp_exon_list)
	####################
	if tmp_strand == "+":
		out_genome_pos = transcriptPosition_to_genomePosition(transcript_position=tmp_transcript_pos,exon_list=tmp_exon_list)
		out_genome_pos = out_genome_pos + 4
		out_motif = tmp_read[out_genome_pos-6:out_genome_pos+7]
	else:
		rev_transcript_pos = transcript_length - tmp_transcript_pos
		out_genome_pos = transcriptPosition_to_genomePosition(transcript_position=rev_transcript_pos,exon_list=tmp_exon_list)
		out_genome_pos = out_genome_pos - 5
		out_motif = tmp_read[out_genome_pos-6:out_genome_pos+7]
		out_motif = reverse(motif=out_motif)
	##################################
	out_genome_pos = str(out_genome_pos)
	outlist = [chrom,out_genome_pos,out_motif,tmp_strand]
	return outlist
	##################################

def reverse(motif):
	motif.upper()
	convert_dict = {"A":"T","T":"A","C":"G","G":"C","N":"N"}
	outlist = []
	for i in range(len(motif)):
		x = motif[i]
		y = convert_dict[x]
		outlist.append(y)
	outlist.reverse()
	outstr = "".join(outlist)
	return outstr

# questiuon !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
def successive_judge(position_list):
	for i in range(len(position_list)-1):
		x = position_list[i]
		y = position_list[i+1]
		if y-x == 1:
			pass
		else:
			return False
	return True

def signal_split_by_pos(position_list):
	outlist = [[position_list[0]]]
	for i in range(1,len(position_list)):
		x = position_list[i]
		y = outlist[-1][-1]
		if x - y == 1:
			outlist[-1].append(x)
		else:
			outlist.append([x])
	return outlist

def sequecne_make(motif_list):
	tmplist = []
	for i in range(len(motif_list)):
		base_i = motif_list[i][0]
		tmplist.append(base_i)
	tmpstr = "".join(tmplist) + motif_list[-1][1:]
	return tmpstr

def signal_merge_and_result_write(tmp_signal_list,length,writefile):
	#print(len(tmp_signal_list))
	#print(tmp_signal_list[0])
	strand = tmp_signal_list[0][3]
	for i in range(length,len(tmp_signal_list)-length):
		tmp_value_list = []
		for j in range(i-length,i+length+1):
			tmp_value_list.append(tmp_signal_list[j])
		#####################
		#print(len(tmp_value_list))
		#####################
		tmp_value_list.sort(key=lambda x:int(x[5]))
		#tmp_value_list.sort()
		#####################
		tmp_pos_list = [int(x[5]) for x in tmp_value_list]
		tmp_motif_list = [x[6] for x in tmp_value_list]
		#normalize_value_i = float('%.4g' % normalize_value_i)
		tmp_mean_list = [str(float('%.4g' % float(x[10]))) for x in tmp_value_list]
		tmp_stdv_list = [str(float('%.4g' % float(x[11]))) for x in tmp_value_list]
		tmp_length_list = [str(float('%.4g' % float(x[12]))) for x in tmp_value_list]
		###################################################
		tmp_skew_list = [str(float('%.4g' % float(x[17]))) for x in tmp_value_list]
		tmp_kurtosis_list = [str(float('%.4g' % float(x[18]))) for x in tmp_value_list]
		tmp_median_list = [str(float('%.4g' % float(x[19]))) for x in tmp_value_list]
		tmp_quantile_25_list = [str(float('%.4g' % float(x[20]))) for x in tmp_value_list]
		tmp_quantile_75_list = [str(float('%.4g' % float(x[21]))) for x in tmp_value_list]
		tmp_max_list = [str(float('%.4g' % float(x[22]))) for x in tmp_value_list]
		tmp_min_list = [str(float('%.4g' % float(x[23]))) for x in tmp_value_list]
		####################################################
		#####################
		current_list = tmp_signal_list[i].copy()
		#####################
		#####################
		if successive_judge(position_list=tmp_pos_list):
			if strand == "+":
				tmp_mean_str = "|".join(tmp_mean_list)
				tmp_stdv_str = "|".join(tmp_stdv_list)
				tmp_length_str = "|".join(tmp_length_list)
				tmp_skew_str = "|".join(tmp_skew_list)
				tmp_kurtosis_str = "|".join(tmp_kurtosis_list)
				tmp_median_str = "|".join(tmp_median_list)
				tmp_quantile_25_str = "|".join(tmp_quantile_25_list)
				tmp_quantile_75_str = "|".join(tmp_quantile_75_list)
				tmp_max_str = "|".join(tmp_max_list)
				tmp_min_str = "|".join(tmp_min_list)
				#################################
				tmp_sequecne_str = sequecne_make(motif_list=tmp_motif_list)
			else:
				tmp_mean_str = "|".join(tmp_mean_list)
				tmp_stdv_str = "|".join(tmp_stdv_list)
				tmp_length_str = "|".join(tmp_length_list)
				tmp_skew_str = "|".join(tmp_skew_list)
				tmp_kurtosis_str = "|".join(tmp_kurtosis_list)
				tmp_median_str = "|".join(tmp_median_list)
				tmp_quantile_25_str = "|".join(tmp_quantile_25_list)
				tmp_quantile_75_str = "|".join(tmp_quantile_75_list)
				tmp_max_str = "|".join(tmp_max_list)
				tmp_min_str = "|".join(tmp_min_list)
				##################################
				tmp_sequecne_str = sequecne_make(motif_list=tmp_motif_list)
			####################
			out_chrom = current_list[0]
			out_genome_pos = current_list[1]
			out_genome_motif = current_list[2]
			out_strand = current_list[3]
			out_transid = current_list[4]
			out_trans_pos = current_list[5]
			out_read_name = current_list[7]
			####################
			output_list = [
			out_chrom,out_genome_pos,out_genome_motif,out_strand,out_transid,out_trans_pos,out_read_name,
			tmp_sequecne_str,tmp_mean_str,tmp_stdv_str,tmp_length_str,
			tmp_skew_str,tmp_kurtosis_str,tmp_median_str,tmp_quantile_25_str,tmp_quantile_75_str,tmp_max_str,tmp_min_str
			]
			output_str = "\t".join(output_list) + "\n"
			####################
			writefile.write(output_str.encode())
			####################
			####################
		else:
			pass

def args_make():
	import argparse
	parser = argparse.ArgumentParser(description='m6A peak quality control by motif and miCLIP dataset overlap')
	parser.add_argument('--nanopolish_filename', required=True, help="nanopolish result filename")
	parser.add_argument('--signal_filename', required=True, help="output signal filename")
	parser.add_argument('--fasta_filename', required=True, help="fasta file name")
	parser.add_argument('--gtf_filename', required=True, help="gtf filename")
	parser.add_argument('--sample_select_length', required=True, help="sample_select_length")
	parser.add_argument('--signal_extend_length', required=True, help="signal_extend_length")
	args = parser.parse_args()
	return args


def main():
	args = args_make()
	nanopolish_filename = args.nanopolish_filename
	signal_filename = args.signal_filename
	fasta_filename = args.fasta_filename
	gtf_filename = args.gtf_filename
	sample_select_length = int(args.sample_select_length)
	signal_extend_length = int(args.signal_extend_length)
	#########################################
	out_transcript_to_exon_dict,out_transcript_to_strand_dict,out_transcript_to_chrom_dict = gtf_fileread(filename=gtf_filename)
	print("trainscrip result loaded")
	out_chrom_to_read_dict = fasta_fileread(filename=fasta_filename)
	print("genome result loaded")
	#########################################
	nanopolish_fileread(
		nanopolish_filename=nanopolish_filename,
		signal_filename=signal_filename,
		chrom2read_dict=out_chrom_to_read_dict,
		transcript_to_chrom_dict=out_transcript_to_chrom_dict,
		transcript_to_exon_dict=out_transcript_to_exon_dict,
		transcript_to_strand_dict=out_transcript_to_strand_dict,
		sample_select_length=sample_select_length,
		signal_extend_length=signal_extend_length
		)
if __name__=="__main__":
	main()