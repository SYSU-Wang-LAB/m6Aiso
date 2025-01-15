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


def signal_file_annotation(filename,transcript_to_exon_dict,transcript_to_strand_dict,transcript_to_chrom_dict,chrom_to_read_dict):
	"""
	kmer    kmer_contig_readindex_tranpos   P1_mean P1_std  P1_length       P0_mean P0_std  P0_length       N1_mean      N1_std  N1_length       baseflank
	GAACT   GAACT_ENST00000533540_655672_551        89.0    1.5970000000000002      0.0067224390243902435   98.32.7369999999999997       0.011667499999999999    123.8   9.194666666666667       0.002055        AGAACTG
	TGACC   TGACC_ENST00000533540_655672_562        78.3    2.8004000000000002      0.0052  117.5   10.837  0.006106.1   3.989   0.01967 CTGACCT
	"""
	f = open(filename,'r')
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "kmer":
			print("\t".join(['chrom','genome_pos',"strand",'motif']+list_x))
			continue
		sitename = list_x[1]
		transcript_id = sitename.split("_")[1]
		transcript_pos = int(sitename.split("_")[3])
		##########################
		try:
			exon_list = transcript_to_exon_dict[transcript_id]
			strand = transcript_to_strand_dict[transcript_id]
			chrom = transcript_to_chrom_dict[transcript_id]
			##############################################
			read = chrom_to_read_dict[chrom]
		except:
			continue
		if strand == "+":
			genome_pos = transcriptPosition_to_genomePosition(transcript_position=transcript_pos,exon_list=exon_list)
			genome_pos = genome_pos + 1
			motif = read[genome_pos-3:genome_pos+2]
			motif.upper()
		else:
			transcript_length = exon_length_calculate(exon_list=exon_list)
			transcript_pos = transcript_length - transcript_pos
			genome_pos = transcriptPosition_to_genomePosition(transcript_position=transcript_pos,exon_list=exon_list)
			motif = read[genome_pos-3:genome_pos+2]
			motif = reverse(motif=motif)
		print("\t".join([chrom,str(genome_pos),strand,motif]+list_x))
		##############################################



def main():
	import sys
	gtf_filename = "/home/ZJRen/index.file/hg38_index.file/longread/Homo_sapiens.GRCh38.91.GenePred"
	fasta_filename = "/home/ZJRen/index.file/hg38_index.file/10X_resource/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
	(transcript_to_exon_dict,transcript_to_strand_dict,transcript_to_chrom_dict) = gtf_fileread(filename=gtf_filename)
	(chrom_to_read_dict) = fasta_fileread(filename=fasta_filename)
	signal_file_annotation(
		filename=sys.argv[1],
		transcript_to_exon_dict=transcript_to_exon_dict,
		transcript_to_strand_dict=transcript_to_strand_dict,
		transcript_to_chrom_dict=transcript_to_chrom_dict,
		chrom_to_read_dict=chrom_to_read_dict
		)

if __name__ == "__main__":
	main()
