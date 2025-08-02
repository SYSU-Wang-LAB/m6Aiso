

import gzip
#########################################

def signal_readname_select(filename,writefile,site2value_dict):
	f = gzip.open(filename,'rb')
	d = gzip.open(writefile,'ab')
	readname_to_count = {}
	for str_x in f:
		str_x = str_x.decode()
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		#############################
		if list_x[0] == "chrom":
			index_chrom = list_x.index("chrom")
			index_genome_position = list_x.index("genome_position")
			index_genome_motif = list_x.index("genome_motif")
			index_strand = list_x.index("strand")
			index_transcript_id = list_x.index("transcript_id")
			index_transcript_position = list_x.index("transcript_position")
			index_up_and_down_sequecne = list_x.index("up_and_down_sequecne")
			index_readname = list_x.index("readname")
			str2write = str_x + "\n"
			d.write(str2write.encode())
			continue
		#############################
		chrom = list_x[index_chrom]
		genome_position = list_x[index_genome_position]
		genome_motif = list_x[index_genome_motif]
		strand = list_x[index_strand]
		transcript_id = list_x[index_transcript_id]
		transcript_position = list_x[index_transcript_position]
		up_and_down_sequecne = list_x[index_up_and_down_sequecne]
		readname = list_x[index_readname]
		#############################
		try:
			value = site2value_dict[(chrom,genome_position,strand)]
			str2write = str_x + "\n"
			d.write(str2write.encode())
		except:
			pass
	f.close()
	return readname_to_count

def bed_format_fileread(filename):
	"""
	DRACH_oligo_0   28      29      -       0       +
	DRACH_oligo_0   46      47      -       0       +
	DRACH_oligo_0   64      65      -       0       +
	DRACH_oligo_0   82      83      -       0       +
	"""
	f = open(filename,'r')
	outdict = {}
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		#############################
		chrom = list_x[0]
		position = list_x[1]
		strand = list_x[5]
		outdict[(chrom,position,strand)] = 1
	return outdict


def args_make():
	import argparse
	parser = argparse.ArgumentParser(description='m6A peak quality control by motif and miCLIP dataset overlap')
	parser.add_argument('--signal_filename', required=True, help="m6A peak file")
	parser.add_argument('--bed_format_filename', required=True, help="m6A site file")
	parser.add_argument('--output_filename', required=True, help="m6A site file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	import sys
	args = args_make()
	signal_filename = args.signal_filename
	bed_format_filename = args.bed_format_filename
	output_filename = args.output_filename
	#############################
	site2value_dict = bed_format_fileread(filename=bed_format_filename)
	#############################
	signal_readname_select(
		filename=signal_filename,
		writefile=output_filename,
		site2value_dict=site2value_dict
		)
