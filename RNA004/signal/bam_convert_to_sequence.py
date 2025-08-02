import gzip

def bam_tsv_fileread(filename,writefile,transcript2chrom_dict,input_chrom):
	"""
	#Read-Name      Flag    MAPQ    CHROM   READ-POS0       READ-BASE       READ-QUAL       REF-POS1	REF-BASE        CIGAR-OP
	b8b8a1e0-decd-4d21-88eb-260780a7c9d3    0       2       ENST00000390447.3       0       C      $	-18     .       S
	b8b8a1e0-decd-4d21-88eb-260780a7c9d3    0       2       ENST00000390447.3       1       T      &	-17     .       S
	b8b8a1e0-decd-4d21-88eb-260780a7c9d3    0       2       ENST00000390447.3       2       G      '	-16     .       S
	b8b8a1e0-decd-4d21-88eb-260780a7c9d3    0       2       ENST00000390447.3       3       G      '	-15     .       S
	b8b8a1e0-decd-4d21-88eb-260780a7c9d3    0       2       ENST00000390447.3       4       G      &	-14     .       S
	#READ_NAME      FLAG    CHROM   READ_POS        BASE    QUAL    REF_POS REF     OP
	0fd40ad9-b70c-4dc5-9d4e-ec207212b9b4    0       ENST00000456328.2       0       A       &       -332    .       S
	0fd40ad9-b70c-4dc5-9d4e-ec207212b9b4    0       ENST00000456328.2       1       A       (       -331    .       S
	0fd40ad9-b70c-4dc5-9d4e-ec207212b9b4    0       ENST00000456328.2       2       G       (       -330    .       S
	"""
	d = gzip.open(writefile,'ab')
	str2write = "\t".join(['readname','transcript_id','baseread','baseref','refposi','quality']) + "\n"
	str2write = str2write.encode()
	d.write(str2write)
	#######################################################
	#f = open(filename,'r')
	f = gzip.open(filename,'rb')
	tmplist = []
	for str_x in f:
		str_x = str_x.decode()
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "#READ_NAME":
			index_dict = dict(zip(list_x,list(range(len(list_x)))))
			continue
		#print(index_dict)
		readname = list_x[index_dict['#READ_NAME']]
		transcript_id = list_x[index_dict['CHROM']]
		baseread = list_x[index_dict['BASE']]
		baseref = list_x[index_dict['REF']]
		refposi = list_x[index_dict['REF_POS']]
		quality = list_x[index_dict['QUAL']]
		try:
			chrom = transcript2chrom_dict[transcript_id]
		except:
			continue
		if not chrom in input_chrom:
			continue
		if len(tmplist) == 0:
			tmplist.append([readname,transcript_id,baseread,baseref,refposi,quality])
		else:
			if readname == tmplist[-1][0]:
				tmplist.append([readname,transcript_id,baseread,baseref,refposi,quality])
			else:
				resultwrite(tmplist=tmplist,d=d)
				tmplist = []
				tmplist.append([readname,transcript_id,baseread,baseref,refposi,quality])

def annotation_fileread(filename):
	"""
	#trans_id       chrom   strand  trans_start     trans_end       cds_start       cds_end exon_num        exon_start      exon_end gene_id gene_name       gene_type
	ENST00000456328.2       chr1    +       11869   14409   11869   14409   3       11869,12613,13221,      12227,12721,14409,       ENSG00000223972.5       DDX11L1 lncRNA
	ENST00000450305.2       chr1    +       12010   13670   12010   13670   6       12010,12179,12613,12975,13221,13453,    12057,12227,12697,13052,13374,13670,     ENSG00000223972.5       DDX11L1 transcribed_unprocessed_pseudogene
	"""
	f = open(filename,'r')
	outdict = {}
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		############################
		if list_x[0] == "#trans_id":
			continue
		transcript_id = list_x[0]
		chrom = list_x[1]
		############################
		outdict[transcript_id] = chrom
	return outdict

def resultwrite(tmplist,d):
	readname = tmplist[0][0]
	transcript_id = tmplist[0][1]
	baseread = "".join([x[2] for x in tmplist])
	baseref = "".join([x[3] for x in tmplist])
	refposi = ";".join([x[4] for x in tmplist])
	quality = "".join([x[5] for x in tmplist])
	str2write = "\t".join([readname,transcript_id,baseread,baseref,refposi,quality]) + "\n"
	str2write = str2write.encode()
	d.write(str2write)

def args_make():
	import argparse
	parser = argparse.ArgumentParser(description='m6A peak quality control by motif and miCLIP dataset overlap')
	parser.add_argument('--bam_filename', required=True, help="bam_filename")
	parser.add_argument('--sequence_filename', required=True, help="sequence_filename")
	parser.add_argument('--input_chrom', required=True, help="input_chrom")
	parser.add_argument('--annotation_filename', required=True, help="output_filename")
	args = parser.parse_args()
	return args


if __name__ == "__main__":
	import sys
	args = args_make()
	bam_filename = args.bam_filename
	sequence_filename = args.sequence_filename
	input_chrom = args.input_chrom
	annotation_filename = args.annotation_filename
	##############################
	input_chrom = [x for x in input_chrom.split(",")]
	##############################
	transcript2chrom_dict = annotation_fileread(filename=annotation_filename)
	bam_tsv_fileread(
		filename=bam_filename,
		writefile=sequence_filename,
		transcript2chrom_dict=transcript2chrom_dict,
		input_chrom=input_chrom
		)