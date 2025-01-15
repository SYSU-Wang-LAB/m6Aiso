

def signal_file_split(filename,site2value_dict,max_value,min_value):
	"""
	chrom   genome_pos      strand  motif   kmer    kmer_contig_readindex_tranpos   P1_mean P1_std  P1_length   P0_mean  P0_std  P0_length       N1_mean N1_std  N1_length       baseflank
	chr11   123059681       -       GAACT   GAACT   GAACT_ENST00000533540_655672_551        89.0    1.5970000000000002   0.0067224390243902435   98.3    2.7369999999999997      0.011667499999999999    123.8   9.194666666666667    0.002055        AGAACTG
	chr11   123059670       -       TGACC   TGACC   TGACC_ENST00000533540_655672_562        78.3    2.8004000000000002   0.0052  117.5   10.837  0.006   106.1   3.989   0.01967 CTGACCT
	chr11   123059649       -       GGACC   GGACC   GGACC_ENST00000533540_655672_583        75.5    2.464513274336283    0.013920265486725664    111.6   13.009666666666668      0.006666666666666666    119.1   4.812884615384615    0.01053923076923077     TGGACCC
	"""
	f = open(filename,'r')
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "chrom":
			print("\t".join(list_x[4:]))
			continue
		chrom = list_x[0]
		start = list_x[1]
		strand = list_x[2]
		#print(chrom,start,strand)
		try:
			value = site2value_dict[(chrom,start,strand)]
		except:
			value = -1
		if value >= min_value and value < max_value:
			print("\t".join(list_x[4:]))
		else:
			pass

def glori_fileread(filename):
	"""
	Chr     Start   End     sitename        Strand  motif   Gene    DRACH   Ratio
	chr16   68860162        68860163        chr16;68860162;+;TANGO6;GGACT   +       GGACT   TANGO6  True    0.85351
	chr22   42693245        42693246        chr22;42693245;-;A4GALT;GGACT   -       GGACT   A4GALT  True    0.953365
	chr22   46534022        46534023        chr22;46534022;-;CELSR1;GGACC   -       GGACC   CELSR1  True    0.181235
	"""
	f = open(filename,'r')
	outdict = {}
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "Chr":
			continue
		chrom = list_x[0]
		start = list_x[1]
		strand = list_x[4]
		value = float(list_x[-1])
		#print(chrom,start,strand)
		###############################################
		outdict[(chrom,start,strand)] = value
	return outdict

if __name__ == "__main__":
	import sys
	filename = sys.argv[1]
	min_value = float(sys.argv[2])
	max_value = float(sys.argv[3])
	site2value_dict = glori_fileread(filename="/home/ZJRen/public3/m6Apred/m6A/HEK293T_GLORI_m6A_merged.txt")
	signal_file_split(
		filename=filename,
		site2value_dict=site2value_dict,
		max_value=max_value,
		min_value=min_value
		)