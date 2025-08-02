import gzip
import random
import math
import numpy as np


def signal_fileread(filename,output_filename,m6Adict,negative_motif2value_dict=None):
	"""
	#chrom   genome_position genome_motif    strand  transcript_id   transcript_position     readname        up_and_down_sequecne       mean    stdv    length  samples
	#chr1    633122  AAGTA   -       ENST00000440196.3       4       3a8b3f71-7bf0-4c63-8240-fd36d81dbe77    GGAAGTATG  114.4707142857143|103.11469387755102|111.53|126.33|92.90520000000001    5.291396428571429|4.15738775510204|3.496|9.205|2.95144     0.09297|0.01627|0.00896|0.00266|0.0083  118.078,125.686,119.103,103.74,111.787,106.959,116.323,114.421,108.13,120.273,119.103,103.155,121.882,110.178,113.543,117.932,104.179,123.638,110.178,121.882|109.446,105.789,96.1323,90.1336,104.326,100.522,102.716,100.668,103.155,105.35,107.691,95.547,113.25,101.985,108.568,98.3269,96.4249,100.083,96.4249,94.8155|111.641,119.542,111.787,110.617,116.03,105.642,106.666,110.324,110.032,111.641,109.007,111.934,113.104,106.813,114.421,105.935,110.47,109.3,116.323,114.567|117.64,143.682,120.127,131.539,134.465,115.738,134.465,116.176,115.738,120.127,120.127,143.682,115.738,131.246,116.176,117.64,117.64,116.176,131.246,134.465|94.8155,90.1336,92.9135,92.7672,92.3282,91.5967,89.5484,90.5726,94.8155,92.1819,92.9135,105.789,88.6705,88.6705,92.4746,92.0356,91.3041,94.8155,94.6692,88.8169
	"""
	#print(filename)
	#print(output_filename)
	f = gzip.open(filename,'rb')
	writefile = gzip.open(output_filename,'ab')
	##############################
	for str_x in f:
		str_x = str_x.decode()
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "chrom":
			if list_x[-1] == "base_calling_sequecne":
				output_list = ['chrom','genome_position','genome_motif','strand','transcript_id','transcript_position','readname','up_and_down_sequecne','mean','stdv','length','samples','base_calling_sequecne']
				output_str = "\t".join(output_list) + "\n"
				writefile.write(output_str.encode())
			else:
				output_list = ['chrom','genome_position','genome_motif','strand','transcript_id','transcript_position','readname','up_and_down_sequecne','mean','stdv','length','samples']
				output_str = "\t".join(output_list) + "\n"
				writefile.write(output_str.encode())
			continue
		#######################
		chrom = list_x[0]
		genenome_pos = int(list_x[1])
		motif = list_x[2]
		strand = list_x[3]
		if negative_motif2value_dict == None:
			motif_count = 10000
		else:
			try:
				motif_count = negative_motif2value_dict[motif]
			except:
				motif_count = 0
		#######################
		#print(chrom,genenome_pos,strand)
		#######################
		try:
			v = m6Adict[(chrom,genenome_pos,strand)]
			#print(chrom,genenome_pos,strand)
			#print("#######################################")
			#print("#######################################")
		except:
			v = 0
		if v>0 and motif_count>=100:
			output_str = str_x + "\n"
			writefile.write(output_str.encode())
	##############################

def site_fileread(filename,base_type):
	"""
	Chr     Start   End     sitename        Strand  motif   Gene    DRACH   Ratio
	chr16   68860162        68860163        chr16;68860162;+;TANGO6;GGACT   +       GGACT   TANGO6  True    0.85351
	chr2    96186360        96186361        chr2;96186360;-;STARD7;GGATT    -       GGATT   STARD7  False   0.06545000000000001
	
	chr14   73219165        73219166        PSEN1   +       ATTTT   0.09437999999999999
	chr21   43774289        43774290        CSTB    -       GTTCC   0.10696900000000001
	chr15   76285711        76285712        ETFA    -       TATCA   0.1054245
	"""
	f = open(filename,'r')
	outdict = {}
	#outlist = []
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		##################
		if list_x[0] == "Chr":
			continue
		chrom = list_x[0]
		position = int(list_x[1])
		strand = list_x[4]
		motif = list_x[5]
		try:
			DRACH = list_x[7]
			############
			if DRACH == "False":
				pass
		except:
			pass
		############
		if base_type == "1base":
			if strand == "+":
				position = position-1 # m6A
			else:
				position = position-1 # m6A
		if base_type == "0base":
			if strand == "+":
				position = position # Pseu
			else:
				position = position # Pseu
		##################
		##################
		#print(chrom,position,strand)
		try:
			outdict[(chrom,position,strand)] += 1
		except:
			outdict[(chrom,position,strand)] = 1
	return outdict


def negative_motif2value_dict_make(negative_motif_filename):
	f = open(negative_motif_filename,'r')
	out_negative_motif2value_dict = {}
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		#############################
		motif = list_x[0]
		value = int(list_x[1])
		if value >= 200:
			out_negative_motif2value_dict[motif] = value
		else:
			pass
	return out_negative_motif2value_dict
	


def args_make():
	import argparse
	parser = argparse.ArgumentParser(description='m6A peak quality control by motif and miCLIP dataset overlap')
	parser.add_argument('--modification_site_format_filename', required=True, help="m6A site file")
	parser.add_argument('--signal_filename', required=True, help="signal file")
	parser.add_argument('--result_filename', required=True, help="result file")
	parser.add_argument('--base_type', required=True,default="1base", help="result file")
	parser.add_argument('--negative_motif_filename',default=None, help="result file")
	args = parser.parse_args()
	return args


if __name__=="__main__":
	import sys
	args = args_make()
	modification_site_format_filename = args.modification_site_format_filename
	signal_filename = args.signal_filename
	result_filename = args.result_filename
	base_type = args.base_type
	negative_motif_filename = args.negative_motif_filename
	###############################
	outdict = site_fileread(filename=modification_site_format_filename,base_type=base_type)
	if negative_motif_filename == None:
		signal_fileread(filename=signal_filename,output_filename=result_filename,m6Adict=outdict,negative_motif2value_dict=None)
	else:
		out_negative_motif2value_dict = negative_motif2value_dict_make(negative_motif_filename=negative_motif_filename)
		signal_fileread(filename=signal_filename,output_filename=result_filename,m6Adict=outdict,negative_motif2value_dict=out_negative_motif2value_dict)