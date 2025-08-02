import gzip
import random


def label_fileread(filename,result_label_filename,noise_label_filename,noise_percentage):
	f = gzip.open(filename,'rb')
	d1 = gzip.open(result_label_filename,'ab')
	d2 = gzip.open(noise_label_filename,'ab')
	for str_x in f:
		str_x = str_x.decode()
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		#######################################
		sitename = list_x[0]
		v1 = float(list_x[1])
		v2 = float(list_x[2])
		#######################################
		v3 = random.random()
		#######################################
		if v3 <= noise_percentage:
			adj_v1 = v2
			adj_v2 = v1
			adj_v3 = 1.0
		else:
			adj_v1 = v1
			adj_v2 = v2
			adj_v3 = 0.0
		str2write_1 = "\t".join([sitename,str(adj_v1),str(adj_v2)]) + "\n"
		str2write_2 = "\t".join([sitename,str(adj_v3)]) + "\n"
		d1.write(str2write_1.encode())
		d2.write(str2write_2.encode())

def args_make():
	import argparse
	parser = argparse.ArgumentParser(description='label file make based on motif')
	parser.add_argument('--org_label_filename', required=True, help="pos_filename")
	parser.add_argument('--new_label_filename', required=True, help="pos_filename")
	parser.add_argument('--noise_label_filename', required=True, help="neg_filename")
	parser.add_argument('--noise_percentage', required=True, help="noise_percentage")
	args = parser.parse_args()
	return args


if __name__=="__main__":
	import sys
	args = args_make()
	org_label_filename = args.org_label_filename
	new_label_filename = args.new_label_filename
	noise_label_filename = args.noise_label_filename
	noise_percentage = float(args.noise_percentage)
	label_fileread(
		filename=org_label_filename,
		result_label_filename=new_label_filename,
		noise_label_filename=noise_label_filename,
		noise_percentage=noise_percentage
		)