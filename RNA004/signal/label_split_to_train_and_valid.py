

import gzip
import random


def label_fileread(filename,train_label_filename,valid_label_filename,valid_percentage):
	f = gzip.open(filename,'rb')
	d1 = gzip.open(valid_label_filename,'ab')
	d2 = gzip.open(train_label_filename,'ab')
	for str_x in f:
		str_x = str_x.decode()
		#######################################
		v = random.random()
		#######################################
		if v <= valid_percentage:
			str2write = str_x
			d1.write(str2write.encode())
		else:
			str2write = str_x
			d2.write(str2write.encode())

def args_make():
	import argparse
	parser = argparse.ArgumentParser(description='label file make based on motif')
	parser.add_argument('--org_label_filename', required=True, help="pos_filename")
	parser.add_argument('--train_label_filename', required=True, help="pos_filename")
	parser.add_argument('--valid_label_filename', required=True, help="neg_filename")
	parser.add_argument('--valid_percentage', required=True, help="noise_percentage")
	args = parser.parse_args()
	return args


if __name__=="__main__":
	import sys
	args = args_make()
	org_label_filename = args.org_label_filename
	train_label_filename = args.train_label_filename
	valid_label_filename = args.valid_label_filename
	valid_percentage = float(args.valid_percentage)
	label_fileread(
		filename=org_label_filename,
		train_label_filename=train_label_filename,
		valid_label_filename=valid_label_filename,
		valid_percentage=valid_percentage
		)