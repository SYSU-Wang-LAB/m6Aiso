

import sys
import os
import random
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
from utils.dataload_utils import file2dict
from utils.dataload_utils import row2valuelist
from utils.dataload_utils import oneHot_encoding_sequence

def extrem_value_make(filename):
	"""
	max and min value abstract from dataset
	"""
	df = pd.read_csv(filename,sep='\t')
	max_value = df.max()
	min_value = df.min()
	max_value_list,_ = row2valuelist(row=max_value)
	min_value_list,_ = row2valuelist(row=min_value)
	###########################
	return max_value_list,min_value_list

def value_list_write(input_value_list,out_filename):
	d = open(out_filename,'a')
	tmpstr = "\t".join(["mean","std","length"]) + "\n"
	d.write(tmpstr)
	for i in range(len(input_value_list)):
		tmplist = [str(x) for x in input_value_list[i]]
		tmpstr = "\t".join(tmplist) + "\n"
		d.write(tmpstr)
	d.close()

def main():
	import sys
	merge_total_filename = "/home/ZJRen/pytorch/attention/data/merge_total.sorted.baseflank.tsv"
	merge_max_value_list,merge_min_value_list = extrem_value_make(filename=merge_total_filename)
	try:
		os.remove("/home/ZJRen/pytorch/attention/data/merge_max_value_list.txt")
		os.remove("/home/ZJRen/pytorch/attention/data/merge_min_value_list.txt")
	except:
		pass
	value_list_write(input_value_list=merge_max_value_list,out_filename="/home/ZJRen/pytorch/attention/data/merge_max_value_list.txt")
	value_list_write(input_value_list=merge_min_value_list,out_filename="/home/ZJRen/pytorch/attention/data/merge_min_value_list.txt")


if __name__ == "__main__":
	main()