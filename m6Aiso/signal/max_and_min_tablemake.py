

import sys
sys.path.insert(0,"../../")
sys.path.insert(0,"../../m6Aiso/")
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
	df = pd.read_csv(filename,sep='\t',low_memory=False)
	colnames_list = df.columns.tolist()
	##########################################
	max_value = df.max()
	min_value = df.min()
	max_value_list,_ = row2valuelist(row=max_value,colnames_list=colnames_list)
	min_value_list,_ = row2valuelist(row=min_value,colnames_list=colnames_list)
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
	work_dir = sys.argv[1]
	merge_total_filename = work_dir + sys.argv[2]
	merge_max_filename = work_dir + sys.argv[3]
	merge_min_filename = work_dir + sys.argv[4]
	##########################################
	##########################################
	merge_max_value_list,merge_min_value_list = extrem_value_make(filename=merge_total_filename)
	try:
		os.remove(merge_max_filename)
		os.remove(merge_min_filename)
	except:
		pass
	value_list_write(input_value_list=merge_max_value_list,out_filename=merge_max_filename)
	value_list_write(input_value_list=merge_min_value_list,out_filename=merge_min_filename)


if __name__ == "__main__":
	main()