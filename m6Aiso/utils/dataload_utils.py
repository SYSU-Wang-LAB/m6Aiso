import sys
import os
import random
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold



def normalize_value_list_method_by_max_and_min(max_value_list,min_value_list,input_value_list):
	"""
	normalize value by giving max and min value
	normal_value = (value - min_value)/(max_value-min_value)
	"""
	a = len(max_value_list)
	b = len(max_value_list[0])
	out_value_list = []
	for i in range(a):
		tmplist = []
		for j in range(b):
			max_ij = max_value_list[i][j]
			min_ij = min_value_list[i][j]
			tmp_ij = input_value_list[i][j]
			normalize_ij = (tmp_ij - min_ij)/(max_ij - min_ij)
			tmplist.append(normalize_ij)
		out_value_list.append(tmplist)
	return out_value_list

def normalize_value_list_method_by_mean_and_std(mean_value_list,std_value_list,input_value_list):
	"""
	normalize value by giving max and min value
	normal_value = (value - mean)/(std)
	"""
	a = len(mean_value_list)
	b = len(mean_value_list[0])
	out_value_list = []
	for i in range(a):
		tmplist = []
		for j in range(b):
			mean_ij = mean_value_list[i][j]
			std_ij = std_value_list[i][j]
			tmp_ij = input_value_list[i][j]
			normalize_ij = (tmp_ij - mean_ij)/(std_ij)
			tmplist.append(normalize_ij)
		out_value_list.append(tmplist)
	return out_value_list


def file2dict(filename,max_value_list,min_value_list):
	################
	emdedding_dict = embedding_dict_make()
	################
	df = pd.read_csv(filename,sep='\t')
	out_name_list = []
	out_name2value_dict = {}
	out_name2motif_dict = {}
	out_sitename_list = []
	out_sitename2name_dict = {}
	for index,row in df.iterrows():
		value_list,motif_list = row2valuelist(row=row,emdedding_dict=emdedding_dict)
		name = row["kmer_contig_readindex_tranpos"]
		motif = row["kmer"]
		try:
			chrom = row["chrom"]
			genepos = row["genepos"]
			sitename = chrom+";"+str(genepos)+";"+motif
		except:
			sitename = "_".join(name.split("_")[0:2]+name.split("_")[3:])
		###############
		###############
		normal_value_list = normalize_value_list_method_by_max_and_min(
			max_value_list=max_value_list,
			min_value_list=min_value_list,
			input_value_list=value_list
			)
		out_name_list.append(name)
		##### normal_value_list
		out_name2value_dict[name] = normal_value_list
		out_name2motif_dict[name] = motif_list
		try:
			out_sitename2name_dict[sitename].append(name)
		except:
			out_sitename2name_dict[sitename] = [name]
			out_sitename_list.append(sitename)
	return out_name_list,out_name2value_dict,out_name2motif_dict,out_sitename_list,out_sitename2name_dict


def row2valuelist(row,emdedding_dict=None):
	tmplist = []
	###################################
	kmer = row["kmer"]
	baseflank = row["baseflank"]
	name = row["kmer_contig_readindex_tranpos"]
	###################################
	N1_mean = float(row["N1_mean"])
	P0_mean = float(row["P0_mean"])
	P1_mean = float(row["P1_mean"])
	N1_std= float(row["N1_std"])
	P0_std= float(row["P0_std"])
	P1_std= float(row["P1_std"])
	N1_length = float(row["N1_length"])
	P0_length = float(row["P0_length"])
	P1_length = float(row["P1_length"])
	#####################################
	motif_list = oneHot_encoding_sequence(sequence=baseflank)
	#####################################
	#motif_list = embedding_sequence(sequence=baseflank,embedding_dict=emdedding_dict)
	value_list = []
	value_list.append([N1_mean,N1_std,N1_length])
	value_list.append([P0_mean,P0_std,P0_length])
	value_list.append([P1_mean,P1_std,P1_length])
	return value_list,motif_list
	#####################################

def oneHot_encoding_sequence(sequence):
	sequence = sequence
	outlist = []
	for i in range(len(sequence)):
		base_i = sequence[i]
		if base_i == "A":
			outlist.append([1.0,0.0,0.0,0.0])
		elif base_i == "T":
			outlist.append([0.0,1.0,0.0,0.0])
		elif base_i == "C":
			outlist.append([0.0,0.0,1.0,0.0])
		elif base_i == "G":
			outlist.append([0.0,0.0,0.0,1.0])
		elif base_i == "N":
			outlist.append([0.25,0.25,0.25,0.25])
		else:
			raise Exception("oneHot_encoding_sequence error")
	return outlist

def label_dict_merge(pos_name_list,neg_name_list):
	valuedict = {}
	for i in range(len(pos_name_list)):
		name = pos_name_list[i]
		valuedict[name] = 1.0
	for j in range(len(neg_name_list)):
		name = neg_name_list[j]
		valuedict[name] = 0.0
	return valuedict

def value_dict_merge(pos_name2value_dict,neg_name2value_dict):
	merge_name2value_dict = {}
	for key,value in pos_name2value_dict.items():
		merge_name2value_dict[key] = value
	for key,value in neg_name2value_dict.items():
		merge_name2value_dict[key] = value
	return merge_name2value_dict


def random_disorganize_train_namelist(pos_name_list,neg_name_list):
	merge_name_list = pos_name_list + neg_name_list
	random.seed(1234)
	random.shuffle(merge_name_list)
	random.seed(1234)
	random.shuffle(merge_name_list)
	random.seed(1234)
	random.shuffle(merge_name_list)
	return merge_name_list

def random_disorganize_using_namelist(namelist):
	random.seed(1234)
	random.shuffle(namelist)
	random.seed(1234)
	random.shuffle(namelist)
	random.seed(1234)
	random.shuffle(namelist)
	return namelist

def train_numpy_array_make(namelist,merge_value_x_dict,merge_value_y_dict,merge_value_m_dict):
	x_list = []
	y_list = []
	m_list = []
	for i in range(len(namelist)):
		name_i = namelist[i]
		y_value_i = [merge_value_y_dict[name_i]]
		x_value_i = merge_value_x_dict[name_i]
		m_value_i = merge_value_m_dict[name_i]
		x_list.append(x_value_i)
		y_list.append(y_value_i)
		m_list.append(m_value_i)
	xlist = np.array(x_list,dtype = np.float32)
	ylist = np.array(y_list,dtype = np.float32)
	mlist = np.array(m_list,dtype = np.float32)
	return xlist,ylist,mlist

def using_numpy_array_make(namelist,merge_value_x_dict,merge_motif_m_dict):
	x_list = []
	m_list = []
	for i in range(len(namelist)):
		name_i = namelist[i]
		x_value_i = merge_value_x_dict[name_i]
		m_value_i = merge_motif_m_dict[name_i]
		x_list.append(x_value_i)
		m_list.append(m_value_i)
	xlist = np.array(x_list,dtype = np.float32)
	mlist = np.array(m_list,dtype = np.float32)
	return xlist,mlist

def namelist_abstract_by_index(namelist,select_index_list):
	outlist = []
	for i in select_index_list:
		x = namelist[i]
		outlist.append(x)
	return outlist

def kfold_split_namelist(namelist,split_number):
	outlist = []
	kf = KFold(n_splits=split_number,shuffle=True,random_state=116)
	for k,(res_index_list,select_index_list) in enumerate(kf.split(namelist)):
		select_name_list_k = namelist_abstract_by_index(namelist=namelist,select_index_list=select_index_list)
		outlist.append(select_name_list_k)
	return outlist

def splited_array_make(splited_train_neg_name_list,train_pos_name_list,merge_value_x_dict,merge_value_y_dict,merge_value_m_dict):
	splited_train_x_list = []
	splited_train_y_list = []
	splited_train_m_list = []
	splited_train_namelist = []
	for k in range(len(splited_train_neg_name_list)):
		k_train_neg_namelist = splited_train_neg_name_list[k]
		k_train_pos_namelist = train_pos_name_list
		k_merge_namelist = random_disorganize_train_namelist(pos_name_list=k_train_pos_namelist,neg_name_list=k_train_neg_namelist)
		k_xlist,k_ylist,k_mlist = train_numpy_array_make(
			namelist=k_merge_namelist,
			merge_value_x_dict=merge_value_x_dict,
			merge_value_y_dict=merge_value_y_dict,
			merge_value_m_dict=merge_value_m_dict
			)
		splited_train_x_list.append(k_xlist)
		splited_train_y_list.append(k_ylist)
		splited_train_m_list.append(k_mlist)
		splited_train_namelist.append(k_merge_namelist)
	return splited_train_x_list,splited_train_y_list,splited_train_m_list,splited_train_namelist

def sitename_splited_array_make(splited_merge_sitename_list,merge_value_x_dict,merge_value_y_dict,merge_motif_m_dict,merge_sitename2name_dict):
	out_array_x_list = []
	out_value_y_list = []
	out_motif_m_list = []
	#################################################
	for k in range(len(splited_merge_sitename_list)):
		k_merge_sitename_list = splited_merge_sitename_list[k]
		k_array_x_list = []
		k_value_y_list = []
		k_motif_m_list = []
		for sitename_i in k_merge_sitename_list:
			merge_namelist_i = merge_sitename2name_dict[sitename_i]
			merge_namelist_i = random_disorganize_using_namelist(namelist=merge_namelist_i)
			#########################################
			x_i_array,m_i_array = using_numpy_array_make(namelist=merge_namelist_i,merge_value_x_dict=merge_value_x_dict,merge_motif_m_dict=merge_motif_m_dict)
			y_i_value = merge_value_y_dict[sitename_i]
			if len(merge_namelist_i) < 50:
				#k_array_x_list.append(x_i_array)
				#k_value_y_list.append(y_i_value)
				pass
			else:
				k_array_x_list.append(x_i_array)
				k_value_y_list.append(y_i_value)
				k_motif_m_list.append(m_i_array)
		out_array_x_list.append(k_array_x_list)
		out_value_y_list.append(k_value_y_list)
		out_motif_m_list.append(k_motif_m_list)
	#################################################
	return out_array_x_list,out_value_y_list,out_motif_m_list

def motif2readname_dict_make(read_name_list):
	outdict = {}
	outlist = []
	for name in read_name_list:
		motif = name.split("_")[0]
		try:
			outdict[motif].append(name)
		except:
			outdict[motif] = [name]
			outlist.append(motif)
	return outdict,outlist

def motif_splited_array_make(pos_motif2readname_dict,neg_motif2readname_dict,pos_motif_list,neg_motif_list,merge_value_x_dict,merge_value_y_dict,multiN):
	splited_train_x_list = []
	splited_train_y_list = []
	splited_train_namelist = []
	for motif in pos_motif_list:
		k_pos_name_list = pos_motif2readname_dict[motif]
		k_neg_name_list = neg_motif2readname_dict[motif]
		n_pos_name_list = k_pos_name_list * multiN
		n_neg_name_list = k_neg_name_list
		k_merge_namelist = random_disorganize_train_namelist(pos_name_list=n_pos_name_list,neg_name_list=n_neg_name_list)
		k_xlist,k_ylist = train_numpy_array_make(namelist=k_merge_namelist,merge_value_x_dict=merge_value_x_dict,merge_value_y_dict=merge_value_y_dict)
		splited_train_x_list.append(k_xlist)
		splited_train_y_list.append(k_ylist)
		splited_train_namelist.append(k_merge_namelist)
	return splited_train_x_list,splited_train_y_list,splited_train_namelist


def value_list_write(input_value_list,out_filename):
	d = open(out_filename,'a')
	for i in range(len(input_value_list)):
		tmplist = [str(x) for x in input_value_list[i]]
		tmpstr = "\t".join(tmplist) + "\n"
		d.write(tmpstr)
	d.close()

def extrem_value_load_form_file(filename):
	f = open(filename,'r')
	outlist = []
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "mean":
			continue
		tmplist = [float(x) for x in list_x]
		outlist.append(tmplist)
	return outlist

def file2namelist(filename):
	f = open(filename,'r')
	namelist = []
	for str_x in f:
		str_x = str_x.strip("\n")
		namelist.append(str_x)
	return namelist

def embedding_sequence(sequence,embedding_dict):
	#############
	outlist = []
	for i in range(len(sequence)-4):
		motif_i = sequence[i:i+5]
		value = embedding_dict[motif_i]
		outlist.append(value)
	return outlist

def embedding_dict_make():
	codingdict = {
	'A':['A'],
	'T':['T'],
	'C':['C'],
	'G':['G'],
	'N':["A","T","G","G"],
	'R':["A","G"],
	'D':["A","T","G"],
	'H':["A","T","C"]
	}
	value = 1
	outdict = {}
	############################
	outlist = []
	for x in "NDRAC":
		tmplist = codingdict[x]
		outlist.append(tmplist)
	for i0 in range(len(outlist[0])):
		for i1 in range(len(outlist[1])):
			for i2 in range(len(outlist[2])):
				for i3 in range(len(outlist[3])):
					for i4 in range(len(outlist[4])):
						tmpmotif = "".join([outlist[0][i0],outlist[1][i1],outlist[2][i2],outlist[3][i3],outlist[4][i4]])
						outdict[tmpmotif] = value
						value += 1
	############################
	outlist = []
	for x in "DRACH":
		tmplist = codingdict[x]
		outlist.append(tmplist)
	for i0 in range(len(outlist[0])):
		for i1 in range(len(outlist[1])):
			for i2 in range(len(outlist[2])):
				for i3 in range(len(outlist[3])):
					for i4 in range(len(outlist[4])):
						tmpmotif = "".join([outlist[0][i0],outlist[1][i1],outlist[2][i2],outlist[3][i3],outlist[4][i4]])
						outdict[tmpmotif] = value
						value += 1
	############################
	outlist = []
	for x in "RACHN":
		tmplist = codingdict[x]
		outlist.append(tmplist)
	for i0 in range(len(outlist[0])):
		for i1 in range(len(outlist[1])):
			for i2 in range(len(outlist[2])):
				for i3 in range(len(outlist[3])):
					for i4 in range(len(outlist[4])):
						tmpmotif = "".join([outlist[0][i0],outlist[1][i1],outlist[2][i2],outlist[3][i3],outlist[4][i4]])
						outdict[tmpmotif] = value
						value += 1
	return outdict