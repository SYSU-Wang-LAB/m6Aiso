import torch
import numpy as np
from torch.utils.data import Dataset,DataLoader
##################################### extrem_value_load_form_file
from .dataload_utils import extrem_value_load_form_file
from .dataload_utils import file2dict
from .dataload_utils import label_dict_merge
from .dataload_utils import value_dict_merge
from .dataload_utils import kfold_split_namelist
from .dataload_utils import splited_array_make
from .dataload_utils import random_disorganize_train_namelist
from .dataload_utils import random_disorganize_using_namelist
from .dataload_utils import sitename_splited_array_make
from .dataload_utils import motif_splited_array_make
from .dataload_utils import motif2readname_dict_make
from .dataload_utils import using_numpy_array_make

def dataload_for_read_level_m6A_pred(
	pos_filename,
	neg_filename,
	split_number=10,
	max_value_filename="/home/ZJRen/pytorch/NGS/data/merge_max_value_list.txt",
	min_value_filename="/home/ZJRen/pytorch/NGS/data/merge_min_value_list.txt"):
	merge_max_value_list = extrem_value_load_form_file(filename=max_value_filename)
	merge_min_value_list = extrem_value_load_form_file(filename=min_value_filename)
	pos_name_list,pos_name2value_dict,pos_name2motif_dict,_,_ = file2dict(
		filename=pos_filename,
		max_value_list=merge_max_value_list,
		min_value_list=merge_min_value_list
		)
	############################################
	############################################
	neg_name_list,neg_name2value_dict,neg_name2motif_dict,_,_ = file2dict(
		filename=neg_filename,
		max_value_list=merge_max_value_list,
		min_value_list=merge_min_value_list
		)
	############################################
	merge_value_y_dict = label_dict_merge(pos_name_list=pos_name_list,neg_name_list=neg_name_list)
	merge_value_x_dict = value_dict_merge(pos_name2value_dict=pos_name2value_dict,neg_name2value_dict=neg_name2value_dict)
	merge_value_m_dict = value_dict_merge(pos_name2value_dict=pos_name2motif_dict,neg_name2value_dict=neg_name2motif_dict)
	#############################################
	#############################################
	if split_number == 1:
		splited_neg_name_list = [neg_name_list]
	else:
		splited_neg_name_list = kfold_split_namelist(namelist=neg_name_list,split_number=split_number)
	splited_train_x_list,splited_train_y_list,splited_train_m_list,splited_train_name_list = splited_array_make(
		splited_train_neg_name_list=splited_neg_name_list,
		train_pos_name_list=pos_name_list,
		merge_value_x_dict=merge_value_x_dict,
		merge_value_y_dict=merge_value_y_dict,
		merge_value_m_dict=merge_value_m_dict
		)
	##############################################
	return splited_train_x_list,splited_train_y_list,splited_train_m_list,splited_train_name_list

def dataload_for_site_level_m6A_pred(
	infor_table,
	level_table,
	max_value_filename="/home/ZJRen/pytorch/NGS/data/merge_max_value_list.txt",
	min_value_filename="/home/ZJRen/pytorch/NGS/data/merge_min_value_list.txt"):
	merge_max_value_list = extrem_value_load_form_file(filename=max_value_filename)
	merge_min_value_list = extrem_value_load_form_file(filename=min_value_filename)
	###############################################
	_,merge_name2value_dict,merge_name2motif_dict,_,_ = file2dict(
		filename=infor_table,
		max_value_list=merge_max_value_list,
		min_value_list=merge_min_value_list,
		)
	###############################################
	(m6Aname_to_value_dict,
	m6Aname_to_motif_dict,
	m6Aname_to_level_dict,
	m6Aname_list) = m6A_level_table_read(
		filename=level_table,
		merge_value_dict=merge_name2value_dict,
		merge_motif_dict=merge_name2motif_dict,
		)
	shuffed_m6Aname_list = random_disorganize_using_namelist(namelist=m6Aname_list)
	###############################################
	return shuffed_m6Aname_list,m6Aname_to_value_dict,m6Aname_to_motif_dict,m6Aname_to_level_dict


class Datamake_for_read_level_m6A_pred(Dataset):
	def __init__(self,data_feature_1,data_feature_2,data_target):
		self.len = len(data_feature_1)
		self.features_1 = torch.from_numpy(data_feature_1)
		self.features_2 = torch.from_numpy(data_feature_2)
		self.target = torch.from_numpy(data_target)
	def __getitem__(self,index):
		return self.features_1[index],self.features_2[index],self.target[index]
	def __len__(self):
		return self.len