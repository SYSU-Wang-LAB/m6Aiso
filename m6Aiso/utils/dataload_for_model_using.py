import torch
import numpy as np
#####################################
from .dataload_utils import extrem_value_load_form_file
from .dataload_utils import file2dict
from .dataload_utils import label_dict_merge
from .dataload_utils import value_dict_merge
from .dataload_utils import kfold_split_namelist
from .dataload_utils import splited_array_make
from .dataload_utils import random_disorganize_train_namelist
from .dataload_utils import random_disorganize_using_namelist
from .dataload_utils import sitename_splited_array_make
from .dataload_utils import using_numpy_array_make
from torch.utils.data import Dataset,DataLoader

def dataload_for_read_level_m6A_using(use_filename,model_type,target_label,max_value_filename,min_value_filename):
	merge_max_value_list = extrem_value_load_form_file(filename=max_value_filename)
	merge_min_value_list = extrem_value_load_form_file(filename=min_value_filename)
	use_name_list,use_name2value_dict,use_name2motif_dict,_,_ = file2dict(
		filename=use_filename,
		max_value_list=merge_max_value_list,
		min_value_list=merge_min_value_list
		)
	out_index_list,use_index2readname_dict,out_label_list = readname_list_to_index(use_name_list=use_name_list,target_label=target_label)
	out_label_array = np.array(out_label_list,dtype = np.float32)
	out_index_array = np.array(out_index_list,dtype = np.float32)
	use_x_value_list,use_m_motif_list = using_numpy_array_make(namelist=use_name_list,merge_value_x_dict=use_name2value_dict,merge_motif_m_dict=use_name2motif_dict)
	############
	if model_type == "2d":
		use_x_value_array = use_x_value_list.reshape(use_x_value_list.shape[0],1,use_x_value_list.shape[1],use_x_value_list.shape[2])
		use_m_motif_array = use_m_motif_list.reshape(use_m_motif_list.shape[0],use_m_motif_list.shape[1],use_m_motif_list.shape[2])
	else:
		use_x_value_array = use_x_value_list
		use_m_motif_array = use_m_motif_list
	############
	using_dataset = Datamake_for_read_level_m6A_using(
		data_feature_1=use_x_value_array,
		data_feature_2=use_m_motif_array,
		data_feature_3=out_index_array,
		data_target=out_label_array
		)
	return using_dataset,use_index2readname_dict

def readname_list_to_index(use_name_list,target_label):
	out_index2readname_dict = {}
	out_index_list = []
	out_label_list = []
	for index_i in range(len(use_name_list)):
		out_index_list.append([index_i])
		out_index2readname_dict[index_i] = use_name_list[index_i]
		out_label_list.append([target_label])
	return out_index_list,out_index2readname_dict,out_label_list
		
class Datamake_for_read_level_m6A_using(Dataset):
	def __init__(self,data_feature_1,data_feature_2,data_feature_3,data_target):
		self.len = len(data_feature_1)
		self.features_1 = torch.from_numpy(data_feature_1)
		self.features_2 = torch.from_numpy(data_feature_2)
		self.features_3 = torch.from_numpy(data_feature_3)
		self.target = torch.from_numpy(data_target)
	def __getitem__(self,index):
		return self.features_1[index],self.features_2[index],self.features_3[index],self.target[index]
	def __len__(self):
		return self.len


