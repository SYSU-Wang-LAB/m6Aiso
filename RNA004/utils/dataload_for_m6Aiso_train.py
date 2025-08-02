import torch
import time
from torch.utils.data import Dataset,DataLoader
from .dataload_utils_for_m6Aiso import signal_data_read,sitename_to_label_dictmake,sitename_id_make

def Train_data_load(pos_signal_filename,neg_signal_filename,label_filename):
	pos_name2sequence_dict,pos_name2signal_dict,_ = signal_data_read(filename=pos_signal_filename,using_start=0,using_end=3)

	#print(pos_name2sequence_dict["ATGGACCGG_ENST00000335791.10_e840684e-b5aa-42a9-a942-88a7102401c7_759"])
	#print(pos_name2signal_dict["ATGGACCGG_ENST00000335791.10_e840684e-b5aa-42a9-a942-88a7102401c7_759"])

	neg_name2sequence_dict,neg_name2signal_dict,_ = signal_data_read(filename=neg_signal_filename,using_start=0,using_end=3)
	meg_sitename2label_dict,meg_sitename_list = sitename_to_label_dictmake(filename=label_filename)
	meg_id2sitename_dict,meg_id_list = sitename_id_make(sitename_list=meg_sitename_list)
	################################################################################
	meg_sitename2sequence_dict = value_dict_merge(sitename_list=meg_sitename_list,pos_name2value_dict=pos_name2sequence_dict,neg_name2value_dict=neg_name2sequence_dict)
	meg_sitename2signal_dict = value_dict_merge(sitename_list=meg_sitename_list,pos_name2value_dict=pos_name2signal_dict,neg_name2value_dict=neg_name2signal_dict)
	################################################################################
	meg_sequence_list = value_list_make(sitename_list=meg_sitename_list,name2value_dict=meg_sitename2sequence_dict)
	meg_signal_list = value_list_make(sitename_list=meg_sitename_list,name2value_dict=meg_sitename2signal_dict)
	meg_label_list = value_list_make(sitename_list=meg_sitename_list,name2value_dict=meg_sitename2label_dict)
	dataloader = Datamake_for_m6Aiso_model_training(
		data_feature_1 = meg_sequence_list,
		data_feature_2 = meg_signal_list,
		data_feature_3 = meg_id_list,
		data_target = meg_label_list
		)
	return dataloader,meg_id2sitename_dict

def value_dict_merge(sitename_list,pos_name2value_dict,neg_name2value_dict):
	out_sitename2value_dict = {}
	for sitename in sitename_list:
		try:
			value = pos_name2value_dict[sitename]
			#out_sitename2value_dict[sitename] = value
		except:
			pass
		try:
			value = neg_name2value_dict[sitename]
			#out_sitename2value_dict[sitename] = value
		except:
			pass
		out_sitename2value_dict[sitename] = value
	return out_sitename2value_dict

def value_list_make(sitename_list,name2value_dict):
	out_value_list = []
	for sitename in sitename_list:
		value = name2value_dict[sitename]
		out_value_list.append(value)
	return out_value_list

class Datamake_for_m6Aiso_model_training(Dataset):
	def __init__(self,data_feature_1,data_feature_2,data_feature_3,data_target):
		self.len = len(data_feature_1)
		self.features_1 = torch.tensor(data_feature_1,dtype=torch.float32)
		self.features_2 = torch.tensor(data_feature_2,dtype=torch.float32)
		self.features_3 = torch.tensor(data_feature_3,dtype=torch.int)
		self.target = torch.tensor(data_target,dtype=torch.float32)
	def __getitem__(self,index):
		return self.features_1[index],self.features_2[index],self.features_3[index],self.target[index]
	def __len__(self):
		return self.len


if __name__ == "__main__":
	import sys
	dataloader = Train_data_load(
		pos_signal_filename="/home/ZJRen/pytorch/semi_suprevise_model/bin/test/test_positive.flanking.tsv.gz",
		neg_signal_filename="/home/ZJRen/pytorch/semi_suprevise_model/bin/test/test_negative.flanking.tsv.gz",
		label_filename="/home/ZJRen/pytorch/semi_suprevise_model/bin/test/test_label.txt.gz"
		)
