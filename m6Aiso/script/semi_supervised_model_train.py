import torch
import time
import sys
import os
import random
import numpy as np
import torch.nn as nn
import torch.nn.functional as F
import torch.utils
import torch.backends.cudnn as cudnn
from torch.utils.data import Dataset,DataLoader
from torch.autograd import Variable
###############################
sys.path.insert(0,"../../")
sys.path.insert(0,"../../m6Aiso/")
sys.path.insert(0,"../")
sys.path.insert(0,"./m6Aiso/")
###############################
from m6Aiso.utils.dataload_for_model_train import dataload_for_read_level_m6A_pred
from m6Aiso.utils.dataload_for_model_train import Datamake_for_read_level_m6A_pred
from m6Aiso.utils.dataload_for_model_using import dataload_for_read_level_m6A_using
from m6Aiso.utils.dataload_for_model_using import Datamake_for_read_level_m6A_using
from m6Aiso.utils.get_lossfunction import self_binary_cross_entropy_loss
###################################################################
from m6Aiso.blocks.model import AttentionNet
from m6Aiso.blocks.model import Res1dNet
from m6Aiso.blocks.model import Res2dNet
###############################
from m6Aiso.semi.semi_model_train import Predict_model_train
from m6Aiso.semi.semi_model_using import Predict_model_using
from m6Aiso.semi.semi_positive_data_clean import positive_data_clean_and_FPR_count
from m6Aiso.semi.semi_negative_data_make import motif_balance_neg_data_make
###############################


class semi_model_train():
	def __init__(
		self,
		model_name,
		orginal_pos_filename,
		orginal_neg_filename,
		balance_neg_filename,
		cleaned_pos_filename,
		pred_pos_save_filename,
		max_value_filename,
		min_value_filename,
		model_save_filename,
		negative_splited_number=5,
		negative_over_sample_number=12,
		min_negative_sample_number=100,
		positive_cutoff=0.5,
		percentage_cutoff=0.95,
		learning_rate = 0.001
		):
		self.model_name = model_name
		self.orginal_pos_filename = orginal_pos_filename
		self.orginal_neg_filename = orginal_neg_filename
		self.balance_neg_filename = balance_neg_filename
		self.cleaned_pos_filename = cleaned_pos_filename
		####################################################################
		self.max_value_filename = max_value_filename
		self.min_value_filename = min_value_filename
		self.pred_pos_save_filename = pred_pos_save_filename
		self.model_save_filename = model_save_filename
		####################################################################
		self.negative_splited_number = negative_splited_number
		self.negative_over_sample_number = negative_over_sample_number
		self.min_negative_sample_number = min_negative_sample_number
		self.positive_cutoff = positive_cutoff
		self.percentage_cutoff = percentage_cutoff
		self.learning_rate = learning_rate
		self.epochs = 5
		self.batch_size = 512
		####################################################################

	def model_init(self):
		if self.model_name == "AttentionNet":
			model = AttentionNet()
			model = model.cuda()
			self.model_type = "1d"
		elif self.model_name == "Res1dNet":
			model = Res1dNet()
			model = model.cuda()
			self.model_type = "1d"
		elif self.model_name == "Res2dNet":
			model = Res2dNet()
			model = model.cuda()
			self.model_type = "2d"
		else:
			print(self.model_name)
			raise Exception("model make error")
		return model

	def optimizer_init(self,model):
		optimizer = torch.optim.Adam(model.parameters(),self.learning_rate)
		return optimizer

	def criterion_init(self,model):
		criterion = self_binary_cross_entropy_loss()
		criterion = criterion.cuda()
		return criterion

	def model_save(self,model,model_save_path):
		torch.save(model,model_save_path)

	def model_weight_save(self,model,model_save_path):
		torch.save(model.state_dict(),model_save_path)

	def train_data_loader_make(self,pos_filename,neg_filename,negative_splited_number=None):
		if negative_splited_number == None:
			negative_splited_number = self.negative_splited_number
		out_data_loader_list = []
		#################
		splited_train_x_list,splited_train_y_list,splited_train_m_list,_ = dataload_for_read_level_m6A_pred(
			pos_filename=pos_filename,
			neg_filename=neg_filename,
			max_value_filename=self.max_value_filename,
			min_value_filename=self.min_value_filename,
			split_number=negative_splited_number
		)
		#################
		for k in range(negative_splited_number):
			k_train_x_list = splited_train_x_list[k]
			k_train_m_list = splited_train_m_list[k]
			k_train_y_list = splited_train_y_list[k]
			if self.model_type == "2d":
				k_train_x_array = k_train_x_list.reshape(k_train_x_list.shape[0],1,k_train_x_list.shape[1],k_train_x_list.shape[2])
				k_train_m_array = k_train_m_list.reshape(k_train_m_list.shape[0],k_train_m_list.shape[1],k_train_m_list.shape[2])
			if self.model_type == "1d":
				k_train_x_array = k_train_x_list
				k_train_m_array = k_train_m_list
			k_train_y_array = k_train_y_list
			k_train_dataset = Datamake_for_read_level_m6A_pred(data_feature_1=k_train_x_array,data_feature_2=k_train_m_array,data_target=k_train_y_array)
			k_train_loader = DataLoader(dataset=k_train_dataset,batch_size=self.batch_size,shuffle=True,drop_last=False)
			out_data_loader_list.append(k_train_loader)
		return out_data_loader_list

	def using_data_loader_make(self,use_filename,target_label):
		using_dataset,use_index2readname_dict = dataload_for_read_level_m6A_using(
			use_filename = use_filename,
			model_type = self.model_type,
			max_value_filename = self.max_value_filename,
			min_value_filename = self.min_value_filename,
			target_label = target_label
			)
		using_data_loader = DataLoader(dataset=using_dataset,batch_size=self.batch_size,shuffle=True,drop_last=False)
		return using_data_loader,use_index2readname_dict

	def model_train(self,pos_filename,neg_filename,times):
		model = self.model_init()
		optimizer = self.optimizer_init(model=model)
		criterion = self.criterion_init(model=model)
		data_loader_list = self.train_data_loader_make(
			pos_filename=pos_filename,
			neg_filename=neg_filename
			)
		for k in range(self.negative_splited_number):
			k_train_loader = data_loader_list[k]
			model = Predict_model_train(
				model_type = self.model_type,
				model = model,
				epochs = self.epochs,
				k_train_loader = k_train_loader,
				optimizer = optimizer,
				criterion = criterion
				)
			model_save_path_k = self.model_save_filename % (str(times),str(k))
			self.model_save(model=model,model_save_path=model_save_path_k)

	def output_model_train(self,pos_filename,neg_filename,neg_select_times):
		model = self.model_init()
		optimizer = self.optimizer_init(model=model)
		criterion = self.criterion_init(model=model)
		###########################################################
		data_loader_list = self.train_data_loader_make(
			pos_filename=pos_filename,
			neg_filename=neg_filename,
			negative_splited_number=neg_select_times
			)
		###########################################################
		####################################
		for k in range(neg_select_times):
			k_train_loader = data_loader_list[k]
			model = Predict_model_train(
				model_type = self.model_type,
				model = model,
				epochs = self.epochs,
				k_train_loader = k_train_loader,
				optimizer = optimizer,
				criterion = criterion
				)
			model_save_path_k = self.model_save_filename % ("over",str(k))
			self.model_save(model=model,model_save_path=model_save_path_k)
		####################################

	def model_load_for_using(self,times):
		out_model_list = []
		for k in range(self.negative_splited_number):
			model_save_path_k = self.model_save_filename % (str(times),str(k))
			model_loaded_k = torch.load(model_save_path_k)
			model_loaded_k.eval()
			out_model_list.append(model_loaded_k)
		return out_model_list


	def model_using(self,pos_filename,pred_result_save_filename,times):
		#pred_result_save_filename = self.pred_pos_save_filename % (str(times))
		#pred_result_save_filename = self.pred_result_save_path + "semi_m6A_model_predicate_result_%s.txt" % (str(times))
		model_list = self.model_load_for_using(times=times)
		using_data_loader,use_index2readname_dict = self.using_data_loader_make(use_filename=pos_filename,target_label=1)
		###################
		Predict_model_using(
			model_type = self.model_type,
			using_model_list=model_list,
			using_loader=using_data_loader,
			use_index2sitename_dict=use_index2readname_dict,
			pred_result_save_filename=pred_result_save_filename
			)
		###################

	def positive_data_clean(self,pred_result_save_filename,uncleaned_positive_filename,cleaned_positive_filename,times):
		FPR = positive_data_clean_and_FPR_count(
			pred_result_filename=pred_result_save_filename,
			uncleaned_positive_filename=uncleaned_positive_filename,
			cleaned_positive_filename=cleaned_positive_filename,
			cutoff=self.positive_cutoff,
			cutoff_perventage=self.percentage_cutoff
			)
		return FPR

	def negative_data_make(self,cleaned_pos_filename,cleaned_neg_filename,motif_balanced_neg_filename,times):
		########################################
		########################################
		motif_balance_neg_data_make(
			apobe_neg_filename=cleaned_neg_filename,
			apobe_pos_filename=cleaned_pos_filename,
			balan_neg_filename=motif_balanced_neg_filename,
			min_negative_sample_number=self.min_negative_sample_number,
			negative_over_sample_number=self.negative_over_sample_number
			)

	def step_semi_predict_model_train(self,max_train_times,min_false_positive_rate):
		np.random.seed(123)
		torch.manual_seed(123)
		torch.cuda.manual_seed(123)
		torch.cuda.set_device(0)
		cudnn.benchmark = True
		cudnn.enabled=True
		#####################################
		#####################################
		times = 0
		while times<=max_train_times:
			###########################################################
			# motif balanced negative data make
			###########################################################
			if times == 0:
				cleaned_neg_filename=self.orginal_neg_filename
				cleaned_pos_filename=self.orginal_pos_filename
				motif_balanced_neg_filename = self.balance_neg_filename % (str(times))
			else:
				cleaned_neg_filename=self.orginal_neg_filename
				cleaned_pos_filename=self.cleaned_pos_filename % (str(times))
				motif_balanced_neg_filename = self.balance_neg_filename % (str(times))
			self.negative_data_make(
				cleaned_pos_filename=cleaned_pos_filename,
				cleaned_neg_filename=cleaned_neg_filename,
				motif_balanced_neg_filename=motif_balanced_neg_filename,
				times=times
				)
			print("motif balanced negative data make at times %s"%(times))
			###########################################################
			# semi model training
			###########################################################
			if times == 0:
				tmp_pos_filename = self.orginal_pos_filename
				tmp_neg_filename = self.balance_neg_filename % (str(times))
			else:
				tmp_pos_filename = self.cleaned_pos_filename % (str(times))
				tmp_neg_filename = self.balance_neg_filename % (str(times))
			self.model_train(
				pos_filename=tmp_pos_filename,
				neg_filename=tmp_neg_filename,
				times=times
				)
			print("semi model training at times %s"%(times))
			############################################################
			# semi model using
			############################################################
			pred_pos_result_save_filename = self.pred_pos_save_filename % (str(times))
			self.model_using(
				pos_filename=tmp_pos_filename,
				pred_result_save_filename=pred_pos_result_save_filename,
				times=times
				)
			print("semi model using at times %s"%(times))
			############################################################
			# cleaned positive data make
			############################################################
			if times == 0:
				uncleaned_positive_filename = self.orginal_pos_filename
				cleaned_positive_filename = self.cleaned_pos_filename % (str(times+1))
			else:
				uncleaned_positive_filename = self.cleaned_pos_filename % (str(times))
				cleaned_positive_filename = self.cleaned_pos_filename % (str(times+1))
			FPR = self.positive_data_clean(
				pred_result_save_filename=pred_pos_result_save_filename,
				uncleaned_positive_filename=uncleaned_positive_filename,
				cleaned_positive_filename=cleaned_positive_filename,
				times=times
				)
			print("cleaned positive data make at times %s,FPR:%s"%(times,FPR))
			if FPR < min_false_positive_rate or times == max_train_times:
				cleaned_neg_filename = self.orginal_neg_filename
				cleaned_pos_filename = cleaned_positive_filename
				motif_balanced_neg_filename = self.balance_neg_filename % ("over")
				motif_balance_neg_data_make(
					apobe_neg_filename=cleaned_neg_filename,
					apobe_pos_filename=cleaned_pos_filename,
					balan_neg_filename=motif_balanced_neg_filename,
					min_negative_sample_number=100,
					negative_over_sample_number=10
				)
				self.output_model_train(
					pos_filename=cleaned_pos_filename,
					neg_filename=motif_balanced_neg_filename,
					neg_select_times=10
					)
				break
			#self.percentage_cutoff += 0.0015
			#if self.percentage_cutoff >= 1:
			#	self.percentage_cutoff = 0.99
			times += 1


def args_make(parser):
	import argparse
	parser.add_argument('--model_name', required=True, help="model was selected form training: AttentionNet,Res1dNet,Res2dNet")
	parser.add_argument('--orginal_pos_filename', required=True, help="m6A positive dataset generated form data perpare")
	parser.add_argument('--orginal_neg_filename', required=True, help="m6A negative dataset generated form data perpare")
	parser.add_argument('--balance_neg_filename', required=True, help="m6A positive dataset generated form data perpare")
	parser.add_argument('--cleaned_pos_filename', required=True, help="m6A negative dataset generated form data perpare")
	parser.add_argument('--pred_pos_save_filename', required=True, help="m6A negative dataset generated form data perpare")
	parser.add_argument('--model_save_filename', required=True, help="m6A negative dataset generated form data perpare")
	parser.add_argument('--max_value_filename', required=True, help="max value of current signal mean,std and dwelltime")
	parser.add_argument('--min_value_filename', required=True, help="min value of current signal mean,std and dwelltime")
	parser.add_argument('--out_dir', required=True, help="output directory. /home/ZJRen/pytorch/m6Aiso/m6Aiso/data/")
	parser.add_argument('--learning_rate', required=False,default=0.001, help="learning_rate")
	parser.add_argument('--percentage_cutoff', required=False,default=0.95, help="learning_rate")
	return parser

def main(args):
	###################################
	if not os.path.exists(args.out_dir):
		os.makedirs(args.out_dir)
	###################################
	out_dir = args.out_dir
	###################################
	model_name = args.model_name
	balance_neg_filename = out_dir + args.balance_neg_filename                #balance_m6A_model_train_negative.times_%s.tsv
	cleaned_pos_filename = out_dir + args.cleaned_pos_filename                #cleaned_m6A_model_train_positive.times_%s.tsv
	pred_pos_save_filename = out_dir + args.pred_pos_save_filename        #pos_m6A_predicate_result.times_%s.tsv
	model_save_filename = out_dir + args.model_save_filename              #semi_m6A_predicate_model.times_%s.number_%s.pt
	###################################
	orginal_pos_filename = args.orginal_pos_filename
	orginal_neg_filename = args.orginal_neg_filename
	max_value_filename = args.max_value_filename
	min_value_filename = args.min_value_filename
	#######################################
	learning_rate = float(args.learning_rate)
	percentage_cutoff = float(args.percentage_cutoff)
	#######################################
	np.random.seed(123)
	torch.manual_seed(123)
	torch.cuda.manual_seed(123)
	torch.cuda.set_device(0)
	cudnn.benchmark = True
	cudnn.enabled=True
	#######################################
	SEMI = semi_model_train(
		model_name=model_name,
		orginal_pos_filename=orginal_pos_filename,
		orginal_neg_filename=orginal_neg_filename,
		balance_neg_filename=balance_neg_filename,
		cleaned_pos_filename=cleaned_pos_filename,
		pred_pos_save_filename=pred_pos_save_filename,
		max_value_filename=max_value_filename,
		min_value_filename=min_value_filename,
		model_save_filename=model_save_filename,
		negative_splited_number=2,
		negative_over_sample_number=5,
		min_negative_sample_number=10,
		positive_cutoff=0.5,
		percentage_cutoff=percentage_cutoff,
		learning_rate = learning_rate
		)
	SEMI.step_semi_predict_model_train(
		max_train_times=50,
		min_false_positive_rate=0.05
		)
	#######################################
	#######################################

if __name__=="__main__":
	import argparse
	parser = argparse.ArgumentParser(description='model training')
	subparsers = parser.add_subparsers(help="semi_model_train,m6A_predication")
	parser_sub_1 = subparsers.add_parser('semi_supervised_model_train',help="semi supervised model training.")
	parser_sub_1 = args_make(parser=parser_sub_1)
	args = args_make(parser_sub_1)
	#parser_sub_1 = subparsers.add_parser('semi_supervised_model_train',help="semi supervised model training.")
	#parser_sub_1 = semi_supervised_model_train.args_make(parser_sub_1)
	#parser_sub_1.set_defaults(func=semi_supervised_model_train.main)
	main(args)