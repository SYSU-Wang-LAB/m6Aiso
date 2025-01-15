import os
import sys
import torch
import time
import random
import argparse
import numpy as np
import torch.nn as nn
import torch.nn.functional as F
import torch.utils
import torch.backends.cudnn as cudnn
from torch.utils.data import Dataset,DataLoader
from torch.autograd import Variable
###################
sys.path.insert(0,"../")
sys.path.insert(0,"./m6Aiso/")
from m6Aiso.utils.dataload_for_model_using import dataload_for_read_level_m6A_using_by_step
from m6Aiso.utils.dataload_for_model_using import Datamake_for_read_level_m6A_using
from m6Aiso.blocks.model import AttentionNet
from m6Aiso.blocks.model import Res1dNet
from m6Aiso.blocks.model import Res2dNet
###################
class m6A_predication():
	def __init__(self,model_path,model_type,max_value_filename,min_value_filename):
		self.model_path = model_path
		self.model_type = model_type
		self.max_value_filename = max_value_filename
		self.min_value_filename = min_value_filename
		self.batch_size = 256
		self.colnames_list = []
		self.max_lines = 1000000

	def model_load_for_using(self):
		model = torch.load(self.model_path)
		model.eval()
		############################
		return model

	def data_load_for_using(self,signalfile):
		using_dataset,use_index2readname_dict,colnames_list,endType = dataload_for_read_level_m6A_using_by_step(
			usefile=signalfile,
			colnames_list=self.colnames_list,
			model_type=self.model_type,
			target_label=1,
			max_lines=self.max_lines,
			max_value_filename=self.max_value_filename,
			min_value_filename=self.min_value_filename
			)
		self.colnames_list = colnames_list
		using_data_loader = DataLoader(dataset=using_dataset,batch_size=self.batch_size,shuffle=False,drop_last=False)
		return using_data_loader,use_index2readname_dict,endType

	def pred_result_write(self,y_pred_list,index_using,use_index2sitename_dict,writefile):
		index_list = index_using.tolist()
		for i in range(len(index_list)):
			index_i = index_list[i][0]
			sitename_i = use_index2sitename_dict[index_i]
			y_pred_i = y_pred_list[i]
			str2write = "\t".join([sitename_i,str(y_pred_i)]) + "\n"
			writefile.write(str2write)
	
	def Predict_2d_model_using(self,model,using_loader,use_index2sitename_dict,pred_result_save_filename):
		writefile = open(pred_result_save_filename,'a')
		for step,(input_using,motif_using,index_using,target_using) in enumerate(using_loader):
			input_using = input_using.transpose(3,2)
			motif_using = motif_using.transpose(2,1)
			########################
			input_using = Variable(input_using,requires_grad=False).cuda()
			motif_using = Variable(motif_using,requires_grad=False).cuda()
			target_using = Variable(target_using,requires_grad=False).cuda()
			########################
			y_pred_read,_ = model(value=input_using,motif=motif_using)
			y_pred_read = y_pred_read.view(-1)
			y_pred_list = y_pred_read.tolist()
			self.pred_result_write(
				y_pred_list=y_pred_list,
				index_using=index_using,
				use_index2sitename_dict=use_index2sitename_dict,
				writefile=writefile
				)

	def Predict_1d_model_using(self,model,using_loader,use_index2sitename_dict,pred_result_save_filename):
		writefile = open(pred_result_save_filename,'a')
		for step,(input_using,motif_using,index_using,target_using) in enumerate(using_loader):
			input_using = input_using.transpose(2,1)
			motif_using = motif_using.transpose(2,1)
			########################
			input_using = Variable(input_using,requires_grad=False).cuda()
			motif_using = Variable(motif_using,requires_grad=False).cuda()
			target_using = Variable(target_using,requires_grad=False).cuda()
			########################
			y_pred_read,_ = model(value=input_using,motif=motif_using)
			y_pred_read = y_pred_read.view(-1)
			y_pred_list = y_pred_read.tolist()
			self.pred_result_write(
				y_pred_list=y_pred_list,
				index_using=index_using,
				use_index2sitename_dict=use_index2sitename_dict,
				writefile=writefile
				)


	def model_using(self,using_signal_filename,predict_result_filename):
		####################
		model = self.model_load_for_using()
		####################
		signalfile = open(using_signal_filename,'r')
		endType = True
		while endType == True:
			using_data_loader,use_index2readname_dict,endType_new = self.data_load_for_using(signalfile=signalfile)
			####################
			if self.model_type == "1d":
				self.Predict_1d_model_using(
					model=model,
					using_loader=using_data_loader,
					use_index2sitename_dict=use_index2readname_dict,
					pred_result_save_filename=predict_result_filename
					)
			if self.model_type == "2d":
				self.Predict_2d_model_using(
					model=model,
					using_loader=using_data_loader,
					use_index2sitename_dict=use_index2readname_dict,
					pred_result_save_filename=predict_result_filename
					)
			endType = endType_new


def args_make(parser):
	print(sys.argv[0])
	print(os.getcwd())
	parser.add_argument('--model_path',help='m6A predict model name for using',default="./m6Aiso/module/semi_model_7mer.times_over.epoch_0.pt")
	parser.add_argument('--model_type',help="1d or 2d",default="2d")
	parser.add_argument('--using_signal_filename',required=True,help="Data perpare generated signal.tsv")
	parser.add_argument('--predict_result_filename', required=True, help="model predicated molecular m6A level filename")
	parser.add_argument('--max_value_filename',help="max value of current signal mean,std and dwelltime",default="./m6Aiso/data/merge_max_value_list.txt")
	parser.add_argument('--min_value_filename',help="min value of current signal mean,std and dwelltime",default="./m6Aiso/data/merge_min_value_list.txt")
	return parser

def main(args):

	np.random.seed(123)
	torch.manual_seed(123)
	torch.cuda.manual_seed(123)
	torch.cuda.set_device(0)
	cudnn.benchmark = True
	cudnn.enabled=True
	#################
	#################
	PRED = m6A_predication(
		model_path=args.model_path,
		model_type=args.model_type,
		max_value_filename=args.max_value_filename,
		min_value_filename=args.min_value_filename
		)
	PRED.model_using(
		using_signal_filename=args.using_signal_filename,
		predict_result_filename=args.predict_result_filename
		)
