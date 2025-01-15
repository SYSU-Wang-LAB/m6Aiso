import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.backends.cudnn as cudnn
import torch.utils
import numpy as np
from torch.autograd import Variable

def Predict_2d_model_using(using_model_list,using_loader,use_index2sitename_dict,pred_result_save_filename):
	writefile = open(pred_result_save_filename,'a')
	for step,(input_using,motif_using,index_using,target_using) in enumerate(using_loader):
		input_using = input_using.transpose(3,2)
		motif_using = motif_using.transpose(2,1)
		########################
		input_using = Variable(input_using,requires_grad=False).cuda()
		motif_using = Variable(motif_using,requires_grad=False).cuda()
		target_using = Variable(target_using,requires_grad=False).cuda()
		########################
		y_pred_list = []
		for model in using_model_list:
			y_pred_read,_ = model(value=input_using,motif=motif_using)
			y_pred_read = y_pred_read.view(-1)
			y_pred_list.append(y_pred_read.tolist())
		########################
		pred_result_write(
			y_pred_list=y_pred_list,
			index_using=index_using,
			use_index2sitename_dict=use_index2sitename_dict,
			writefile=writefile
			)
	writefile.close()

def Predict_1d_model_using(using_model_list,using_loader,use_index2sitename_dict,pred_result_save_filename):
	writefile = open(pred_result_save_filename,'a')
	for step,(input_using,motif_using,index_using,target_using) in enumerate(using_loader):
		input_using = input_using.transpose(2,1)
		motif_using = motif_using.transpose(2,1)
		########################
		input_using = Variable(input_using,requires_grad=False).cuda()
		motif_using = Variable(motif_using,requires_grad=False).cuda()
		target_using = Variable(target_using,requires_grad=False).cuda()
		########################
		y_pred_list = []
		for model in using_model_list:
			y_pred_read,_ = model(value=input_using,motif=motif_using)
			y_pred_read = y_pred_read.view(-1)
			y_pred_list.append(y_pred_read.tolist())
		########################
		pred_result_write(
			y_pred_list=y_pred_list,
			index_using=index_using,
			use_index2sitename_dict=use_index2sitename_dict,
			writefile=writefile
			)
	writefile.close()

def pred_result_write(y_pred_list,index_using,use_index2sitename_dict,writefile):
	index_list = index_using.tolist()
	for i in range(len(index_list)):
		index_i = index_list[i][0]
		sitename_i = use_index2sitename_dict[index_i]
		y_pred_i = [x[i] for x in y_pred_list]
		y_pred_i_mean = np.mean(y_pred_i)
		str2write = "\t".join([sitename_i,str(y_pred_i_mean)]) + "\n"
		writefile.write(str2write)

def Predict_model_using(model_type,using_model_list,using_loader,use_index2sitename_dict,pred_result_save_filename):
	if model_type == "2d":
		Predict_2d_model_using(
			using_model_list=using_model_list,
			using_loader=using_loader,
			use_index2sitename_dict=use_index2sitename_dict,
			pred_result_save_filename=pred_result_save_filename
			)
	if model_type == "1d":
		Predict_1d_model_using(
			using_model_list=using_model_list,
			using_loader=using_loader,
			use_index2sitename_dict=use_index2sitename_dict,
			pred_result_save_filename=pred_result_save_filename
			)
