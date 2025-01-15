
import torch
from torch.nn import BCELoss
import torch.nn as nn
import torch.backends.cudnn as cudnn

class self_defined_BCEloss(nn.Module):
	def __init__(self):
		super(self_defined_BCEloss,self).__init__()
	def forward(self,y_pred,y_true):
		mean = torch.mean(y_pred)
		loss = torch.pow((mean - y_true),2)
		return loss

class self_defined_siteloss(nn.Module):
	def __init__(self):
		super(self_defined_siteloss,self).__init__()
	def forward(self,y_pred,y_true):
		y_pred = y_pred.reshape(-1)
		y_pred_top5 = y_pred.topk(5).values
		mean = torch.mean(y_pred)
		#################
		#################
		#y_site = 1 - torch.prod((1 - y_pred),axis=-1)
		y_site = 1 - torch.prod((1 - y_pred_top5),axis=-1)
		loss = torch.pow((y_site-y_true),2)
		return loss

class self_binary_cross_entropy_loss(nn.Module):
	def __init__(self):
		super(self_binary_cross_entropy_loss,self).__init__()
	def forward(self,y_pred,y_true):
		loss = BCELoss()(y_pred,y_true)
		return loss

class self_weighted_binary_cross_entropy_loss(nn.Module):
	def __init__(self):
		super(self_weighted_binary_cross_entropy_loss,self).__init__()
	def forward(self,y_pred,y_true):
		labels,counts = torch.unique(y_true,return_counts=True)
		pos_weight,neg_weight = self.label_count(labels=labels,counts=counts)
		sample_weights = torch.where(y_true==0,neg_weight,pos_weight)
		loss = BCELoss(reduction='none')(y_pred,y_true) * sample_weights
		loss = torch.mean(loss)
		return loss
	def label_count(self,labels,counts):
		if labels.size()[0] == 1:
			if labels[0] == 1:
				pos_weight = counts[0]
				neg_weight = torch.tensor(0,device='cuda:0')
			if labels[0] == 0:
				pos_weight = torch.tensor(0,device='cuda:0')
				neg_weight = counts[0]
		else:
			if labels[0] == 1:
				pos_weight = counts[0]
				neg_weight = counts[1]
			if labels[0] == 0:
				pos_weight = counts[1]
				neg_weight = counts[0]
		return pos_weight,neg_weight