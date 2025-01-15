import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from m6Aiso.blocks.attention_blocks import ScaledDotProduct_self_attention_block
from m6Aiso.blocks.transformer_blocks import Positionnal_encoding
from m6Aiso.blocks.transformer_blocks import PositionwiseFeedForward
from m6Aiso.blocks.transformer_blocks import MultiHead_attention


class AttentionNet(nn.Module):
	def __init__(self,d_model=32):
		super(AttentionNet,self).__init__()
		self.conv_motif = nn.Sequential(
			nn.Conv1d(in_channels=4,out_channels=2,kernel_size=5,padding=0,bias=False),
			nn.Tanh(),
			)
		self.conv_layer = nn.Sequential(
			nn.Conv1d(in_channels=5,out_channels=d_model,kernel_size=3,padding=1,bias=False),
			nn.BatchNorm1d(d_model),
			nn.ReLU()
			)
		self.attention_layer = MultiHead_attention(d_model=d_model,header_num=2,dropout_ratio=0.2)
		self.ffn = PositionwiseFeedForward(d_model=d_model,hidden_size=8,drop_ratio=0.2)
		self.global_pooling = nn.MaxPool1d(kernel_size=3,stride=1,padding=1)
		self.full_connect = nn.Sequential(
			nn.Linear(d_model*3,16),
			nn.ReLU(),
			nn.Linear(16,1),
			nn.Sigmoid()
			)

	def forward(self,value,motif):
		mot1 = self.conv_motif(motif)
		x = torch.concat((mot1,value),dim=1)
		out = self.conv_layer(x)
		out = out.transpose(2,1)
		att = self.attention_layer(q=out,k=out,v=out,mask=None)
		att = att + out
		out = self.ffn(att)
		out = out.transpose(2,1)
		out = self.global_pooling(out)
		out = out.reshape(out.size(0),-1)
		y_read = self.full_connect(out)
		y_site = self.site_m6A_pred(y_read)
		return y_read,y_site

	def site_m6A_pred(self,y_read):
		y_pred = y_read.reshape(-1)
		y_pred_top5 = y_pred.topk(5).values
		y_site = 1 - torch.prod((1 - y_pred_top5),axis=-1)
		return y_site


class Res2dNet(nn.Module):
	def __init__(self,d_model=32):
		super(Res2dNet,self).__init__()

		self.conv_motif = nn.Sequential(
			nn.Conv1d(in_channels=4,out_channels=2,kernel_size=5,padding=0,bias=False),
			nn.Tanh(),
			)		

		self.conv_layer_1 = nn.Sequential(
			nn.Conv2d(in_channels=1,out_channels=d_model,kernel_size=3,padding=1,bias=False),
			nn.BatchNorm2d(d_model),
			nn.Tanh(),
			nn.Dropout(p=0.2,inplace=False)
			)
		self.conv_layer_2 = nn.Sequential(
			nn.Conv2d(in_channels=d_model,out_channels=d_model*4,kernel_size=3,padding=1,bias=False),
			nn.BatchNorm2d(d_model*4),
			nn.Tanh(),
			nn.Dropout(p=0.2,inplace=False)
			)
		self.conv_layer_3 = nn.Sequential(
			nn.Conv2d(in_channels=d_model*4,out_channels=d_model,kernel_size=3,padding=1,bias=False),
			nn.BatchNorm2d(d_model),
			nn.Tanh(),
			nn.Dropout(p=0.2,inplace=False)
			)
		self.global_pooling = nn.AdaptiveMaxPool2d(1)
		self.full_connect = nn.Sequential(
			nn.Linear(d_model,16),
			nn.ReLU(),
			nn.Linear(16,1),
			nn.Sigmoid()
			)
	def forward(self,value,motif):
		mot1 = self.conv_motif(motif)
		###############
		mot2 = mot1.view(mot1.size()[0],1,mot1.size()[1],mot1.size()[2])
		x = torch.concat((mot2,value),dim=2)
		###############
		res1 = self.conv_layer_1(x)
		res2 = self.conv_layer_2(res1)
		res3 = self.conv_layer_3(res2)
		out = res1 + res3
		out = self.global_pooling(out)
		out = out.reshape(out.size(0),-1)
		y_read = self.full_connect(out)
		y_site = self.site_m6A_pred(y_read)
		return y_read,y_site

	def site_m6A_pred(self,y_read):
		y_pred = y_read.reshape(-1)
		y_pred_top5 = y_pred.topk(5).values
		y_site = 1 - torch.prod((1 - y_pred_top5),axis=-1)
		return y_site


class Res1dNet(nn.Module):
	def __init__(self,d_model=32):
		super(Res1dNet,self).__init__()
		self.conv_motif = nn.Sequential(
			nn.Conv1d(in_channels=4,out_channels=2,kernel_size=5,padding=0,bias=False),
			nn.Tanh(),
			)
		self.conv_layer_1 = nn.Sequential(
			nn.Conv1d(in_channels=5,out_channels=d_model,kernel_size=3,padding=1,bias=False),
			nn.BatchNorm1d(d_model),
			nn.ReLU(),
			nn.Dropout(p=0.2,inplace=False)
			)
		self.conv_layer_2 = nn.Sequential(
			nn.Conv1d(in_channels=d_model,out_channels=4*d_model,kernel_size=3,padding=1,bias=False),
			nn.BatchNorm1d(4*d_model),
			nn.ReLU(),
			nn.Dropout(p=0.2,inplace=False)
			)
		self.conv_layer_3 = nn.Sequential(
			nn.Conv1d(in_channels=4*d_model,out_channels=d_model,kernel_size=3,padding=1,bias=False),
			nn.BatchNorm1d(d_model),
			nn.ReLU(),
			nn.Dropout(p=0.2,inplace=False)
			)
		self.global_pooling = nn.MaxPool1d(kernel_size=3,stride=1,padding=1)
		self.full_connect = nn.Sequential(
			nn.Linear(d_model*3,16),
			nn.ReLU(),
			nn.Linear(16,1),
			nn.Sigmoid()
			)

	def forward(self,value,motif):
		mot1 = self.conv_motif(motif)
		x = torch.concat((mot1,value),dim=1)
		res1 = self.conv_layer_1(x)
		res2 = self.conv_layer_2(res1)
		res3 = self.conv_layer_3(res2)
		out = res1 + res3
		out = self.global_pooling(out)
		out = out.reshape(out.size(0),-1)
		y_read = self.full_connect(out)
		y_site = self.site_m6A_pred(y_read)
		return y_read,y_site

	def site_m6A_pred(self,y_read):
		y_pred = y_read.reshape(-1)
		y_pred_top5 = y_pred.topk(5).values
		y_site = 1 - torch.prod((1 - y_pred_top5),axis=-1)
		return y_site

