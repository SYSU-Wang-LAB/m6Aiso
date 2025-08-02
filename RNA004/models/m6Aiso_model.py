import math
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from torch.autograd import Variable


class m6AisoNet(nn.Module):
	def __init__(self,dropout_probility,d_model,number_of_classes):
		super(m6AisoNet,self).__init__()
		#################################################
		self.sequence_conv = nn.Sequential(
			nn.Conv1d(in_channels=5,out_channels=2,kernel_size=9,padding=0,bias=False),
			nn.Tanh()
			)
		#################################################
		self.conv_layer_1 = nn.Sequential(
			nn.Conv2d(in_channels=1,out_channels=d_model,kernel_size=3,padding=1,bias=False),
			nn.BatchNorm2d(d_model),
			nn.Tanh(),
			nn.Dropout(p=dropout_probility,inplace=False)
			)
		self.conv_layer_2 = nn.Sequential(
			nn.Conv2d(in_channels=d_model,out_channels=d_model*4,kernel_size=3,padding=1,bias=False),
			nn.BatchNorm2d(d_model*4),
			nn.Tanh(),
			nn.Dropout(p=dropout_probility,inplace=False)
			)
		self.conv_layer_3 = nn.Sequential(
			nn.Conv2d(in_channels=d_model*4,out_channels=d_model,kernel_size=3,padding=1,bias=False),
			nn.BatchNorm2d(d_model),
			nn.Tanh(),
			nn.Dropout(p=dropout_probility,inplace=False)
			)
		self.global_pooling = nn.AdaptiveMaxPool2d(1)
		self.full_connect = nn.Sequential(
			nn.Linear(d_model,16),
			nn.ReLU(),
			nn.Linear(16,number_of_classes),
			nn.Sigmoid()
			)

	def forward(self,sequence_batch,signal_batch):
		###############
		print(sequence_batch.shape)
		print(signal_batch.shape)

		embeding_sequence_batch = self.sequence_conv(sequence_batch)
		###############
		embeding_sequence_batch = embeding_sequence_batch.view(embeding_sequence_batch.size()[0],1,embeding_sequence_batch.size()[1],embeding_sequence_batch.size()[2])
		signal_batch = signal_batch.view(signal_batch.size()[0],1,signal_batch.size()[1],signal_batch.size()[2])
		x = torch.concat((embeding_sequence_batch,signal_batch),dim=2)
		###############
		res1 = self.conv_layer_1(x)
		res2 = self.conv_layer_2(res1)
		res3 = self.conv_layer_3(res2)
		out = res1 + res3
		out = self.global_pooling(out)
		out = out.reshape(out.size(0),-1)
		y_read = self.full_connect(out)
		return y_read


if __name__=="__main__":
	sequence = torch.zeros((255,5,11))
	signal = torch.zeros((255,3,3))
	model = m6AisoNet(dropout_probility=0.2,d_model=32,number_of_classes=3)
	out = model(
		sequence_batch=sequence,
		signal_batch=signal
		)
	print(out.shape)