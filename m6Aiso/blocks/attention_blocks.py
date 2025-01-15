######################
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
######################

##############################################
##############################################
#### SE attention
class SE_attention_block(nn.Module):
	def __init__(self,channel,reduction,kernel_size=3):
		super(SE_attention_block,self).__init__()
		self.avg_pool = nn.AdaptiveAvgPool2d(1)
		self.fc = nn.Sequential(
			nn.Linear(channel,channel//reduction,bias=False),
			nn.ReLU(inplace=True),
			nn.Linear(channel//reduction,channel,bias=False),
			nn.Sigmoid()
			)
	def forward(self,x):
		b,c,_,_ = x.size()
		y = self.avg_pool(x)
		y = y.view(b,c)
		y = self.fc(y)
		y = y.view(b,c,1,1)
		return x * y.expand_as(x)

##############################################
##############################################
#### CNAM attention
class Channel_attention(nn.Module):
	def __init__(self,channel,reduction=4):
		super(Channel_attention,self).__init__()
		self.maxpool = nn.AdaptiveMaxPool2d(1)
		self.avgpool = nn.AdaptiveAvgPool2d(1)
		self.Squeeze = nn.Sequential(
			nn.Conv2d(channel,channel//reduction,1,bias=False),
			nn.ReLU(),
			nn.Conv2d(channel//reduction,channel,1,bias=False)
			)
		self.sigmoid = nn.Sigmoid()

	def forward(self,x):
		max_result = self.maxpool(x)
		avg_result = self.avgpool(x)
		max_out = self.Squeeze(max_result)
		avg_out = self.Squeeze(avg_result)
		output = self.sigmoid(max_out+avg_out)
		return output

class Spatial_attention(nn.Module):
	def __init__(self,kernel_size=7):
		super(Spatial_attention,self).__init__()
		self.conv = nn.Conv2d(2,1,kernel_size=kernel_size,padding=kernel_size//2)
		self.sigmoid = nn.Sigmoid()
	def forward(self,x):
		max_result = torch.max(x,dim=1,keepdim=True)
		print(max_result.values.size())
		avg_result = torch.mean(x,dim=1,keepdim=True)
		print(avg_result.size())
		result = torch.cat([max_result.values,avg_result],1)
		output = self.conv(result)
		output = self.sigmoid(output)
		return output


class CBAM_attention_block(nn.Module):
	def __init__(self,channel,reduction,kernel_size=7):
		super(CBAM_attention_block,self).__init__()
		self.channel_attention = Channel_attention(channel=channel,reduction=reduction)
		self.spatial_attention = Spatial_attention(kernel_size=kernel_size)

	def forward(self,x):
		b,c,_,_ = x.size()
		residual = x
		out = x*self.channel_attention(x)
		out = out*self.spatial_attention(out)
		return out + residual

##############################################
##############################################
#### ECA attention
class ECA_attention_block(nn.Module):
	def __init__(self,kernel_size=3):
		super(ECA_attention_block,self).__init__()
		self.gap = nn.AdaptiveAvgPool2d(1)
		self.conv = nn.Conv1d(1,1,kernel_size=kernel_size)
		self.sigmoid = nn.Sigmoid()
	def forward(self,x):
		out = self.gap(x)
		out = out.seqeeze(-1)
		out = out.permute(0,2,1)
		out = self.conv(out)
		out = self.sigmoid(out)
		y = out.permute(0,2,1).unsqueeze(-1)
		return x*y.expand_as(x)


##############################################
##############################################
#### ScaledDotProduct attention

class ScaledDotProduct_self_attention_block(nn.Module):
	def __init__(self,d_model,d_key,d_value,h_number,dropout_ratio=0.1):
		super(ScaledDotProduct_self_attention_block,self).__init__()
		self.W_q = nn.Linear(d_model,h_number * d_key)
		self.W_k = nn.Linear(d_model,h_number * d_key)
		self.W_v = nn.Linear(d_model,h_number * d_value)
		self.W_o = nn.Linear(h_number * d_value,d_model)
		self.dropout = nn.Dropout(dropout_ratio)

		self.d_model = d_model
		self.d_key = d_key
		self.d_value = d_value
		self.h_number = h_number

	def forward(self,queries,keys,values,attention_mask=None,attention_weights=None):
		print(queries.size())
		bs,nq = queries.size()[:2]
		nk = keys.size()[1]
		q = self.W_q(queries).view(bs,nq,self.h_number,self.d_key).permute(0, 2, 1, 3)
		k = self.W_k(keys).view(bs,nk,self.h_number,self.d_key).permute(0, 2, 3, 1)
		v = self.W_v(values).view(bs,nk,self.h_number,self.d_value).permute(0, 2, 1, 3)
		#######################
		# att = quare * keys
		att = torch.matmul(q,k)
		att = att / np.sqrt(self.d_key)
		#######################
		if attention_weights is not None:
			att = att * attention_weights
		if attention_mask is not None:
			att = att.masked_fill(attention_mask,-np.inf)
		att = torch.softmax(att,-1)
		att = self.dropout(att)
		out = torch.matmul(att,v).permute(0,2,1,3)
		out = out.contiguous().view(bs,nq,self.h_number*self.d_value)
		y = self.W_o(out)
		return y


##############################################
##############################################
#### CNAM attention

class Multihead_self_attention_block(nn.Module):
	def __init__(self,d_model,d_key,d_value,h_number,header_num,H,W,ApplyTransform,dropout_ratio=0.1):
		super(Multihead_self_attention_block,self).__init__()
		self.W_q = nn.Linear(d_model,h_number * d_key)
		self.W_k = nn.Linear(d_model,h_number * d_key)
		self.W_v = nn.Linear(d_model,h_number * d_value)
		self.W_o = nn.Linear(h_number * d_value,d_model)
		self.dropout = nn.Dropout(dropout_ratio)
		self.d_model = d_model
		self.d_key = d_key
		self.d_value = d_value
		self.h_number = h_number
		self.header_num = header_num
		self.H = H
		self.W = W
		self.ApplyTransform = ApplyTransform
		self.attention_block = nn.Sequential()
		self.attention_block.add_module('conv_1',nn.Conv2d(d_model,d_model,kernel_size=header_num+1,stride=header_num,padding=header_num//2,groups=d_model))
		self.attention_block.add_module('liner_1',nn.LayerNorm(d_model))
		self.transform_block = nn.Sequential()
		self.transform_block.add_module("conv_1",nn.Conv2d(h_number,h_number,kernel_size=1,stride=1))
		self.transform_block.add_module("softmax",nn.Softmax(-1))
		self.transform_block.add_module("normal",nn.InstanceNorm2d(h_number))


	def transpose_for_multihead_classical(self,inputX):
		X = inputX.reshape(inputX.size()[0],inputX.size()[1],self.header_num,-1)
		X = X.permute(0,2,1,3)
		Y = X.reshape(-1,X.size()[2],X.size()[3])
		return Y


	def forward(self,queries,keys,values,attention_mask=None,attention_weights=None):
		print(queries.size())
		bs,nq,channel = queries.size()[:3]
		nk = keys.size()[1]
		q = self.W_q(queries).view(bs,nq,self.h_number,self.d_key).permute(0, 2, 1, 3)
		############################################################
		print(queries.size())
		x = queries.permute(0,2,1).view(bs,channel,self.H,self.W)
		print(x.size())
		x = self.attention_block.conv_1(x)
		print(x.size())
		x = x.contiguous().view(bs,channel,-1).permute(0,2,1)
		print(x.size())
		x = self.attention_block.liner_1(x)
		print(x.size())
		############################################################
		k = self.W_k(x).view(bs,-1,self.h_number,self.d_key).permute(0, 2, 3, 1)
		v = self.W_v(x).view(bs,-1,self.h_number,self.d_value).permute(0, 2, 1, 3)
		#######################
		# att = quare * keys
		att = torch.matmul(q,k)
		print(att.size())
		att = att / np.sqrt(self.d_key)
		#######################
		if attention_weights is not None:
			att = att * attention_weights
		if attention_mask is not None:
			att = att.masked_fill(attention_mask,-np.inf)
		if self.ApplyTransform == True:
			att = self.transform_block(att)
		else:
			att = torch.softmax(att,-1)
		att = self.dropout(att)
		out = torch.matmul(att,v).permute(0,2,1,3)
		out = out.contiguous().view(bs,nq,self.h_number*self.d_value)
		y = self.W_o(out)
		return y

##############################################
##############################################
#### CNAM attention


##############################################
##############################################
#### CNAM attention


if __name__ == "__main__":
	inputs = torch.randn((50,49,512))
	#sa = ScaledDotProduct_attention_block(d_model=512,d_key=512,d_value=512,h_number=8)
	#output = sa(queries=inputs,keys=inputs,values=inputs)
	#print(output.size())

	emsa = Multihead_self_attention_block(d_model=512,d_key=512,d_value=512,h_number=2,H=7,W=7,header_num=2,ApplyTransform=True)
	output = emsa(queries=inputs,keys=inputs,values=inputs)
	print(output)
	print(output.size())