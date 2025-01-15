######################
import math
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from torch.autograd import Variable
######################

class Positionnal_encoding(nn.Module):
	def __init__(self,d_model,max_len,dropout_ratio=0.1):
		"""
		d_model: dimension of model
		max_len: max sequence length
		arange(0,10) --> tensor(0,1,2,3,4,5,6,7,8,9)
		torch.arange(0,10,step=2) --> tensor(0,2,4,6,8)
		unsqueeze(dim=1) --> convert tensor from size=(100) to site=(100,1) 
		"""
		super(Positionnal_encoding,self).__init__()
		self.posi_encoding = torch.zeros(max_len,d_model)
		self.posi_encoding.requires_grad = False
		position = torch.arange(0,max_len).float().unsqueeze(1)
		divterm = torch.exp(torch.arange(0,d_model,step=2) * (-(np.log(10.0)/d_model)))
		self.posi_encoding[:,0::2] = torch.sin(position*divterm)
		self.posi_encoding[:,1::2] = torch.cos(position*divterm)
		self.dropout = nn.Dropout(p=dropout_ratio)
	def forward(self,x):
		batch_size = x.size(0)
		sequence_length = x.size(-2)
		print(self.posi_encoding.size())
		return x+self.posi_encoding[:sequence_length,:]

class TokenEmbedding(nn.Embedding):
	"""
	vocab_size: size of vocabulary
	d_model:dimensions of model
	"""
	def __init__(self,vocab_size,d_model):
		super(TokenEmbedding,self).__init__(vocab_size,d_model)

class TransformerEmbedding(nn.Module):
	def __init__(self,vocab_size,d_model,max_len,drop_ratio):
		super(TransformerEmbedding,self).__init__()
		self.token_embedding = TokenEmbedding(vocab_size,d_model)
		self.posit_embedding = Positionnal_encoding(d_model,max_len)
	def forward(self,x):
		t = self.token_embedding(x)
		p = self.posit_embedding(x)
		return self.dropout(t + p)

class DotProduct_attention(nn.Module):
	def __init__(self):
		super(DotProduct_attention,self).__init__()
		self.softmax = nn.Softmax(dim=-1)
	def forward(self,q,k,v,mask=None,dropout=None,e=-10000):
		batch_size,header_num,length,d_tensor = k.size()
		d_model = header_num * d_tensor
		k_t = k.transpose(-2,-1)
		alpha_score = torch.matmul(q,k_t) / math.sqrt(d_tensor)
		if mask is not None:
			alpha_score = alpha_score.masked_fill(mask==0,e)
		alpha_score = self.softmax(alpha_score)
		if dropout is not None:
			alpha_score = dropout(alpha_score)
		v = torch.matmul(alpha_score,v)
		return v,alpha_score

class MultiHead_attention(nn.Module):
	def __init__(self,d_model,header_num,dropout_ratio=None):
		super(MultiHead_attention,self).__init__()
		self.header_num = header_num
		self.attention_block = DotProduct_attention()
		self.W_q = nn.Linear(d_model,d_model)
		self.W_k = nn.Linear(d_model,d_model)
		self.W_v = nn.Linear(d_model,d_model)
		self.W_concat = nn.Linear(d_model,d_model)
		if not dropout_ratio == None:
			self.dropout = nn.Dropout(p=dropout_ratio)
		else:
			self.dropout = None
	def forward(self,q,k,v,mask=None):
		q = self.W_q(q)
		k = self.W_q(k)
		v = self.W_q(v)
		q = self.split(q)
		k = self.split(k)
		v = self.split(v)
		out,att = self.attention_block(q=q,k=k,v=v,mask=mask,dropout=self.dropout,e=-10000)
		out = self.concat(out)
		out = self.W_concat(out)
		return out

	def split(self,input_tensor):
		"""
		input tensor:[batch_size,length,d_model]
		output tensor:[batch_size,head_number,length,d_tensor]
		"""
		batch_size,length,d_model = input_tensor.size()
		d_tensor = d_model//self.header_num
		output_tensor = input_tensor.view(batch_size,length,self.header_num,d_tensor).transpose(1, 2)
		return output_tensor

	def concat(self,input_tensor):
		"""
		input tensor:[batch_size,head_number,length,d_tensor]
		output tensor:
		"""
		batch_size,header_number,length,d_tensor = input_tensor.size()
		d_model = header_number * d_tensor
		output_tensor = input_tensor.transpose(1,2).contiguous().view(batch_size,length,d_model)
		return output_tensor

class LayerNorm(nn.Module):
	def __init__(self,d_model,eps=1e-12):
		super(LayerNorm,self).__init__()
		self.gamma = nn.Parameter(torch.ones(d_model))
		self.beta = nn.Parameter(torch.zeros(d_model))
		self.eps = eps
	def forward(self,x):
		# -1 mean last dimension
		mean = x.mean(-1,keepdim=True)
		var = x.var(-1,unbiased=False,keepdim=True)
		out = (x-mean)/torch.sqrt(var+self.eps)
		out = self.gamma * out + self.beta
		return out

class PositionwiseFeedForward(nn.Module):
	def __init__(self,d_model,hidden_size,drop_ratio=0.1):
		super(PositionwiseFeedForward,self).__init__()
		self.linear1 = nn.Linear(d_model,hidden_size)
		self.linear2 = nn.Linear(hidden_size,d_model)
		self.relu = nn.ReLU()
		self.dropout = nn.Dropout(p=drop_ratio)
	def forward(self,x):
		out = self.linear1(x)
		out = self.relu(out)
		out = self.dropout(out)
		out = self.linear2(out)
		return out

class EncoderLayer(nn.Module):
	def __init__(self):
		super(EncoderLayer,self).__init__()
		self.attention_block = MultiHead_attention(d_model=d_model,hidden_size=hidden_size,dropout_ratio=dropout_ratio)
		self.normal_layer_1 = LayerNorm(d_model=d_model)
		self.dropout_layer_1 = nn.Dropout(p=dropout_ratio)
		self.ffn = PositionwiseFeedForward(d_model=d_model,hidden_size=hidden_size,dropout_ratio=dropout_ratio)
		self.normal_layer_2 = LayerNorm(d_model=d_model)
		self.dropout_layer_2 = nn.Dropout(p=dropout_ratio)

	def forward(self,x):
		# 1. compute self attention
		x_old_1 = x
		att = self.attention_block(q=x,k=x,v=x,mask=mask)
		# 2. add and normalize
		att = self.dropout_layer_1(att)
		att = self.normal_layer_1(att + x_old_1)
		# 3. positionwise feed forward network
		x_old_2 = att
		ffn_value = self.ffn(att)
		# 4. add and normal
		y = self.dropout_layer_2(ffn_value)
		y = self.normal_layer_2(y+x_old_2)
		return x

class Encoder(nn.Module):
	def __init__(self):
		super(Encoder,self).__init__()
		self.emd = TransformerEmbedding(
			d_model=d_model,
			max_len=max_len,
			vocab_size=vocab_size,
			dropout_ratio=dropout_ratio)
		layers = []
		for i in range(num_layers):
			layers.append(EncoderLayer(
				d_model=d_model,
				hidden_size=hidden_size,
				header_num=header_num,
				dropout_ratio=dropout_ratio
				)
			)

		self.layers = nn.ModuleList(layers)
	def forward(self,x):
		x = self.emd(x)
		for layer in self.layers:
			x = layer(x)
		return x

class DecoderLayer(nn.Module):
	def __init__(self,d_model,hidden_size,header_num,dropout_ratio):
		super(DecoderLayer,self).__init__()
		self.self_attention = MultiHead_attention(d_model=d_model,header_num=header_num)
		self.normal_layer_1 = LayerNorm(d_model=d_model)
		self.dropout_layer_1 = nn.Dropout(p=dropout_ratio)
		self.enc_dec_attention = MultiHead_attention(d_model=d_model,header_num=header_num)
		self.normal_layer_2 = LayerNorm(d_model=d_model)
		self.dropout_layer_2 = nn.Dropout(p=dropout_ratio)
		self.ffn = PositionwiseFeedForward(d_model=d_model,hidden_size=hidden_size,dropout_ratio=dropout_ratio)
		self.normal_layer_3 = LayerNorm(d_model=d_model)
		self.dropout_layer_3 = nn.Dropout(p=dropout_ratio)
	def forward(self,x):
		pass

		


	
def main():
	pe = Positionnal_encoding(d_model=8,max_len=5)
	x = Variable(torch.zeros(100,32,5,8))
	ma = MultiHead_attention(d_model=32,header_num=2,dropout_ratio=None)
	y1 = pe(x)
	print(y1.size())
	y1 = y1.reshape((100,32,40))
	print(y1[0,0,:])
	y1 = y1.transpose(1, 2)
	print(y1[0,:,0])
	print(y1.size())
	ma(q=y1,k=y1,v=y1)
	a = torch.arange(0,10,step=2)
	print(a.size())

if __name__=="__main__":
	main()