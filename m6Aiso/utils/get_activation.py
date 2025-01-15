import torch
from torch import nn

def get_activate_function(activation_name):

	"""
	this function is using to get as loss function by give activation function name
	"""

	allowed_activation = ('tanh', 'sigmoid', 'relu', 'softmax')
	out_activation_function = None
	#####################
	if activation_name == "relu":
		out_activation_function = nn.ReLU()
	elif activation_name == "tanh":
		out_activation_function = nn.Tanh()
	elif activation_name == "sigmoid":
		out_activation_function = nn.Sigmoid()
	elif activation_name == "softmax":
		out_activation_function = nn.Softmax(dim=1)
	else:
		raise ValueError("Invalid activation")
	return out_activation_function