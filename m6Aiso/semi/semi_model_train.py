import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.backends.cudnn as cudnn
import torch.utils
import numpy as np
from torch.autograd import Variable

def Predict_2d_model_train(model,epochs,k_train_loader,optimizer,criterion,grad_clip=5):
	for epoch in range(epochs):
		for step,(input_train,motif_train,target_train) in enumerate(k_train_loader):
			input_train = input_train.transpose(3,2)
			motif_train = motif_train.transpose(2,1)
			########################
			input_train = Variable(input_train,requires_grad=False).cuda()
			motif_train = Variable(motif_train,requires_grad=False).cuda()
			input_target = Variable(target_train,requires_grad=False).cuda()
			########################
			optimizer.zero_grad()
			y_pred_read,y_pred_site = model(value=input_train,motif=motif_train)
			y_true = input_target
			loss = criterion(y_pred=y_pred_read,y_true=y_true)
			loss.backward()
			nn.utils.clip_grad_norm_(model.parameters(),grad_clip)
			optimizer.step()
			########################
			########################
	return model

def Predict_1d_model_train(model,epochs,k_train_loader,optimizer,criterion,grad_clip=5):
	for epoch in range(epochs):
		for step,(input_train,motif_train,target_train) in enumerate(k_train_loader):
			input_train = input_train.transpose(2,1)
			motif_train = motif_train.transpose(2,1)
			########################
			input_train = Variable(input_train,requires_grad=False).cuda()
			motif_train = Variable(motif_train,requires_grad=False).cuda()
			input_target = Variable(target_train,requires_grad=False).cuda()
			########################
			optimizer.zero_grad()
			y_pred_read,y_pred_site = model(value=input_train,motif=motif_train)
			y_true = input_target
			loss = criterion(y_pred=y_pred_read,y_true=y_true)
			loss.backward()
			nn.utils.clip_grad_norm_(model.parameters(),grad_clip)
			optimizer.step()
			########################
			########################
	return model

def Predict_model_train(model_type,model,epochs,k_train_loader,optimizer,criterion):
	if model_type == "1d":
		out_model = Predict_1d_model_train(model,epochs,k_train_loader,optimizer,criterion,grad_clip=5)
	if model_type == "2d":
		out_model = Predict_2d_model_train(model,epochs,k_train_loader,optimizer,criterion,grad_clip=5)
	return model

