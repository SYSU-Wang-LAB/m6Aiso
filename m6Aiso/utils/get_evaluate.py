import os
import numpy as np


class model_evaluate():
	def __init__(self,y_pred,y_true):
		self.y_pred = y_pred
		self.y_true = y_true

	def statistic_evaluate(self,cutoff):
		TP,TN,FP,FN = self.TP_TN_FP_FN_count(cutoff=cutoff)
		acc = self.accurracy(TP,TN,FP,FN)
		sen = self.sensitivity(TP,TN,FP,FN)
		spe = self.specificity(TP,TN,FP,FN)
		pre = self.precision(TP,TN,FP,FN)
		rec = self.recall(TP,TN,FP,FN)
		F1 = self.F1_score(TP,TN,FP,FN)
		str_result_1 = "cutoff is:%f,result is:%f\t%f\t%f\t%f"%(cutoff,TP,TN,FP,FN)
		str_result_2 = "cutoff is:%f,acc is:%f\t%f\t%f\t%f\t%f\t%f"%(cutoff,acc,sen,spe,pre,rec,F1)
		return str_result_1,str_result_2

	def TP_TN_FP_FN_count(self,cutoff):
		TP_count = 0
		TN_count = 0
		FP_count = 0
		FN_count = 0
		for i in range(len(self.y_pred)):
			if self.y_pred[i] >= cutoff:
				if self.y_true[i] == 1:
					TP_count+=1
				else:
					FP_count+=1
			else:
				if self.y_true[i] == 0:
					TN_count+=1
				else:
					FN_count+=1
		return TP_count,TN_count,FP_count,FN_count

	def accurracy(self,TP,TN,FP,FN):
		return (TP+TN)/(TP+TN+FP+FN)

	def sensitivity(self,TP,TN,FP,FN):
		return (TP)/(TP+FN)

	def specificity(self,TP,TN,FP,FN):
		return (TN)/(TN+FP)

	def precision(self,TP,TN,FP,FN):
		if TP+FP == 0:
			return 1.0
		return (TP)/(TP+FP)

	def recall(self,TP,TN,FP,FN):
		return (TP)/(TP+FN)

	def F1_score(self,TP,TN,FP,FN):
		return (2*TP)/(2*TP+FP+FN)

	def MCC(self,TP,TN,FP,FN):
		z = (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)
		e = np.sqrt(z)
		return (TP*TN - FP*FN)/e
