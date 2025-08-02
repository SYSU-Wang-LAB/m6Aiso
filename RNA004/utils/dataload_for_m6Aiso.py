


class Datamake_for_read_level_mod_using(Dataset):
	def __init__(self,data_feature_1,data_feature_2,data_feature_3,data_feature_4,data_target):
		self.len = len(data_feature_1)
		self.features_1 = torch.tensor(data_feature_1,dtype=torch.float32)
		self.features_2 = torch.tensor(data_feature_2,dtype=torch.float32)
		self.features_3 = torch.tensor(data_feature_3,dtype=torch.float32)
		self.features_4 = torch.tensor(data_feature_4,dtype=torch.int)
		self.target = torch.tensor(data_target,dtype=torch.float32)
	def __getitem__(self,index):
		return self.features_1[index],self.features_2[index],self.features_3[index],self.features_4[index],self.target[index]
	def __len__(self):
		return self.len