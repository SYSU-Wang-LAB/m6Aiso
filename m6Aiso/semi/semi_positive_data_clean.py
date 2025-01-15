


def pred_result_to_dict(filename):
	f = open(filename,'r')
	out_sitename_to_value_dict = {}
	out_sitename_list = []
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		sitename = list_x[0]
		pred_value = float(list_x[1])
		out_sitename_to_value_dict[sitename] = pred_value
		out_sitename_list.append(sitename)
	f.close()
	return out_sitename_to_value_dict,out_sitename_list

def positive_result_write(old_pos_filename,new_pos_filename,valuedict):
	f = open(old_pos_filename,'r')
	d = open(new_pos_filename,'a')
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "kmer":
			d.write(str_x+"\n")
			continue
		sitename = list_x[1]
		try:
			value = valuedict[sitename]
			d.write(str_x+"\n")
		except:
			pass
	f.close()
	d.close()


def site_select_by_percentage(valuedict,sitename_list,cutoff_perventage,reverse=True):
	unsort_value_list = []
	for sitename in sitename_list:
		value = valuedict[sitename]
		unsort_value_list.append([sitename,value])
	# reverse = False --> [[1, 2], [2, 3], [3, 7], [4, 5]]
	# reverse = True  --> [[4, 5], [3, 7], [2, 3], [1, 2]]
	unsort_value_list.sort(key=lambda x:x[1],reverse=reverse)
	cutoff_number = int(cutoff_perventage * len(unsort_value_list))
	out_dict = {}
	for i in range(cutoff_number):
		sitename = unsort_value_list[i][0]
		value = unsort_value_list[i][1]
		out_dict[sitename] = value
	return out_dict

def false_positive_rate_count(sitename2value_dict,sitename_list,cutoff):
	false_positive_count = 0
	true_positive_count = 0
	for sitename_i in sitename_list:
		value_i = sitename2value_dict[sitename_i]
		if value_i >= cutoff:
			true_positive_count += 1
		else:
			false_positive_count += 1
	FPR = false_positive_count / (false_positive_count + true_positive_count)
	return FPR

def positive_data_clean_and_FPR_count(pred_result_filename,uncleaned_positive_filename,cleaned_positive_filename,cutoff=0.5,cutoff_perventage=0.95):
	sitename2value_dict,sitename_list = pred_result_to_dict(filename=pred_result_filename)
	FPR = false_positive_rate_count(
		sitename2value_dict=sitename2value_dict,
		sitename_list=sitename_list,
		cutoff = cutoff
		)
	######################################
	selected_sitename2value_dict = site_select_by_percentage(
		valuedict=sitename2value_dict,
		sitename_list=sitename_list,
		cutoff_perventage=cutoff_perventage,
		reverse=True
		)
	######################################
	positive_result_write(
		old_pos_filename = uncleaned_positive_filename,
		new_pos_filename = cleaned_positive_filename,
		valuedict = selected_sitename2value_dict
		)
	return FPR
	################################
