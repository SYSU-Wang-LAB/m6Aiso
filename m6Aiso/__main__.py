#!/bin/python3
import argparse
import sys
sys.path.insert(0,"/home/ZJRen/pytorch/m6Aiso/m6Aiso")
sys.path.insert(0,"/home/ZJRen/pytorch/m6Aiso/")
sys.path.insert(0,"./")
sys.path.insert(0,"./m6Aiso/")
sys.path.insert(0,"./m6Aiso/blocks/")
print("\n".join(sys.path))
from m6Aiso.script import semi_supervised_model_train
from m6Aiso.script import molecular_m6A_predication
from m6Aiso.script import current_signal_abstract_for_m6A_pred

def get_version():
	version = "m6Aiso version:1.0.0"
	return version

def args_make():
	import argparse
	parser = argparse.ArgumentParser(description='model training')
	parser.add_argument('--version',action='version',version=get_version(),help="get version")
	subparsers = parser.add_subparsers(help="semi_model_train,m6A_predication")
	#################################
	parser_sub_1 = subparsers.add_parser('semi_supervised_model_train',help="semi supervised model training.")
	parser_sub_1 = semi_supervised_model_train.args_make(parser_sub_1)
	parser_sub_1.set_defaults(func=semi_supervised_model_train.main)
	#################################
	parser_sub_2 = subparsers.add_parser('molecular_m6A_predication',help="molecular level m6A predication.")
	parser_sub_2 = molecular_m6A_predication.args_make(parser_sub_2)
	parser_sub_2.set_defaults(func=molecular_m6A_predication.main)
	#################################
	parser_sub_3 = subparsers.add_parser('current_signal_abstract_for_m6A_pred',help="Split nanopolish and normalized the current intensity.")
	parser_sub_3 = current_signal_abstract_for_m6A_pred.args_make(parser_sub_3)
	parser_sub_3.set_defaults(func=current_signal_abstract_for_m6A_pred.main)

	args = parser.parse_args()
	return args

def main():
    args = args_make()
    args.func(args)
    
if __name__=="__main__":
	main()
