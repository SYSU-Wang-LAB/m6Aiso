import gzip
#########################################

def signal_readname_select(filename):
	f = gzip.open(filename,'rb')
	readname_to_count = {}
	for str_x in f:
		str_x = str_x.decode()
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		#############################
		if list_x[0] == "chrom":
			index_chrom = list_x.index("chrom")
			index_genome_position = list_x.index("genome_position")
			index_genome_motif = list_x.index("genome_motif")
			index_strand = list_x.index("strand")
			index_transcript_id = list_x.index("transcript_id")
			index_transcript_position = list_x.index("transcript_position")
			index_up_and_down_sequecne = list_x.index("up_and_down_sequecne")
			index_readname = list_x.index("readname")
			continue
		#############################
		chrom = list_x[index_chrom]
		genome_position = list_x[index_genome_position]
		genome_motif = list_x[index_genome_motif]
		strand = list_x[index_strand]
		transcript_id = list_x[index_transcript_id]
		transcript_position = list_x[index_transcript_position]
		up_and_down_sequecne = list_x[index_up_and_down_sequecne]
		readname = list_x[index_readname]
		#############################
		try:
			readname_to_count[readname] +=1
		except:
			readname_to_count[readname] =1
	f.close()
	return readname_to_count


def bam_tsv_fileread(filename):
	"""
	readname        transcript_id   baseread        baseref refposi quality
	6715ca76-9fb7-4169-b4ad-7e1371b86e9c    ENST00000440196.3|ENSG00000230021.10|OTTHUMG00000191652.4|OTTHUMT00000493605.1|ENST00000440196|ENSG00000230021|1022|processed_transcript|       GGACAGCTCATGAGTGCAAGACGTCTTGTGATGTAATTATTATACGAATAGGGGCTCTCAATCGGGAGTACTACTCGATTGTCAACGTACAAGGAGTCGCAGGTCGCCTGGTTCTAGGAATAATGGGGAAGTATGTAGGAGTTGAAGATTAGTCCGCCGTAGTCGGTGTACTCGTAGGTTCAGTTGTCATCGGCGGCCAATTGATTTGATAGTTGGAGAGGGATCGTTGACCTCGTCTGTTATGTAAAGGATGCGTAGGGATGGGAGGGCGATGAGGACTAGGATGATGGCGGGCAGGATAGTTCAGACAGGTTTCTATTCCTGAGCGTCTGAGATGTTAGTATTAGTTAGTTTGTTGTGAGTGTTAGGAAAAGGGCATACAGGACTGGGAAGCAGATAAGGAAATGATTATGAGGGCGTGATCATGAAAGGTGATAAGCTCTTTTATGATAGGGGAAGTAGCGTCTTGTAGACCTACTTGCGCTGCATGTGCCATTAAGATATAGGATTTAGCCTATAATTTAACTTTGACAAGTTATGAATGGTTTTCTAATACCTTTTAAAAAAACAG     ............................................................................................................................GGGGAAGTATGTAGGAGTTGAAGATTAGTCCGCCGTAGTCGGTGTATTCGTAGGTTCAG.TACCATTGATGGCCAATTGATTTGATGACCTTAGTTTAGGTATTGGGGCCAAAGGATGGATGACCATTTCAAACGATCCAGGCTAAGCCAGGAGGAGAGCTCAAAGTCTGATCTGCTCTGCTGCCCCCTGCCCCATACACGTGATGGAGCAGAAAACGTGCTGTGTGAACCTGTGACTTCAGGGCCTGTTGACGTGGTCGTGCTTGCATACTCTCTGGACTGGACCTCACTGTGGGAACAACAAGATCAACAAGAGGAGCAAGAACAACATCAAGAGTCAGGGCCCGGGGGTCCTGACGGGTACAGGATGGGTACAGACCCACACAGGAATCCCAGAGTGTGTTCCACAGCAGGACACGCCTGCGCTGAAAGAGTGGGCAGAAAGGA     -123;-122;-121;-120;-119;-118;-117;-116;-115;-114;-113;-112;-111;-110;-109;-108;-107;-106;-105;-104;-103;-102;-101;-100;-99;-98;-97;-96;-95;-94;-93;-92;-91;-90;-89;-88;-87;-86;-85;-84;-83;-82;-81;-80;-79;-78;-77;-76;-75;-74;-73;-72;-71;-70;-69;-68;-67;-66;-65;-64;-63;-62;-61;-60;-59;-58;-57;-56;-55;-54;-53;-52;-51;-50;-49;-48;-47;-46;-45;-44;-43;-42;-41;-40;-39;-38;-37;-36;-35;-34;-33;-32;-31;-30;-29;-28;-27;-26;-25;-24;-23;-22;-21;-20;-19;-18;-17;-16;-15;-14;-13;-12;-11;-10;-9;-8;-7;-6;-5;-4;-3;-2;.;0;1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;.;60;61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;84;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;115;116;117;118;119;120;121;122;123;124;125;126;127;128;129;130;131;132;133;134;135;136;137;138;139;140;141;142;143;144;145;146;147;148;149;150;151;152;153;154;155;156;157;158;159;160;161;162;163;164;165;166;167;168;169;170;171;172;173;174;175;176;177;178;179;180;181;182;183;184;185;186;187;188;189;190;191;192;193;194;195;196;197;198;199;200;201;202;203;204;205;206;207;208;209;210;211;212;213;214;215;216;217;218;219;220;221;222;223;224;225;226;227;228;229;230;231;232;233;234;235;236;237;238;239;240;241;242;243;244;245;246;247;248;249;250;251;252;253;254;255;256;257;258;259;260;261;262;263;264;265;266;267;268;269;270;271;272;273;274;275;276;277;278;279;280;281;282;283;284;285;286;287;288;289;290;291;292;293;294;295;296;297;298;299;300;301;302;303;304;305;306;307;308;309;310;311;312;313;314;315;316;317;318;319;320;321;322;323;324;325;326;327;328;329;330;331;332;333;334;335;336;337;338;339;340;341;342;343;344;345;346;347;348;349;350;351;352;353;354;355;356;357;358;359;360;361;362;363;364;365;366;367;368;369;370;371;372;373;374;375;376;377;378;379;380;381;382;383;384;385;386;387;388;389;390;391;392;393;394;395;396;397;398;399;400;401;402;403;404;405;406;407;408;409;410;411;412;413;414;415;416;417;418;419;420;421;422;423;424;425;426;427;428;429;430;431;432;433;434;435;436;437;438;439;440;441;442;443;444;445;446       .78C3/*028<D>629=DI>&$8747A93@1-7;:740)5038<81.;6+479>6*$399:133A;A;B67##%91726;6,+1304$(,;<D?4<;96>;@A569<-+:)'&)05:64404)-8:CA--?<>%5./BG;45/,7@7DE@?=:1?2:6*'$9)%&#,:7/.-*41344-<85<6&)$'%''*))31868B37<*368?B8)9+'*-&/6/>=AA7>1=<A=665957;82<==@<BHLEDPI85:@0<<=:%&233.<;71=GA@.1(&-/ACFGA>9$"$+5432:>D?84A8D@7AG,7641&)52109/853744/'%%94BB:0-718?587<:;<63/31099*#)$2;;96>?@3JI@>4.*-.2558C8@02/=;5;<:88D=4=<7)1;=%(<%5*)3257??@:<698AJRDA<CDEBDA/+4$2-2*(2:5,4CHG986/(),,,,(;71,7>BA5>@;B9<C::9:69<@86%$$23===:6,+8;AB3?<>:3:9:<9:''&<%0+9:;:8<=-=65=798DC56*/+331/1-./-$//...-,,+)+
	"""
	f = gzip.open(filename,'rb')
	out_readname_to_bam_dict = {}
	for str_x in f:
		str_x = str_x.decode()
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "readname":
			continue
		readname = list_x[0]
		transcript_id = list_x[1]
		baseread = list_x[2]
		baseref = list_x[3]
		refposi = list_x[4]
		quality = list_x[5]
		##################################
		out_readname_to_bam_dict[(readname,transcript_id)] = [baseread,refposi]
		##################################
	f.close()
	return out_readname_to_bam_dict


def up_and_down_sequecne_make(baseread,refposi,transcript_position):
	transcript_position = int(transcript_position)
	refposi = [x for x in refposi.split(";")]
	outlist = []
	for i in range(len(refposi)):
		refposi_i = refposi[i]
		if refposi_i == ".":
			continue
		refposi_i = int(refposi[i])
		baseread_i = baseread[i]
		############################
		if refposi_i == transcript_position-1:
			outlist.append([baseread_i,refposi_i])
		if refposi_i == transcript_position+0:
			outlist.append([baseread_i,refposi_i])
		if refposi_i == transcript_position+1:
			outlist.append([baseread_i,refposi_i])
		if refposi_i == transcript_position+2:
			outlist.append([baseread_i,refposi_i])
		if refposi_i == transcript_position+3:
			outlist.append([baseread_i,refposi_i])
		if refposi_i == transcript_position+4:
			outlist.append([baseread_i,refposi_i])
		if refposi_i == transcript_position+5:
			outlist.append([baseread_i,refposi_i])
		if refposi_i == transcript_position+6:
			outlist.append([baseread_i,refposi_i])
		if refposi_i == transcript_position+7:
			outlist.append([baseread_i,refposi_i])
		if refposi_i == transcript_position+8:
			outlist.append([baseread_i,refposi_i])
		if refposi_i == transcript_position+9:
			outlist.append([baseread_i,refposi_i])
		if refposi_i == transcript_position+10:
			outlist.append([baseread_i,refposi_i])
		if refposi_i == transcript_position+11:
			outlist.append([baseread_i,refposi_i])
		############################
	outlist.sort(key=lambda x:int(x[1]))
	outstr = "".join([x[0] for x in outlist])
	############################
	return outstr


def signal_fileread(filename,out_readname_to_bam_dict,writefile):
	"""
	chrom   genome_position genome_motif    strand  transcript_id   transcript_position     readname	up_and_down_sequecne    mean    stdv    length  samples
	chr14   22008033        TTCCTTATT       +       ENST00000390447.3       317     b8b8a1e0-decd-4d21-88eb-260780a7c9d3    TTCCTTATT       72.16000000000001|73.14|74.30133333333333|87.36|96.56  1.5198461538461538|1.613|1.0717333333333334|2.491|3.198  0.026000000000000002|0.00500|0.005|0.00900|0.00700      71.5026,74.4007,70.923,70.0535,75.8498,71.7924,70.6331,72.0822,72.9517,73.8211,68.8943,73.5313,73.5313,72.6618,73.5313,73.5313,72.9517,70.6331,73.2415,72.372|75.56,75.56,73.5313,73.5313,71.5026,71.5026,76.1396,76.1396,70.923,70.923,71.7924,70.923,73.5313,75.56,73.5313,76.4294,71.7924,72.372,73.2415,71.7924|72.9517,72.9517,77.009,77.009,74.6905,74.6905,75.8498,75.8498,74.1109,74.1109,76.4294,75.8498,80.197,76.1396,73.8211,72.9517,73.5313,71.2128,72.0822,72.372|94.6877,86.5729,91.2099,83.9646,88.6016,84.2544,90.0507,90.6303,90.6303,88.8914,88.022,90.6303,88.022,90.0507,89.1812,88.8914,86.8627,87.7322,88.022,88.8914|98.7451,96.1368,93.8183,100.194,92.3692,99.6145,99.0349,100.194,98.4553,92.0794,90.3405,94.3979,98.4553,97.0062,88.8914,94.9775,98.4553,97.8757,100.774,100.774
	"""
	f = gzip.open(filename,'rb')
	d = gzip.open(writefile,'ab')
	readname_to_count = {}
	for str_x in f:
		str_x = str_x.decode()
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		#############################
		if list_x[0] == "chrom":
			index_chrom = list_x.index("chrom")
			index_genome_position = list_x.index("genome_position")
			index_genome_motif = list_x.index("genome_motif")
			index_strand = list_x.index("strand")
			index_transcript_id = list_x.index("transcript_id")
			index_transcript_position = list_x.index("transcript_position")
			index_up_and_down_sequecne = list_x.index("up_and_down_sequecne")
			index_readname = list_x.index("readname")
			str2write = "\t".join(list_x+["base_calling_sequecne"]) + "\n"
			d.write(str2write.encode())
			continue
		#############################
		chrom = list_x[index_chrom]
		genome_position = list_x[index_genome_position]
		genome_motif = list_x[index_genome_motif]
		strand = list_x[index_strand]
		transcript_id = list_x[index_transcript_id]
		transcript_position = list_x[index_transcript_position]
		up_and_down_sequecne = list_x[index_up_and_down_sequecne]
		readname = list_x[index_readname]
		#############################
		#print(out_readname_to_bam_dict)
		#############################
		try:
			baseread,refposi = out_readname_to_bam_dict[(readname,transcript_id)]
		except:
			continue
		base_sequence = up_and_down_sequecne_make(baseread=baseread,refposi=refposi,transcript_position=transcript_position)
		#print(base_sequence)

		if len(base_sequence) == 13:
			str2write = "\t".join(list_x+[base_sequence]) + "\n"
			d.write(str2write.encode())
		#############################
		#############################

def args_make():
	import argparse
	parser = argparse.ArgumentParser(description='m6A peak quality control by motif and miCLIP dataset overlap')
	parser.add_argument('--signal_filename', required=True, help="m6A peak file")
	parser.add_argument('--bam_convert_filename', required=True, help="m6A site file")
	parser.add_argument('--output_filename', required=True, help="m6A site file")
	args = parser.parse_args()
	return args

if __name__ == "__main__":
	import sys
	args = args_make()
	signal_filename = args.signal_filename
	bam_convert_filename = args.bam_convert_filename
	output_filename = args.output_filename
	###################################
	###################################
	#readname_to_count_dict = signal_readname_select(filename=signal_filename)
	###################################
	#print(readname_to_count_dict)
	###################################
	readname_to_bam_dict = bam_tsv_fileread(
		filename=bam_convert_filename
		)
	#print(readname_to_bam_dict)
	signal_fileread(
		filename=signal_filename,
		out_readname_to_bam_dict=readname_to_bam_dict,
		writefile=output_filename
		)

		

	