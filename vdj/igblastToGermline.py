# Replaces VDJ sequence with germline VDJ
# Input: changeo MakeDb IgBLAST output


'''
# decodeCIGAR breaks MakeDb cigar string into tuple: [(),()]
import importlib.util
spec = importlib.util.spec_from_file_location("changeo", "/users/PAS0472/osu8725/anaconda3/lib/python3.6/site-packages/changeo/Alignment.py")
changeo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(changeo)
decodeCIGAR=changeo.decodeCIGAR
'''
# 

def imgtToVDJ(string):
	string = string.split('.')
	return ''.join(string)

from sys import argv

script, igblast = argv

with open(igblast, 'r') as f:
	for line in f:
		try:
			line = line.strip().split('\t')
			NAME = line[0]
			SEQUENCE = line[1]
			VDJ_CALL = ' | '.join([line[7], line[8], line[9]])

			# Where V,D,J starts/ends relative to query sequence
			V_SEQ_START = int(line[12]) - 1
			J_SEQ_END = int(line[24]) + int(line[25]) - 1
			SEQUENCE_VDJ = line[10]

			# Where V,D,J starts/ends on query relative to germline
			GERMLINE_IMGT = line[30]
			GERMLINE_VDJ = imgtToVDJ(GERMLINE_IMGT)
			V_GERM_START_VDJ = int(line[14])-1
			J_GERM_END_VDJ = V_GERM_START_VDJ + int(line[17]) + int(line[18]) + int(line[22]) + int(line[23]) + int(line[27]) - 1


			# Replace just the portion of the query that is VDJ associated with the germline
			# Retain other sequences on 5' and 3' ends of query (if they exist)

			fiveprime_keep = SEQUENCE[0:V_SEQ_START]
			threeprime_keep = SEQUENCE[J_SEQ_END:]
			germline_keep = GERMLINE_VDJ[V_GERM_START_VDJ:J_GERM_END_VDJ+1]
			print(V_GERM_START_VDJ,J_GERM_END_VDJ+1)
			new_germline = fiveprime_keep + '|' + germline_keep + '|' + threeprime_keep
			print('>%s %s\n%s' % (NAME, VDJ_CALL, new_germline))


			# Print 
			FWR1_END = V_SEQ_START + len(imgtToVDJ(line[43]))
			CDR1_END = FWR1_END + len(imgtToVDJ(line[47]))
			FWR2_END = CDR1_END + len(imgtToVDJ(line[44]))
			CDR2_END = FWR2_END + len(imgtToVDJ(line[48]))
			FWR3_END = CDR2_END + len(imgtToVDJ(line[45]))
			CDR3_END = FWR3_END + len(imgtToVDJ(line[49]))
			FWR4_END = CDR3_END + len(imgtToVDJ(line[46]))

			regions = [0, V_SEQ_START, FWR1_END, CDR1_END, FWR2_END, CDR2_END, FWR3_END, CDR3_END, FWR4_END, len(SEQUENCE_VDJ)]
			region_length = len(regions)-1
			region_names = ['-', 'FWR1', 'CDR1', 'FWR2', 'CDR2', 'FWR3', 'CDR3', 'FWR4', '-']

			for n, r in enumerate(regions):
				if n != region_length:
					for i in range(regions[n], regions[n+1]):
						#print(region_names[n])
						hi = region_names[n]
				else:
					break


		except:
			pass