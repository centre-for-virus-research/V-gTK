import re
import os

'''
This class creates the dictionary for a given gff3 file
'''
class GffDictionary:

	def __init__(self, gff_file):
		self.gff_file = gff_file
		self.gff_dict = self.create_gff_dict()

	def create_gff_dict(self):
		gff_dict = {}
		with open(self.gff_file) as f:
			for each_line in f:
				if not each_line.startswith('#'):
					seqid,source,feature,start,end,score,strand,phase,attributes = each_line.strip().split('\t')
					match = re.search(r'product=([^;]+)', attributes)
					product = match.group(1) if match else None

					if feature not in gff_dict:
						gff_dict[feature] = []
						gff_dict[feature].append({'start': start, 'end': end, 'product': product})
					else:
						gff_dict[feature].append({'start': start, 'end': end, 'product':product})
		
		return gff_dict


