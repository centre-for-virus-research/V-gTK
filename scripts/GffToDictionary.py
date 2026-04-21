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
		raw_cds = []
		mature_regions = []
		with open(self.gff_file) as f:
			for each_line in f:
				if not each_line.startswith('#'):
					seqid,source,feature,start,end,score,strand,phase,attributes = each_line.strip().split('\t')
					match = re.search(r'product=([^;]+)', attributes)
					product = match.group(1) if match else None
					entry = {'start': start, 'end': end, 'product': product}

					if feature not in gff_dict:
						gff_dict[feature] = []
						gff_dict[feature].append(entry)
					else:
						gff_dict[feature].append(entry)

					if feature == 'CDS':
						raw_cds.append(entry)
					elif feature == 'mature_protein_region_of_CDS':
						mature_regions.append(entry)

		if mature_regions:
			gff_dict['CDS'] = mature_regions
		elif raw_cds:
			gff_dict['CDS'] = raw_cds
		
		return gff_dict


