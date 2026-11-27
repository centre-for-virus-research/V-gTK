import re
import os
from urllib.parse import unquote

'''
This class creates the dictionary for a given gff3 file
'''
class GffDictionary:

	def __init__(self, gff_file):
		self.gff_file = gff_file
		self.gff_dict = self.create_gff_dict()

	@staticmethod
	def _parse_attrs(attr_str):
		attrs = {}
		for part in attr_str.split(';'):
			if '=' in part:
				k, v = part.split('=', 1)
				attrs[k.strip()] = unquote(v.strip())
		return attrs

	def create_gff_dict(self):
		gff_dict = {}
		raw_cds = []
		mature_regions = []
		with open(self.gff_file, encoding='utf-8', errors='replace') as f:
			for each_line in f:
				stripped = each_line.strip()
				if stripped == '##FASTA':
					break
				if not stripped or stripped.startswith('#'):
					continue
				fields = stripped.split('\t')
				if len(fields) < 9:
					continue
				seqid, _source, feature, start, end, _score, _strand, _phase, attributes = fields[:9]
				attrs = self._parse_attrs(attributes)
				product = attrs.get('product')
				if product is None:
					match = re.search(r'product=([^;]+)', attributes)
					product = unquote(match.group(1)) if match else None
				entry = {'start': start, 'end': end, 'product': product}
				meta = {
					'seqid': seqid,
					'start_int': int(start),
					'end_int': int(end),
					'id': attrs.get('ID', ''),
					'parent': attrs.get('Parent', ''),
					'entry': entry,
				}

				gff_dict.setdefault(feature, []).append(entry)
				if feature == 'CDS':
					raw_cds.append(meta)
				elif feature == 'mature_protein_region_of_CDS':
					mature_regions.append(meta)

		if mature_regions:
			mature_parent_cds = [
				cds
				for cds in raw_cds
				if any(
					(cds['id'] and m['parent'] and cds['id'] in m['parent'].split(','))
					or (
						m['seqid'] == cds['seqid']
						and m['start_int'] >= cds['start_int']
						and m['end_int'] <= cds['end_int']
					)
					for m in mature_regions
				)
			]
			standalone_cds = [
				cds
				for cds in raw_cds
				if not any(
					cds['seqid'] == p['seqid']
					and cds['start_int'] >= p['start_int']
					and cds['end_int'] <= p['end_int']
					for p in mature_parent_cds
				)
			]
			combined = sorted(
				standalone_cds + mature_regions,
				key=lambda item: (item['start_int'], item['end_int']),
			)
			gff_dict['CDS'] = [item['entry'] for item in combined]
		elif raw_cds:
			gff_dict['CDS'] = [cds['entry'] for cds in raw_cds]

		return gff_dict


