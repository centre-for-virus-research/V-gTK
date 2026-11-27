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

	@staticmethod
	def _extract_attr(attributes, key):
		match = re.search(rf'(?:^|;){re.escape(key)}=([^;]+)', attributes)
		return unquote(match.group(1)) if match else None

	@staticmethod
	def _spans_overlap(a_start, a_end, b_start, b_end):
		try:
			s1, e1 = int(a_start), int(a_end)
			s2, e2 = int(b_start), int(b_end)
		except (TypeError, ValueError):
			return False
		if s1 > e1 or s2 > e2:
			return False
		return s1 <= e2 and s2 <= e1

	def create_gff_dict(self):
		gff_dict = {}
		raw_cds = []
		mature_regions = []
		mature_parents = set()
		seen_mature = set()
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
					product = self._extract_attr(attributes, 'product')
				feat_id = attrs.get('ID', '')
				parent_id = attrs.get('Parent', '')
				entry = {'start': start, 'end': end, 'product': product}
				start_int = int(start) if str(start).isdigit() else 0
				end_int = int(end) if str(end).isdigit() else 0
				meta = {
					'seqid': seqid,
					'start_int': start_int,
					'end_int': end_int,
					'id': feat_id,
					'parent': parent_id,
					'entry': entry,
				}

				gff_dict.setdefault(feature, []).append(entry)
				if feature == 'CDS':
					raw_cds.append(meta)
				elif feature == 'mature_protein_region_of_CDS':
					if parent_id:
						for p in parent_id.split(','):
							if p.strip():
								mature_parents.add(p.strip())
					dedup_key = (seqid, start, end, product)
					if dedup_key not in seen_mature:
						seen_mature.add(dedup_key)
						mature_regions.append(meta)

		if mature_regions:
			# Replace only polyprotein CDSs (and nested frameshift CDSs inside them)
			# covered by mature peptide annotations, while retaining standalone
			# downstream CDSs (e.g. SARS-CoV-2 S, ORF3a, E, M, ORF6-ORF10).
			mature_parent_cds = [
				cds
				for cds in raw_cds
				if (cds['id'] and cds['id'] in mature_parents)
				or any(
					m['seqid'] == cds['seqid']
					and self._spans_overlap(cds['start_int'], cds['end_int'], m['start_int'], m['end_int'])
					for m in mature_regions
				)
			]
			standalone_cds = [
				cds
				for cds in raw_cds
				if not any(
					cds['seqid'] == p['seqid']
					and self._spans_overlap(cds['start_int'], cds['end_int'], p['start_int'], p['end_int'])
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


