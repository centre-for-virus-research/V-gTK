import re
import os

'''
This class creates the dictionary for a given gff3 file
'''
class GffDictionary:

	def __init__(self, gff_file):
		self.gff_file = gff_file
		self.gff_dict = self.create_gff_dict()

	@staticmethod
	def _extract_attr(attributes, key):
		match = re.search(rf'(?:^|;){re.escape(key)}=([^;]+)', attributes)
		return match.group(1) if match else None

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
		raw_cds_meta = []
		mature_regions = []
		mature_parents = set()
		seen_mature = set()
		with open(self.gff_file) as f:
			for each_line in f:
				if not each_line.startswith('#') and each_line.strip():
					seqid,source,feature,start,end,score,strand,phase,attributes = each_line.strip().split('\t')
					product = self._extract_attr(attributes, 'product')
					feat_id = self._extract_attr(attributes, 'ID')
					parent_id = self._extract_attr(attributes, 'Parent')
					entry = {'start': start, 'end': end, 'product': product}

					if feature not in gff_dict:
						gff_dict[feature] = []
						gff_dict[feature].append(entry)
					else:
						gff_dict[feature].append(entry)

					if feature == 'CDS':
						raw_cds.append(entry)
						raw_cds_meta.append((entry, feat_id))
					elif feature == 'mature_protein_region_of_CDS':
						if parent_id:
							for p in parent_id.split(','):
								if p.strip():
									mature_parents.add(p.strip())
						dedup_key = (start, end, product)
						if dedup_key not in seen_mature:
							seen_mature.add(dedup_key)
							mature_regions.append(entry)

		if mature_regions:
			# Replace only polyprotein CDSs (and nested frameshift CDSs inside them)
			# covered by mature peptide annotations, while retaining standalone
			# downstream CDSs (e.g. SARS-CoV-2 S, ORF3a, E, M, ORF6-ORF10).
			polyprotein_spans = []
			for entry, feat_id in raw_cds_meta:
				is_parent = bool(feat_id and feat_id in mature_parents)
				overlaps_mature = any(
					self._spans_overlap(entry['start'], entry['end'], m['start'], m['end'])
					for m in mature_regions
				)
				if is_parent or overlaps_mature:
					polyprotein_spans.append((entry['start'], entry['end']))

			standalone_cds = [
				entry for entry, _ in raw_cds_meta
				if not any(
					self._spans_overlap(entry['start'], entry['end'], p_start, p_end)
					for p_start, p_end in polyprotein_spans
				)
			]
			combined = mature_regions + standalone_cds
			combined.sort(key=lambda e: (int(e['start']), int(e['end'])) if str(e['start']).isdigit() and str(e['end']).isdigit() else (0, 0))
			gff_dict['CDS'] = combined
		elif raw_cds:
			gff_dict['CDS'] = raw_cds
		
		return gff_dict


