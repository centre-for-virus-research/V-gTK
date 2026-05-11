import os
from Bio import SeqIO

class RemoveRedundantSequence:
	def __init__(self, fasta_file, output_file):
		self.fasta_file = fasta_file
		self.output_file = output_file 

	def remove_redundant_fasta(self):
		seen_sequences = set()
		unique_records = []

		for record in SeqIO.parse(self.fasta_file, "fasta"):
			seq = str(record.seq)
			if seq not in seen_sequences:
				seen_sequences.add(seq)
				unique_records.append(record)

		temp_file = self.output_file + ".tmp"

		SeqIO.write(unique_records, temp_file, "fasta")
		os.replace(temp_file, self.fasta_file)

		print(f"Removed redundancy: {len(unique_records)} unique sequences saved to '{self.fasta_file}'")
