import os
import csv
import time
import subprocess
from os.path import join
from argparse import ArgumentParser

class SoftwareVersionChecker:
	def __init__(self, tmp_dir, output_dir, table_name):
		self.tmp_dir = tmp_dir
		self.output_dir = output_dir
		self.table_name = table_name
		self.software_commands = {
			"Nextalign": ["nextalign", "--version"],
			"BLAST": ["blastn", "-version"],
			"MAFFT": ["mafft", "--version"],
			"IQ-TREE": ["iqtree2", "--version"],
			"FastTree": ["FastTree", "-help"], 
			"Python": ["python", "--version"] }
		self.software_versions = {}
        
		self.check_software_versions()
		self.save_versions_to_tsv()
    
	def get_software_version(self, software, command):
		try:
			result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
			output = result.stdout if result.stdout else result.stderr
            
			lines = output.split("\n")
			for line in lines:
				if "version" in line.lower() or software.lower() in line.lower():
					return line.strip()
			return output.strip().split("\n")[0]  # Fallback: Return first line of output
		except FileNotFoundError:
			return "Not Installed"
    
	def check_software_versions(self):
		self.software_versions = {
		software: self.get_software_version(software, command)
		for software, command in self.software_commands.items()
        }
    
	def save_versions_to_tsv(self):
		os.makedirs(self.tmp_dir, exist_ok=True)
		os.makedirs(join(self.tmp_dir, self.output_dir), exist_ok=True)
		output_path = os.path.join(self.tmp_dir, self.output_dir, self.table_name)
        
		with open(output_path, "w", newline="") as file:
			writer = csv.writer(file, delimiter="\t")
			writer.writerow(["Software", "Version"])
			for software, version in self.software_versions.items():
				writer.writerow([software, version])
			writer.writerow(["Time of creation", time.time()])
			writer.writerow(["vgtk version", "v.1.0.0"])
        
if __name__ == "__main__":
    parser = ArgumentParser(description='Create a TSV with available software version information')
    parser.add_argument('-d', '--tmp_dir', help='Temp directory', default="tmp")
    parser.add_argument('-o', '--output_dir', help='Output directory', default="Software_info")
    parser.add_argument('-f', '--table_name', help='Name of the TSV to be saved', default="software_info.tsv")
    args = parser.parse_args()
    
    checker = SoftwareVersionChecker(args.tmp_dir, args.output_dir, args.table_name)

