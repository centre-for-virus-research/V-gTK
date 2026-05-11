import os
import pandas as pd

'''
Load text file 
'''
class TextFileLoader:

	def __init__(self, filepath, delimiter):
		self.gb_matrix = gb_matrix
		self.delimiter= delimiter

	def load_gb_matrix(self):
		df = pd.read_csv(filepath, delimiter=delimiter)
		return df




