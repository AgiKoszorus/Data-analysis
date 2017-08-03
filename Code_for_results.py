# code for going through the analysis 
# using classes 
# should be universal

import os
import numpy as np
import matplotlib.pyplot as plt 
from sys import argv
import satlas as s 
import pandas as pd 


directory_1 = 'D:\\Agi\\Documents\\PhD\\CRIS\\K\\K data\\K_Online_Results\\'

class Analysis(object):
	def __init__ (self, directory):
		self.all_data=dict()
		for mass in os.listdir(directory):
			
			self.all_data[mass]={}
			subdir = os.path.join(directory,'{}\\Results_Diode_Correction'.format(mass))
			for filename in os.listdir(subdir):
				if filename.endswith(".txt"):
					scan_no = filename.strip('.txt')
					data = pd.read_csv(os.path.join(subdir,filename), delimiter = '\t',  header=None, names = ['parameter', 'value', 'error'])
					self.all_data[mass][scan_no]= {}
					for index, row in data.iterrows():
						self.all_data[mass][scan_no][row['parameter']]=[row ['value'], row ['error']]


	def get_value(self, mass, scan_no,parameter):
		
		param = self.all_data[mass][scan_no][parameter][0]
		error = self.all_data[mass][scan_no][parameter][1]

		return param,error


	def get_list(self, mass, parameter):
		list_of_params=dict()
		

		for scan,par_dict in self.all_data[mass].items():
			print(scan, par_dict)
			param = par_dict[parameter]
			list_of_params [scan] = [param[0], param[1]]

		return list_of_params

	def get_scan(self,scan):
		
		for mass_no, scan_dict in self.all_data.items():
			try:
				scan_params = scan_dict[scan]
				break
			except:
				pass
			
		return scan_params



a = Analysis(directory_1)
b=a.get_list('49', 'Centroid')

# print(b)
masses=[]
for mass, cent in b.items():
	agi=cent[0]
	#print(mass, agi, cent[1])
	plt.plot(mass, agi, 'bo', alpha=0.5)
	plt.errorbar(x=[float(mass)], y=[float(agi)], yerr=[float(cent[1])], fmt='r')
	# plt.axhline(y=411.8, color='k', linestyle='-', alpha=0.5)
	# plt.axhline(y=411.8+0.2, color='k', linestyle='--', alpha=0.5)
	# plt.axhline(y=411.8-0.2, color='k', linestyle='--', alpha=0.5)

plt.show()

# plt.axhline(y=0.0, color='r', linestyle='-')
# plt.show()



