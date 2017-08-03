import os
import numpy as np
import matplotlib.pyplot as plt 
from sys import argv
import satlas as s 
import pandas as pd 

#masses = [ 39, 40, 41, 42, 43, 44, 45, 46, 47, 48]
#for i in range(len(masses)):
	#mass = masses[i]
c = 299792458

mass =38
scan_no_1 = 169
OFFSET = 389286074.578447-384228152.36085065


def doppler(V,m):
	beta = np.sqrt(1-((m* 931.494095*10**6)**2/(V+m* 931.494095*10**6)**2))
	dop = np.sqrt((1-beta)/(1+beta))
	return dop

scan_info= pd.read_csv ('C:\\Users\\MyStuff\\OneDrive\\Documents\\CRIS\\K\\K_{}_fullscans.txt'.format(mass))
diode_info= pd.read_csv ('C:\\Users\\MyStuff\\OneDrive\\Documents\\CRIS\\K\\No_diode_correction.txt'.format(mass))
low_res_scans_list =  pd.read_csv ('C:\\Users\\MyStuff\\OneDrive\\Documents\\CRIS\\K\\low_resolution_scans.txt')
new_start_end =  pd.read_csv ('C:\\Users\\MyStuff\\OneDrive\\Documents\\CRIS\\K\\special_start_end.txt', delimiter = '\t')


Info = pd.read_csv ('C:\\Users\\MyStuff\\OneDrive\\Documents\\CRIS\\K\\Data_exl.csv')

directory= 'D:\\Agi\\Documents\\PhD\\CRIS\\K\\K data\\K_online_binned\\{}\\'.format(mass) 
directory_save = 'D:\\Agi\\Documents\\PhD\\CRIS\\K\\K data\\K_online_results\\{}\\Results_corrected_diode\\'.format(mass) 

list_no_diode_correction = diode_info['no_correction'].values
list_high_res =scan_info['full_scans'].values
low_resolution_scans = low_res_scans_list['Low_res'].values

special_scan = new_start_end['scan'].values

if not os.path.exists(directory_save):
    os.makedirs(directory_save)

Spins = Info.SPIN

for i in range (len(Spins)):
	isotope= int(Info.K[i]) 
	if isotope==mass:
		I_info = Info.SPIN[i]
		m= float(Info.MASS[i])

		Isotope_shif = float(Info.Isotope_shift[i])+550# +2900

		A_L =float(Info['A_Lower'][i])
		B_L =float(0)
		A_U =float( Info['A_upper'][i])
		B_U =0

for filename in os.listdir(directory):
	if filename.endswith(".txt"): 


		data = np.loadtxt(directory + filename, skiprows=1, delimiter = ',')
		scan_no = int(filename.strip('.txt'))
		if  not scan_no in list_no_diode_correction and scan_no==scan_no_1: 

			for_length=data[:,2]


			if scan_no in special_scan:
				start = new_start_end[new_start_end['scan'].values==scan_no]['start'].values
				print(start)
				end = new_start_end[new_start_end['scan'].values==scan_no]['stop'].values
				print(end)
				if end ==0:
					end =len(for_length)+1
					print('new end', end)

				start = 0
				stop = len(for_length)+1
				
				wavelength=data[start:stop,2]
				wavelength_err=data[start:stop, 3]
				diode_laser=data[start:stop,4]
				diode_laser_err=data[start:stop,5]
				counts=data[start:stop,13]
				n = data[start:stop,6]
				rate= counts/n
				v_iscool = data[start:stop,7]
				yerr = np.sqrt(counts)
				yerr[yerr == 0] = 1
				yerr = yerr / n*2.3
				freq= c*wavelength*10**(-4)
				frequency  = freq * doppler(v_iscool, m)
				frequency= frequency-(diode_laser*c*10**(-4))
		
				I=I_info
				J=[0.5,0.5]
				ABC=[A_L, A_U, B_L, B_U,0, 0]

				centroid = Isotope_shif
				frequency -= OFFSET

				scale=np.max(rate)
				bkg=[rate[0]]
				scale -= bkg[0]

				if scan_no in low_resolution_scans:
					fwhm = 60 	
					spectrum_asym = s.HFSModel(I=I,J=J,ABC=ABC,centroid=centroid,fwhm = fwhm, 
					                      background_params =[0.001, 0.00000005], scale = scale, 
					                      use_racah=False
					                      ,shape='asymmlorentzian', 
					                      asymmetryparams = {'a':-0.009})
					                      #, shape = 'voigt')
					spectrum = spectrum_asym
					print('asym')
				elif int(scan_no) in list_high_res:
					fwhm=[30,30]
					spectrum_sym = s.HFSModel(I=I,J=J,ABC=ABC,centroid=centroid,fwhm = fwhm, 
				                      background_params =[0.001, 0.00000005], scale = scale, 
				                      use_racah=False
				                      #,shape='asymmlorentzian', 
				                      #asymmetryparams = {'a':-0.009})
				                      , shape = 'voigt')
					spectrum = spectrum_sym
					print('sym')
				else:
					print('nothing')
					pass

				# spectrum.set_variation({'Cl':False,'Cu':False, 'FWHMG':True, 'FWHML':True})
				# spectrum.plot(x=frequency,y=rate,yerr=yerr)
				# #do the fit!
				# s.chisquare_fit(spectrum,x=frequency,y=rate,yerr=yerr)

				# ## show the fit
				# spectrum.display_chisquare_fit()
				# print(spectrum.get_result_frame(vary=True))


				# fig, (ax1, ax2) = plt.subplots(2,1, sharex=True)
				# ax1.errorbar(x=frequency,y=rate,yerr=yerr, fmt='bo')
				# ax1.plot(frequency,spectrum(frequency), 'g-', linewidth=2.5)
				# ax1.set_xlabel('MHz')
				# ax1.set_ylabel('Rate')
				# ax2.errorbar(x=frequency,y=(rate-spectrum(frequency))/yerr,yerr=yerr, fmt='ro')
				# ax1.set_title('{}'.format(filename))
				
				#fig.savefig(directory_save+'{}.png'.format(scan_no))	

#


				# scan_no = filename.strip('.txt')
				# with open(directory_save+'{}.txt'.format(scan_no), 'w') as f:
				# 	for p in spectrum.params.values():
				# 		f.write('{}\t{}\t{}\n'.format(p.name,p.value,p.stderr))
				# 	f.write('{}\t{}\n'.format('dof',spectrum.ndof_chi))
				# 	f.write('{}\t{}\n'.format('chisqr',spectrum.chisqr_chi))
				# 	f.write('{}\t{}\n'.format('redchisqr',spectrum.redchi_chi))