import os
import numpy as np
import matplotlib.pyplot as plt 
from sys import argv
import satlas as s 

'''Fitting hyperfine structures
   CRIS data
   binned
   satlas
   no text file with HFS information''' 


# offset is the NIST transition frequency
OFFSET = 389286074.578447#-384230484.4685+572000

#speed of light
c = 299792458


# doppler shifting the frequency point by point
def doppler(V,m):
	beta = np.sqrt(1-((m* 931.494095*10**6)**2/(V+m* 931.494095*10**6)**2))
	dop = np.sqrt((1-beta)/(1+beta))
	return dop

mass = 49
directory= 'D:\\Agi\\Documents\\PhD\\CRIS\\K\\K data\\K_online_binned\\{}\\'.format(mass) 
directory_save = 'D:\\Agi\\Documents\\PhD\\CRIS\\K\\K data\\K_Online_Results\\{}\\Results\\'.format(mass) 
if not os.path.exists(directory_save):
    os.makedirs(directory_save)


m=48.9682108

# number added to the offset is the isotope shift
OFFSET +=710


scan_no_1 = '29'



for filename in os.listdir(directory):
	print(filename)
	if filename.endswith(".txt") and  scan_no_1 in filename:


		data = np.loadtxt(directory + filename, skiprows=1, delimiter = ',')
		scan_no = filename.strip('.txt')
		#if int(scan_no) in list_of_full_scans: 




		wavelength=data[:,2]
		wavelength_err=data[:, 3]
		diode_laser=data[:,4]
		diode_laser_err=data[:,5]


		#wavelength=wavelength-diode_laser
		counts=data[:,13]
		n = data[:,6]
		rate= counts/n
	
		v_iscool = data[:,7]
		#diodes = data[:,8]
		#diodes_err = data[:,9]


		yerr = np.sqrt(counts)
		yerr[yerr == 0] = 1
		yerr = yerr / n


		freq= c*wavelength*10**(-4)
		frequency  = freq * doppler(v_iscool, m)


		I=0.5
		J=[0.5,0.5]
		ABC=[2368, 285, 0, 0,0, 0]


		centroid = 0
		frequency -= OFFSET

		fwhm=[20,20]
		scale=np.max(rate)
		bkg=[rate[0]]
		scale -= bkg[0]

		
		spectrum = s.HFSModel(I=I,J=J,ABC=ABC,centroid=centroid,fwhm = fwhm, background_params = bkg, scale = scale, use_racah=False)#, shape = 'lorentzian')
		spectrum.set_variation({'Cl':False,'Cu':False, 'FWHMG':True, 'FWHML':True})


		## Plot the initial guess to check
		spectrum.plot(x=frequency,y=rate,yerr=yerr)

		#do the fit!
		s.chisquare_fit(spectrum,x=frequency,y=rate,yerr=yerr)

		## show the fit
		spectrum.display_chisquare_fit()
		print(spectrum.get_result_frame(vary=True))

		plt.errorbar(x=frequency,y=rate,yerr=yerr, fmt='bo')
		plt.plot(frequency,spectrum(frequency), 'g-', linewidth=2.5)

		fig.savefig(directory_save+'{}.png'.format(scan_no))	

		



		scan_no = filename.strip('.txt')
		with open(directory_save+'{}.txt'.format(scan_no), 'w') as f:
			for p in spectrum.params.values():
				if p.name == 'Centroid':
					p.value += OFFSET
				
				f.write('{}\t{}\t{}\n'.format(p.name,p.value,p.stderr))