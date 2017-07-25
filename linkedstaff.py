#Analysis of even_even isotopes for the bachelor project
# focus on isotope shifts

OFFSET = 389286074.578447#-384228152.36085065
import os
import numpy as np
import matplotlib.pyplot as plt 
from sys import argv
import satlas as s 
import pandas as pd 

mass =38
c = 299792458

list_to_merge=[68,69]

def doppler(V,m):
	beta = np.sqrt(1-((m* 931.494095*10**6)**2/(V+m* 931.494095*10**6)**2))
	dop = np.sqrt((1-beta)/(1+beta))
	return dop

directory= 'D:\\Agi\\Documents\\PhD\\CRIS\\K\\K data\\K_online_binned\\{}\\'.format(mass) 
directory_save = 'D:\\Agi\\Documents\\PhD\\CRIS\\K\\K data\\K_Online_Results\\{}\\Results\\'.format(mass) 

if not os.path.exists(directory_save):
    os.makedirs(directory_save)

Info = pd.read_csv ('C:\\Users\\MyStuff\\OneDrive\\Documents\\CRIS\\K\\Data_exl.csv')
Spins = Info.SPIN

for i in range (len(Spins)):
	isotope= int(Info.K[i]) 
	if isotope==mass:
		I_info = Info.SPIN[i]
		m= float(Info.MASS[i])

		Isotope_shif = float(Info.Isotope_shift[i])+600
		A_L =float(Info['A_Lower'][i])
		print(A_L)
		B_L =float(0)
		A_U =float( Info['A_upper'][i])
		print(A_U)
		B_U =0



spectra = []
xs = []
ys = []
yerrs = []

for filename in os.listdir(directory):
	if filename.endswith(".txt"):#  and  scan in filename:


		data = np.loadtxt(directory + filename, skiprows=1, delimiter = ',')
		scan_no = filename.strip('.txt')
		
		if int(scan_no) in list_to_merge:



			wavelength=data[5:,2]
			wavelength_err=data[5:, 3]
			diode_laser=data[5:,4]
			diode_laser_err=data[5:,5]


			counts=data[5:,13]
			n = data[5:,6]
			rate= counts/n
			#v_iscool = 9960*np.ones(len(wavelength))
		
			v_iscool = data[5:,7]
			#diodes = data[5:,8]
			#diodes_err = data[5:,9]


			yerr = np.sqrt(counts)
			yerr[yerr == 0] = 1
			yerr = yerr / n


			freq= c*wavelength*10**(-4)
			frequency  = freq * doppler(v_iscool, m)
			#frequency=frequency-(diode_laser*c*10**(-4))

			#print(frequency-OFFSET)

			I=I_info
			J=[0.5,0.5]
			ABC=[A_L, A_U, B_L, B_U,0, 0]

			#centroid fro 226
			centroid = Isotope_shif
			frequency -= OFFSET

			#for 226
			#centroid=4.19696217*10**8
			#fwhm=100
			scale=np.max(rate)
			bkg=[rate[2]]
			scale -= bkg[0]


			fwhm=100
		
			spectrum = s.HFSModel(I=I,J=J,ABC=ABC,centroid=centroid,fwhm = fwhm, 
			background_params = bkg,
			scale = scale, use_racah=False

			,shape='asymmlorentzian', 
			asymmetryparams = {'a':-0.009})#, shape = 'lorentzian')

			spectrum.set_variation({'Al':True,'Au':True,'Bl':False,'Bu':False,'Cl':False,'Cu':False, 'Background':True, 'FWHML':True, 'FWHMG':True})
			#spectrum.set_value({'Amp5_2__3_2':0.01})


			spectra.append(spectrum)
			xs.append(frequency)
			ys.append(rate)
			yerrs.append(yerr)

scan_no=list_to_merge[0]
linked_spectrum = s.LinkedModel(spectra)
# linked_spectrum.models[0].set_variation({'Scale':True})
linked_spectrum.shared = ['Al','Au','Centroid', 'Amp5_2__5_2', 'Amp5_2__3_2','Amp3_2__5_2' ,'Amp3_2__3_2' ,'FWHMG', 'FWHML', 'Background']
#'Amp3_2__3_2','Amp5_2__3_2','Amp3_2__5_2' ,'Amp5_2__5_2'
## Plot the initial guess to check
fig,ax=linked_spectrum.plot(show=False)
linked_spectrum.plot(x=xs,y=ys,yerr=yerrs,ax=ax)

#do the fit!
s.chisquare_fit(linked_spectrum,x=xs,y=ys,yerr=yerrs)

## show the fit
linked_spectrum.display_chisquare_fit()
print(linked_spectrum.get_result_frame(vary=True))

linked_spectrum.plot(x=xs,y=ys,yerr=yerrs)
#plt.show()
fig.savefig(directory_save+'{}.png'.format(scan_no))
plt.close()

with open(directory_save+'{}.txt'.format(scan_no), 'w') as f:
	for p in linked_spectrum.params.values():
		p.name=p.name.strip('s0_')

		f.write('{}\t{}\t{}\n'.format(p.name,p.value,p.stderr))

# print(filename)