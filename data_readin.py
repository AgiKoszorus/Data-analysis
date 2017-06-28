import numpy as np
import h5py
import pandas as pd
import sys,os
import satlas as sat
import matplotlib.pyplot as plt



mass=49

def get_format(filepath):
    with open(filepath,'r') as reader:
        _ = reader.readline() ## mass
        _ = reader.readline() ## scan
        frmt = eval(reader.readline()) ## format

    return frmt

def calchist(df,binsize=1,binpar=['x', 'time']):
    
    bins = np.arange(df[binpar].min(),
             df[binpar].max()+binsize,binsize)

    return hist(df,bins,binpar)
    
def hist(df,bins,binpar):
    bin_selection = np.digitize(df[binpar], bins)
    df['bin_selection'] = bin_selection
    groups = df.groupby('bin_selection')
    
    binned = pd.DataFrame()
    binned['time'] = groups['time'].agg(np.mean)
    binned[['x','xerr']] = groups['x'].agg([np.mean,np.std])
    binned[['x2','xerr2']] = groups['x2'].agg([np.mean,np.std])
    binned['n'] = groups['bin_selection'].sum()/binned.index.values
    binned[['v','verr']] = groups['v'].agg([np.mean,np.std])
    binned[['diodes_diode','diodes_diodeerr']] = groups['diodes_diode'].agg([np.mean,np.std])
    binned[['diodes_m2','diodes_m2err']] = groups['diodes_m2'].agg([np.mean,np.std])

    binned['y'] = groups['y'].sum()
    binned['yerr'] = np.sqrt(binned['y'])
    
    return binned

def extract_and_bin(file_name = None, scans = [], x_bin_size = 1, wave1 = 10000, wave2=10000, wavevar = 50 ):
    if file_name is None:
        file_name = 'D:\\Agi\\Documents\\PhD\\CRIS\\K\\K data\\K_online_data\\{}\\'.format(mass) 

    y_par = 'Counts'

    x_origin,x_par_name = ('wavemeter','wavenumber_1')
    x2_origin,x2_par_name = ('wavemeter','wavenumber_2')
    y_origin,y_par_name = ('cris','Counts')
    v_origin,v_par_name = ('iscool','voltage')
    d_origin,d_par_name = ('diodes','M2_FPI')
    d2_origin,d2_par_name = ('diodes','M2_inj')

    origins = (x_origin,x2_origin,y_origin,v_origin, d_origin, d2_origin)   # removed  ---> v_origin,
    par_names = (x_par_name, x2_par_name,y_par_name, v_par_name, d_par_name, d2_par_name)  # removed ----> v_par_name,
    aliases = ('x','x2','y','v',  'diodes_diode', "diodes_m2")  # removed ----->  'v',


    for scan_number in scans:  

        frames = []
        print(scan_number)
        if int(scan_number) > -1 :#and not int(scan_number) in useless:
            scan_path = "scan_{0:04d}".format(scan_number)
          
            data_sets = []
            for i,o in enumerate(origins):
                _ = "{}\\{}_ds.csv".format(scan_path,o)
                path = os.path.join(file_name,_)
                data_set = np.loadtxt(path,delimiter = ';')

                
                _ = "{}\\metadata_{}_ds.txt".format(scan_path,o)
                path = os.path.join(file_name,_)
                frmt = get_format(path)
                col_index = frmt.index(par_names[i])

                frames.append(pd.DataFrame({'time':data_set[:,0],
                                aliases[i]:data_set[:,col_index]}) )

            data_frame = pd.concat(frames)

            data_frame.sort_values(by = 'time', inplace=True)
            
            for par_name in aliases:
                if not par_name == y_par:
                    data_frame[par_name].fillna(method='ffill', inplace=True)
                    data_frame[par_name].fillna(method='bfill', inplace=True)

            data_frame.dropna(inplace=True)
            data_frame.reset_index(inplace=True)


            if len(data_frame)>0: 
                median = np.median(data_frame['x'].values)
                data_frame = data_frame[np.logical_and(data_frame['x'] > median-wavevar, data_frame['x'] < median + wavevar)]
                data_frame_binned_time = calchist(data_frame, binsize = 1, binpar = 'time')
                data_frame_binned = calchist(data_frame, binsize = x_bin_size, binpar = 'x')
                print(data_frame_binned)
                directory_save = 'D:\\Agi\\Documents\\PhD\\CRIS\\K\\K data\\K_online_binned\\{}\\'.format(mass)
                if not os.path.exists(directory_save):
                    os.makedirs(directory_save)

                plt.plot(data_frame_binned['x'],data_frame_binned['y'], 'ro')
                plt.savefig(directory_save+'{}.png'.format(scan_number))
                #plt.show()
           
                plt.close()          

        # directory = 'C:\\Users\\User\\OneDrive\\Documents\\CRIS\\K_experiment\\'
        # if not os.path.exists(directory):
        #     os.makedirs(directory)


        data_frame_binned.to_csv(directory_save+'{}.txt'.format(scan_number))


        # data_frame['time'].to_csv(directory+'{}.txt'.format(scan_number))


dsats = extract_and_bin(file_name = None, scans = [ 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 53, 54, 55, 56, 57], x_bin_size =0.0004, wave1 = 13003 ,  wave2 =12816, wavevar = 10)
#54, 53, 52, 46, 45, 43,41, 40,39, 38,37,36,35