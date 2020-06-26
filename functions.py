
# coding: utf-8

# In[ ]:

import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt


def E_or_PA(E_PA):
    return E_PA

def O_or_S(O_S):
    return O_S

def Orbitload(Orbit):
    return Orbit

def load_yTag_data(file, yTag, yTag_to_intens):
    # load data:
    data = open('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/%s.d2s'%file, 'r')
    
    for line in range(yTag):
        line = data.readline()
    
    line = line.strip()
    columns = line.split(',')
    N = len(columns)
    
    if N == 30:
        end = 6
    else:
        end = 7
    
    vec = np.zeros(N)
    for i,val in enumerate(columns[1:N-1]):
        vec[i+1] = float(val)
    vec[0] = float(columns[0][7::])
    vec[-1] = float(columns[N-1][:end])
    
    for headerline in range(yTag_to_intens):
        headerline = data.readline()
    
    intens = []
    for time in data:
        time = time.strip()
        columns2 = time.split()
        for j in columns2[1:]:
            intens.append(float(j))
    intens = np.reshape(intens,(int(len(intens)/N),N))  
    data.close()
    return vec, intens


def load_data(file, header):
    data = open('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/%s.d2s'%file, 'r')

    # drop header
    for line in range(header):
        line = data.readline()

    # split columns, float-values, units
    arr = []
    for line in data:
        line = line.strip()
        columns = line.split()
        for j in columns[1:]:
            arr.append(float(j))
    arr = np.reshape(arr,(int(len(arr)/(len(columns)-1)),len(columns)-1))
    data.close()
    return arr


def load_Bfield(file, arr_data, header):
    data = open('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/%smag'%file, 'r')

    # drop header
    for line in range(header):
        line = data.readline()

    # split columns, float-values, units

    for line in data:
        line = line.strip()
        columns = line.split()
        N = len(columns)
        for i in range(N):
            for xfunc in columns[i:i+1]:
                arr_data[i].append(float(xfunc))
    data.close()
    return arr_data


def find(arr, arr_reference, operator, value):
    if operator == '>':
        ind = [i for i,x in enumerate(arr_reference) if x >= value]
    else:
        ind = [i for i,x in enumerate(arr_reference) if x <= value]
    arr_new = []
    for i in ind:
        arr_new.append(arr[i])
    return arr_new, ind


def enlargeVector(vector, size):
    old_indices = np.arange(0,len(vector))
    new_length = size
    new_indices = np.linspace(0,len(vector)-1,new_length)
    spl = interp1d(old_indices,vector)
    new_array = spl(new_indices)
    return spl(new_indices)


def butter_highpass(cutoff, fs, order):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, fs, order):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y


def fillVal(arr):
    arr = list(np.transpose(arr))
    for i in range(len(arr)):
    	arr[i] = [x for x in arr[i] if x != -1.00E38]
    return arr
    
    
def mean(arr):
    arr_mean = []
    for i in arr:
        arr_mean.append(np.mean(i))
    return arr_mean
 
    
def distribution_function_E(I_mean,E,m):
    f = np.divide(np.array(I_mean)*(m**2),2*np.array(E))
    return f

def distribution_function_PA(I_mean,E,m):
    f = np.array(I_mean)*(m**2)/(2*E)
    return f
