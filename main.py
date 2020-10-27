# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 10:59:25 2020

@author: michi
"""

from __future__ import print_function
from ipywidgets import interact, interactive, fixed, interact_manual
from scipy.signal import argrelextrema

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.image as mpimg

from functions import Orbitload
from functions import load_yTag_data
from functions import load_data
from functions import load_Bfield
from functions import enlargeVector



############## user input ###################

# This section asks the user for the path to the store data sets, 
# the data sets they would like to analyze and the type of proton data set and
# saves the answers in variables

# filepath = interactive(Orbitload, Orbit = '')
# print('Please give the local path to where the data set you would like to analyze 
# 	  'is stored, in the format eg C:/Users/Username/Documents/Git/Jupiters_energetic_particles/data')
# display(filepath)

timespan = interactive(Orbitload, Orbit = '')
print('Which timespan would you like to analyze? Please type in format "Orbitxx/dd-dd_mm_')
display(timespan)

ProtonData = interactive(Orbitload, Orbit = '')
print('Proton dataset')
display(ProtonData)

answer_Orbit = timespan.kwargs['Orbit']
answer_Proton = ProtonData.kwargs['Orbit']
answer_filepath = filepath.kwargs['Orbit']

fileIntE_O = '%sEn_O'%answer_Orbit
fileIntE_S = '%sEn_S'%answer_Orbit
fileIntE_He = '%sEn_He'%answer_Orbit
fileIntE_H = '%sEn_%s_H'%(answer_Orbit, answer_Proton)
fileIntPA_O = '%sPA_O'%answer_Orbit
fileIntPA_S = '%sPA_S'%answer_Orbit
fileIntPA_He = '%sPA_He'%answer_Orbit
fileIntPA_H = '%sPA_%s_H'%(answer_Orbit, answer_Proton)
filedensO = '%sdens_O'%answer_Orbit
filedensS = '%sdens_S'%answer_Orbit
filedensHe = '%sdens_He'%answer_Orbit
filedensH = '%sdens_%s_H'%(answer_Orbit, answer_Proton)


############## parameters ###################

# ion masses
m_O = 16*1.672*(10**(-27)); #[kg]
m_S = 32*1.67*(10**(-27)); #[kg]
m_He = 4*1.67*(10**(-27)); #[kg]
m_H = 1*1.67*(10**(-27)); #[kg]
m = [m_O, m_S, m_He, m_H] #kg

# the data sets will be imported and saved in a list. The number of file header lines
# to drop varies for each data set and is declared here.
header_mag = 1
header_dens = 14
yTag_En = 13
if answer_Proton == 'HiResIon' or answer_Proton == 'LoTOFxE':
    yTag_En_H = 13
if answer_Proton == 'HiTOFxE':
    yTag_En_H = 14
yTag_PA = 14
yTag_to_intens_En_S = 2
yTag_to_intens = 3



############## data import ###################

# data sets are imported using the function load_yTag_data() in functions.py
# the list I_E[i] contains all intensity vs energy data for the 4 species H+, He+, O+ and S+. 
# The list I_PA all intensity vs pitch angle data for the 4 species H+, He+, O+ and S+.
# The list n[i] contains the densities. The first column of each list is the timestamp.

files_E = [[fileIntE_O, yTag_En, yTag_to_intens], [fileIntE_S, yTag_En, yTag_to_intens_En_S], 
         [fileIntE_He, yTag_En, yTag_to_intens], [fileIntE_H, yTag_En_H, yTag_to_intens]]
files_PA = [[fileIntPA_O, yTag_PA, yTag_to_intens], [fileIntPA_S, yTag_PA, yTag_to_intens],
           [fileIntPA_He, yTag_PA, yTag_to_intens], [fileIntPA_H, yTag_PA, yTag_to_intens]]

​
n_spec = len(files_E)
I_E = [[]]*n_spec
I_PA = [[]]*n_spec
for i,species in enumerate(files_E):
    I_E[i] = [load_yTag_data(species[0],species[1],species[2])][0]
for i,species in enumerate(files_PA):
    I_PA[i] = [load_yTag_data(species[0],species[1],species[2])][0]

# import density:
files_n = [filedensO, filedensS, filedensHe, filedensH]
n = [[]]*n_spec
for i,val in enumerate(files_n):
    n[i] = load_data(val, header_dens)
	
# import B-field
B = [[] for k in range(9)]
B = load_Bfield(answer_Orbit, B, header_mag)


############## convert time to seconds ###################

# round minute & second values for allocation in 2.1.2

second_B = [int(round(val,0)) for val in B[4]]
minute_B = [int(round(val,0)) for val in B[3]]
min_sec_vec_B = np.transpose([minute_B, second_B])
min_sec_B = [float(f'{i[0]}.{i[1]}') for i in min_sec_vec_B]


############## Linearisiere B-feld zum Plotten ###################

day_B_vec = [B[1][0]+B[2][0]/24, B[1][-1]+B[2][-1]/24]
day_B_vec = enlargeVector(day_B_vec,len(B[-1]))


n_arr = [[]]*n_spec
temp = []
for i,spec in enumerate(n):
    for val in spec:
        temp.append(val[0])
    n_arr[i] = temp
    temp = []
n_interp = [[]]*n_spec
for i,specarr in enumerate(n_arr):
    n_interp[i] = enlargeVector(specarr,len(B[8]))
	
x_PA = []
I_PA_nofill = [[]]*n_spec
for spec,specArr in enumerate(I_PA):
    for PAval in specArr[1]:
        x_PA.append([x for x in PAval if x != -1.0E38])
    I_PA_nofill[spec] = x_PA
    x_PA = []       
	
PA_full = [[]]*n_spec
PA_temp = []  
for spec,specArr in enumerate(I_PA_nofill):
    for PA,PAarr in enumerate(specArr):
        if len(PAarr) == 30:
            PA_temp.append(PAarr)
    PA_full[spec] = PA_temp
    PA_temp = []
print(len(PA_full[1]),'/',len(I_PA[1][1]), ' Intensity-Daten sind voll PA-abgedeckt.')



############## calculate number of B_above & B_below: ###################

​
maglatExtr = argrelextrema(np.array(B[6]), np.greater_equal)
num = []
i_now = 0
for i,magval in enumerate(maglatExtr[0]):
    if abs(magval - i_now) > 1000:
        num.append(magval)
        i_now = magval
    else:
        continue
n_pass = len(num)

# Anzahl der rel. Extrema plus der erste Wert (index 0) --> +1
num_below = int(np.floor((len(num)+1)/2))
num_above = int(num_below +1)

​

# --> take vector 'num' as passages to to find changepoints in mean with matlab function findchangepts
print(num)

#np.savetxt('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/%s/num.out'%answer_Orbit[0:7], num, delimiter=',')


# changepoints by variation in mean (matlab):

iptload = open('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/%s/ipt.txt'%answer_Orbit[0:7],'r')
ipt = []
for line in iptload:
    line = line.strip()
    columns = line.split()
    for j in columns[0:]:
        ipt.append(float(j))
for i,val in enumerate(ipt):
    ipt[i] = int(val)

ipt_lo = [val for i,val in enumerate(ipt[0:len(ipt):2])]
ipt_hi = [val for i,val in enumerate(ipt[1:len(ipt):2])]


############## Determine B_inside/B_above/B_below ###################
# select percentage(s) by which current sheet thickness will be reduced on each side of the drop:

perc = [0,5,10,20] #in percent
n_perc = len(perc)


# get all 'inside-indices':
Bins_perc_temp = [[[]]*n_pass]*n_perc
Bins_perc_temp_edgeabo = [[[]]*n_pass]*n_perc
Bins_perc_temp_edgebel = [[[]]*n_pass]*n_perc
Bin = [[]]*n_pass
Bins_perc = [] # percentage-Variationen x current-sheet Durchgänge
Bins_perc_edgeabo = []# transition zone from above and inside the current sheet
Bins_perc_edgebel = [] # transition zone from below and inside the current sheet

​
#mean Variation:
for i in range(0,len(ipt_lo),1):
    Bin[i] = [k for k in range(ipt_lo[i],ipt_hi[i]+1,1)]
   

#convolution:
#for i,val in enumerate(passage_inds):
#    Bin[i] = [k for k in range(val[-1],val[-2]+1,1)]


for j,percentage in enumerate(Bins_perc_temp):
    for i,crossing in enumerate(percentage):
        perc_val = round(len(Bin[i])*(perc[j]/100))
        Bins_perc_temp[j][i] = (Bin[i][perc_val:len(Bin[i])-perc_val])
        if i%2 == 0:
            Bins_perc_temp_edgeabo[j][i] = Bin[i][0:perc_val]
            Bins_perc_temp_edgebel[j][i] = Bin[i][len(Bin[i])-perc_val:len(Bin[i])]
        if i%2 != 0:
            Bins_perc_temp_edgebel[j][i] = Bin[i][0:perc_val]
            Bins_perc_temp_edgeabo[j][i] = Bin[i][len(Bin[i])-perc_val:len(Bin[i])]
        Bins_perc.append(Bins_perc_temp[j][i])
        Bins_perc_edgebel.append(Bins_perc_temp_edgebel[j][i])
        Bins_perc_edgeabo.append(Bins_perc_temp_edgeabo[j][i])
Bins_perc = np.reshape(Bins_perc,(n_perc,n_pass))
Bins_perc_edgeabo = np.reshape(Bins_perc_edgeabo,(n_perc,n_pass))
Bins_perc_edgebel = np.reshape(Bins_perc_edgebel,(n_perc,n_pass))


B_inside = []
for arr in Bins_perc:
    B_inside.append(np.concatenate(arr, axis=None))
    
B_edgeabo = []
for arr in Bins_perc_edgeabo:
    B_edgeabo.append(np.concatenate(arr, axis=None))
    
B_edgebel = []
for arr in Bins_perc_edgebel:
    B_edgebel.append(np.concatenate(arr, axis=None))
    
B_edge = []
for i in range(0,n_perc,1):
    B_edge.append(np.sort(list(B_edgeabo[i]) + list(B_edgebel[i])))



# calculate indices of above and below current sheet for all percentage variations 
#and current sheet crossings and concatenate to B_above & B_below:


#Above:

B_abs_temp = []
B_abs = [[[]]*(num_above-2) for i in range(n_perc)]
B_abs1 = [[]]*n_perc
B_abs2 = [[]]*n_perc
B_above = [[]]*n_perc

# calculate before 1st and after last crossing:
	
for i,val in enumerate(perc):
    B_abs1[i] = [k for k in range(0,Bins_perc[i][0][0],1)]
    B_abs2[i] = [k for k in range(Bins_perc[i][len(num)-1][-1]+1,len(B[-1]),1)]

# calculate above-indices between all other crossings:    

for i,val in enumerate(perc):    
    for j,arr in enumerate(B_abs[i]):
        for k in range(Bins_perc[i][j*2+1][-1]+1,Bins_perc[i][j*2+2][0],1):
            B_abs_temp.append(k)
        B_abs[i][j] = B_abs_temp
        B_abs_temp = []

#concatenate to 1 vector
#all indices above current sheet for i = percentage variations

for i,val in enumerate(perc):    
    B_above[i] = B_abs1[i] + list(np.concatenate(B_abs[i])) + B_abs2[i]    

#Below:

B_bes_temp = []
B_bes = [[[]]*(num_below) for i in range(n_perc)]
B_below = [[]]*n_perc


# calculate below-indices inbetween all crossings:    
for i,val in enumerate(perc):    
    for j,arr in enumerate(B_bes[i]):
        for k in range(Bins_perc[i][j*2][-1]+1,Bins_perc[i][j*2+1][0],1):
            B_bes_temp.append(k)
        B_bes[i][j] = B_bes_temp
        B_bes_temp = []

​
#concatenate to 1 vector
for i,val in enumerate(perc):
    B_below[i] = list(np.concatenate(B_bes[i]))

#concatenate regions to 1 vector

Bclass = [B_inside, B_above, B_below, B_edgeabo, B_edgebel, B_edge]
Nclass = len(Bclass)


# calculate percentages of each region:
N_in = [[]]*n_perc
N_ab = [[]]*n_perc
N_be = [[]]*n_perc
N_edgeabo = [[]]*n_perc
N_edgebel = [[]]*n_perc
N_edge = [[]]*n_perc

​
for i in range(n_perc):
    N_in[i] = len(B_inside[i])
    N_ab[i] = len(B_above[i])
    N_be[i] = len(B_below[i])
    N_edgeabo[i] = len(B_edgeabo[i])
    N_edgebel[i] = len(B_edgebel[i])
    N_edge[i] = len(B_edge[i])

print('B_inside: ',round(N_in[0]/len(B[8]),4)*100,'%')
print('B_above: ',round(N_ab[0]/len(B[8]),4)*100,'%')
print('B_below: ',round(N_be[0]/len(B[8]),4)*100,'%')
print('B_edge_abo: ',round(N_edgeabo[3]/len(B[8]),4)*100,'%')
print('B_edge_bel: ',round(N_edgebel[3]/len(B[8]),4)*100,'%')


############## Calculation of Intensity average for inside, above/below (Quasi-Interpolation) ###################


# get time vector in intensity matrix

time_vec = open('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/%s.d2s'%fileIntPA_O, 'r')

for line in range(17):
    line = time_vec.readline()
    
# split columns, float-values, units
time_Int=[]
for line in time_vec:
    line = line.strip()
    columns = line.split()
    for i in columns[0:1]:
        time_Int.append(i)
time_vec.close()



# separate time vector in days, hours, minutes, seconds and find time in B-field vector closest 
# to each time of Intensity matrix

date_Int = [dt.datetime(int(i[4:8]),int(i[9:11]),int(i[12:14])) for i in time_Int]
day_Int = [int(i.strftime("%j")) for i in date_Int]
hourInt = [int(i[15:17]) for i in time_Int]
minInt = [int(i[18:20]) for i in time_Int]
secInt = [int(i[21:23]) for i in time_Int]
min_sec_vec_Int = np.transpose([minInt,secInt])
min_sec_Int = [float(f'{i[0]}.{i[1]}') for i in min_sec_vec_Int]
Int_time = np.transpose([day_Int,hourInt,min_sec_Int])


# allocate closest times of both arrays by hours, minutes and seconds
# the array 'seconds_allocate_Int_B' gives list of B-field-time Indices that are closest to each intensity-time value

seconds_allocate_Int_B = []
index = 0
for time in Int_time:
    index +=1
    print('{} %'.format(index/20),end='\r')
    dates = []
    dates.append([lines for lines,val in enumerate(B[1]) if val == time[0]])
    hours = []
    hours.append([lines for lines,val in enumerate(B[2][dates[0][0]:dates[0][-1]]) if val == time[1]])
    hours[0] = [x+dates[0][0] for x in hours[0]]
    seconds_allocate_Int_B.append(min(range(len(min_sec_B[hours[0][0]:hours[0][-1]])), 
                                      key=lambda i:abs(min_sec_B[hours[0][0]:hours[0][-1]][i]-time[2])))
    seconds_allocate_Int_B[-1] = seconds_allocate_Int_B[-1] + hours[0][0]




##############################################################################
# Allocate B-field
# Erklärung der Datenstruktur:

# die Listen sind auf folgender Weise verschachtelt:
# Liste[Region][Verschiebungs-Prozentsatz für Ränder][Species][EnergieKanal/PitchWinkel][Vektor]
# Die Nummerierung ist:

#         Region: 0 = inside, 1 = above, 2 = below, 3 = Randregion_{above},
#		  4 = Randregion_{below}, 5 = Randregion_{gesamt}
#         [Verschiebungs-Prozentsatz für Ränder]: Prozentsatz, um den die Zeit 
#		  innerhalb der current sheet an beiden Rändern jeweils gekürzt wurde um die 
#		  Randregionen genauer zu betrachten. Die Randregionen sind also 2Prozentsatzcurrentsheet_crossingtime breit
#         [Species]: 0 = Oxygen, 1 = Sulphur, 2 = Helium, 3 = Protons
#         [EnergieKanal/PitchWinkel]: Energie je nach Species: Oxygen(0-5), 
#		  Sulphur(0-4), Helium(0-2), Protons depending on dataset. PitchWinkel: 0-29



# Get matrix of intensity indices that belong to the times where the spacecraft 
# is inside/above/below the current sheet:

Int = [[[]]*n_perc]*Nclass
Int_class = [[]]*n_perc
Int_temp = []

for clas,Bclassif in enumerate(Bclass):
    print(clas)
    for ind,p in enumerate(Bclassif):
        index1 = 0
        for valB in p:
            index1 +=1
            print('{:1%}'.format(index1/N_in[ind]),end ='\r')
            for pos,val in enumerate(seconds_allocate_Int_B):
                if valB == val:
                    Int_temp.append(pos)
        Int_class[ind] = Int_temp
        Int_temp = []
    Int[clas] = Int_class
    Int_class = [[]]*n_perc


Int_E = [[[[]]*n_spec for i in range(n_perc)] for i in range(Nclass)]
Int_PA = [[[[]]*n_spec for i in range(n_perc)] for i in range(Nclass)]

x_E = [[[]]*n_spec for i in range(n_perc)]
x_PA = [[[]]*n_spec for i in range(n_perc)]
x_E_temp = []
x_PA_temp = []

for clas,classif in enumerate(Int_E):
    for percind,percArr in enumerate(classif):
        for spec,specArr in enumerate(percArr):
            for j in Int[clas][percind]:
                x_E_temp.append(I_E[spec][1][j])
                x_PA_temp.append(I_PA[spec][1][j])
            x_E[percind][spec] = x_E_temp
            x_PA[percind][spec] = x_PA_temp
            x_E_temp = []
            x_PA_temp = []
    Int_E[clas] = x_E
    Int_PA[clas] = x_PA
    x_E = [[[]]*n_spec for i in range(n_perc)]
    x_PA = [[[]]*n_spec for i in range(n_perc)]


# transpose to delete fillValues:
for clas,classif in enumerate(Int_E):
    for p,percarr in enumerate(classif):
        for i,val in enumerate(percarr):
            Int_E[clas][p][i] = np.transpose(Int_E[clas][p][i])
            Int_PA[clas][p][i] = np.transpose(Int_PA[clas][p][i])


# delete fillValues
Int_E_nofill = [[[[]]*n_spec for i in range(n_perc)] for i in range(Nclass)]
Int_PA_nofill = [[[[]]*n_spec for i in range(n_perc)] for i in range(Nclass)]

​
x_E = [[[]]*n_spec for i in range(n_perc)]
x_PA = [[[]]*n_spec for i in range(n_perc)]
x_E_temp = []
x_PA_temp = []

for clas,classif in enumerate(Int_E):
    for percind,percArr in enumerate(classif):
        for spec,specArr in enumerate(percArr):
            for Eval in specArr:
                x_E_temp.append([x for x in Eval if x != -1.0E38])
            x_E[percind][spec] = x_E_temp
            x_E_temp = []
    Int_E_nofill[clas] = x_E
    x_E = [[[]]*n_spec for i in range(n_perc)]


for clas,classif in enumerate(Int_PA):
    for percind,percArr in enumerate(classif):
        for spec,specArr in enumerate(percArr):
            for PAval in specArr:
                x_PA_temp.append([x for x in PAval if x != -1.0E38])
            x_PA[percind][spec] = x_PA_temp
            x_PA_temp = []
    Int_PA_nofill[clas] = x_PA
    x_PA = [[[]]*n_spec for i in range(n_perc)]



# take average

Int_E_mean = [[[[]]*n_spec for i in range(n_perc)] for i in range(Nclass)]
Int_PA_mean = [[[[]]*n_spec for i in range(n_perc)] for i in range(Nclass)]

x_E = [[[]]*n_spec for i in range(n_perc)]
x_PA = [[[]]*n_spec for i in range(n_perc)]
x_E_temp = []
x_PA_temp = []

for clas,classif in enumerate(Int_E_nofill):
    for percind,percArr in enumerate(classif):
        for spec,val in enumerate(percArr):
            x_E_temp.append([np.mean(x) for x in val])
        x_E[percind] = x_E_temp
        x_E_temp = []
    Int_E_mean[clas] = x_E
    x_E = [[[]]*n_spec for i in range(n_perc)]
    
for clas,classif in enumerate(Int_PA_nofill):
    for percind,percArr in enumerate(classif):
        for spec,val in enumerate(percArr):
            x_PA_temp.append([np.mean(x) for x in val])
        x_PA[percind] = x_PA_temp
        x_PA_temp = []
    Int_PA_mean[clas] = x_PA
    x_PA = [[[]]*n_spec for i in range(n_perc)]

I_E_mean = []
for i in range(4):
    #print(i)
    I_E_mean.append(np.mean(I_E[i][0]))
I_E_mean


############## distribution function ###################

#number density
n_cbm = np.array(n)*1E6 #1/m^3

​
#phase-space density
# müsste intergriert n ergeben, aber nur bei integration über gesamte verteilungsfunktion. 
#Hier wird nicht die gesamte Verteilungsfunktion über E observed

f_E = [[[[]]*n_spec for i in range(n_perc)] for i in range(Nclass)] # in s^3/m^6
f_PA = [[[[]]*n_spec for i in range(n_perc)] for i in range(Nclass)] # in s^3/m^6


x_E = [[[]]*n_spec for i in range(n_perc)]
x_PA = [[[]]*n_spec for i in range(n_perc)]
x_E_temp = []
x_PA_temp = []

for clas,classif in enumerate(Int_E_mean):
    for percind,percArr in enumerate(classif):
        for spec,val in enumerate(percArr):
            if val == []:
                break
            else:
                x_E_temp.append(np.divide(np.array(val)*10000/(1.6*(10**(-19))*1000)*(m[spec]**2),
                                      2*I_E[spec][0]*1.6*(10**(-19))*1000))

        x_E[percind] = x_E_temp
        x_E_temp = []
    f_E[clas] = x_E
    x_E = [[[]]*n_spec for i in range(n_perc)]

       
for clas,classif in enumerate(Int_PA_mean):
    for percind,percArr in enumerate(classif):
        for spec,val in enumerate(percArr):
            if val == []:
                break
            else:
                x_PA_temp.append(np.divide(np.array(val)*10000/(1.6*(10**(-19))*1000)*(m[spec]**2),
                                      2*I_E_mean[spec]*1.6*(10**(-19))*1000))
        x_PA[percind] = x_PA_temp
        x_PA_temp = []
    f_PA[clas] = x_PA
    x_PA = [[[]]*n_spec for i in range(n_perc)]

for region in range(3):
    for species in range(4):
        np.savetxt('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/f_E/f_E%s_%s_0_%s.txt'
				   %(answer_Orbit[5:7],region,species),f_E[region][0][species], delimiter = ',')

for species in range(4):
    np.savetxt('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/f_E/n%s_%s.txt'
			   %(answer_Orbit[5:7],species),n_cbm[species], delimiter = ',')
    np.savetxt('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/f_E/I_E%s_%s_0.txt'
			   %(answer_Orbit[5:7],species),I_E[species][0], delimiter = ',')



############## PA-Count for S/N ###################

# load in counts for PA
yTag_PA_counts = 13

​
fileIntPA_O_counts = '%sPA_O_counts'%answer_Orbit
fileIntPA_S_counts = '%sPA_S_counts'%answer_Orbit
fileIntPA_He_counts = '%sPA_He_counts'%answer_Orbit
fileIntPA_H_counts = '%sPA_%s_H_counts'%(answer_Orbit, answer_Proton)

files_PA_counts = [[fileIntPA_O_counts, yTag_PA_counts, yTag_to_intens], [fileIntPA_S_counts, yTag_PA_counts, yTag_to_intens],
           [fileIntPA_He_counts, yTag_PA_counts, yTag_to_intens], [fileIntPA_H_counts, yTag_PA_counts, yTag_to_intens]]

​
counts_PA = [0,0,0,0]

for i,species in enumerate(files_PA_counts):
    counts_PA[i] = [load_yTag_data(species[0],species[1],species[2])][0]
PA_counts = [[[[]]*n_spec for i in range(n_perc)] for i in range(Nclass)]


x_PA = [[[]]*n_spec for i in range(n_perc)]
x_PA_temp = []
  
for clas,classif in enumerate(PA_counts):
    for percind,percArr in enumerate(classif):
        for spec,specArr in enumerate(percArr):
            for j in Int[clas][percind]:
                x_PA_temp.append(counts_PA[spec][1][j])
            x_PA[percind][spec] = x_PA_temp
            x_PA_temp = []
    PA_counts[clas] = x_PA
    x_PA = [[[]]*n_spec for i in range(n_perc)]

for clas,classif in enumerate(PA_counts):
    for p,percArr in enumerate(classif):
        for i,val in enumerate(percArr):
            PA_counts[clas][p][i] = np.transpose(PA_counts[clas][p][i])
PA_nofill_counts = [[[[]]*n_spec for i in range(n_perc)] for i in range(Nclass)]


x_PA = [[[]]*n_spec for i in range(n_perc)]
x_PA_temp = []    

for clas,classif in enumerate(PA_counts):
    for percind,percArr in enumerate(classif):
        for spec,specArr in enumerate(percArr):
            for PAval in specArr:
                x_PA_temp.append([x for x in PAval if x != -1.0E38])
            x_PA[percind][spec] = x_PA_temp
            x_PA_temp = []
    PA_nofill_counts[clas] = x_PA
    x_PA = [[[]]*n_spec for i in range(n_perc)]



#############################################################################

############## standard deviation = (⎯⎯√∑(𝑐𝑜𝑢𝑛𝑡𝑠)) & 
############## variance = ∑(𝑐𝑜𝑢𝑛𝑡𝑠) (Poisson statistic) ###################

PA_sum_counts = [[[[]]*n_spec for i in range(n_perc)] for i in range(Nclass)]

x_PA = [[[]]*n_spec for i in range(n_perc)]
x_PA_temp = []
   
for clas,classif in enumerate(PA_nofill_counts):
    for percind,percArr in enumerate(classif):
        for spec,val in enumerate(percArr):
            x_PA_temp.append([np.sum(x) for x in val])
        x_PA[percind] = x_PA_temp
        x_PA_temp = []
    PA_sum_counts[clas] = x_PA
    x_PA = [[[]]*n_spec for i in range(n_perc)]
PA_sqrt_counts = [[[[]]*n_spec for i in range(n_perc)] for i in range(Nclass)]​


x_PA = [[[]]*n_spec for i in range(n_perc)]
x_PA_temp = []
  
for clas,classif in enumerate(PA_sum_counts):
    for percind,percArr in enumerate(classif):
        for spec,val in enumerate(percArr):
            x_PA_temp.append([np.sqrt(x) for x in val])
        x_PA[percind] = x_PA_temp
        x_PA_temp = []
    PA_sqrt_counts[clas] = x_PA
    x_PA = [[[]]*n_spec for i in range(n_perc)]



############## Rose criterion ###################

PA_vec = np.arange(3,180,6)
Int_PA_mean_gr5 = [[[[]]*n_spec for i in range(n_perc)] for i in range(Nclass)]
I_PA_gr5 = [[[[]]*n_spec for i in range(n_perc)] for i in range(Nclass)]


x_PA = [[[]]*n_spec for i in range(n_perc)]
i_PA = [[[]]*n_spec for i in range(n_perc)]
x_PA_temp = []
i_PA_temp = []

   
for clas,classif in enumerate(PA_sqrt_counts):
    for percind,percArr in enumerate(classif):
        for spec,specArr in enumerate(percArr):
            for i,PAval in enumerate(specArr):
                if PAval >= 5:
                    x_PA_temp.append(Int_PA_mean[clas][percind][spec][i])
                    i_PA_temp.append(PA_vec[i])
            x_PA[percind][spec] = x_PA_temp
            i_PA[percind][spec] = i_PA_temp
            x_PA_temp = []
            i_PA_temp = []
    Int_PA_mean_gr5[clas] = x_PA
    I_PA_gr5[clas] = i_PA
    x_PA = [[[]]*n_spec for i in range(n_perc)]
    i_PA = [[[]]*n_spec for i in range(n_perc)]


############## Save variables ###################

# save Intensity (energy):
regions = [0,1,2]

for region in regions:
    for species in range(4):
        np.savetxt('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/Orbit%s/save/Int_E_mean%s_%s_0_%s.txt'%
                   (answer_Orbit[5:7],answer_Orbit[5:7],region,species),Int_E_mean[region][0][species], delimiter = ',')

for species in range(4):
    np.savetxt('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/Orbit%s/save/Int_E_mean%s_5_3_%s.txt'%
               (answer_Orbit[5:7],answer_Orbit[5:7],species),Int_E_mean[5][3][species], delimiter = ',')        


# save Intensity (PA):
for region in regions:
    for species in range(4):
        np.savetxt('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/Orbit%s/save/Int_PA_mean_gr5_%s_%s_0_%s.txt'%
                   (answer_Orbit[5:7],answer_Orbit[5:7],region,species),Int_PA_mean_gr5[region][0][species], delimiter = ',')
        np.savetxt('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/Orbit%s/save/I_PA_gr5_%s_%s_0_%s.txt'%
                   (answer_Orbit[5:7],answer_Orbit[5:7],region,species),I_PA_gr5[region][0][species], delimiter = ',')

        
regions_edge = [3,4,5]
for region in regions_edge:
    for species in range(4):
        np.savetxt('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/Orbit%s/save/Int_PA_mean_gr5_%s_%s_3_%s.txt'%
                   (answer_Orbit[5:7],answer_Orbit[5:7],region,species),Int_PA_mean_gr5[region][3][species], delimiter = ',')
        np.savetxt('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/Orbit%s/save/I_PA_gr5_%s_%s_3_%s.txt'%
                   (answer_Orbit[5:7],answer_Orbit[5:7],region,species),I_PA_gr5[region][3][species], delimiter = ',')



############## Einlesen aller Orbits & Plotten ###################

n_orb = 12
n_spec = 3
regions = [0,1,2]
regions_edge = [3,4,5]


# Einladen Intensity (E)
Int_E_mean_orb = [[[[]]*n_spec for i in range(3)] for i in range(n_orb)]
Int_temp = [[[]]*n_spec for i in range(3)]
temp = []
for orb in range(n_orb):
    for region in regions:
        for species in range(n_spec):
            if orb <= 8:
                load = open('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/Orbit0%s/save/Int_E_mean0%s_%s_0_%s.txt'%
                            (orb+1,orb+1,region,species),'r')
            else:
                load = open('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/Orbit%s/save/Int_E_mean%s_%s_0_%s.txt'%
                            (orb+1,orb+1,region,species),'r')
            for line in load:
                temp.append(float(line))
            Int_temp[region][species] = temp
            temp = []
    Int_E_mean_orb[orb] = Int_temp
    Int_temp = [[[]]*n_spec for i in range(3)]

    
    
Int_E_mean_perc_orb = [[[]]*n_spec for i in range(n_orb)]
Int_temp = [[]]*n_spec
temp2 = []
for orb in range (n_orb):
    for species in range(n_spec):
        if orb <= 8:
            load = open('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/Orbit0%s/save/Int_E_mean0%s_5_3_%s.txt'%
                            (orb+1,orb+1,species),'r')
        else:
            load = open('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/Orbit%s/save/Int_E_mean%s_5_3_%s.txt'%
                            (orb+1,orb+1,species),'r')
        for line in load:
            temp2.append(float(line))
        Int_temp[species] = temp2
        temp2 = []
    Int_E_mean_perc_orb[orb] = Int_temp
    Int_temp = [[]]*n_spec


# Einladen Intensity (PA)

Int_PA_mean_gr5_orb = [[[[]]*n_spec for i in range(3)] for i in range(n_orb)]
I_PA_gr5_orb = [[[[]]*n_spec for i in range(3)] for i in range(n_orb)]
Int_temp = [[[]]*n_spec for i in range(3)]
I_temp = [[[]]*n_spec for i in range(3)]
temp_Int = []
temp_I = []
for orb in range(n_orb):
    for region in regions:
        for species in range(n_spec):
            if orb <= 8:
                load_Int = open('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/Orbit0%s/save/Int_PA_mean_gr5_0%s_%s_0_%s.txt'%
                            (orb+1,orb+1,region,species),'r')
                load_I = open('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/Orbit0%s/save/I_PA_gr5_0%s_%s_0_%s.txt'%
                            (orb+1,orb+1,region,species),'r')
            else:
                load_Int = open('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/Orbit%s/save/Int_PA_mean_gr5_%s_%s_0_%s.txt'%
                            (orb+1,orb+1,region,species),'r')
                load_I = open('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/Orbit%s/save/I_PA_gr5_%s_%s_0_%s.txt'%
                            (orb+1,orb+1,region,species),'r')
            for line in load_Int:
                temp_Int.append(float(line))
            for line in load_I:
                temp_I.append(float(line))    
            Int_temp[region][species] = temp_Int
            temp_Int = []
            I_temp[region][species] = temp_I
            temp_I = []
    Int_PA_mean_gr5_orb[orb] = Int_temp
    Int_temp = [[[]]*n_spec for i in range(3)]
    I_PA_gr5_orb[orb] = I_temp
    I_temp = [[[]]*n_spec for i in range(3)]

    
    
Int_PA_mean_gr5_perc_orb = [[[[]]*n_spec for i in range(3)] for i in range(n_orb)]
I_PA_gr5_perc_orb = [[[[]]*n_spec for i in range(3)] for i in range(n_orb)]
Int_temp = [[[]]*n_spec for i in range(3)]
I_temp = [[[]]*n_spec for i in range(3)]
temp_Int = []
temp_I = []
for orb in range(n_orb):
    for r,region in enumerate(regions_edge):
        for species in range(n_spec):
            if orb <= 8:
                load_Int = open('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/Orbit0%s/save/Int_PA_mean_gr5_0%s_%s_3_%s.txt'%
                            (orb+1,orb+1,region,species),'r')
                load_I = open('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/Orbit0%s/save/I_PA_gr5_0%s_%s_3_%s.txt'%
                            (orb+1,orb+1,region,species),'r')
            else:
                load_Int = open('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/Orbit%s/save/Int_PA_mean_gr5_%s_%s_3_%s.txt'%
                            (orb+1,orb+1,region,species),'r')
                load_I = open('C:/Users/michi/Documents/GitHub/Jupiters_energetic_particles/Orbit%s/save/I_PA_gr5_%s_%s_3_%s.txt'%
                            (orb+1,orb+1,region,species),'r')
            for line in load_Int:
                temp_Int.append(float(line))
            for line in load_I:
                temp_I.append(float(line))    
            Int_temp[r][species] = temp_Int
            temp_Int = []
            I_temp[r][species] = temp_I
            temp_I = []
    Int_PA_mean_gr5_perc_orb[orb] = Int_temp
    Int_temp = [[[]]*n_spec for i in range(3)]
    I_PA_gr5_perc_orb[orb] = I_temp
    I_temp = [[[]]*n_spec for i in range(3)]
	