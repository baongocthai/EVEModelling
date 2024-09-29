# -*- coding: utf-8 -*-
"""
Thai Bao Ngoc
9 September 2024
Environmental Modelling, Assignment 1 - Task 2.1
"""
# =============================================================================
# Simple Advection Transport Model
# The model calculates the concentration of specific pollutant in a river, 
# assuming the river as a series of 0-D, completely mixed tank reactors (CSTR)
# =============================================================================
"""
INPUT
- Qi [m3/s] river flow for river node i
STATE VARIABLES
- Ci [g/m3] pollutant concentration in river node i
- Vcso [m3] volume discharged during CSO event
OUTPUT
- EQS_exc [TRUE/FALSE] exceedance of EQS at river node i
PARAMETERS
- EQS [g/m3] Environmental Quality Standard
- N [-] number of river nodes
- C0 [g/m3] pollutant concentration upstream (first node = Lyngby Lake = node Novana_Model_MOELLEAA_DKI_3687 in HIP)
- tsco [s] duration of an archetypical CSO event (tcso = 4.3hr = 15480s, Table 2, pg. 11)
- theta [-] fraction of the total yearly volume discharged by an arcetypical CSO event (theta = 2.31% = 0.0231, Table 2, pg. 11)
"""
#%% Main block
# clean the workspace
from IPython import get_ipython
get_ipython().magic('reset -sf')
# close all figures
import matplotlib.pyplot as plt
plt.close('all')

# import all libraries
import pandas as pd
import numpy as np
import os

# set working directory
# directory = os.getcwd()
directory = r'C:\Users\ngoc2\Danmarks Tekniske Universitet\Environmental Modelling - General\Step 2\code'
os.chdir(directory)

# import data
CSOdata_file = 'input/Distance_fromCSO_toRiverNode.csv'
CSOdata = pd.read_csv(CSOdata_file) #hubName = River Node ID closest to CSO
RiverData_file = 'input/river_input_data.csv'
RiverData = pd.read_csv(RiverData_file) #monthly data (from vandfoering.maanedsresultat layer in Task 1.1)

# data preparation
idxMoelleAA = RiverData['beregningspunktlokalid'].str.contains('MOELLEAA') #only Moelleaa river
RiverData = RiverData[idxMoelleAA]
RiverDataYY = RiverData[RiverData['aar']==2019].reset_index(drop=True) #only 2019
RiverDataMM = RiverDataYY[RiverDataYY['maaned'] == 'september'].reset_index(drop=True) #need to discuss
idxUp = RiverDataMM['beregningspunktlokalid'].str.contains('3687')
idxDown = RiverDataMM['beregningspunktlokalid'].str.contains('13500')
Up = RiverDataMM[idxUp].index.astype(int)[0]
Down = RiverDataMM[idxDown].index.astype(int)[0]
RiverMM_UP_DOWN = RiverDataMM.iloc[Up:Down+1].reset_index(drop=True)

# initialize variables
N = len(RiverMM_UP_DOWN) #number of river nodes, need to confirm Task 1.1, 44?
RiverQ = pd.DataFrame(np.zeros([N,6]))
RiverQ.columns = ['flow','ID','x_coord','y_coord','distance','Qadded']
RiverQ[['flow','ID','x_coord','y_coord']] = RiverMM_UP_DOWN[['vandfoering','beregningspunktlokalid','X', 'Y']]
RiverQ.index = RiverQ.pop('ID')

RiverC = pd.DataFrame(np.zeros([N,4],dtype=float))
RiverC.columns = ['conc','ID','x_coord','y_coord']
EQS_exc = RiverC.copy()
EQS_exc.columns = ['Exc','ID','x_coord','y_coord']

# calculate distance from beginning of river
RiverQ['distance'] = [float(RiverQ.index[i].split("_")[-1]) - 3687 for i in range(len(RiverQ))]

# simple model advection-dilution model
C0 = 0.67 #need to discuss, 
RiverC['conc'][0] = C0
theta = 0.0231
t_CSO = 4.3*60*60
CSO_conc = 5 # pollutant concentration 50th percentile, 95th percentile = 60
EQS = 4.9 #regulation for copper
# RiverQ['Qadded'][0] = RiverQ['flow'][0] #need to discuss

hubname = CSOdata['HubName'].to_list()

for i in range(1, N, 1):
    # check ì there is an overflow connected to this point
    isThereCSO = RiverQ.index[i] in hubname
    
    # find these CSO structures are located in the CSO data matrix
    idxCSO = CSOdata[CSOdata['HubName'] == RiverQ.index[i]].index
    
    if len(idxCSO) > 0: # if connected to CSO
        # create an empty value for MP flux from CSO
        CSO_flux = 0
        
        #create an empty value for total CSO flow added to river
        CSO_Qtot = 0 
        
        for j in range(len(idxCSO)): # loop over the connected CSOs
            # calculate the volume of the single discharge event
            V_CSO = CSOdata['Vol_m3_yea'][idxCSO[j]] * theta
            
            #calculate the flow fo the single discharge event
            Q_CSO = V_CSO/t_CSO
            CSO_flux = CSO_flux + Q_CSO * CSO_conc
            CSO_Qtot = CSO_Qtot + Q_CSO
            print(CSO_Qtot, CSO_flux)
            
        # calculate total flow added to the river by CSO
        RiverQ['Qadded'][i] = RiverQ['Qadded'][i-1]+CSO_Qtot
            
        # calculate concentration for the node
        RiverC['conc'][i] = (RiverC['conc'][i-1] * (RiverQ['flow'][i-1] + RiverQ['Qadded'][i-1]) 
                                 + CSO_flux) / (RiverQ['flow'][i] + RiverQ['Qadded'][i])
    else:
        # calculate flow added to the river (no CSO)
        RiverQ['Qadded'][i] = RiverQ['Qadded'][i-1]
        
        # calculate concentration for the node
        # RiverC['conc'][i] = (RiverC['conc'][i-1] * RiverQ['flow'][i-1] + RiverQ['Qadded'][i-1]) / (RiverQ['flow'][i] + RiverQ['Qadded'][i])
        RiverC['conc'][i] = RiverC['conc'][i-1] * (RiverQ['flow'][i-1] + RiverQ['Qadded'][i-1]) / (RiverQ['flow'][i] + RiverQ['Qadded'][i])

EQS_exc['Exc'] = RiverC['conc'] > EQS

# export and visualize results

# plot concentration profile along river
fig, ax = plt.subplots()
fig.set_size_inches(6, 3)
plt.plot(RiverQ['distance'],RiverC['conc'], color='orange')
plt.ylabel('Cu concentration [ug/l]') #add y-label
plt.xlabel('Distance (m)') #add x-label
plt.savefig("output/RiverC_5.png", bbox_inches='tight')
plt.close()

# Export RiverC & EQS_exc as csv to import to QGIS

RiverC['ID'], EQS_exc['ID'] = RiverQ.index, RiverQ.index
RiverC['x_coord'] = list(RiverQ['x_coord'])
RiverC['y_coord'] = list(RiverQ['y_coord'])
RiverC.to_csv('output/RiverC_5.csv')
RiverQ.to_csv('output/RiverQ_5.csv')

# plot concentration profile along river with exceedance
# plt.plot(RiverQ['distance'], EQS_exc['Exc'])
# plt.close()

        