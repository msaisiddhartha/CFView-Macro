CFViewBackward(912)
"""This macro provides CFView post-processing for multi-stage compressor.
It computes:
    1. Contour Plots at different spanwise sections
    2. Cartesian plots of X-Gradient of Entropy vs X-direction
    3. Cp distribution at 0.5 span
    4. Individual Losses through correlations and the quantities are
    mass averaged at the surface planes.
    5. Domain numbers are 3, 14 and 7 respectively or  """

import sys
import math
import numpy as np
import pylab as py
from tabulate import tabulate

# Case Data
# Case Data
WorkSplitR1 = 35
dalpha = 25

project_name = 'MSD_5sect'
case_name = '4kgs_FR'
file_dir = 'C:/Users/msais/Box Sync/Thesis Work/Multi-Stage_data/DiffuserConstArea/WorkSplitRotor1=' + \
    str(WorkSplitR1) + '/Stator' + str(dalpha) + 'deg/' + \
    project_name + '/' + project_name + '_' + case_name + '/'
RunFile = str(project_name + '_' + case_name + '.me.cfv')

# OutFile = "Scalar Data_Azimuthal averaged.txt"

Quant = ['Magnitude of V', 'Magnitude of W', 'Vm', 'Vt', 'Wt']

nsect =7
# ------------------------------Geometry Data----------------------------------
r = np.zeros((nsect, 2))
x = np.zeros((nsect, 2))
x[0] = 0
r[0] = [0.0449341, 0.11274362]

[X01, R01] = [-0.065005893244729426, 0.21920467299175289]
R = 0.186                         # Radius of the circle
[X04, R04] = [0, 0.209]  # (x,r)-coordinates for center of ellipse
[a, b] = [0.105761, R04 - r[0][1]]

# Radius at LE and TE of blade

r[1] = [0.075, 0.12]
gap_rs = np.array([0.0025, 0.00125])
r[2] = r[1] + gap_rs
s1_len = np.array([0.04, 0.02])
r[3] = r[2] + s1_len
gap_sr = np.array([0.005, 0.003])
r[4] = r[3] + gap_sr
r[5] = 0.2
r[6] = 0.3

for i in range(nsect - 2):
    x[i + 1][0] = X01 + (R**2 - (r[i + 1][0] - R01)**2)**0.5
    x[i + 1][1] = X04 + (a / b) * (b**2 - (r[i + 1][1] - R04)**2)**0.5

FileOpenProject(file_dir + RunFile,)

dm = [3, 14, 7]
Viewnum = 2
c=0
for j in Quant:
    QntFieldScalar(j)
    for i in range(nsect):
        if i%2==0:
            ViewActivate(RunFile + ':1')
            RprSection(x[i][0],r[i][0],0,x[i][1],r[i][1],0,0,0,1 ,'Section '+str(i+1),0 ,'',0)

# ViewActivate(RunFile + ':6')
#
#
# ActivePlotCurveOutput(file_dir + j + '_1_Az_avg' + '.dat' ,'Section 1 on domain3')
# ActivePlotCurveOutput(file_dir + j + '_3_Az_avg' + '.dat' ,'Section 3 on domain14')
# ActivePlotCurveOutput(file_dir + j + '_5_Az_avg' + '.dat' ,'Section 5 on domain7')
