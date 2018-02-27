CFViewBackward(912)
"""This macro provides CFView post-processing for 1.5-stage centrifugal compressor.
It computes:
    1. Contour Plots at different spanwise sections
    2. Cartesian plots of X-Gradient of Entropy vs X-direction
    3. Cp distribution at 0.5 span
    4. Individual Losses through correlations and the quantities are
    mass averaged at the surface planes."""

import sys
import math
import numpy as np
import pylab as py
from tabulate import tabulate

# Case Data
WorkSplitR1 = 15
dalpha = 25

project_name = 'MSD_Bckswp45'
case_name = '4kgs_111_mav'
file_dir = 'C:/Users/msais/Box Sync/Thesis Work/Multi-Stage_data/DiffuserConstArea/WorkSplitRotor1=' + \
    str(WorkSplitR1) + '/Stator' + str(dalpha) + 'deg/' + \
    project_name + '/' + project_name + '_' + case_name + '/'
RunFile = str(project_name + '_' + case_name + '.run')

# 3D-View Contour Plot
Qnt = ['Entropy', 'Magnitude of W', 'Relative Mach Number']
span = [0.25, 0.5, 0.95]
EntropyRange = [[-10, 120], [-10, 120], [-10, 170]]
WRange = [[0, 400], [0, 400], [0, 400]]
RMach = [[0, 1], [0, 1], [0, 1]]

# R1CamPos = [0.261069,-0.0804571,0.0386409,0.158219,-0.0247637,0.0239268,-0.474188,-0.880254,-0.0172853,0.0471529,0.0703255]
# S1CamPos = [0.214698, -0.0509236, 0.00635897, 0.16529, -0.00694214, 0.0328556, -0.26149, -0.69639, 0.66833, 0.0285029, 0.0425103]
# R2CamPos = [0.360388,0.0561327,-0.0619866,0.215378,0.0653739,0.0752378,0.285215,-0.887805,0.361185,0.0799439,0.119231]
FileOpenProject(file_dir + RunFile,)
ViewActivate(RunFile + ':1')
SetNumecaLogo(0, 0)
for i in range(len(span)):
    CutPlaneSave(span[i], 0, 0, 1, 0, 0, 2)
GmtRepetitionToggle()
GmtRepetitionNumber(3, 3, 3)
scm = [0.222479, -0.0683883, -0.00803093, 0.158146, -0.0111208,
       0.0264699, -0.26149, -0.69639, 0.66833, 0.0371132, 0.0553519]
SetCamera(scm[0], scm[1], scm[2], scm[3], scm[4], scm[5],
          scm[6], scm[7], scm[8], scm[9], scm[10],)
GmtToggleBoundary()

for j in range(len(Qnt)):
    for i in range(len(span)):
        SelectFromProject('CUT' + str(i + 1),)

        # Entropy
        QntFieldScalar(Qnt[j],)
        SclContourStrip()
        ColormapStripesOnly()
        ColormapNumOfSmallTicks(3)
        ColormapTicksNumberTextType(10, 12, 2, 0, 1, 0, 0, 1, 0, 0, 0, 0)
        ColormapLabelTextType(10, 12, 2, 2, 1, 0, 0, 1, 0, 0, 0, 0)
        if j == 0:
            RprRangeIn(EntropyRange[i][0], EntropyRange[i][1])
        if j == 1:
            RprRangeIn(WRange[i][0], WRange[i][1])
        if j == 2:
            RprRangeIn(RMach[i][0], RMach[i][1])
        SetCamera(0.261069, -0.0804571, 0.0386409, 0.158219, -0.0247637,
                  0.0239268, -0.474188, -0.880254, -0.0172853, 0.0471529, 0.0703255)
        Print(8, 0, 0, 1, 100, 1317, 704, 0, file_dir + Qnt[j] + '_B2b_' +
              str(span[i]) + '_R1_' + str(WorkSplitR1) + '_' + str(dalpha) + '.png', '', 1, 1, 1)
        SetCamera(0.214698, -0.0509236, 0.00635897, 0.16529, -0.00694214,
                  0.0328556, -0.26149, -0.69639, 0.66833, 0.0285029, 0.0425103)
        Print(8, 0, 0, 1, 100, 1920, 1080, 0, file_dir + Qnt[j] + '_B2b_' +
              str(span[i]) + '_S1_' + str(WorkSplitR1) + '_' + str(dalpha) + '.png', '', 1, 1, 1)
        SetCamera(0.360388, 0.0561327, -0.0619866, 0.215378, 0.0653739,
                  0.0752378, 0.285215, -0.887805, 0.361185, 0.0799439, 0.119231)
        Print(8, 0, 0, 1, 100, 1701, 873, 0, file_dir + Qnt[j] + '_B2b_' +
              str(span[i]) + '_R2_' + str(WorkSplitR1) + '_' + str(dalpha) + '.png', '', 1, 1, 1)
        DeleteAll()
