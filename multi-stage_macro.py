CFViewBackward(912)
"""This macro provides CFView post-processing for multi-stage compressor.
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
dalpha = 30
WorkSplitR1 = 35

project_name = 'MSD_48bladesR2'
case_name = '4kgs_FR'
file_dir = 'C:/Users/msais/Box Sync/Thesis Work/Multi-Stage_data/DiffuserConstArea/WorkSplitRotor1=' + \
    str(WorkSplitR1) + '/Stator' + str(dalpha) + 'deg/' + \
    project_name + '/' + project_name + '_' + case_name + '/'
RunFile = str(project_name + '_' + case_name + '.run')
LossFile = "LossData_cyl.txt"

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -------------------------------------Input data------------------------------
Cp = 1006  # Specific Heat
Rgas = 287  # Universal gas constant
N = 22363  # Shaft Speed
m = 4  # mass flow rate
r6 = 0.2  # Outlet radius
Lb = [0.05247181, 0.025, 0.0775]  # Hydraulic Length
Dh = [0.05640476, 0.03375, 0.010875]  # Hydraulic Diameter
cl = [0.0005, 0.0004, 0.0003]  # Clearance for each row
Z = [24, 24, 24]    # Number of blades for each row
nu = 3.8e-5  # Kinematic viscosity

# Reference values
Vref = 100.74  # velocity (Inlet velocity)
rho = 1.2  # density
Pref = 101325  # Pressure

w = 2 * np.pi * N / 60  # Angular Velocity
U6 = r6 * w
phi = (m / rho) / (U6 * ((2 * r6)**2))

# ---------------------------Input for Post-processing ------------------------
# Loss file data
nrows = 3  # Number of blade rows
nsect = 7  # Number of planes to create

# 3D-View Contour Plot
span = [0.5, 0.9]
EntropyRange = [[-10, 120], [-10, 170]]

# Quantities
Qnt_Cart = {'ds/dx': 'Grad (Entropy)_X',
            'Cp': '(Static Pressure-101325)/(0.5*1.2*(Magnitude of V*Magnitude of V+(Vt-Wt)*(Vt-Wt)))'}

# Losses
Qnt_LossParm = ["dH_inc", "dH_bl", "dH_sf",
                "dH_cl", "dH_vld", "dH_rc", "dH_df", "dH_lk"]
f_inc = 0.5

# Averaged Quantities

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

b = np.zeros((2 * nrows, 2))
for i in range(2 * nrows):
    b = ((r[i][1] - r[i][0])**2 + (x[i][1] - x[i][0])**2)**0.5

# -------------------------Initializing Parameters-----------------------------
# Scalar data
Vt = np.zeros(nsect)  # Tangential Velocity
T0 = np.zeros(nsect)  # Total Temperature
P0 = np.zeros(nsect)  # Total Pressure
P0rel = np.zeros(nsect)  # Relative Total Pressure
sw = np.zeros(nsect)  # Swirl
r = np.zeros(nsect)  # Radius at mid-span of section plane
s = np.zeros(nsect)  # Entropy
cf = np.zeros(nrows)  # Skin Friction
W = np.zeros(nsect)
W1tip = np.zeros(nsect)
P = np.zeros(nsect)
rho = np.zeros(nsect)

FileOpenProject(file_dir + RunFile,)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ==============================================================================
# 3D View
# ==============================================================================

# -----------------------Contour Plots at Differnet Span-----------------------

# R1CamPos = [0.261069,-0.0804571,0.0386409,0.158219,-0.0247637,0.0239268,-0.474188,-0.880254,-0.0172853,0.0471529,0.0703255]
# S1CamPos = [0.214698, -0.0509236, 0.00635897, 0.16529, -0.00694214, 0.0328556, -0.26149, -0.69639, 0.66833, 0.0285029, 0.0425103]
# R2CamPos = [0.360388,0.0561327,-0.0619866,0.215378,0.0653739,0.0752378,0.285215,-0.887805,0.361185,0.0799439,0.119231]

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


for i in range(len(span)):
    SelectFromProject('CUT' + str(i + 1),)

    # Entropy
    QntFieldScalar('Entropy',)
    SclContourStrip()
    ColormapStripesOnly()
    ColormapNumOfSmallTicks(3)
    ColormapTicksNumberTextType(10, 12, 2, 0, 1, 0, 0, 1, 0, 0, 0, 0)
    ColormapLabelTextType(10, 12, 2, 2, 1, 0, 0, 1, 0, 0, 0, 0)
    RprRangeIn(EntropyRange[i][0], EntropyRange[i][1])
    SetCamera(0.261069, -0.0804571, 0.0386409, 0.158219, -0.0247637,
              0.0239268, -0.474188, -0.880254, -0.0172853, 0.0471529, 0.0703255)
    Print(8, 0, 0, 1, 100, 1317, 704, 0, file_dir + 'EntropyB2b' +
          str(span[i]) + 'R1_' + str(WorkSplitR1) + '_' + str(dalpha) + '.png', '', 1, 1, 1)
    SetCamera(0.214698, -0.0509236, 0.00635897, 0.16529, -0.00694214,
              0.0328556, -0.26149, -0.69639, 0.66833, 0.0285029, 0.0425103)
    Print(8, 0, 0, 1, 100, 1920, 1080, 0, file_dir + 'EntropyB2b' +
          str(span[i]) + 'S1_' + str(WorkSplitR1) + '_' + str(dalpha) + '.png', '', 1, 1, 1)
    SetCamera(0.360388, 0.0561327, -0.0619866, 0.215378, 0.0653739,
              0.0752378, 0.285215, -0.887805, 0.361185, 0.0799439, 0.119231)
    Print(8, 0, 0, 1, 100, 1701, 873, 0, file_dir + 'EntropyB2b' +
          str(span[i]) + 'R2_' + str(WorkSplitR1) + '_' + str(dalpha) + '.png', '', 1, 1, 1)

    # Relative Velocity
    QntFieldScalar('Magnitude of W',)
    SclContourStrip()
    ColormapStripesOnly()
    ColormapLabelTextType(10, 12, 2, 2, 1, 0, 0, 1, 0, 0, 0, 0)
    SclContourStrip()
    RprRangeIn(0, 400)
    SetCamera(0.261069, -0.0804571, 0.0386409, 0.158219, -0.0247637,
              0.0239268, -0.474188, -0.880254, -0.0172853, 0.0471529, 0.0703255)
    Print(8, 0, 0, 1, 100, 1317, 704, 0, file_dir + 'RelVelB2b' +
          str(span[i]) + 'R1_' + str(WorkSplitR1) + '_' + str(dalpha) + '.png', '', 1, 1, 1)
    SetCamera(0.214698, -0.0509236, 0.00635897, 0.16529, -0.00694214,
              0.0328556, -0.26149, -0.69639, 0.66833, 0.0285029, 0.0425103)
    Print(8, 0, 0, 1, 100, 1920, 1080, 0, file_dir + 'RelVelB2b' +
          str(span[i]) + 'S1_' + str(WorkSplitR1) + '_' + str(dalpha) + '.png', '', 1, 1, 1)
    SetCamera(0.360388, 0.0561327, -0.0619866, 0.215378, 0.0653739,
              0.0752378, 0.285215, -0.887805, 0.361185, 0.0799439, 0.119231)
    Print(8, 0, 0, 1, 100, 1701, 873, 0, file_dir + 'RelVelB2b' +
          str(span[i]) + 'R2_' + str(WorkSplitR1) + '_' + str(dalpha) + '.png', '', 1, 1, 1)
    DeleteAll()

# ------------------------------Cartesian Plots----------------------------------

ViewNum = 2
Curvenum = 1
QntFieldScalar('Entropy')
FieldGradient('Entropy')
QntFieldVector('Grad (Entropy)')
for key in Qnt_Cart.keys():
    QntFieldDerived(0, key, Qnt_Cart[key], '', '0')

for key in Qnt_Cart.keys():
    ViewActivate(RunFile + ':1')
    SelectFromProject('row 1_blade_(r.p.m. -22363)',
                      'row 2_blade', 'row 3_blade_(r.p.m. -22363)')
    QntFieldScalar(key)
    SclPlotNormalizedGridLine(0, 0.5, 0, 1, 'row_1_flux_1_Main_Blade_skin.Imin blade_(aap-ps)_rotating',
                              0, 0.5, 0, 1, 'row_1_flux_1_Main_Blade_skin.Imin blade_(aap-ss)_rotating', 0)
    ViewActivate(RunFile + ':' + str(ViewNum))
    PlotPlaneX()
    SelectPlotCurves('Gridline I=0.5 on row_1_flux_1_Main_Blade_skin.Imin blade_(aap-ps)_rotating',
                     'Gridline I=0.5 on row_1_flux_1_Main_Blade_skin.Imin blade_(aap-ss)_rotating')
    PlotCurvesMerge('Gridline I=0.5 on row_1_flux_1_Main_Blade_skin.Imin blade_(aap-ps)_rotating',
                    'Gridline I=0.5 on row_1_flux_1_Main_Blade_skin.Imin blade_(aap-ss)_rotating')
    ActivePlotCurveOutput(file_dir + '3DB2b0.5R1_' + str(WorkSplitR1) +
                          '_' + str(dalpha) + '.dat', 'merged curves ' + str(Curvenum))
    Curvenum += 1
    SelectPlotCurves('Gridline I=0.5 on row_2_flux_1_Main_Blade_skin.Imin blade_(aap-ps)',
                     'Gridline I=0.5 on row_2_flux_1_Main_Blade_skin.Imin blade_(aap-ss)')
    PlotCurvesMerge('Gridline I=0.5 on row_2_flux_1_Main_Blade_skin.Imin blade_(aap-ps)',
                    'Gridline I=0.5 on row_2_flux_1_Main_Blade_skin.Imin blade_(aap-ss)')
    ActivePlotCurveOutput(file_dir + '3DB2b0.5S1_' + str(WorkSplitR1) +
                          '_' + str(dalpha) + '.dat', 'merged curves ' + str(Curvenum))
    Curvenum += 1
    SelectPlotCurves('Gridline I=0.5 on row_3_flux_1_Main_Blade_skin.Imin blade_(aap-ps)_rotating',
                     'Gridline I=0.5 on row_3_flux_1_Main_Blade_skin.Imin blade_(aap-ss)_rotating')
    PlotCurvesMerge('Gridline I=0.5 on row_3_flux_1_Main_Blade_skin.Imin blade_(aap-ps)_rotating',
                    'Gridline I=0.5 on row_3_flux_1_Main_Blade_skin.Imin blade_(aap-ss)_rotating')
    ActivePlotCurveOutput(file_dir + '3DB2b0.5R2_' + str(WorkSplitR1) +
                          '_' + str(dalpha) + '.dat', 'merged curves ' + str(Curvenum))
    Curvenum += 1
    DeletePlot()
    ViewNum += 1

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ==============================================================================
# Cylindrical view - Calculating Loss and Exporting to File(Surface Area Average Values)
# ==============================================================================
ViewActivate(RunFile + ':1')
sys.stdout = open(file_dir + LossFile, "w")

#---------------------Evaluating Parameters for Loss---------------------------
ViewOpenRTZ(-0.783962, -0.683962, 0.508295, 0.608295)
QntFieldDerived(0, 'swirl', 'sqrt((x*x)+(y*y))*Vt', '', '0')
QntFieldDerived(0, 'r', 'sqrt((x*x)+(y*y))', '', '0')
LimitsFull()
ViewOriginal(1,)
CutPlaneSave(0.11, 0, 0, 0, 0, 1, 1)  # Rotor 1 Inlet
CutPlaneSave(0.0760, 0, 0.05247182, 0.27081969, 0, 1, 1)  # Rotor 1 Outlet
CutPlaneSave(0.0765, 0, 0.05547558, 0.27442423, 0, 1, 1)  # Stator 1 Inlet
CutPlaneSave(0.1185, 0, 0.09072516, 0.65675436, 0, 1, 1)  # Stator 1 Outlet
CutPlaneSave(0.1215, 0, 0.09387836, 0.71827631, 0, 1, 1)  # Rotor 2 Inlet
CutPlaneSave(0.2001, 0, 0.12, 1, 0, 0, 1)  # Rotor 2 Outlet
CutPlaneSave(0.3, 0, 0.12, 1, 0, 0, 1)  # Diffuser Outlet
for c in range(nsect):
    plne_nme = 'CUT' + str(c + 3)
    SelectFromProject(plne_nme)
    QntFieldScalar('Absolute Total Temperature')
    T0[c] = WeightedIntegral()
    QntFieldScalar('Vt')
    Vt[c] = WeightedIntegral()
    QntFieldScalar('Relative Total Pressure')
    P0rel[c] = WeightedIntegral()
    QntFieldScalar('Absolute Total Pressure')
    P0[c] = WeightedIntegral()
    QntFieldScalar('swirl')
    sw[c] = WeightedIntegral()
    QntFieldScalar('r')
    r[c] = WeightedIntegral()
    QntFieldScalar('Entropy')
    s[c] = WeightedIntegral()
    QntFieldScalar('Magnitude of W')
    W[c] = WeightedIntegral()
    QntFieldScalar('Static Pressure')
    P[c] = WeightedIntegral()
    QntFieldScalar('Density')
    rho[c] = WeightedIntegral()

SelectFromProject('CUT9.D9')
QntFieldScalar('Absolute Total Temperature')
T0[-1] = WeightedIntegral()
QntFieldScalar('Vt')
Vt[-1] = WeightedIntegral()
QntFieldScalar('Relative Total Pressure')
P0rel[-1] = WeightedIntegral()
QntFieldScalar('Absolute Total Pressure')
P0[-1] = WeightedIntegral()
QntFieldScalar('swirl')
sw[-1] = WeightedIntegral()
QntFieldScalar('r')
r[-1] = WeightedIntegral()
QntFieldScalar('Entropy')
s[-1] = WeightedIntegral()
QntFieldScalar('Magnitude of W')
W[-1] = WeightedIntegral()
QntFieldScalar('Static Pressure')
P[-1] = WeightedIntegral()
QntFieldScalar('Density')
rho[-1] = WeightedIntegral()

print('\nJ    Swirl      T0          P0      P0rel       r     Entropy')

for c in range(nsect):
    print(str(c + 1) + '    ' + '%.4f' % (sw[c]) + '    ' + '%.4f' % (T0[c]) + '    ' + '%.4f' % (
        P0[c]) + '    ' + '%.4f' % (P0rel[c]) + '    ' + '%.4f' % (r[c]) + '    ' + '%.4f' % (s[c]))

# Evaluating blade properties averaged
QntSolidScalar('Cf')
for i in range(nrows):
    if i % 2 == 0:  # Check for rotor row number
        SelectFromProject('row ' + str(i + 1) + '_blade_(r.p.m. -22363)')
    else:  # Else a stator
        SelectFromProject('row ' + str(i + 1) + '_blade')
    cf[i] = SclAverage()

# ---------------------Calculating Individual Losses---------------------------
loss = np.zeros((len(Qnt_LossParm), nrows))
Df = np.zeros(nrows)
dH_Aero = np.zeros(nrows)
c = 0
for i in range(nrows):
    # Aerodynamic Enthalpy
    dH_Euler = Cp * (T0[c + 1] - T0[c])
#     Df = 1 - (W[c + 1] / W1tip[c + 1]) + ((0.75 * dH_Euler) / U6 **
#               2) / ((W1tip[c + 1] / W[c + 1]) * [(Z[c] / np.pi) *
#               (1 - (r[c][1] / np.average(r[c + 1]))) + 2 * (r[c][1] / np.average(r[c + 1]))])
    c += 2
print(dH_Euler)

# c = 0
# for i in range(nrows):
#
#     # Skin Friction loss
#     loss[2][i] = 2 * cf[i] * Lb[i] * Wavg[i]**2 / Dh[i]
#
#     # Clearance Loss
#     loss[3][i] = 0.6 * (cl / 0.0675) * (Vt[c + 1]) * ())((4 * np.pi) / (0.0675 * 24)) * (0.1125**2)
    # c += 2
# print(loss[0])

# Vaneless diffuser loss
dH_vld = Cp * T0[-2] * ((P[-1] / P0[-1])**(7 / 2) - (P[-1] / P0[-2])**(7 / 2))
print(dH_vld)

# Disk Friction only for rotor 2
Re_df = U6 * r6 / nu
f_df = 0.0622 / Re_df**0.2
rho_avg = (rho[4] + rho[5]) / 2
dH_df = f_df * (rho_avg * r6**2 * U6**3) / (4 * m)
print(dH_df)

#Leakage Loss for Rotor 2
b_avg = (b[-2] + b[-1]) / 2
r_avg = 0.5 * (r[5][0] + 0.5 * (r[4][0] + r[4][1]))
dP_cl = m * (r[5][0] * Vt[5] - r_avg * Vtm[4]) / (Z[2] * r_avg * b_avg)
U_cl = 0.816 * (2 * dP_cl * rho[5])**0.5
m_cl = rho[5] * Z[2] * cl[2] * L_theta * U_cl
dH_lk = m_cl * U_cl * U6 / (2 * m)
print(dH_lk)

# Rowheaders = ['Row','Cf','SkinFricLoss','IncLoss','MixLoss','IntLoss','ExtLoss','Overall(Int+Ext)']
# mainlist = [[] for i in range(len(Rowheaders))]
# for i  in range(nrows+1):
#     if i==0:
#         for j in range(len(Rowheaders)):
#             mainlist[i].append(Rowheaders[j])
#     else:
#         mainlist[i].append(i)
#         mainlist[i].append('%.4f'%float(cf[i-1]))
#         mainlist[i].append('%.4f'%float(dssf[i-1]))
#         mainlist[i].append('%.4f'%float(dsinc[i-1]))
#         mainlist[i].append('%.4f'%float(dsmix[i-1]))
#         mainlist[i].append('%.4f'%float(dsint[i-1]))
#         mainlist[i].append('%.4f'%float(dsext[i-1]))
#         mainlist[i].append('%.4f'%float(ds[i-1]))
# print('\n'+tabulate(mainlist))

sys.stdout.close()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # ==============================================================================
# # Open Turbomachinery Mode
# # ==============================================================================
# SetTurboMode()
# OpenTurboModeStandard3DView()
# OpenTurboModeBladeToBladeView()
# OpenTurboModeBladeView()
# OpenTurboModeMeridionalView(0, 0)
#
# # ==============================================================================
# # Blade to blade plots
# # ==============================================================================
# ViewActivate(RunFile + ':' + str(ViewNum))
# GmtToggleBoundary()
# QntFieldScalar('Absolute Total Pressure')
# SclContourStrip()
# ViewNum += 1
# # ==============================================================================
# # Blade view
# # ==============================================================================
# ViewActivate(RunFile + ':' + str(ViewNum))
# GmtToggleBoundary()
# QntFieldScalar('Absolute Total Pressure')
# SclContourStrip()
# ViewNum += 1
# # ==============================================================================
# # Azimuthal Plots
# # ==============================================================================
# ViewActivate(RunFile + ':1')
# # ViewActivate(':'+str(ViewNum))
# ViewActivate(RunFile + ' meridional view:' + str(ViewNum))
# GmtToggleBoundary()
# QntFieldScalar('Static Pressure')
# SclContourStrip()
# ViewNum += 1
