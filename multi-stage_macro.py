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
WorkSplitR1 = 35
dalpha = 25

project_name = 'MSD_5sect'
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

chord_r1 = 0.0534744
chord_s1 = 0.07272

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
WRange = [[0, 400], [0, 400]]

# Losses
Qnt_LossParm = ["Inc_Loss", "BldeLoading_Loss", "SkinFriction_Loss",
                "Clearance_Loss", "Recir_Loss", "DiskFriction_Loss", "Leakage_Loss", "VanelessDiff_Loss"]
f_inc = 0.5
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

bw = np.zeros(nsect)
for i in range(nsect):
    bw[i] = ((r[i][1] - r[i][0])**2 + (x[i][1] - x[i][0])**2)**0.5


# -------------------------Initializing Parameters-----------------------------
# Scalar data
Vt = np.zeros(nsect)  # Tangential Velocity
T0 = np.zeros(nsect)  # Total Temperature
P0 = np.zeros(nsect)  # Total Pressure
P0rel = np.zeros(nsect)  # Relative Total Pressure
sw = np.zeros(nsect)  # Swirl
rm = np.zeros(nsect)  # Radius at mid-span of section plane
s = np.zeros(nsect)  # Entropy
W = np.zeros(nsect)
W_tip = np.zeros(nsect)
W_hub = np.zeros(nsect)
V_tip = np.zeros(nsect)
P = np.zeros(nsect)
rho = np.zeros(nsect)
V = np.zeros(nsect)
Vm_m = np.zeros(nsect)
Wt = np.zeros(nsect)
alpha = np.zeros(nsect)

cf = np.zeros(nrows)  # Skin Friction

W_hub = [114.073, 1, 287.376, 1, 243.4279, 1, 1]
W_tip = [296.548, 1, 230.909, 1, 337.775, 1, 1]
V_tip = [5.449, 1, 213.587, 1, 122.495, 1, 1]
Vm_m = [124.380, 1, 130.725, 1, 135.77, 1, 1]

FileOpenProject(file_dir + RunFile,)
ViewNum = 2

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ==============================================================================
# 3D View
# ==============================================================================

# -----------------------Contour Plots at Differnet Span-----------------------

# R1CamPos = [0.261069,-0.0804571,0.0386409,0.158219,-0.0247637,0.0239268,-0.474188,-0.880254,-0.0172853,0.0471529,0.0703255]
# S1CamPos = [0.214698, -0.0509236, 0.00635897, 0.16529, -0.00694214, 0.0328556, -0.26149, -0.69639, 0.66833, 0.0285029, 0.0425103]
# R2CamPos = [0.360388,0.0561327,-0.0619866,0.215378,0.0653739,0.0752378,0.285215,-0.887805,0.361185,0.0799439,0.119231]
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
    RprRangeIn(WRange[i][0], WRange[i][1])
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
    rm[c] = WeightedIntegral()
    QntFieldScalar('Entropy')
    s[c] = WeightedIntegral()
    QntFieldScalar('Magnitude of W')
    W[c] = WeightedIntegral()
    QntFieldScalar('Static Pressure')
    P[c] = WeightedIntegral()
    QntFieldScalar('Density')
    rho[c] = WeightedIntegral()
    QntFieldScalar('Magnitude of V')
    V[c] = WeightedIntegral()
    QntFieldScalar('Wt')
    Wt[c] = WeightedIntegral()
    QntFieldScalar('atan(Vt/Vm)')
    alpha[c] = WeightedIntegral()

SelectFromProject('CUT9.D9')
QntFieldScalar('Absolute Total Temperature')
T0[-1] = WeightedIntegral()
QntFieldScalar('Relative Total Pressure')
P0rel[-1] = WeightedIntegral()
QntFieldScalar('Absolute Total Pressure')
P0[-1] = WeightedIntegral()
QntFieldScalar('swirl')
sw[-1] = WeightedIntegral()
QntFieldScalar('r')
rm[-1] = WeightedIntegral()
QntFieldScalar('Entropy')
s[-1] = WeightedIntegral()
QntFieldScalar('Static Pressure')
P[-1] = WeightedIntegral()
QntFieldScalar('atan(Vt/Vm)')
alpha[c] = WeightedIntegral()

# print('\nJ    Swirl      T0          P0      P0rel       r     Entropy')
Rowheaders = ['J', 'Swirl', 'T0', 'P0', 'P0rel', 'alpha', 'Entropy']
mainlist = [[] for i in range(nsect + 1)]
for i in range(nsect + 1):
    if i == 0:
        for j in range(len(Rowheaders)):
            mainlist[i].append(Rowheaders[j])
    else:
        mainlist[i].append(i)
        mainlist[i].append('%.4f' % float(sw[i - 1]))
        mainlist[i].append('%.4f' % float(T0[i - 1]))
        mainlist[i].append('%.4f' % float(P0[i - 1]))
        mainlist[i].append('%.4f' % float(P0rel[i - 1]))
        mainlist[i].append('%.4f' % float(alpha[i - 1]))
        mainlist[i].append('%.4f' % float(s[i - 1]))
print('\n' + tabulate(mainlist) + '\n')

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
dH_Euler = np.zeros(nrows)
dH_inc = np.zeros(nrows)
dH_bld = np.zeros(nrows)
dH_sf = np.zeros(nrows)
dH_cl = np.zeros(nrows)
dH_rc = np.zeros(nrows)
dH_df = np.zeros(nrows)
dH_lk = np.zeros(nrows)

U = w * rm
f__inc = 0.5

c = 0
for i in range(nrows):
    # Aerodynamic Enthalpy
    dH_Euler[i] = Cp * (T0[c + 1] - T0[c])

    # Diffusion Factor for rotors
    if i == 0:
        ss_r1 = 2 * np.pi * rm[c] / Z[i]
        sol_r1 = chord_r1 / ss_r1
        DR = W[c + 1] / W[c]
        Df[i] = 1 - DR + (DR - 1) / (2 * sol_r1)
    elif i == 2:
        Wratio = W[c + 1] / W_tip[c]
        sol = Z[i] / np.pi
        rad_ratio = r[c][1] / (0.5 * (r[c + 1][0] + r[c + 1][1]))
        num = 0.75 * dH_Euler[i] / U[c + 1]**2
        Df[i] = 1 - Wratio + (Wratio * num) / \
            (sol * (1 - rad_ratio) + 2 * rad_ratio)
    else:
        ss_s1 = 2 * np.pi * rm[c] / Z[i]
        sol = chord_s1 / ss_s1
        DR = V[c + 1] / V[c]
        Df[i] = 1 - DR + (DR - 1) / (2 * sol)

    # Incidence Loss
    dH_inc[i] = f_inc * Wt[c]**2 / 2

    """Losses calculated assuming all loss is similar to rotor 2 loss s=correlations"""
    # Blade Loading Loss
    dH_bld[i] = 0.05 * Df[i] * U[c + 1]**2

    # Skin Friction Loss
    W_avg = (V_tip[c] + V[c + 1] + W_tip[c] + 2 * W_hub[c] + 3 * W[c + 1]) / 8
    dH_sf[i] = 2 * cf[i] * Lb[i] * W_avg**2 / Dh[i]

    # Clearance Loss for rotor 2
    if i == 2:
        sold = (4 * np.pi) / (bw[c + 1] * Z[i])
        frac = ((r[c][1]**2 - r[c][0]**2)) / \
            ((r[c + 1][0] - r[c][1]) * (1 + rho[c + 1] / rho[c]))
        dH_cl[i] = 0.6 * (cl[i] / bw[c + 1]) * abs(Vt[c + 1]) * \
            (sold * frac * abs(Vt[c + 1]) * Vm_m[c])**0.5

    # Disk Friction
    Re_df = U[c + 1] * rm[c + 1] / nu
    f_df = 0.0622 / Re_df**0.2
    rho_avg = (rho[c] + rho[c + 1]) / 2
    dH_df[i] = f_df * (rho_avg * rm[c + 1]**2 * U[c + 1]**3) / (4 * m)

    # Recirculation Loss
    dH_rc[i] = abs(8e-5 * math.sinh(3.5 * alpha[c + 1]**3)
                   * Df[i]**2 * U[c + 1]**2)

    # Leakage Loss for Rotors only
    if i % 2 == 0:
        b_avg = (bw[c] + bw[c + 1]) / 2
        r_avg = 0.5 * (r[c + 1][0] + 0.5 * (r[c][0] + r[c][1]))
        r1_m = 0.5 * (r[c][0] + r[c][1])
        dP_cl = m * (r[c + 1][0] * abs(Vt[c + 1]) - r_avg *
                     abs(Vt[c])) / (Z[i] * r_avg * b_avg)
        U_cl = 0.816 * (2 * dP_cl * rho[c + 1])**0.5
        m_cl = rho[c + 1] * Z[i] * cl[i] * Lb[i] * U_cl
        dH_lk[i] = m_cl * U_cl * U6 / (2 * m)

    c += 2

# Vaneless diffuser loss
dH_vld = Cp * T0[5] * ((P[6] / P0[6])**(1 / 3.5) -
                       (P[6] / P0[5])**(1 / 3.5))

dH_int = dH_inc + dH_bld + dH_sf + dH_cl  # Internal Loss
dH_par = dH_rc + dH_df + dH_lk  # Exernal Loss
dH_act = dH_Euler + dH_par  # Actual Work Done
Eff_isen = (dH_Euler - dH_int) / dH_act
Eff_isen[1] = 0
Rowheaders = ['Row', 'Df', 'Euler_Work'] + Qnt_LossParm[0:7] + ['Efficiency']
mainlist = [[] for i in range(nrows + 1)]
for i in range(nrows + 1):
    if i == 0:
        for j in range(len(Rowheaders)):
            mainlist[i].append(Rowheaders[j])
    else:
        mainlist[i].append(i)
        mainlist[i].append('%.4f' % float(Df[i - 1]))
        mainlist[i].append('%.4f' % float(dH_Euler[i - 1]))
        mainlist[i].append('%.4f' % float(dH_inc[i - 1]))
        mainlist[i].append('%.4f' % float(dH_bld[i - 1]))
        mainlist[i].append('%.4f' % float(dH_sf[i - 1]))
        mainlist[i].append('%.4f' % float(dH_cl[i - 1]))
        mainlist[i].append('%.4f' % float(dH_rc[i - 1]))
        mainlist[i].append('%.4f' % float(dH_df[i - 1]))
        mainlist[i].append('%.4f' % float(dH_lk[i - 1]))
        mainlist[i].append('%.4f' % float(Eff_isen[i - 1]))
print('\n' + tabulate(mainlist))

dH_inc_ov = np.sum(dH_inc)
dH_bld_ov = np.sum(dH_bld)
dH_sf_ov = np.sum(dH_sf)
dH_cl_ov = np.sum(dH_cl)
dH_rc_ov = np.sum(dH_rc)
dH_df_ov = np.sum(dH_df)
dH_lk_ov = np.sum(dH_lk)
dH_Euler_ov = np.sum(dH_Euler)

dH_int = dH_inc_ov + dH_bld_ov + dH_sf_ov + dH_cl_ov + dH_vld  # Internal Loss
dH_par = dH_rc_ov + dH_df_ov + dH_lk_ov  # Exernal Loss
dH_act = dH_Euler_ov + dH_par  # Actual Work Done
Eff_isen = (dH_Euler_ov - dH_int) / dH_act

Rowheaders = ['Euler_Work'] + Qnt_LossParm + ['Efficiency']
mainlist = [[] for i in range(2)]
for i in range(2):
    if i == 0:
        for j in range(len(Rowheaders)):
            mainlist[i].append(Rowheaders[j])
    else:
        mainlist[i].append('%.4f' % float(dH_Euler_ov))
        mainlist[i].append('%.4f' % float(dH_inc_ov))
        mainlist[i].append('%.4f' % float(dH_bld_ov))
        mainlist[i].append('%.4f' % float(dH_sf_ov))
        mainlist[i].append('%.4f' % float(dH_cl_ov))
        mainlist[i].append('%.4f' % float(dH_rc_ov))
        mainlist[i].append('%.4f' % float(dH_df_ov))
        mainlist[i].append('%.4f' % float(dH_lk_ov))
        mainlist[i].append('%.4f' % float(dH_vld))
        mainlist[i].append('%.4f' % float(Eff_isen))
print('\n' + tabulate(mainlist))

sys.stdout.close()

# # Blade Loading Loss for rotor 2
# dH_bld[2] = 0.05 * Df[2] * U6**2
#
# # Skin Friction Loss for rotor 2
# W_avg = (V_tip[4] + V[5] + W_tip[4] + 2 * W_hub[4] + 3 * W[5]) / 8
# dH_sf[2] = 2 * cf[2] * Lb[2] * W_avg**2 / Dh[2]
#
#
# # Clearance Loss for rotor 2
# sold = (4 * np.pi) / (bw[5] * Z[2])
# frac = ((r[4][1]**2 - r[4][0]**2)) / \
#     ((r[5][0] - r[4][1]) * (1 + rho[5] / rho[4]))
# dH_cl[2] = 0.6 * (cl[2] / bw[5]) * abs(Vt[5]) * \
#     (sold * frac * abs(Vt[5]) * Vm_m[4])**0.5
#
# # Vaneless diffuser loss
# dH_vld[2] = Cp * T0[5] * ((P[6] / P0[6])**(1 / 3.5) -
#                           (P[6] / P0[5])**(1 / 3.5))
#
#
# # Disk Friction only for rotor 2
# Re_df = U6 * r6 / nu
# f_df = 0.0622 / Re_df**0.2
# rho_avg = (rho[4] + rho[5]) / 2
# dH_df[2] = f_df * (rho_avg * r6**2 * U6**3) / (4 * m)
#
# # Recirculation Loss for Rotor 2
# dH_rc[2] = abs(8e-5 * math.sinh(3.5 * alpha[5]**3) * Df[2]**2 * U6**2)
#
# # Leakage Loss for Rotor 2
# b_avg = (bw[4] + bw[5]) / 2
# r_avg = 0.5 * (r[5][0] + 0.5 * (r[4][0] + r[4][1]))
# r1_m = 0.5 * (r[4][0] + r[4][1])
# dP_cl = m * (r[5][0] * abs(Vt[5]) - r_avg *
#              abs(Vt[4])) / (Z[2] * r_avg * b_avg)
# U_cl = 0.816 * (2 * dP_cl * rho[5])**0.5
# m_cl = rho[5] * Z[2] * cl[2] * Lb[2] * U_cl
# dH_lk[2] = m_cl * U_cl * U6 / (2 * m)
