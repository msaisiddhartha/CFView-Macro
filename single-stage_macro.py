CFViewBackward(1210)
"""This macro provides CFView post-processing for single-stage centrifugal compressor.
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

project_name = 'baseline4'
case_name = '4kgs_SA_no_boost'
file_dir = 'C:/Users/msais/Box Sync/Thesis Work/Baseline/' + \
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
Lb = [0.15497081]
cl = [0.0005]  # Clearance for each row
Z = [24]    # Number of blades for each row
nu = 3.8e-5  # Kinematic viscosity
Wt_ac = [-184.6286]
betab = [-30]
betabi = [-62]

# Reference values
Vref = 100.74  # velocity (Inlet velocity)
rho = 1.2  # density
Pref = 101325  # Pressure

w = 2 * np.pi * N / 60  # Angular Velocity
U6 = r6 * w

#=======R1 = 0.35, Dalpha = 35====================
W_hub = [114.073]
W_tip = [271.222]
V_tip = [5.449]
Vm_m = [116.003]
# chord = 0.0534744
#=================================================

# ---------------------------Input for Post-processing ------------------------
# Loss file data
nrows = 1  # Number of blade rows
nsect = 7  # Number of planes to create

# 3D-View Contour Plot
Qnt = ['Entropy', 'Magnitude of W', 'Relative Mach Number']
span = [0.25, 0.5, 0.75, 0.9, 0.95]
EntropyRange = [[-10, 120],[-10, 120],[-10, 120],[-10, 120], [-10, 170]]
WRange = [[0, 400],[0, 400],[0, 400],[0, 400], [0, 400]]
RMach = [[0, 1],[0, 1],[0, 1],[0, 1], [0, 1], [0, 1]]

# Losses
Qnt_LossParm = ["Incidence_Loss", "BldeLoading_Loss", "SkinFriction_Loss",
                "Clearance_Loss", "Recir_Loss", "DiskFriction_Loss", "Leakage_Loss", "VanelessDiff_Loss"]

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

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ==============================================================================
# 3D View
# ==============================================================================
FileOpenProject(file_dir + RunFile,)
ViewNum = 2

# -----------------------Contour Plots at Differnet Span-----------------------

# R1CamPos = [0.261069,-0.0804571,0.0386409,0.158219,-0.0247637,0.0239268,-0.474188,-0.880254,-0.0172853,0.0471529,0.0703255]
# S1CamPos = [0.214698, -0.0509236, 0.00635897, 0.16529, -0.00694214, 0.0328556, -0.26149, -0.69639, 0.66833, 0.0285029, 0.0425103]
# R2CamPos = [0.360388,0.0561327,-0.0619866,0.215378,0.0653739,0.0752378,0.285215,-0.887805,0.361185,0.0799439,0.119231]
ViewActivate(RunFile + ':1')
SetNumecaLogo(0, 0)
for i in range(len(span)):
    CutPlaneSave(span[i], 0, 0, 1, 0, 0, 2)
GmtRepetitionToggle()
GmtRepetitionNumber(3)
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
              str(span[i]) + '_R1_' + '.png', '', 1, 1, 1)
        SetCamera(0.214698, -0.0509236, 0.00635897, 0.16529, -0.00694214,
                  0.0328556, -0.26149, -0.69639, 0.66833, 0.0285029, 0.0425103)
        Print(8, 0, 0, 1, 100, 1920, 1080, 0, file_dir + Qnt[j] + '_B2b_' +
              str(span[i]) + '_S1_' + '.png', '', 1, 1, 1)
        SetCamera(0.360388, 0.0561327, -0.0619866, 0.215378, 0.0653739,
                  0.0752378, 0.285215, -0.887805, 0.361185, 0.0799439, 0.119231)
        Print(8, 0, 0, 1, 100, 1701, 873, 0, file_dir + Qnt[j] + '_B2b_' +
              str(span[i]) + '_R2_' '.png', '', 1, 1, 1)
        DeleteAll()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ==============================================================================
# Cylindrical view - Exporting to File(Surface Area Average Values)
# ==============================================================================
ViewActivate(RunFile + ':1')

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
P = np.zeros(nsect)
rho = np.zeros(nsect)
V = np.zeros(nsect)
Wt = np.zeros(nsect)
alpha = np.zeros(nsect)
beta = np.zeros(nsect)
Wui = np.zeros(nsect)
Vm = np.zeros(nsect)
M = np.zeros(nsect)
Mrel = np.zeros(nsect)
T = np.zeros(nsect)
cf = np.zeros(nrows)  # Skin Friction

#---------------------Evaluating Parameters for Loss---------------------------
ViewOpenRTZ(-0.783962, -0.683962, 0.508295, 0.608295)
QntFieldDerived(0, 'swirl', 'sqrt((x*x)+(y*y))*Vt', '', '0')
QntFieldDerived(0, 'r', 'sqrt((x*x)+(y*y))', '', '0')
QntFieldDerived(0, 'Wui', 'Vt - r*2340 - Wt', '', '0')
QntFieldDerived(0, 'Vui', 'Wt + r*2340 - Vt', '', '0')
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
    if c < nsect - 1:
        plne_nme = 'CUT' + str(c + len(span) + 1)
    else:
        plne_nme = 'CUT' + str(c + len(span) + 1) + '.D1'
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
    QntFieldScalar('atan(Wt/Wm)')
    beta[c] = WeightedIntegral()
    QntFieldScalar('Wui')
    Wui[c] = WeightedIntegral()
    QntFieldScalar('Vm')
    Vm[c] = WeightedIntegral()
    QntFieldScalar('Absolute Mach Number')
    M[c] = WeightedIntegral()
    QntFieldScalar('Relative Mach Number')
    Mrel[c] = WeightedIntegral()
    QntFieldScalar('Static Temperature')
    T[c] = WeightedIntegral()

# ==============================================================================
# Cartesian view - Calculating Loss and Exporting to File(Surface Area Average Values)
# ==============================================================================
ViewOpen(-0.504382, -0.404382, 0.154613, 0.254613)
LimitsFull()


# Evaluating blade properties averaged
QntSolidScalar('Cf')
for i in range(nrows):
    if i % 2 == 0:  # Check for rotor row number
        SelectFromProject('row ' + str(i + 1) + '_blade_(r.p.m. 22363)')
    else:  # Else a stator
        SelectFromProject('row ' + str(i + 1) + '_blade')
    cf[i] = SclAverage()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
DH = np.zeros(nrows)
Rx = np.zeros(nrows)
U = w * rm
r_inc = np.linspace(r[0][0], r[0][1], 10)
U_inc = np.average(w * r_inc)
f_inc = 0.5

c = 0
for i in range(nrows):
    # Aerodynamic Enthalpy
    dH_Euler[i] = Cp * (T0[c + 5] - T0[c])
    c += 2

# Diffusion Factor

DH[0] = W[-2] / W[0]
Wratio = W[5] / W_tip[0]
sol_r2 = Z[0] / np.pi
rad_ratio = r[0][1] / (0.5 * (r[5][0] + r[5][1]))
num = 0.75 * dH_Euler[0] / U[5]**2
Df[0] = 1 - Wratio + (Wratio * num) / \
    (sol_r2 * (1 - rad_ratio) + 2 * rad_ratio)

c = 0
Rx[0] = (P[-2] - P[0])/(P[-1] - P[0])
for i in range(nrows):
    # Incidence Loss
    dH_inc[i] = f_inc * (Wt_ac[i] + Wt[c])**2 / 2

    """Losses calculated assuming all loss is similar to rotor 2 loss correlations"""
    # Blade Loading Loss
    dH_bld[i] = 0.05 * Df[i] * U[c + 5]**2

    # Skin Friction Loss
    nu2 = r[c][0] / r[c][1]
    L = 0.5 * (1 - rm[c] / rm[c + 5]) / np.cos(np.deg2rad(betab[c]))
    p = Z[i] / (np.pi * np.cos(np.deg2rad(betab[c])))
    bi = 2 * rm[c + 5] * Z[i] / 2 * np.pi * \
        r[i][1] * np.cos(np.deg2rad(betabi[c]))
    dia_rat = r[i][1] / rm[c + 5]
    denom = (2 / (1 - nu2)) + (2 * Z[i] * np.sqrt(1 + 0.5 *
                                                 (1 + nu2**2) * np.tan(np.radians(betabi[i]))**2))/(np.pi*(1 + nu2))
    dhyd = (1 / (p + bi)) + dia_rat/denom
    W_avg = (W[c]**2 + W[c + 5]**2) / 2
    dH_sf[i] = 2 * cf[i] * L * W_avg / dhyd

    # Disk Friction
    Re_df = U[c + 5] * rm[c + 5] / nu
    f_df = 0.0622 / Re_df**0.2
    rho_avg = (rho[c] + rho[c + 5]) / 2
    dH_df[i] = f_df * (rho_avg * rm[c + 5]**2 * U[c + 5]**3) / (4 * m)

    # Recirculation Loss
    dH_rc[i] = abs(8e-5 * math.sinh(3.5 * alpha[c + 5]**3)
                   * Df[i]**2 * U[c + 5]**2)

    # Clearance Loss for rotor 2
    sold = (4 * np.pi) / (bw[c + 5] * Z[i])
    frac = ((r[c][1]**2 - r[c][0]**2)) / \
        ((r[c + 5][0] - r[c][1]) * (1 + rho[c + 5] / rho[c]))
    dH_cl[i] = 0.6 * (cl[i] / bw[c + 5]) * abs(Vt[c + 5]) * \
        (sold * frac * abs(Vt[c + 5]) * Vm_m[c])**0.5

    # Leakage Loss for Rotors only
    if i % 2 == 0:
        b_avg = (bw[c] + bw[c + 5]) / 2
        r_avg = 0.5 * (r[c + 5][0] + 0.5 * (r[c][0] + r[c][1]))
        r1_m = 0.5 * (r[c][0] + r[c][1])
        dP_cl = m * (r[c + 5][0] * abs(Vt[c + 5]) - r_avg *
                     abs(Vt[c])) / (Z[i] * r_avg * b_avg)
        U_cl = 0.816 * (2 * dP_cl * rho[c + 5])**0.5
        m_cl = rho[c + 5] * Z[i] * cl[i] * Lb[i] * U_cl
        dH_lk[i] = m_cl * U_cl * U6 / (2 * m)

    c += 2

# Vaneless diffuser loss
dH_vld = Cp * T0[5] * ((P[6] / P0[6])**(1 / 3.5) -
                       (P[6] / P0[5])**(1 / 3.5))

dH_int = dH_inc + dH_bld + dH_sf + dH_cl + dH_vld  # Internal Loss
dH_par = dH_rc + dH_df + dH_lk  # Exernal Loss
dH_act = dH_Euler + dH_par  # Actual Work Done
Eff_isen = (dH_Euler - dH_int) / dH_act

# ==============================================================================
# Exporting data to file
# ==============================================================================

sys.stdout = open(file_dir + LossFile, "w")

Rowheaders = ['J', 'r', 'Swirl', 'Vt', 'Vm', 'T', 'P', 'M', 'Mrel',
              'T0', 'P0', 'P0rel', 'alpha_m', 'beta_m', 'Entropy']
mainlist = [[] for i in range(nsect + 1)]
for i in range(nsect + 1):
    if i == 0:
        for j in range(len(Rowheaders)):
            mainlist[i].append(Rowheaders[j])
    else:
        mainlist[i].append(i)
        mainlist[i].append('%.4f' % float(rm[i - 1]))
        mainlist[i].append('%.4f' % float(sw[i - 1]))
        mainlist[i].append('%.4f' % float(Vt[i - 1]))
        mainlist[i].append('%.4f' % float(Vm[i - 1]))
        mainlist[i].append('%.4f' % float(T[i - 1]))
        mainlist[i].append('%.4f' % float(P[i - 1]))
        mainlist[i].append('%.4f' % float(M[i - 1]))
        mainlist[i].append('%.4f' % float(Mrel[i - 1]))
        mainlist[i].append('%.4f' % float(T0[i - 1]))
        mainlist[i].append('%.4f' % float(P0[i - 1]))
        mainlist[i].append('%.4f' % float(P0rel[i - 1]))
        mainlist[i].append('%.4f' % float(alpha[i - 1] * 180 / np.pi))
        mainlist[i].append('%.4f' % float(beta[i - 1] * 180 / np.pi))
        mainlist[i].append('%.4f' % float(s[i - 1]))
print('\n' + tabulate(mainlist) + '\n')

Rowheaders = ['Row', 'Cf', 'Df', 'Rx', 'DeHaller', 'Euler_Work'] + \
    Qnt_LossParm + ['Efficiency']
mainlist = [[] for i in range(nrows + 1)]
for i in range(nrows + 1):
    if i == 0:
        for j in range(len(Rowheaders)):
            mainlist[i].append(Rowheaders[j])
    else:
        mainlist[i].append(i)
        mainlist[i].append('%.4f' % float(cf[i - 1]))
        mainlist[i].append('%.4f' % float(Df[i - 1]))
        mainlist[i].append('%.4f' % float(Rx[i - 1]))
        mainlist[i].append('%.4f' % float(DH[i - 1]))
        mainlist[i].append('%.4f' % float(dH_Euler[i - 1]))
        mainlist[i].append('%.4f' % float(dH_inc[i - 1]))
        mainlist[i].append('%.4f' % float(dH_bld[i - 1]))
        mainlist[i].append('%.4f' % float(dH_sf[i - 1]))
        mainlist[i].append('%.4f' % float(dH_cl[i - 1]))
        mainlist[i].append('%.4f' % float(dH_rc[i - 1]))
        mainlist[i].append('%.4f' % float(dH_df[i - 1]))
        mainlist[i].append('%.4f' % float(dH_lk[i - 1]))
        mainlist[i].append('%.4f' % float(dH_vld))
        mainlist[i].append('%.4f' % float(Eff_isen[i - 1]))
print('\n' + tabulate(mainlist))

phi = (m / rho[5]) / (U[5] * ((2 * r[5][0])**2))
print("\n" + "Flow Coefficient")
print(phi)
sys.stdout.close()
