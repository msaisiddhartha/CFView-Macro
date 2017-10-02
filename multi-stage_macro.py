CFViewBackward(912)
"""This macro provides CFView post-processing for multi-stage compressor.
It computes:
    1. Contour Plots at different spanwise sections
    2. Cartesian plots of X-Gradient of Entropy vs X-direction
    3. Cp distribution at 0.5 span
    4. Loss through correlations and the quantities are mass averaged."""

import sys
import math
import numpy as np
import pylab as py
from tabulate import tabulate

#Case Data
dalpha = 30
WorkSplitR1 = 35

project_name = 'MSD_48bladesR2'
case_name = '4kgs_FR'
file_dir = 'C:/Users/msais/Box Sync/Thesis Work/Multi-Stage_data/DiffuserConstArea/WorkSplitRotor1='+str(WorkSplitR1)+'/Stator'+str(dalpha)+'deg/'+project_name+'/'+project_name+'_'+case_name+'/'
RunFile = str(project_name+'_'+case_name+'.run')
LossFile = "LossData_cyl.txt"

#----Input data----
Cp = 1006							#Specific Heat
R = 287								#Universal gas constant
N = 22363							#Shaft Speed
m = 4								#mass flow rate
r2 = 0.2							#Outlet radius
L = [0.05247181, 0.025, 0.0775]       #Hydraulic Length
Dh = [0.05640476, 0.03375, 0.010875]	#Hydraulic Diameter
cl = [0.0005, 0.0004, 0.0003]         #Clearance for each row

#Reference values
Vref = 100.74						#velocity (Inlet velocity)
rho = 1.2							#density
Pref = 101325						#Pressure

w = 2*np.pi*N/60					#Angular Velocity
U2 = r2*w
phi = (m/rho)/(U2*((2*r2)**2))

#------Input for Post-processing ------
#Loss file data
nrows = 3							#Number of blade rows
nsect = 7           				#Number of planes to create

#3D-View Contour Plot
span = [0.5,0.9]
EntropyRange =[[-10,120],[-10,170]]

#----Initializing Parameters----
#Scalar data
Vt = np.zeros(nsect)				#Tangential Velocity
T0 = np.zeros(nsect)				#Total Temperature
P0 = np.zeros(nsect)				#Total Pressure
P0rel = np.zeros(nsect)				#Relative Total Pressure
sw = np.zeros(nsect)				#Swirl
r = np.zeros(nsect)					#Radius at mid-span of section plane
s =  np.zeros(nsect)				#Entropy
cf =  np.zeros(nrows)				#Skin Friction

#For Loss calulation
dsext = np.zeros(nrows)				#Change in Entropy For External Loss
dsint = np.zeros(nrows)				#Change in Entropy For Internal Loss
dssf = np.zeros(nrows)				#Change in Entropy For Skin Friction Loss
dsmix = np.zeros(nrows)				#Change in Entropy For Mixing Loss
lambdaext = np.zeros(nrows)			#Head Loss Coefficient For External Loss
lambdaint = np.zeros(nrows)			#Head Loss Coefficient For Internal Loss
lambdasf = np.zeros(nrows)			#Head Loss Coefficient For Skin Friction Loss
lambdamix = np.zeros(nrows)			#Head Loss Coefficient For Mixing Loss
fc = np.zeros(nrows)				#Correction Factor

Qnt_Cart = {'ds/dx':'Grad (Entropy)_X',
       'Cp':'(Static Pressure-101325)/(0.5*1.2*(Magnitude of V*Magnitude of V+(Vt-Wt)*(Vt-Wt)))'}


FileOpenProject(file_dir+RunFile,)

#%%
# ==============================================================================
# 3D View
# ==============================================================================

#-----------------------Contour Plots at Differnet Span-------------------------

# R1CamPos = [0.261069,-0.0804571,0.0386409,0.158219,-0.0247637,0.0239268,-0.474188,-0.880254,-0.0172853,0.0471529,0.0703255]
# S1CamPos = [0.214698, -0.0509236, 0.00635897, 0.16529, -0.00694214, 0.0328556, -0.26149, -0.69639, 0.66833, 0.0285029, 0.0425103]
# R2CamPos = [0.360388,0.0561327,-0.0619866,0.215378,0.0653739,0.0752378,0.285215,-0.887805,0.361185,0.0799439,0.119231]

SetNumecaLogo(0,0)
for i in range(len(span)):
    CutPlaneSave(span[i], 0, 0, 1, 0, 0, 2)
GmtRepetitionToggle()
GmtRepetitionNumber(3, 3, 3)
scm =[0.222479, -0.0683883, -0.00803093, 0.158146, -0.0111208, 0.0264699, -0.26149, -0.69639, 0.66833, 0.0371132, 0.0553519]
SetCamera(scm[0],scm[1],scm[2],scm[3],scm[4],scm[5],scm[6],scm[7],scm[8],scm[9],scm[10],)
GmtToggleBoundary()


for i in range(len(span)):
    SelectFromProject('CUT'+str(i+1),)

    # Entropy
    QntFieldScalar('Entropy',)
    SclContourStrip()
    ColormapStripesOnly()
    ColormapNumOfSmallTicks(3)
    ColormapTicksNumberTextType(10,12,2,0,1,0,0,1,0,0,0,0)
    ColormapLabelTextType(10,12,2,2,1,0,0,1,0,0,0,0)
    RprRangeIn(EntropyRange[i][0],EntropyRange[i][1])
    SetCamera(0.261069,-0.0804571,0.0386409,0.158219,-0.0247637,0.0239268,-0.474188,-0.880254,-0.0172853,0.0471529,0.0703255)
    Print(8,0,0,1,100,1317,704,0 ,file_dir+'EntropyB2b'+str(span[i])+'R1_'+str(WorkSplitR1)+'_'+str(dalpha)+'.png' ,'',1,1,1)
    SetCamera(0.214698, -0.0509236, 0.00635897, 0.16529, -0.00694214, 0.0328556, -0.26149, -0.69639, 0.66833, 0.0285029, 0.0425103)
    Print(8, 0, 0, 1, 100, 1920, 1080, 0, file_dir+'EntropyB2b'+str(span[i])+'S1_'+str(WorkSplitR1)+'_'+str(dalpha)+'.png', '', 1, 1, 1)
    SetCamera(0.360388,0.0561327,-0.0619866,0.215378,0.0653739,0.0752378,0.285215,-0.887805,0.361185,0.0799439,0.119231)
    Print(8,0,0,1,100,1701,873,0 ,file_dir+'EntropyB2b'+str(span[i])+'R2_'+str(WorkSplitR1)+'_'+str(dalpha)+'.png' ,'',1,1,1)

    # Relative Velocity
    QntFieldScalar('Magnitude of W',)
    SclContourStrip()
    ColormapStripesOnly()
    ColormapLabelTextType(10,12,2,2,1,0,0,1,0,0,0,0)
    SclContourStrip()
    RprRangeIn(0,400)
    SetCamera(0.261069,-0.0804571,0.0386409,0.158219,-0.0247637,0.0239268,-0.474188,-0.880254,-0.0172853,0.0471529,0.0703255)
    Print(8,0,0,1,100,1317,704,0 ,file_dir+'RelVelB2b'+str(span[i])+'R1_'+str(WorkSplitR1)+'_'+str(dalpha)+'.png' ,'',1,1,1)
    SetCamera(0.214698, -0.0509236, 0.00635897, 0.16529, -0.00694214, 0.0328556, -0.26149, -0.69639, 0.66833, 0.0285029, 0.0425103)
    Print(8, 0, 0, 1, 100, 1920, 1080, 0, file_dir+'RelVelB2b'+str(span[i])+'S1_'+str(WorkSplitR1)+'_'+str(dalpha)+'.png', '', 1, 1, 1)
    SetCamera(0.360388,0.0561327,-0.0619866,0.215378,0.0653739,0.0752378,0.285215,-0.887805,0.361185,0.0799439,0.119231)
    Print(8,0,0,1,100,1701,873,0 ,file_dir+'RelVelB2b'+str(span[i])+'R2_'+str(WorkSplitR1)+'_'+str(dalpha)+'.png' ,'',1,1,1)
    DeleteAll()

#------------------------------Cartesian Plots----------------------------------

ViewNum = 2
Curvenum = 1
QntFieldScalar('Entropy')
FieldGradient('Entropy')
QntFieldVector('Grad (Entropy)')
for key in Qnt_Cart.keys():
    QntFieldDerived(0 ,key ,Qnt_Cart[key] ,'' ,'0')

for key in Qnt_Cart.keys():
    ViewActivate(RunFile+':1')
    SelectFromProject('row 1_blade_(r.p.m. -22363)','row 2_blade','row 3_blade_(r.p.m. -22363)')
    QntFieldScalar(key)
    SclPlotNormalizedGridLine(0,0.5,0,1 ,'row_1_flux_1_Main_Blade_skin.Imin blade_(aap-ps)_rotating',0,0.5,0,1 ,'row_1_flux_1_Main_Blade_skin.Imin blade_(aap-ss)_rotating',0)
    ViewActivate(RunFile+':'+str(ViewNum))
    PlotPlaneX()
    SelectPlotCurves('Gridline I=0.5 on row_1_flux_1_Main_Blade_skin.Imin blade_(aap-ps)_rotating' ,'Gridline I=0.5 on row_1_flux_1_Main_Blade_skin.Imin blade_(aap-ss)_rotating')
    PlotCurvesMerge('Gridline I=0.5 on row_1_flux_1_Main_Blade_skin.Imin blade_(aap-ps)_rotating' ,'Gridline I=0.5 on row_1_flux_1_Main_Blade_skin.Imin blade_(aap-ss)_rotating')
    ActivePlotCurveOutput(file_dir+'3DB2b0.5R1_'+str(WorkSplitR1)+'_'+str(dalpha)+'.dat' ,'merged curves '+str(Curvenum))
    Curvenum+=1
    SelectPlotCurves('Gridline I=0.5 on row_2_flux_1_Main_Blade_skin.Imin blade_(aap-ps)' ,'Gridline I=0.5 on row_2_flux_1_Main_Blade_skin.Imin blade_(aap-ss)')
    PlotCurvesMerge('Gridline I=0.5 on row_2_flux_1_Main_Blade_skin.Imin blade_(aap-ps)' ,'Gridline I=0.5 on row_2_flux_1_Main_Blade_skin.Imin blade_(aap-ss)')
    ActivePlotCurveOutput(file_dir+'3DB2b0.5S1_'+str(WorkSplitR1)+'_'+str(dalpha)+'.dat' ,'merged curves '+str(Curvenum))
    Curvenum+=1
    SelectPlotCurves('Gridline I=0.5 on row_3_flux_1_Main_Blade_skin.Imin blade_(aap-ps)_rotating' ,'Gridline I=0.5 on row_3_flux_1_Main_Blade_skin.Imin blade_(aap-ss)_rotating')
    PlotCurvesMerge('Gridline I=0.5 on row_3_flux_1_Main_Blade_skin.Imin blade_(aap-ps)_rotating' ,'Gridline I=0.5 on row_3_flux_1_Main_Blade_skin.Imin blade_(aap-ss)_rotating')
    ActivePlotCurveOutput(file_dir+'3DB2b0.5R2_'+str(WorkSplitR1)+'_'+str(dalpha)+'.dat' ,'merged curves '+str(Curvenum))
    Curvenum+=1
    DeletePlot()
    ViewNum+=1


#==============================================================================
#Cylindrical view - Calculating Loss and Exporting to File(Surface Area Average Values)
#==============================================================================
ViewActivate(RunFile+':1')
sys.stdout = open(file_dir+LossFile,"w")

#---------------------Evaluating Parameters for Loss---------------------------
ViewOpenRTZ(-0.783962, -0.683962, 0.508295, 0.608295)
QntFieldDerived(0 ,'swirl' ,'sqrt((x*x)+(y*y))*Vt' ,'' ,'0')
QntFieldDerived(0 ,'r' ,'sqrt((x*x)+(y*y))' ,'' ,'0')
LimitsFull()
ViewOriginal(1,)
CutPlaneSave(0.11,0,0,0,0,1,1)    #Rotor 1 Inlet
CutPlaneSave(0.0760,0,0.05247182,0.27081969,0,1,1)   #Rotor 1 Outlet
CutPlaneSave(0.0765,0,0.05547558,0.27442423,0,1,1)  #Stator 1 Inlet
CutPlaneSave(0.1185,0,0.09072516,0.65675436,0,1,1)  #Stator 1 Outlet
CutPlaneSave(0.1215,0,0.09387836,0.71827631,0,1,1)  #Rotor 2 Inlet
CutPlaneSave(0.2001,0,0.12,1,0,0,1)    #Rotor 2 Outlet
CutPlaneSave(0.3,0,0.12,1,0,0,1)    #Diffuser Outlet
for c in range(nsect):
    plne_nme = 'CUT'+str(c+3)
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

SelectFromProject('CUT9.D17')
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

print('\nJ    Swirl      T0          P0      P0rel       r     Entropy')

for c in range(nsect):
    print(str(c+1)+'    '+'%.4f'%(sw[c])+'    '+'%.4f'%(T0[c])+'    '+'%.4f'%(P0[c])+'    '+'%.4f'%(P0rel[c])+'    '+'%.4f'%(r[c])+'    '+'%.4f'%(s[c]))

#Evaluating blade properties averaged
QntSolidScalar('Cf')
for i in range(nrows):
    if i%2==0:                  #Check for rotor row number
        SelectFromProject('row '+str(i+1)+'_blade_(r.p.m. -22363)')
    else:                       #Else a stator
        SelectFromProject('row '+str(i+1)+'_blade')
#    cf[i]=WeightedIntegral()
    cf[i]=SclAverage()

#---------------------Calculating Individual Losses---------------------------
c =0
for i in range(nrows):
    #External loss
    dsext[i] = abs(Cp*np.log(T0[c+1]/(T0[c]+((sw[c+1]-sw[c])*w)/Cp)))
    lambdaext[i] = T0[-1]*dsext[i]/(U2)**2
    #Inlernal loss
    dsint[i] = abs(Cp*np.log(1+(sw[c+1]-sw[c])*w/(Cp*T0[c]))-R*np.log(P0[c+1]/P0[c]))
    lambdaint[i] = T0[-1]*dsint[i]/(U2)**2
    #Skin friction loss
    fc[i] = P0rel[c+1]/P0rel[c]    #Correction factor
    dssf[i] = fc[i]*(0.5*Vref**2)*(4*cf[i]*L[i])/(Dh[i]*T0[-3])
    lambdasf[i] = T0[-1]*dssf[i]/(U2)**2
    #Mixing loss
    dsmix[i] = -R*np.log(P0[c+2]/P0[c+1])
    lambdamix[i] = T0[-1]*dsmix[i]/(U2)**2
    c+=2

#Incidence Loss
dsinc = dsint-abs(dssf)
#Diffuser loss
dsdiff = -R*np.log(P0[-1]/P0[-2])
lambdadiff = T0[-1]*dsdiff/(U2)**2

ds = dsext+dsint
lambd = T0[-1]*ds/(U2)**2

Rowheaders = ['Row','Cf','SkinFricLoss','IncLoss','MixLoss','IntLoss','ExtLoss','Overall(Int+Ext)']
mainlist = [[] for i in range(len(Rowheaders))]
for i  in range(nrows+1):
    if i==0:
        for j in range(len(Rowheaders)):
            mainlist[i].append(Rowheaders[j])
    else:
        mainlist[i].append(i)
        mainlist[i].append('%.4f'%float(cf[i-1]))
        mainlist[i].append('%.4f'%float(dssf[i-1]))
        mainlist[i].append('%.4f'%float(dsinc[i-1]))
        mainlist[i].append('%.4f'%float(dsmix[i-1]))
        mainlist[i].append('%.4f'%float(dsint[i-1]))
        mainlist[i].append('%.4f'%float(dsext[i-1]))
        mainlist[i].append('%.4f'%float(ds[i-1]))
print('\n'+tabulate(mainlist))

sys.stdout.close()

#%%
#==============================================================================
#Open Turbomachinery Mode
#==============================================================================
SetTurboMode()
OpenTurboModeStandard3DView()
OpenTurboModeBladeToBladeView()
OpenTurboModeBladeView()
OpenTurboModeMeridionalView(0,0)

#==============================================================================
#Blade to blade plots
#==============================================================================
ViewActivate(RunFile+':'+str(ViewNum))
GmtToggleBoundary()
QntFieldScalar('Absolute Total Pressure')
SclContourStrip()
ViewNum+=1
#==============================================================================
#Blade view
#==============================================================================
ViewActivate(RunFile+':'+str(ViewNum))
GmtToggleBoundary()
QntFieldScalar('Absolute Total Pressure')
SclContourStrip()
ViewNum+=1
#==============================================================================
#Azimuthal Plots
#==============================================================================
ViewActivate(RunFile+':1')
ViewActivate(':'+str(ViewNum))
#ViewActivate(RunFile+'.run meridional view:4')
GmtToggleBoundary()
QntFieldScalar('Static Pressure')
SclContourStrip()
ViewNum+=1
