# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 6.13-2 replay file
# Internal Version: 2013_07_18-07.22.41 126428
# Run by HPMSL on Tue Feb 02 13:06:29 2016
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from time import *
from os import *
from numpy import *
import numpy as np
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import xlwt
import xlrd
import odbAccess



# input vriables
#name_of_job = '2D-1-{0}'.format(0) 

i_start=0

Loading_day = 
Load_amount = 
amount_day = 

n_callus_elem=3394
type_of_elem = 'CAX4R'

Em = 1e6
Eib = 1e9
Emb = 6e9
Ec = 15e6
Ef = 5e6

E1_c = Em*ones(n_callus_elem)
E2_c = Em*ones(n_callus_elem)
E3_c = Em*ones(n_callus_elem)
E4_c = Em*ones(n_callus_elem)
E5_c = Em*ones(n_callus_elem)
E6_c = Em*ones(n_callus_elem)
E7_c = Em*ones(n_callus_elem)
E8_c = Em*ones(n_callus_elem)
E9_c = Em*ones(n_callus_elem)
E10_c = Em*ones(n_callus_elem)
Enew_c = Em*ones(n_callus_elem)
E_c = Em*ones(n_callus_elem)
nu_c = 0.17*ones(n_callus_elem)
D_c = 0.5*(2.37e-6)*ones(n_callus_elem)
cm_c = 0*ones(n_callus_elem)
cm_p = 0*ones(n_callus_elem)
cib_c = 0*ones(n_callus_elem)
cmb_c = 0*ones(n_callus_elem)
cc_c = 0*ones(n_callus_elem)
cf_c = 0*ones(n_callus_elem)
ctot = 0*ones(n_callus_elem)

path = 'C:/Masters Project/Full/Day1/Load1e5/'

if i_start == 0:
    print('lets start')
else:

    
    Ex = xlrd.open_workbook('E_elem_{0}.xls'.format(i_start))
    for sheet in Ex.sheets():
        n_row=sheet.nrows
        n_col=sheet.ncols
        
        for i in range(n_row):
            
            E1_c[i] = (sheet.cell(i,5).value)
            E2_c[i] = (sheet.cell(i,6).value)
            E3_c[i] = (sheet.cell(i,7).value)
            E4_c[i] = (sheet.cell(i,8).value)
            E5_c[i] = (sheet.cell(i,9).value)
            E6_c[i] = (sheet.cell(i,10).value)
            E7_c[i] = (sheet.cell(i,11).value)
            E8_c[i] = (sheet.cell(i,12).value)
            E9_c[i] = (sheet.cell(i,13).value)
            E10_c[i] = (sheet.cell(i,14).value)
            
            
            Enew_c[i] = (sheet.cell(i,16).value)
            E_c[i] = (sheet.cell(i,17).value)
            nu_c[i] = (sheet.cell(i,18).value)
            D_c[i] =(sheet.cell(i,19).value)
            cm_c[i] = (sheet.cell(i,20).value)
            cib_c[i] = (sheet.cell(i,21).value)
            cmb_c[i] = (sheet.cell(i,22).value)
            cc_c[i] = (sheet.cell(i,23).value)
            cf_c[i] = (sheet.cell(i,24).value)
            ctot[i] = (sheet.cell(i,25).value)
        
    

#E_for_excel = 0*ones(n_callus_elem)



f_bone=0.15
f_rbone=0.1
f_cbone=0.1
f_cart=0.1
f_fiber=0.1
f_excd=0.5


#book = xlwt.Workbook()

for iter in range(amount_day):
    book = xlwt.Workbook()
    i_run = iter + i_start
    print('number of iterations',i_run+1)
    name_of_job = '2D-stress-d=6-{0}'.format(i_run+1)
# openMdb('Model 1.cae')
        
                
    mdb.ModelFromInputFile(name=name_of_job,inputFileName='{1}Job-2D-stress-d=6-{0}.inp'.format(i_run,path))
#: The model "2D-5" has been imported from an input file. 
#: Please scroll up to check for error and warning messages.
    
    a = mdb.models[name_of_job].rootAssembly
    p1 = mdb.models[name_of_job].parts['PART-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)
    
#%% THIS IS THE PART I NEED TO EDIT

    if i_run == Loading_day:
        print('Load applied')
        mdb.models[name_of_job].Pressure(amplitude=UNSET, createStepName='Step-1', 
            distributionType=UNIFORM, field='', magnitude=Load_amount, name='AppliedLoad'
            , region=mdb.models[name_of_job].rootAssembly.surfaces['SURF-1']) 

#%%
# Change the material properties for each element
    
    
    for ii in range(n_callus_elem):
        EE = E_c[ii]
        nuu = nu_c[ii]
        mdb.models[name_of_job].materials['MAT-CALLUS-{0}'.format(ii+1)].elastic.setValues(table=((EE, nuu), ))
    
# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 6.13-2 replay file
# Internal Version: 2013_07_18-07.22.41 126428
# Run by HPMSL on Fri Jun 03 13:08:36 2016
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
    mdb.models[name_of_job].fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'EE', 'LE', 'FLVEL', 'U', 'RF', 'COORD'))

# Creat of job
    mdb.Job(name='Job-{0}'.format(name_of_job), model=name_of_job, description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
    
    
    
# Submission of Job

    mdb.jobs['Job-{0}'.format(name_of_job)].submit(consistencyChecking=OFF)
    mdb.jobs['Job-{0}'.format(name_of_job)].waitForCompletion()

# openning odb file

    o1 = session.openOdb(name='{1}Job-{0}.odb'.format(name_of_job,path))
    session.viewports['Viewport: 1'].setValues(displayedObject=o1)
    session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
    session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(deformationScaling=UNIFORM, uniformScaleFactor=1)
    session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(maxAutoCompute=ON)
    session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(minAutoCompute=ON)
# Read from odb

# Define some handy variables

    odb = openOdb(path='{1}Job-{0}.odb'.format(name_of_job,path))

    step1 = odb.steps['Step-1']

# Last frame in step

    frame = step1.frames[-1]

# all Strain in model

    Strain = frame.fieldOutputs['EE']
    Stress = frame.fieldOutputs['S']
    FlVel = frame.fieldOutputs['FLVEL']
    
# Defines elements = elements in elem set of CALLUS

    elems = odb.rootAssembly.instances['PART-1-1'].elementSets['WHOLE_CALLUS']

    
# get the coordinates
    
# Sets Strain_at_elems = Strains in elements in set 'WHOLE_CALLUS'

    Strn_at_elems = Strain.getSubset(region=elems)
    Strs_at_elems = Stress.getSubset(region=elems)
    FlVel_at_elems = FlVel.getSubset(region=elems)
    
    n_elem = len(Strn_at_elems.values)
    n_elem = n_callus_elem
    E11=zeros((n_elem))
    E22=zeros((n_elem))
    E33=zeros((n_elem))
    E12=zeros((n_elem))
    E1=zeros((n_elem))
    E2=zeros((n_elem))
    E3=zeros((n_elem))
    E_oct=zeros((n_elem))
    S11=zeros((n_elem))
    S22=zeros((n_elem))
    S33=zeros((n_elem))
    S12=zeros((n_elem))
    S1=zeros((n_elem))
    S2=zeros((n_elem))
    S3=zeros((n_elem))
    FlVel1=zeros((n_elem))
    FlVel2=zeros((n_elem))
    FF=zeros((n_elem))
    Hydro_P=zeros((n_elem))
    Sai=zeros((n_elem))
    ie = 0
    iie= -1
    iie2 = -1
    for StrnVal in Strn_at_elems.values:
        iie = iie + 1
        iie2 = iie2 + 1
        if iie2 == 4:
            ie = ie + 1
            iie2 = 0
        E11[ie] = E11[ie] + 0.25*StrnVal.data[0]
        E22[ie] = E22[ie] + 0.25*StrnVal.data[1]
        E33[ie] = E33[ie] + 0.25*StrnVal.data[2]
        E12[ie] = E12[ie] + 0.25*StrnVal.data[3]
        
        E1[ie] = E1[ie] + 0.25*StrnVal.minPrincipal
        E2[ie] = E2[ie] + 0.25*StrnVal.midPrincipal
        E3[ie] = E3[ie] + 0.25*StrnVal.maxPrincipal
        
        E_oct[ie] = 2.0/3*(((E1[ie] - E2[ie])**2 + (E1[ie] - E3[ie])**2 + (E3[ie] - E2[ie])**2)**0.5)
        
        Sai[ie] = max(abs(E1[ie]),abs(E2[ie]),abs(E3[ie]))
        #print(StrnVal)
        #print(ali)
#        Sai[ie] = ( (E1[ie] - E_oct[ie])**2 + (E2[ie] - E_oct[ie])**2 + (E3[ie] - E_oct[ie])**2 )**0.5
    
    
    ie = 0
    iie= -1
    iie2 = -1
    for StrsVal in Strs_at_elems.values:
        iie = iie + 1
        iie2 = iie2 + 1
        if iie2 == 4:
            ie = ie + 1
            iie2 = 0
        S11[ie] = S11[ie] + 0.25*StrsVal.data[0]
        S22[ie] = S22[ie] + 0.25*StrsVal.data[1]
        S33[ie] = S33[ie] + 0.25*StrsVal.data[2]
        S12[ie] = S12[ie] + 0.25*StrsVal.data[3]
        S1[ie] = S1[ie] + 0.25*StrsVal.minPrincipal
        S2[ie] = S2[ie] + 0.25*StrsVal.midPrincipal
        S3[ie] = S3[ie] + 0.25*StrsVal.maxPrincipal
        
        Hydro_P[ie] = Hydro_P[ie] + 0.25*StrsVal.press
    
    
        #S11[ie] = CoordVal.data[0]
#    odb.close()
    
    ie = 0
    iie= -1
    iie2 = -1
    for FlVelVal in FlVel_at_elems.values:
        iie = iie + 1
        iie2 = iie2 + 1
        if iie2 == 4:
            ie = ie + 1
            iie2 = 0
        FlVel1[ie] = FlVel1[ie] + 0.25*FlVelVal.data[0]
        FlVel2[ie] = FlVel2[ie] + 0.25*FlVelVal.data[1]
        FF[ie] = ((FlVel1[ie])**2 + (FlVel2[ie])**2)**0.5



    
# %%%%%%%%%%%%%%%%   Diffusion analysis  %%%%%%%%%%%%%%%%%%%%%%%%%     
 # %%%%%%%%%%%%%%%%   Diffusion analysis  %%%%%%%%%%%%%%%%%%%%%%%%%     
    name_of_job = '2D-diffusion-d=6-{0}'.format(i_run+1)
# openMdb('Model 1.cae')
    mdb.ModelFromInputFile(name=name_of_job,inputFileName='{1}Job-2D-diffusion-d=6-{0}.inp'.format(i_run,path))
    a = mdb.models[name_of_job].rootAssembly
    p1 = mdb.models[name_of_job].parts['PART-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)
    
# Change the material properties for each element
    
    
    for ii in range(n_callus_elem):
        DD = D_c[ii]
        cc = cm_c[ii]
        del mdb.models[name_of_job].materials['MAT-CALLUS-{0}'.format(ii+1)].elastic
        mdb.models[name_of_job].materials['MAT-CALLUS-{0}'.format(ii+1)].Diffusivity(table=((DD,cc), ))
        mdb.models[name_of_job].materials['MAT-CALLUS-{0}'.format(ii+1)].Solubility(table=((1.0, ),))
        
    mdb.models[name_of_job].fieldOutputRequests['F-Output-1'].setValues(variables=('CONC', 'NNC', 'COORD'))
# Creat of job
    mdb.Job(name='Job-{0}'.format(name_of_job), model=name_of_job, description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
    
    
# application of initial Etration

    nodes = mdb.models[name_of_job].rootAssembly.sets['WHOLE_NODE']
    
    
    
    elem = mdb.models[name_of_job].rootAssembly.sets['PART-1-1.WHOLE_CALLUS']
    #print(elem.elements[0].connectivity[0])
    
    
    n_callus_elem = len(elem.elements)
    print(n_callus_elem)
    n_callus_node = len(nodes.nodes)
    print(n_callus_node)
    s_concen_nodes=0*ones((n_callus_node,2))
    concen_nodes=0*ones(n_callus_node)
    
    
    for jj in range(n_callus_elem):
        for kk in range(4):
            i_n = elem.elements[jj].connectivity[kk]
            #print(i_n)
            #print(jj)
            s_concen_nodes[i_n,0]= s_concen_nodes[i_n,0] + cm_c[jj]
            s_concen_nodes[i_n,1]= s_concen_nodes[i_n,1] + 1
            
        
        
    
    for ii in range(n_callus_node):
        concen_nodes[ii]=s_concen_nodes[ii,0]/s_concen_nodes[ii,1]
        
    #print(concen_nodes)
    
    iee=0;
    if i_run==0:
        for ii in range(4400):
            
            try:
                region = mdb.models[name_of_job].rootAssembly.sets['ND_{0}'.format(ii+1)]
                mdb.models[name_of_job].ConcentrationBC(name='BC-{0}'.format(ii+1), createStepName='Step-1', region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', magnitude=concen_nodes[iee], amplitude=UNSET)
                mdb.models[name_of_job].boundaryConditions['BC-{0}'.format(ii+1)].deactivate('Step-2')
                iee=iee+1;
                
            except:
                pass
    else:
        for ii in range(n_callus_node):
            mdb.models[name_of_job].boundaryConditions['Con-BC-{0}'.format(ii+1)].setValues(magnitude=concen_nodes[ii])

			
			
    if i_run==0:
        mdb.models[name_of_job].steps['Step-2'].setValues(timePeriod=7.0, maxInc=7.0)
    else:
        mdb.models[name_of_job].steps['Step-2'].setValues(timePeriod=1.0, maxInc=1.0)

            
# Submission of Job

    mdb.jobs['Job-{0}'.format(name_of_job)].submit(consistencyChecking=OFF)
    mdb.jobs['Job-{0}'.format(name_of_job)].waitForCompletion()

# openning odb file

    o1 = session.openOdb(name='{1}Job-{0}.odb'.format(name_of_job,path))
    session.viewports['Viewport: 1'].setValues(displayedObject=o1)
    session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
    session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(deformationScaling=UNIFORM, uniformScaleFactor=1)
    session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(maxAutoCompute=ON)
    
# Define some handy variables

    odb = openOdb(path='{1}Job-{0}.odb'.format(name_of_job,path))

    step1 = odb.steps['Step-2']

# Last frame in step

    frame = step1.frames[-1]

# all Strain in model

    Concentration = frame.fieldOutputs['CONC']

# Defines elements = elements in elem set of CALLUS

    elems = odb.rootAssembly.instances['PART-1-1'].elementSets['WHOLE_CALLUS']


# Sets Strain_at_elems = Strains in elements in set 'WHOLE_CALLUS'

    Concentration_at_elems = Concentration.getSubset(region=elems,position=CENTROID)
    

    n_elem = len(Concentration_at_elems.values)
    
    
    system('cls')
    ie=-1
    for ConcntrnVal in Concentration_at_elems.values:
        ie = ie + 1
        #print(ie)
        #print(len(ConcntrnVal.data))
        #print(ConcntrnVal.data)
        cm_c[ie] = ConcntrnVal.data
        
    Density_max = max(cm_c);
    
#get the coordinates of nodes in elements
    
    Coord1=0*ones(n_callus_elem)
    Coord2=0*ones(n_callus_elem)
    
    nnn=len(odb.steps['Step-2'].frames[1].fieldOutputs['COORD'].values[1].instance.nodes)
    
    
    Coordin = odb.steps['Step-2'].frames[-1].fieldOutputs['COORD'].getSubset
    
    
    ie=-1
    for ii in range(n_callus_elem):
        ie = ie + 1
         
        Coord1[ie]=odb.steps['Step-2'].frames[1].fieldOutputs['COORD'].values[ie].instance.nodes[ie].coordinates[0]
        Coord2[ie]=odb.steps['Step-2'].frames[1].fieldOutputs['COORD'].values[ie].instance.nodes[ie].coordinates[1]
        

# Mechanobiological Regulation
# Claes et al


    

# differentiation
#    print(max(E_p))
    nu_p = nu_c
    mean_Sai = mean(Sai)
    ih = zeros((n_elem))
    for ii in range(n_callus_elem):
        ih[ii] = (E_oct[ii]*100/3.75 + FF[ii]*1e6/3)*1
        ctot[ii] = cm_c[ii] + cib_c[ii] + cmb_c[ii] + cc_c[ii] + cf_c[ii]
        if ctot[ii]>1:
            fft = 1/(0.25*ctot[ii]+1)
        else:
            fft = 1
            
        
        ymb=1/(1+1*2.71828**(15*(ih[ii]-0.533)))
        yib=1/(1+1*2.71828**(-15*(ih[ii]-0.533))+1*2.71828**(15*(ih[ii]-1)))
        yc=1/(1+1*2.71828**(-15*(ih[ii]-1))+1*2.71828**(15*(ih[ii]-3)))
        yf=1/(1+1*2.71828**(-15*(ih[ii]-3)))
        
        
                
        formed_mbone_cell = ymb*fft*f_bone * cm_c[ii]
        formed_ibone_cell = yib*fft*f_bone * cm_c[ii]
        formed_cart_cell = yc*fft*f_cart * cm_c[ii]
        formed_fiber_cell = yf*fft*f_fiber * cm_c[ii]
        
        cmb_c[ii] = cmb_c[ii] + formed_mbone_cell + ymb*f_rbone*cib_c[ii]
        cib_c[ii] = cib_c[ii] + formed_ibone_cell - ymb*f_rbone*cib_c[ii] + yc*f_cbone*cc_c[ii]
        cc_c[ii] = cc_c[ii] + formed_cart_cell - yc*f_cbone*cc_c[ii]
        cf_c[ii] = cf_c[ii] + formed_fiber_cell
        cm_c[ii] = cm_c[ii] - formed_mbone_cell - formed_ibone_cell - formed_cart_cell - formed_fiber_cell    
         
        #
# Update of Material Properties
    
    
    
    #E10_c = 0*ones(n_callus_elem)
    
     
    for ii in range(n_callus_elem):
        
        E1_c[ii] = E2_c[ii]
        E2_c[ii] = E3_c[ii]
        E3_c[ii] = E4_c[ii]
        E4_c[ii] = E5_c[ii]
        E5_c[ii] = E6_c[ii]
        E6_c[ii] = E7_c[ii]
        E7_c[ii] = E8_c[ii]
        E8_c[ii] = E9_c[ii]
        E9_c[ii] = E10_c[ii]
        
        Density_max = max(1,cib_c[ii] + cmb_c[ii] + cc_c[ii] + cf_c[ii] + cm_c[ii])
        
        if (cib_c[ii] + cmb_c[ii] + cc_c[ii] + cf_c[ii] + cm_c[ii])<1:
            Enew_c[ii] = (cib_c[ii]*Eib + cmb_c[ii]*Emb + cc_c[ii]*Ec + cf_c[ii]*Ef + (1-(cib_c[ii] + cmb_c[ii] + cc_c[ii] + cf_c[ii]))*Em)/1
        else:
            Enew_c[ii] = (cib_c[ii]*Eib + cmb_c[ii]*Emb + cc_c[ii]*Ec + cf_c[ii]*Ef + cm_c[ii]*Em)/Density_max
            
    E10_c = Enew_c 
    E_c = (E1_c + E2_c + E3_c + E4_c + E5_c + E6_c + E7_c + E8_c + E9_c + E10_c)/10
    
    
    print('you are all set for iteration {0} and you can set your i_start = {0}'.format(i_run))
    print('you can set your i_start = {0}'.format(i_run+1))
    # Definition of new field output 
    
    #odb = Odb(name='myData', analysisTitle='derived data', description='test problem', path='testWrite.odb')
    
    #nodb = odb
    #ie=-1
    #for ii in range(n_callus_elem):
    #    ie = ie + 1
        #nodb.steps['Step-1'].frames[-1].fieldOutputs['NNC11'].values[ie].data = E_c[ie]
        

    
    
    
    # Save in Excel
    sh = book.add_sheet('sheet_{0}'.format(0))
    for i,e in enumerate(E_c):
        sh.write(i,0,i+1)
        sh.write(i,1,Coord1[i])
        sh.write(i,2,Coord2[i])
        
        sh.write(i,5,E1_c[i])
        sh.write(i,6,E2_c[i])
        sh.write(i,7,E3_c[i])
        sh.write(i,8,E4_c[i])
        sh.write(i,9,E5_c[i])
        sh.write(i,10,E6_c[i])
        sh.write(i,11,E7_c[i])
        sh.write(i,12,E8_c[i])
        sh.write(i,13,E9_c[i])
        sh.write(i,14,E10_c[i])
        
        
        sh.write(i,16,Enew_c[i])
        sh.write(i,17,E_c[i])
        sh.write(i,18,nu_c[i])
        sh.write(i,19,D_c[i])
        sh.write(i,20,cm_c[i])
        sh.write(i,21,cib_c[i])
        sh.write(i,22,cmb_c[i])
        sh.write(i,23,cc_c[i])
        sh.write(i,24,cf_c[i])
        sh.write(i,25,ctot[i])
        


    name_of_job = 'final_E_{0}'.format(i_run+1)
# openMdb('Model 1.cae')
    mdb.ModelFromInputFile(name=name_of_job,inputFileName='{0}Job-2D-diffusion-d=6-E.inp'.format(path))
    a = mdb.models[name_of_job].rootAssembly
    p1 = mdb.models[name_of_job].parts['PART-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)
    
# Change the material properties for each element
    
    
    for ii in range(n_callus_elem):
        DD = 1e-20
        cc = 0
        del mdb.models[name_of_job].materials['MAT-CALLUS-{0}'.format(ii+1)].elastic
        mdb.models[name_of_job].materials['MAT-CALLUS-{0}'.format(ii+1)].Diffusivity(table=((DD,cc), ))
        mdb.models[name_of_job].materials['MAT-CALLUS-{0}'.format(ii+1)].Solubility(table=((1.0, ),))
        
    mdb.models[name_of_job].fieldOutputRequests['F-Output-1'].setValues(variables=('CONC', 'NNC', 'COORD'))
# Creat of job
    mdb.Job(name='Job-{0}'.format(name_of_job), model=name_of_job, description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
    
    
# application of initial concentration

    nodes = mdb.models[name_of_job].rootAssembly.sets['WHOLE_NODE']
    
    
    
    elem = mdb.models[name_of_job].rootAssembly.sets['PART-1-1.WHOLE_CALLUS']
    #print(elem.elements[0].connectivity[0])
    
    
    n_callus_elem = len(elem.elements)
    n_callus_node = len(nodes.nodes)
    s_E_nodes=0*ones((n_callus_node,2))
    E_nodes=0*ones(n_callus_node)
    
    
    for jj in range(n_callus_elem):
        for kk in range(4):
            i_n = elem.elements[jj].connectivity[kk]
            s_E_nodes[i_n,0]= s_E_nodes[i_n,0] + E_c[jj]
            s_E_nodes[i_n,1]= s_E_nodes[i_n,1] + 1
            
        
        
    
    for ii in range(n_callus_node):
        E_nodes[ii]=s_E_nodes[ii,0]/s_E_nodes[ii,1]
        
    #print(E_nodes)
  
    iee=0;

    for ii in range(4400):
           
        try:
            region = mdb.models[name_of_job].rootAssembly.sets['ND_{0}'.format(ii+1)]
            mdb.models[name_of_job].ConcentrationBC(name='BC-{0}'.format(ii+1), createStepName='Step-1', region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', magnitude=E_nodes[iee], amplitude=UNSET)
            iee=iee+1;
                
        except:
            pass

# Submission of Job

    mdb.jobs['Job-{0}'.format(name_of_job)].submit(consistencyChecking=OFF)
    mdb.jobs['Job-{0}'.format(name_of_job)].waitForCompletion()

# openning odb file

    o1 = session.openOdb(name='{1}Job-{0}.odb'.format(name_of_job,path))
    session.viewports['Viewport: 1'].setValues(displayedObject=o1)
    session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
    session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(deformationScaling=UNIFORM, uniformScaleFactor=1)
    #session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(maxAutoCompute=ON)
    session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(maxAutoCompute=OFF, maxValue=2e9)
    session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(minAutoCompute=OFF, minValue=1e6)
    book.save('E_elem_{0}.xls'.format(i_run+1))# -*- coding: mbcs -*-
