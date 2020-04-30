# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 6.13-1 replay file
# Internal Version: 2013_05_16-10.28.56 126354
# Run by WXD on Fri Apr 26 15:20:25 2019
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=246.69841003418, 
    height=185.862945556641)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
mdb.ModelFromInputFile(name='3Dmesh_Solid', 
    inputFileName='C:\Users\WXD\Desktop\ABHmodeling\Model/3Dmesh_Solid.inp')
#: The model "3Dmesh_Solid" has been created.
#: The part "PART-1" has been imported from the input file.
#: The model "3Dmesh_Solid" has been imported from an input file. 
#: Please scroll up to check for error and warning messages.
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
a = mdb.models['3Dmesh_Solid'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
p = mdb.models['3Dmesh_Solid'].parts['PART-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
session.viewports['Viewport: 1'].view.setValues(nearPlane=713.672, 
    farPlane=1215.24, width=488.143, height=348.267, viewOffsetX=-32.1857, 
    viewOffsetY=-28.6379)
elemType1 = mesh.ElemType(elemCode=C3D20, elemLibrary=STANDARD)
p = mdb.models['3Dmesh_Solid'].parts['PART-1']
z1 = p.elements
#elems1 = z1[0:-1] #Select all cells equally
pickedRegions =(z1, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, ))
a = mdb.models['3Dmesh_Solid'].rootAssembly
a.regenerate()
a = mdb.models['3Dmesh_Solid'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
mdb.Job(name='3Dmesh_Solid', model='3Dmesh_Solid', description='', 
    type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
    memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
mdb.jobs['3Dmesh_Solid'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "3Dmesh_Solid.inp".