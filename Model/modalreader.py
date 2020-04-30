#-*-coding:UTF-8-*-
from abaqus import *
from abaqusConstants import *
from visualization import *
# from driverUtils import executeOnCaeStartup
# import math
# import os
from textRepr import *
from odbAccess import *
from abaqusConstants import *
# import string
from jobMessage import *
# import displayGroupOdbToolset as dgo
# import sys
#
# Parameters
# 
FileName='3Dmesh_SC'
StepName='Step-1'
OutputFileName='dampInfor.txt'
#
# Open file
#
odb=openOdb(FileName+'.odb')
StepHere=odb.steps[StepName]
OutputRegion=StepHere.historyRegions['Assembly ASSEMBLY']
f=open(OutputFileName,'w+')

CDarray=OutputRegion.historyOutputs['CD'].data
ModalNumtotal=len(CDarray)
Efarray=OutputRegion.historyOutputs['EIGFREQ'].data
ModalInfor=[]
for modalNum in range(0,ModalNumtotal):
	ModalInfor.append([Efarray[modalNum][1],CDarray[modalNum][1]])
	f.write(str(Efarray[modalNum][1])+' '+str(CDarray[modalNum][1])+'\n')
f.close()
odb.close()
print 'FINISHED !'