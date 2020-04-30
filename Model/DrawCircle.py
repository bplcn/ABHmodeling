# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2016 replay file
# Internal Version: 2015_09_25-04.31.09 126547
# Run by bpl on Today
#
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
# Parameter 
from Parameter import *
Mdb()
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=100.0, height=100.0)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
# Create Plate
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(referenceRepresentation=ON)
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=(-xBound, -yBound), point2=(xBound, yBound))
p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-1']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
#Grid type
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#f ]', ), )
p.setMeshControls(regions=pickedRegions, elemShape=QUAD, algorithm=ADVANCING_FRONT)
#ADVANCING_FRONT， ADVANCING_FRONT
# Partitation the Plate
p = mdb.models['Model-1'].parts['Part-1']
f1, e1, d2 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f1[0], sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',  sheetSize=500.0, gridSpacing=10.0, transform=t)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Part-1']
p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
s.ArcByCenterEnds(center=(gapX, gapY), point1=(gapX+R2, gapY), point2=(gapX+R2, gapY), direction=CLOCKWISE)
s.ArcByCenterEnds(center=(gapX, gapY), point1=(gapX+R1, gapY), point2=(gapX+R1, gapY), direction=CLOCKWISE)
s.ArcByCenterEnds(center=(gapX, gapY), point1=(gapX+Rd, gapY), point2=(gapX+Rd, gapY), direction=CLOCKWISE)


s.ArcByCenterEnds(center=(gapX, -gapY), point1=(gapX+R2, -gapY), point2=(gapX+R2, -gapY), direction=CLOCKWISE)
s.ArcByCenterEnds(center=(gapX, -gapY), point1=(gapX+R1, -gapY), point2=(gapX+R1, -gapY), direction=CLOCKWISE)
s.ArcByCenterEnds(center=(gapX, -gapY), point1=(gapX+Rd, -gapY), point2=(gapX+Rd, -gapY), direction=CLOCKWISE)

s.ArcByCenterEnds(center=(-gapX, -gapY), point1=(-gapX+R2, -gapY), point2=(R2, -gapY), direction=CLOCKWISE)
s.ArcByCenterEnds(center=(-gapX, -gapY), point1=(-gapX+R1, -gapY), point2=(R1, -gapY), direction=CLOCKWISE)
s.ArcByCenterEnds(center=(-gapX, -gapY), point1=(-gapX+Rd, -gapY), point2=(Rd, -gapY), direction=CLOCKWISE)

s.ArcByCenterEnds(center=(-gapX, gapY), point1=(-gapX+R2, gapY), point2=(R2, gapY), direction=CLOCKWISE)
s.ArcByCenterEnds(center=(-gapX, gapY), point1=(-gapX+R1, gapY), point2=(R1, gapY), direction=CLOCKWISE)
s.ArcByCenterEnds(center=(-gapX, gapY), point1=(-gapX+Rd, gapY), point2=(Rd, gapY), direction=CLOCKWISE)


p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedFaces = f.findAt(((0.0,0.0,0.0),),)
e, d1 = p.edges, p.datums
p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
# Mesh the Part
p = mdb.models['Model-1'].parts['Part-1']
p.seedPart(size=lengElem1, deviationFactor=0.1, minSizeFactor=0.1)
#1
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
PickedgeDic=e.getClosest(coordinates=((gapX, gapY+R2, 0.0),),) 
pickedEdges=(PickedgeDic[0][0],) 
p.seedEdgeBySize(edges=pickedEdges, size=lengElem2, deviationFactor=0.1, minSizeFactor=0.1, constraint=FINER)
#2
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
PickedgeDic=e.getClosest(coordinates=((gapX, -gapY+R2, 0.0),),) 
pickedEdges=(PickedgeDic[0][0],) 
p.seedEdgeBySize(edges=pickedEdges, size=lengElem2, deviationFactor=0.1, minSizeFactor=0.1, constraint=FINER)
#3
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
PickedgeDic=e.getClosest(coordinates=((-gapX, -gapY+R2, 0.0),),) 
pickedEdges=(PickedgeDic[0][0],) 
p.seedEdgeBySize(edges=pickedEdges, size=lengElem2, deviationFactor=0.1, minSizeFactor=0.1, constraint=FINER)
#4
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
PickedgeDic=e.getClosest(coordinates=((-gapX, gapY+R2, 0.0),),) 
pickedEdges=(PickedgeDic[0][0],) 
p.seedEdgeBySize(edges=pickedEdges, size=lengElem2, deviationFactor=0.1, minSizeFactor=0.1, constraint=FINER)


p = mdb.models['Model-1'].parts['Part-1']
p.generateMesh()
# Assembly
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Part-1']
a.Instance(name='Part-1-1', part=p, dependent=ON)
# Output the mesh
mdb.Job(name='2Dmesh', model='Model-1') 
mdb.jobs['2Dmesh'].writeInput(consistencyChecking=OFF) 