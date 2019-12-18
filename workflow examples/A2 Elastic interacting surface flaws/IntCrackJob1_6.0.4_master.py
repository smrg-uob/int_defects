#Abaqus CAE script for interaction between surface and sub-surface defects
#Harry Coules 2016
#Abaqus 6.12
#
## INCLUDE MODULES
from abaqusConstants import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
## CREATE AND PARTITION PARTS
#Plate
#Note Abaqus' maximum allowed value of sheetSize is 1e4.
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=1e4)
mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(#plEdgeB#, 0.0), 
	point2=(#plEdgeA#, #b#))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Plate', type=
	DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Plate'].BaseSolidExtrude(depth=#plDepth#, sketch=
	mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']
#Partition the plate
mdb.models['Model-1'].ConstrainedSketch(gridSpacing=70.71, name='__profile__', 
	sheetSize=1e4, transform=
	mdb.models['Model-1'].parts['Plate'].MakeSketchTransform(
	sketchPlane=mdb.models['Model-1'].parts['Plate'].faces[1], 
	sketchPlaneSide=SIDE1, 
	sketchUpEdge=mdb.models['Model-1'].parts['Plate'].edges[1], 
	sketchOrientation=RIGHT, origin=(0.0, #b#, 0.0)))
mdb.models['Model-1'].parts['Plate'].projectReferencesOntoSketch(filter=
	COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(#plEdgeB#, #meshTransZ#), 
	point2=(#plEdgeA#, #meshTransZ#))
mdb.models['Model-1'].sketches['__profile__'].HorizontalConstraint(
	addUndoState=False, entity=
	mdb.models['Model-1'].sketches['__profile__'].geometry[6])
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(#meshTransSketchB#, 0.0), point2=
	(#meshTransSketchB#, #meshTransZ#))
mdb.models['Model-1'].sketches['__profile__'].VerticalConstraint(addUndoState=
	False, entity=mdb.models['Model-1'].sketches['__profile__'].geometry[7])
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(#meshTransSketchA#, 0.0), point2=(
	#meshTransSketchA#, #meshTransZ#))
mdb.models['Model-1'].sketches['__profile__'].VerticalConstraint(addUndoState=
	False, entity=mdb.models['Model-1'].sketches['__profile__'].geometry[8])
mdb.models['Model-1'].parts['Plate'].PartitionFaceBySketch(faces=
	mdb.models['Model-1'].parts['Plate'].faces.getSequenceFromMask(('[#2 ]', ), 
	), sketch=mdb.models['Model-1'].sketches['__profile__'], sketchUpEdge=
	mdb.models['Model-1'].parts['Plate'].edges[1])
del mdb.models['Model-1'].sketches['__profile__']
mdb.models['Model-1'].parts['Plate'].PartitionCellByExtrudeEdge(cells=
	mdb.models['Model-1'].parts['Plate'].cells.getSequenceFromMask(('[#1 ]', ), 
	), edges=(mdb.models['Model-1'].parts['Plate'].edges[0], 
	mdb.models['Model-1'].parts['Plate'].edges[4], 
	mdb.models['Model-1'].parts['Plate'].edges[8]), line=
	mdb.models['Model-1'].parts['Plate'].edges[18], sense=FORWARD)
mdb.models['Model-1'].parts['Plate'].PartitionCellByExtrudeEdge(cells=
	mdb.models['Model-1'].parts['Plate'].cells.getSequenceFromMask(('[#2 ]', ), 
	), edges=(mdb.models['Model-1'].parts['Plate'].edges[17], ), line=
	mdb.models['Model-1'].parts['Plate'].edges[24], sense=FORWARD)
mdb.models['Model-1'].parts['Plate'].PartitionCellByExtrudeEdge(cells=
	mdb.models['Model-1'].parts['Plate'].cells.getSequenceFromMask(('[#1 ]', ), 
	), edges=(mdb.models['Model-1'].parts['Plate'].edges[25], ), line=
	mdb.models['Model-1'].parts['Plate'].edges[6], sense=FORWARD)
#TubeShellSmallA
#Sweep sketch
mdb.models['Model-1'].ConstrainedSketch(name='__sweep__', sheetSize=1e4)
mdb.models['Model-1'].sketches['__sweep__'].EllipseByCenterPerimeter(
	axisPoint1=(#aA2#, 0.0), axisPoint2=(0.0, #aA1#), center=(0.0, 0.0))
#Auto-trim the Crack A sweep sketch if necessary.
	#autoTrimA#
#Profile sketch
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=#b#, 
	transform=(1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.0, -1.0, -0.0, 1.0, 0.0, 0.0))
mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
	0.0, 0.0), point1=(#rpA1#, 0.0))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='TubeShellSmallA', 
	type=DEFORMABLE_BODY)
mdb.models['Model-1'].parts['TubeShellSmallA'].BaseShellSweep(path=
	mdb.models['Model-1'].sketches['__sweep__'], sketch=
	mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']
del mdb.models['Model-1'].sketches['__sweep__']
#TubeShellLargeA
mdb.models['Model-1'].Part(name='TubeShellLargeA', objectToCopy=
	mdb.models['Model-1'].parts['TubeShellSmallA'])
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=
	mdb.models['Model-1'].parts['TubeShellLargeA'].features['Shell sweep-1'].sketch)
mdb.models['Model-1'].parts['TubeShellLargeA'].projectReferencesOntoSketch(
	filter=COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__edit__'], 
	upToFeature=
	mdb.models['Model-1'].parts['TubeShellLargeA'].features['Shell sweep-1'])
mdb.models['Model-1'].sketches['__edit__'].RadialDimension(curve=
	mdb.models['Model-1'].sketches['__edit__'].geometry[2], radius=#rpA2#, 
	textPoint=(0.0292775630950928, 0.0293945372104645))
mdb.models['Model-1'].parts['TubeShellLargeA'].features['Shell sweep-1'].setValues(
	sketch=mdb.models['Model-1'].sketches['__edit__'])
del mdb.models['Model-1'].sketches['__edit__']
mdb.models['Model-1'].parts['TubeShellLargeA'].regenerate()
#TubeShellSmallB
#Sweep sketch
mdb.models['Model-1'].ConstrainedSketch(name='__sweep__', sheetSize=1e4)
mdb.models['Model-1'].sketches['__sweep__'].EllipseByCenterPerimeter(
	axisPoint1=(#aB2#, 0.0), axisPoint2=(0.0, #aB1#), center=(0.0, 0.0))
#Auto-trim the Crack B sweep sketch if necessary.
	#autoTrimB#
#Profile sketch
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=#b#, 
	transform=(1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.0, -1.0, -0.0, 0.5, 0.0, 0.0))
mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
	0.0, 0.0), point1=(#rpB1#, 0.0))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='TubeShellSmallB', 
	type=DEFORMABLE_BODY)
mdb.models['Model-1'].parts['TubeShellSmallB'].BaseShellSweep(path=
	mdb.models['Model-1'].sketches['__sweep__'], sketch=
	mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']
del mdb.models['Model-1'].sketches['__sweep__']
#TubeShellLargeB
mdb.models['Model-1'].Part(name='TubeShellLargeB', objectToCopy=
	mdb.models['Model-1'].parts['TubeShellSmallB'])
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=
	mdb.models['Model-1'].parts['TubeShellLargeB'].features['Shell sweep-1'].sketch)
mdb.models['Model-1'].parts['TubeShellLargeB'].projectReferencesOntoSketch(
	filter=COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__edit__'], 
	upToFeature=
	mdb.models['Model-1'].parts['TubeShellLargeB'].features['Shell sweep-1'])
mdb.models['Model-1'].sketches['__edit__'].RadialDimension(curve=
	mdb.models['Model-1'].sketches['__edit__'].geometry[2], radius=#rpB2#, 
	textPoint=(0.050321102142334, 0.056938249617815))
mdb.models['Model-1'].parts['TubeShellLargeB'].features['Shell sweep-1'].setValues(
	sketch=mdb.models['Model-1'].sketches['__edit__'])
del mdb.models['Model-1'].sketches['__edit__']
mdb.models['Model-1'].parts['TubeShellLargeB'].regenerate()
# Create EdgeShellA
mdb.models['Model-1'].Part(name='EdgeShellA', objectToCopy=
	mdb.models['Model-1'].parts['TubeShellLargeA'])
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=
	mdb.models['Model-1'].parts['EdgeShellA'].features['Shell sweep-1'].sketch)
mdb.models['Model-1'].parts['EdgeShellA'].projectReferencesOntoSketch(filter=
	COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__edit__'], 
	upToFeature=
	mdb.models['Model-1'].parts['EdgeShellA'].features['Shell sweep-1'])
mdb.models['Model-1'].sketches['__edit__'].delete(objectList=(
	mdb.models['Model-1'].sketches['__edit__'].geometry[2], ))
mdb.models['Model-1'].sketches['__edit__'].Line(point1=(0.0, 0.0), point2=(0.0, 
	#rpA2#))
mdb.models['Model-1'].sketches['__edit__'].VerticalConstraint(addUndoState=
	False, entity=mdb.models['Model-1'].sketches['__edit__'].geometry[3])
mdb.models['Model-1'].parts['EdgeShellA'].features['Shell sweep-1'].setValues(
	sketch=mdb.models['Model-1'].sketches['__edit__'])
del mdb.models['Model-1'].sketches['__edit__']
mdb.models['Model-1'].parts['EdgeShellA'].regenerate()
# Create EdgeShellB
mdb.models['Model-1'].Part(name='EdgeShellB', objectToCopy=
	mdb.models['Model-1'].parts['TubeShellLargeB'])
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=
	mdb.models['Model-1'].parts['EdgeShellB'].features['Shell sweep-1'].sketch)
mdb.models['Model-1'].parts['EdgeShellB'].projectReferencesOntoSketch(filter=
	COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__edit__'], 
	upToFeature=
	mdb.models['Model-1'].parts['EdgeShellB'].features['Shell sweep-1'])
mdb.models['Model-1'].sketches['__edit__'].delete(objectList=(
	mdb.models['Model-1'].sketches['__edit__'].geometry[2], ))
mdb.models['Model-1'].sketches['__edit__'].Line(point1=(0.0, 0.0), point2=(0.0, 
	#rpB2#))
mdb.models['Model-1'].sketches['__edit__'].VerticalConstraint(addUndoState=
	False, entity=mdb.models['Model-1'].sketches['__edit__'].geometry[3])
mdb.models['Model-1'].parts['EdgeShellB'].features['Shell sweep-1'].setValues(
	sketch=mdb.models['Model-1'].sketches['__edit__'])
del mdb.models['Model-1'].sketches['__edit__']
mdb.models['Model-1'].parts['EdgeShellB'].regenerate()
#Create assembly and plate instance
mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Plate-1', part=
	mdb.models['Model-1'].parts['Plate'])
#Create Ellipse A instances in the assembly
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=
	'TubeShellSmallA-1', part=mdb.models['Model-1'].parts['TubeShellSmallA'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=
	'TubeShellLargeA-1', part=mdb.models['Model-1'].parts['TubeShellLargeA'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='EdgeShellA-1', 
	part=mdb.models['Model-1'].parts['EdgeShellA'])
#Translate Ellipse A instances
mdb.models['Model-1'].rootAssembly.translate(instanceList=('TubeShellSmallA-1', 
	), vector=(#xA#, 0.0, 0.0))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('TubeShellLargeA-1', 
	), vector=(#xA#, 0.0, 0.0))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('EdgeShellA-1', ), 
	vector=(#xA#, 0.0, 0.0))
#Create Ellipse B instances in the assembly
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=
	'TubeShellSmallB-1', part=mdb.models['Model-1'].parts['TubeShellSmallB'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=
	'TubeShellLargeB-1', part=mdb.models['Model-1'].parts['TubeShellLargeB'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='EdgeShellB-1', 
	part=mdb.models['Model-1'].parts['EdgeShellB'])
#Rotate Ellipse B instances
	#rotateEllipseB#
#Translate Ellipse B instances
mdb.models['Model-1'].rootAssembly.translate(instanceList=('TubeShellSmallB-1', 
	), vector=(#xB#, #yB#, 0.0))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('TubeShellLargeB-1', 
	), vector=(#xB#, #yB#, 0.0))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('EdgeShellB-1', ), 
	vector=(#xB#, #yB#, 0.0))
#Partition using shell instances, create PlatePartitioned
mdb.models['Model-1'].rootAssembly.InstanceFromBooleanCut(cuttingInstances=(
	mdb.models['Model-1'].rootAssembly.instances['EdgeShellA-1'], 
	mdb.models['Model-1'].rootAssembly.instances['TubeShellSmallA-1'], 
	mdb.models['Model-1'].rootAssembly.instances['TubeShellLargeA-1'], 
	mdb.models['Model-1'].rootAssembly.instances['EdgeShellB-1'], 
	mdb.models['Model-1'].rootAssembly.instances['TubeShellSmallB-1'], 
	mdb.models['Model-1'].rootAssembly.instances['TubeShellLargeB-1']), 
	instanceToBeCut=mdb.models['Model-1'].rootAssembly.instances['Plate-1'], 
	name='PlatePartitioned', originalInstances=SUPPRESS)
mdb.models['Model-1'].parts['PlatePartitioned'].setValues(geometryRefinement=
    FINE)
#Further partitioning if divideFlag = true
	#crackZoneDivide#
#Partition Crack B if it is a subsurface crack
	#partitionSubsurfaceB#
#Further partitioning if divideFlag = true
#NOTE - No code has been written here. It was not necessary to divide the crack region.
#Make instance independent
mdb.models['Model-1'].rootAssembly.makeIndependent(instances=(
	mdb.models['Model-1'].rootAssembly.instances['PlatePartitioned-1'], ))
#Add crack features
	#edgeCoords.A.tip#
	#edgeCoords.B.tip#
## NON-GEOMETRIC MODEL FEATURES (material, step, BCs etc.)
#Material properties
	#material#
mdb.models['Model-1'].HomogeneousSolidSection(material='Material-1', name=
	'Section-1', thickness=None)
mdb.models['Model-1'].parts['PlatePartitioned'].Set(cells=
	mdb.models['Model-1'].parts['PlatePartitioned'].cells.getByBoundingBox(#boundBoxMag#)
	, name='Set-1')
mdb.models['Model-1'].parts['PlatePartitioned'].SectionAssignment(offset=0.0, 
	offsetField='', offsetType=MIDDLE_SURFACE, region=
	mdb.models['Model-1'].parts['PlatePartitioned'].sets['Set-1'], sectionName=
	'Section-1', thicknessAssignment=FROM_SECTION)
mdb.models['Model-1'].rootAssembly.regenerate()
#Step
mdb.models['Model-1'].StaticStep(name='ApplyLoad', previous='Initial', nlgeom=#nlgeomFlag#)
	#step#
#Boundary conditions
	#edgeCoords.sets.zSymm#
	#edgeCoords.sets.xSymm#
	#edgeCoords.sets.fixed#
#Coupling constraint and concentrated load
	#edgeCoords.sets.coupling#
#Crack face pressure load
	#edgeCoords.sets.crackFacePressure#
#Pipe pressure load
	#edgeCoords.sets.pipePressure#
#Stress in plane of crack (x dir.)
	#edgeCoords.sets.biaxialStress#
#Create History output requests
	#contourIntegralA#
	#contourIntegralB#
## MESHING
## Crack A
	#edgeCoords.A.radial1#
	#edgeCoords.A.radial2#
	#edgeCoords.A.circ#
	#edgeCoords.A.sweep#
	#edgeCoords.A.mid#
## Crack B
	#edgeCoords.B.radial1#
	#edgeCoords.B.radial2#
	#edgeCoords.B.circ#
	#edgeCoords.B.sweep#
	#edgeCoords.B.mid#
## OTHER GEOMETRY - around the cracks
#Edge along x direction on xy symmetry plane - in between cracks
	#edgeCoords.aroundCracks.inBetweenCracks#
#Edge along x direction on xy symmetry plane - opposite side of plate
	#edgeCoords.aroundCracks.oppositeCracks#
#Edge along x direction on xy symmetry plane - either side of cracks
	#edgeCoords.aroundCracks.eitherSideCracks#
## OTHER GEOMETRY - aligned with x-axis
#Edges along x direction in crack zone at z = 3
	#edgeCoords.otherXdir.crackZone#
#Edges along x direction at z = 0 and z = 3 away from crack zone
	#edgeCoords.otherXdir.outer#
##Edges along x direction at the other side of the plate (z = 1000)
	#edgeCoords.otherXdir.lower#
## OTHER GEOMETRY - aligned with y-axis
#All edges along y direction
	#edgeCoords.otherYdir.all#
## OTHER GEOMETRY - aligned with z-axis
#Edges along the z direction at the same height in z as the crack zone
	#edgeCoords.otherZdir.crackZone#
#Edges along the z direction lower down in z
	#edgeCoords.otherZdir.lower#
## Assign mesh controls
	#edgeCoords.meshControls#
#Generate mesh
mdb.models['Model-1'].rootAssembly.generateMesh(regions=(
	mdb.models['Model-1'].rootAssembly.instances['PlatePartitioned-1'], ))
##WRITE INPUT
#Create job
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
	explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
	memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
	multiprocessingMode=DEFAULT, name='IntCrackJob1', nodalOutputPrecision=
	SINGLE, numCpus=1, numGPUs=0, queue=None, scratch='', type=ANALYSIS, 
	userSubroutine='', waitHours=0, waitMinutes=0)
#Write input
mdb.jobs['IntCrackJob1'].writeInput()