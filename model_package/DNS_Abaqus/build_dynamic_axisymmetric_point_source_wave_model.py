# -*- coding: mbcs -*-
import argparse
import inspect
import os
import sys

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


def main(model_name, job_name, seed_size=0.2, domain_width=20., domain_height=20.,
         source_depth=1., receiver_radius=5., receiver_depth=5.,
         density=1000., modulus=2.5e6, poisson=0.25,
         duration=0.5, increment=0.0025, load=-100.):
    """Create an Abaqus model approximating an axisymmetric, elastic halfspace subjected a point load

    :param str model_name: The name of the Abaqus model
    :param str job_name: The name of the job (input file) to create
    :param float seed_size: The approximate global seed size for meshing
    :param float domain_width: The width of the domain to simulate (m)
    :param float domain_height: The height of the domain to simulate (m)
    :param float source_depth: The depth on the axis of symmetry to apply the point load (m)
    :param float receiver_radius: The radius from the axis of symmetry to measure quantities of interest (m)
    :param float receiver_depth: The depth to measure quantities of interest (m)
    :param float density: The density of the material (kg/m^3)
    :param float modulus: The elastic modulus of the material (Pa)
    :param float poisson: The Poisson ratio of the material
    :param float duration: The duration to simulate (s)
    :param float increment: The fixed time step (s)
    :param float load: The value of the concentrated point load (N)

    :returns: write ``{model_name}.cae`` and ``{job_name}.inp``

    **Node sets:**

    * ``HALFSPACE-1.ALL`` - all nodes of the meshed domain
    * ``HALFSPACE-1.RADIAL_FIX`` - the nodes on the axis of symmetry
    * ``HALFSPACE-1.RECEIVER`` - the node located at (receiver_radius, receiver_depth)
    * ``HALFSPACE-1.RECEIVER_SURFACE`` - the node located at (receiver_radius, 0)
    * ``HALFSPACE-1.SOURCE`` - the node located at (0, receiver_depth) for applying the point load
    """

    model = mdb.Model(name=model_name)

    # sketch geometry
    model.ConstrainedSketch(name='__profile__', sheetSize=20.0)
    sketch = model.sketches['__profile__']
    sketch.sketchOptions.setValues(viewStyle=AXISYM)
    sketch.ConstructionLine(point1=(0.0, -10.0), point2=(0.0, 10.0))
    sketch.FixedConstraint(entity=sketch.geometry[2])
    sketch.rectangle(point1=(0.0, 0.0), point2=(domain_width, -1*domain_height))

    # make part
    model.Part(dimensionality=AXISYMMETRIC, name='halfspace', type=DEFORMABLE_BODY)
    part = model.parts['halfspace']
    part.BaseShell(sketch=model.sketches['__profile__'])
    del model.sketches['__profile__']

    # partitions
    half_width = domain_width / 2
    half_height = domain_height / 2
    model.ConstrainedSketch(
        gridSpacing=1.41, name='__profile__',
        sheetSize=56.56, transform=
        part.MakeSketchTransform(sketchPlane=part.faces[0],
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT,
        origin=(half_width, -1*half_height, 0.0)))
    part.projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=model.sketches['__profile__'])
    sketch = model.sketches['__profile__']
    if source_depth > 1.e-3:
        sketch.Line(
            point1=(-1*half_width, half_height - source_depth),
            point2=(half_width, half_height - source_depth))
    sketch.HorizontalConstraint(addUndoState=False, entity=sketch.geometry[7])
    sketch.Line(
        point1=(receiver_radius - half_width, half_height),
        point2=(receiver_radius - half_width, -1*half_height))
    sketch.VerticalConstraint(addUndoState=False, entity=sketch.geometry[8])
    sketch.Line(
        point1=(-1*half_width, half_height - receiver_depth),
        point2=(half_width, half_height - receiver_depth))
    sketch.HorizontalConstraint(addUndoState=False, entity=sketch.geometry[9])
    part.PartitionFaceBySketch(
        faces=part.faces.getSequenceFromMask(('[#1 ]', ), ),
        sketch=model.sketches['__profile__'])
    del model.sketches['__profile__']

    # sets
    part.Set(name='source', vertices=part.vertices.getSequenceFromMask(('[#10 ]', ), ))
    part.Set(name='receiver', vertices=part.vertices.getSequenceFromMask(('[#40 ]', ), ))
    part.Set(name='reciever_surface', vertices=part.vertices.getSequenceFromMask(('[#8 ]', ), ))
    part.Set(name='radial_fix', edges=part.edges.getSequenceFromMask(('[#c040 ]', ), ))

    # mesh
    part.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=seed_size)
    part.setMeshControls(
        elemShape=QUAD,
        regions=part.faces.getSequenceFromMask(
        ('[#3f ]', ), ), technique=STRUCTURED)
    part.setElementType(
        elemTypes=(ElemType(elemCode=CAX8, elemLibrary=STANDARD), ElemType(elemCode=CAX6M, elemLibrary=STANDARD)),
        regions=(part.faces.getSequenceFromMask(('[#3f ]', ), ), ))
    part.generateMesh()

    # material and section
    model.Material(name='elastic')
    model.materials['elastic'].Density(table=((density, ), ))
    model.materials['elastic'].Elastic(table=((modulus, poisson), ))
    model.HomogeneousSolidSection(material='elastic', name='Section-1', thickness=None)
    part.Set(faces=part.faces.getSequenceFromMask(('[#3f ]', ), ), name='all')
    part.SectionAssignment(
        offset=0.0,
        offsetField='', offsetType=MIDDLE_SURFACE, region=
        part.sets['all'], sectionName=
        'Section-1', thicknessAssignment=FROM_SECTION)

    # assembly
    model.rootAssembly.DatumCsysByThreePoints(
        coordSysType=CYLINDRICAL, origin=(0.0, 0.0, 0.0), point1=(1.0, 0.0, 0.0), point2=(0.0, 0.0, -1.0))
    model.rootAssembly.Instance(
        dependent=ON, name='halfspace-1', part=model.parts['halfspace'])
    model.rootAssembly.regenerate()

    # BC
    model.DisplacementBC(
        amplitude=UNSET, createStepName='Initial',
        distributionType=UNIFORM, fieldName='', localCsys=None, name='fix_radial',
        region=model.rootAssembly.instances['halfspace-1'].sets['radial_fix'],
        u1=SET, u2=UNSET, ur3=SET)

    # Step
    model.ImplicitDynamicsStep(
        initialInc=increment, maxNumInc=1000,
        name='dynamic_step', noStop=OFF, nohaf=OFF, previous='Initial',
        timeIncrementationMethod=FIXED, timePeriod=duration)
    model.rootAssembly.regenerate()

    # Load
    model.ConcentratedForce(
        cf2=load, createStepName='dynamic_step',
        distributionType=UNIFORM, field='', localCsys=None,
        name='Load-1',
        region=model.rootAssembly.instances['halfspace-1'].sets['source'])

    # History request
    model.HistoryOutputRequest(
        createStepName='dynamic_step', frequency=1, name='H-Output-2', rebar=EXCLUDE, region=
        model.rootAssembly.allInstances['halfspace-1'].sets['receiver'],
        sectionPoints=DEFAULT, variables=('U1', 'U2'))
    model.HistoryOutputRequest(
        createStepName='dynamic_step', frequency=1, name='H-Output-3', rebar=EXCLUDE,
        region=model.rootAssembly.allInstances['halfspace-1'].sets['reciever_surface'],
        sectionPoints=DEFAULT, variables=('U1', 'U2'))

    # Save cae file
    mdb.saveAs(pathName='{}.cae'.format(model_name))

    # Write input file
    mdb.Job(model=model_name, name=job_name)
    mdb.jobs[job_name].writeInput()

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)

    # set description and program
    prog = "abaqus cae -noGui {} -- ".format(basename)
    cli_description = "Create an Abaqus model approximating an axisymmetric, elastic halfspace subjected a point load"

    # add parser arguments
    parser = argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('--model-name', type=str, required=True,
        help='The name of the model')
    parser.add_argument('--job-name', type=str, required=True,
        help='The name of the job (input file) to create')
    parser.add_argument('--seed-size', type=float, required=False, default=0.2,
        help='The approximate global seed size for meshing')
    parser.add_argument('--domain-width', type=float, required=False, default=20.,
        help='The width of the domain to simulate (m)')
    parser.add_argument('--domain-height', type=float, required=False, default=20.,
        help='The height of the domain to simulate (m)')
    parser.add_argument('--source-depth', type=float, required=False, default=1.,
        help='The depth on the axis of symmetry to apply the point load (m)')
    parser.add_argument('--receiver-radius', type=float, required=False, default=5.,
        help='The radius from the axis of symmetry to measure quantities of interest (m)')
    parser.add_argument('--receiver-depth', type=float, required=False, default=5.,
        help='The depth to measure quantities of interest (m)')
    parser.add_argument('--density', type=float, required=False, default=1000.,
        help='The density of the material (kg/m^3)')
    parser.add_argument('--modulus', type=float, required=False, default=2.5e6,
        help='The elastic modulus of the material (Pa)')
    parser.add_argument('--poisson', type=float, required=False, default=0.25,
        help='The Poisson ratio of the material')
    parser.add_argument('--duration', type=float, required=False, default=0.5,
        help='The duration to simulate (s)')
    parser.add_argument('--increment', type=float, required=False, default=0.0025,
        help='The fixed time step (s)')
    parser.add_argument('--load', type=float, required=False, default=-100,
        help='The value of the concentrated point load (N)')

    return parser


if __name__ == '__main__':
    parser = get_parser()
    # Abaqus does not strip the CAE options, so we have to skip the unknown options related to the CAE CLI.
    args, unknown = parser.parse_known_args()
    sys.exit(main(model_name=args.model_name,
                  job_name=args.job_name,
                  seed_size=args.seed_size,
                  domain_width=args.domain_width,
                  source_depth=args.source_depth,
                  receiver_radius=args.receiver_radius,
                  density=args.density,
                  modulus=args.modulus,
                  poisson=args.poisson,
                  duration=args.duration,
                  increment=args.increment,
                  load=args.load,
                  ))
