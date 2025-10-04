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


def get_nodes(part):
    '''Collect all nodes and nodal coordinates for a given Abaqus part.

    :param object part: the Abaqus model database (mdb) part to operate on.

    :returns: dictionary of nodes and nodal coordinates.
    '''

    dict = {}
    nodes = part.nodes[:]
    for node in nodes:
        x = node.coordinates[0]
        y = node.coordinates[1]
        dict[node.label]=[x,y]

    return(dict)


def main(model_name='halfspace_model', job_name='halfpsace_job',
         seed_size=0.2, domain_width=20., domain_height=20.,
         source_depth=1., receiver_radius=5., receiver_depth=5.,
         density=1000., modulus=2.5e6, poisson=0.25,
         duration=0.5, increment=0.0025, load=100., load_direction='y'):
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
    :param float load: The magnitude of the concentrated point load (N)
    :param str load_direction: The direction to apply concetrated load, either 'y' for down or 'x' for right

    :returns: write ``{model_name}.cae`` and ``{job_name}.inp``

    **Node sets:**

    * ``HALFSPACE-1.ALL`` - all nodes of the meshed domain
    * ``HALFSPACE-1.RADIAL_FIX`` - the nodes on the axis of symmetry
    * ``HALFSPACE-1.RECEIVER`` - the node located at (receiver_radius, receiver_depth)
    * ``HALFSPACE-1.RECEIVER_SURFACE`` - the node located at (receiver_radius, 0)
    * ``HALFSPACE-1.SOURCE`` - the node located at (0, receiver_depth) for applying the point load
    """

    model = mdb.Model(name=model_name)
    tol = 1.e-8

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

    # mesh
    part.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=seed_size)
    part.setMeshControls(
        elemShape=QUAD,
        regions=part.faces.getByBoundingBox(xMin=0, xMax=domain_width, yMin=-1*domain_height, yMax=0),
        technique=STRUCTURED)
    part.setElementType(
        elemTypes=(ElemType(elemCode=CAX8, elemLibrary=STANDARD), ElemType(elemCode=CAX6M, elemLibrary=STANDARD)),
        regions=(part.faces.getByBoundingBox(xMin=0, xMax=domain_width, yMin=-1*domain_height, yMax=0), ))
    part.generateMesh()

    # sets
    nodes = get_nodes(part)
    print(nodes)
    source_node, receiver_node, receiver_surface_node, radial_fix_nodes = [], [], [], []
    for n in nodes:
        x, y = nodes[n]
        # source node
        if (x < tol) and (abs(y + source_depth) < tol):
            source_node.append(n)
        # reciever_node
        if (abs(x - receiver_radius) < tol) and (abs(y + receiver_depth) < tol):
            receiver_node.append(n)
        # receiver surface node
        if (abs(x - receiver_radius) < tol) and (abs(y) < tol):
            receiver_surface_node.append(n)
        # radial_fix_nodes
        if x < tol:
            radial_fix_nodes.append(n)
    part.SetFromNodeLabels(name='source', nodeLabels=tuple(source_node))
    part.SetFromNodeLabels(name='receiver', nodeLabels=tuple(receiver_node))
    part.SetFromNodeLabels(name='receiver_surface', nodeLabels=tuple(receiver_surface_node))
    part.SetFromNodeLabels(name='radial_fix', nodeLabels=tuple(radial_fix_nodes))

    # material and section
    model.Material(name='elastic')
    model.materials['elastic'].Density(table=((density, ), ))
    model.materials['elastic'].Elastic(table=((modulus, poisson), ))
    model.HomogeneousSolidSection(material='elastic', name='Section-1', thickness=None)
    part.Set(faces=part.faces.getByBoundingBox(xMin=0, xMax=domain_width, yMin=-1*domain_height, yMax=0),
             name='all')
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
    if load_direction == 'y':
        cf1_load = 0.
        cf2_load = -1*load
    elif load_direction == 'x':
        cf1_load = load
        cf2_load = 0.
    else:
        print('specify valid load direction')
    model.ConcentratedForce(
        cf1=cf1_load, cf2=cf2_load,
        createStepName='dynamic_step',
        distributionType=UNIFORM, field='', localCsys=None,
        name='Load-1',
        region=model.rootAssembly.instances['halfspace-1'].sets['source'])

    # History request
    model.HistoryOutputRequest(
        createStepName='dynamic_step', frequency=1, name='H-Output-2',
        region=model.rootAssembly.allInstances['halfspace-1'].sets['receiver'],
        sectionPoints=DEFAULT,
        variables=('U1', 'U2'))
    model.HistoryOutputRequest(
        createStepName='dynamic_step', frequency=1, name='H-Output-3',
        region=model.rootAssembly.allInstances['halfspace-1'].sets['receiver_surface'],
        sectionPoints=DEFAULT,
        variables=('U1', 'U2'))

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
    parser.add_argument('--model-name', type=str, required=False, default='halfspace_model',
        help='The name of the model')
    parser.add_argument('--job-name', type=str, required=False, default='halfspace_job',
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
    parser.add_argument('--load', type=float, required=False, default=100,
        help='The value of the concentrated point load (N)')
    parser.add_argument('--load-direction', type=str, required=False, default='y',
        help='The direction to apply concetrated load, either "y" for down or "x" for right')
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
                  load_direction=args.load_direction,
                  ))
