import sys
import os
import argparse
import inspect
import numpy as np

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
import inspect


def main(model_name, diam, height, seed, material_E, material_nu, material_rho, total_force, duration, num_steps, fix_lateral_dofs=False, finite_rise=None):
    '''Creates, partitions, and meshes an Abaqus model. Material
    properties, boundary conditions, and loads are applied. The complete
    Abaqus job is written as an input file (.inp).

    :param str model_name: The name of the Abaqus model
    :param float diam: The diameter of the cylinder
    :param float height: The height of the cylinder
    :param float seed:  The approximate global seed size for meshing
    :param float material_E: The elastic modulus of the material
    :param float material_nu: The Poisson ratio of the material
    :param float material_rho: The density (g/cm^3) of the material. This value will be multiplied by 1.00e-9 to convert to units of tonne/mm^3
    :param float total_force: The force applied to cylinder
    :param float duration: The duration of the simulation
    :param int num_steps: The number of fixed time increments
    :param bool fix_lateral_dofs: Option to force all x- and y-displacements to be fixed
    :param int finite_rise: Optional extra number of time steps over which to to ramp force

    :returns: write ``model_name``.cae and ``model_name``.inp

    **Node sets:**

    * ``ALL NODES`` - all nodes of the meshed domain.
    * ``BOTTOM`` - nodes of the bottom z-face of the rectangular domain.
    * ``LOAD_HERE`` - fictious node kinematically coupled to the TOP nodeset used for load application and force summation.
    * ``TOP`` - nodes of the top z-face of the rectangular domain.
    * ``X-PLANE`` - nodes of central plane of rectangular domain with normal in the x-direction.
    * ``Y-PLANE`` - nodes of central plane of rectangular domain with normal in the y-direction.
    * ``material_SET`` - nodes associated with cylindrical region.
    * ``VOID_SET`` - nodes associated with outer void region.

    **Element sets:**

    * ``ALL ELEMENTS`` - all elements of the meshed domain.
    * ``BOTTOM`` - elements of the bottom z-face of the rectangular domain.
    * ``TOP`` - elements of the top z-face of the rectangular domain.
    * ``X-PLANE`` - elements of central plane of rectangular domain with normal in the x-direction.
    * ``Y-PLANE`` - elements of central plane of rectangular domain with normal in the y-direction.
    * ``material_SET`` - elements associated with cylindrical region.
    * ``VOID_SET`` - elements associated with outer void region.
    * ``SET_1`` - elements in first octant of rectangular domain. Micro-domain set #1.
    * ``SET_2`` - elements in second octant of rectangular domain. Micro-domain set #2.
    * ``SET_3`` - elements in third octant of rectangular domain. Micro-domain set #3.
    * ``SET_4`` - elements in fourth octant of rectangular domain. Micro-domain set #4.
    * ``SET_5`` - elements in fifth octant of rectangular domain. Micro-domain set #5.
    * ``SET_6`` - elements in sixth octant of rectangular domain. Micro-domain set #6.
    * ``SET_7`` - elements in seventh octant of rectangular domain. Micro-domain set #7.
    * ``SET_8`` - elements in eighth octant of rectangular domain. Micro-domain set #8.
    '''

    # === Parameters ===
    modelName_1 = model_name
    material_rho  = material_rho*(1.0e-9)

    tol = 0.001

    # === ASSERTIONS ===
    assert diam > 0.0
    assert height > 0.0
    assert seed > 0.0
    assert material_E > 0.0
    assert material_rho > 0.0
    assert num_steps > 0

    # === Create model ===
    modelName = modelName_1

    viewportName = session.Viewport(name=modelName)
    viewportName.makeCurrent()
    viewportName.maximize()

    Model = mdb.Model(name=modelName)

    # === Create Initial Parts Before Boolean Merge ===
    ## Cylinder
    rad = diam/2
    sketch = Model.ConstrainedSketch(name='__profile__', sheetSize=200.0)
    sketch.CircleByCenterPerimeter(
        center=(0.0, 0.0), point1=(0.0, rad))
    cylinder = Model.Part(dimensionality=THREE_D, name='cylinder', type=DEFORMABLE_BODY)
    cylinder.BaseSolidExtrude(depth=height, sketch=sketch)
    del sketch

    # === Assemble and Merge ===
    assembly = Model.rootAssembly
    assembly.DatumCsysByDefault(CARTESIAN)

    ## Create instances for cube and cylinder
    assembly.Instance(dependent=ON, name='cylinder-1',
        part=Model.parts['cylinder'])

    part = cylinder

    ## define top and bottom sets and top surface
    assembly.Set(name='top',
        faces=assembly.instances['cylinder-1'].faces.findAt(
            ((tol, tol, height),), ))
    assembly.Set(name='bottom',
        faces=assembly.instances['cylinder-1'].faces.findAt(
            ((tol, tol, 0.),), ))
    part.Surface(name='top',
        side1Faces=part.faces.findAt(
            ((tol, tol, height),),))

    # === MATERIALS and Sections ===
    material = Model.Material(name='material')
    material.Elastic(table=((material_E, material_nu), ))
    material.Density(table=((material_rho, ), ))

    Model.HomogeneousSolidSection(material='material', name='material',
        thickness=None)

    part.Set(name='all',
        cells=part.cells.findAt(((0.0, 0.0, tol),),) )
    part.SectionAssignment(offset=0.0,
        offsetField='', offsetType=MIDDLE_SURFACE,
        region=part.sets['all'],
        sectionName='material',
        thicknessAssignment=FROM_SECTION)

    # === MESH ===

    ## Partition
    W = diam
    H2 = height/2
    width = diam/2
    part.PartitionCellByPlaneThreePoints(
        cells=part.cells.getByBoundingBox(xMin=-W,xMax=W,yMin=-W,yMax=W,zMin=-W,zMax=W),
        point1=part.InterestingPoint(part.edges[0],MIDDLE),
        point2=part.vertices[0],
        point3=part.vertices[1])

    part.PartitionCellByPlanePointNormal(
        cells=part.cells.getByBoundingBox(xMin=-W,xMax=W,yMin=-W,yMax=W,zMin=-W,zMax=W),
        normal=part.edges[0],
        point=part.InterestingPoint(part.edges[4], CENTER))

    part.PartitionCellByPlanePointNormal(
        cells=part.cells.getByBoundingBox(xMin=-W,xMax=W,yMin=-W,yMax=W,zMin=-W,zMax=W),
        normal=part.edges[5],
        point=part.InterestingPoint(part.edges[5], MIDDLE))

    ## Mesh
    part.setElementType(elemTypes=(
        ElemType(elemCode=C3D8, elemLibrary=STANDARD, secondOrderAccuracy=OFF,
        distortionControl=DEFAULT), ElemType(elemCode=C3D6, elemLibrary=STANDARD),
        ElemType(elemCode=C3D4, elemLibrary=STANDARD)),
        regions=(part.cells.getByBoundingBox(xMin=-W,xMax=W,yMin=-W,yMax=W,zMin=-W,zMax=W),))

    part.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=seed)
    part.generateMesh()

    # === MICRODOMAINS ===

    nodes = get_nodes(part)
    centroids = get_centroids(part,nodes)
    A = -width
    B = 0.
    C = width
    uDomain_1 = get_sets(part,B,C,A,B,0.,H2, centroids, 'Set_1')
    uDomain_2 = get_sets(part,B,C,B,C,0.,H2, centroids, 'Set_2')
    uDomain_3 = get_sets(part,A,B,B,C,0.,H2, centroids, 'Set_3')
    uDomain_4 = get_sets(part,A,B,A,B,0.,H2, centroids, 'Set_4')
    uDomain_5 = get_sets(part,B,C,A,B,H2,2*H2, centroids, 'Set_5')
    uDomain_6 = get_sets(part,B,C,B,C,H2,2*H2, centroids, 'Set_6')
    uDomain_7 = get_sets(part,A,B,B,C,H2,2*H2, centroids, 'Set_7')
    uDomain_8 = get_sets(part,A,B,A,B,H2,2*H2, centroids, 'Set_8')

    all = get_sets(part,A,C,A,C,0.,height, centroids, 'ALL')

    ## Get mid planes for symmetry!
    assembly.Set(name='y-plane',
        faces=assembly.instances['cylinder-1'].faces.findAt(
            ((A+tol, 0.0, tol),),
            ((C-tol, 0.0, tol),),
            ((-tol, 0.0, tol),),
            ((tol, 0.0, tol),),
            ((A+tol, 0.0, H2+tol),),
            ((C-tol, 0.0, H2+tol),),
            ((-tol, 0.0, H2+tol),),
            ((tol, 0.0, H2+tol), ),) )
    assembly.Set(name='x-plane',
        faces=assembly.instances['cylinder-1'].faces.findAt(
            ((0.0, A+tol, tol),),
            ((0.0, C-tol, tol),),
            ((0.0, -tol, tol),),
            ((0.0, tol, tol),),
            ((0.0, A+tol, H2+tol),),
            ((0.0, C-tol, H2+tol),),
            ((0.0, -tol, H2+tol),),
            ((0.0, tol, H2+tol), ),) )

    # === STEP, LOADS, BCs, ETC. ===

    # calculate number of time steps and duration if modified
    time_inc = duration/num_steps
    if finite_rise:
        duration = duration + finite_rise*time_inc

    Model.ImplicitDynamicsStep(
        name='load',
        previous='Initial',
        timePeriod=duration,
        initialInc=time_inc,
        nohaf=OFF,
        timeIncrementationMethod=FIXED,
        maxNumInc=1000000,
        nlgeom=ON,
        )
    ## Turn off algorithmic damping
    Model.steps['load'].setValues(alpha=0.0)

    ## Fix
    ### Bottom
    Model.DisplacementBC(
        amplitude=UNSET, createStepName='Initial',
        distributionType=UNIFORM, fieldName='', localCsys=None,
        name='fix',
        region=Model.rootAssembly.sets['bottom'],
        u1=UNSET, u2=UNSET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    if fix_lateral_dofs:
        ### All x- and y-DOFs fixed!
        Model.DisplacementBC(amplitude=UNSET, createStepName='Initial',
            distributionType=UNIFORM, fieldName='', localCsys=None,
            name='fix_all', region=assembly.sets['cylinder-1.all'],
            u1=SET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    else:
        ### x-plane
        Model.DisplacementBC(amplitude=UNSET, createStepName='Initial',
            distributionType=UNIFORM, fieldName='', localCsys=None,
            name='fix_x-plane', region=assembly.sets['x-plane'],
            u1=SET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
        ### y-plane
        Model.DisplacementBC(amplitude=UNSET, createStepName='Initial',
            distributionType=UNIFORM, fieldName='', localCsys=None,
            name='fix_y-plane', region=assembly.sets['y-plane'],
            u1=UNSET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

    ## RP
    assembly.ReferencePoint(point=(0.0, 0.0, 1.1*height))
    assembly.Set(name='load_here', referencePoints=(assembly.referencePoints.values()[-1],))
    Model.Equation(name='Constraint-1',
        terms=((1.0, 'top', 3), (-1.0, 'load_here', 3)))

    # Prescribed Load
    Model.Pressure(
        amplitude=UNSET,
        createStepName='load',
        distributionType=TOTAL_FORCE,
        magnitude=total_force,
        name='P0',
        region=Model.rootAssembly.instances['cylinder-1'].surfaces['top'])

    if finite_rise:
        Model.TabularAmplitude(name='Amp-1',
            smooth=SOLVER_DEFAULT,
            timeSpan=STEP,
            data=((0.0, 0.0),
                  (finite_rise*time_inc, 1.0),
                  (duration, 1.0)),)
        Model.loads['P0'].setValues(amplitude='Amp-1')

    # output requests
    '''
        Make sure that you edit the KEYWORDS so that only 'U' is sampled at the nodes!
        Make extra sure that 'COORD' is sampled at the element level!!!!!
    '''
    Model.fieldOutputRequests['F-Output-1'].setValues(
        timeInterval=time_inc,
        variables=('S', 'U', 'A', 'V', 'EVOL', 'IVOL', 'COORD'))

    # Model.FieldOutputRequest(createStepName='load',
        # name='F-Output-2',
        # timeInterval=(duration/num_steps),
        # variables=('S', 'U', 'A', 'V', 'EVOL', 'IVOL', 'COORD'))

    Model.HistoryOutputRequest(createStepName='load',
        name='LOAD_HERE', rebar=EXCLUDE,
        region=assembly.sets['load_here'],
        sectionPoints=DEFAULT,
        variables=('U3', 'RT'),
        timeInterval=time_inc)

    mdb.saveAs(pathName='{}.cae'.format(modelName))

    # === Create Job ===
    #jobName = modelName + '_job'
    jobName = modelName
    mdb.Job(model=modelName, name=modelName)
    mdb.jobs[jobName].writeInput()

    return 0


def get_nodes(part):
    """
    Utility function to collect all nodes and nodal coordinates for a given Abaqus part.

    :param object part: the Abaqus model database (mdb) part to operate on.

    :returns: dictionary of nodes and nodal coordinates.
    """
    dict = {}
    nodes = part.nodes[:]
    for node in nodes:
        #print(node.label)
        x = node.coordinates[0]
        y = node.coordinates[1]
        z = node.coordinates[2]
        dict[node.label]=[x,y,z]
    return(dict)


def get_centroids(part, nodes):
    """
    Utility function to calculate centroids of all elements for a given Abaqus part.

    :param object part: the Abaqus model database (mdb) part to operate on.
    :param dict nodes: Dictionary of nodes and nodal coordinates.

    :returns: dictionary of elements and element centroid coordinates.
    """

    dict = {}
    elements=part.elements[:]
    for e in elements:
        #print(e.label)
        #print(e.connectivity)
        label = e.label
        c_x, c_y, c_z = [],[],[]
        for c in e.connectivity:
            c = c+1
            #print(c)
            c_x.append(nodes[c][0])
            c_y.append(nodes[c][1])
            c_z.append(nodes[c][2])
        g_x = np.mean(c_x)
        g_y = np.mean(c_y)
        g_z = np.mean(c_z)
        dict[label]=[g_x,g_y,g_z]

    return(dict)


def get_sets(part, xMin, xMax, yMin, yMax, zMin, zMax, centroids, set_name):
    """
    Utility function to generate element set for all elements with centroids
        contained within a minimum and maximum range of x-, y-, and z-coordinates.

    :param object part: the Abaqus model database (mdb) part to operate on.
    :param float xMin: minimum coordinate in x-direction
    :param float xMax: maximum coordinate in x-direction
    :param float yMin: minimum coordinate in y-direction
    :param float yMax: maximum coordinate in y-direction
    :param float zMin: minimum coordinate in z-direction
    :param float zMax: maximum coordinate in z-direction
    :param dict centroids: dictionary of element IDs and associated coordinates of centroid
    :param str set_name: name of output element set

    :returns: tuple of element labels to include in new set. Part level element
        set is generated.
    """

    set_elems=[]
    elements = part.elements[:]
    for e in elements:
        flag = True
        g_x, g_y, g_z = centroids[e.label]
        if (g_x < xMin) or (g_x > xMax):
            flag = False
        if (g_y < yMin) or (g_y > yMax):
            flag = False
        if (g_z < zMin) or (g_z > zMax):
            flag = False

        # final
        if flag == True:
            set_elems.append(e.label)

    set_elem = tuple(set_elems)
    part.SetFromElementLabels(name=set_name, elementLabels=set_elem)
    return(set_elem)


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)

    # set description and program
    prog = "abaqus cae -noGui {} -- ".format(basename)
    cli_description = "Create an Abaqus model of an elastic cylinder under dynamic compression"

    # add parser arguments
    parser = argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('--model-name', type=str, required=True,
        help='Specify the name of the model')
    parser.add_argument('--diam', type=float, required=True,
        help='Specify the diameter (mm) of the cylinder')
    parser.add_argument('--height', type=float, required=True,
        help='Specify the height (mm) of the cylinder')
    parser.add_argument('--seed', type=float, required=True,
        help='Specify the approximate global seed size (mm) for meshing')
    parser.add_argument('--material-E', type=float, required=True,
        help='Specify the elastic modulus (MPa) of the material')
    parser.add_argument('--material-nu', type=float, required=True,
        help='Specify the Poisson ratio of the material')
    parser.add_argument('--material-rho', type=float, required=True,
        help='Specify the density (g/cm^3) of the material. This value will be multiplied by 1.00e-9 to convert to units of tonne/mm^3')
    parser.add_argument('--total-force', type=float, required=True,
        help='Specify the force applied to cylinder.')
    parser.add_argument('--duration', type=float, required=True,
        help='Specify the duration of the simulation.')
    parser.add_argument('--num-steps', type=int, required=True,
        help='Specify the number of fixed time increments.')
    parser.add_argument('--fix-lateral-dofs', action='store_true', required=False,
        help='Option to force all x- and y-displacements to be fixed')
    parser.add_argument('--finite-rise', type=int, required=False, default=None,
        help='Optional extra number of time steps over which to to ramp force')

    return parser


if __name__ == '__main__':
    parser = get_parser()
    # Abaqus does not strip the CAE options, so we have to skip the unknown options related to the CAE CLI.
    args, unknown = parser.parse_known_args()
    sys.exit(main(model_name=args.model_name,
                  diam=args.diam,
                  height=args.height,
                  seed=args.seed,
                  material_E=args.material_E,
                  material_nu=args.material_nu,
                  material_rho=args.material_rho,
                  total_force=args.total_force,
                  duration=args.duration,
                  num_steps=args.num_steps,
                  fix_lateral_dofs=args.fix_lateral_dofs,
                  finite_rise=args.finite_rise,
                  ))
