#!python
import argparse
import pathlib
import sys

import numpy

import file_io.xdmf


def random_points_in_cylinder(npts, radius, height, xoffset=0., yoffset=0., zoffset=0.):
    """Generate random points in a cylinder. The bottom of the cylinder is centered on the original

    :params int npts: The number microfield points
    :params float radius: The radius of the cylinder
    :params float height: The height of the cylinder

    :returns: x, y, and z points in the cylindrical domain    
    """

    r = radius * numpy.sqrt(numpy.random.uniform(0, 1, npts)) #sqrt for uniform distribution over area
    theta = numpy.random.uniform(0, 2 * numpy.pi, npts)
    z = numpy.random.uniform(0, height, npts) + zoffset
    
    # Convert to Cartesian coordinates
    x = r * numpy.cos(theta) + xoffset
    y = r * numpy.sin(theta) + yoffset
    
    return numpy.column_stack((x, y, z))


def random_points_in_cube(npts, height):
    """Generate random points in a cube. The bottom of the cube is centered on the original

    :params int npts: The number microfield points
    :params float height: The side length of the cube

    :returns: x, y, and z points in the cylindrical domain    
    """

    x = numpy.random.uniform(-0.5*height, 0.5*height, npts)
    y = numpy.random.uniform(-0.5*height, 0.5*height, npts)
    z = numpy.random.uniform(0, height, npts)
    
    return numpy.column_stack((x, y, z))


def calculate_cylinder_displacement_uniaxial_stress(coords, disp_z, disp_lat, height, rad, times):
    """Calculate the displacement history for a cylinder subjected to uniaxial stress

    :params array-like coords: The x, y, and z coordinates of points
    :param float disp_z: The final displacement to prescribe in the z-direction
    :param float disp_lat: The final displacement to prescribe in the x- and y-directions
    :params float height: The height of the cylinder
    :params float rad: The radius of the cylinder
    :params array-like times: Array of time steps

    :returns: list of displacement arrays for each micro point and time step
    """

    all_disps = []
    for t in times:
        disps = numpy.zeros_like(coords)
        for i, (x,y,z) in enumerate(coords):
            disps[i,0] = t * (x/rad) * disp_lat
            disps[i,1] = t * (y/rad) * disp_lat
            disps[i,2] = t * (z/height) * disp_z
        all_disps.append(disps)

    return all_disps


def calculate_cube_displacement_uniaxial_stress(coords, disp_z, disp_lat, height, times):
    """Calculate the displacement history for a cube subjected to uniaxial stress

    :params array-like coords: The x, y, and z coordinates of points
    :param float disp_z: The final displacement to prescribe in the z-direction
    :param float disp_lat: The final displacement to prescribe in the x- and y-directions
    :params float height: The side length of the cube
    :params array-like times: Array of time steps

    :returns: list of displacement arrays for each micro point and time step
    """

    all_disps = []
    for t in times:
        disps = numpy.zeros_like(coords)
        for i, (x,y,z) in enumerate(coords):
            disps[i,0] = t * (x/height) * disp_lat
            disps[i,1] = t * (y/height) * disp_lat
            disps[i,2] = t * (z/height) * disp_z
        all_disps.append(disps)

    return all_disps


def calculate_uniaxial_stress(npts, final_stress, times):
    """Calculate the uniaxial stress history

    :params int npts: The number microfield points
    :params float final_stress: The final stress to prescribe for the zz component
    :params array-like times: Array of time steps

    :returns: list of stress arrays for each micro point and time step
    """

    all_stresses = []
    for t in times:
        stress = numpy.array([[1., 0., 0., 0., 1., 0., 0., 0., final_stress*t] for _ in range(npts)])
        all_stresses.append(stress)

    return all_stresses


def generate_microfield_file(data_filename, times, npts, coords,
                             all_disps, all_stresses, volumes, densities):
    """Generate the microfield .xdmf and .h5 files for a Micromorphic Filter study

    :params str data_filename: The output filename for the h5 + XDMF file pair (no suffix)
    :params array-like times: Array of time steps
    :params int npts: The number microfield points
    :params array-like coords: The x, y, and z coordinates of points
    :params list all_disps: list of displacement arrays for each micro point and time step
    :params list stresses: list of stress arrays for each micro point and time step
    :params array-like volumes: The volumes associated with each micro point
    :params array-like densities: The densities associated with each micro point

    :returns: ``{output_file}.xdmf`` and ``{output_file}.h5``
    """

    xdmf = file_io.xdmf.XDMF(output_filename=data_filename)

    point_name = 'points'
    conn_name  = 'connectivity'

    for i, t in enumerate(times):
        grid = xdmf.addGrid(xdmf.output_timegrid, {})
        xdmf.addTime(grid, t)
        xdmf.addPoints(grid, coords, duplicate=point_name)
        xdmf.addConnectivity(grid, "POLYVERTEX", numpy.array([v for v in range(npts)]).reshape((-1,1)), duplicate=conn_name)

        xdmf.addData(grid, "disp", all_disps[i], "Node", dtype='d')
        xdmf.addData(grid, "stress", all_stresses[i], "Cell", dtype='d')
        xdmf.addData(grid, "volume", volumes, "Cell", dtype='d')
        xdmf.addData(grid, "density", densities, "Cell", dtype='d')
        print(f'shape of disps = {numpy.shape(all_disps[i])}')

    xdmf.write()
    return 0


def build_microfield_study(output_file, npts, n_steps,
                           disp_z, disp_lat, final_stress,
                           geometry='cylinder', radius=None, height=None,
                           ):
    """Generate a microfield data file of either a cylinder or cube

    :params str output_file: The output filename for the h5 + XDMF file pair (no suffix)
    :params int npts: The number microfield points
    :params int n_steps: The number of timesteps of data to generate
    :params float disp_z: The final displacement to prescribe in the z-direction
    :params float disp_lat: The final displacement to prescribe in the x- and y-directions
    :params float final_stress: The final stress to prescribe for the zz component
    :params str geometry: The geometry type to generate: either "cylinder" or "cube"
    :params float radius: The radius of the cylinder if geometry type is "cylinder"
    :params float height: The height of the cylinder if geometry type is "cylinder" or the side length of the cube if geometry type is "cube
    """

    times = numpy.linspace(0, 1, n_steps)
    import os
    print(f'{os.getcwd()}')
    if geometry == 'cylinder':
        assert radius != None, "radius must be specified!"
        assert height != None, "height must be specified!"
        coords = random_points_in_cylinder(npts, radius, height)
        all_disps = calculate_cylinder_displacement_uniaxial_stress(coords, disp_z, disp_lat, height, radius, times)
        total_volume = numpy.pi * radius * radius * height
    elif geometry == 'cube':
        assert height != None, "height must be specified!"
        coords = random_points_in_cube(npts, height)
        all_disps = calculate_cube_displacement_uniaxial_stress(coords, disp_z, disp_lat, height, times)
        total_volume = height * height * height
    else:
        print('Valid geometry type must be specified!')

    all_stresses = calculate_uniaxial_stress(npts, final_stress, times)

    # evenly split cylindrical volume for each point, keep constant over time
    volumes = numpy.array([(total_volume / npts) for _ in range(npts)]).reshape((-1,1))
    densities = numpy.array([1. for _ in range(npts)]).reshape((-1,1))

    # generate the microfield data file
    generate_microfield_file(output_file, times, npts, coords,
                             all_disps, all_stresses, volumes, densities)

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Generate a microfield data file of either a cylinder or cube"
    parser = argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help='Specify the output filename for the h5 + XDMF file pair (no suffix)')
    parser.add_argument('--npts', type=int, required=True,
        help='The number microfield points')
    parser.add_argument('--n-steps', type=int, required=True,
        help='The number of timesteps of data to generate')
    parser.add_argument('--disp-z', type=float, required=True,
        help='The final displacement to prescribe in the z-direction')
    parser.add_argument('--disp-lat', type=float, required=True,
        help='The final displacement to prescribe in the x- and y-directions')
    parser.add_argument('--final-stress', type=float, required=True,
        help='The final stress to prescribe for the zz component')
    parser.add_argument('--geometry', type=str, required=False, default='cylinder',
        help='The geometry type to generate: either "cylinder" or "cube"')
    parser.add_argument('--radius', type=float, required=False, default=None,
        help='The radius of the cylinder if geometry type is "cylinder"')
    parser.add_argument('--height', type=float, required=False, default=None,
        help='The height of the cylinder if geometry type is "cylinder"\
              or the side length of the cube if geometry type is "cube"')

    return parser


if __name__ == '__main__':
    parser = get_parser()

    args, unknown = parser.parse_known_args()
    sys.exit(build_microfield_study(output_file=args.output_file,
                                    npts=args.npts,
                                    n_steps=args.n_steps,
                                    disp_z=args.disp_z,
                                    disp_lat=args.disp_lat,
                                    final_stress=args.final_stress,
                                    geometry=args.geometry,
                                    radius=args.radius,
                                    height=args.height,
                                    ))