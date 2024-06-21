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
#from visualization import *
import visualization
#from connectorBehavior import *
import numpy as np
import argparse
#import visualization
from abaqus import session
from odbAccess import *
#session.viewports
import inspect
import re
import os


def extract_frames(input_file, output_file, field, viewVec, camVec, fram=None, refine_seq=None):
    '''Extracts 3D field output from a completed Abaqus simulation to save as 2D image

    :param str input_file: Relative or absolute path to the output database (odb)
        to operate on.
    :param str output_file: Relative or absolute path the output png file.
    :param str field: Field output to collect (e.g. 'S' for stress)
    :param tuple viewVec: tuple of vector components that define the camera view.
    :param tuple camVec: tuple of vector components that define the vertical view
        orientation.
    :param int fram: Simulation frame number to extract field output. Final frame will be plotted if nothing is specified.
    :param tuple refine_seq: additional options required by some fields to
        further define what field quantity is extracted.

    :returns: ``output_file``
    '''

    odb = input_file
    out_name = output_file

    print('input file = {}'.format(input_file))
    print('output file = {}'.format(output_file))

    # initialize odb session
    myViewport = session.Viewport(name='output', width=300, height=200, border=OFF)
    session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
    myOdb = visualization.openOdb(path=odb, readOnly=TRUE)
    myViewport.setValues(displayedObject=myOdb)

    # set frame and define field
    if fram:
        frame = myOdb.steps['load'].frames[fram]
    else:
        frame = myOdb.steps['load'].frames[-1]
    plotField = frame.fieldOutputs[field]

    # set viewpoint
    myViewport.view.setViewpoint(viewVector=viewVec, cameraUpVector=camVec)

    # set background color --> doesn't work yet :(
    #session.graphicsOptions.setValues(backgroundStyle=SOLID, backgroundColor='White')
    #session.graphicsOptions.setValues()

    # set annotation options
    myViewport.viewportAnnotationOptions.setValues(
        triad=ON, compass=OFF, title=OFF, state=OFF,
        legendBackgroundStyle=OTHER, legendBackgroundColor='White')

    # set primary variable
    if refine_seq:
        myViewport.odbDisplay.setPrimaryVariable(
            field=plotField, outputPosition=INTEGRATION_POINT, refinement=refine_seq)
    else:
        myViewport.odbDisplay.setPrimaryVariable(
            field=plotField, outputPosition=INTEGRATION_POINT)

    # display and output
    myViewport.odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF,))

    session.printToFile(fileName=out_name, format=PNG, canvasObjects=(myViewport, ))

    # close odb
    myOdb.close()

    return 0


def get_parser():

    # The global '__file__' variable doesn't appear to be set when executing from Abaqus CAE
    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)

    prog = "abaqus cae -noGui {} -- ".format(basename)
    cli_description = "Extracts 3D field output from a completed Abaqus simulation to save as 2D image"

    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=prog)
    parser.add_argument('-i', '--input-file', type=str, required=True,
                        help="The Abaqus input file created by ``build_model.py``. ")
    parser.add_argument('-o', '--output-file', type=str, required=True,
                        help="The modified Abaqus input file")
    parser.add_argument('--frame', type=int, required=False, default=None,
                        help="Simulation frame number to extract field output. Final frame will be plotted if nothing is specified.")
    parser.add_argument('--field', type=str, required=True,
                        help='Field to extract')

    return parser


if __name__ == '__main__':

    parser = get_parser()

    # set view options
    viewVec = (5,4,2.5)
    camVec = (0,0,1)

    # set field and refinement seqeuence
    #field = 'S'
    refine_seq = (INVARIANT, 'Mises')

    args, unknown = parser.parse_known_args()
    sys.exit(extract_frames(input_file=args.input_file,
                            output_file=args.output_file,
                            field=args.field,
                            viewVec=viewVec,
                            camVec=camVec,
                            fram=args.frame,
                            refine_seq=refine_seq,
                            ))
