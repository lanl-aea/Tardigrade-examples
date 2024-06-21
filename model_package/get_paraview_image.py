import sys
import os
import argparse
import inspect
import pathlib

import numpy
import pandas
from paraview.simple import *

def set_defaults(display):
    display.Representation = 'Surface'
    display.ColorArrayName = [None, '']
    display.LookupTable = None
    display.MapScalars = 1
    display.MultiComponentsMapping = 0
    display.InterpolateScalarsBeforeMapping = 1
    display.Opacity = 1.0
    display.PointSize = 2.0
    display.LineWidth = 1.0
    display.RenderLinesAsTubes = 0
    display.RenderPointsAsSpheres = 0
    display.Interpolation = 'Gouraud'
    display.Specular = 0.0
    display.SpecularColor = [1.0, 1.0, 1.0]
    display.SpecularPower = 100.0
    display.Luminosity = 0.0
    display.Ambient = 0.0
    display.Diffuse = 1.0
    display.Roughness = 0.3
    display.Metallic = 0.0
    display.Texture = None
    display.RepeatTextures = 1
    display.InterpolateTextures = 0
    display.SeamlessU = 0
    display.SeamlessV = 0
    display.UseMipmapTextures = 0
    display.BaseColorTexture = None
    display.NormalTexture = None
    display.NormalScale = 1.0
    display.MaterialTexture = None
    display.OcclusionStrength = 1.0
    display.EmissiveTexture = None
    display.EmissiveFactor = [1.0, 1.0, 1.0]
    display.FlipTextures = 0
    display.BackfaceRepresentation = 'Follow Frontface'
    display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
    display.BackfaceOpacity = 1.0
    display.Position = [0.0, 0.0, 0.0]
    display.Scale = [1.0, 1.0, 1.0]
    display.Orientation = [0.0, 0.0, 0.0]
    display.Origin = [0.0, 0.0, 0.0]
    display.Pickable = 1
    display.Triangulate = 0
    display.UseShaderReplacements = 0
    display.ShaderReplacements = ''
    display.NonlinearSubdivisionLevel = 1
    display.UseDataPartitions = 0
    display.OSPRayUseScaleArray = 0
    #display.OSPRayScaleArray = 'diagnostic_quantitiesdual.nodal_volume'
    display.OSPRayScaleFunction = 'PiecewiseFunction'
    display.OSPRayMaterial = 'None'
    display.Orient = 0
    display.OrientationMode = 'Direction'
    display.SelectOrientationVectors = 'None'
    display.Scaling = 0
    display.ScaleMode = 'No Data Scaling Off'
    display.ScaleFactor = 600.0
    display.SelectScaleArray = 'None'
    display.GlyphType = 'Arrow'
    display.UseGlyphTable = 0
    display.GlyphTableIndexArray = 'None'
    display.UseCompositeGlyphTable = 0
    display.UseGlyphCullingAndLOD = 0
    display.LODValues = []
    display.ColorByLODIndex = 0
    display.GaussianRadius = 30.0
    display.ShaderPreset = 'Sphere'
    display.CustomTriangleScale = 3
    display.CustomShader = """ // This custom shader code define a gaussian blur
     // Please take a look into vtkSMPointGaussianRepresentation.cxx
     // for other custom shader examples
     //VTK::Color::Impl
       float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
       float gaussian = exp(-0.5*dist2);
       opacity = opacity*gaussian;
    """
    display.Emissive = 0
    display.ScaleByArray = 0
    #display.SetScaleArray = ['POINTS', 'diagnostic_quantitiesdual.nodal_volume']
    #display.ScaleArrayComponent = ''
    display.UseScaleFunction = 1
    display.ScaleTransferFunction = 'PiecewiseFunction'
    display.OpacityByArray = 0
    #display.OpacityArray = ['POINTS', 'diagnostic_quantitiesdual.nodal_volume']
    #display.OpacityArrayComponent = ''
    display.OpacityTransferFunction = 'PiecewiseFunction'
    display.DataAxesGrid = 'GridAxesRepresentation'
    display.SelectionCellLabelBold = 0
    display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
    display.SelectionCellLabelFontFamily = 'Arial'
    display.SelectionCellLabelFontFile = ''
    display.SelectionCellLabelFontSize = 18
    display.SelectionCellLabelItalic = 0
    display.SelectionCellLabelJustification = 'Left'
    display.SelectionCellLabelOpacity = 1.0
    display.SelectionCellLabelShadow = 0
    display.SelectionPointLabelBold = 0
    display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
    display.SelectionPointLabelFontFamily = 'Arial'
    display.SelectionPointLabelFontFile = ''
    display.SelectionPointLabelFontSize = 18
    display.SelectionPointLabelItalic = 0
    display.SelectionPointLabelJustification = 'Left'
    display.SelectionPointLabelOpacity = 1.0
    display.SelectionPointLabelShadow = 0
    display.PolarAxes = 'PolarAxesRepresentation'
    display.ScalarOpacityFunction = None
    display.ScalarOpacityUnitDistance = 67.25568450249996
    #display.ExtractedBlockIndex = 0
    display.SelectMapper = 'Projected tetra'
    display.SamplingDimensions = [128, 128, 128]
    display.UseFloatingPointFrameBuffer = 1

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
    display.OSPRayScaleFunction.UseLogScale = 0

    # init the 'Arrow' selected for 'GlyphType'
    display.GlyphType.TipResolution = 6
    display.GlyphType.TipRadius = 0.1
    display.GlyphType.TipLength = 0.35
    display.GlyphType.ShaftResolution = 6
    display.GlyphType.ShaftRadius = 0.03
    display.GlyphType.Invert = 0

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    display.ScaleTransferFunction.Points = [2535.0086363001424, 0.0, 0.5, 0.0, 7276674.992940602, 1.0, 0.5, 0.0]
    display.ScaleTransferFunction.UseLogScale = 0

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    display.OpacityTransferFunction.Points = [2535.0086363001424, 0.0, 0.5, 0.0, 7276674.992940602, 1.0, 0.5, 0.0]
    display.OpacityTransferFunction.UseLogScale = 0

    # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
    display.DataAxesGrid.XTitle = 'X Axis'
    display.DataAxesGrid.YTitle = 'Y Axis'
    display.DataAxesGrid.ZTitle = 'Z Axis'
    display.DataAxesGrid.XTitleFontFamily = 'Arial'
    display.DataAxesGrid.XTitleFontFile = ''
    display.DataAxesGrid.XTitleBold = 0
    display.DataAxesGrid.XTitleItalic = 0
    display.DataAxesGrid.XTitleFontSize = 12
    display.DataAxesGrid.XTitleShadow = 0
    display.DataAxesGrid.XTitleOpacity = 1.0
    display.DataAxesGrid.YTitleFontFamily = 'Arial'
    display.DataAxesGrid.YTitleFontFile = ''
    display.DataAxesGrid.YTitleBold = 0
    display.DataAxesGrid.YTitleItalic = 0
    display.DataAxesGrid.YTitleFontSize = 12
    display.DataAxesGrid.YTitleShadow = 0
    display.DataAxesGrid.YTitleOpacity = 1.0
    display.DataAxesGrid.ZTitleFontFamily = 'Arial'
    display.DataAxesGrid.ZTitleFontFile = ''
    display.DataAxesGrid.ZTitleBold = 0
    display.DataAxesGrid.ZTitleItalic = 0
    display.DataAxesGrid.ZTitleFontSize = 12
    display.DataAxesGrid.ZTitleShadow = 0
    display.DataAxesGrid.ZTitleOpacity = 1.0
    display.DataAxesGrid.FacesToRender = 63
    display.DataAxesGrid.CullBackface = 0
    display.DataAxesGrid.CullFrontface = 1
    display.DataAxesGrid.ShowGrid = 0
    display.DataAxesGrid.ShowEdges = 1
    display.DataAxesGrid.ShowTicks = 1
    display.DataAxesGrid.LabelUniqueEdgesOnly = 1
    display.DataAxesGrid.AxesToLabel = 63
    display.DataAxesGrid.XLabelFontFamily = 'Arial'
    display.DataAxesGrid.XLabelFontFile = ''
    display.DataAxesGrid.XLabelBold = 0
    display.DataAxesGrid.XLabelItalic = 0
    display.DataAxesGrid.XLabelFontSize = 12
    display.DataAxesGrid.XLabelShadow = 0
    display.DataAxesGrid.XLabelOpacity = 1.0
    display.DataAxesGrid.YLabelFontFamily = 'Arial'
    display.DataAxesGrid.YLabelFontFile = ''
    display.DataAxesGrid.YLabelBold = 0
    display.DataAxesGrid.YLabelItalic = 0
    display.DataAxesGrid.YLabelFontSize = 12
    display.DataAxesGrid.YLabelShadow = 0
    display.DataAxesGrid.YLabelOpacity = 1.0
    display.DataAxesGrid.ZLabelFontFamily = 'Arial'
    display.DataAxesGrid.ZLabelFontFile = ''
    display.DataAxesGrid.ZLabelBold = 0
    display.DataAxesGrid.ZLabelItalic = 0
    display.DataAxesGrid.ZLabelFontSize = 12
    display.DataAxesGrid.ZLabelShadow = 0
    display.DataAxesGrid.ZLabelOpacity = 1.0
    display.DataAxesGrid.XAxisNotation = 'Mixed'
    display.DataAxesGrid.XAxisPrecision = 2
    display.DataAxesGrid.XAxisUseCustomLabels = 0
    display.DataAxesGrid.XAxisLabels = []
    display.DataAxesGrid.YAxisNotation = 'Mixed'
    display.DataAxesGrid.YAxisPrecision = 2
    display.DataAxesGrid.YAxisUseCustomLabels = 0
    display.DataAxesGrid.YAxisLabels = []
    display.DataAxesGrid.ZAxisNotation = 'Mixed'
    display.DataAxesGrid.ZAxisPrecision = 2
    display.DataAxesGrid.ZAxisUseCustomLabels = 0
    display.DataAxesGrid.ZAxisLabels = []
    display.DataAxesGrid.UseCustomBounds = 0
    display.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

    # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
    display.PolarAxes.Visibility = 0
    display.PolarAxes.Translation = [0.0, 0.0, 0.0]
    display.PolarAxes.Scale = [1.0, 1.0, 1.0]
    display.PolarAxes.Orientation = [0.0, 0.0, 0.0]
    display.PolarAxes.EnableCustomBounds = [0, 0, 0]
    display.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    display.PolarAxes.EnableCustomRange = 0
    display.PolarAxes.CustomRange = [0.0, 1.0]
    display.PolarAxes.PolarAxisVisibility = 1
    display.PolarAxes.RadialAxesVisibility = 1
    display.PolarAxes.DrawRadialGridlines = 1
    display.PolarAxes.PolarArcsVisibility = 1
    display.PolarAxes.DrawPolarArcsGridlines = 1
    display.PolarAxes.NumberOfRadialAxes = 0
    display.PolarAxes.AutoSubdividePolarAxis = 1
    display.PolarAxes.NumberOfPolarAxis = 0
    display.PolarAxes.MinimumRadius = 0.0
    display.PolarAxes.MinimumAngle = 0.0
    display.PolarAxes.MaximumAngle = 90.0
    display.PolarAxes.RadialAxesOriginToPolarAxis = 1
    display.PolarAxes.Ratio = 1.0
    display.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
    display.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
    display.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
    display.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
    display.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
    display.PolarAxes.PolarAxisTitleVisibility = 1
    display.PolarAxes.PolarAxisTitle = 'Radial Distance'
    display.PolarAxes.PolarAxisTitleLocation = 'Bottom'
    display.PolarAxes.PolarLabelVisibility = 1
    display.PolarAxes.PolarLabelFormat = '%-#6.3g'
    display.PolarAxes.PolarLabelExponentLocation = 'Labels'
    display.PolarAxes.RadialLabelVisibility = 1
    display.PolarAxes.RadialLabelFormat = '%-#3.1f'
    display.PolarAxes.RadialLabelLocation = 'Bottom'
    display.PolarAxes.RadialUnitsVisibility = 1
    display.PolarAxes.ScreenSize = 10.0
    display.PolarAxes.PolarAxisTitleOpacity = 1.0
    display.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
    display.PolarAxes.PolarAxisTitleFontFile = ''
    display.PolarAxes.PolarAxisTitleBold = 0
    display.PolarAxes.PolarAxisTitleItalic = 0
    display.PolarAxes.PolarAxisTitleShadow = 0
    display.PolarAxes.PolarAxisTitleFontSize = 12
    display.PolarAxes.PolarAxisLabelOpacity = 1.0
    display.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
    display.PolarAxes.PolarAxisLabelFontFile = ''
    display.PolarAxes.PolarAxisLabelBold = 0
    display.PolarAxes.PolarAxisLabelItalic = 0
    display.PolarAxes.PolarAxisLabelShadow = 0
    display.PolarAxes.PolarAxisLabelFontSize = 12
    display.PolarAxes.LastRadialAxisTextOpacity = 1.0
    display.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
    display.PolarAxes.LastRadialAxisTextFontFile = ''
    display.PolarAxes.LastRadialAxisTextBold = 0
    display.PolarAxes.LastRadialAxisTextItalic = 0
    display.PolarAxes.LastRadialAxisTextShadow = 0
    display.PolarAxes.LastRadialAxisTextFontSize = 12
    display.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
    display.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
    display.PolarAxes.SecondaryRadialAxesTextFontFile = ''
    display.PolarAxes.SecondaryRadialAxesTextBold = 0
    display.PolarAxes.SecondaryRadialAxesTextItalic = 0
    display.PolarAxes.SecondaryRadialAxesTextShadow = 0
    display.PolarAxes.SecondaryRadialAxesTextFontSize = 12
    display.PolarAxes.EnableDistanceLOD = 1
    display.PolarAxes.DistanceLODThreshold = 0.7
    display.PolarAxes.EnableViewAngleLOD = 1
    display.PolarAxes.ViewAngleLODThreshold = 0.7
    display.PolarAxes.SmallestVisiblePolarAngle = 0.5
    display.PolarAxes.PolarTicksVisibility = 1
    display.PolarAxes.ArcTicksOriginToPolarAxis = 1
    display.PolarAxes.TickLocation = 'Both'
    display.PolarAxes.AxisTickVisibility = 1
    display.PolarAxes.AxisMinorTickVisibility = 0
    display.PolarAxes.ArcTickVisibility = 1
    display.PolarAxes.ArcMinorTickVisibility = 0
    display.PolarAxes.DeltaAngleMajor = 10.0
    display.PolarAxes.DeltaAngleMinor = 5.0
    display.PolarAxes.PolarAxisMajorTickSize = 0.0
    display.PolarAxes.PolarAxisTickRatioSize = 0.3
    display.PolarAxes.PolarAxisMajorTickThickness = 1.0
    display.PolarAxes.PolarAxisTickRatioThickness = 0.5
    display.PolarAxes.LastRadialAxisMajorTickSize = 0.0
    display.PolarAxes.LastRadialAxisTickRatioSize = 0.3
    display.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
    display.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
    display.PolarAxes.ArcMajorTickSize = 0.0
    display.PolarAxes.ArcTickRatioSize = 0.3
    display.PolarAxes.ArcMajorTickThickness = 1.0
    display.PolarAxes.ArcTickRatioThickness = 0.5
    display.PolarAxes.Use2DMode = 0
    display.PolarAxes.UseLogAxis = 0

    return 0


def paraview_image(input_file, output_file, field, field_min, field_max, legend_title):
    '''
    '''

    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'XML Unstructured Grid Reader'
    simulation = XMLUnstructuredGridReader(FileName=[input_file])
    simulation.CellArrayStatus = ['Rank']
    simulation.PointArrayStatus = [field]

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    renderView1.ViewSize = [770, 780]

    # get layout
    layout1 = GetLayout()

    # show data in view
    Display = Show(simulation, renderView1, 'UnstructuredGridRepresentation')

    #set_defaults(Display)

    # reset view to fit data
    #renderView1.ResetCamera()

    # update the view to ensure updated data information
    renderView1.Update()
    renderView1.CameraPosition = [0.0, 0.0, 3.2036135254332487]
    renderView1.CameraParallelScale = 0.8291561935301535
    # set scalar coloring
    #ColorBy(Display, ('POINTS', field))

    # # rescale color and/or opacity maps used to include current data range
    # Display.RescaleTransferFunctionToDataRange(True, False)

    # # show color bar/color legend
    # Display.SetScalarBarVisibility(renderView1, True)

    # # get color transfer function/color map and opacity transfer function/opacity map
    # field_without_periods = field.replace(".", "")
    # field_LUT = GetColorTransferFunction(field_without_periods)
    # field_PWF = GetOpacityTransferFunction(field_without_periods)

    # # Rescale and invert
    # field_LUT.RescaleTransferFunction(field_min, field_max)
    # field_PWF.RescaleTransferFunction(field_min, field_max)
    # field_LUT.InvertTransferFunction()

    # # get color legend/bar for diagnostic_quantitiesprojectedCauchy_stress_zzLUT in view renderView1
    # ColorBar = GetScalarBar(field_LUT, renderView1)

    # # Properties modified on diagnostic_quantitiesprojectedCauchy_stress_zzLUTColorBar
    # ColorBar.ScalarBarLength = 0.5

    # # Properties modified on diagnostic_quantitiesprojectedCauchy_stress_zzLUTColorBar
    # ColorBar.Title = legend_title

    # # current camera placement for renderView1
    # renderView1.CameraPosition = [16280.24330465352, 16277.759211197848, 8875.0]
    # renderView1.CameraFocalPoint = [3000.0, 3000.0, 2742.5000000000005]
    # renderView1.CameraViewUp = [-0.22046996075558448, -0.2185805334116402, 0.9505869485838283]
    # renderView1.CameraParallelScale = 5117.470254813933

    # # save screenshot
    # SaveScreenshot(output_file, renderView1, ImageResolution=[770, 780],
        # FontScaling='Scale fonts proportionally',
        # OverrideColorPalette='',
        # StereoMode='No change',
        # TransparentBackground=0,
        # CompressionLevel='5')
    # SaveScreenshot(filename=output_file, viewOrLayout=renderView1)

    return 0


def get_parser():

    basename = pathlib.Path(__file__).name
    cli_description = "Create a csv containing the extents of a DNS file"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=basename)
    parser.add_argument('-i', '--input-file', type=str, required=True,
        help='The input file to get an image from')
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help='The name of the output image file')
    parser.add_argument('--field', type=str, required=True,
        help='The field to plot')
    parser.add_argument('--field-min', type=float, required=True,
        help='The minimum value of the field to plot')
    parser.add_argument('--field-max', type=float, required=True,
        help='The maximum value of the field to plot')
    parser.add_argument('--legend-title', type=str, required=True,
        help='The string to show in the legend')
    parser.add_argument('--camera-position', nargs="+", required=False,
        help='The position of the camera')

    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    print(args)
    print('we parsed something!')
    sys.exit(paraview_image(input_file=args.input_file,
                            output_file=args.output_file,
                            field=args.field,
                            field_min=args.field_min,
                            field_max=args.field_max,
                            legend_title=args.legend_title,
                            ))