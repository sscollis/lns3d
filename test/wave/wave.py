# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

import os
import glob

def sortKeyFunc(s):
  return int(os.path.basename(s).split('.')[2]) 

# create the list of q-files
#
# Requires the LNS3D_DIR be set in the environment
#
LNS3D_DIR = os.getenv("LNS3D_DIR")
os.chdir(LNS3D_DIR+'/test/wave')
files=glob.glob("output.q.*")
files.sort(key=sortKeyFunc)

# create a new 'PLOT3D Reader'
gridxyz = PLOT3DReader(FileName='grid.xyz',
    QFileName=[''],
    FunctionFileName='')
gridxyz.Functions = ['Scalar - Pressure']

# Properties modified on gridxyz
gridxyz.QFileName = files
gridxyz.FunctionFileName = ''

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [2026, 1546]

# show data in view
gridxyzDisplay = Show(gridxyz, renderView1)

# get color transfer function/color map for 'Density'
densityLUT = GetColorTransferFunction('Density')

# get opacity transfer function/opacity map for 'Density'
densityPWF = GetOpacityTransferFunction('Density')

# trace defaults for the display properties.
gridxyzDisplay.Representation = 'Surface'
gridxyzDisplay.ColorArrayName = ['POINTS', 'Density']
gridxyzDisplay.LookupTable = densityLUT
gridxyzDisplay.OSPRayScaleArray = 'Density'
gridxyzDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
gridxyzDisplay.SelectOrientationVectors = 'Momentum'
gridxyzDisplay.ScaleFactor = 50.0
gridxyzDisplay.SelectScaleArray = 'Density'
gridxyzDisplay.GlyphType = 'Arrow'
gridxyzDisplay.GlyphTableIndexArray = 'Density'
gridxyzDisplay.GaussianRadius = 2.5
gridxyzDisplay.SetScaleArray = ['POINTS', 'Density']
gridxyzDisplay.ScaleTransferFunction = 'PiecewiseFunction'
gridxyzDisplay.OpacityArray = ['POINTS', 'Density']
gridxyzDisplay.OpacityTransferFunction = 'PiecewiseFunction'
gridxyzDisplay.DataAxesGrid = 'GridAxesRepresentation'
gridxyzDisplay.PolarAxes = 'PolarAxesRepresentation'
gridxyzDisplay.ScalarOpacityFunction = densityPWF
gridxyzDisplay.ScalarOpacityUnitDistance = 58.291548897025436

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
gridxyzDisplay.OSPRayScaleFunction.Points = [1.980318699596129, 0.0, 0.5, 0.0, 2.698034067905081, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
gridxyzDisplay.ScaleTransferFunction.Points = [-0.12499051376435383, 0.0, 0.5, 0.0, 0.125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
gridxyzDisplay.OpacityTransferFunction.Points = [-0.12499051376435383, 0.0, 0.5, 0.0, 0.125, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# show color bar/color legend
gridxyzDisplay.SetScalarBarVisibility(renderView1, True)

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# update the view to ensure updated data information
renderView1.Update()

# rename source object
RenameSource('Wave', gridxyz)

animationScene1.Play()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [250.0, 250.0, 1366.0254037844388]
renderView1.CameraFocalPoint = [250.0, 250.0, 0.0]
renderView1.CameraParallelScale = 353.5533905932738

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
