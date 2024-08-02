import numpy as np
import matplotlib.pyplot as plt
from mantid.simpleapi import *
import calibrationObject as calObj
import sys

crystalCalibrant = "sapphire"
outputDetcalName = "test"
runs = range(61325,61346+1)

#calibration calculation parameters (for fiddling)
wavelength = [0.5, 3.5] #TODO: calculate in calObj
max_peaks = 1000
density_threshold = 10000
peak_radii = [0.1, 0.12, 0.15]
sig_noise = 50


#############################################################
# DON' EDIT BELOW
#############################################################

#beamline axes TODO: (Generalise...create a goniometer set-up file to store? This 
#the goniometer could be automated via an input to calibrationObject 

beamLineAxis0 ='BL3:Mot:omega,0,1,0,1'
beamLineAxis1 ='BL3:Mot:phi,0.707,0.707,0,1'

#create calibration object and use it to set parameters
cal = calObj.create(runs[0],crystalCalibrant,verbose=True)

outdir = cal.calDir #output directory
calibration_file = f"{cal.calDir}Default.DetCal"

a = cal.crystal.a
b = cal.crystal.b
c = cal.crystal.a
alpha = cal.crystal.alpha
beta = cal.crystal.beta
gamma = cal.crystal.gamma


Q_max = 4*np.pi/wavelength[0]
d_max = max([a,b,c])

inst=cal.Inst #TODO: needed?

#check that runs are all from the same state

allID = []
prog = 0
for i,run in enumerate(runs):
    prog = 100*((i+1)/len(runs))
    print(f"checking all {len(runs)} requested runs are from same state: {prog:.1f}% complete", end = "\r")
    tempCal = calObj.create(run,crystalCalibrant)
    allID.append(tempCal.stateID)

stateIDSet = set(allID)

if len(stateIDSet) != 1:
    print("ERROR all runs must be from same state!")
    print(f"There are {len(stateIDSet)} states in run list:")
    for state in stateIDSet:
        print(state)
    sys.exit()

prog = 0
print("Finding and integrating peaks")
for i, run in enumerate(runs):

    print(f"Working on run {run} ({(i+1)}/{len(runs)} runs)", end="\r")
    cal = calObj.create(run,crystalCalibrant)

    filename = cal.nxsFile

    LoadEventNexus(Filename=filename,
                   OutputWorkspace='data')

    LoadIsawDetCal(InputWorkspace='data',
                   Filename=calibration_file)

    SetGoniometer(Workspace='data',
                  Axis0=beamLineAxis0,
                  Axis1=beamLineAxis1,
                #   Axis2=beamLineAxis0, #TODO automatically accommodate 1,2 or 3 axes
                  Average=True)

    ConvertToMD(InputWorkspace='data',
                QDimensions='Q3D',
                dEAnalysisMode='Elastic',
                Q3DFrames='Q_sample',
                MinValues=[-Q_max,-Q_max,-Q_max],
                MaxValues=[+Q_max,+Q_max,+Q_max],
                OutputWorkspace='md')

    FindPeaksMD(InputWorkspace='md',
                PeakDistanceThreshold=2*np.pi/d_max,
                MaxPeaks=max_peaks,
                DensityThresholdFactor=density_threshold,
                OutputWorkspace='peaks_ws')

    IntegratePeaksMD(InputWorkspace='md',
                     PeakRadius=peak_radii[0],
                     BackgroundInnerRadius=peak_radii[1],
                     BackgroundOuterRadius=peak_radii[2],
                     PeaksWorkspace='peaks_ws',
                     OutputWorkspace='peaks_ws',
                     Ellipsoid=True,
                     FixQAxis=True,
                     FixMajorAxisLength=False,
                     UseCentroid=True,
                     MaxIterations=3)

    FilterPeaks(InputWorkspace='peaks_ws', 
                FilterVariable='Signal/Noise', 
                FilterValue=sig_noise, 
                Operator='>',
                OutputWorkspace='peaks_ws')

    if i == 0:
        CloneWorkspace(InputWorkspace='peaks_ws',
                       OutputWorkspace='peaks')
    else:
        CombinePeaksWorkspaces(LHSWorkspace='peaks_ws', 
                               RHSWorkspace='peaks', 
                               OutputWorkspace='peaks')

FindUBUsingLatticeParameters(PeaksWorkspace='peaks',
                             a=a,
                             b=b,
                             c=c,
                             alpha=alpha,
                             beta=beta,
                             gamma=gamma,
                             NumInitial=50,
                             Iterations=10)

IndexPeaks(PeaksWorkspace='peaks', Tolerance=0.05)

FilterPeaks(InputWorkspace='peaks', 
            FilterVariable='h^2+k^2+l^2', 
            FilterValue=0, 
            Operator='!=',
            OutputWorkspace='peaks')

FilterPeaks(InputWorkspace='peaks', 
            FilterVariable='QMod', 
            FilterValue=0, 
            Operator='>',
            OutputWorkspace='peaks')

FilterPeaks(InputWorkspace='peaks',
            BankName='',
            Criterion='!=',
            OutputWorkspace='peaks')

SCDCalibratePanels(PeakWorkspace='peaks',
                   RecalculateUB=True,
                   a=a,
                   b=b,
                   c=c,
                   alpha=alpha,
                   beta=beta,
                   gamma=gamma,
                   OutputWorkspace='calibration_table',
                   DetCalFilename=os.path.join(outdir, 'calibration.DetCal'),
                   CSVFilename=os.path.join(outdir, 'calibration.csv'),
                   XmlFilename=os.path.join(outdir, 'calibration.xml'),
                   CalibrateT0=False,
                   SearchRadiusT0=10,
                   CalibrateL1=True,
                   SearchRadiusL1=0.2,
                   CalibrateBanks=True,
                   SearchRadiusTransBank=0.5,
                   SearchRadiusRotXBank=5,
                   SearchRadiusRotYBank=5,
                   SearchRadiusRotZBank=5,
                   VerboseOutput=True,
                   SearchRadiusSamplePos=0.1,
                   TuneSamplePosition=True,
                   CalibrateSize=False,
                   SearchRadiusSize=0.1,
                   FixAspectRatio=True)

CloneWorkspace(InputWorkspace='peaks',
               OutputWorkspace='calibration_ws')

LoadEmptyInstrument(InstrumentName='inst',
                    OutputWorkspace='inst')

LoadParameterFile(Workspace='inst', 
                  Filename=os.path.join(outdir, 'calibration.xml'))

sample_pos = mtd['inst'].getInstrument().getComponentByName('sample-position').getPos()

for bank in np.unique(mtd['calibration_ws'].column(13)):
    MoveInstrumentComponent(Workspace='inst', 
                            ComponentName=bank, 
                            X=-sample_pos[0], Y=-sample_pos[1], Z=-sample_pos[2], 
                            RelativePosition=True)

MoveInstrumentComponent(Workspace='inst', 
                        ComponentName='sample-position', 
                        X=-sample_pos[0], Y=-sample_pos[1], Z=-sample_pos[2], 
                        RelativePosition=True)

MoveInstrumentComponent(Workspace='inst', 
                        ComponentName='moderator', 
                        X=0, Y=0, Z=-sample_pos[2], 
                        RelativePosition=True)

ApplyInstrumentToPeaks(InputWorkspace='calibration_ws', 
                       InstrumentWorkspace='inst',
                       OutputWorkspace='calibration_ws')

SCDCalibratePanels(PeakWorkspace='calibration_ws',
                   RecalculateUB=False,
                   a=a,
                   b=b,
                   c=c,
                   alpha=alpha,
                   beta=beta,
                   gamma=gamma,
                   OutputWorkspace='tmp',
                   DetCalFilename=os.path.join(outdir, 'cal.DetCal'),
                   CSVFilename=os.path.join(outdir, 'cal.csv'),
                   XmlFilename=os.path.join(outdir, 'cal.xml'),
                   CalibrateT0=False,
                   SearchRadiusT0=10,
                   CalibrateL1=False,
                   SearchRadiusL1=0.2,
                   CalibrateBanks=False,
                   SearchRadiusTransBank=0.5,
                   SearchRadiusRotXBank=5,
                   SearchRadiusRotYBank=5,
                   SearchRadiusRotZBank=5,
                   VerboseOutput=False,
                   SearchRadiusSamplePos=0.1,
                   TuneSamplePosition=False,
                   CalibrateSize=False,
                   SearchRadiusSize=0.1,
                   FixAspectRatio=False)

# ---

CloneWorkspace(InputWorkspace='peaks',
               OutputWorkspace='detcal_ws')

LoadEmptyInstrument(InstrumentName='inst',
                    OutputWorkspace='inst')

LoadParameterFile(Workspace='inst', 
                  Filename=os.path.join(outdir, 'calibration.xml'))

sample_pos = mtd['inst'].getInstrument().getComponentByName('sample-position').getPos()

for bank in np.unique(mtd['detcal_ws'].column(13)):
    MoveInstrumentComponent(Workspace='inst', 
                            ComponentName=bank, 
                            X=-sample_pos[0], Y=-sample_pos[1], Z=-sample_pos[2], 
                            RelativePosition=True)

MoveInstrumentComponent(Workspace='inst', 
                        ComponentName='sample-position', 
                        X=-sample_pos[0], Y=-sample_pos[1], Z=-sample_pos[2], 
                        RelativePosition=True)

MoveInstrumentComponent(Workspace='inst', 
                        ComponentName='moderator', 
                        X=0, Y=0, Z=-sample_pos[2], 
                        RelativePosition=True)

ApplyInstrumentToPeaks(InputWorkspace='detcal_ws', 
                       InstrumentWorkspace='inst',
                       OutputWorkspace='detcal_ws')

# ---

CloneWorkspace(InputWorkspace='peaks',
               OutputWorkspace='cal')

LoadEmptyInstrument(InstrumentName='inst',
                    OutputWorkspace='inst')

LoadParameterFile(Workspace='inst', 
                  Filename=os.path.join(outdir, 'cal.xml'))

ApplyInstrumentToPeaks(InputWorkspace='cal', 
                       InstrumentWorkspace='inst',
                       OutputWorkspace='cal')
                      