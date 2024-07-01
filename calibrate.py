import numpy as np
import matplotlib.pyplot as plt
from mantid.simpleapi import *

outdir = '/SNS/TOPAZ/IPTS-31189/shared/YAG/calibration_2024'
calibration_file = '/SNS/TOPAZ/shared/calibration/2022B/TOPAZ_2022B.DetCal'

a = 11.9386
b = 11.9386
c = 11.9386
alpha = 90
beta = 90
gamma = 90

IPTS = 31189
runs = range(46643, 46692+1)

wavelength = [0.5, 3.5]
Q_max = 4*np.pi/wavelength[0]

d_max = max([a,b,c])

max_peaks = 1000
density_threshold = 10000

peak_radii = [0.1, 0.12, 0.15]
sig_noise = 50

for i, run in enumerate(runs):

    filename = '/SNS/TOPAZ/IPTS-{}/nexus/TOPAZ_{}.nxs.h5'.format(IPTS, run)

    LoadEventNexus(Filename=filename,
                   OutputWorkspace='data')

    LoadIsawDetCal(InputWorkspace='data',
                   Filename=calibration_file)

    SetGoniometer(Workspace='data',
                  Axis0='BL12:Mot:goniokm:omega,0,1,0,1',
                  Axis1='BL12:Mot:goniokm:chi,0,0,1,1',
                  Axis2='BL12:Mot:goniokm:phi,0,1,0,1',
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

LoadEmptyInstrument(InstrumentName='TOPAZ',
                    OutputWorkspace='TOPAZ')

LoadParameterFile(Workspace='TOPAZ', 
                  Filename=os.path.join(outdir, 'calibration.xml'))

sample_pos = mtd['TOPAZ'].getInstrument().getComponentByName('sample-position').getPos()

for bank in np.unique(mtd['calibration_ws'].column(13)):
    MoveInstrumentComponent(Workspace='TOPAZ', 
                            ComponentName=bank, 
                            X=-sample_pos[0], Y=-sample_pos[1], Z=-sample_pos[2], 
                            RelativePosition=True)

MoveInstrumentComponent(Workspace='TOPAZ', 
                        ComponentName='sample-position', 
                        X=-sample_pos[0], Y=-sample_pos[1], Z=-sample_pos[2], 
                        RelativePosition=True)

MoveInstrumentComponent(Workspace='TOPAZ', 
                        ComponentName='moderator', 
                        X=0, Y=0, Z=-sample_pos[2], 
                        RelativePosition=True)

ApplyInstrumentToPeaks(InputWorkspace='calibration_ws', 
                       InstrumentWorkspace='TOPAZ',
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

LoadEmptyInstrument(InstrumentName='TOPAZ',
                    OutputWorkspace='TOPAZ')

LoadParameterFile(Workspace='TOPAZ', 
                  Filename=os.path.join(outdir, 'calibration.xml'))

sample_pos = mtd['TOPAZ'].getInstrument().getComponentByName('sample-position').getPos()

for bank in np.unique(mtd['detcal_ws'].column(13)):
    MoveInstrumentComponent(Workspace='TOPAZ', 
                            ComponentName=bank, 
                            X=-sample_pos[0], Y=-sample_pos[1], Z=-sample_pos[2], 
                            RelativePosition=True)

MoveInstrumentComponent(Workspace='TOPAZ', 
                        ComponentName='sample-position', 
                        X=-sample_pos[0], Y=-sample_pos[1], Z=-sample_pos[2], 
                        RelativePosition=True)

MoveInstrumentComponent(Workspace='TOPAZ', 
                        ComponentName='moderator', 
                        X=0, Y=0, Z=-sample_pos[2], 
                        RelativePosition=True)

ApplyInstrumentToPeaks(InputWorkspace='detcal_ws', 
                       InstrumentWorkspace='TOPAZ',
                       OutputWorkspace='detcal_ws')

# ---

CloneWorkspace(InputWorkspace='peaks',
               OutputWorkspace='cal')

LoadEmptyInstrument(InstrumentName='TOPAZ',
                    OutputWorkspace='TOPAZ')

LoadParameterFile(Workspace='TOPAZ', 
                  Filename=os.path.join(outdir, 'cal.xml'))

ApplyInstrumentToPeaks(InputWorkspace='cal', 
                       InstrumentWorkspace='TOPAZ',
                       OutputWorkspace='cal')
                      