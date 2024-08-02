import numpy as np
import matplotlib.pyplot as plt
from mantid.simpleapi import *
import calibrationObject

#create calibration object called "cal" using the run number TODO: and a calibrant
cal = calibrationObject.create(61325)

#this is autogenerate the calibration directory for us
outdir = cal.calDir
