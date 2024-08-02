# class to hold the details of SNAP calibration
from mantid.simpleapi import *
import os

#add this folder to the path
import sys
sys.path.append("/SNS/SNAP/shared/Malcolm/code/SNAPRedScripted/")

import json

#import SNAPTools some useful functions for snap data
import SNAPTools as snp
import importlib
importlib.reload(snp)

sxlCalibHome = '/SNS/SNAP/shared/Calibration/SingleCrystal/'
calibrantLibrary = f"{sxlCalibHome}/CalibrantSamples/"

class create():
    
    def initState(self):

    #     #checks if state exists and that it contains a default DetCal file
    #     #if neither exists, they'll be created

        if not os.path.exists(self.calDir):
            # attempt to make directory
            try:
                os.makedirs(self.calDir)
                print(f"Created new state directory:{self.calDir}")
            except:
                print(f"ERROR couldn/'t initialise directory:{self.calDir}")
                print("Check write priviledges")
                return

        self.detCalPath = f"{self.calDir}Default.DetCal"

        if not os.path.exists(self.detCalPath):
            self.makeDetCal()
            print(f"Initialised state")
            return
        else:
            print(f"State available")

        return
    
    def __init__(self,runNumber):

        self.isLite = False #not using lite mode for now.
        instDict = snp.loadSNAPInstPrm()

        self.Inst = instDict["name"]
        self.calDirectory=sxlCalibHome

        self.stateID,self.stateDict,errorState= snp.StateFromRunFunction(runNumber)
        if errorState['value'] != 0: #something went wrong
            print(f"Error in {errorState['function']}")
            return

        self.nxsFile = self.stateDict["nxsFile"]
        #buildPath to output folder
        self.outdir = sxlCalibHome + self.stateID + '/' #don't like this name to making an alias
        self.calDir = sxlCalibHome + self.stateID + '/'

        #initialise state if necessary
        self.initState()

        outputFileType='detcal'


    # def setCalibrant(self,jsonName):
    #     #TODO

    

    def makeDetCal(self):

        print("generating instrument geometry")
        #create a DetCal file and write to calDir
        SNAP = LoadEmptyInstrument(InstrumentName="SNAP")
        

        pars = {
        "det_arc1":str(self.stateDict["det_arc1"]),
        "det_lin1":str(self.stateDict["det_lin1"]),    
        "det_arc2":str(self.stateDict["det_arc2"]),
        "det_lin2":str(self.stateDict["det_lin2"])
        }

        for key in pars:
            AddSampleLog(Workspace="SNAP",
            LogName=key,LogText=pars[key],
            LogType='Number Series')
            
        LoadInstrument(Workspace="SNAP",MonitorList='-1,1179648', RewriteSpectraMap='False',InstrumentName='SNAP')
        
        try:
            SaveIsawDetCal(InputWorkspace="SNAP",
                    Filename=self.detCalPath)
        except:
            print(f"Error: couldn/'t create DetCal at :{self.detCalPath}")
            print("Do you have write priviledges there?")
        
        DeleteWorkspace(Workspace="SNAP")

        return



