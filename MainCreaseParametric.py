"""
Main Procedure to conduct the parametric analysis
including CreaseArrayInput.py, RunABAQSU.py and odbWrite.py
The load is subjected to the part middle edge
"""
from CreaseArrayInput import CreaseArray
from RunABAQUS import RunABAQUS
from DataProcess import DataPostProcess
import os
import shutil
import pandas as pd
from math import *
import numpy as np
from scipy import stats

"""
Working directory settings, and copy the post process file (.py)
to the calculation directory
"""
PathNow = os.getcwd()
print("The Current program code path", PathNow)
Path = 'C:\ABAQUSTemp\OrigamiArray\Simulation1003'
# shutil.rmtree(Path)
# os.mkdir(Path)
shutil.copyfile("OdbWriteStep.py", Path+"\OdbWriteStep.py")
os.chdir(Path)
PathNew = os.getcwd()
print("The Calculation analysis save path", PathNew)

"""
Input Modulus for parametric analysis
parameter Stiffness: the crease stiffness per length
Parameter SectorAngle: The plate geometry, [width, height, Sector Angle]
the crease length is the height, and  Sector angle between line 1/4 to 1/2, default = 90
Parameter InitialAngle: Define the initial rest angle of the single crease
Parameter dx: mesh size
Parameter foldingstate: 1 is unfolding and 0 is the folding process,default = 1
Parameter Loading region: for the ratio of load region to the edge set. It is the middle region. 0.5 means the middle 50% are loaded
CreaseRegion: defined to expression the cut region. [0.2 0.5] expresses that there is no connection during the region [0.2 0.5]
ArrayNumber: The count of basic origami unit
"""


def Crease(jobnumber, Stiffness, height, SectorAngle, InitialAngle, width, thickness, Loadregion, ArrayNumber):

    jobname = "CreaseArray"+str(jobnumber)
    with open('readmejob.txt',"w") as f:
        f.write(jobname)
    dx = 5.0
    Model = CreaseArray(jobname, Stiffness, dx=dx, width=width, height=height,
                        SectorAngle=SectorAngle, InitialAngle=InitialAngle, thickness=thickness,
                        foldingState=1, Loadregion=Loadregion, Creaseregion=[0.8, 0.2], ArrayNumber = ArrayNumber)
    Model.Inpwrite()
    ABAQUSModel = RunABAQUS(jobname)
    ABAQUSModel.Simulate()
    # for Step in range(1, 2):
    #     stepname = "Step-"+str(Step)
    #     with open('readmeStep.txt',"w") as f:
    #         f.write(stepname)
    #     os.system('abaqus cae noGUI=OdbWriteStep.py')
    #     DataModel = DataPostProcess(jobname, jobnumber, Stiffness, dx, InitialAngle)
    #     DataModel.Process()

"""
Main procedure for Crease Array analysis
"""

Allseresult = pd.DataFrame(columns=['jobnumber', 'Stiffness', 'height', 'width', 'SectorAngle',
                                    'InitialAngle', 'thickness', 'LoadRegion', 'ArrayNumber'])
jobnumber = 1
# STi = [0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1]
STi = [0.05]
for i in range(len(STi)):
    Stiffness = STi[i]
    for HEj in range(8, 9):
        height = 100 * HEj
        for Wm in range(2, 3):
            width = 40.0 * Wm
            for Sn in range(6, 7):
                SectorAngle = 15 / 180 * pi * Sn
                for Agk in range(15, 16, 1):
                    InitialAngle = 10 / 180 * pi * Agk
                    for Thp in range(2, 3, 1):
                        thickness = 0.05 * Thp
                        for Lp in range(16, 18, 2):
                            Loadregion = 0.05 * Lp
                            for AN in range(10, 11, 1):
                                ArrayNumber = AN
                                # Crease Procedure
                                jobnumber = jobnumber + 1
                                print(jobnumber)
                                Crease(jobnumber, Stiffness, height, SectorAngle, InitialAngle, width, thickness,
                                       Loadregion, ArrayNumber)
                                Allseresult.loc[jobnumber] = [jobnumber, Stiffness, height, width, SectorAngle*180/pi,
                                                              InitialAngle*180/pi, thickness, Loadregion, ArrayNumber]

Allseresult.to_excel("LastALLSE.xlsx", index=False)
