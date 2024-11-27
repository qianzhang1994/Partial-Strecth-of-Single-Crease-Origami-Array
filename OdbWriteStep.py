from abaqus import *
from odbAccess import *

"""
Read the job Information
and the Step Information
Set-1 corresponds to the unfolding to the plane
Set-2 corresponds to the external tension of the plane
"""
f = open('readmejob.txt','r')
req = f.readline()
f.close()
path = req + '.odb'
an_odb_object = openOdb(path=path)

f = open('readmeStep.txt','r')
step = f.readline()
f.close()
StepNumber = step

"""
The ALLSE of the plate and the total systems
as PlateASSE and TotalALLSE
"""
Platedata = an_odb_object.steps[StepNumber].historyRegions["ElementSet SET-PLATE"].historyOutputs['ALLSE'].data
Totaldata = an_odb_object.steps[StepNumber].historyRegions["Assembly ASSEMBLY"].historyOutputs['ALLSE'].data
Framenumber = len(Platedata)
DataFile=open('Data_PlateALLSE.txt', 'w')
DataFile.write('time,PlateALLSE\n')
for time, Data in Platedata:
    DataFile.write('%10.4E, %10.4E\n' % (time, Data))
DataFile.close()

DataFile=open('Data_TotalALLSEl.txt','w')
DataFile.write('TotalALLSE\n')
for time, Data in Totaldata:
    DataFile.write('%10.4E\n' % (Data))
DataFile.close()

"""
Record all the output result documents name in a file
Resultfile_Label
"""
ResultLabel = open('Resultfile_Label.txt','w')
ResultLabel.write("Data_PlateALLSE\n")
ResultLabel.write("Data_TotalALLSEl\n")
ResultLabel.close()


"""
Read the OUTPUTX and OUTPUTY from the SIngleCreaseInput.py
to obtain the configuration of system from The X direction and Y direction.
"""
countX = len(open("OutputX.txt").readlines())
countY = len(open("OutputY.txt").readlines())

OutPutXList = []
for line in open('OutputX.txt').readlines():
    OutPutXList.append(line.strip("\n"))

OutPutYList = []
for line in open('OutputY.txt').readlines():
    OutPutYList.append(line.strip("\n"))

"""
OutputX and OUTPUTY for the node set coordinates in field output
"""
for i in range(countX):
    ResultLabel = open('Resultfile_Label.txt', 'a+')
    ResultLabel.write("Data_CoordX" + str(i+1) + "\n")
    ResultLabel.close()

    DataFile = open("Data_CoordX" + str(i + 1) + ".txt", 'w')
    NodeLabel = OutPutXList[i]
    NodeLabel = int(NodeLabel)
    line = "XCoordX." +str(NodeLabel) + ",XCoordY." + str(NodeLabel) + ",XCoordZ." + str(NodeLabel)
    DataFile.write(line)
    DataFile.write('\n')
    for k in range(Framenumber):
        Creasedata = an_odb_object.steps[StepNumber].frames[k].fieldOutputs['COORD'].values[int(NodeLabel) - 1].data
        line = str(Creasedata[0]) + "," + str(Creasedata[1]) + "," + str(Creasedata[2])
        DataFile.write(line)
        DataFile.write('\n')
    DataFile.close()

for i in range(countY):
    ResultLabel = open('Resultfile_Label.txt', 'a+')
    ResultLabel.write("Data_CoordY" + str(i+1) + "\n")
    ResultLabel.close()

    DataFile = open("Data_CoordY" + str(i + 1) + ".txt", 'w')
    NodeLabel = OutPutYList[i]
    NodeLabel = int(NodeLabel)
    line = "YCoordX." + str(NodeLabel) + ",YCoordY." + str(NodeLabel) + ",YCoordZ." + str(NodeLabel)
    DataFile.write(line)
    DataFile.write('\n')
    for k in range(Framenumber):
        Creasedata = an_odb_object.steps[StepNumber].frames[k].fieldOutputs['COORD'].values[int(NodeLabel) - 1].data
        line = str(Creasedata[0]) + "," + str(Creasedata[1]) + "," + str(Creasedata[2])
        DataFile.write(line)
        DataFile.write('\n')
    DataFile.close()

"""
Connection Pair moment output (CM1) in history output
"""
for i in range(countY):
    ResultLabel = open('Resultfile_Label.txt', 'a+')
    ResultLabel.write("Data_CreasePair"+str(i+1) + "\n")
    ResultLabel.close()

    DataFile = open("Data_CreasePair"+str(i+1)+".txt", 'w')
    CTMdata = an_odb_object.steps[StepNumber].historyRegions["Element ASSEMBLY." + str(i+1)].historyOutputs[
        'CTM1'].data
    DataFile.write("Element ASSEMBLY." + str(i+1)+"\n")
    for time, Data in CTMdata:
        DataFile.write('%10.4E\n' % Data)
    DataFile.close()


