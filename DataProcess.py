"""
Data Processing Modulus
Combine all result in Pandas.DataFrame
df: columns
time/plateAllSE/TotalALLSE/XCOORD/YCOORD/Conection Pair Moment
XCOORD: the coordinates of the nodes in middle plate 0 along the width directin
YCOORD: the coordinate of the nodes in the connection edge of plate 0
XCOORD and YCOORD are calculated to determine the configuration.
counX and CountY are the number of nodes for XCOORD and YCOORD
The number of the connection pair is the same as the YCOORD set
The total series of the dataframe is
1 + 1 + 1 + countX *3 +countY * 3 + CountY
"""

import pandas as pd
from math import *
import matplotlib.pyplot as plt


class DataPostProcess(object):
    def __init__(self, jobname, jobnumber, Stiffness, dx, InitialAngle):
        self.jobname = jobname
        self.jobnumber = jobnumber
        self.stiffness = Stiffness
        self.dx = dx
        self.InitialAngle = InitialAngle

    def Process(self):
        jobname = self.jobname
        jobnumber = self.jobnumber

        # countX = len(open("OutputX.txt").readlines())
        # countY = len(open("OutputY.txt").readlines())

        f = open('readmeStep.txt', 'r')
        step = f.readline()
        f.close()
        StepNumber = step

        Resultfilelist = []
        f = open('Resultfile_Label.txt', 'r')
        for line in f.readlines():
            Resultfilelist.append(line.strip("\n"))

        fpath = Resultfilelist[0] + ".txt"
        df = pd.read_csv(fpath)
        for i in range(1, len(Resultfilelist)):
            fpath_temp = Resultfilelist[i] + ".txt"
            df_temp = pd.read_csv(fpath_temp)
            df = pd.concat([df, df_temp], axis=1)

        resultname = "CreaseArray" + str(jobnumber) + str(StepNumber) + ".xlsx"
        df.to_excel(resultname, index=False)

        """
        Post process for the results from the excel file
        (1) The strain energy curves for plate and total;
        (2) The X configuration for the specified frame;
        (3) The Y Configuration for the specified frame;
        (4) The Crease spring Angle results for the specified frame
        """
        fpath = "CreaseArray" + str(jobnumber) + str(StepNumber) + ".xlsx"
        df = pd.read_excel(fpath)

        """
        (1) The strain energy curves for plate and total;
        """
        time = df.loc[:, "time"]
        print(time)
        PlateALLSE = df.loc[:, "PlateALLSE"]
        TotalALLSE = df.loc[:, "TotalALLSE"]
        CreaseALLSE = TotalALLSE - PlateALLSE
        fig, ax = plt.subplots()
        plt.xlabel('Time')
        plt.ylabel('ALLSE')
        line1, = ax.plot(time, PlateALLSE, linewidth=2.0)
        line2, = ax.plot(time, TotalALLSE, linewidth=2.0)
        line3, = ax.plot(time, CreaseALLSE, linewidth=2.0)
        plt.legend([line1, line2, line3], ["PlateALLSE", "TotalALLSE", "CreaseALLSE"])
        figname = "CreaseArray" + str(jobnumber) + str(StepNumber) + "ALLSE"
        fig.suptitle(figname)
        fig.savefig(figname + '.jpg')
        plt.show()

        """
        (2) The X configuration for the specified frame;
        Specified frame list Frames[] and the simulation process percent Prcess[]
        """
        CreaseXData = pd.DataFrame()
        timelist = time.tolist()
        timeFram = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
        Frames = []
        for nu in range(len(timeFram)):
            answer = []
            for ni in timelist:
                answer.append(abs(timeFram[nu] - ni))
            Frames.append(answer.index(min(answer)))

        fig, ax = plt.subplots()
        figname = "CreaseArray" + str(jobnumber) + str(StepNumber) + "XCoord"
        fig.suptitle(figname)
        plt.xlabel('Coordinate X')
        plt.ylabel('Coordinate Z')
        for i in range(len(Frames)):
            XListX = list(filter(lambda x: "XCoordX" in x, df.columns.values.tolist()))
            XListZ = list(filter(lambda x: "XCoordZ" in x, df.columns.values.tolist()))
            dfXListX = df[XListX].loc[Frames[i], :]
            dfXListZ = df[XListZ].loc[Frames[i], :]
            label = str(floor(timeFram[i] * 100)) + "%"
            CreaseXData['dfXListX' + str(i)] = dfXListX.tolist()
            CreaseXData['dfXListZ' + str(i)] = dfXListZ.tolist()
            ax.plot(dfXListX, dfXListZ, linewidth=2.0, label=label)
        ax.legend()
        fig.savefig(figname + '.jpg')
        plt.show()
        figname = "CreaseXData" + str(jobnumber) + str(StepNumber) + '.xlsx'
        CreaseXData.to_excel(figname, index=False)

        """
        (3) The Y configuration for the specified frame;
        """
        CreaseYData = pd.DataFrame()
        fig, ax = plt.subplots()
        figname = "CreaseArray" + str(jobnumber)+ str(StepNumber) + "YCoord"
        fig.suptitle(figname)
        plt.xlabel('Coordinate Y')
        plt.ylabel('Coordinate Z')
        for i in range(len(Frames)):
            YListY = list(filter(lambda x: "YCoordY" in x, df.columns.values.tolist()))
            YListZ = list(filter(lambda x: "YCoordZ" in x, df.columns.values.tolist()))
            dfYListY = df[YListY].loc[Frames[i], :]
            dfYListZ = df[YListZ].loc[Frames[i], :]
            CreaseYData['dfYListY' + str(i)] = dfYListY.tolist()
            CreaseYData['dfYListZ' + str(i)] = dfYListZ.tolist()
            label = str(floor(timeFram[i] * 100)) + "%"
            ax.plot(dfYListY, dfYListZ, linewidth=2.0, label=label)
        ax.legend()
        fig.savefig(figname + '.jpg')
        plt.show()

        """
        (4) The Crease Spring Angle configuration for the specified frame;
        """
        fig, ax = plt.subplots()
        figname = "CreaseArraye" + str(jobnumber) + str(StepNumber) + "CreaseMomentAngle"
        fig.suptitle(figname)
        plt.xlabel('Coordinate Y')
        plt.ylabel('Crease Angle')
        for i in range(len(Frames)):
            YListY = list(filter(lambda x: "YCoordY" in x, df.columns.values.tolist()))
            CMoment = list(filter(lambda x: "Element ASSEMBLY" in x, df.columns.values.tolist()))
            dfYListY = df[YListY].loc[Frames[i], :]
            dfCMoment = df[CMoment].loc[Frames[i], :]
            dfCangle = - dfCMoment / (self.stiffness * self.dx) + self.InitialAngle
            label = str(floor(timeFram[i] * 100)) + "%"
            CreaseYData['dfCangle' + str(i)] = dfCangle.tolist()
            ax.plot(dfYListY, dfCangle, linewidth=2.0, label=label)
        ax.legend()
        fig.savefig(figname + '.jpg')
        plt.show()
        figname = "CreaseYData" + str(jobnumber) + str(StepNumber) + '.xlsx'
        CreaseYData.to_excel(figname, index=False)

        """
        (5) The X configuration for the last frame
        Specified frame list Frames[] and the simulation process percent Prcess[]
        """
        fig, ax = plt.subplots()
        figname = "CreaseArray" + str(jobnumber) + str(StepNumber) +"XCoordLast"
        fig.suptitle(figname)
        plt.xlabel('Coordinate X')
        plt.ylabel('Coordinate Z')
        XListX = list(filter(lambda x: "XCoordX" in x, df.columns.values.tolist()))
        XListZ = list(filter(lambda x: "XCoordZ" in x, df.columns.values.tolist()))
        dfXListX = df[XListX].loc[df.index[-1], :]
        dfXListZ = df[XListZ].loc[df.index[-1], :]
        ax.plot(dfXListX, dfXListZ, linewidth=2.0)
        fig.savefig(figname + '.jpg')
        plt.show()

        # """
        # (6) The Y configuration parameter for the crease (average);
        # """
        # YListZ = list(filter(lambda x: "YCoordZ" in x, df.columns.values.tolist()))
        # dfYListZ = df[YListZ]
        # dfYListZNew = abs(dfYListZ)
        # NlistZpara = dfYListZNew.mean(axis=1)
        #
        # fig, ax = plt.subplots()
        # figname = "SingleCrease" + str(jobnumber) + "CreaseZ"
        # fig.suptitle(figname)
        # plt.xlabel('Time')
        # plt.ylabel('Coordinate Z Average')
        # ax.plot(time, NlistZpara, linewidth=2.0)
        # fig.savefig(figname + '.jpg')
        # plt.show()
        #
        # CreaseZ = pd.DataFrame()
        # CreaseZ['time'] = time.tolist()
        # CreaseZ['CreaseZ'] = NlistZpara.tolist()
        # figname = "CreaseZ" + str(jobnumber)+'.xlsx'
        # CreaseZ.to_excel(figname, index=False)


        # return AllSERatioResult
