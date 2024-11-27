"""
Created on June 21 2022
Modified On July 7 2022
@author: Zhang Qian
"""

from math import *
import numpy as np


class CreaseArray(object):
    """
    A crease array class
    This procedure is used for the crease array simulation for Edge Displacement actuator.
    Tension and Compression.
    Crease Stiffness
            *              *
           * *           *   *            ---
          *   *        *       *          ---
         *     *     *           *        ---
        *       *  *               *      ---
    D<-*         *                  *->D  ---
    There are many plates connected by creases.
    """

    def __init__(self, jobname, Stiffness, dx, width, height, SectorAngle, InitialAngle, thickness, foldingState,
                 Loadregion, Creaseregion, ArrayNumber):
        self.jobname = jobname
        self.stiffness = Stiffness
        self.dx = dx
        self.width = width
        self.height = height
        self.SectorAngle = SectorAngle
        self.InitialAngle = InitialAngle
        self.thickness = thickness
        self.foldingstate = foldingState
        self.loadregion = Loadregion
        self.Creaseregion = Creaseregion
        self.ArrayNumber = ArrayNumber

    def _coordinate(self, geometry, InitialAngle, PlateNumber):
        """
        The vertex number of the plate starts from the left-up node,
        then increases count-clockwise, like
        1, 4
        2, 3
        The node 3 of left plate is selected as the origin point.
        the line 3/4 of left plate is the y-axis
        the 3/4 of left is connected by 2/1 of right plate
        :arg geometry [width, Height and sector angle] of the plate.
        :arg float InitialAngle: the initial angle of the two plate.
        :arg int PlateNumber: 0 for left plate and 1 for right plate
        """
        Rotation = (pi - InitialAngle) / 2
        T = np.array([[cos(Rotation), 0, sin(Rotation)],
                      [0, 1, 0],
                      [-sin(Rotation), 0, cos(Rotation)]])
        width = geometry[0]
        height = geometry[1]
        SecAngle = geometry[2]
        CoordR = np.ones((4, 3))

        CoordR[0] = np.array([0, height, 0])
        CoordR[1] = np.array([0, 0, 0])
        temp3 = np.array([width * sin(SecAngle), -width * cos(SecAngle), 0])
        CoordR[2] = np.dot(T, temp3)
        CoordR[3] = CoordR[2] + CoordR[0]

        if (PlateNumber % 2) == 1:
            Coord = CoordR
        else:
            CoordL = np.ones((4, 3))
            M = np.array([[-1, 0, 0],
                          [0, 1, 0],
                          [0, 0, 1]])
            CoordL[0] = np.dot(M, CoordR[3])
            CoordL[1] = np.dot(M, CoordR[2])
            CoordL[2] = CoordR[1]
            CoordL[3] = CoordR[0]
            Coord = CoordL

        # Translation Matrix in the Horizon Direction
        InitialWidth = width * sin(SecAngle) * cos(Rotation) * 2
        TM = np.array([InitialWidth, 0, 0])
        # The tension plane is defined as the coordinate plane
        InitialHeight = width * sin(SecAngle) * sin(Rotation)
        TZM = np.array([0, 0, InitialHeight])
        CoordArray = Coord + PlateNumber // 2 * TM + TZM

        return CoordArray

    def _Mesh(self, nw, nh, Coord, identifier):
        """
        Only support Quadrilateral mesh
        Calculate the node list and element list in global system
        Aspect_Ratio means the mesh size in the width and height direction.
        nw: width/dx/Aspect_Ratio
        nh: height/dx
        identifier: the plate number, int
        :return: NodeCoord Node list A(,4), (Node Number, x, y, z)
        :return: Element element lsit B(,5),
                (Element Number, Node1, Node2, Node 3, Node4 )
        :return: EdgeSet --Edge Node set for Connection design.
        Edge Number: 0 is for line 1-2, 1 is for line 2-3. 2 is for line 4-3
                3 is for line 1-4. The order of the line is considered.
        """
        NodeCoord = np.zeros(((nw + 1) * (nh + 1), 4))
        Element = np.zeros((nw * nh, 5), dtype=np.int32)
        EdgeSet = []
        # Node Information
        Number = 1
        for i in range(nw + 1):
            start = Coord[0] + (i / nw) * (Coord[3] - Coord[0])
            end = Coord[1] + (i / nw) * (Coord[2] - Coord[1])
            for j in range(nh + 1):
                temp = start + (j / nh) * (end - start)
                Global_Number = Number + identifier * (nw + 1) * (nh + 1)
                NodeCoord[Number - 1] = np.array([Global_Number, temp[0], temp[1], temp[2]])
                Number += 1
        # Element Information
        Number = 1
        for i in range(nw):
            for j in range(nh):
                Global_symbol = identifier * (nw + 1) * (nh + 1)
                Element[Number - 1] = np.array([Number + identifier * nw * nh,
                                                Global_symbol + i * (nh + 1) + j + 1,
                                                Global_symbol + i * (nh + 1) + j + 2,
                                                Global_symbol + (i + 1) * (nh + 1) + j + 2,
                                                Global_symbol + (i + 1) * (nh + 1) + j + 1])
                Number += 1
        # Edge Information for node set using the connection and loading
        for i in range(4):
            plateIncrease = identifier * (nw + 1) * (nh + 1)
            if (i % 2) == 0:
                EdgeIncrease = 1 if i == 0 else nw * (nh + 1) + 1
                EdgeNode = [itme for itme in range(plateIncrease + EdgeIncrease, plateIncrease + EdgeIncrease + nh + 1)]
                EdgeSet.append(EdgeNode)
            else:
                EdgeIncrease = nh + 1 if i == 1 else 1
                EdgeNode = [itme for itme in
                            range(plateIncrease + EdgeIncrease, plateIncrease + EdgeIncrease + (nh + 1) * nw + 1,
                                  nh + 1)]
                EdgeSet.append(EdgeNode)

        return NodeCoord, Element, EdgeSet

    def Inpwrite(self):
        """
        Basic Information Input Module
        The plate geometry, [width, height, Sector Angle]
        the crease length is the height, and  Sector angle between line 1/4 to 1/2, default = 90
        InitialAngle: Define the initial rest angle of the one crease
        For the crease array, there are two kinds of plates, left and right.
        """
        width = self.width
        height = self.height
        SectorAngle = self.SectorAngle
        geometry = np.array([width, height, SectorAngle])
        InitialAngle = self.InitialAngle
        PlateNumber = 2 * self.ArrayNumber

        """
        Mesh Information
        dx is the size of mesh, nw, nh are the number of seed of
        edge with the width and height.
        NodeNumber and ElementNumber are the number of the nodes and 
        elements for each plate.
        dx > 10 * thickness for shell element
        Aspect Ratio is used to design the seed distribution in the width and height direction.
        If Aspect ratio is larger than 1.0, there are larger mesh size in the Height direction. 
        """
        dx = self.dx
        Aspect_Ratio = 1.0
        nw = floor(width / dx)
        nh = floor(height / dx / Aspect_Ratio)
        NodeNumber = (nw + 1) * (nh + 1)
        ElementNumber = nw * nh

        """
        Node List and Element List Modulus
        Assembly.
        Including NodeCoord and Element Information for all plate.
        The geometry points for each plate is defined as
        1  4
        2  3
        Plate Number: 0 is for left plate, while 1 is for right plate
        Edge Number: 0 is for line 1-2, 1 is for line 2-3. 2 is for line 4-3 
                    3 is for line 1-4. The order of the line is considered.
        Position of the edge node in the mesh, Plate Number, Edge Number, Node Number 
        such as Plate Number = 0, Edge Number = 1, Node Number = 5 means the nodes is located
        at the left plate, line 2-3, the fifth node. 
        """
        NodeCoord = []
        Element = []
        EdgeSet = []
        for i in range(PlateNumber):
            Coord_temp = self._coordinate(geometry, InitialAngle, i)
            NodeCoord_temp, Element_temp, Edge_temp = self._Mesh(nw, nh, Coord_temp, i)
            if i == 0:
                NodeCoord = NodeCoord_temp
                Element = Element_temp
                EdgeSet = Edge_temp
            else:
                NodeCoord = np.concatenate((NodeCoord, NodeCoord_temp))
                Element = np.concatenate((Element, Element_temp))
                EdgeSet = EdgeSet + Edge_temp

        """
        Material Information of FEM Model
        The plate is elastic material
        Including Thickness, Elastic Module, Poisson ratio 
        unit [mm], [kg] [N] [rad]
        PET material except for the thickness
        ref: Embedded Actuation for Shape-Adaptive Origami, 2021, k=0.05
        B = Et^3/12(1-v^2)
        """
        thickness = self.thickness
        Elastic_Module = 3200.0
        Possion_ratio = 0.43
        Stiffness = self.stiffness
        CreaseStiffness = Stiffness * dx
        B = Elastic_Module * thickness ** 3 / 12 / (1 - Possion_ratio ** 2)
        ratio = Stiffness * self.width / B
        print("The scale parameter kL/B is ", ratio)

        Total_NodeNumber = PlateNumber * NodeNumber
        Total_ElementNumber = PlateNumber * ElementNumber

        """
        Loading Information
        According to the edge number in the Total system,
        Transform: global = 4 * plate + local 
        For example (1,2) corresponds to the the second plate and the third edge.
        The number of the global system is 4 * 1 + 2 = 6
        In self-folding:
        Load_Edge_left U1 U2 U3 are constrained
        Load_Edge_right U2 and U3 are constrained
        load region: 1.0(default). A number expresses the middle region is loaded. For example, 0.2 means middle 20% 
        boundary region is selected for the loading 
        """
        # loading region edge number
        Load_Edge_left = 4 * 0 + 0
        Load_Edge_right = 4 * (PlateNumber - 1) + 2

        # Loading region
        lowerbound = ceil((0.5 - self.loadregion / 2.0) * len(EdgeSet[Load_Edge_left]))
        uplowerbound = floor((0.5 + self.loadregion / 2.0) * len(EdgeSet[Load_Edge_left]))
        # Node set for load
        LoadLeft = EdgeSet[Load_Edge_left][lowerbound:uplowerbound]
        LoadRight = EdgeSet[Load_Edge_right][lowerbound:uplowerbound]

        # Loading Amplitude
        # The first step is for unfolding or folding, similar as the hinge without stiffness
        if self.foldingstate == 1:
            Dis = width * sin(SectorAngle) * (1 - sin(InitialAngle / 2)) * (PlateNumber / 2)
        elif self.foldingstate == 0:
            Dis = -width * sin(SectorAngle) * (sin(InitialAngle / 2)) * (PlateNumber / 2)
        else:
            print("error for the folding state input")

        # The second step is for Over tensioning
        Dis2 = (width * sin(SectorAngle) * (1 - sin(InitialAngle / 2)) + width * sin(SectorAngle) * 0.005) * (
                    PlateNumber / 2)

        """
        Connection Information
        Use the specified method instead of the query method
        For the two plate, the connection Information can be given 
        as (0,2) to (1,0) between the third edge of first plate and 
        the first edge of the second plate
        The local coordinate system is also defined through the node information
        There is only one connection edge
        ConectionEdge=[[2，4],[6，8],……]
        """
        ConnectionEdge = []
        ConnectionNumber = PlateNumber - 1
        for num in range(ConnectionNumber):
            PlatePair = [num, num + 1]
            EdgePair = [2, 0]
            EdgePair_global = [4 * x + y for x, y in zip(PlatePair, EdgePair)]
            ConnectionEdge.append(EdgePair_global)
        print("The Connection Pair (Edge Number) is ", ConnectionEdge)

        """
        Output Information for analysis
        Configuration Label
        OutputX : Plate0 the middle line along the width direction.
        OutputY : Plate0 the connection edge (Number: 2)
        """
        EdgePair_global = ConnectionEdge[self.ArrayNumber - 1]
        Label = EdgePair_global[0]
        OutputY = EdgeSet[Label]
        with open('OutputY.txt', "w") as f:
            for i in range(len(OutputY)):
                f.write(str(OutputY[i]) + "\n")

        OutputX = []
        for i in range(PlateNumber):
            for x in range(nw + 1):
                OutputX.append(ceil((nh + 1) / 2) + x * (nh + 1) + i * (nw + 1) * (nh + 1))
        with open('OutputX.txt', "w") as f:
            for i in range(len(OutputX)):
                f.write(str(OutputX[i]) + "\n")

        """
        Inp file Module
        Label information: Width, Height, Sector Angle, Initial Angle, PlateNumber
        Define the name of the inp file as "CreaseArray"
        Static calculation in General/static solver. 
        """

        inp_file = open(self.jobname + ".inp", "w")
        inp_file.write("*Heading\n")
        inp_file.write("**Job Name and Model Name:CreaseArray\n")
        inp_file.write("*Preprint, echo=NO, model=NO, history=NO, contact=NO\n")
        inp_file.write("**\n")

        inp_file.write("**PARTS\n")
        inp_file.write("*Part,name=Crease\n")
        inp_file.write("*Node\n")
        # for the node list
        for i in range(Total_NodeNumber):
            inp_file.write(str(i + 1) + "," + str(NodeCoord[i][1]) + ","
                           + str(NodeCoord[i][2]) + "," + str(NodeCoord[i][3]) + "\n")
        inp_file.write("**\n")

        inp_file.write("*Element,type=S4R\n")
        for i in range(Total_ElementNumber):
            inp_file.write(str(i + 1) + "," + str(Element[i][1]) + "," + str(Element[i][2]) + ","
                           + str(Element[i][3]) + "," + str(Element[i][4]) + "\n")
        inp_file.write("**\n")

        inp_file.write("*Nset, nset=SET-All, generate\n")
        inp_file.write(" 1," + str(Total_NodeNumber) + ",1\n")
        inp_file.write("*Elset, elset=SET-All, generate\n")
        inp_file.write("1," + str(Total_ElementNumber) + ",1\n")

        inp_file.write("** Section: \n")
        inp_file.write("*Shell Section, elset=SET-All, material=Self-define\n")
        inp_file.write(str(thickness) + ", 5\n")
        inp_file.write("*End Part\n")

        """
        Assembly Modulus, including the node set definition
        element set definition, and connection definition.
        """

        inp_file.write("** ASSEMBLY\n")
        inp_file.write("*Assembly, name=Assembly\n")
        inp_file.write("*Instance, name=CreaseArrayModel, part=Crease\n")
        inp_file.write("*End Instance\n")

        inp_file.write("*Nset, nset=Set-Left, instance=CreaseArrayModel\n")
        length = len(LoadLeft)
        # New Line set
        for item in range(length - 1):
            if item % 10 == 9:
                inp_file.write("\n")
            inp_file.write(str(LoadLeft[item]) + ",")
        inp_file.write(str(LoadLeft[-1]) + "\n")

        inp_file.write("*Nset, nset=Set-Right, instance=CreaseArrayModel\n")
        length = len(LoadRight)
        for item in range(length - 1):
            if item % 10 == 9:
                inp_file.write("\n")
            inp_file.write(str(LoadRight[item]) + ",")
        inp_file.write(str(LoadRight[-1]) + "\n")

        inp_file.write("*Elset, elset=SET-Plate, instance=CreaseArrayModel, generate\n")
        inp_file.write("1," + str(Total_ElementNumber) + ",1\n")

        # there are many connection Pairs in this analysis
        # Each connection pair can have different definition.
        inp_file.write("*Element, type=CONN3D2\n")
        for num in range(ConnectionNumber):
            EdgePair_global = ConnectionEdge[num]
            Pair1 = EdgePair_global[0]
            Pair2 = EdgePair_global[1]
            PairLen = len(EdgeSet[Pair1])
            for i in range(PairLen):
                if (i / PairLen < self.Creaseregion[0]) | (i / PairLen > self.Creaseregion[1]):
                    inp_file.write(str(i + (num * PairLen) + 1) + ", CreaseArrayModel." + str(
                        EdgeSet[Pair1][i]) + ", CreaseArrayModel." + str(EdgeSet[Pair2][i]) + "\n")

        for num in range(ConnectionNumber):
            EdgePair_global = ConnectionEdge[num]
            Pair1 = EdgePair_global[0]
            Pair2 = EdgePair_global[1]
            PairLen = len(EdgeSet[Pair1])
            inp_file.write("*Nset, nset=Set-Connection" + str(num + 1) + ", instance=CreaseArrayModel\n")
            for item in range(PairLen - 1):
                if (item / PairLen < self.Creaseregion[0]) | (item / PairLen > self.Creaseregion[1]):
                    if item % 6 == 5:
                        inp_file.write("\n")
                    inp_file.write(str(EdgeSet[Pair1][item]) + "," + str(EdgeSet[Pair2][item]) + ",")
            inp_file.write(str(EdgeSet[Pair1][-1]) + "," + str(EdgeSet[Pair2][-1]) + "\n")
            inp_file.write("*Elset, elset=Set-Connection" + str(num + 1) + ", generate\n")
            inp_file.write(str(num * PairLen + 1) + "," + str((num + 1) * PairLen) + ",1 \n")
            # inp_file.write("*Orientation, name=csys-Con, DEFINITION=NODES\n")
            # inp_file.write("CreaseArrayModel." + str(EdgeSet[Pair1][0]) + ", CreaseArrayModel.1, CreaseArrayModel." + str(EdgeSet[Pair1][-1]) + "\n")
            inp_file.write("*Orientation, name=csys-Con" + str(num + 1) + "\n")
            inp_file.write("0.,1.0,0.,-1.0,0.,0.\n")
            inp_file.write("1,0.\n")
            inp_file.write("*Connector Section, elset=Set-Connection" + str(num + 1) + ", behavior=ConnProp-1\n")
            inp_file.write("Hinge,\n")
            inp_file.write("csys-Con" + str(num + 1) + ",\n")
        # There are the same behavior for the hinge stiffness.
        inp_file.write("*End Assembly\n")
        inp_file.write("*Connector Behavior, name=ConnProp-1\n")
        inp_file.write("*Connector Elasticity, component=4\n")
        inp_file.write(str(CreaseStiffness) + ",\n")

        """
        Standard Module
        Including Material Information, Step Information, 
        Boundary Information
        """
        inp_file.write("** MATERIALS\n")
        inp_file.write("*Material, name=Self-define\n")
        inp_file.write("*Elastic\n")
        inp_file.write(str(Elastic_Module) + "," + str(Possion_ratio) + "\n")

        inp_file.write("** STEP: Step-1\n")
        inp_file.write("*Step, name=Step-1, nlgeom=YES, inc=1000\n")
        inp_file.write("*Static, stabilize, factor = 1e-12, allsdtol = 0, continue=NO\n")
        inp_file.write("0.01, 1., 1e-15, 0.01\n")

        inp_file.write("** BOUNDARY CONDITIONS\n")
        inp_file.write("** Name: BC-1 Type: Displacement/Rotation\n")
        inp_file.write("*Boundary\n")
        inp_file.write("Set-left, 1, 1, " + str(-Dis * 1.0) + "\n")
        inp_file.write("** Name: BC-2 Type: Displacement/Rotation\n")
        inp_file.write("*Boundary\n")
        inp_file.write("Set-left, 2, 2\n")
        inp_file.write("** Name: BC-3 Type: Displacement/Rotation\n")
        inp_file.write("*Boundary\n")
        inp_file.write("Set-left, 3, 3\n")
        inp_file.write("** Name: BC-4 Type: Displacement/Rotation\n")
        inp_file.write("*Boundary\n")
        inp_file.write("Set-right, 1, 1, " + str(Dis * 1.0) + "\n")
        inp_file.write("** Name: BC-5 Type: Displacement/Rotation\n")
        inp_file.write("*Boundary\n")
        inp_file.write("Set-right, 2, 2\n")
        inp_file.write("** Name: BC-6 Type: Displacement/Rotation\n")
        inp_file.write("*Boundary\n")
        inp_file.write("Set-right, 3, 3\n")

        inp_file.write("** CONTROLS\n")
        inp_file.write("*Controls, reset\n")
        inp_file.write("*Controls, parameters=time incrementation\n")
        inp_file.write("8, 10, , , , , , 50, , , \n")

        inp_file.write("** OUTPUT REQUESTS\n")
        inp_file.write("*Restart, write, number interval=10, time marks=YES\n")
        inp_file.write("** FIELD OUTPUT: F-Output-1\n")
        inp_file.write("*Output, field\n")
        inp_file.write("*Node Output\n")
        inp_file.write("CF, COORD, RF, U\n")
        inp_file.write("*Element Output, directions=YES\n")
        inp_file.write("LE, PE, PEEQ, PEMAG, S\n")
        inp_file.write("** FIELD OUTPUT: F-Output-2\n")
        inp_file.write("**\n")
        inp_file.write("*Output, field\n")
        inp_file.write("*Element Output, elset=Set-Connection1, directions=YES\n")
        inp_file.write("CTF,\n")
        inp_file.write("** HISTORY OUTPUT: H-Output-1\n")
        inp_file.write("*Output, history, variable=PRESELECT\n")
        inp_file.write("** HISTORY OUTPUT: H-Output-2\n")
        inp_file.write("*Output, history\n")
        for num in range(ConnectionNumber):
            inp_file.write("*Element Output, elset=Set-Connection"+str(num+1)+", directions=YES\n")
            inp_file.write("CTF,\n")
        inp_file.write("** HISTORY OUTPUT: H-Output-3\n")
        inp_file.write("*Output, history\n")
        inp_file.write("*Energy Output, elset=SET-Plate\n")
        inp_file.write("ALLSE, \n")
        inp_file.write("*End Step\n")

        # inp_file.write("** STEP: Step-2\n")
        # inp_file.write("*Step, name=Step-2, nlgeom=YES, inc=300\n")
        # inp_file.write("*Static, stabilize, factor = 1e-10, allsdtol = 0, continue=NO\n")
        # inp_file.write("0.01, 1., 1e-15, 0.01\n")
        #
        # inp_file.write("** BOUNDARY CONDITIONS\n")
        # inp_file.write("** Name: BC-1 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary, op = NEW\n")
        # inp_file.write("** Name: BC-2 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary, op = NEW\n")
        # inp_file.write("Set-left, 2, 2\n")
        # inp_file.write("** Name: BC-3 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary, op = NEW\n")
        # inp_file.write("Set-left, 3, 3\n")
        # inp_file.write("** Name: BC-4 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary, op = NEW\n")
        # inp_file.write("** Name: BC-5 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary, op = NEW\n")
        # inp_file.write("Set-right, 2, 2\n")
        # inp_file.write("** Name: BC-6 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary, op = NEW\n")
        # inp_file.write("Set-right, 3, 3\n")
        #
        #
        # inp_file.write("** CONTROLS\n")
        # inp_file.write("*Controls, reset\n")
        # inp_file.write("*Controls, parameters=time incrementation\n")
        # inp_file.write("8, 10, , , , , , 50, , , \n")
        # inp_file.write("** OUTPUT REQUESTS\n")
        # inp_file.write("*Restart, write, frequency=0\n")
        # inp_file.write("** FIELD OUTPUT: F-Output-1\n")
        # inp_file.write("*Output, field\n")
        # inp_file.write("*Node Output\n")
        # inp_file.write("CF, COORD, RF, U\n")
        # inp_file.write("*Element Output, directions=YES\n")
        # inp_file.write("LE, PE, PEEQ, PEMAG, S\n")
        # inp_file.write("** FIELD OUTPUT: F-Output-2\n")
        # inp_file.write("**\n")
        # inp_file.write("*Output, field\n")
        # for num in range(ConnectionNumber):
        #     inp_file.write("*Element Output, elset=Set-Connection"+str(num+1)+", directions=YES\n")
        #     inp_file.write("CTF,\n")
        # inp_file.write("** HISTORY OUTPUT: H-Output-1\n")
        # inp_file.write("*Output, history, variable=PRESELECT\n")
        # inp_file.write("** HISTORY OUTPUT: H-Output-2\n")
        # inp_file.write("*Output, history\n")
        # inp_file.write("*Element Output, elset=Set-Connection1\n")
        # inp_file.write("CTM1, \n")
        # inp_file.write("** HISTORY OUTPUT: H-Output-3\n")
        # inp_file.write("*Output, history\n")
        # inp_file.write("*Energy Output, elset=SET-Plate\n")
        # inp_file.write("ALLSE, \n")
        # inp_file.write("*End Step\n")
        inp_file.close()


