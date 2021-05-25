from enum import Enum
from abc import ABCMeta, abstractmethod

class Indicators(Enum):
    NOTHING = 0

class Lines(Enum):
    NOLINE = 0
    SINGLELINE = 1
    DOUBLELINE = 2

class Modes(Enum):
    NOMODE = 0
    INT_FLOAT = 1
    INT_FLOAT_FLOAT = 2
    INT_INT_INT_INT = 3

class Parameters(Enum):
    THERMAL_CONDUCTIVITY = 0
    HEAT_SOURCE = 1

class Sizes(Enum):
    NODES = 0
    ELEMENTS = 1
    DIRICHLET = 2
    NEUMANN = 3

class Item:
    def setId(self, identifier):
        self._id = identifier
    
    def setX(self, x_coord):
        self._x = x_coord

    def setY(self, y_coord):
        self._y = y_coord
    
    def setNode1(self, node_1):
        self._node1 = node_1

    def setNode2(self, node_2):
        self._node2 = node_2
    
    def setNode3(self, node_3):
        self._node3 = node_3
    
    def setValue(self, value_to_assign):
        self._value = value_to_assign
    
    def getId(self):
        return self._id
    
    def getX(self):
        return self._x
    
    def getY(self):
        return self._y
    
    def getNode1(self):
        return self._node1
    
    def getNode2(self):
        return self._node2
    
    def getNode3(self):
        return self._node3
    
    def getValue(self):
        return self._value
    
    @abstractmethod
    def setValues(self, a, b, c, d, e, f, g):
        pass

class Node(Item):
    def setValues(self, a, b, c, d, e, f, g):
        self._id = a
        self._node1 = d
        self._node2 = e
        self._node3 = f

class Condition(Item):
    def setValues(self, a, b, c, d, e, f, g):
        self._node1 = d
        self._value = g

class Mesh:
    parameters = []
    sizes = []

    def setParameters(self, k, Q):
        self.parameters.insert(Parameters.THERMAL_CONDUCTIVITY,k)
        self.parameters.insert(Parameters.HEAT_SOURCE,Q)
    
    def setSizes(self, nnodes, neltos, ndirich, nneu):
        self.sizes.insert(Sizes.NODES, nnodes)
        self.sizes.insert(Sizes.ELEMENTS, neltos)
        self.sizes.insert(Sizes.DIRICHLET, ndirich)
        self.sizes.insert(Sizes.NEUMANN, nneu)
    
    def getSize(self, s):
        return self.sizes.index(s)
    
    def getParameter(self, p):
        return self.parameters.index(p)
    
    def createData(self):
        self.node_list = []
        self.element_list = []
        self.indices_dirich = []
        self.dirichlet_list = []
        self.neuman_list = []
    
    def getNodes(self):
        return self.node_list
    
    def getElements(self):
        return self.element_list
    
    def getDirichetIndices(self):
        return self.indices_dirich
    
    def getDirichlet(self):
        return self.dirichlet_list
    
    def getNeumann(self):
        return self.neuman_list
    
    def getNode(self, i):
        return self.node_list[i]
    
    def getElement(self, i):
        return self.element_list[i]
    
    def getCondition(self, i, type):
        if type == Sizes.DIRICHLET : return self.dirichlet_list[i]
        else : return self.neuman_list[i]