import math_tools as mt
import classes as c
from math import sqrt

def showMatrix(K):
    for i in range(len(K[0])):
        print("[\t")
        for j in range(len(K)):
            print(str(K[i][j]) + "\t")
        print("]\n")

def showKs(Ks):
    for i in range(len(Ks)):
        print("K del elemento " + str(i + 1) + ":\n")
        showMatrix(Ks[i])
        print("*************************************\n")

def showVector(b):
    print("[\t")
    for i in range(len(b)):
        print(str(b[i]) + "\t")
    print("]\t")

def showbs(bs):
    for i in range(len(bs)):
        print("b del elemento " + str(i+1) + ":\n")
        showVector(bs[i])
        print("*************************************\n")

def calculateLocalD(i, m):
    e = m.getElement(i)
    n1 = m.getNode(e.getNode1()-1)
    n2 = m.getNode(e.getNode2()-1)
    n3 = m.getNode(e.getNode3()-1)

    a = n2.getX() - n1.getX()
    b = n2.getY() - n1.getY()
    c = n3.getX() - n1.getX()
    d = n3.getY() - n1.getY()

    D = a*d - b*c

    return D

def calculateMagnitude(v1, v2):
    return sqrt(pow(v1,v2) + pow(v2,2))

def calculateLocalArea(i, m):
    e = m.getElement(i)
    n1 = m.getNode(e.getNode1()-1)
    n2 = m.getNode(e.getNode2()-1)
    n3 = m.getNode(e.getNode3()-1)
    
    a = calculateMagnitude(n2.getX() - n1.getX(), n2.getY() - n1.getY())
    b = calculateMagnitude(n3.getX() - n2.getX(), n3.getY() - n2.getY())
    c = calculateMagnitude(n3.getX() - n1.getX(), n3.getY() - n1.getY())

    s = (a+b+c)/2

    A = sqrt(s*(s-a)*(s-b)*(s-c))
    return A

def calculateLocalA(i, A, m):
    e = m.getElement(i)
    n1 = m.getNode(e.getNode1() - 1)
    n2 = m.getNode(e.getNode2() - 1)
    n3 = m.getNode(e.getNode3() - 1)

    A[0][0] = n3.getY() - n1.getY()
    A[0][1] = n1.getY() - n2.getY()
    A[1][0] = n1.getX() - n3.getX()
    A[1][1] = n2.getX() - n1.getX()

def calculateB(B):
    B[0][0] = -1
    B[0][1] = 1
    B[0][2] = 0
    B[1][0] = -1
    B[1][1] = 0
    B[1][2] = 1

def createLocalK(element, m):
    k = m.getParameter(c.Parameters.THERMAL_CONDUCTIVITY)
    A = []
    B = []
    K = []
    Bt = []
    At = []

    D = calculateLocalD(element, m)
    Ae = calculateLocalArea(element, m)

    mt.ZeroesFirst(A,2)
    mt.ZeroesSecond(B,2,3)
    calculateLocalA(element, A, m)
    calculateB(B)
    mt.transpose(A, At)
    mt.transpose(B, Bt)

    mt.productRealMatrix(k*Ae/(D*D),mt.productMatrixMatrix(Bt,mt.productMatrixMatrix(At,mt.productMatrixMatrix(A,B,2,2,3),2,2,3),3,2,3),K)

    return K

def calculateLocalJ(i, m):
    e = m.getElement(i)
    n1 = m.getNode(e.getNode1() - 1)
    n2 = m.getNode(e.getNode2() - 1)
    n3 = m.getNode(e.getNode3() - 1)

    a = n2.getX() - n1.getX()
    b = n3.getX() - n1.getX()
    c = n2.getY() - n1.getY()
    d = n3.getY() - n1.getY()

    J = a*d -b*c

    return J

def createLocalb(element, m):
    b = []

    q = m.getParameter(c.Parameters.HEAT_SOURCE)
    J = calculateLocalJ(element, m)

    b_i = ( q * J ) / 6

    b.append(b_i)
    b.append(b_i)
    b.append(b_i)

    return b

def crearSistemasLocales(m, localKs, localbs):
    for i in range(m.getSize(c.Sizes.ELEMENTS)):
        localKs.append(createLocalK(i, m))
        localbs.append(createLocalb(i, m))

def assemblyK(e, localK, K):
    index1 = e.getNode1() - 1
    index2 = e.getNode2() - 1
    index3 = e.getNode3() - 1

    K[index1][index1] += localK[0][0]
    K[index1][index2] += localK[0][1]
    K[index1][index3] += localK[0][2]
    K[index2][index1] += localK[1][0]
    K[index2][index2] += localK[1][1]
    K[index2][index3] += localK[1][2]
    K[index3][index1] += localK[2][0]
    K[index3][index2] += localK[2][1]
    K[index3][index3] += localK[3][1]

def assemblyb(e, localb, b):
    index1 = e.getNode1() - 1
    index2 = e.getNode2() - 1
    index3 = e.getNode3() - 1

    b[index1] += localb[0]
    b[index2] += localb[1]
    b[index3] += localb[2]

def ensamblaje(m, localKs, localbs, K, b):
    for i in range(m.getSize(c.Sizes.ELEMENTS)):
        e = m.getElement(i)
        assemblyK(e, localKs[i], K)
        assemblyb(e, localbs[i], b)

def applyNeumann(m, b):
    for i in range(m.getSize(c.Sizes.NEUMANN)):
        n = m.getCondition(i, c.Sizes.NEUMANN)
        b[n.getNode1() - 1] += n.getValue()

def applyDirichlet(m, K, b):
    for i in range(m.getSize(c.Sizes.DIRICHLET)):
        d = m.getCondition(i, c.Sizes.DIRICHLET)
        index = d.getNode1() - 1

        del K[index]
        del b[index]

        for row in range(len(K)):
            cell = K[row][index]
            del K[row][index]
            b[row] += -1*c.getValue()*cell 

def calculate(K, b, T):
    print("Iniciando calculo de respuesta...\n")
    Kinv = []
    print("Calculo de inversa...\n")
    mt.inversaMatrix(K, Kinv)
    print("Calculo de respuesta...\n")
    mt.productMatrixVector(Kinv, b, T)

