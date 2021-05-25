import numpy as np

def ZeroesFirst(M, n):
    for i in range(n):
        row = [0.0]*n
        M.append(row)

def ZeroesSecond(M, n, m):
    for i in range(n):
        row = [0.0]*m 
        M.append(row)

def ZeroesThird(v, n):
    for i in range(n):
        v.append(0.0)

def copyMatrix(A, copy):
    ZeroesFirst(copy,len(A))
    for i in range(len(A)):
        for j in range(len(A)):
            copy[i][j] = A[i][j]

def calculateMember(i, j, r, A, B):
    member = 0
    for k in range (r):
        member += A[i][k] * B[k][j]
    return member

def productMatrixMatrix(A, B, n, r, m):
    R = []
    ZeroesSecond(R, n , m)
    for i in range(n):
        for j in range(m):
            R[i][j] = calculateMember(i, j, r, A, B)
    return R

def productMatrixVector(A, v, R):
    for f in range(len(A)):
        cell = 0.0
        for c in range(len(v)):
            cell += A[f][c] * v[c]
    R[f] += cell

def productRealMatrix(real, M, R):
    ZeroesFirst(R, len(M))
    for i in range(len(M)):
        for j in range(len(M[0])):
            R[i][j] = real * M[i][j]

def getMinor( M, i, j):
    del M[i]
    for i in range(len(M)):
        del M[i][j]

def determinant(M):
    if len(M) == 1 : return M[0][0]
    else:
        det = 0.0
        for i in range(len(M[0])):
            minor = []
            copyMatrix(M, minor)
            getMinor(minor, 0, i)
            det +=  pow(-1, i) * M[0][i] * determinant(minor)
        return det

def cofactors(M, Cof):
    ZeroesFirst(Cof, len(M))
    for i in range(len(M)):
        for j in range(len(M[0])):
            minor = []
            copyMatrix(M, minor)
            getMinor(minor, i, j)
            Cof[i][j] = pow(-1, i+j) * determinant(minor)

def transpose(M, T):
    ZeroesSecond(T,len(M[0]), len(M))
    for i in range(len(M)):
        for j in range(len(M[0])):
            T[j][i]  = M[i][j]

def inversaMatrix(M, Minv):
    print("Iniciando calculo de inversa...\n")
    Cof = []
    Adj = []
    print("Calculo de determinante...\n")
    det = determinant(M)
    if det == 0 : exit()
    print("Iniciando calculo de cofactores...\n")
    cofactors(M,Cof)
    print("Calculo de adjunta...\n")
    transpose(Cof,Adj)
    print("Calculo de inversa...\n")
    productRealMatrix(1/det, Adj, Minv)