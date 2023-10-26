import numpy
from math import comb
from collections import Counter


##
#Input: A set (V) of vertices and a list (E) of edges
#Output: A dict (edgeSet) of edge sets and a dict (deg) of degrees
##
def vertexData(V,E):
    #Initalize edgeSet and deg
    edgeSet = {}
    deg = {}
    for v in V:
        edgeSet[v] = []
        deg[v] = 0
    #Update edgeSet and deg by iterating through E
    for e in E:
        for v in e:
            edgeSet[v].append(e)
            deg[v] = deg[v]+1
    #Return the two dicts
    return edgeSet,deg

##
#Input: A set (S) of vertices and a dict (deg) of degrees
#Output: (int) Volume of S
##
def vol(S,deg):
    total = 0
    for v in S:
        total = total + deg[v]
    return total

##
#Input: A list (E) of edges
#Output: A dict (edgePartition) sorting edges by size
##
def sortBySize(E):
    #Get edge sizes
    edgeSizes = set()
    for e in E:
        edgeSizes.add(len(e))
    #Initialize edgePartition
    edgePartition = {}
    for k in edgeSizes:
        edgePartition[k] = []
    #Update edgePartition by iterating through E
    for e in E:
        edgePartition[len(e)].append(e)
    #Return dict
    return edgePartition    

#The next group of functions are different ways of computing the expected simplicial matrix and normalizing factors.

##
#Input: Number of vertices (n) and dict (numEdges) giving the number of edges of each size
#Output: Expected matrix of containment pairs based on Edros-Renyi model
##
def expectedER(n,numEdges):
    #Initialize matrix
    edgeSizes = set(numEdges.keys())
    maxEdgeSize = max(edgeSizes)
    expected = [[0]*(maxEdgeSize+1) for i in range(maxEdgeSize+1)]
    #Iterate through all pairs of edge sizes k<l and compute the expected number of containment pairs
    for k in edgeSizes:
        if (k<maxEdgeSize):
            numk = numEdges[k]
            for l in edgeSizes:
                if (k<l):
                    numl = numEdges[l]
                    expected[k][l] = comb(l,k)*numl*numk/comb(n,k)
    #Return matrix 
    return expected

##
#Input: A dict (deg) of degrees and a dict (edgePartition) partitioning the edges by size
#Output: Expected matrix of containment pairs based on "pseudo"-Chung-Lu model
##
def expectedCL(deg,edgePartition):
    #Initialize matrix
    edgeSizes = set(edgePartition.keys())
    maxEdgeSize = max(edgeSizes)
    expected = [[0]*(maxEdgeSize+1) for i in range(maxEdgeSize+1)]
    #Compute things applied to all computations in the coming loop
    V = set(deg.keys())
    n = len(V)
    L = vol(V,deg)
    #Iterate through all pairs of edge sizes k<l and compute the expected number of containment pairs
    for k in edgeSizes:
        if (k<maxEdgeSize):
            numk = len(edgePartition[k])
            for l in edgeSizes:
                if (k<l):
                    for e in edgePartition[l]:
                        expected[k][l] = expected[k][l] + numk*((vol(e)/L)**k)
    #Return matrix
    return expected

##
#Input: A dict giving the number of edges of each size
#Output: A matrix (norm) of normalizing constants given by the product of the number of edges of each size
##
def normalizeProd(numEdges):
    #Initialize matrix
    edgeSizes = set(numEdges.keys())
    maxEdgeSize = max(edgeSizes)
    norm = [[1]*(maxEdgeSize+1) for i in range(maxEdgeSize+1)]
    #Iterate through all pairs of edge sizes k<l and compute the normalizing constant
    for k in edgeSizes:
        if (k<maxEdgeSize):
            numk = numEdges[k]
            for l in edgeSizes:
                if (k<l):
                    numl = numEdges[l]
                    norm[k][l] = numk*numl
    #Return matrix
    return norm

##
#Input: A dict giving the number of edges of each size
#Output: A matrix (norm) of normalizing constants given by fixing the number of large edges and assuming all small edges are contained
##
def normalizeMax(numEdges):
    edgeSizes = set(numEdges.keys())
    maxEdgeSize = max(edgeSizes)
    norm = [[1]*(maxEdgeSize+1) for i in range(maxEdgeSize+1)]
    #Iterate through all pairs of edge sizes k<l and compute the normalizing constant
    for k in edgeSizes:
        if (k<maxEdgeSize):
            for l in edgeSizes:
                if (k<l):
                    numl = numEdges[l]
                    norm[k][l] = numl*comb(l,k)
    #Return matrix
    return norm

##
#Input: A dict giving the number of edges of each size
#Output: A matrix (norm) of normalizing constants given by the point-wise minimum of the previous two functions
##
def normalizeMin(numEdges):
    #Initialize matrix
    edgeSizes = set(numEdges.keys())
    maxEdgeSize = max(edgeSizes)
    norm = [[1]*(maxEdgeSize+1) for i in range(maxEdgeSize+1)]
    #Iterate through all pairs of edge sizes k<l and compute the normalizing constant
    for k in edgeSizes:
        if (k<maxEdgeSize):
            numk = numEdges[k]
            for l in edgeSizes:
                if (k<l):
                    numl = numEdges[l]
                    norm[k][l] = min(numk*numl,numl*comb(l,k))
    #Return matrix
    return norm

##
#Input: A set (vertices) of vertices and a list (edges) of edges
#Output: A matrix M measuring the simpliciality of each k<l pair of edge sizes
##
def simplicialMatrix(vertices,edges):
    #Organizing data
    V = vertices
    E = edges
    n = len(V)
    m = len(E)
    edgeSet,deg = vertexData(V,E)
    edgePartition = sortBySize(E)
    edgeSizes = set(edgePartition.keys())
    maxEdgeSize = max(edgeSizes)
    numEdges = {}
    for k in edgeSizes:
        numEdges[k] = len(edgePartition[k])
    #Initialize matrix
    M = [[0]*(maxEdgeSize+1) for i in range(maxEdgeSize+1)]
    #Iterate through all pairs of edge sizes k<l and count containment pairs
    for e in E:
        if (len(e)<maxEdgeSize):
            lene = len(e)
            #For each edge e, to count edges f containing e, we need only consider edges containing v for any vertex v in e.
            #So we find v in e that minimizes deg[v] and check all f in edgeSet[v].
            minv = next(iter(e))
            for v in e:
                if (deg[v]<deg[minv]):
                    minv = v
            for f in edgeSet[minv]:
                lenf = len(f)
                if (lene<lenf and e.issubset(f)):
                    M[lene][lenf] = M[lene][lenf] + 1
    #At this point, M is the matrix of num containment pairs.
    #Now we apply normalization.
                    
    #Options for expected number of containment pairs
    #################################################
    expected = expectedER(n,numEdges)
    #expected = expectedCL(deg,edgePartition)
    #################################################

    #Options for normalization
    #################################################
    #norm = normalizeProd(numEdges)
    norm = normalizeMax(numEdges)
    #norm = normalizeMin(numEdges)
    #################################################
    
    for k in edgeSizes:
        if (k<maxEdgeSize):
            for l in edgeSizes:
                if (k<l):
                    #Options for applying expected and norm
                    #################################################
                    M[k][l] = (M[k][l]-expected[k][l])/norm[k][l]
                    #if (expected[k][l] != 0):
                    #    M[k][l] = M[k][l]/expected[k][l]
                    #################################################
    #Return matrix
    return M
  