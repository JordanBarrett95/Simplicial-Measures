import numpy
from math import comb
from collections import Counter
import matplotlib.pyplot as plt


#################
####FUNCTIONS####
#################


##
#Different options for normalizing
#1) product of number of edges
#2) Number of possible small edges in big edges
#3) Smaller of 1 and 2
##
def normalize(k,l,numk,numl,option):
    if (option==1):
        return numk*numl
    if (option==2):
        return numl*comb(l,k)
    if (option==3):
        return min(numk*numl,numl*comb(l,k))
    return 1




#################
####MAIN CODE####
#################

filenames = {'tags-math-sx',
             'tags-ask-ubuntu',
             'email-eu',
             'email-enron',
             'contact-primary-school',
             'contact-high-school',
             'hospital-lyon'
            }
for filename in filenames:

    ##
    #load data
    ##
    f = open('{0}-simple.txt'.format(filename),'r')
    file = f.read().splitlines()
    f.close()


    ##
    #extract nodes and edges
    ##
    V = set()
    E = []
    Vstart = file.index('nodes')
    Estart = file.index('edges')
    end = len(file)
    for i in range(Vstart+1,Estart):
        V.add(file[i])
    for i in range(Estart+1,end):
        E.append(eval(file[i]))
        
    n = len(V)
    m = len(E)


    ##
    #Get degrees and list of edges for each vertex
    ##
    vertexEdgeSets = {}
    deg = {}
    for v in V:
        vertexEdgeSets[v] = []
        deg[v] = 0
        
    for e in E:
        for v in e:
            vertexEdgeSets[v].append(e)
            deg[v] = deg[v]+1

    #Volume of a set of vertices
    def vol(S):
        total = 0
        for v in S:
            total = total + deg[v]
        return total

    #Declaring volume and average degree of whole graph
    L = vol(V)
    D = L/n

    ##
    #Partition edges by size
    ##
    edgeSizes = set()
    for e in E:
        edgeSizes.add(len(e))

    edgePartition = {}
    for k in edgeSizes:
        edgePartition[k] = []

    for e in E:
        edgePartition[len(e)].append(e)

    numEdges = {}
    for k in edgeSizes:
        numEdges[k] = len(edgePartition[k])

    minEdgeSize = min(edgeSizes)
    maxEdgeSize = max(edgeSizes)


    ##
    #Computing q_H(i,j) and expected_H(i,j) for each pair of edge sizes (i,j)
    ##

    #Initializing matrices
    q_H = [[0]*(maxEdgeSize+1) for i in range(maxEdgeSize+1)]
    expected_H = [[0]*(maxEdgeSize+1) for i in range(maxEdgeSize+1)]
    a_H = [[0]*(maxEdgeSize+1) for i in range(maxEdgeSize+1)]

    #Computing q_H
    for e in E:
        lene = len(e)
        minv = next(iter(e))
        for v in e:
            if (deg[v]<deg[minv]):
                minv = v
        for f in vertexEdgeSets[minv]:
            lenf = len(f)
            if (lene<lenf and e.issubset(f)):
                q_H[lene][lenf] = q_H[lene][lenf] + 1

    ###############################################################################################################################
    #Different options for computing expected_H
    #
    #Erdos-Renyi basic version
    #for k in edgeSizes:
    #    numk = numEdges[k]
    #    for l in edgeSizes:
    #        if (k<l):
    #            numl = numEdges[l]
    #            expected = comb(l,k)*numl*numk/comb(n,k)
    #            expected_H[k][l] = expected
    #
    #Loosely based on Chung-Lu model
    expected = 0
    for k in edgeSizes:
        numk = numEdges[k]
        for l in edgeSizes:
            if (k<l):
                #numl = numEdges[l]
                for e in edgePartition[l]:
                    expected = numk*((vol(e)/L)**k)
                    expected_H[k][l] = expected_H[k][l] + expected
    ###############################################################################################################################

    ###############################################################################################################################
    #Different options for computing a_H
    #
    #Based on modularity: (q_H - expected_H)/normalization
    #option = 2
    #for k in edgeSizes:
    #    for l in edgeSizes:
    #        if (k<l):
    #            a_H[k][l] = (q_H[k][l] - expected_H[k][l])/normalize(k,l,numEdges[k],numEdges[l],option)
    #
    #Based on "boosting": q_H/expected_H
    for k in edgeSizes:
        for l in edgeSizes:
            if (k<l and (expected_H[k][l] != 0)):
                a_H[k][l] = q_H[k][l]/expected_H[k][l]
    ###############################################################################################################################


    #################
    ######PLOTS######
    #################

    a_H = numpy.array(a_H)

    plt.imshow(a_H, cmap='hot', interpolation='nearest')
    plt.colorbar()
    plt.title('{0} simplicial heat-map'.format(filename))
    plt.savefig('{0}-heat.png'.format(filename))
    plt.close()

    x = []
    y = []
    for k in range(minEdgeSize,maxEdgeSize):
        x.append(k)
        numNonZero = 0
        yk = 0
        for l in edgeSizes:
            if (l>k and a_H[k][l] != 0):
                numNonZero = numNonZero + 1
                yk = yk + a_H[k][l]
        yk = yk/max(numNonZero,1)
        y.append(yk)

    plt.plot(x,y)
    plt.title('{0} bottom profile'.format(filename))
    plt.savefig('{0}-bottom.png'.format(filename))
    plt.close()

    x = []
    y = []
    for l in range(minEdgeSize+1,maxEdgeSize+1):
        x.append(l)
        numNonZero = 0
        yl = 0
        for k in edgeSizes:
            if (k<l and a_H[k][l] != 0):
                numNonZero = numNonZero + 1
                yl = yl + a_H[k][l]
        yl = yl/max(numNonZero,1)
        y.append(yl)

    plt.plot(x,y)
    plt.title('{0} top profile'.format(filename))
    plt.savefig('{0}-top.png'.format(filename))
    plt.close()


