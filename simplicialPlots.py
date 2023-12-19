import numpy
from math import comb
from collections import Counter
import matplotlib.pyplot as plt
from HypergraphFunctions import simplicialMatrix as sM, sortBySize
from HypergraphFunctions import ChungLu as CL, ChungLu2 as CL2, ChungLu3 as CL3

filenames = {'hospital-lyon',
             'NDC-substances',
             'email-eu',
             #'tags-math-sx',
             #'tags-ask-ubuntu,
             #'email-enron',
             #'contact-primary-school',
             #'contact-high-school'
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
        if (len(eval(file[i]))<=5):
            E.append(eval(file[i]))

    ##
    #Collect data
    ##
    edgePartition = sortBySize(E)
    edgeSizes = set(edgePartition.keys())
    numEdges = {}
    for i in edgeSizes:
        numEdges[i] = len(edgePartition[i])
    
    ##
    #Getting simplicial matrix
    ##
    M = sM(V,E)
    length = len(M)

    
    ##
    #Getting simplicial matrix of 20 Chung-Lu samples
    ##
    Vcl,Ecl = CL3(V,E)
    Mcl = sM(Vcl,Ecl)
    for i in range(20):
        Vcl,Ecl = CL3(V,E)
        Madd = sM(Vcl,Ecl)
        for j in range(len(Mcl)):
            for k in range(len(Mcl)):
                Mcl[j][k] = ((i+1)*Mcl[j][k] + Madd[j][k])/(i+2)

    ##
    #Method 1: (M-Mcl)/M
    ##
    M1 = [row[:] for row in M]
    for i in edgeSizes:
        for j in edgeSizes:
            if (M[i][j] != 0):
                M1[i][j] = (M[i][j]-Mcl[i][j])/M[i][j]

    ##
    #Method 2: Normalize by product of number of edges
    ##
    M2 = [row[:] for row in M]
    for i in edgeSizes:
        for j in edgeSizes:
            if (M[i][j] != 0):
                M2[i][j] = (M[i][j]-Mcl[i][j])/(numEdges[i]*numEdges[j])

    ##
    #Method 3: Normalize by max small per big
    ##
    M3 = [row[:] for row in M]
    for i in edgeSizes:
        for j in edgeSizes:
            if (M[i][j] != 0):
                M3[i][j] = (M[i][j]-Mcl[i][j])/(numEdges[j]*comb(j,i))

    ##
    #Method 4: Normalize by number of small edges
    ##
    M4 = [row[:] for row in M]
    for i in edgeSizes:
        for j in edgeSizes:
            if (M[i][j] != 0):
                M4[i][j] = (M[i][j]-Mcl[i][j])/(numEdges[i])

    ##
    #Method 5: Normalize by number of big edges
    ##
    M5 = [row[:] for row in M]
    for i in edgeSizes:
        for j in edgeSizes:
            if (M[i][j] != 0):
                M5[i][j] = (M[i][j]-Mcl[i][j])/(numEdges[j])

    ##
    #Method 6: +1 and div
    ##
    M6 = [row[:] for row in M]
    for i in edgeSizes:
        for j in edgeSizes:
            if (M[i][j] != 0):
                M6[i][j] = (M[i][j]+1)/(Mcl[i][j]+1)

    ##
    #Method 7: +epsilon and div
    ##
    eps = 2**(-10)
    M7 = [row[:] for row in M]
    for i in edgeSizes:
        for j in edgeSizes:
            if (M[i][j] != 0):
                M7[i][j] = (M[i][j]+eps)/(Mcl[i][j]+eps)

    ##
    #Method 8: div or infinity
    ##
    M8 = [row[:] for row in M]
    for i in edgeSizes:
        for j in edgeSizes:
            if (Mcl[i][j] != 0):
                M8[i][j] = M[i][j]/Mcl[i][j]
            elif (M[i][j] != 0):
                M8[i][j] = -1
  
    #################
    ######PLOTS######
    #################

    
    plt.imshow(M, cmap='cool', interpolation='nearest')
    plt.colorbar()
    for (j,i),label in numpy.ndenumerate(M):
        if (label != 0):
            plt.text(i,j,round(label,2),ha='center',va='center')
    plt.title('{0} simplicial pairs'.format(filename))
    plt.savefig('{0}-simplicial-pairs-new.png'.format(filename))
    #plt.show()
    plt.close()


    plt.imshow(Mcl, cmap='cool', interpolation='nearest')
    plt.colorbar()
    for (j,i),label in numpy.ndenumerate(Mcl):
        if (label != 0):
            plt.text(i,j,round(label,2),ha='center',va='center')
    plt.title('{0} resampled pairs'.format(filename))
    plt.savefig('{0}-resampled-pairs-new.png'.format(filename))
    #plt.show()
    plt.close()


    methods = {
        1: M1,
        2: M2,
        3: M3,
        4: M4,
        5: M5,
        6: M6,
        7: M7,
        8: M8
    }
    
    for n in methods.keys():
        plt.imshow(methods[n], cmap='cool', interpolation='nearest')
        plt.colorbar()
        for (j,i),label in numpy.ndenumerate(methods[n]):
            if (label != 0):
                plt.text(i,j,round(label,2),ha='center',va='center')
        plt.title('{0}-method {1}'.format(filename,n))
        #plt.savefig('{0}-method {1}.png'.format(filename,n))
        plt.show()
        plt.close()


    '''
    x = []
    y1 = []
    y2 = []
    m = len(M)
    for k in range(m):
        x.append(k)
        numNonZero1 = 0
        numNonZero2 = 0
        y1k = 0
        y2k = 0
        for l in range(m):
            if (l>k and M[k][l] != 0):
                if (M[k][l] != 0):
                    numNonZero1 = numNonZero1 + 1
                    y1k = y1k + M[k][l]
                if (Mcl[k][l] != 0):
                    numNonZero2 = numNonZero2 + 1
                    y2k = y2k + Mcl[k][l]
        y1k = y1k/max(numNonZero1,1)
        y2k = y2k/max(numNonZero2,1)
        y1.append(y1k)
        y2.append(y2k)

    plt.plot(x,y1,label = 'Original')
    plt.plot(x,y2,label = 'Chung-Lu')
    plt.title('{0} bottom profile'.format(filename))
    plt.legend()
    plt.savefig('{0}-bottom.png'.format(filename))
    #plt.show()
    plt.close()

    x = []
    y1 = []
    y2 = []
    for l in range(m):
        x.append(l)
        numNonZero1 = 0
        numNonZero2 = 0
        y1l = 0
        y2l = 0
        for k in range(m):
            if (k<l):
                if (M[k][l] != 0):
                    numNonZero1 = numNonZero1 + 1
                    y1l = y1l + M[k][l]
                if (Mcl[k][l] != 0):
                    numNonZero2 = numNonZero2 + 1
                    y2l = y2l + Mcl[k][l]
        y1l = y1l/max(numNonZero1,1)
        y2l = y2l/max(numNonZero2,1)
        y1.append(y1l)
        y2.append(y2l)

    plt.plot(x,y1,label = 'Origianl')
    plt.plot(x,y2,label = 'Chung-Lu')
    plt.title('{0} top profile'.format(filename))
    plt.legend()
    plt.savefig('{0}-top.png'.format(filename))
    #plt.show()
    plt.close()

    for k in range(1,4):
        x = []
        y1 = []
        y2 = []
        for l in range(k+1,m):
            x.append(l)
            y1.append(M[k][l])
            y2.append(Mcl[k][l])

        plt.plot(x,y1,label = 'Origianl')
        plt.plot(x,y2,label = 'Chung-Lu')
        plt.title('{0} {1}-in-x profile'.format(filename,k))
        plt.legend()
        plt.savefig('{0}-{1}-in-x.png'.format(filename,k))
        #plt.show()
        plt.close()

    for l in range(m-3,m):
        x = []
        y1 = []
        y2 = []
        for k in range(l):
            x.append(k)
            y1.append(M[k][l])
            y2.append(Mcl[k][l])

        plt.plot(x,y1,label = 'Origianl')
        plt.plot(x,y2,label = 'Chung-Lu')
        plt.title('{0} x-in-{1} profile'.format(filename,l))
        plt.legend()
        plt.savefig('{0}-x-in-{1}.png'.format(filename,l))
        #plt.show()
        plt.close()
    '''

