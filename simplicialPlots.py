import numpy
from math import comb
from collections import Counter
import matplotlib.pyplot as plt
from HypergraphFunctions import simplicialMatrix as sM


filenames = {'hospital-lyon',
             #'tags-math-sx',
             #'tags-ask-ubuntu',
             #'email-eu',
             'email-enron',
             'contact-primary-school',
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
        E.append(eval(file[i]))
        
    ##
    #Getting simplicial matrix
    ##
    M = sM(V,E)
        
    
    #################
    ######PLOTS######
    #################

    plt.imshow(M, cmap='hot', interpolation='nearest')
    plt.colorbar()
    plt.title('{0} simplicial heat-map'.format(filename))
    #plt.savefig('{0}-heat.png'.format(filename))
    plt.show()
    plt.close()

    x = []
    y = []
    m = len(M)
    for k in range(m):
        x.append(k)
        numNonZero = 0
        yk = 0
        for l in range(m):
            if (l>k and M[k][l] != 0):
                numNonZero = numNonZero + 1
                yk = yk + M[k][l]
        yk = yk/max(numNonZero,1)
        y.append(yk)

    plt.plot(x,y)
    plt.title('{0} bottom profile'.format(filename))
    #plt.savefig('{0}-bottom.png'.format(filename))
    plt.show()
    plt.close()

    x = []
    y = []
    for l in range(m):
        x.append(l)
        numNonZero = 0
        yl = 0
        for k in range(m):
            if (k<l and M[k][l] != 0):
                numNonZero = numNonZero + 1
                yl = yl + M[k][l]
        yl = yl/max(numNonZero,1)
        y.append(yl)

    plt.plot(x,y)
    plt.title('{0} top profile'.format(filename))
    #plt.savefig('{0}-top.png'.format(filename))
    plt.show()
    plt.close()

