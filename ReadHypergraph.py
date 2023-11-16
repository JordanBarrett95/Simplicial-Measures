
# Sample filenames: 
# 'hospital-lyon',
#'tags-math-sx',
#'tags-ask-ubuntu',
#'email-eu',
#'email-enron',
#'contact-primary-school',
#'contact-high-school'

def readHG(filename):
    #load data
    f = open(filename,'r')
    file = f.read().splitlines()
    f.close()
    
    #extract nodes and edges
    V = set()
    E = []
    Vstart = file.index('nodes')
    Estart = file.index('edges')
    end = len(file)
    for i in range(Vstart+1,Estart):
        V.add(file[i])
    for i in range(Estart+1,end):
        if (len(eval(file[i]))<=10):
            E.append(eval(file[i]))
        
    return(V,E)
