from numpy import *

def vdinv(i,v):
    for k in range(i):
       v[k] = 1/float(v[k])
    return v

def vdln(i,v):
    for k in range(i):
        v[k] = math.log(v[k])
    return v

def vdsqrtinv(i,v):
    for k in range(i):
        v[k] = 1/math.sqrt(v[k])
    return v
# end vector sections

def sq(i):
    return i*i

def inv(i):
    try:
        return 1/i
    except:
        pass

def dist2(a,b):
    """ Distance between atom a and b"""
    return (sq(a.x - b.x) + sq(a.y - b.y) + sq(a.z - b.z))

def dist(a,b):
    return math.sqrt(dist2(a,b))

def fgb(a,b):
    f = dist2(a,b) + a.R * b.R * exp(-dist2(a,b)/(4 * a.R * b.R))
    return sqrt(f)

def linenum(natoms, linesize):
    ''' use this function to parse top and crd files for molecule method.
        linesize is number of items of line.split()
    '''
    if (natoms % linesize) == 0:
        return int(natoms / linesize)
    else:
        return int(natoms / linesize + 1)

