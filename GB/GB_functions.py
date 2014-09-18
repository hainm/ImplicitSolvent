import numpy as np
from numpy import *
from math import radians, degrees, acos

def sq(i):
    return i*i

def dist2(a,b):
    """ Distance between atom a and b"""
    return (sq(a.x - b.x) + sq(a.y - b.y) + sq(a.z - b.z))

def dist(a,b):
    return np.sqrt(dist2(a,b))

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

def getNormedV(v):
    '''get normalized vector'''

def vlen(v):
    '''return len of vector'''
    return np.sqrt(sum(v**2))

def vdot(v1,v2):
    '''return dot product of 2 vectors'''
    return np.sum(v1*v2)

def vangle(v1,v2):
    ''' return angle between two vectors
        unit: rad
        Status: validated
    '''

def vdihedral(v1,v2,v3,v4):
    """calculate dihedral angle of four given vectors"""

if __name__ == '__main__':
    v1 = array([2.0,3.0,4.0])
    v2 = array([1.0,-2.0,3.0])
    v3 = ''
    v4 = ''

    print vlen(v1)
    print dot(v1,v2)
    cosv1v2 = vdot(v1,v2)/(vlen(v1)*vlen(v2))
    print degrees(acos(cosv1v2))

