__doc__ ='''
Calculate ASA for a molecule
Learn from
https://github.com/boscoh/pdbremix/blob/master/pdbremix/asa.py
'''

import numpy as np
from math import sqrt
from math import pi as PI
from GBClass import Molecule

class ASA():
    def __init__(self):
        pass

    def get_points(self,npoints=100):
        ''' get points in spherical
            npoints = number of points in the sphere
            return narray [len = 3*npoints]

        '''
        indexarr = np.arange(npoints)
        inc = PI * (3. - sqrt(5.))
        offset = 2./npoints
        Y = indexarr*offset - 1. + (offset/2.)
        R = np.sqrt(1. - Y**2)
        phi = inc * indexarr
        X = np.cos(phi)* R
        Z = np.sin(phi) * R
        return np.array([X, Y, Z])

    def get_neighbor_list(self,atom,mol,cut=6.0):
        '''get list of neighbor atoms
            An atom j is considered be neighbor to "atom" if their distance <= cut
        '''

    def get_ASA(self,mol):
        ''' calculate ASA for mol instance of Molecule
            Return ASA (A^2)
        '''


if __name__ == '__main__':
    #status: not finish this code yet
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    asa = ASA()
    X0, Y0, Z0 = asa.get_points(npoints= 1000)
    at_radius = 1.5
    X, Y, Z = X0*at_radius, Y0*at_radius, Z0*at_radius
    top = "../Examples/Trp-cage/Tc5b.ff99SB.mb3.top"
    crd = "../Examples/Trp-cage/Tc5b.nat.crd"
    mol = Molecule(top=top,crd=crd)
    print mol.xx.__len__
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X, Y, Z, color='b')
#    plt.show()


