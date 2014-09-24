"""Force field class"""
from GB.GBClass import Molecule

class ForceFieldForMolecule(Molecule):
    def __init__(self):
        pass
    def getBond(self):
        pass
    def getAngle(self):
        pass
    def getDih(self):
        pass
    def setBondK(self):
        pass
    def setDihK(self):
        pass
    def setDihPhi(self):
        pass
    def setDihPsi(self):
        pass
    def setDihPhiPrime(self):
        pass
    def setDihPsiPrime(self):
        pass
    def get_EBOND(self):
        pass
    def get_EANGLE(self):
        pass
    def get_EDIHED(self):
        pass
    def get_E_1_4NB(self):
        pass
    def get_E1_4EEL(self):
        pass
    def get_EVDWAALS(self):
        pass
    def get_EEELEC(self):
        pass
    def get_EGB(self):
        '''same as Molecule.get_solv()'''
        pass
    def get_ERESTRAINT(self):
        pass
