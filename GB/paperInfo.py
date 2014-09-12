class PaperInfo():
    HELPIGB = '''igb = 1, 2, 5, 7, or 8'''

    DESCRITION_IGB='''calculate solvation energy using igb=1,2,5,7,8 in AMBER
    '''

    GBHCT = '''Hawkins, G. D.; Cramer, C. J.; Truhlar, D. G., Pairwise solute descreening of solute charges from a dielectric medium. Chem. Phys. Lett. 1995, 246 (1-2), 122-129.'''
    GBOBC = '''Onufriev, A.; Bashford, D.; Case, D. A., Exploring protein native states and large-scale conformational changes with a modified generalized born model.Proteins: Struct., Funct., Bioinf. 2004, 55 (2), 383-394.'''
    GBNeck = '''Mongan, J.; Simmerling, C.; McCammon, J. A.; Case, D. A.; Onufriev, A., Generalized Born Model with a Simple, Robust Molecular Volume Correction. J. Chem. Theory Comput. 2007, 3 (1), 156-169.'''
    GBNeck2_GBNecknu2='''
    1. Protein
    Nguyen, H.; Roe, D. R.; Simmerling, C., Improved Generalized Born Solvent Model Parameters for Protein Simulations. J. Chem. Theory Comput. 2013, 9 (4), 2020-2034.

    2. Nucleic acid
    Nguyen, H.; Perez, A.; Simmerling, C., Refinement of Generalized-Born Neck Parameters for Nucleic Acid and Their Complex with Protein. (in preperation).
    '''
    def __init__(self,igb=''):
        self.igb = igb
        if self.igb:
            self.printInfo(self.igb)
    def printInfo(self,igb):
        if igb == 1:
            print self.GBHCT
        elif igb in [2,5]:
            print self.GBOBC
        elif igb == 7:
            print self.GBNeck
        elif igb == 8:
            print self.GBNeck2_GBNecknu2

if __name__ == '__main__':
    info = PaperInfo(8)
#    info.printInfo(8)
