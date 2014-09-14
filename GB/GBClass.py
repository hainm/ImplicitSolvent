import sys
import re
from GB_functions import *

FSMAX  = 5.0
RGBMAX = 25.0
AMBCONST = 18.222 #amber constant to convert amber charge to normal charge

class Atom:

    def __init__(self,x=0,y=0,z=0,r=0,R=0,q=0, atype='',offset=0.09,ro=0,rs=0,sf=0):
        """ x,y,z: atom coordinates
        r, R: bondi and effective radius
        q: partial charge
        offset: offset
        ro: ro = r - offest
        rs = r * sf
        sf: Scaling Factor (Sh, Sc, Sn, Ss, So)
        """
        self.x = x
        self.y = y
        self.z = z
        self.r = r
        self.R = R
        self.q = q
        self.sf = sf
        self.ro = ro
        self.rs = rs
        self.offset = offset
        self.atype = atype

    def setCoord(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z

    def setEffRad(self,R):
        self.R = R

    def setTypeChargeRadiusFs(self,atype,charge,radius,sf):
        self.atype = atype
        self.q = charge
        self.r = radius
        self.sf = sf

    def setRoRs(self):
        self.ro = self.r - self.offset
        self.rs = self.ro * self.sf #change self.r to self.ro

    def showCoord(self):
        print self.x, self.y, self.z

    def showAll(self):
        print 'type: ', (self.atype),' r= ',self.r,' ro= ',self.ro,' sf= ',self.sf, ' rs= ', self.rs, \
                            'coor = (',self.x, self.y, self.z,')',' R= ', self.R

class Molecule:

    def __init__(self,ein=1.0,eout=78.5,offset=0.09,crd='',top='',pdb='',molType=''):
        '''GBNeck: scalingDict = {'H':1.091, 'C':0.484, 'N':0.700, 'O':1.066}
         GBOBC: scalingDict = {}
         self.fss: array of atom.rs
         self.rborn: array of atom.r
         self.Rbonrs contain is array of self.rborn
         self.effRad_PB: array of "perfect radii" from PB caclulation
         sefl.selfPB: PB array of self term
         self.atypes: array of atom types
            type: string, len = 2
         self.mod_atypes = array(self.atypes,'c')
         molType: molecule type, if mol get from top, molType = top
         self.typekeys: array of atom type key
         self.AMBCONST: constant to convert amber charge to normal charge
       '''
        #initialize data
        self.atoms = []
        self.xx = []
        self.fss = []
        self.rborn = []
        self.Rborns = []
        self.atypes = []
        self.atname_pdb = []
        self.mod_atypes = []
        self.natoms = 0
        self.charges = []
        self.effRads = []
        self.ein = ein
        self.eout = eout
        self.effRad_PB = []
        self.selfEnPB = []
        self.offset = offset
        self.typekeys = []
        self.AMBCONST = 18.222
        self.reslb = []
        self.reslblen4 = []
        self.respt = []
        if top and crd and not pdb:
            self.getAtoms(top,crd)
        elif top and not crd:
            self.getAtomsFromTop(top)
            if pdb:
                self.getCoordFromPDB(pdb)


    def getAtomsFromTop(self,top):
        ''' get atoms from top only (without coordinate file)'''
        fh = open(top)
        lines = fh.readlines() # put all lines of topology file in lines array

        atypes = []
        atname_pdb = []
        charges = []
        radii = []
        fss = []
        respt = []
        reslb = []

        # get atom.q, atom.r, atom.atype, screening parameters
        for line in lines:
            # get atom number
            if line.startswith('%FLAG POINTERS'):
                index = lines.index(line)
                natoms12 = lines[index+2].split()[0]
                natoms = int(natoms12)

            #get atom charge (amber format)
            if line.startswith('%FLAG CHARGE'):
                index = lines.index(line)
                numline = linenum(natoms,5)
                for i in range(index+2,index+numline+2):
                    for charge in lines[i].split():
                        charges.append(float(charge)) # You could use 18.2222 to convert charge in AMBER format to normal charge

            # get atom radii
            if line.startswith('%FLAG RADII'):
                index = lines.index(line)
                numline = linenum(natoms,5)
                for i in range(index+2,index+numline+2):
                    for radius in lines[i].split():
                        radii.append(float(radius))

            # get atom type
            if line.startswith('%FLAG AMBER_ATOM_TYPE'):
                index = lines.index(line)
                numline = linenum(natoms,20)
                for i in range(index+2,index+numline+2):
                    for atype in lines[i].split():
                        if len(atype) == 1:
                            #add space to atype so len(atype) = 2
                            # for example: 'H' should be 'H '
                            atype = atype + ' '
                        atypes.append(atype)

            # get screening parameter to sf
            if line.startswith('%FLAG SCREEN'):
                index = lines.index(line)
                numline = linenum(natoms,5)
                for i in range(index+2,index+numline+2):
                    for s in lines[i].split():
                        fss.append(float(s))

            #get FLAG RESIDUE_LABEL:
            if line.startswith('%FLAG RESIDUE_LABEL'):
                i = lines.index(line) + 2
                while not lines[i].startswith('%FLAG'):
                    reslb = reslb + [st for st in lines[i].split()]
                    i = i +1

            if line.startswith('%FLAG ATOM_NAME'):
                i = lines.index(line) + 2
                while not lines[i].startswith('%FLAG'):
                    #avoid this: HG11HG12HG13CG2, HG21HG22HG23CG1
                    for st in lines[i].split():
                        if len(st) in [5,6,7,8]:
                            atname_pdb.append(st[:4])
                            atname_pdb.append(st[4:])
                        elif len(st) in [9,10,11,12]:
                            atname_pdb.append(st[:4])
                            atname_pdb.append(st[4:8])
                            atname_pdb.append(st[8:])
                        elif len(st) > 12:
                            atname_pdb.append(st[:4])
                            atname_pdb.append(st[4:8])
                            atname_pdb.append(st[8:11])
                            atname_pdb.append(st[12:])
                        else:
                            atname_pdb.append(st)
                    i = i+1
                for i,st in enumerate(atname_pdb):
                    if len(st) == 1:
                        atname_pdb[i] +=  '   '
                    elif len(st) == 2:
                        atname_pdb[i] += '  '
                    elif len(st) == 3:
                        atname_pdb[i] += ' '

            #get FLAG RESIDUE_POINTER
            if line.startswith('%FLAG RESIDUE_POINTER'):
                i = lines.index(line) + 2
                while not lines[i].startswith('%FLAG'):
                    respt = respt + [int(st) for st in lines[i].split()]
                    i = i +1

        # set atom types, charges, radii for atoms of molecule and add atoms to array stored in molecule.atoms
        for i in range(len(charges)):
            atom = Atom()
            atom.setTypeChargeRadiusFs(atypes[i],charges[i],radii[i],fss[i])
            atom.setRoRs()
            fss[i] = atom.rs
            self.atoms.append(atom)
        fh.close()

        self.molType = top
        self.fss = fss
        self.rborn = radii
        self.atypes = atypes #work with numpy array when using f2py
        self.atname_pdb = atname_pdb
        #self.mod_atypes = array(atypes,'c')
        self.natoms = natoms
        self.charges = charges
        self.typekeys = self.getinfo().keys()
        self.reslb = reslb
        self.respt = respt

    def getAtoms(self,top,crd):
        ''' Set: atom.x, atom.y, atom.z, atom.q, atom.r, atom.atype,
        atom.offset, atom.ro, atom.rs for molecule.
        This is an example of crd file:
        109
              -1.4290652   5.2696260   0.4109236  -0.6377563   5.5230198  -0.2945710
               0.3292836   5.0308183  -0.1912677  -0.4089177   6.5594278  -0.5427559
              -1.0704247   4.9698918  -1.6515413  -2.2444098   4.6086778  -1.7336542
                ...(skip)
      '''
        self.getAtomsFromTop(top)
        xx = []
        # get atom coordinate from restart or crd file
        fh = open(crd)
        lines = fh.readlines()
        try:
            natoms = int(lines[0])
        except:
            lines.pop(0)
            natoms = int(lines[0])
        if natoms != self.natoms:
            print "Are use sure you use the correct coordinate?"
            sys.exit()

        numline = linenum(natoms,2)
        index = 0
        for i in range(1,numline+1):
            line = lines[i].split()
            self.atoms[index].setCoord(float(line[0]),float(line[1]),float(line[2]))
            xx.append(float(line[0]))
            xx.append(float(line[1]))
            xx.append(float(line[2]))
            try:
                self.atoms[index+1].setCoord(float(line[3]),float(line[4]),float(line[5]))
                xx.append(float(line[3]))
                xx.append(float(line[4]))
                xx.append(float(line[5]))
            except:
                pass
            index = index + 2
        fh.close() # close crd file
        # while using module from FORTRAN, I could not use numpy to pass these array to function? why?
        self.xx = xx # coordinate file: the same with the amber file. Should I use numpy.array or just use normal array?

    def getBlock(self,index,crd):
        '''get block array from crd file
        '''
        pass

    def printInfo(self):
        '''
        '''

    def getCoordFromTraj(self,top,crd):
        """get coordinate from top and long trajectories file
        Return array contain coordinate for each snapshot.
        NOT DONE YET
        """
        index = 0
        while 1:
            index += 1
            block  = self.getBlock(index,crd)
            self.Rborns.append(block)
        pass

    def getCoordFromPDB(self,pdb):
        '''get coordinate from pdb file generated by ptraj
        For the safety, these pdb files must be genereated by ptraj with "nowrap"
        "REMARK 1 PDB file generated by ptraj (set  1296)"
        Look at this folder for reference:
        PDB generated by ptraj: /mnt/raidc/haichit/GB/PB/ala10/clusterAnalysis/pdbFromCluster
        PDB generated by vmd (from top and crd files): /mnt/raidc/haichit/GB/PB/trp_nmr/trp_nmr_vmd.pdb
        NOTE: Please check coordinates of atoms from different pdb file format
        VMD format:
            CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1

            ATOM      1  N   MET  X   1      27.340  24.430   2.614  1.00  0.00

        ~/scripts/pdbProcess.sh format: (no header)
            ATOM      1  N   ASN  A   1      -8.901   4.127  -0.555  0.00  0.00   N

        ptraj format:
            REMARK   1                     PDB file generated by ptraj (set     1)
            ATOM      1 N    PRO     1      50.614  59.061  -8.731  0.00  0.00

        ambpdb format: the same as ptraj but make sure to remove "TER" near final "END" line
        '''
        lines = open(pdb,'r').readlines()
        if lines[0].startswith("CRYST1"):    #generated from vmd.
            try:
                #there is no problem with coordinate (see Molecule.write2NewPdb())
                #skip header "CRYST1" in PDB file generated by vmd, remove "END" line
                arrpdb = loadtxt(pdb,skiprows=1,comments="END",usecols=(6,7,8))
            except:                          #re-write pdb file
                self.write2NewPdb(pdb)
                arrpdb = loadtxt(pdb+".mod",skiprows=1,comments="END",usecols=(6,7,8))
        elif lines[0].startswith("REMARK"):   #generated from ptraj
#            print 'test'
            try:
                #no problem with coordinate (see Molecule.write2NewPdb())
                #skip header "REMARK" in PDB file generated by ptraj, remove "END" line
                arrpdb = loadtxt(pdb,skiprows=1,comments="END",usecols=(5,6,7))
#                print arrpdb.__len__()
            except:                           #re-write pdb file
                self.write2NewPdb(pdb)
                arrpdb = loadtxt(pdb+".mod",skiprows=1,comments="END",usecols=(5,6,7))
        else:
            #no header: generated by ambpdb and then removing header
            #pdb file generated by tleap have coordinates in colum 6,7,8 (index: 5,6,7)
            try:
                #there is no problem with coordinate (see Molecule.write2NewPdb())
                #remove "END" line, do not skip first line. Make sure to remove TER
                arrpdb = loadtxt(pdb,skiprows=0,comments="END",usecols=(5,6,7))
            except: #re-write pdb file
                self.write2NewPdb(pdb)
                arrpdb = loadtxt(pdb+".mod",skiprows=0,comments="END",usecols=(5,6,7))

        for i,arrtmp in enumerate(arrpdb):
            atomi = self.atoms[i]
            atomi.x,atomi.y,atomi.z = arrtmp
        self.xx = [coord for coord in ravel(arrpdb)]

    def updateFs(self,scalingDict,offset=0.09):
        '''when using GBNeck, you have to change default scaling parameter set
         (in topology file). This function will update self.fss
        scalingDict = { ... }
        '''
        self.offset = offset
        sf5 = [scalingDict[at] for at in ["H","C","N","O","S"]]
        self.fss = fss.fss.fsscal(sf5,self.rborn,self.atypes,self.offset)

    def updateFs_sf(self,sf):
        '''update Fss from scaling factor array
        '''
        self.fss = fss.fss.fsscal(sf,mol.rborn,mol.mod_atypes,mol.offset)

    def updatebondi(self,bkey):
        '''bkey = 'bondi','mbondi','mbondi2',
        TODO: NOT DONE YET
        '''

        if bkey == 'bondi':
            for i,atom in enumerate(self.atoms):
                if atom.atype == 'H': self.rborn[i] = ''
        pass

    def getGBArray_3Pars(self, gbArray):
        '''gbArray = [alpha, beta, gamma]
        '''
        alpha, beta, gamma = gbArray
        gbvalpha,gbvbeta,gbvgamma = [],[],[]
        for atom in self.atoms:
            gbvalpha.append(alpha)
            gbvbeta.append(beta)
            gbvgamma.append(gamma)
        return gbvalpha,gbvbeta,gbvgamma

    def getGBArray_HO(self,gbArray):
        ''' gbArray = [gbalpha,gbbeta,gbgamma,
                        gbalphaH,gbbetaH,gbgammaH,
                        gbalphaO,gbbetaO,gbgammaO]
        '''
        gbvalpha,gbvbeta,gbvgamma = [],[],[]
        for atom in self.atoms:
            if atom.atype == 'H':
                gbvalpha.append(gbArray[3])
                gbvbeta.append(gbArray[4])
                gbvgamma.append(gbArray[5])
            elif atom.atype == 'O':
                gbvalpha.append(gbArray[6])
                gbvbeta.append(gbArray[7])
                gbvgamma.append(gbArray[8])
            else:
                gbvalpha.append(gbArray[0])
                gbvbeta.append(gbArray[1])
                gbvgamma.append(gbArray[2])
        return gbvalpha,gbvbeta,gbvgamma

    def getGBArray_NC(self,gbArray):
        ''' gbArray = [gbalpha,gbbeta,gbgamma,
                        gbalphaN,gbbetaN,gbgammaN,
                        gbalphaC,gbbetaC,gbgammaC]
            N is amide nitrogen (backbone)
            C is carboxyl carbon (backbone)
        '''
        gbvalpha,gbvbeta,gbvgamma = [],[],[]
        for atom in self.atoms:
            if atom.atype == 'N':
                gbvalpha.append(gbArray[3])
                gbvbeta.append(gbArray[4])
                gbvgamma.append(gbArray[5])
            elif atom.atype == 'C':
                gbvalpha.append(gbArray[6])
                gbvbeta.append(gbArray[7])
                gbvgamma.append(gbArray[8])
            else:
                gbvalpha.append(gbArray[0])
                gbvbeta.append(gbArray[1])
                gbvgamma.append(gbArray[2])
        return gbvalpha,gbvbeta,gbvgamma

    def getGBArray_NC_4pars(self,gbArray):
        '''getGBArray_NC_4pars: use 4 parameters
        exp: tanh((gbalpha - gbbeta*psi + gbgamma*psi**2 + gblambda*psi**3)*psi)
         gbArray = [gbalpha,gbbeta,gbgamma,gblambda
                        gbalphaN,gbbetaN,gbgammaN,gblambdaN
                        gbalphaC,gbbetaC,gbgammaC,gblambdaC]
            N is amide nitrogen (backbone)
            C is carboxyl carbon (backbone)
        use:
        '''
        gbvalpha,gbvbeta,gbvgamma,gbvlambda = [],[],[],[]
        for atom in self.atoms:
            if atom.atype == 'N':
                gbvalpha.append(gbArray[4])
                gbvbeta.append(gbArray[5])
                gbvgamma.append(gbArray[6])
                gbvlambda.append(gbArray[7])
            elif atom.atype == 'C':
                gbvalpha.append(gbArray[8])
                gbvbeta.append(gbArray[9])
                gbvgamma.append(gbArray[10])
                gbvlambda.append(gbArray[11])
            else:
                gbvalpha.append(gbArray[0])
                gbvbeta.append(gbArray[1])
                gbvgamma.append(gbArray[2])
                gbvlambda.append(gbArray[3])
        return gbvalpha,gbvbeta,gbvgamma,gbvlambda

    def getGBArray_NCHO(self,gbArray):
        ''' For BackBone atoms
            gbArray = [gbalpha,gbbeta,gbgamma,
                        gbalphaN,gbbetaN,gbgammaN,
                        gbalphaC,gbbetaC,gbgammaC,
                        gbalphaH,gbbetaH,gbgammaH,
                        gbalphaO,gbbetaO,gbgammaO]
            H,N is amide hydrogen and nitrogen (backbone)
            C,O is carboxyl carbon  and Oxygen(backbone)
        '''
        gbvalpha,gbvbeta,gbvgamma = [],[],[]
        for atom in self.atoms:
            if atom.atype == 'N':
                gbvalpha.append(gbArray[3])
                gbvbeta.append(gbArray[4])
                gbvgamma.append(gbArray[5])
            elif atom.atype == 'C':
                gbvalpha.append(gbArray[6])
                gbvbeta.append(gbArray[7])
                gbvgamma.append(gbArray[8])
            elif atom.atype == 'H':
                gbvalpha.append(gbArray[9])
                gbvbeta.append(gbArray[10])
                gbvgamma.append(gbArray[11])
            elif atom.atype == 'O':
                gbvalpha.append(gbArray[12])
                gbvbeta.append(gbArray[13])
                gbvgamma.append(gbArray[14])
            else:
                gbvalpha.append(gbArray[0])
                gbvbeta.append(gbArray[1])
                gbvgamma.append(gbArray[2])
        return gbvalpha,gbvbeta,gbvgamma

    def getGBArray_HCNO_allatom(self,gbArray):
        '''For all atom
           gbArray =   [gbalphaH,gbbetaH,gbgammaH,
                        gbalphaC,gbbetaC,gbgammaC,
                        gbalphaN,gbbetaN,gbgammaN,
                        gbalphaO,gbbetaO,gbgammaO]
           return gbalphav,gbbetav,gbgammav VECTOR
           (each atom has it own (gbalpha,gbbeta,gbgamma) set)
        '''
        gbvalpha,gbvbeta,gbvgamma = [],[],[]
        for atom in self.atoms:
            if atom.atype[0] == 'H':
                gbvalpha.append(gbArray[0])
                gbvbeta.append(gbArray[1])
                gbvgamma.append(gbArray[2])
            elif atom.atype[0] == 'C':
                gbvalpha.append(gbArray[3])
                gbvbeta.append(gbArray[4])
                gbvgamma.append(gbArray[5])
            elif atom.atype[0] == 'N':
                gbvalpha.append(gbArray[6])
                gbvbeta.append(gbArray[7])
                gbvgamma.append(gbArray[8])
            elif atom.atype[0] == 'O' or atom.atype[0] == 'S':
                gbvalpha.append(gbArray[9])
                gbvbeta.append(gbArray[10])
                gbvgamma.append(gbArray[11])
        return gbvalpha,gbvbeta,gbvgamma

    def getGBArray_HCNOP_allatom(self,gbArray):
        '''For all atom
           gbArray =   [gbalphaH,gbbetaH,gbgammaH,
                        gbalphaC,gbbetaC,gbgammaC,
                        gbalphaN,gbbetaN,gbgammaN,
                        gbalphaO,gbbetaO,gbgammaO,
                        gbalphaP,gbbetaP,gbgammaP,
                        gbalphaOp,gbbetaOp,gbgammaOp]
           return gbalphav,gbbetav,gbgammav VECTOR
           (each atom has it own (gbalpha,gbbeta,gbgamma) set)
        '''
        gbvalpha,gbvbeta,gbvgamma = [],[],[]
        for atom in self.atoms:
            if atom.atype[0] == 'H':
                gbvalpha.append(gbArray[0])
                gbvbeta.append(gbArray[1])
                gbvgamma.append(gbArray[2])
            elif atom.atype[0] == 'C':
                gbvalpha.append(gbArray[3])
                gbvbeta.append(gbArray[4])
                gbvgamma.append(gbArray[5])
            elif atom.atype[0] == 'N':
                gbvalpha.append(gbArray[6])
                gbvbeta.append(gbArray[7])
                gbvgamma.append(gbArray[8])
            elif atom.atype[0] == 'O' or atom.atype[0] == 'S':
                gbvalpha.append(gbArray[9])
                gbvbeta.append(gbArray[10])
                gbvgamma.append(gbArray[11])
            elif atom.atype[0] == 'P':
                gbvalpha.append(gbArray[12])
                gbvbeta.append(gbArray[13])
                gbvgamma.append(gbArray[14])
        return gbvalpha,gbvbeta,gbvgamma

    def getGBArray_HCNOPOp_allatom(self,gbArray):
        '''STATUS: NOT DONE YET
           For all atom
           gbArray =   [gbalphaH,gbbetaH,gbgammaH,
                        gbalphaC,gbbetaC,gbgammaC,
                        gbalphaN,gbbetaN,gbgammaN,
                        gbalphaO,gbbetaO,gbgammaO,
                        gbalphaP,gbbetaP,gbgammaP,
                        gbalphaOp,gbbetaOp,gbgammaOp]
           return gbalphav,gbbetav,gbgammav VECTOR
           (each atom has it own (gbalpha,gbbeta,gbgamma) set)
        '''
        gbvalpha,gbvbeta,gbvgamma = [],[],[]
        for atom in self.atoms:
            if atom.atype[0] == 'H':
                gbvalpha.append(gbArray[0])
                gbvbeta.append(gbArray[1])
                gbvgamma.append(gbArray[2])
            elif atom.atype[0] == 'C':
                gbvalpha.append(gbArray[3])
                gbvbeta.append(gbArray[4])
                gbvgamma.append(gbArray[5])
            elif atom.atype[0] == 'N':
                gbvalpha.append(gbArray[6])
                gbvbeta.append(gbArray[7])
                gbvgamma.append(gbArray[8])
            elif atom.atype[0] == 'S':
                gbvalpha.append(gbArray[9])
                gbvbeta.append(gbArray[10])
                gbvgamma.append(gbArray[11])
            elif atom.atype[0] == 'P':
                gbvalpha.append(gbArray[12])
                gbvbeta.append(gbArray[13])
                gbvgamma.append(gbArray[14])
            elif atom.atype[0] == 'O':
                at = atom.atype[:2]
                if at == 'O2' or at == 'OH' or at == 'OS':
                    #Op atoms
                    gbvalpha.append(gbArray[15])
                    gbvbeta.append(gbArray[16])
                    gbvgamma.append(gbArray[17])
                else:
                    #not Op atoms
                    gbvalpha.append(gbArray[9])
                    gbvbeta.append(gbArray[10])
                    gbvgamma.append(gbArray[11])
        return gbvalpha,gbvbeta,gbvgamma

    def get_offset(self,offset_pro,offset_P):
        '''Add desctiption here
        '''
        N = len(self.atoms)
        offset = zeros(N)
        for i in range(N):
            atom = self.atoms[i]
            if atom.atype[0] == 'O':
                at = atom.atype[:2]
                if at == 'O2' or at == 'OH' or at == 'OS':
                   offset[i] = offset_P
                else:
                    offset[i] = offset_pro
            elif atom.atype[0] == 'P':
                offset[i] = offset_P
            else:
                offset[i] = offset_pro
        return offset

    def getGBArray_NCHO_4pars(self,gbArray):
        ''' gbArray = [gbalpha,gbbeta,gbgamma,gblambda
                        gbalphaN,gbbetaN,gbgammaN,gblambdaN,
                        gbalphaC,gbbetaC,gbgammaC,gblambdaC,
                        gbalphaH,gbbetaH,gbgammaH,gblambdaH,
                        gbalphaO,gbbetaO,gbgammaO,gblambdaO]
            H,N is amide hydrogen and nitrogen (backbone)
            C,O is carboxyl carbon  and Oxygen(backbone)
        '''
        gbvalpha,gbvbeta,gbvgamma,gbvlambda = [],[],[],[]
        for atom in self.atoms:
            if atom.atype == 'N':
                gbvalpha.append(gbArray[4])
                gbvbeta.append(gbArray[5])
                gbvgamma.append(gbArray[6])
                gbvlambda.append(gbArray[7])
            elif atom.atype == 'C':
                gbvalpha.append(gbArray[8])
                gbvbeta.append(gbArray[8])
                gbvgamma.append(gbArray[10])
                gbvlambda.append(gbArray[11])
            elif atom.atype == 'H':
                gbvalpha.append(gbArray[12])
                gbvbeta.append(gbArray[13])
                gbvgamma.append(gbArray[14])
                gbvlambda.append(gbArray[15])
            elif atom.atype == 'O':
                gbvalpha.append(gbArray[16])
                gbvbeta.append(gbArray[17])
                gbvgamma.append(gbArray[18])
                gbvlambda.append(gbArray[19])
            else:
                gbvalpha.append(gbArray[0])
                gbvbeta.append(gbArray[1])
                gbvgamma.append(gbArray[2])
                gbvlambda.append(gbArray[3])
        return gbvalpha,gbvbeta,gbvgamma,gbvlambda

    def getGBArray_NCHO_4pars_Hall(self,gbArray):
        ''' gbArray = [gbalpha,gbbeta,gbgamma,gblambda
                        gbalphaN,gbbetaN,gbgammaN,gblambdaN,
                        gbalphaC,gbbetaC,gbgammaC,gblambdaC,
                        gbalphaH,gbbetaH,gbgammaH,gblambdaH,
                        gbalphaO,gbbetaO,gbgammaO,gblambdaO]
            N is amide nitrogen (backbone)
            C,O is carboxyl carbon  and Oxygen(backbone)
            H is ALL TYPE OF H (HO,HC,H1,HA...)
        '''
        gbvalpha,gbvbeta,gbvgamma,gbvlambda = [],[],[],[]
        for atom in self.atoms:
            if atom.atype == 'N':
                gbvalpha.append(gbArray[4])
                gbvbeta.append(gbArray[5])
                gbvgamma.append(gbArray[6])
                gbvlambda.append(gbArray[7])
            elif atom.atype == 'C':
                gbvalpha.append(gbArray[8])
                gbvbeta.append(gbArray[8])
                gbvgamma.append(gbArray[10])
                gbvlambda.append(gbArray[11])
            elif atom.atype[0] == 'H':
                gbvalpha.append(gbArray[12])
                gbvbeta.append(gbArray[13])
                gbvgamma.append(gbArray[14])
                gbvlambda.append(gbArray[15])
            elif atom.atype == 'O':
                gbvalpha.append(gbArray[16])
                gbvbeta.append(gbArray[17])
                gbvgamma.append(gbArray[18])
                gbvlambda.append(gbArray[19])
            else:
                gbvalpha.append(gbArray[0])
                gbvbeta.append(gbArray[1])
                gbvgamma.append(gbArray[2])
                gbvlambda.append(gbArray[3])
        return gbvalpha,gbvbeta,gbvgamma,gbvlambda

    def getGBArray_NCHCT_4pars_Hall(self,gbArray):
        ''' gbArray = [gbalpha,gbbeta,gbgamma,gblambda
                        gbalphaN,gbbetaN,gbgammaN,gblambdaN,
                        gbalphaC,gbbetaC,gbgammaC,gblambdaC,
                        gbalphaH,gbbetaH,gbgammaH,gblambdaH,
                        gbalphaO,gbbetaO,gbgammaO,gblambdaO]
            N is amide nitrogen (backbone)
            C is carboxyl carbon  and Oxygen(backbone)
            H is ALL TYPE OF H (HO,HC,H1,HA...)
            CT is carbon with CT type
        '''
        gbvalpha,gbvbeta,gbvgamma,gbvlambda = [],[],[],[]
        for atom in self.atoms:
            if atom.atype == 'N':
                gbvalpha.append(gbArray[4])
                gbvbeta.append(gbArray[5])
                gbvgamma.append(gbArray[6])
                gbvlambda.append(gbArray[7])
            elif atom.atype == 'C':
                gbvalpha.append(gbArray[8])
                gbvbeta.append(gbArray[8])
                gbvgamma.append(gbArray[10])
                gbvlambda.append(gbArray[11])
            elif atom.atype[0] == 'H':
                gbvalpha.append(gbArray[12])
                gbvbeta.append(gbArray[13])
                gbvgamma.append(gbArray[14])
                gbvlambda.append(gbArray[15])
            elif atom.atype == 'CT':
                gbvalpha.append(gbArray[16])
                gbvbeta.append(gbArray[17])
                gbvgamma.append(gbArray[18])
                gbvlambda.append(gbArray[19])
            else:
                gbvalpha.append(gbArray[0])
                gbvbeta.append(gbArray[1])
                gbvgamma.append(gbArray[2])
                gbvlambda.append(gbArray[3])
        return gbvalpha,gbvbeta,gbvgamma,gbvlambda

    def getGBArray_alltypes_4pars(self,gbArray):
        ''' gbArray = [ gbalpha,gbbeta,gbgamma,gblambda   #0-3
                        gbalphaN,gbbetaN,gbgammaN,gblambdaN, #4-7
                        gbalphaC,gbbetaC,gbgammaC,gblambdaC, #8-11
                        gbalphaH,gbbetaH,gbgammaH,gblambdaH, #12-15
                        gbalphaCT,gbbetaCT,gbgammaCT,gblambdaCT, #16-19
                        gbalphaH1,gbbetaH1,gbgammaH1,gblambdaH1, #20-23
                        gbalphaHC,gbbetaHC,gbgammaHC,gblambdaHC, #24-27
                        gbalphaHA,gbbetaHA,gbgammaHA,gblambdaHA, #28-31
                        gbalphaCA,gbbetaCA,gbgammaCA,gblambdaCA, #32-35
                        gbalphaO,gbbetaO,gbgammaO,gblambdaO]     #36-39

            N is amide nitrogen (backbone)
            C is carboxyl carbon  and Oxygen(backbone)
            H is ALL TYPE OF H (HO,HC,H1,HA...)
            CT is carbon with CT type
        '''
        gbvalpha,gbvbeta,gbvgamma,gbvlambda = [],[],[],[]
        for atom in self.atoms:
            if atom.atype == 'N':
                gbvalpha.append(gbArray[4])
                gbvbeta.append(gbArray[5])
                gbvgamma.append(gbArray[6])
                gbvlambda.append(gbArray[7])
            elif atom.atype == 'C':
                gbvalpha.append(gbArray[8])
                gbvbeta.append(gbArray[9])
                gbvgamma.append(gbArray[10])
                gbvlambda.append(gbArray[11])
            elif atom.atype[0] == 'H':
                gbvalpha.append(gbArray[12])
                gbvbeta.append(gbArray[13])
                gbvgamma.append(gbArray[14])
                gbvlambda.append(gbArray[15])
            elif atom.atype == 'CT':
                gbvalpha.append(gbArray[16])
                gbvbeta.append(gbArray[17])
                gbvgamma.append(gbArray[18])
                gbvlambda.append(gbArray[19])
            elif atom.atype == 'H1':
                gbvalpha.append(gbArray[20])
                gbvbeta.append(gbArray[21])
                gbvgamma.append(gbArray[22])
                gbvlambda.append(gbArray[23])
            elif atom.atype == 'HC':
                gbvalpha.append(gbArray[24])
                gbvbeta.append(gbArray[25])
                gbvgamma.append(gbArray[26])
                gbvlambda.append(gbArray[27])
            elif atom.atype == 'HA':
                gbvalpha.append(gbArray[28])
                gbvbeta.append(gbArray[29])
                gbvgamma.append(gbArray[30])
                gbvlambda.append(gbArray[31])
            elif atom.atype == 'CA':
                gbvalpha.append(gbArray[32])
                gbvbeta.append(gbArray[33])
                gbvgamma.append(gbArray[34])
                gbvlambda.append(gbArray[35])
            elif atom.atype == 'O':
                gbvalpha.append(gbArray[36])
                gbvbeta.append(gbArray[37])
                gbvgamma.append(gbArray[38])
                gbvlambda.append(gbArray[39])
            else:
                gbvalpha.append(gbArray[0])
                gbvbeta.append(gbArray[1])
                gbvgamma.append(gbArray[2])
                gbvlambda.append(gbArray[3])

        return gbvalpha,gbvbeta,gbvgamma,gbvlambda

    def getEffRadFromReff(self,reff):
        '''Get effective radii from reff array
       '''
        self.effRads = reff
        for i in range(len(reff)):
            self.atoms[i].R = reff[i]

    def updateRad_PB(self):
        '''update perfect radii'''
        self.getEffRadFromReff(self.effRad_PB)

    def getDeltaG(self):
        ''' Calculate Total, Self, Cross Solvation Energy'''
        sumTot = 0
        sumSelf = 0
        sumCross = 0
        #reff = self.effRads
        #sumSelf = sum(self.charges**2/reff)
        # atom.R == 0?, if yes, atom.R = 0.0000000000001 to avoid 0/0
        for i in range(self.natoms):
            if self.atoms[i].R == 0: self.atoms[i].R = 0.001

        for i in range(self.natoms):
            atomi = self.atoms[i]
            if atomi.R == 0: self.atoms[i]
            else: sumSelf = sumSelf + sq(atomi.q)/atomi.R

        for i in range(self.natoms):
            for j in range(i+1,len(self.atoms)):
                atomi = self.atoms[i]
                atomj = self.atoms[j]
                fgbtmp = fgb(atomi,atomj)
                sumCross = sumCross + atomi.q * atomj.q / fgbtmp
        sumCross = -1.0 * (1/self.ein - 1/self.eout) *  sumCross # 1.0 and 0.5 are BIG different. If you use 1.0 here, you would get the correct answers.
        sumSelf  = -0.5 * (1/self.ein - 1/self.eout) *  sumSelf
        sumTot   =  sumCross + sumSelf
        return (sumTot, sumSelf, sumCross)

    def getDeltaGFromFortran(self):
        '''using FORTRAN function'''
        #return deltaG.deltag.dltg(self.xx,self.effRads,self.rborn,self.charges,self.ein,self.eout)
        return deltaG_testNozero.deltag.dltg(self.xx,self.effRads,self.rborn,self.charges,self.ein,self.eout)

    def getDeltaGFromFortran_for_GivenEffRadSet(self,reff):
        '''use this function for optimztion to save time'''
        return deltaG_testNozero.deltag.dltg(self.xx,reff,self.rborn,self.charges,self.ein,self.eout)

    def getDeltaGFromFortranFull(self,igb,egbVersion):
        """given igb and egbVersion, gb enery will be calculated"""
        if igb == 5:
            scalingDict = {'H':0.85, 'C':0.72, 'N':0.79, 'O':0.85, 'S':0.96}
            gbArray =  [1.0,0.8,4.85]
        elif igb == 7:
            scalingDict = {'H':1.09085413633, 'C':0.484353823306,
                            'N':0.700147318409, 'O':1.06557401132,'S': 0.602256336067}
            gbArray = [1.09511284,1.90792938,2.50798245]
        gbalpha, gbbeta, gbgamma = gbArray
        egb = egbVersion.effrad.egb_calc_radii
#        t0 = time()
        self.updateFs(scalingDict)
#        t01 = time()
        reff = egb(igb,self.xx,self.fss,FSMAX,RGBMAX,self.rborn,offset,
                   gbalpha,gbbeta,gbgamma,gbneckscale)[0]
        return deltaG.deltag.dltg(self.xx,reff,self.rborn,self.charges,self.ein,self.eout)
#        t1 = time()
#        self.getEffRadFromReff(reff)
#        print self.getDeltaGFromFortran()
#        t2 = time()
#        print self.getDeltaGFromFortran_for_GivenEffRadSet(reff)
##        t3 = time()
##        print  t01-t0,t2-t1,t3-t2

    def getSelfEn_perfectR_PB(self,pbfile):
        '''get self energy term from PB calculation'''
        from numpy import loadtxt as lt
        import sys

        selfPB = lt(pbfile)
        if self.natoms != len(selfPB):
            print "PB and GB do not match the number of atoms. Exit"
            sys.exit()

        temp = - 0.5 * (1/self.ein - 1/self.eout)
        self.selfEnPB = selfPB  # which term is correct?

        for i, ene in enumerate(selfPB):
            if ene == 0: selfPB[i] = 0.00000001 # avoid 0/0
        self.effRad_PB = temp*pow(array(self.charges),2)/self.selfEnPB

    def getSelfEne(self,igb,egbVersion):
        ''' get self energy array
        '''
        selfEne = []
        if igb == 5:
            scalingDict = {'H':0.85, 'C':0.72, 'N':0.79, 'O':0.85, 'S':0.96}
            gbArray =  [1.0,0.8,4.85]
        elif igb == 7:
            scalingDict = {'H':1.09085413633, 'C':0.484353823306,
                            'N':0.700147318409, 'O':1.06557401132,'S': 0.602256336067}
            gbArray = [1.09511284,1.90792938,2.50798245]
        self.updateFs(scalingDict)
        gbalpha, gbbeta, gbgamma = gbArray
        egb = egbVersion.effrad.egb_calc_radii
        reff = egb(igb,self.xx,self.fss,FSMAX,RGBMAX,self.rborn,offset,
                   gbalpha,gbbeta,gbgamma,gbneckscale)[0]
        self.getEffRadFromReff(reff)
        tmp = -0.5 * (1/self.ein - 1/self.eout)
        for i in range(len(self.atoms)):
            atomi = self.atoms[i]
            temp = tmp * sq(atomi.q)/atomi.R
            selfEne.append(temp)
        return array(selfEne)

    def getSelfEne_NCHO(self,gbArray_NCHO,egbVersion):
        """Use gbArray_NCHO. Target: modified GBOBC
        """
        selfEne = []
        gbvalpha,gbvbeta,gbvgamma = self.getGBArray_NCHO(gbArray_NCHO)
        egb = egbVersion.effrad.egb_calc_radii
        reff = egb(self.xx,self.fss,FSMAX,RGBMAX,self.rborn,offset,gbvalpha,gbvbeta,gbvgamma)[0]
        self.getEffRadFromReff(reff)
        tmp = -0.5 * (1/self.ein - 1/self.eout)
        for i in range(len(self.atoms)):
            atomi = self.atoms[i]
            temp = tmp * sq(atomi.q)/atomi.R
            selfEne.append(temp)
        return array(selfEne)

    def _converString(self,errStr):
        ''' convert this -94.376-136.032 to -94.376  -136.032 '''
        import re
        #
        errStr = re.split("\d-",errStr)
        return errStr[0], "-"+errStr[-1]

    def write2NewPdb(self,pdb):
        ''' When pdb file has this kind of error: -94.376-136.032
            --> you need to seperate them so that numpy.loadtxt could work properly
        '''
        fh = open(pdb,'r')
        lines = fh.readlines()
        for i,line in enumerate(lines):
            line = re.sub("\d-","  -",line)
            #avoid this: -94.376-136.032-94.376
            lines[i] = re.sub("\d-","  -",line)
        fh.close()
        open(pdb+".mod",'w').writelines(lines)

    def getinfo(self):
        '''print my info
        AtomNum = total number of atom type
        '''
        dict = {}
        #creat dictionary of atoms, initiate atom number.
        N = len(self.atoms)
        for i in range(N):
            key = self.atname_pdb[i]
            atom = self.atoms[i]
            if not dict.has_key(key):
                dict[key] = {"Amber type":key,"Radii":atom.r,"Charge" :atom.q/AMBCONST,"AtomNum":0,
                                   }
#                print dict[atom.atype]
        for i in range(N):
            key = self.atname_pdb[i]
            dict[key]["AtomNum"] = dict[key]["AtomNum"] + 1
        return dict

    def typeExtract(self,atype,reff_gb,reff_pb,psi):
        '''plot psi vs. reff for specific atom type
        mol = Molecule
        '''
        n = self.natoms
        reff_tmp_gb,reff_tmp_pb,psi_tmp = [],[],[]
        if (len(reff_gb) != n) or (len(reff_pb) != n) or (len(psi) != n):
            print "Please check you array lenght. Thanks."
            sys.exit()
        for i in range(self.natoms):
            atom = self.atoms[i]
            if atom.atype == atype:
                reff_tmp_gb.append(reff_gb[i])
                reff_tmp_pb.append(reff_pb[i])
                psi_tmp.append(psi[i])
        return array(reff_tmp_gb),array(reff_tmp_pb),array(psi_tmp)

    def typeExtractAll(self,reff_gb,reff_pb,psi):
        '''mol = Molecule
        Extract all atom types info.
        '''
        #Creating dictionary of reff_tmp_gb,reff_tmp_pb,psi_tmp
        #example: {"H":[],"C":[],...}
        reff_tmp_gb,reff_tmp_pb,psi_tmp = {},{},{}
        for at in self.typekeys:
            reff_tmp_gb[at],reff_tmp_pb[at],psi_tmp[at] = [],[],[]
        for i,atom in enumerate(self.atoms):
            for at in self.typekeys:
                if atom.atype == at:
                    reff_tmp_gb[at].append(reff_gb[i])
                    reff_tmp_pb[at].append(reff_pb[i])
                    psi_tmp[at].append(psi[i])
        for at in self.typekeys:
            reff_tmp_gb[at] = array(reff_tmp_gb[at])
            reff_tmp_pb[at] = array(reff_tmp_pb[at])
            psi_tmp[at] = array(psi_tmp[at])
        return reff_tmp_gb,reff_tmp_pb,psi_tmp

    def typeExtractAll_with_keys(self,reff_gb,reff_pb,atkeys):
        '''mol = Molecule
        Extract all atom types info.
        VALIDATED: YES
        '''
        #Creating dictionary of reff_tmp_gb,reff_tmp_pb,psi_tmp
        #example: {"H   ":[],"C    ":[],...}
        reff_tmp_gb,reff_tmp_pb = {},{}
        for at in atkeys:
            reff_tmp_gb[at],reff_tmp_pb[at] = [],[]
        if 'all' in [at.lower() for at in atkeys]:
            #all atom
            reff_tmp_gb['All'] = reff_gb
            reff_tmp_pb['All'] = reff_pb
        for i,atom in enumerate(self.atoms):
            for at in atkeys:
                atarr = []
                if at == 'BB':
                    atarr = ['CA  ','C   ','N   ','O   ']
                elif at == 'BBH':
                    atarr = ['CA  ','C   ','N   ','O   ','H   ']
                else:
                    atarr = [at]
#                print atarr
#                print self.atname_pdb[i]
                if self.atname_pdb[i] in atarr:
                    reff_tmp_gb[at].append(reff_gb[i])
                    reff_tmp_pb[at].append(reff_pb[i])
#                print self.atname_pdb[i]
        for at in atkeys:
            reff_tmp_gb[at] = array(reff_tmp_gb[at])
            reff_tmp_pb[at] = array(reff_tmp_pb[at])
        return reff_tmp_gb,reff_tmp_pb

    def showme(self):
        print "atom type:", self.typekeys
        for at in self.typekeys:
            print self.getinfo()[at]
        print "Residue label: ",self.reslb
        print "Residue pointer: ",self.respt

    def createAtomNameString(self):
        '''get string of atom name such as HHHCNOS...
        '''
        str = ' '.join([atom.atype[0] for atom in self.atoms])
        return str

    def createAtomNameString_2chars(self):
        '''get string of atom name such as HHHCNOS...
        '''
        str = ' '.join([atom.atype for atom in self.atoms])
        return str

    def createString_from_array(self,arr):
        '''create string of an array --> use join?
        '''
        return ' '.join([str(a) for a in arr])

    def show_res_info(self):
        '''show residue info: NOT DONE YET
        ind_1,ind_2 are begining, end indeces of atoms.
        '''
        N = len(self.respt)
        for i in range(N):
            if i == N - 1:
                ind_1,ind_2 = self.respt[i], self.natoms
            else:
                ind_1,ind_2 = self.respt[i], self.respt[i+1]
            a = [at for at in self.atname_pdb[ind_1:ind_2]]
            print "Residue %s: " %(self.reslb[i]), a

    def is_backbone(self,at):
        '''Check atom is backbone or not
        at = atomname, len(at) = 4
        '''
        if at in ['CA  ','C   ','N   ','O   ']:
            return True
        else:
            return False

    def is_backbone_addH(self,at):
        '''Check atom is backbone or not
        at = atomname, len(at) = 4
        '''
        if at in ['CA  ','C   ','N   ','O   ','H   ']: return True
        else: return False

    def get_back_bone_based_arr(self,reff):
        '''return child-array of reff for each residue, using only backbone atoms (no H)
        '''
        N = len(self.reslb)
        res = []
        #add back bone atom info to res[i]
        for i in range(N):
            if i == N - 1:
                ind_1,ind_2 = self.respt[i], self.natoms
            else:
                ind_1,ind_2 = self.respt[i], self.respt[i+1]
            res.append([reff[ind] for ind in range(ind_1,ind_2) if self.is_backbone(self.atname_pdb[ind])])
        return res

    def get_back_bone_based_arr_addH(self,reff):
        '''return child-array of reff for each residue, using only backbone atoms (add amide H)
        '''
        N = len(self.reslb)
        res = []
        #add back bone atom info to res[i]
        for i in range(N):
            if i == N - 1:
                ind_1,ind_2 = self.respt[i], self.natoms
            else:
                ind_1,ind_2 = self.respt[i], self.respt[i+1]
            res.append([reff[ind] for ind in range(ind_1,ind_2) if self.is_backbone_addH(self.atname_pdb[ind])])
        return res

    def write_coord(self,data,fname):
        '''Write data as Amber coordinate
        atom_number
        x1 y1 z1 x2 y2 z2
        ...
        '''
        fh = open(fname,'w')
        natom = len(data)/3
        fh.write(str(natom)+"\n")
        for i in range(natom/2):
            start = 6*i
            end   = 6*(i+1)
            st    = '  '.join([str(a) for a in data[start:end]]) + '\n'
            fh.write(st)
        fh.close()

    def update_rborn(self,fname,type='siz'):
        '''update rborn from a file'''
        if type == 'siz':
            #topology file is converted to *crg and *siz by Dan's program
            rnew = loadtxt(fname,usecols=(2,),skiprows=1)
            self.rborn = rnew

#Functions
def create_GB_vector(atypes,gbArray):
    '''For all atom
       gbArray =   [gbalphaH,gbbetaH,gbgammaH,
                    gbalphaC,gbbetaC,gbgammaC,
                    gbalphaN,gbbetaN,gbgammaN,
                    gbalphaO,gbbetaO,gbgammaO]
        return vectors: gbvalpha,gbvbeta,gbvgamma
    '''
    gbvalpha,gbvbeta,gbvgamma = [],[],[]
    for atype in atypes:
        if atype[0] == 'H':
            gbvalpha.append(gbArray[0])
            gbvbeta.append(gbArray[1])
            gbvgamma.append(gbArray[2])
        elif atype[0] == 'C':
            gbvalpha.append(gbArray[3])
            gbvbeta.append(gbArray[4])
            gbvgamma.append(gbArray[5])
        elif atype[0] == 'N':
            gbvalpha.append(gbArray[6])
            gbvbeta.append(gbArray[7])
            gbvgamma.append(gbArray[8])
        elif atype[0] == 'O' or atom.atype[0] == 'S':
            gbvalpha.append(gbArray[9])
            gbvbeta.append(gbArray[10])
            gbvgamma.append(gbArray[11])
    return gbvalpha,gbvbeta,gbvgamma

def test_gbEne(top,pdb,igb,egbVersion,sf=''):
    ''' Test GB energy
    '''
    mol = Molecule(top=top,pdb=pdb)
    if not sf:
        if igb in [1,5]:
            scalingDict = {'H':0.85, 'C':0.72, 'N':0.79, 'O':0.85, 'S':0.96}
            gbArray =  [1.0,0.8,4.85]
            mol.updateFs(scalingDict)
            gbneckscale = 0.0
        elif igb == 7:
            scalingDict = {'H':1.09085413633, 'C':0.484353823306,
                                'N':0.700147318409, 'O':1.06557401132,'S': 0.602256336067}
            gbArray = [1.09511284,1.90792938,2.50798245]
            mol.updateFs(scalingDict)
            gbneckscale = 0.361825
        gbalpha, gbbeta, gbgamma = gbArray
    if sf:
        scalingDict = {'H':sf[0], 'C':sf[1], 'N':sf[2], 'O':sf[3], 'S':sf[4]}
        mol.offset = sf[5]
        mol.updateFs(scalingDict,offset=sf[5])
        gbalpha,gbbeta,gbgamma = sf[6:9]
        gbneckscale = sf[-1]
    egb = egbVersion.effrad.egb_calc_radii
    reff = egb(igb,mol.xx,mol.fss,FSMAX,RGBMAX,mol.rborn,mol.offset,
               gbalpha,gbbeta,gbgamma,gbneckscale)[0]
    print deltaG.deltag.dltg(mol.xx,reff,mol.rborn,mol.charges,mol.ein,mol.eout)

def getEffRad_1(top,igb,egbVersion,pdb="",crd=""):
    """ get reff array from top,crd,gbArray and egbVersion"""
    if pdb:
        mol = Molecule(top=top,pdb=pdb)
    elif crd:
        mol = Molecule(top=top,crd=crd)
    if igb in  [1,5]:
        scalingDict = {'H':0.85, 'C':0.72, 'N':0.79, 'O':0.85, 'S':0.96}
        gbArray =  [1.0,0.8,4.85]
        mol.updateFs(scalingDict)
    elif igb == 7:
        scalingDict = {'H':1.09085413633, 'C':0.484353823306,
                    'N':0.700147318409, 'O':1.06557401132,'S': 0.602256336067}
        gbArray = [1.09511284,1.90792938,2.50798245]
        mol.updateFs(scalingDict)
    gbalpha, gbbeta, gbgamma = gbArray
    egb = egbVersion.effrad.egb_calc_radii
    reff = egb(igb,mol.xx,mol.fss,FSMAX,RGBMAX,mol.rborn,offset,
                           gbalpha,gbbeta,gbgamma,gbneckscale)[0]
    return reff

if __name__ == '__main__':

    top,crd = "trp_nmr.top", "trp_nmr.crd"
    mol = Molecule(top=top, crd=crd)
    N = len(mol.atypes)
    for i in range(N):
        print mol.atypes[i], mol.atname_pdb[i]
