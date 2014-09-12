'''Hai Nguyen
I made this work for my GB project. You're free to use for whatever purpose.
It would be great if you can give me credit for my work.
    1. Protein
    Nguyen, H.; Roe, D. R.; Simmerling, C., Improved Generalized Born Solvent Model Parameters for Protein Simulations. J. Chem. Theory Comput. 2013, 9 (4), 2020-2034.

    2. Nucleic acid
    Nguyen, H.; Perez, A.; Simmerling, C., Refinement of Generalized-Born Neck Parameters for Nucleic Acid and Their Complex with Protein. (in preperation).
'''

import sys
import argparse
from GB.GBClass import Molecule
from GB.CheckInput import *
from GB.paperInfo import *

pinfo = PaperInfo()
DESCRITION_IGB = pinfo.DESCRITION_IGB
HELPIGB = pinfo.HELPIGB

parser = argparse.ArgumentParser(description=DESCRITION_IGB)
parser.add_argument('--top',help='AMBER topology file',dest='top')
parser.add_argument('--crd',help='CRD file',dest='crd')
parser.add_argument('--pdb',help='PDB file, preferred ptraj format',dest='pdb')
parser.add_argument('--igb',action="store",type=int,help=HELPIGB,dest='igb')
args = parser.parse_args()

if not (args.igb and args.top):
    print "you must provide igb and top file\n. Use '--help' for further info"
    sys.exit()
if not (args.crd or args.pdb):
    print "You must provide either CRD or PDB files"
    sys.exit()

check = CheckInput()
check.check_igb(igb)
check.check_top_crd(top,crd)

#initialize molecule to get access solvation energy calculating method
mol = Molecule(top=top,crd=crd,pbd=pdb)
print mol.get_solv(igb=igb)
 #print paper's information for different GB model.
 pinfo.printInfo(args.igb)


