import sys
import argparse
from GBClass import Molecule
from CheckInput import *
from paperInfo import *

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
pinfo = PaperInfo(args.igb) #print paper's information for different GB model.


