'''calculate solvation energy using igb=1,2,5,7,8 in AMBER
usage = python this_script --igb=<igb> --top=top --crd=crd [--pdb=pdb]
'''
from CheckInput import *
from GBClass import Molecule

check = CheckInput()
check.check_igb(igb)
check.check_top_crd(top,crd)

mol = Molecule(top=top,crd=crd,pbd=pdb)

print mol.get_solv(igb=igb)


