#!/usr/bin/env python
#STATUS: not validated yet

DESCRITION = """convert new topology file: changing mbondi2 to mbondi3 radii"""

'''Hai Nguyen
I made this work for my GB project. You're free to use for whatever purpose.
It would be great if you can give me credit for my work.
    1. Protein
    Nguyen, H.; Roe, D. R.; Simmerling, C., Improved Generalized Born Solvent Model Parameters for Protein Simulations. J. Chem. Theory Comput. 2013, 9 (4), 2020-2034.
'''
import sys
import argparse
from GB.Utilities import Utilities
from GB.paperInfo import *

parser = argparse.ArgumentParser(description=DESCRITION)
parser.add_argument('--oldtop',help='Old top file with mbondi2 radii',dest='topold')
parser.add_argument('--newtop',help='New top file with mbondi3 radii',dest='topnew')
args = parser.parse_args()
oldtop = args.topold
newtop = args.topnew

if not (args.topold and args.topnew):
    print "need old and new top file. Use 'python %s --help' for help" %(sys.argv[0])
    sys.exit(1)

util = Utilities()
util.convertmb2tomb3(oldtop,newtop)
pinfo = PaperInfo()
pinfo.printInfo(8)
