import sys


class CheckInput:

    def check_igb(self, igb):
        if igb not in [1, 2, 5, 7, 8]:
            print "igb must be 1, 2, 5, 7 or 8"
            sys.exit()

    def check_top_crd(self, mol):
        '''check whether topology and crd files have the same NATOMS'''
        pass
