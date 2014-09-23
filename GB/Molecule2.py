from GBClass import Molecule
import re

#parse_file function was copied and modified from
#https://github.com/pele-python/pele
data_casts = {"a": str,
              "A": str,
              "e": float,
              "E": float,
              "i": int,
              "I": int}


class Molecule2(Molecule):
    def __init__(self,top='',crd=''):
        Molecule.__init__(self,top='',crd='')
        self.topology_data = {}


    def parse_file(self, fname):
        '''this parse_file was copied (and modified) from above link'''
        flag = None
        format = None
        data_length = None
        data_type = None

        lines = open(fname,'r').readlines()
        for line in lines:
            # Remove whitespace and newlines
            line = line.rstrip()
            if line.startswith("%VERSION"): continue
            elif line.startswith("%FLAG"):
                flag = line[6:]
                #AMBER topology file uses upper case for flag
                #I prefer lower case
                self.topology_data[flag.lower()] = []
            elif line.startswith("%FORMAT"):
                format      = line[8:-1]
                data_length = int(re.split("[aEI\.\)]", format)[1])
                data_type   = re.findall("[a-zA-Z]", format)[0]
            else:
                split_line = [line[i:i+data_length] for i in range(0, len(line), data_length)]
                formatted_line = map(data_casts[data_type], split_line)
                self.topology_data[flag.lower()] += formatted_line
        return self.topology_data

if __name__ == '__main__':
    top = "Tc5b.ff99SB.mb3.top"
    crd = "Tc5b.nat.crd"
    mol = Molecule(top=top, crd=crd)
    print mol.atoms[0].xyz

