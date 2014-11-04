from GBClass import *
#to be added more stuff
#numpy was imported from GB_functions.py

class Utilities:

    def convertmb2tomb3(topold,topnew):
        """Convert topology file from mbondi2 to mbondi3 radii set
        """

        mol = Molecule(top=topold)
        word = "mbondi2"
        if not mol.check_word_in_file(word,topold):
            print "I can not find mbondi2 keyword in old toplogy file. Handle with care"

        #change radii for OE(GLU), OD(ASP), OXT and HH,HE(ARG)
        for i in range(mol.natoms):
            at = mol.atname_pdb[i]
            if mol.whichResidue(i) == 'GLU' and at[:2] == 'OE':
                mol.rborn[i] = '1.40000000E+00'
            if mol.whichResidue(i) == 'ASP' and at[:2] == 'OD':
                mol.rborn[i] = '1.40000000E+00'
            if mol.whichResidue(i) == 'ARG' and (at[:2] == 'HH' or at[:2] == 'HE'):
                mol.rborn[i] = '1.17000000E+00'
            if (i == mol.natoms - 1) and at[:3] == 'OXT':
                #change radii of OXT also
                mol.rborn[i] = '1.40000000E+00'
                mol.rborn[i-1] = '1.40000000E+00'

        #write to file
        fh = open(topnew,'w')
        lines = mol.topologylines
        rborn = mol.rborn
        print rborn

        for line in lines:
            if line.startswith('%FLAG RADII'):
                index = lines.index(line)
                numline = linenum(mol.natoms,mol.radii_line_format)
                k = 0
                for i in range(index+2,index+numline+2):
                    lines[i] = '  ' + '  '.join(rborn[k:k+5]) + '\n'
                    k += mol.radii_line_format

        fh.writelines(lines)
        fh.close()

    def calc_histogram(bins=100,lowerb=0,upperb=20,data=''):
        """histogram of data
        Why do I need this?
        """
        hist,bedge  = np.histogram(data,bins=bins,normed=True,range=(lowerb,uppperb))
        mean = hist
        bedge = bedge[:-1]
        return bedge, mean
