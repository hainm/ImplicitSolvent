from GBClass import Molecule

top = "Tc5b.ff99SB.mb3.top"
crd = "Tc5b.nat.crd"
mol = Molecule(top=top, crd=crd)
print mol.atname_pdb
