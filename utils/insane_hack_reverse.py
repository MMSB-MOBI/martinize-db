import argparse
import re

parser = argparse.ArgumentParser(description="Replace hacked coordinates by nan coordinates in gro files. To reverse insane_hack.py modifications. It will put nan as coordinates for atoms that are specified in atoms txt file")
parser.add_argument('input', type=str, help="input.gro")
parser.add_argument('output', type=str, help="output.gro")
parser.add_argument('atoms', type=str, help="txt file with hacked atoms line")
args=parser.parse_args()

atom_line_re = re.compile('^[\d\s]{5}.{5}.{5}[\d\s]{5}[\d\s-]{4}\.[\d]{3}[\d\s-]{4}\.[\d]{3}[\d\s-]{4}\.[\d]{3}')
previous_coords = ""
out = open(args.output, "w")

with open(args.atoms, "r") as atoms:
    hacked_atoms = [int(l.rstrip()) for l in atoms]

atom_nb = 0
with open(args.input) as inp: 
    for l in inp:
        if atom_line_re.match(l):
            if atom_nb in hacked_atoms:
                print(atom_nb)
                l = l[:20] + "     nan     nan     nan" + l[44:]
            atom_nb += 1
        out.write(l)
        

out.close()