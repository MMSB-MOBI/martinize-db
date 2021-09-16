import argparse

def open_and_write(input_file, output_file, atom_file, first_coords=None):
    print("open and write", first_coords)
    sweep=False
    current_coords = None
    atom_nb = 0
    out = open(output_file, "w")
    inp = open(input_file)
    atoms = open(atom_file,"w")
    for l in inp:
        if l.startswith("ATOM"):
            coords = l[30:54]
            if "nan" in coords:
                atoms.write(str(atom_nb) + "\n")
                if not current_coords:
                    if not first_coords:
                        sweep=True
                    else:
                        l = l[:30] + first_coords + l[54:]
                else:    
                    l = l[:30] + current_coords + l[54:]
            else:
                current_coords = coords
                if sweep:
                    break
            atom_nb += 1

        out.write(l)

    inp.close()
    out.close()
    atoms.close()
    
    if sweep:
        open_and_write(input_file, output_file, atom_file, current_coords)

    

parser = argparse.ArgumentParser(description="Replace nan coordinates by hack coordinates in pdb files. Required to make insane works. It will replace nan by previous coordinates (or by the next one if nan is in first line)\n")
parser.add_argument('input', type=str, help="input.pdb")
parser.add_argument('output', type=str, help="output.pdb")
parser.add_argument('atoms', type=str, help="List of atoms that have been hacked (txt file)")

args=parser.parse_args()

open_and_write(args.input, args.output, args.atoms)