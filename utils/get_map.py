import pyproteinsExt.structure.coordinates as PDB
import ccmap
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-f', help='File containing the original structure of the protein in pdb format.')
parser.add_argument('-t', help='Threshold', default=7.0, type=float)
parser.add_argument('-o', help='JSON out file (default: stdout)', default=sys.stdout)
args = parser.parse_args()

parser = PDB.Parser()
pdbObj = parser.load(file=args.f)

atoms = pdbObj.atomDictorize

res = ccmap.compute(
  atoms['x'],
  atoms['y'],
  atoms['z'],
  atoms['chainID'],
  atoms['seqRes'],
  atoms['resName'],
  atoms['name'],
  args.t
)

out = args.o

if type(out) is str:
  out = open(file=out, mode='w')
# else: the argument is a IOWrapper

out.write(res + '\n')
