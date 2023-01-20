import argparse
import MDAnalysis as mda
from MDAnalysis.analysis import align

parser = argparse.ArgumentParser(description="Superimpose a structure against another.")
parser.add_argument("mol1")
parser.add_argument("mol2")
parser.add_argument("output")
args = parser.parse_args()

mobile = mda.Universe(args.mol1)
ref = mda.Universe(args.mol2)
align.alignto(mobile, ref, select='protein and backbone', weights="mass")

mobile.atoms.write(args.output)

