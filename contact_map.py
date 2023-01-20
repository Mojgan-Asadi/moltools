import argparse
import numpy as np
import MDAnalysis as mda

parser = argparse.ArgumentParser(description="Genarate a contact matrix between amino acid side chains")
parser.add_argument("mol")
parser.add_argument("output")
parser.add_argument("--distance", type=float, default=12.0)
args = parser.parse_args()

u = mda.Universe(args.mol)

sidechains = u.select_atoms("protein and ( name CA or not backbone)").split('residue')
resnumbers = np.asarray([sc.residues[0].resid for sc in sidechains])
side_coms = np.asarray([sc.center_of_mass() for sc in sidechains])

dr_array = np.sqrt(np.sum((side_coms[:, np.newaxis, :] - side_coms[ np.newaxis, :, :])**2, axis=-1))
contacts = (dr_array < args.distance)

with open(args.output, 'w') as f:
    for i, res in enumerate(resnumbers):
        f.write(f"{res} ")
        rescontacts = np.nonzero(contacts[i, :])[0]
        for j in rescontacts:
            if i != j:
                f.write(f"{resnumbers[j]} ")
        f.write("\n")


