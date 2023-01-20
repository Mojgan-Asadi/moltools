import argparse
import pathlib as path
import numpy as np
import MDAnalysis as mda


# Auxiliary function to merge
def merge(*grps):
    combined = mda.Merge(*grps)
    combined_chains = np.concatenate([a.chainIDs for a in grps])
    combined_types = np.concatenate([a.types for a in grps])
    combined.add_TopologyAttr("chainID", combined_chains)
    combined.add_TopologyAttr("type", combined_types)
    return combined


# Function to add a Cl ion to a system
def add_ion(pos, *a, ion='Cl'):
    clu = mda.Universe.empty(1,
                         n_residues=1,
                         atom_resindex=[0],
                         residue_segindex=[0],
                         trajectory=True)
    clu.add_TopologyAttr('name', [ion.upper()])
    clu.add_TopologyAttr('type', [ion])
    clu.add_TopologyAttr('resname', [ion.upper()])
    clu.add_TopologyAttr("resid", [1])
    clu.add_TopologyAttr('elements', [ion])
    clu.add_TopologyAttr("chainID", ['Y'])
    coords = np.asarray(pos)
    clu.atoms.positions = coords

    combined = merge(*a, clu.atoms)

    return combined


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Insert an ion in the center of a residue selection")
    parser.add_argument('pdb_file')
    parser.add_argument('ressel_file')
    parser.add_argument('ion')

    args = parser.parse_args()
    pdb_file = path.Path(args.pdb_file)
    pdb_name = str(pdb_file.stem)
    ion = args.ion
    # import the PDB
    u = mda.Universe(pdb_file)
    ressel_file = args.ressel_file
    with open(ressel_file) as f:
        res_groups = f.readlines()

    for i, sel_str in enumerate(res_groups):
        sel = f"({sel_str}) and (name CA or not backbone)"

        atomsel : mda.AtomGroup = u.select_atoms(sel)
        c = atomsel.center_of_geometry()
        u_with_ion = add_ion(c, u.atoms, ion=ion)
        out_file = pdb_file.with_stem(pdb_name+f'_{ion}_{i}')
        print(f"{i}: {c}, {out_file}")
        u_with_ion.atoms.write(out_file)

