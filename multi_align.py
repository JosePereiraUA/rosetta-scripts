# Questions: jose.manuel.pereira@ua.pt

import os
import argparse
from ze_utils.molecule_manipulation import *

#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\


#                             Multi Align Script:
# ______________________________________________________________________________
#  Align all PDB structures in a directory (optionally using a subset of atoms
# selected) and export the resulting aligned RMSD values. This script outputs
# all combinations possible between the files.
#
# > Check align.py script to align 2 single structures among themselves.


def multi_align(elem):
    """
    Read all PDB files in the current directory and perform RMSD analysis.
    Print the RMSD pairs in descending order (Highest RMSD value first).
    """
    s = [f for f in os.listdir('.') if (os.path.isfile(f) and f[-4:] == ".pdb")]
    structures = [Molecule(f) for f in s]
    results = []
    for i in range(len(s) - 1):
        str_i = structures[i]
        for j in range(i+1, len(s)):
            str_j = structures[j]

            # Obtain masks, if provided with an element filter
            if elem is not None:
                str_i_mask = str_i.get_mask_for_elements(elem)
                str_j_mask = str_j.get_mask_for_elements(elem)
            else:
                str_i_mask, str_j_mask = None, None

            # Perform alignment, based on the created masks, if provided
            rms = str_i.align(str_j, str_i_mask, str_j_mask, verbose = False)
            results.append((str_i.name, str_j.name, rms))

    results = sorted(results, key = lambda x: x[2], reverse = True)
    for r in results:
        print("%s <-> %s : %7.3f" % (r[0], r[1], r[2]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Align all PDB files in
        current directory among eachother (ALL ATOMS) and print the
        corresponding RMSD value in descending order""")
    parser.add_argument('-ca', action='store_true',
        help='If set, align only CA vs CA instead of all atoms')
        
    args = parser.parse_args()
    elem = ["CA"] if args.ca else None
    multi_align(elem)