# Questions: jose.manuel.pereira@ua.pt

import argparse
from ze_utils.molecule_manipulation import *

#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\


#                            Align Script:
# ______________________________________________________________________________
#  Align 2 PDB structures (optionally using a subset of atoms selected) and
# export the resulting aligned input structure to a new file. This script
# supports .PDB, .GRO and .XYZ files, both for input and output.
#
# Based on the following script:
# https://github.com/biopython/biopython/blob/master/Bio/SVDSuperimposer/
#
# > Check multi_align.py to align more than 2 structures in 1 on 1 combinations
# between multiple files in a single folder


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Align 2 structures.
        Export the result to a file. This script supports .PDB, .GRO and .XYZ
        files, both for input and output.""")
    parser.add_argument('-m', required=True, type=str,
        help='Structure that will suffer the rotation and be outputed',
        metavar='MOVABLE')
    parser.add_argument('-r', required=True, type=str,
        help='Reference structure', metavar='REFERENCE')
    parser.add_argument('-o', required=True, type=str,
        help='Output file. Format will be infered from file extension.',
        metavar='OUTPUT')
    parser.add_argument('-f', type=str, nargs='+',
        help='Filter alignment by atom element. (Optional)', metavar='ELEM')
    args = parser.parse_args()

    # 1) Load molecules
    movable   = Molecule(args.m)
    reference = Molecule(args.r)
    
    # 2) Obtain masks, if provided with an element filter
    if args.f is not None:
        movable_mask   = movable.get_mask_for_elements(args.f)
        reference_mask = reference.get_mask_for_elements(args.f)
    else:
        movable_mask, reference_mask = None, None

    # 3) Perform alignment, based on the created masks, if provided
    movable.align(reference, movable_mask, reference_mask)

    # 4) Print the alignment to a file
    movable.export(args.o)