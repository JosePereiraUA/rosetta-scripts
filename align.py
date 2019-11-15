import sys
import argparse
import time

# -*- coding: utf-8 -*-
# From: https://github.com/biopython/biopython/blob/master/Bio/SVDSuperimposer/__init__.py
# From: https://github.com/biopython/biopython/blob/master/Bio/PDB/Superimposer.py

#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\


#                            Align Script:
# ______________________________________________________________________________
#  Align 2 PDB structures (optionally using a subset of atoms selected) and
# export the resulting aligned input structure to a new file.

from ze_utils.molecule_manipulation import *

def main(args):
    prediction = Molecule(args.m)
    reference = Molecule(args.r)
    prediction.align(reference, args.f)
    prediction.print_structure(args.o)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Align 2 structures.
        Export the result to a file.""")
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

    sys.exit(main(args))