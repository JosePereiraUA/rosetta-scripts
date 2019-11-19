import argparse
from ze_utils.molecule_manipulation import *

#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\


#                            Blosum62 Script:
# ______________________________________________________________________________
#  Calculate blosum62 score and conserved sequence percentage between a query
# and a reference protein. The script is able to read both .PDB and .GRO files.
# The BLOSUM62 matrix used, by default, is located in
# ~/Desktop/scripts/static/b62.txt and was extracted from
# https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
#
# > Check functions color_by_sequence_blosum62_score and
# color_by_sequence_conservation on ze_utils.pymol_tools.py in order to color a
# PDB file in PyMOL based on BLOSUM62 score or conserved status
# 
# > Check multi_blosum62.py script for calculating blosum62 score and sequence
# conservation percentage on multiple query sequences towards a single
# reference sequence.

def validate_arguments(args):
    """
    Validate the script arguments to prevent known input errors.
    """
    if not args.query_file[-4:] in [".pdb", ".gro"]:
        exit("ERROR: Query file should be in PDB or GRO format")
    if not args.reference_file[-4:] in [".pdb", ".gro"]:
        exit("ERROR: Reference file should be in PDB or GRO format")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Calculate blosum62 score
        and conserved sequence percentage between a query and a reference
        protein. The script is able to read both .PDB and .GRO files.""")
    parser.add_argument('query_file', metavar='QUERY', type=str,
        help='The query PDB/GRO file')
    parser.add_argument('reference_file', metavar='REFERENCE', type=str,
        help='The reference PDB/GRO file')

    args = parser.parse_args()
    validate_arguments(args)

    query     = Molecule(args.query)
    reference = Molecule(args.reference)

    print("%15s | %15s | %8s | %%CONSERVED" % \
        ("REFERENCE", "QUERY", "BLOSUM62"))
    print(" %s" % ("-"*56))
    print("%15s | %15s | %8d | %9.2f%%" % \
        (reference.name, query.name, query.blosum62(reference), \
        query.conserved(reference)))