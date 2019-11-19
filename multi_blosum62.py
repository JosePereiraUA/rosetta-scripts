import os
import argparse
from ze_utils.molecule_manipulation import *

#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\


#                          Multi Blosum62 Script:
# ______________________________________________________________________________
#  Calculate both the Blosum62 and conserved sequence percentage between all
# PDB or GRO structures in a directory and a reference structure.
#
# > Check blosum62.py script to make the same calculations between 2 single
# structures.


class DEFAULT:
    """
    Define defaults for the script. Values can be modified using arguments.
    """
    parameter = 'blosum62'


def validate_arguments(args):
    """
    Validate the script arguments to prevent known input errors.
    """
    if not os.path.isdir(args.query_folder):
        exit("ERROR: Query should be a folder")
    if not args.reference[-4:] in [".pdb", ".gro"]:
        exit("ERROR: Reference file should be in PDB or GRO format")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Calculate blosum62 score
        and conserved sequence percentage between a query and a reference
        protein. The script is able to read both .PDB and .GRO files.""")
    parser.add_argument('query_folder', metavar='QUERY', type=str,
        help='Folder where all PDB/GRO file will be compared to the reference')
    parser.add_argument('reference', metavar='REFERENCE', type=str,
        help='The reference PDB/GRO file')
    parser.add_argument('-s', '--sort_by', metavar='', type=str,
        help='Use blosum62 or conserved when sorting the files (Default: %s)'
        % (DEFAULT.parameter), default = DEFAULT.parameter)
    parser.add_argument('-r', '--reverse', action='store_false',
        help='If set, print the LOWEST values instead')

    args = parser.parse_args()
    validate_arguments(args)

    # Load the reference molecule
    reference = Molecule(args.reference)

    # Print header
    print("%15s | %15s | %8s | %%CONSERVED" % \
        ("REFERENCE", "QUERY", "BLOSUM62"))
    print(" %s" % ("-"*56))

    # Load query molecules
    s = []
    for f in os.listdir(args.query_folder):
        if (os.path.isfile(f) and f[-4:] in [".pdb", ".gro"]):
            s.append(f)
    if len(s) == 0:
        exit("No PDB or GRO files found in the query directory.")
    query_molecules = [Molecule(f) for f in s]

    # Calculate blosum62 and conserved values
    data = []
    for query in query_molecules:
        data.append({
            "name": query.name,
            "blosum62": query.blosum62(reference),
            "conserved": query.conserved(reference)})
    
    # Sort data based on args.sort_by
    data = sorted(data, key=lambda entry: entry[args.sort_by], \
        reverse=args.reverse)

    # Print the results to the user
    for entry in data:
        print("%15s | %15s | %8d | %9.2f%%" % \
            (reference.name, entry["name"], entry["blosum62"], \
            entry["conserved"]))