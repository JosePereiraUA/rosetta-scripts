# Questions: jose.manuel.pereira@ua.pt

import json
import argparse
from sys import argv
from ze_utils.common import load_data_from_fasc_file

#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\


def extract_n_decoys_by_parameter(data, parameter, n, reverse = False):
    """
    Extract the name of the 'n' decoys with the lowest values of 'paramater'.
    If 'reversed = True', extract the highest value first.
    """
    data = sorted(data, key=lambda entry: entry["total_score"], reverse=reverse)
    extracted = []
    for entry in data[0:n]:
        print(str(entry["decoy"]), entry[parameter])
        extracted.append((str(entry["decoy"]), entry[parameter]))
    return extracted


def validate_arguments(args):
    """
    Validate the script arguments to prevent known input errors.
    """
    if args.n_files < 0:
        exit("ERROR: Number of files to print must be a non-negative value")


class DEFAULT:
    """
    Define defaults for the script. Values can be modified using arguments.
    """
    n_files = 1
    parameter = 'total_score'


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Print the names of the
        'n' files with the LOWEST given parameter 'p', from the input .fasc
        file""")
    parser.add_argument('input_file', metavar='FASC', type=str,
        help='The input FASC file')
    parser.add_argument('-n', '--n_files', metavar='', type=int,
        help='Print the first N_FILES (Default: %d)' % (DEFAULT.n_files),
        default = DEFAULT.n_files)
    parser.add_argument('-p', '--parameter', metavar='', type=str,
        help='Use this parameter when sorting the files (Default: %s)'
        % (DEFAULT.parameter), default = DEFAULT.parameter)
    parser.add_argument('-r', '--reverse', action='store_true',
        help='If set, print the HIGHEST values instead')
    args = parser.parse_args()
    validate_arguments(args)

    data = load_data_from_fasc_file(args.input_file)
    extract_n_decoys_by_parameter(data,
        parameter = args.parameter,
        n = args.n_files,
        reverse = args.reverse)