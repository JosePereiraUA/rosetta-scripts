import json
import argparse
from sys import argv

#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\

def load_data_from_fasc_file(file_path):
    """
    Returns an array with all Dictionary entries found in the .fasc file.
    """
    json_path = "["
    with open(file_path, "r") as file_in:
        for line in file_in:
            json_path += line[:-1]
            json_path += ","
    json_path = json_path[:-1] + "]"
    return json.loads(json_path)


def extract_n_decoys_by_parameter(data, parameter, n, reverse):
    """
    Extract the name of the 'n' decoys with the lowest values of 'paramater'.
    If 'reversed = True', extract the highest value first.
    """
    data = sorted(data, key=lambda entry: entry["total_score"], reverse=reverse)
    for entry in data[0:n]:
        print(str(entry["decoy"]), entry[parameter])


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