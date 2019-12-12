import json
import argparse
import matplotlib.ticker as ticker
from matplotlib import pyplot as plt
from ze_utils.molecule_manipulation import *


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


class DEFAULT:
    """
    Define defaults for the script. Values can be modified using arguments.
    """
    original = -1


def validate_arguments(args):
    """
    Validate the script arguments to prevent known input errors.
    """
    if args.original != None and args.original[-4:] != ".pdb":
        exit("ERROR: Original file should be in PDB format")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Plots the energy/RMSD
        variation between each structure defined in the .FASC file and the best
        overall structure, based on the total score. If an highlign index is 
        provided, that index on the ordered list will be highligthed in red.""")
    parser.add_argument('input_file', metavar='FASC', type=str,
        help='The input FASC file')
    parser.add_argument('-o', '--original', metavar='', type=str,
        help='Original structure for comparison purposes.')
    args = parser.parse_args()
    validate_arguments(args)

    data            = load_data_from_fasc_file(args.input_file)
    data            = sorted(data, key = lambda entry: entry["total_score"])
    structure_names = [entry['decoy'] for entry in data]
    structures      = [Molecule(f) for f in structure_names]

    zero_energy     = data[0]['total_score']

    rmsd, energies = [0], [0]
    for index, structure in enumerate(structures[1:], start = 1):
        rms = structures[0].align(structure, verbose = False)
        rmsd.append(rms)
        energies.append(data[index]['total_score'] - zero_energy)

    fig, ax = plt.subplots()
    ax.scatter(rmsd, energies, c = 'lightseagreen')
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('+%d'))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%5.2f'))
    ax.set(xlabel='RMSD (Angstrom)', ylabel='Energy variation (REU)')

    if args.original != None:
        from align_pyrosetta import *
        from pyrosetta import *
        init()
        original_pose  = pose_from_pdb(args.original)
        best_pose      = pose_from_pdb(structure_names[0])
        rmsd           = align(best_pose, original_pose)
        ax.axvline(rmsd, c = 'orchid')
    
    plt.show()