# Questions: jose.manuel.pereira@ua.pt

import argparse
import matplotlib.ticker as ticker
from matplotlib import pyplot as plt
from ze_utils.molecule_manipulation import *
from ze_utils.common import load_data_from_fasc_file, progress_bar


#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\


#                        Convergence Funnel Script:
# ______________________________________________________________________________
#  The objective of this script is to plot the convergence RMSD vs Energy funnel
# plot. This is achieved in several steps:
#
#  1) Read an input FASC file;
#  2) Obtain the best structure overall in terms of total score;
#  3) Calculate RMSD and energy variantion values for each structure vs the best
# structure;
#  4) If defined, calculate the RMSD value between the original and best
# structures;
#  5) Plot the results.


class DEFAULT:
    """
    Define defaults for the script. Values can be modified using arguments.
    """
    output = "convergence_funnel.png"


def validate_arguments(args):
    """
    Validate the script arguments to prevent known input errors.
    """
    if args.input_file != None and args.input_file[-5:] != ".fasc":
        exit("ERROR: Input file should be in .FASC format")
    if args.reference != None and args.reference[-4:] != ".pdb":
        exit("ERROR: Original file should be in PDB format")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Plots the energy/RMSD
        variation between each structure defined in the .FASC file and the best
        overall structure, based on the total score. If a 'reference' structure
        file is provided, show the difference between that strucutre and the
        best structure as a line in the plot.""")

    parser.add_argument('input_file', metavar='FASC', type=str,
        help='The input FASC file')
    parser.add_argument('-r', '--reference', metavar='', type=str,
        help='Reference structure for comparison purposes.')
    parser.add_argument('-o', '--output', metavar='', type=str,
        help='Filename for the output image of the plot (Default: %s).' % \
            (DEFAULT.output), default = DEFAULT.output)
    args = parser.parse_args()
    validate_arguments(args)

    # Read .FASC file and extract data
    data        = load_data_from_fasc_file(args.input_file)

    # Order extracted data based on the total score parameter
    data        = sorted(data, key = lambda entry: entry["total_score"])

    # Obtain the best energy data point
    zero_energy  = data[0]['total_score']

    # Extract filenames only from the data
    filenames   = [entry['decoy'] for entry in data]
    n_filenames = len(filenames)

    # Load each extracted filename as a Molecule object
    structures  = []
    print(" : Loading all PDB files found")
    for (index, filename) in enumerate(filenames):
        structures.append(Molecule(filename))

        # Print loading progress in a progress bar
        pb = progress_bar(index, n_filenames, 55)
        pb_per = (index / n_filenames) * 100
        sys.stdout.write("\r" + " %s %5.2f %-s" % (pb, pb_per, filename))
        sys.stdout.flush()
    print("\r" + " %s %5.2f" % (progress_bar(1, 1, 55), 100.0) + " " * 20)

    # Using the loaded Molecule objects, calculate RMSD values vs best structure
    # Additionally, extract each structure energy difference vs best structure
    rmsds, energies = [0], [0]
    n_structures = len(structures)
    print(" : Measuring RMSD and energy values")
    for index, structure in enumerate(structures[1:], start = 1):

        # Calculate RMSD
        rms = structures[0].align(structure, verbose = False)
        rmsds.append(rms)

        # Calulate energy variation
        energies.append(data[index]['total_score'] - zero_energy)

        # Print calculation progress in progress bar
        pb = progress_bar(index, n_structures, 55)
        pb_per = (index / n_structures) * 100
        sys.stdout.write("\r" + " %s %5.2f %-s" % (pb, pb_per, structure.name))
        sys.stdout.flush()
    print("\r" + " %s %5.2f" % (progress_bar(1, 1, 55), 100.0) + " " * 20)

    # Create the plot figure
    fig, ax = plt.subplots()

    # If defined, calculate the RMSD value of the original vs best structures
    # and show the result as a vertical line in the plot
    if args.reference != None:
        print(" : Calculating RMSD value for original structure")
        from align_pyrosetta import *
        from pyrosetta import *
        init()
        
        original_pose  = pose_from_pdb(args.original)
        best_pose      = pose_from_pdb(filenames[0])
        original_rmsd  = align(best_pose, original_pose)
        ax.axvline(original_rmsd, c = 'orchid')

    # Scatter plot the RMSD vs energy values
    print(" : Saving results to file %s" % (args.output))
    ax.scatter(rmsds, energies, c = 'lightseagreen')
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('+%d'))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%5.2f'))
    ax.set(xlabel='RMSD (Angstrom)', ylabel='Energy variation (REU)')
    plt.savefig(args.output, bbox_inches='tight')

    print(" : Plotting results")
    plt.show()