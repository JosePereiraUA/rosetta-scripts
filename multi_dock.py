# Questions: jose.manuel.pereira@ua.pt

import os
import glob
import shutil
import argparse
import subprocess
import numpy as np
from pyrosetta import *
from ze_utils.common import overwrite_dir
from ze_utils.pyrosetta_classes import DockingGrid
from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from single_dock import \
    deploy_decoys_on_slurm, deploy_decoys_on_pyjobdistributor
from ze_utils.pyrosetta_tools import \
    get_centroid_coordinates_from_selector, load_pre_filter

#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\

#               Multi Starting Position Dock-Design Script:
# ______________________________________________________________________________
#  Performs docking + design algorithm from multiple starting positions.
#
#  The docking + design protocol is performed by calling the single_dock.py 
# script. Read its docstrings and comments for more detailed information.
# Furthermore, and in addition to the single_dock.py script, this script will
# attempt to initiallize multiple docking + design protocols from different
# starting positions defined in a docking grid. A default docking grid is
# provided in the script, and can be changed by editing this file. If 'slurm'
# option the arguments is set to True, instead of calculating each of the new
# starting positions simulations in sequential order, individual slurm bash
# scripts will be automatically created and run, spawning 'n_decoys' processes
# that will each be in charge of one of the decoys for each simulation.
#
#  The multi_dock.py script will create or overwrite N dock folders, where N is
# the number of points on the DockingGrid, containing:
#
#  1) For each dock 'dock_n/', where n is the current docking point:
#   start_n.pdb         : This dock starting position, on the grid
#
#  2) For every decoy on every dock, where d is the current decoy:
#   n_d.sh              : The slurm launching bash script*
#   n_d_pre_filter.json : The JSON file with the pre filter parameters*
#   n_d.err             : The SLURM standard error receiver*
#   n_d.out             : The SLURM standard output receiver*
#   n_d_status.txt      : Holds current simulation step
#   n_d_energy.dat      : Holds accepted steps energy data
#   n_d_data.json       : Holds all steps pre filter acceptance/rejection data
#
#  3) For every accepted step on every decoy, where s is the current step:
#   n_d_s.pdb           : The new proposed and accepted structure
#   n_d_s.txt           : The PyMOL selection for the designed region


class DEFAULT:
    """
    Define defaults for the script. Values can be modified using arguments.
    """
    n_decoys   = 100
    n_steps    = 2000


def validate_arguments(args):
    """
    Validate the script arguments to prevent known input errors.
    """
    if not args.input_file[-4:] == ".pdb":
        exit("ERROR: Input file should be in PDB format")
    if args.n_decoys < 0:
        exit("ERROR: Number of decoys must be a non-negative value")
    if args.n_steps < 0:
        exit("ERROR: Number of PASSO steps must be a non-negative value")


# ______________________________________________________________________________
#                                 M A I N
# ______________________________________________________________________________

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = """Perform a docking+design
        algorithm from multiple starting positions defined in a custom grid.""")
    parser.add_argument('input_file', metavar = 'INPUT', type = str,
        help = 'The input PDB file')
    parser.add_argument('-nd', '--n_decoys', metavar = '', type = int,
        help = 'Number of decoys (Default: %d)' % (DEFAULT.n_decoys),
        default = DEFAULT.n_decoys)
    parser.add_argument('-ns', '--n_steps', metavar = '', type = int,
        help = 'Number of PASSO protocol steps (Default: %d)' % \
        (DEFAULT.n_steps), default = DEFAULT.n_steps)
    parser.add_argument('-s', '--slurm', action = 'store_true',
        help = "Use SLURM to launch parallel decoys (Default: False)")
    parser.add_argument('-d', '--dry_run', action = 'store_true',
        help = """Dry Run: Export only the starting positions to the correct
        folder, but don't run simulations on them (Default: False)""")
    parser.add_argument('-pf', '--pre_filter', metavar='', type=str,
        help='Input pre filter JSON file (Default: auto)',
        default = "auto")

    args = parser.parse_args()
    validate_arguments(args)
    init()

    # Remove old slurm file that may exist in the current folder
    if args.slurm:
        slurms = glob.glob("slurm.*")
        for slurm in slurms:
            os.remove(slurm)

    # Load the ABC model pose
    pose           = pose_from_pdb(args.input_file)
    score_function = get_fa_scorefxn()

    # 'load_pre_filter' returns a default PreFilter if no JSON file is provided
    # Any changes to single default values can be made after the loading of the
    # pre filter.
    pre_filter = load_pre_filter(args.pre_filter)
    # Ex. pre_filter.contact_min_count = 6

    # Define the center for the grid. By default, and when considering the ABC
    # model, this is the position in the vector AB at the same distance from
    # the centroid of the ligand as the initial chain C position.
    centroidA = get_centroid_coordinates_from_selector(ChainSelector("A"), pose)
    centroidB = get_centroid_coordinates_from_selector(ChainSelector("B"), pose)
    centroidC = get_centroid_coordinates_from_selector(ChainSelector("C"), pose)
    cb_distance = np.linalg.norm(centroidB - centroidC)
    ac_normal = (centroidC - centroidA) / np.linalg.norm(centroidC - centroidA)
    grid_center = centroidC + cb_distance * ac_normal

    # Create the grid around the calculated grid center. The hardcoded values on
    # this script constitute a default docking grid, but can be freely modified.
    # Check the DockingGrid class docstring for more detailed information on how
    # to do this.
    dg = DockingGrid(grid_center,
        (10, 40), (2, 8), (20, 20), 4, 3, 4, [centroidB])

    # Print to a configuration file the initial conditions of the simulation:
    #  Line 1 : Number of docks, decoys and steps
    #  Line 2 : Initial total score of the pose
    #  Line 3 : Initial conformation file name
    with open("init.conf", "w") as conf:
        conf.write("%d %d %d\n" % (len(dg.points), args.n_decoys, args.n_steps))
        conf.write("%f\n" % (score_function(pose)))
        conf.write("%s\n" % (args.input_file))
    
    # Orient the grid. In this case, the grid's X [1, 0, 0] axis should match
    # the orientation of the length of the ligand (this is the DE vector, where
    # D and E are the alpha carbon atoms of the first and last residues on the
    # ligand), therefore locking 2 of the 3 rotation degrees of freedom of the
    # grid. The third rotation degree of freedom will be left at the original
    # state (matching the natural axis). 
    lig_res = get_residues_from_subset(ChainSelector("C").apply(pose))
    D = np.array(pose.residue(list(lig_res)[0]).atom(2).xyz())
    E = np.array((pose.residue(list(lig_res)[len(lig_res) - 1]).atom(2).xyz()))
    de_normal = (E - D) / np.linalg.norm(D - E)
    dg.orient(np.array([1, 0, 0]), de_normal)

    # Ex. Print the current docking grid starting points in PDB format.
    dg.as_pdb("grid.pdb")

    # Looping over the number of points in the docking grid, move the Chain B   
    # of the target protein to the new positions and run the simultaneous
    # position and design optimization protocol on each one, asynchronously.
    p = Pose()
    for point_id in range(len(dg.points)):
        # Create a copy of the original position of the Chain B of the target
        # protien.
        p.assign(pose)
        # Translate the chain B from the original position to the new position
        # defined the the grid point.
        dg.translate_selector_to_point(point_id, ChainSelector("B"), p)

        # Sabe this new starting position in it's own folder. Will overwrite any
        # old data in such folder.
        overwrite_dir("dock_%d" % (point_id))
        p.dump_pdb("dock_%d/start_%d.pdb" % (point_id, point_id))

        # A dry run ONLY creates the starting positions in each new directory
        if args.dry_run:
            continue

        current_work_directory = os.getcwd()
        os.chdir("dock_%d" % (point_id))

        if args.slurm:
            deploy_decoys_on_slurm(
                "start_%d.pdb" % (point_id),
                "%d" % (point_id),
                args.n_decoys,
                args.n_steps,
                pre_filter)
        else:
            from hb_design import get_designer_mover
            design_mover = get_designer_mover(192, 10, 1)
            deploy_decoys_on_pyjobdistributor(
                "start_%d.pdb" % (point_id),
                "%d" % (point_id),
                args.n_decoys,
                args.n_steps,
                score_function,
                pre_filter,
                design_mover)

        os.chdir(current_work_directory)