# Questions: jose.manuel.pereira@ua.pt

import json
import argparse
from pyrosetta import *
from ze_utils.pyrosetta_tools import set_ABC_model_fold_tree, load_pre_filter
from ze_utils.pyrosetta_classes import PASSO

#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\

#                      Single Decoy Dock-Design Script:
# ______________________________________________________________________________
#  Performs docking + design algorithm (PASSO protocol). The objetive of this
# script is to sample the position of a movable chain on the input structure,
# and for each provable new position, sample the sequence in the interface in
# order to design, at the same time, both the position and the sequence of part
# of a protein. The script assumes an ABC model for the protein. This model
# essentially assumes that the protein is split in three different chains on the
# input PDB file. Those parts should be defined as follows:
#  - Chain A: (Fixed) Fixed part of the protein (PBZ, in this example);
#  - Chain B: (Movable) Movable part of the protein, which will be subject to
# the pseudo docking process (random movement + design);
#  - Chain C: (Fixed) Target peptide ligand around which the interface will be 
# defined;
#
#  Furthermore, in order to recreate the existing loop (or a new one), the order
# of the chains is important, and should be respected: A -> B -> C. In order to
# maintain the fold_tree coherence and minimize unexpected errors, all chains
# must be continuous (not have intra-chain breaks such as AAAA-B-AAAA).
#  The linker loop between any of this regions can (and should) be removed. This
# model can be set-up using the 'ze_utils.molecule_manipulation' module, by
# manually setting the chains between a range of consecutive residues and
# deleting the loop residues, as follows:
#
# pdb = Molecule("input.pdb")             # Load structure
# pdb.set_chain_from_range('A', 376, 473) # Set Chain A
# pdb.set_chain_from_range('B', 483, 572) # Set Chain B
# pdb.set_chain_from_range('C',   1,   7) # Set Chain C
# pdb.remove_residues_in_range(474, 482)  # Remove atoms in the loop region
# pdb.sort_residues_by_chain()            # (Necessary) Renumber all residues
# pdb.export("output.pdb")                # Save edited structure
#
#  Or by automatically identifying chains from the connection graphs, as long as
# the chains are deterministically separated in the CONECT records of the input
# PDB file:
#
# pdb = Molecule("input.pdb")            # Load structure
# pdb.remove_residues_in_range(474, 482) # Remove atoms in the loop region
# pdb.define_chains_from_connections()   # Automatically identify chains
# pdb.renumber_residues()                # (Necessary) Renumber all residues
# pdb.export("output.pdb")               # Save edited structure
#
# Note: By renumbering the residues, the missing loop anchors will have
# a different index.
# 
# Note: Pyrosetta automatically caps these new terminals where the loop was. In
# order to uncap the terminals, check 'uncap' and 'uncap_selection' functions
# in ze_utils.pyrosetta_tools module.
#
#    THIS SCRIPT PERFORMS ONLY ONE DECOY ON THE PASSO PROTOCOL.
#  > Check single_dock.py script for multiple decoy PASSO simulation
#  > Check multi_dock.py script for multiple starting positions PASSO simulation
#  > Check loop.py for loop rebuilding/design
#
# For more information and similar scripts, please read:
# https://graylab.jhu.edu/pyrosetta/downloads/scripts/demo/D100_Docking.py


class DEFAULT:
    """
    Define defaults for the script. Values can be modified using arguments.
    """
    n_steps    = 2000
    output     = "output"
    pre_filter = "auto"

def validate_arguments(args):
    """
    Validate the script arguments to prevent known input errors.
    """
    if not args.input_file[-4:] == ".pdb":
        exit("ERROR: Input file should be in PDB format")
    if args.n_steps < 0:
        exit("ERROR: Number of PASSO steps must be a non-negative value")


def single_dock_decoy(input_file, output_prefix, n_steps, docker = "auto"):
    """
    Launch a new PASSO protocol from the 'input_file' pose (must be a PDB file).
    All output files from the PASSO protocol will have the 'output_prefix'. The
    protocol will run for 'n_steps', using the given 'docker' (when set to
    "auto" will use the default PASSO protocol).
    """

    assert input_file[-4:] == ".pdb", \
        "Input file for PASSO protocol must be in PDB format (in in %s)." % \
            (input_file[-4:])

    if docker == "auto":
        docker = PASSO(n_steps)
    else:
        assert type(docker) == PASSO, "Docker object must be of type PASSO."

    # Load the pose and the score function
    pose = pose_from_pdb(input_file)
    score_function = get_fa_scorefxn()

    set_ABC_model_fold_tree(pose)

    # PASSO is pre-defined to print the best results during the simulation,
    # since it is based on Monte-Carlo and the last position does not
    # necessarily correspond to the best overall result. Therefore, the name
    # of the PDB needs to point towards the corrent directory.
    pose.pdb_info().name("%s.pdb" % (output_prefix))

    # Apply the PASSO protocol to the new starting position pose
    docker.apply(pose)

    return pose


# ______________________________________________________________________________
#                                 M A I N
# ______________________________________________________________________________
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Perform a docking+design
        algorithm.""")
    parser.add_argument('input_file', metavar='INPUT', type=str,
        help='The input PDB file')
    parser.add_argument('-ns', '--n_steps', metavar='', type=int,
        help='Number of PASSO protocol steps (Default: %d)' % (DEFAULT.n_steps),
        default = DEFAULT.n_steps)
    parser.add_argument('-o', '--output', metavar='', type=str,
        help='Output name prefix (Default: %s)' % (DEFAULT.output),
        default = DEFAULT.output)
    parser.add_argument('-pf', '--pre_filter', metavar='', type=str,
        help='Input pre filter JSON file (Default: %s)' % (DEFAULT.pre_filter),
        default = DEFAULT.pre_filter)

    args = parser.parse_args()
    validate_arguments(args)
    init()

    # load_pre_filter returns a default PreFilter if no JSON file is provided
    # Any changes to single default values can be made after the loading of the
    # pre filter.
    pre_filter = load_pre_filter(args.pre_filter)
    # Ex. pre_filter.contact_min_count = 6
    
    from hb_design import get_designer_mover
    test = get_designer_mover(192, 10, 1, "../rotlib.dat")
    docker = PASSO(args.n_steps, pre_filter = pre_filter, design_mover = test)
    single_dock_decoy(args.input_file, args.output, args.n_steps, docker)