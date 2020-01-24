# Questions: jose.manuel.pereira@ua.pt

import argparse
from relax_decoy import relax
from pyrosetta import *
from ze_utils.pyrosetta_classes import Fragment
from ze_utils.pyrosetta_tools import uncap_selection
from pyrosetta.rosetta.core.kinematics import FoldTree
from pyrosetta.rosetta.core.pack.task import TaskFactory
from ze_utils.pyrosetta_tools import get_residues_from_selector
from pyrosetta.rosetta.protocols.loop_modeler import LoopModeler
from pyrosetta.rosetta.protocols.simple_task_operations import \
    RestrictToLoopsAndNeighbors
from pyrosetta.rosetta.protocols.loops import \
    Loop, Loops, fold_tree_from_loops
from pyrosetta.rosetta.core.select.residue_selector import \
    ChainSelector, ResidueIndexSelector

#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\

#                        Loop Close-Design Script:
# ______________________________________________________________________________
#  Performs loop closure (using KIC protocol) with automatic design of the loop
# aminoacids using the LoopModeler mover. The loop to add can be obtained from
# an existing loop (from another pose), or generated from a string. Furthermore,
# existing loops can be extended with additional residues at the C/N terminal.
#  This script assumes the ABC model, so the loop will be added between the end
# of chain A and the beginning of chain B. Read the dosctrings and documentation
# of the single_dock_decoy.py script for more information on how to set up this
# model.
#
#  How to append a new loop from a sequence:
#
#  This strategy can be employed to create a de-novo loop from scratch OR to
# append extra residues at the end of a pre-existing reference loop in order to
# accomodate the necessity for a bigger loop.
#
#  Build a dummy pose from the requested sequence. The pose should have one
# extra residue in the beginning as the pose. This anchor aminoacid nature is
# irrelevant, as only the backbone phi, psi and omega angles are necessary for
# correct positioning. Moreover, if the loop is intended to be closed to another
# part of the protein (i.e: is not a terminal loop), and in order to prevent
# pyrosetta from automatically adding the terminal caps, thus
# allowing for a correct loop closing, use the uncap_selection function in
# ze_utils.pyrosetta_tools module.
#
# Example:
# loop_pose = pose_from_sequence("A" * 7)
# ids = list(range(2, len(loop_pose.sequence()) + 1))
# loop = Fragment(loop_pose, ids)
# loop.append_to(pose, anchor)
#
# For more information and similar scripts, please read:
# https://graylab.jhu.edu/pyrosetta/downloads/scripts/test/T660_LoopBuilding.py
# https://graylab.jhu.edu/pyrosetta/downloads/scripts/demo/D080_Loop_modeling.py
# https://drive.google.com/drive/folders/1eRpTwwS7rHPoLYSLYIPVac6SvuiTVCNX


class DEFAULT:
    """
    Define defaults for the script. Values can be modified using arguments.
    """
    cut_off       = 9.0
    n_decoys      = 10
    output_prefix = "relooped"


def validate_arguments(args):
    """
    Validate the script arguments to prevent known input errors.
    """
    if not args.input_file[-4:] == ".pdb":
        exit("ERROR: Input file should be in PDB format")
    if args.cutoff < 0:
        exit("ERROR: Cut-off must be a non-negative value")
    if args.n_decoys < 0:
        exit("ERROR: Number of decoys must be a non-negative value")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Rebuild the loop between two
    parts of an input protein, between the last residue of chain A and the first
    of chain B. The new loop residues can be built from a sequence or taken from
    a reference structure. In any case, the secondary structure can be set to be
    alpha helix, beta sheet, stretched and original conformations. Currently,
    this options can only be changed in the script body itself.""")
    parser.add_argument('input_file', metavar='INPUT', type=str,
        help='The input PDB file')
    parser.add_argument('-co', '--cutoff', metavar='', type=float,
        help='Repack around loop in cut-off (in Å) (Default: %5.1f)' \
            % (DEFAULT.cut_off), default = DEFAULT.cut_off)
    parser.add_argument('-nd', '--n_decoys', metavar='', type=int,
        help='Number of decoys (Default: %d)' % (DEFAULT.n_decoys),
        default = DEFAULT.n_decoys)
    parser.add_argument('-o', '--output', metavar='', type=str,
        help='Output prefix (Default: %s)' % (DEFAULT.output_prefix),
        default = DEFAULT.output_prefix)

    args = parser.parse_args()
    validate_arguments(args)
    init()


    pose = pose_from_pdb(args.input_file)

    # Set the anchor residue at which the new loop will be appended in chain A
    anchor_A = pose.chain_end(1)

    # Rebuild loop from original loop section. Currently loop definition needs
    # to be performed manually by defining each residue ID and the original PDB
    # file, in this script. A prototype automatic identification function has
    # been developed in ze_utils.other.identify_loops, but without success.
    ref_pose = pose_from_pdb("relax_570.pdb")          # !
    ids = [99, 100, 101, 102, 103, 104, 105, 106, 107] # !
    old_loop = Fragment(ref_pose, ids)

    # When appending, the conformation can be set to "alpha", "beta" or
    # "stretched", and the fragment will adopt the new conformation. Default is
    # "auto", using the original phi and psi angles.
    # Ex: old_loop.append_to(pose, anchor_A, backbone = "stretched")
    # Ex: old_loop.append_to(pose, anchor_A, backbone = "alpha")
    # Ex: old_loop.append_to(pose, anchor_A, backbone = "beta")
    old_loop.append_to(pose, anchor_A)

    # Create the Loop object and add it to a Loops instance. A Loops instance is
    # basically a list of Loop objects, but is necessary for the CCD Mover. A
    # loop is automatically defined between the start - 1 and the end + 1.
    loop =  Loop(anchor_A - 1, anchor_A + len(ids) + 1, anchor_A + len(ids))
    loops = Loops()
    loops.add_loop(loop)

    # Set-up the pose fold tree based on the Loop description
    fold_tree = FoldTree()
    fold_tree_from_loops(pose, loops, fold_tree)
    pose.fold_tree(fold_tree)

    # Uncap the anchor residue at which the new loop will be appended in chain B
    # Uncapping means removing any atoms that wouldn't exist in a closed loop
    # (such as extra oxygens and hydrogens at C and N terminals, respectively),
    # and setting the correct pyrosetta variant. Variants are tags appended to
    # residues. By default, the LoopModeler can't perform the KIC protocol if
    # the variant is set to NTermProteinFull or CTermProteinFull (marking
    # terminals), for example.
    anchor_B = ResidueIndexSelector(pose.chain_begin(2))
    uncap_selection(pose, anchor_B)

    # LoopModeler can receive a task_factory instructing the mover to perform
    # design efforts while closing the looping. The task factory needs to be
    # populated with operations, such as the RestrictToLoopsAndNeighbors. We can
    # set this operation to design the loop and repack the neighbours in a 9
    # angstrom cutoff.
    task_factory = TaskFactory()
    loop_op = RestrictToLoopsAndNeighbors()
    loop_op.set_include_neighbors(True)
    loop_op.set_design_neighbors(False)
    loop_op.set_design_loop(True)
    loop_op.set_cutoff_distance(args.cutoff)
    loop_op.set_loops(loops)
    task_factory.push_back(loop_op)

    # The LoopModeler mover is comprised of 3 stages:
    #  - Build Stage: attempts to close the loop as is
    #  - Centroid Stage: uses centroid mode to sample the conformational space
    # and find the lowest energy position for each aminoacid on the loop
    #  - Fullatom Stage: returns the pose to full atom for high resolution
    # refinement. During this stage, the sidechains can be repacked and designed
    # Any one of this stages can be disabled:
    # Ex. loop_modeler.disable_fullatom_stage()
    # Or limited to only a few steps (for example, for test purposes)
    # Ex. loop_modeler.fullatom_stage().mark_as_test_run()
    loop_modeler = LoopModeler()
    loop_modeler.centroid_stage()
    loop_modeler.fullatom_stage()
    loop_modeler.set_loops(loops)
    loop_modeler.set_task_factory(task_factory)

    score_function = get_fa_scorefxn()
    job_man = PyJobDistributor(args.output, args.n_decoys, score_function)
    while not job_man.job_complete:
        # This 'while' loop will run 'n_decoys' times, and the second time
        # it runs it should start from the original pose, not the previously
        # designed structure.
        p = Pose(pose)
        loop_modeler.apply(p)
        relax(p)

        # Change the reconnected chain B to be continuous with chain A. This
        # prevents pyrosetta from trying to change the connection residue
        # between the two chains to a terminal variant.
        # NOTE: DOESN'T WORK! NEEDS TO BE SET MANUALLY!
        for index in range(1, len(p.residues)):
            if p.pdb_info().chain(index) == 'B':
                p.pdb_info().chain(index, 'A')

        job_man.output_decoy(p)