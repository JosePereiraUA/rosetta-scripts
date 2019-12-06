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
# existing loops can be extended with additional residues at the C terminal.
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

init()


pose = pose_from_pdb("candidate.pdb")

# Define the anchor residue at which the new loop will be appended in chain A
anchor_A = pose.chain_end(1)

# Rebuild loop from original loop section. Currently loop definition needs to be
# performed manually by defining each residue ID. A prototype automatic
# identification function has been developed in ze_utils.other.identify_loops,
# but without success.
ref_pose = pose_from_pdb("clamshell.pdb")
ids = [474, 475, 476, 477, 478, 479, 480, 481, 482]
old_loop = Fragment(ref_pose, ids)
old_loop.append_to(pose, anchor_A)

# Create the Loop object and add it to a Loops instance. A Loops instance is
# basically a list of Loop objects, but is necessary for the CCD Mover. A loop
# is automatically defined between the start - 1 and the end + 1 residues.
loop =  Loop(anchor_A - 1, anchor_A + len(ids) + 1, anchor_A + len(ids))
loops = Loops()
loops.add_loop(loop)

# Set-up the pose fold tree based on the Loop description
fold_tree = FoldTree()
fold_tree_from_loops(pose, loops, fold_tree)
pose.fold_tree(fold_tree)

# Uncap the anchor residue at which the new loop will be appended in chain B
# Uncapping means removing any atoms that wouldn't exist in a closed loop (such
# as extra oxygens and hydrogens at C and N terminals, respectively), and
# setting the correct pyrosetta variant. Variants are tags appended to residues.
# By default, the LoopModeler can't perform the KIC protocol if the variant is
# set to NTermProteinFull or CTermProteinFull (marking terminals), for example.
anchor_B = ResidueIndexSelector(pose.chain_begin(2))
uncap_selection(pose, anchor_B)

# LoopModeler can receive a task_factory instructing the mover to perform design
# efforts while closing the looping. The task factory needs to be populated with
# operations, such as the RestrictToLoopsAndNeighbors. We can set this operation
# to design the loop and repack the neighbours in a 9 angstrom cutoff.
task_factory = TaskFactory()
loop_op = RestrictToLoopsAndNeighbors()
loop_op.set_include_neighbors(True)
loop_op.set_design_neighbors(False)
loop_op.set_design_loop(True)
loop_op.set_cutoff_distance(9)
loop_op.set_loops(loops)
task_factory.push_back(loop_op)

# The LoopModeler mover is comprised of 3 stages:
#  - Build Stage: attempts to close the loop as is;
#  - Centroid Stage: uses centroid mode to sample the conformational space and
# find the lowest energy position for each aminoacid on the loop;
#  - Fullatom Stage: returns the pose to full atom for high resolution
# refinement. During this stage, the sidechains can be repacked and designed.
loop_modeler = LoopModeler()
loop_modeler.set_loops(loops)
loop_modeler.set_task_factory(task_factory)

loop_modeler.apply(pose)

pose.dump_pdb("design.pdb")