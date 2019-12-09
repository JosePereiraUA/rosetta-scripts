from pyrosetta import *
from pyrosetta.rosetta.protocols.moves import *
from pyrosetta.rosetta.std import map_int_int
from pyrosetta.rosetta.core.pack.task.operation import *
from pyrosetta.rosetta.core.select.residue_selector import *
from pyrosetta.rosetta.protocols.minimization_packing import MinMover, \
    PackRotamersMover

from ze_utils.pyrosetta_tools import *
import argparse
import sys
import os

#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\


#                            Design Script:
# ______________________________________________________________________________
#  Perform a sequence design by mutating the nature of the aminoacids in the
# structure. PyRosetta does this by changing the sidechain and its rotamers.
# The design effort will be localized in the surrounding aminoacids of a ligand
# peptide, defined as the chain B of the input PDB structure.
#
#  Necessary parts:
# - Starting Pose;
# - Score function;
# - Rotamer Packer (samples both the conformation and nature of sidechain
#   rotamers for each residue in a simulated-annealing Monte Carlo algorithm);
#   Is usually comprised of two components:
#    . PackerTask: Defines what residues are allowed to move and what new 
#      rotamers are allowed to be tried on;
#    . PackRotamersMover: The actual Monte Carlo algorithm object;
#   However, the PackerTask information is static and only relevant for the Pose
#   it was constructed for. A TaskFactory produces new PackerTasks every step of
#   PackRotamersMover, adapting as things change during the simulation. The list
#   of instructions the define new PackerTasks are called TaskOperations. For
#   more infromation, read Unit 6 of PyRosetta User's Manual;
# - Minimizer (run after the Monte Carlo to eliminate small clashes). The
#   default type is set to "linmin".
#   It requires a MoveMap (which defined the conformational degrees of freedom
#   allowed). A MoveMap stores information per residue and contains two switches
#   "is backbone movement allowed" and is "chi movement allowed", both set to
#   False by default. For more information, read Unit 5.5 of PyRosetta User's
#   Manual and the Rosetta Commons documentation -
# (/rosetta_basics/structural_concepts/minimization-overview);
#
#  Optinal parts:
# - Selectors (perform design on only a subset of residues);
#   TaskOperations can optinally be done at the level a single residue. This 
#   type of TaskOperations are called ResLevelTaskOperations. The application or
#   not of a TaskOperation of this type can be controlled using an
#   OperateOnResidueSubset command that requires the target TaskOperation and a
#   ResidueSubset (this is an array of booleans of 'n_atoms' size, where each
#   entry states whether the target TaskOperation is performed or not on the
#   i-th residue of the Pose. It can be obtained from ResidueSelectors).
#   For more information, read the Rosetta Commons documentation on the subject:
# (/scripting_documentation/RosettaScripts/ResidueSelectors/ResidueSelectors)
# - SequenceMover (groups individual movers and performs one after the other);
#
#                          PyRosetta User's Manual: 
# https://graylab.jhu.edu/pyrosetta/downloads/documentation/PyRosetta_Manual.pdf
#
#                  Rosetta Common Latest Documentation:
# https://www.rosettacommons.org/docs/latest (*add the topic here*)
# ______________________________________________________________________________

class DEFAULT:
    """
    Define defaults for the script. Values can be modified using arguments.
    """
    chain         = "B"
    cut_off       = 9.0
    n_decoys      = 100
    n_cycles      = 4
    output_prefix = "design"
    export_input  = None


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
    if args.n_cycles < 0:
        exit("ERROR: Number of cycles must be a non-negative value")


def export_input_file(args,
    pose,
    selections = None,
    score_function = None,
    sequence_mover = None,
    f_out = sys.stdout,
    clean_print = True):

    if clean_print:
        export_input_file(args, pose, selections, score_function,
        sequence_mover, f_out = open(os.devnull, "w"), clean_print = False)
    print("Pose:", file = f_out)
    print(pose.pdb_info(), file = f_out)
    print("\nRuntime arguments:", file = f_out)
    print(" Chain: %s" % (args.chain), file = f_out)
    print(" Cut off: %9.3f" % (args.cutoff), file = f_out)
    print(" Number of decoys: %d" % (args.n_decoys), file = f_out)
    print(" Number of cycles: %d" % (args.n_cycles), file = f_out)
    print(" Output name prefix: %s" % (args.output), file = f_out)
    if not selections == None:
        print("\nSelections:")
        for i, sel in enumerate(selections):
            res = str(list(get_residues_from_subset(sel.apply(pose))))
            pymol = get_pymol_selection_from_selector(sel, pose)
            print("--- Selection %d ---" % (i), file = f_out)
            print(" Type: %s" % (sel.class_name()), file = f_out)
            print(" Residues: %s" % (res), file = f_out)
            print(" PyMOL: %s" % (pymol), file = f_out)
    if not score_function == None:
        print("\nScore Function: %s" % (score_function.get_name()), file = f_out)
        print(score_function.show(pose), file = f_out)
    if not sequence_mover == None:
        print("\nSequence Mover:", file = f_out)
        for i in range(sequence_mover.nr_moves()):
            print(" %s" % (sequence_mover.get_mover(i)), file = f_out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Design the active-site of a
        protein surrousing a ligand defined in another chain, on the input
        PDB.""")
    parser.add_argument('input_file', metavar='INPUT', type=str,
        help='The input PDB file')
    parser.add_argument('-c', '--chain', metavar='', type=str,
        help='The ligand chain identifier (Default: %s)' % (DEFAULT.chain),
        default = DEFAULT.chain)
    parser.add_argument('-co', '--cutoff', metavar='', type=float,
        help='Active site identification cut-off in Å (Default: %5.1f)' \
            % (DEFAULT.cut_off), default = DEFAULT.cut_off)
    parser.add_argument('-nd', '--n_decoys', metavar='', type=int,
        help='Number of decoys (Default: %d)' % (DEFAULT.n_decoys),
        default = DEFAULT.n_decoys)
    parser.add_argument('-nc', '--n_cycles', metavar='', type=int,
        help='Number of pack/min cycles (Default: %d)' % (DEFAULT.n_cycles),
        default = DEFAULT.n_cycles)
    parser.add_argument('-o', '--output', metavar='', type=str,
        help='Output prefix (Default: %s)' % (DEFAULT.output_prefix),
        default = DEFAULT.output_prefix)
    parser.add_argument('-e', '--export', metavar='', type=str,
        help='Export input details (Default: %s)' % (DEFAULT.export_input),
        default = DEFAULT.export_input)
    parser.add_argument("-d", "--dry", action='store_true',
        help="Don't actually run the simulation, just export the input")

    args = parser.parse_args()
    validate_arguments(args)
    init()

    # Define the Starting Pose:
    pose = pose_from_pdb(args.input_file)
    # Ex. Change the name of the Pose
    # Ex. pose.pdb_info().name("Clamshell")

    # Ex. Visualize the Pose on PyMOL
    # Ex. pmm = PyMOLMover()
    # Ex. pmm.apply(pose)

    # Select the movable residues:
    # For this specific case, the selection for design should contain all the
    # residues surrounding the ligand peptide (defined in the chain B of the
    # PDB). The surrounding residues are considered as the residues with a C⍺-C⍺
    # distance lower than the defined 'cutoff'. NOTE: In
    # NeighborhoodResidueSelector, setting 'include_focus_in_subset' to False
    # removes the target ligand itself from the selected residue subset. In
    # PyRosetta, by default, both the repacking and design settings of a
    # PackerTask are set to True. Therefore, the selection should contain all
    # the residues that are intended to be frozen in place.
    lig =  ChainSelector(args.chain)
    interface_and_lig = NeighborhoodResidueSelector(lig, args.cutoff)
    not_interface_not_lig = NotResidueSelector(interface_and_lig)
    # Ex. Get the actual interface and ligand residue numbers in the Pose:
    # Ex. selection = get_residues_from_subset(interface_and_lig.apply(pose))
    # Ex. Visualize the interface residues on PyMOL (in red, over a black Pose):
    # Ex. color_map = map_int_int()
    # Ex. for index in selection: color_map[index] = XC_red
    # Ex. pmm.send_colors(pose, color_map, XC_black)
    # Ex. Get a PyMOL selection string for the given 'subset' (lig):
    # Ex. get_pymol_selection_from_subset(lig, pose)

    # Define the Score Function:
    score_function = get_fa_scorefxn()

    # Define the Rotamer Packer:
    # Ex. Define the Packer Task:
    # Ex. packer_task = standard_packer_task(pose)

    # Define the Packer Tasks through a TaskFactory:
    task_factory = standard_task_factory()
    # An optinal but recommended Task is to include extra rotamers (away from
    # the regular optimums) for both chi1 (ex1) and chi2 (ex2)
    task_factory.push_back(ExtraRotamers(0, 1, 1)) # ex1
    task_factory.push_back(ExtraRotamers(0, 2, 1)) # ex2
    # Restrict the design to the residues on the interface and ligand
    stop = PreventRepackingRLT()
    task_factory.push_back(OperateOnResidueSubset(stop, not_interface_not_lig))
    # However, on the ligand, perform repacking only
    repack_only = RestrictToRepackingRLT()
    task_factory.push_back(OperateOnResidueSubset(repack_only, lig))

    # Define the PackRotamersMover (with the defined TaskFactory):
    # This is essentially a MonteCarlo object. 
    pack_mover = PackRotamersMover(score_function)
    pack_mover.task_factory(task_factory)

    # Define the Minimizer:
    #   Define the MoveMap (allowing movement in the sidechains only):
    move_map = MoveMap()
    move_map.set_chi(True)
    # Ex. Set a single residue (1) Chi movements to False:
    # Ex. move_map.set_chi(1, False)
    # Ex. Set a range of residues (2 to 10) Chi movements to True:
    # Ex. move_map.set_chi_true_range(2, 10)
    #    Note: Using the range, there's currently no function to set to False.
    # Ex. Visualize the current MoveMap:
    # Ex. move_map.show()
    #    Define the MinMover (with the current MoveMap, Score Function and set
    #    the minimization type to "lbfgs_armijo_nonmonotone"):
    min_mover = MinMover()
    min_mover.movemap(move_map)
    min_mover.score_function(score_function)
    min_mover.min_type('lbfgs_armijo_nonmonotone')
    # Ex. Visualize the current MinMover:
    # Ex. min_mover.show()

    # Define the SequenceMover:
    sequence_mover = SequenceMover()
    sequence_mover.add_mover(pack_mover)
    sequence_mover.add_mover(min_mover)

    # Export all the defined input information to an external file
    if not args.export == None:
        with open(args.export, "w") as export_input:
            export_input_file(args,
                pose,
                [lig, interface_and_lig, not_interface_not_lig],
                score_function,
                sequence_mover,
                f_out = export_input)

    # Run the script, deploying decoys with PyJobDistributor:
    #   (Possible BUG) PyJobDistributor doesn't start when it exists a .pdb or
    #   .in_progress file with the same name as the future output from this
    #   script. Relevant when re-starting the same script without cleaning the
    #   folder.
    if not args.dry:
        job_man = PyJobDistributor(args.output, args.n_decoys, score_function)
        while not job_man.job_complete:
            # This 'while' loop will run 'n_decoys' times, and the second time it
            # runs it should start from the original pose, not the previously
            # designed structure.
            p = Pose(pose)
            # A single run of the script is comprised of 'n_cycles' loops of
            # sidechain packaging > minimization.
            for i in range(args.n_cycles):
                # Ex. Run each step on SequenceMover individually
                # Ex. pack_mover.apply(p)
                # Ex. min_mover.apply(p)
                sequence_mover.apply(p)
            job_man.output_decoy(p)


#             A U X I L I A R Y   F U N C T I O N S
# ______________________________________________________________________________
# 
# Minimalistic version of the above script.
# Aimed to be called from other scripts.

def design(pose, score_function, designable, repackable):

    designer = get_design_mover(pose, score_function, designable, repackable)
    p = Pose(pose)
    designer.apply(p)
    return p

def get_design_mover(score_function, designable, repackable):
    """
    Returns a design Mover.
    The design mover used in this example is a SequenceMover comprised of 2\
    movers: a PackRotamersMover and a MinMover, in this order. For the
    PackRotamersMover, extra rotamers on chi1 and chi2 are enabled. During this
    step, residues on the 'designable' ResidueSelector will be subject to design
    efforts, while residues on the 'repackable' ResidueSelector will only change
    conformation to rotamers of the same aminoacid. During the minimization step
    only the sidechains are allowed to relax.
    """

    blockable    = NotResidueSelector(OrResidueSelector(designable, repackable))
    task_factory = standard_task_factory()
    task_factory.push_back(ExtraRotamers(0, 1, 1)) # ex1
    task_factory.push_back(ExtraRotamers(0, 2, 1)) # ex2
    block        = PreventRepackingRLT()
    repack       = RestrictToRepackingRLT()
    task_factory.push_back(OperateOnResidueSubset(block, blockable))
    task_factory.push_back(OperateOnResidueSubset(repack, repackable))
    
    pack_mover   = PackRotamersMover(score_function)
    pack_mover.task_factory(task_factory)

    move_map     = MoveMap()
    move_map.set_chi(True)
    min_mover    = MinMover()
    min_mover.movemap(move_map)
    min_mover.score_function(score_function)
    min_mover.min_type('lbfgs_armijo_nonmonotone')

    designer     = SequenceMover()
    designer.add_mover(pack_mover)
    designer.add_mover(min_mover)

    return designer