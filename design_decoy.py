# Questions: jose.manuel.pereira@ua.pt

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


#                         Single Decoy Design Script:
# ______________________________________________________________________________
#  Perform a sequence design by mutating the nature of the aminoacids in the
# structure. PyRosetta does this by changing the sidechain and its rotamers.
# The design effort will be localized in the surrounding aminoacids of a ligand
# peptide, defined as the chain C (by default) of the input PDB structure. If a
# blocked chain ID is provided, no design or repack will be performed on
# residues of that chain.
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
#
#    THIS SCRIPT PERFORMS ONLY ONE DECOY OF THE DESIGN PROTOCOL.
#  > Check design.py script for multiple decoy design simulation
# ______________________________________________________________________________

class DEFAULT:
    """
    Define defaults for the script. Values can be modified using arguments.
    """
    ligand_chain   = "C"
    blocked_region = None
    cutoff         = 9.0
    n_cycles       = 4
    output_prefix  = "design"
    score_function = "auto"
    designable     = "auto"
    repackable     = "auto"
    sele_file      = None


def design(pose, score_function = DEFAULT.score_function,
    designable = DEFAULT.designable, repackable = DEFAULT.repackable,
    cutoff = DEFAULT.cutoff, ligand_chain = DEFAULT.ligand_chain,
    blocked_region = DEFAULT.blocked_region, n_cycles = DEFAULT.n_cycles,
    sele_file = DEFAULT.sele_file):

    designer = get_designer_mover(score_function, designable, repackable,
        cutoff, ligand_chain, blocked_region)

    # Print the designed region to a file, if sele_file is defined and the file
    # doesn't exist
    if sele_file != None and not os.path.exists(sele_file):
        with open(sele_file, "w") as file_out:
            design_reg = get_designable_region(designable, cutoff, ligand_chain,
                blocked_region)
            file_out.write(get_pymol_selection_from_selector(design_reg, pose))

    p = Pose(pose)
    for i in range(n_cycles):
        designer.apply(p)
    return p


def get_designable_region(designable = DEFAULT.designable,
    cutoff = DEFAULT.cutoff, ligand_chain = DEFAULT.ligand_chain,
    blocked_region = DEFAULT.blocked_region):
    """
    Return the designable region.
    """

    if designable == "auto":
        if blocked_region == None:
            designable = NeighborhoodResidueSelector(
                ChainSelector(ligand_chain),
                cutoff,
                include_focus_in_subset = False)
        else:
            print(blocked_region)
            designable = AndResidueSelector(
                NeighborhoodResidueSelector(
                    ChainSelector(ligand_chain),
                    cutoff,
                    include_focus_in_subset = False),
                NotResidueSelector(ResidueIndexSelector(blocked_region)))
    else:
        assert isinstance(designable, ResidueSelector) or designable == "auto",\
        "Designable selection must be a ResidueSelector or set to 'auto'"
    
    return designable


def get_designer_mover(score_function = DEFAULT.score_function,
    designable = DEFAULT.designable, repackable = DEFAULT.repackable,
    cutoff = DEFAULT.cutoff, ligand_chain = DEFAULT.ligand_chain,
    blocked_region = DEFAULT.blocked_region):
    """
    Returns a design Mover.
    The design mover used in this example is a SequenceMover comprised of 2
    movers: a PackRotamersMover and a MinMover, in this order. For the
    PackRotamersMover, extra rotamers on chi1 and chi2 are enabled. During this
    step, residues on the 'designable' ResidueSelector will be subject to design
    efforts, while residues on the 'repackable' ResidueSelector will only change
    conformation to rotamers of the same aminoacid. During the minimization step
    only the sidechains are allowed to relax. By default, when no custom
    designable or repackable ResidueSelector's are provided, both designable and
    repackable regions are set to 'auto', where the repackable region is defined
    as all the residues of the ligand_chain (C, by default) and the designable
    region is defined as the residues within a cutoff (9.0 Angstrom, by default)
    from the ligand_chain. If a blocked region is provided, no design nor repack
    will be performed on those residues. A blocked region is defined with the
    following syntax: "start-end", where start and end are residue index
    numbers. For example, blocking the region from residue 1 to 56, one would
    use "1-56".
    """

    # --- SCORE FUNCTION
    if score_function == "auto":
        try:
            score_function = get_fa_scorefxn()
        except:
            score_function = get_fa_scorefxn()
    else:
        from pyrosetta.rosetta.core.scoring import ScoreFunction
        assert type(score_function) == ScoreFunction, \
            "Score function for relaxer mover must be of type ScoreFunction."

    # --- DESIGNABLE REGION
    designable = get_designable_region(designable, cutoff, ligand_chain,
        blocked_region)

    # --- REPACKABLE REGION
    if repackable == "auto":
        repackable = ChainSelector(ligand_chain)
    else:
        assert isinstance(repackable, ResidueSelector) or repackable == "auto",\
            "Repackable selection must be a ResidueSelector or set to 'auto'"

    task_factory = standard_task_factory()
    blockable    = NotResidueSelector(OrResidueSelector(designable, repackable))
    task_factory.push_back(ExtraRotamers(0, 1, 1)) # ex1
    task_factory.push_back(ExtraRotamers(0, 2, 1)) # ex2
    block        = PreventRepackingRLT() # NO design, NO repacking
    repack       = RestrictToRepackingRLT() # NO design, ONLY repacking
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


def validate_arguments(args):
    """
    Validate the script arguments to prevent known input errors.
    """
    if not args.input_file[-4:] == ".pdb":
        exit("ERROR: Input file should be in PDB format")
    if args.cutoff < 0:
        exit("ERROR: Cut-off must be a non-negative value")
    if args.n_cycles < 0:
        exit("ERROR: Number of cycles must be a non-negative value")
    if args.blocked_region != None:
        assert len(args.blocked_region.split("-")) == 2, \
            "Blocked region should be in format 'int-int'"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Design the active-site of a
        protein surrousing a ligand defined in another chain, on the input
        PDB. If a blocked region is provided, no design nor repack will be
        performed on those residues. A blocked region is defined with the
        following syntax: "start-end", where start and end are residue index
        numbers. For example, blocking the region from residue 1 to 56, one
        would use "-bl 1-56". A total of n_cycles will be attempted back to back
        to increase convergence. If a sele_file is provided, the region
        determined for design will be printed to that file.""")
    parser.add_argument('input_file', metavar='INPUT', type=str,
        help='The input PDB file name')
    parser.add_argument('-lc', '--ligand_chain', metavar='', type=str,
        help='The ligand chain identifier (Default: %s)' % \
        (DEFAULT.ligand_chain), default = DEFAULT.ligand_chain)
    parser.add_argument('-b', '--blocked_region', metavar='', type=str,
        help='The blocked (no design) region (Default: %s)' % \
        (DEFAULT.blocked_region), default = DEFAULT.blocked_region)
    parser.add_argument('-co', '--cutoff', metavar='', type=float,
        help='Interface identification cut-off in Å (Default: %5.1f)' \
            % (DEFAULT.cutoff), default = DEFAULT.cutoff)
    parser.add_argument('-nc', '--n_cycles', metavar='', type=int,
        help='Number of pack/min cycles (Default: %d)' % (DEFAULT.n_cycles),
        default = DEFAULT.n_cycles)
    parser.add_argument('-o', '--output', metavar='', type=str,
        help='Output prefix (Default: %s)' % (DEFAULT.output_prefix),
        default = DEFAULT.output_prefix)
    parser.add_argument('-sele', '--sele_file', metavar='', type=str,
        help='Selection output file name (Default: %s)' % (DEFAULT.sele_file),
        default = DEFAULT.sele_file)

    args = parser.parse_args()
    validate_arguments(args)
    init()

    # Define the Starting Pose:
    pose = pose_from_pdb(args.input_file)
    # Optional but recommended for looking for NeighborhoodResidueSelector
    pose.update_residue_neighbors() 
    # Ex. Change the name of the Pose
    # Ex. pose.pdb_info().name("Clamshell")

    # Ex. Visualize the Pose on PyMOL
    # Ex. pmm = PyMOLMover()
    # Ex. pmm.apply(pose)

    # Select the non movable residues:
    # For this specific case, the selection for design should contain all the
    # residues surrounding the ligand peptide (defined in the chain C of the
    # PDB, by default). The surrounding residues are considered as the residues
    # with a C⍺-C⍺ distance lower than the defined 'cutoff'. Note: In
    # NeighborhoodResidueSelector, setting 'include_focus_in_subset' to False
    # removes the target ligand itself from the selected residue subset. In
    # PyRosetta, by default, both the repacking and design settings of a
    # PackerTask are set to True. Therefore, the selection should contain all
    # the residues that are intended to be frozen in place. During this process,
    # the designable (no alteration, design is the default in PyRosetta) and
    # repackable (use RestrictToRepackingRLT operation) regions will be defined.
    # If a blocked_region is provided, all residues within that region are
    # considered blocked and no design nor repack will be performed on them.
    designable = get_designable_region("auto", args.cutoff, args.ligand_chain,
        args.blocked_region)

    repackable = ChainSelector(args.ligand_chain)

    # Print the designed region to a file, if sele_file is defined
    if args.sele_file != None:
        with open(args.sele_file, "w") as file_out:
            file_out.write(get_pymol_selection_from_selector(designable, pose))

    # --- Extras: --------------------------------------------------------------
    # Note: Functions marked with an asterisk * require ze_utils.pyrosetta tools
    # Ex. Get the actual designable and repackable residue numbers in the Pose:
    # Ex. selection = get_residues_from_subset(designable.apply(pose))
    # Ex. Get the actual designable and repackable residue objects in the Pose:
    # Ex. selection = get_residues_from_selector(designable, pose) *
    # Ex. Visualize the interface residues on PyMOL (in red, over a black Pose):
    # Ex. color_map = map_int_int()
    # Ex. for index in selection: color_map[index] = XC_red
    # Ex. pmm.send_colors(pose, color_map, XC_black)
    # Ex. Get a PyMOL selection string for the given 'subset' (lig):
    # Ex. get_pymol_selection_from_subset(lig, pose) *

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
    stop         = PreventRepackingRLT()
    blockable    = NotResidueSelector(OrResidueSelector(designable, repackable))
    task_factory.push_back(OperateOnResidueSubset(stop, blockable))
    # However, on the ligand, perform repacking only
    repack_only  = RestrictToRepackingRLT()
    task_factory.push_back(OperateOnResidueSubset(repack_only, repackable))

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

    # Run the script
    for i in range(args.n_cycles):
        # Ex. Run each step on SequenceMover individually
        # Ex. pack_mover.apply(p)
        # Ex. min_mover.apply(p)
        sequence_mover.apply(pose)
    pose.dump_pdb(args.output + ".pdb")