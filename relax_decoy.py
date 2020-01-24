# Questions: jose.manuel.pereira@ua.pt

import argparse
from pyrosetta import *
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.pack.task.operation import *
from ze_utils.pyrosetta_tools import activate_constraints
from pyrosetta.rosetta.protocols.constraint_generator import \
	AddConstraints, CoordinateConstraintGenerator


#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\

#                       Single Decoy Relax Script:
# ______________________________________________________________________________
#  Performs an energy minimization decoy.
#  
#  Necessary parts:
# - Starting Pose;
# - Score function;
# - FastRelax (sucessive rounds of backbone and side-chain packing and
#   minimization, where the repulsive weight in the scoring function is
#   gradually increased. )
#
#  Optinal parts:
# - Constraints (Add energy components that penalize certain conformations based
#   on the defiden constraints. These can be, for example, restraints that lock 
#   the backbone to the original position, try to approximate two residues, etc)
#   In order for a Constraint to be correctly applied in a simulation, two steps
#   are necessary:
#   . Set the constraints in the pose (see bellow);
#   . Enable all the constraint weights in the score function;
#   There are two distinct constraint objects used to set constraints in a pose:
#   . ConstraintMovers: Are the actual constraints. Can be applied to a pose.
#   . ConstraintGenerators: Automatically generate ConstraintMovers, regardless
#     of the pose they are applied to. Are usually more useful. In order to
#     apply automatically generated constraints using a ConstraintGenerator, an
#     AddConstraint object must be initialized, to which an arbitrary number of
#     generators can be added and then all applied at once to all the desired
#     poses. 
# - TaskFactory (Used to enable extra rotamers to be sampled in the packing
#   process. Is applied to the FastRelax object).
#
#    THIS SCRIPT PERFORMS ONLY ONE DECOY ON THE RELAX PROTOCOL.
#  > Check relax.py script for multiple decoy relax simulation
# ______________________________________________________________________________

class DEFAULT:
    """
    Define defaults for the script. Values can be modified using arguments.
    """
    c_weight        = 1.0
    no_constraints = False
    output_prefix  = "relax"
    score_function = "auto"


def validate_arguments(args):
    """
    Validate the script arguments to prevent known input errors.
    """
    if not args.input_file[-4:] == ".pdb":
        exit("ERROR: Input file should be in PDB format")
    if args.c_weight < 0 or args.c_weight > 1.0:
        exit("ERROR: Constraints weight must be a non-negative value (< 1.0)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Relax a PDB structure.""")
    parser.add_argument('input_file', metavar='INPUT', type=str,
        help='The input PDB file')
    parser.add_argument('-nc', '--no_constraints',
        action='store_true', help='Turn off constraints (Default: False)')
    parser.add_argument('-w', '--c_weight', metavar='', type=float  ,
        help='Constraints weight (Default: %5.2f)' % (DEFAULT.c_weight),
        default = DEFAULT.c_weight)
    parser.add_argument('-o', '--output', metavar='', type=str,
        help='Output prefix (Default: %s)' % (DEFAULT.output_prefix),
        default = DEFAULT.output_prefix)

    args = parser.parse_args()
    validate_arguments(args)
    init()

    # Define the starting Pose:
    pose = pose_from_pdb(args.input_file)

    # Define the Score Function:
    score_function = get_fa_scorefxn()

    if args.no_constraints == False:
        # Create the Generator
        coordinates_constraint = CoordinateConstraintGenerator()
        # Add the Generate to a new AddConstraints object
        constraints = AddConstraints()
        constraints.add_generator(coordinates_constraint)
        # Apply the defined AddConstraints object to a Pose
        constraints.apply(pose)
        # Enable all constraint weights on the score function
        activate_constraints(score_function, args.c_weight)

    # An optinal but recommended Task is to include extra rotamers (away from
    # the regular optimums) for both chi1 (ex1) and chi2 (ex2)
    task_factory = standard_task_factory()
    # It seems that Ex1 is enabled by default.
    # Ex. task_factory.push_back(ExtraRotamers(0, 1, 1)) # ex1
    task_factory.push_back(ExtraRotamers(0, 2, 1)) # ex2
    # FastRelax, by default, blocks 'design' on the regular rotamers. However,
    # when adding the ExtraRotamers, new 'designable' rotamers are included
    # and considered during the sampling. Therefore, the TaskFactory must
    # re-restrict the rotamers to repackaging.
    task_factory.push_back(RestrictToRepacking())
    # Altough the PyRosetta manual states that the current rotamer is included
    # by default, this does not seem to be true and therefore needs to be
    # activated.
    task_factory.push_back(IncludeCurrent())

    # Define the FastRelax object with the correct score function and extra
    # rotamers enabled
    fast_relax = FastRelax()
    fast_relax.set_scorefxn(score_function)
    fast_relax.set_task_factory(task_factory)

    fast_relax.apply(pose)
    pose.dump_pdb(args.output + ".pdb")


#             A U X I L I A R Y   F U N C T I O N S
# ______________________________________________________________________________
# 
# Minimalistic version of the above script.
# Aimed to be called from other scripts.

def relax(pose, no_constraints = DEFAULT.no_constraints,
    c_weight = DEFAULT.c_weight):
    """
    TO DO
    """

    relaxer = get_relaxer_mover(pose, no_constraints = no_constraints,
        c_weight = c_weight)
    p = Pose(pose)
    relaxer.apply(p)
    return p


def get_relaxer_mover(pose, score_function = DEFAULT.score_function,
    no_constraints = DEFAULT.no_constraints, c_weight = DEFAULT.c_weight):
    """
    TO DO
    """

    if score_function == "auto":
        try:
            score_function = get_fa_scorefxn()
        except:
            import pyrosetta
            init()
            score_function = get_fa_scorefxn()
    else:
        from pyrosetta.rosetta.core.scoring import ScoreFunction
        assert type(score_function) == ScoreFunction, \
            "Score function for relaxer mover must be of type ScoreFunction."

    if no_constraints == False:
        coordinates_constraint = CoordinateConstraintGenerator()
        constraints = AddConstraints()
        constraints.add_generator(coordinates_constraint)
        constraints.apply(pose)
        activate_constraints(score_function, c_weight)
    task_factory = standard_task_factory()
    task_factory.push_back(ExtraRotamers(0, 2, 1))
    task_factory.push_back(RestrictToRepacking())
    task_factory.push_back(IncludeCurrent())
    fast_relax = FastRelax()
    fast_relax.set_scorefxn(score_function)
    fast_relax.set_task_factory(task_factory)

    return fast_relax