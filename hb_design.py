# Questions: jose.manuel.pereira@ua.pt

import numpy as np
from pyrosetta import *
from operator import mul
from functools import reduce
from ze_utils.pyrosetta_classes import \
    RotamerLibrary, RotamerCycler, DesignHydrogenBondMover
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.select.residue_selector import *
from ze_utils.pyrosetta_tools import \
    calc_d_table_between_selectors_atoms, calc_d_table_between_selectors, \
    ID2Rank, measure_affinity, measure_hbonds_number, get_pymol_selection_from_selector
from pyrosetta.rosetta.core.scoring.hbonds import HBondSet, fill_hbond_set
from ze_utils.common import progress_bar
from design_decoy import get_designer_mover


#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2020 \\

#                  Hydrogen bond focused design Script:
# ______________________________________________________________________________
#  Performs a sequence design. Besides regular design (as described in design.py
# script), performs an initial pre-step where all the rotamers from the 11
# hydrogen bond forming aminoacids are exhaustively explored by enumeration to
# find at least 1 hydrogen bond in order to "seed" the regular design effort, as
# the hydrogen bond forming rotamers are locked in place and all designable
# residues close-by are designed to adapt and estabilize this conformation.
#
#  Hydrogen bonds are identified using HBondSet. HBondSet stores information
# regarding the number, distances, angles and energy values of hydrogens bonds
# in a pose. This object needs to be populated from a pose with:
#  1. It's energy calculated.
#  2. It's neighbours residue list updated.
# Ex. hbond_set = HBondSet()
# Ex. score_function(pose)
# Ex. pose.update_residue_neighbors()
# Ex. fill_hbond_set(pose, False, hbond_set)
# Show the full HBondSet
# Ex. hbond_set.show(pose)
# Show the HBond objects list for a single residue (Ex: residue 10)
# Ex. hbonds_res_10 = hbond_set.residue_hbonds(10)
# A single HBond object contains the distance, angle and energy value of that 
# specific hydrogens bond. This can be visualized:
# Ex. print(hbonds_res_10[1])
# Ex. 49 don: protein 11 23 acc: protein 10 8 -1.36506 0.976471
# Note: In the hbonds objects list, indexing starts at 1, not 0.
# The last number is the weight given to the hydrogen bond (based on distance
# and angle), and the second to last number is the energy value. The don
# (i.e. donor) and acc (i.e. acceptor) identification is given by the residue
# number, followed by the atom number intra-residue.
# This data can be accessed on the object itself.
# Ex. hbonds_res_10[1].energy()
# Ex. hbonds_res_10[1].weight()
#
# This script uses a custom created RotamerCycler to read a rotamer library and
# cycle through all the rotamers. However, a native mover exists for this
# purpose, altough at this time it doesn't seem to be working. Nevertheless,
# some of the gathered information about this mover is documented bellow.
# TryRotamers Mover allows the enumeration of all rotamers of a single
# aminoacid.
#
# TryRotamers arguments:
#  1 - res_num         - (Str ) Residue number to suffer rotamer change
#  2 - score_function  - (SFXN) Score Function employed /?
#  3 - explosion       - (Int ) Until what chi should the rotamer change
#  4 - jump_num        - (Int ) Jump number considered for energy calculation /?
#  5 - clach_check     - (Bool) Filter out clashing rotamers
#  6 - solor_res       - (Bool) Include terminal rotamers
#  7 - include_current - (Bool) Include current rotamers on the search
#
# Each 'apply' call of a TryRotamers object will cycle the res_num to the next
# rotamer on that pose, until all rotamers have been visited. Current status of
# the object can be checked using get_last_move_status function. It will return:
# . MoveStatus.MS_SUCCESS if additional rotamers are available.
# . MoveStatus.FAILT_DO_NOT_RETRY if no additional rotamers are available.
# Ex: next_rotamer = TryRotamers("10", score_function, 4, 0, False, False, True)

if __name__ == "__main__":

    init()

    pmm   = PyMOLMover()
    pmm.keep_history(True)
    pose   = pose_from_pdb("relax_570.pdb")
    pmm.apply(pose)
    rotlib = RotamerLibrary("~/Desktop/script/static/rotlib.dat")
    des    = AndResidueSelector(
                NeighborhoodResidueSelector(
                    ResidueIndexSelector("201"), 9.0, False), 
                NotResidueSelector(
                    OrResidueSelector(
                        ResidueIndexSelector("1-100"),
                        ChainSelector("B"))))
    print(get_pymol_selection_from_selector(des, pose))
    test   = DesignHydrogenBondMover(201, rotlib, 10, designable = des)
    pose   = test.apply(pose)
    pmm.apply(pose)

def get_designer_mover(target_aa, max_length = 5, n_cycles = 4, rl = "auto"):
    if rl == "auto":
        rl = os.path.join(os.path.expanduser("~"),
            "Desktop/scripts/static/rotlib.dat")
    rotlib = RotamerLibrary(rl)
    return DesignHydrogenBondMover(target_aa, rotlib, max_length, n_cycles)