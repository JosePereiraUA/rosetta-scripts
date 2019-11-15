import json
import numpy as np
from pyrosetta.rosetta.core.select import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.core.conformation import Residue
from pyrosetta.rosetta.core.select.residue_selector import *


def get_residues_from_selector(selector, pose):
    """
    Returns a list of Residue instances obtained from the application of the
    given 'selector' on the 'pose'.
    """
    residues = []
    residue_indexes = get_residues_from_subset(selector.apply(pose))
    for residue_index in residue_indexes:
        residues.append(pose.residue(residue_index))
    return residues


def get_pymol_selection_from_selector(selector, pose):
    """
    Returns a selection string to be used on PyMOL with the residues that
    comprise the given 'selector' when applied to the 'pose'.
    """
    # Get pose residue numbers from subset
    residues = get_residues_from_subset(selector.apply(pose))
    # Initiate selection syntax
    selection = "select sele, resi "
    # Get the actual residue index in the PDB
    for index in residues:
        residue_index = pose.pdb_info().pose2pdb(index).split()[0]
        selection += "%s+" % (residue_index)
    return selection[:-1] + " and %s" % (pose.pdb_info().name())


# def get_residues_pdb_format_from_selector(selector, pose):
#     """
#     Returns a list comprised of tuples each with the PDB Index and Chain
#     Identifier (in this order) of residues from a given 'selector' applied to
#     the 'pose'.
#     """
#     residues = get_residues_from_subset(selector.apply(pose))
#     r = []
#     for i in residues:
#         v = pose.pdb_info().pose2pdb(i).split()
#         r.append((int(v[0]), v[1]))
#     return r


def get_residue_energies_from_selector(selector, pose, verbose = True):
    """
    Extract total energy values for each residue in the given 'subset' from the
    pose (with energy previosuly calculated). The subset must be a
    ResidueSelector instance. If verbose is set to True (by
    default), results will be printed.
    """

    # Verify that the pose energy has been calculated
    assert len(pose.scores) != 0,\
        "It seems the given pose energy has not been calculated."
        
    assert issubclass(type(selector), ResidueSelector),\
        "'Selector' should be a ResidueSelector."

    # Get pose residue numbers from subset:
    residues = get_residues_from_subset(selector.apply(pose))
    # Extract energies for the given residues from the pose:
    energies = []
    for index in residues:
        residue_energy = pose.energies().residue_total_energies(index)
        residue_index = int(pose.pdb_info().pose2pdb(index).split()[0])
        energies.append((index, residue_index, residue_energy[total_score]))
    if verbose:
        for _, i, e in energies:
            print("%6d %9.4f" % (i, e))
    return energies


def get_cartesian_coordinates_from_pose(pose):
    """
    Return an (n, 3) matrix, where n is the number of atoms in the pose and 3
    refers to the 3 dimensional X, Y and Z degrees of freedom of each atom.
    """
    coords = np.zeros((pose.total_atoms(), 3))
    index  = 0
    for res_index in range(1, pose.total_residue() + 1):
        residue = pose.residue(res_index)
        for atom_index in range(1, residue.natoms() + 1):
            coords[index] = residue.xyz(atom_index)
            index += 1
    return coords


def get_cartesian_coordinates_from_selector(selector, pose):
    """
    Return an (n, 3) matrix, where n is the number of atoms in the given
    'selector' when applied to the 'pose' and 3 refers to the 3 dimensional X, Y
    and Z degrees of freedom of each atom.
    """

    def count_atoms(_residues):
        count = 0
        for residue in residues:
            count += len(residue.atoms())
        return count
        
    residues   = get_residues_from_selector(selector, pose)
    atom_count = count_atoms(residues)
    coords     = np.zeros((atom_count, 3))
    index      = 0
    for residue in residues:
        for atom_index in range(1, residue.natoms() + 1):
            coords[index] = residue.xyz(atom_index)
            index += 1
    return coords


def get_centroid_coordinates_from_selector(selector, pose):
    coords = get_cartesian_coordinates_from_selector(selector, pose)
    return np.mean(coords, axis=0)


def activate_constraints(score_function, c_weight = 1.0):
    """
    Turn score function weights for constraint related components to the given
    'c_weight' value.
    """
    score_function.set_weight(atom_pair_constraint,  c_weight)
    score_function.set_weight(angle_constraint,      c_weight)
    score_function.set_weight(coordinate_constraint, c_weight)
    score_function.set_weight(dihedral_constraint,   c_weight)
    score_function.set_weight(res_type_constraint,   c_weight)



def load_selections_from_json_file(filename, pose):
    """
    Reads a selections.JSON file and transforms the selection into actual
    ResidueSelectors, if the residues exist. Otherwise, save that selection as
    None.
    """

    with open(filename, "r") as file_in:
        selections = json.load(file_in)
    for selection in selections:
        sel = ResidueIndexSelector(selections[selection])
        try:
            sel.apply(pose)
            selections[selection] = sel
        except RuntimeError:
            selections[selection] = None
    return selections


def save_selections_to_json_file(filename, pose, selectors):
    """
    """
    
    selections = {}
    for selector_index, selector in enumerate(selectors):
        selection = []
        res_index_lst = get_residues_from_subset(selector.apply(pose))
        for res_index in res_index_lst:
            res_pdb_index = int(pose.pdb_info().pose2pdb(res_index).split()[0])
            selection.append(res_pdb_index)
        selections[selector.get_name() + str(selector_index)] = selection
    with open(filename, "w") as file_out:
        json.dump(selections, file_out)


def calc_d_table_between_selectors(selector1, selector2, pose):
    """
    Returns an (n x m) matrix (where n is the length of 'selector1' and m is the
    length of 'selector2'), containing the euclidean distances between all pairs
    of residues from both the selections in the given 'pose'. The atoms that are
    considered for calculating the distances are the "CEN" atoms (if the pose is
    in centroid mode) or the last atom of the sidechain (if the pose is in
    normal mode).
    """

    from ze_utils.pyrosetta_classes import ResidueToAtomMapping

    sc = ResidueToAtomMapping()
    sc_map = sc.cen if pose.is_centroid() else sc.reg
    residues1 = get_residues_from_selector(selector1, pose)
    residues2 = get_residues_from_selector(selector2, pose)
    s1 = len(residues1)
    s2 = len(residues2)

    d_table = np.zeros((s1, s2))
    for i in range(s1):
        res_i = residues1[i]
        atm_i = np.array(res_i.atom(sc_map[res_i.name1()]).xyz())
        for j in range(s2):
            res_j = residues2[j]
            atm_j = np.array(res_j.atom(sc_map[res_j.name1()]).xyz())
            d_table[i, j] = np.linalg.norm(atm_j-atm_i)
    return d_table


def set_ABC_model_fold_tree(pose):
        """
        Sets the ABC model fold tree on a given 'pose'. This model states that:
         1) chains A and B of the pose are a single EDGE on the fold tree;
         2) there is one JUMP between chain B and C;
         3) chain C is a single EDGE on the fold tree;
        """
        
        from pyrosetta import FoldTree
        from pyrosetta.rosetta.core.select.residue_selector import \
            ChainSelector, OrResidueSelector

        chainA = ChainSelector("A")
        chainB = ChainSelector("B")
        chainC = ChainSelector("C")
        chainAB = OrResidueSelector(chainA, chainB)

        # First atom of chain A
        fA = 1
        
        atoms_chainAB = chainAB.apply(pose)
        lB = len([x for x in atoms_chainAB if x == 1])
        atoms_chainC = chainC.apply(pose)
        atom_indexes_chainC = [i for i, x in enumerate(atoms_chainC) if x == 1]
        fC = atom_indexes_chainC[0] + 1
        lC = atom_indexes_chainC[len(atom_indexes_chainC) - 1] + 1
        
        fold_tree = FoldTree()
        fold_tree.add_edge(fA, lB, -1)
        fold_tree.add_edge(lB, fC,  1)
        fold_tree.add_edge(fC, lC, -1)
        if not fold_tree.check_fold_tree(): exit(0)
        pose.fold_tree(fold_tree)