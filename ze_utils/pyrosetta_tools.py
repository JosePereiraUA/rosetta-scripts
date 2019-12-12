import json
import numpy as np
from pyrosetta.rosetta.core.select import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.core.conformation import Residue
from pyrosetta.rosetta.core.select.residue_selector import *

# TODO: Documentation

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


def ID2Rank(pose, ids, as_string = True):
    """
    Returns a selection string (by default) with the rank positions in the pose
    for the given list of residue IDs (from the PDB). This string can be used,
    for example, on instances of ResidueIndexSelector. If 'as_string' is set to
    False, return a list of rank positions instead.
    """

    id2rank = {}
    for rank in range(1, len(pose.residues) + 1):
        id2rank[int(pose.pdb_info().pose2pdb(rank).split()[0])] = rank
    if as_string:
        selection = ""
        for i in ids:
            selection += "%d," % (id2rank[i])
        return selection[:-1]
    else:
        return [id2rank[i] for i in ids]


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


def get_residue_energies_from_selector(selector, pose, verbose = False):
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


def uncap(pose):
    """
    Remove all termini conformations of a pose, switching the atoms to normal
    in-chain conformations.
    """

    from pyrosetta.rosetta.core.select.residue_selector import \
        ResidueIndexSelector

    # Identify all termini in the given pose
    cap_indexes = ""
    for residue in pose.residues:
        if residue.is_terminus():
            cap_indexes += "%d," % (residue.seqpos())
    caps = ResidueIndexSelector(cap_indexes[:-1])

    # Remove termini conformations for the selection in the given pose
    uncap_selection(pose, caps)


def uncap_selection(pose, selection):
    """
    Remove all termini conformations of a selection in the given pose, switching
    the atoms to normal in-chain conformations.
    """

    from pyrosetta.rosetta.protocols.simple_moves import ModifyVariantTypeMover

    # Create the ModifyVariantTypeMover, setting it up to remove the terminal
    # variants in the caps defined in the ResidueIndexSelector
    uncapper = ModifyVariantTypeMover()
    uncapper.set_additional_type_to_remove("LOWER_TERMINUS_VARIANT")
    uncapper.set_additional_type_to_remove("UPPER_TERMINUS_VARIANT")
    uncapper.set_residue_selector(selection)
    uncapper.set_update_polymer_bond_dependent_atoms(True)

    # Apply the ModifyVariantTypeMover to the pose
    uncapper.apply(pose)


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
        for residue in _residues:
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
    """
    """
    
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


def travel_connectivity_map(pose, start):
    """
    Recursively travel the connectivity map, starting on residue 'start', and
    return all residues visited. Assumes a non-branched linear chain.
    """

    def travel(index):
        if index == 0:
            return
        chain.append(index)
        travel(pose.residue(index).connected_residue_at_upper())

    chain = []
    travel(start)

    return chain
        

def smart_fold_tree_cut(fold_tree, pos):
    """
    Creates a cut on the 'fold_tree' as the given 'pos'. This separates the
    existing edge into two unconnected ones, with a jump from residue 1 to the
    residue at 'pos' + 1.
    """

    # Obtain old boundaries to define the new pose
    bl = fold_tree.boundary_left(pos)
    br = fold_tree.boundary_right(pos)

    # Obtain the number of existing jumps. The new jump will be set between
    # residue 1 and the given position + 1 residue
    n  = fold_tree.num_jump() + 1

    # Delete the old edge where the new cut will be placed
    fold_tree.delete_edge(fold_tree.get_residue_edge(pos))

    # Create the new 2 edges (one on each side of the cut) and the jump
    fold_tree.add_edge(    bl,      pos, -1)
    fold_tree.add_edge(      1, pos + 1,  n)
    fold_tree.add_edge(pos + 1,      br, -1)


def set_ABC_model_fold_tree(pose):
        """
        Sets the ABC model fold tree on a given 'pose'. This model states that:
         1) All chains of the pose are a single EDGE on the fold tree;
         2) there is two JUMPs between residue 1 and the first residue of both
         chain B and C;

        Note: If the model PDB has been correctly produced, this should be the
        default fold_tree inferred by pyrosetta. This function additionally
        performs several checks to validate the model. 
        """
        
        from pyrosetta import FoldTree
        from pyrosetta.rosetta.core.select.residue_selector import \
            ChainSelector, OrResidueSelector

        # Verify is the ABC Model is being respected:
        # 1) Should have 3 and only 3 chains
        assert pose.num_chains() == 3, \
            "ABC Model in pose should have 3 chains: A, B and C."

        # 2) All chains should be continuous (have no intra-chain breaks)
        def verify_chain_coherence(chain_index):
            """
            Asserts if it is possible to travel the connectivity map of a chain
            from the first residue and reach the end of the chain, therefore
            lacking breaks.
            """

            chain_res_indexes = travel_connectivity_map(
                pose, pose.chain_begin(chain_index))
            for idx in chain_res_indexes:
                pdb_idx = int(pose.pdb_info().pose2pdb(idx).split()[0])
                assert pdb_idx == idx, \
                    "Chain %s seems to have intra-chain break at residue %d" % \
                        (pose.pdb_info().chain(chain_index), pdb_idx)

        for chain_index in range(1, pose.num_chains() + 1):
            verify_chain_coherence(chain_index)
        
        # 3) The 3 chains should be A, B and C, in order
        assert pose.pdb_info().chain(pose.chain_begin(1)) == 'A', \
            "First chain should be A."
        assert pose.pdb_info().chain(pose.chain_begin(2)) == 'B', \
            "Second chain should be B."
        assert pose.pdb_info().chain(pose.chain_begin(3)) == 'C', \
            "Third chain should be C."

        # Create and empty fold tree
        fold_tree = FoldTree()

        # Set the fold tree to be one constinuous edge across the whole pose
        fold_tree.simple_tree(len(pose.residues))

        # Cut the single edge separating each of the defined chains
        smart_fold_tree_cut(fold_tree, pose.chain_end(1))
        smart_fold_tree_cut(fold_tree, pose.chain_end(2))

        # Apply the new fold tree to the pose
        pose.fold_tree(fold_tree)


def load_pre_filter(filename = "auto"):
    """
    Reads an input JSON 'filename' with all the parameters of a PreFilter and
    returns an instance of PreFilter. In the 'filename' is set to "auto",
    returns a default PreFilter instead.

    To create a JSON file with the parameters of a customly defined PreFilter
    check 'save_pre_filter' function.
    """

    from ze_utils.pyrosetta_classes import PreFilter

    if filename != "auto":
        with open(filename, "r") as pre_filter_config_json:
            pf = json.load(pre_filter_config_json)
        return PreFilter(pf["contact_cutoff"], pf["contact_min_count"],
            pf["anchors_cutoff"], pf["clash_cutoff"], pf["clashes_max_count"],
            pf["sel_A"], pf["sel_B"], pf["sel_C"],
            pf["prevent_block_C_terminal"], pf["prevent_block_N_terminal"],
            pf["max_c_terminal_interaction"], pf["max_n_terminal_interaction"],
            pf["score_function"])
    else:
        return PreFilter()


def save_pre_filter(pre_filter, filename):
    """
    Saves the input 'pre_filter' as a JSON dictionary in 'filename'. However,
    the PreFilter sel_A, sel_B, selC and score_function parameters are ignored
    and set as "auto".

    To load a JSON pre filter file as a PreFilter instance, check
    load_pre_filter function.
    """

    from ze_utils.pyrosetta_classes import PreFilter

    assert type(pre_filter) == PreFilter, \
        "Failed to dump PreFilter. 'pre_filter' needs to be of type PreFilter."

    with open(filename, "w") as pf_json:
        data = pre_filter.__dict__
        data["sel_A"]          = "auto"
        data["sel_B"]          = "auto"
        data["sel_C"]          = "auto"
        data["score_function"] = "auto"
        json.dump(data, pf_json)