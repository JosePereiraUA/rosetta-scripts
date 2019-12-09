#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\

#  The objetive of this script is to provide a way for structure alignment
# inside the PyRosetta framework. This alignment can be done over all the atoms
# in the structure or a selected subset, given by a mask. Read bellow for more
# information.

import numpy as np

def as_matrix(pose, mask = None):
    """
    Return an (n, 3) matrix, where n is the number of atoms in the pose and 3
    refers to the 3 dimensional X, Y and Z degrees of freedom of each atom.
    """
    
    coords = []
    index = 0
    for res_index in range(1, pose.total_residue() + 1):
        residue = pose.residue(res_index)
        for atom_index in range(1, residue.natoms() + 1):
            if residue.type().is_virtual(atom_index):
                continue
            if mask == None or mask[index] == 1:
                coords.append(residue.xyz(atom_index))
            index += 1
    return np.array(coords)


def apply_coordinates(pose, coords):
    """
    Apply a matrix of coordinates 'coords' to this pose.
    """

    from pyrosetta.rosetta.numeric import xyzVector_double_t as c_vector

    i = 0
    for res_index in range(1, pose.total_residue() + 1):
        residue = pose.residue(res_index)
        for atom_index in range(1, residue.natoms() + 1):
            if residue.type().is_virtual(atom_index):
                continue
            atom = residue.atom(atom_index)
            atom.xyz(c_vector(coords[i][0], coords[i][1], coords[i][2]))
            i += 1


def rmsd(pose, reference_pose, movable_mask = None, reference_mask = None):
    """
    Calculate the RMSD between this pose and a reference pose. Also
    returns the number of atoms considered for this calculation. If movable
    and references masks are provided, only consider the atoms whose index
    in the corresponding mask is set to True.
    """

    movable_coords   = as_matrix(pose, movable_mask)
    reference_coords = as_matrix(reference_pose, reference_mask)

    assert len(movable_coords) == len(reference_coords), \
        "Movable (%d) and Reference (%d) number of atoms do not match" % \
        (len(movable_coords), len(reference_coords))

    n_atoms          = len(reference_coords)
    distance         = reference_coords - movable_coords

    return n_atoms, np.sqrt(np.sum(distance * distance) / n_atoms)


def align(pose, reference_pose, movable_mask = None, reference_mask = None,
    verbose = True):
    """
    Align this pose with a reference pose. If movable and reference
    masks are provided (Default: None), only the atoms whose index in the
    corresponding mask is set to True will be considered for alignment. Both
    sets of masked atoms must be on the same size. Verbose flag determines
    if the RMSD value is computed and printed (Default: True).
    """
        
    # 1) Get coordinates as a matrix
    movable_coords    = as_matrix(pose, movable_mask)
    reference_coords  = as_matrix(reference_pose, reference_mask)

    assert len(movable_coords) == len(reference_coords), \
        "Movable (%d) and Reference (%d) number of atoms do not match" % \
        (len(movable_coords), len(reference_coords))
    
    # 2) Center on centroid
    # n                = len(pose.sequence())
    c_movable        = np.mean(movable_coords,   axis = 0)
    c_reference      = np.mean(reference_coords, axis = 0)
    movable_coords   = movable_coords   - c_movable
    reference_coords = reference_coords - c_reference
    
    # 3) Obtain correlation matrix
    cm = np.dot(np.transpose(movable_coords), reference_coords)
    u, d, vt = np.linalg.svd(cm)
    
    # 4) Obtain rotation
    rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
    
    # 5) Apply rotation + translation
    transformed_coords  = np.dot(as_matrix(pose) - c_movable, rot)
    transformed_coords += c_reference
    
    # 6) Apply transformed coordinates
    apply_coordinates(pose, transformed_coords)
    
    # 7) Calculate RMSD (Optional)
    count, rms = rmsd(pose, reference_pose, movable_mask, reference_mask)
    if verbose:
        print("RMSD: %6.3f angstrom (%3d atoms)" % (rms, count))
    
    return rms


def get_mask_from_c_alphas(pose):
    """
    Return a Boolean list with an entry for each non virtual atom in the pose.
    Set the entries whose index correspond to an alpha carbon atom in the pose
    to True.
    """

    mask = []
    for residue in pose.residues:
        for atom_index in range(1, residue.natoms() + 1):
            if residue.type().is_virtual(atom_index):
                continue
            if residue.type().atom_name(atom_index).strip() == 'CA':
                mask.append(1)
            else:
                mask.append(0)

    return mask


def get_empty_mask(pose):
    """
    Return a Boolean list with an entry for each non virtual atom in the pose.
    Set all entries to False.
    """

    mask = []
    for residue in pose.residues:
        for atom_index in range(1, residue.natoms() + 1):
            if not residue.type().is_virtual(atom_index):
                mask.append(0)

    return mask


def set_mask_true_start_stop(pose, mask, start, stop):
    """
    Set the continuous region of atoms between residue 'start' (inclusive) and
    residue 'stop' (inclusive) to True, in the provided mask. 
    """

    def set_after_residue_index(res_index, target):

        index = 0
        for residue in pose.residues:
            for atom_index in range(1, residue.natoms() + 1):
                if residue.type().is_virtual(atom_index):
                    continue
                if residue.seqpos() >= res_index:
                    mask[index] = target
                index += 1
        
    set_after_residue_index(start, 1)
    set_after_residue_index(stop + 1, 0)
    
    return mask


# ----------------------------------------------------
from pyrosetta import *
init()

# A mask is an array of 0 and 1's. Only the atoms whose index in the mask is set
# to 1 will be considered for the alignment. Moreover, by not providing with a
# mask, the alignment function will consider all atoms, as long as both poses
# share the same number of atoms.
#
# Create a mask for 2 given regions on the pose (between residues A and B):
# mask = get_empty_mask(pose)
# mask = set_mask_true_start_stop(pose, mask, 3,  5)
# mask = set_mask_true_start_stop(pose, mask, 56, 178)

pose           = pose_from_pdb("candidate.pdb")
mask           = get_mask_from_c_alphas(pose)

reference_pose = pose_from_pdb("3ch8_3p_rlx.pdb")
reference_mask = get_mask_from_c_alphas(reference_pose)

align(pose, reference_pose, mask, reference_mask)
pose.dump_pdb("alignment_test.pdb")
# ----------------------------------------------------