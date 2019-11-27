# import numpy as np
# from ze_utils.molecule_manipulation import *
# from ze_utils.algebra import \
#     get_rotation_matrix_from_axis_angle, rotate_coords_from_rotation_matrix

# reference = Molecule("3ch8_3p_rlx.pdb")
# loop      = Molecule("placeholder10.pdb")

# idx = [atom for atom in reference.atoms if atom.chain_name == "A"][-1].res_index
# reference_mask = reference.get_mask_for_residue_indexes([idx])
# loop_mask = loop.get_mask_for_residue_indexes([1])

# masked_atoms = []
# first = True
# atoms, atom_indexes = loop.get_residue_atoms_from_indexes([1])
# for index, atom in enumerate(atoms):
#     if atom.elem == "H":
#         if first:
#             first = False
#             continue
#         masked_atoms.append(index)
# for masked_atom in masked_atoms:
#     loop_mask[atom_indexes[masked_atom]] = False

# loop.align(reference, loop_mask, reference_mask)
# loop.export("teste.pdb")

from pyrosetta import *
from pyrosetta.rosetta.core.chemical import *
# from pyrosetta.rosetta.core.io.pdb import dump_pdb
# from pyrosetta.rosetta.core.io import StructFileRepOptions
from pyrosetta.rosetta.core.select.residue_selector import *
# from pyrosetta.rosetta.protocols.grafting import delete_region
from ze_utils.pyrosetta_tools import get_residues_from_selector
# from pyrosetta.rosetta.core.pose import declare_cutpoint_chemical_bond
# from pyrosetta.rosetta.protocols.loops import Loop, add_single_cutpoint_variant
init()

# def reload_pose_from_file(pose):

#     import os

#     temporary_number = 0
#     while True:
#         filename = "tmp%d.pdb" % (temporary_number)
#         if not os.path.exists(filename):
#             pose.dump_pdb(filename)
#             new_pose = pose_from_pdb(filename)
#             with open(filename, "a") as file_out:
#                 file_out.write("CONECT   90   99\n")
#             os.system("rm -f %s" % (filename))
#             return new_pose
#         else:
#             temporary_number += 1

from pyrosetta.rosetta.protocols.simple_moves import ModifyVariantTypeMover
from pyrosetta.rosetta.utility import vector1_unsigned_long
pose = pose_from_pdb("candidate.pdb")

cap_indexes = vector1_unsigned_long()
cap_indexes.append(90)
caps = ResidueIndexSelector(cap_indexes)

uncapper = ModifyVariantTypeMover()
uncapper.set_additional_type_to_remove("LOWER_TERMINUS_VARIANT")
uncapper.set_additional_type_to_remove("UPPER_TERMINUS_VARIANT")
uncapper.set_residue_selector(caps)
uncapper.set_update_polymer_bond_dependent_atoms(True)

uncapper.apply(pose)
pose.dump_pdb("test.pdb")


exit(1)
# ----------------------------------------

pose2 = pose_from_pdb("candidate.pdb")
residue_name_to_recover = pose.residue(90).name1()
t = pose.phi(91)
# print("PHI 91: %f" % (pose.phi(91)))

# print("PHI 89: %f" % (pose.phi(89)))
# print("PSI 89: %f" % (pose.psi(89)))
# print("OMEGA 89: %f" % (pose.omega(89)))
phi_to_recover = pose.phi(90)
psi_to_recover = pose.psi(90)
omega_to_recover = pose.omega(90)

delete_region(pose, 90, 90)

residue_provider = pose_from_sequence("A" + residue_name_to_recover + "A")
# residue_provider.set_phi(  2, phi_to_recover)
# residue_provider.set_psi(  2, psi_to_recover)
# residue_provider.set_omega(2, omega_to_recover)

residue_to_recover = residue_provider.residue(2)
pose.prepend_polymer_residue_before_seqpos(residue_to_recover, 90, True)
# pose.set_phi(  91, t)

print(pose2.fold_tree())
# pose.fold_tree(FoldTree(pose2.fold_tree()))
fold_tree = FoldTree()
fold_tree.add_edge(1, 81,  -1)
fold_tree.add_edge(1, 82,   1)
fold_tree.add_edge(82, 89, -1)
fold_tree.add_edge(1, 92,   2)
fold_tree.add_edge(92, 179, -1)
fold_tree.add_edge(92, 90, -1)
fold_tree.check_fold_tree()
if not fold_tree.check_fold_tree(): exit(0)
pose.fold_tree(fold_tree)
print(pose.fold_tree())

pose.set_phi(91, t)
pose.set_phi(  90, phi_to_recover)
pose.set_psi(  90, psi_to_recover)
pose.set_omega(90, omega_to_recover)
# for i in range(1, len(pose.sequence())):
#     pose.set_phi(i, pose2.phi(i))
#     pose.set_psi(i, pose2.psi(i))
#     pose.set_omega(i, pose2.omega(i))

# print("T:", t)
# print("PHI 89: %f" % (pose.phi(89)))
# print("PSI 89: %f" % (pose.psi(89)))
# print("OMEGA 89: %f" % (pose.omega(89)))

# pose.set_phi(91, t)
# print("PHI 91 %f" % (pose.phi(91)))
pose.dump_pdb("teste2.pdb")
exit(1)
# pose.residue(90).type().remove_variant_type(LOWER_TERMINUS_VARIANT)
# pose.residue(91).type().properties().set_property(LOWER_TERMINUS, False)
# print(pose.residue(90).show())
# print(pose.residue(91).show())
# pose.dump_pdb("teste.pdb")
# exit(1)
residue_provider = pose_from_sequence("A" + "A" * 10 + "A")

# Get last residue in chain A
chain_A = ChainSelector("A")
last_chain_A_res_index = get_residues_from_selector(chain_A, pose)[-1].seqpos()
for res_index in range(2, 11):
    res = residue_provider.residue(res_index)
    pose.append_polymer_residue_after_seqpos(res, last_chain_A_res_index, True)

    pose.set_phi(  last_chain_A_res_index, 180.0)
    pose.set_psi(  last_chain_A_res_index, 180.0)
    pose.set_omega(last_chain_A_res_index, 180.0)

    last_chain_A_res_index += 1

pose = reload_pose_from_file(pose)

# chain_C = ChainSelector("C")
# my_loop = Loop(80, 100, 90)
# add_single_cutpoint_variant(pose, my_loop)

# connect_A = get_residues_from_selector(chain_A, pose)[-1].seqpos()
# pose.residue(connect_A).type().remove_variant_type(UPPER_TERMINUS_VARIANT)

# connect_B = get_residues_from_selector(chain_C, pose)[ 0].seqpos()
# pose.residue(connect_B).type().remove_variant_type(LOWER_TERMINUS_VARIANT)

# pose = reload_pose_from_file(pose)

# declare_cutpoint_chemical_bond(pose, connect_B, connect_A)
# correctly_add_cutpoint_variants(pose, connect_A, connect_B)
# pose.dump_pdb("test.pdb")



# options = StructFileRepOptions()
# options.set_max_bond_length(9999.9)
# options.set_write_all_connect_info(True)
# dump_pdb(pose, "test.pdb", options)