import numpy as np
from design import *
from pyrosetta import *
from ze_utils.molecule_manipulation import *
from ze_utils.pyrosetta_tools import *
from ze_utils.pyrosetta_classes import *
from pyrosetta.rosetta.numeric import *
from pyrosetta.rosetta.protocols.rigid import *
from pyrosetta.rosetta.protocols.docking import *
from pyrosetta.rosetta.core.conformation import Atom
from pyrosetta.rosetta.core.select.residue_selector import *
from pyrosetta.rosetta.core.simple_metrics.metrics import \
	InteractionEnergyMetric

init()

#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\

#                          Dock-Design Script:
# ______________________________________________________________________________
#  Performs docking + design algorithm. The objetive of this script is to sample
# the position of a movable chain on the input structure, and for each provable
# new position, sample the sequence in the interface in order to design, at the
# same time, both the position and the sequence of part of a protein. The script
# assumes an ABC model for the protein. This model essentially assumes that the
# protein is split in three different chains on the input PDB file. Those parts
# should be defined as follows:
#  - Chain A: (Fixed) Fixed part of the protein (PBZ, in this example);
#  - Chain B: (Fixed) Target peptide ligand around which the interface will be 
# defined;
#  - Chain C: (Movable) Movable part of the protein, which will be subject to
# the pseudo docking process (random movement + design);
#  This model can be set-up using the 'ze_utils.molecule_manipulation' module,
# as follows:
#
# pdb = Molecule("clamshell_relax.pdb")       # Load structure
# pdb.set_chain_from_range('A', 376, 473)     # Set Chain A
# pdb.set_chain_from_range('B',   1,   7)     # Set Chain B
# pdb.set_chain_from_range('C', 483, 572)     # Set Chain C
# pdb.remove_residues_in_range(474, 482)         # Remove atoms in the loop region
# pdb.sort_residues_by_chain()                # (Optional) Renumber all residues
# pdb.print_structure("clamshell_3p_rlx.pdb") # Save edited structure
#
# For more information and similar scripts, please read:
# https://graylab.jhu.edu/pyrosetta/downloads/scripts/demo/D100_Docking.py


# !
# os.system("rm -f pos* sop*")

# pose = pose_from_pdb("starting_structures/clamshell_3p_rlx.pdb")
# # chainA =  ChainSelector("A")
# chainB =  ChainSelector("B")
# chainC =  ChainSelector("C")

# centroidA = get_centroid_coordinates_from_selector(ChainSelector("A"), pose)
# centroidB = get_centroid_coordinates_from_selector(ChainSelector("B"), pose)
# centroidC = get_centroid_coordinates_from_selector(ChainSelector("C"), pose)

# bc_distance = np.linalg.norm(centroidC - centroidB)
# ab_normalized = (centroidB - centroidA) / np.linalg.norm(centroidB - centroidA)
# grid_center = centroidB + bc_distance * ab_normalized


# with open("centroids.pdb", "w") as file_out:
#     for p in [centroidA, centroidB, centroidC, grid_center]:
#         a = ze_utils.common.Atom(chain_name="D", x = p[0], y = p[1], z = p[2])
#         file_out.write(a.as_pdb())

# dg = DockingGrid(grid_center, (20, 40), (10, 20), (30, 40), 5, 3, 5)
# dg = DockingGrid(grid_center, (-2, 40), (2, 10), (12, 38), 2, 2, 2) # <-
# dg.generate_grid_positions()
# dg.orient(np.array([0, 0, 1]), ab_normalized) # Orient Z with AB # <-

# lig_res = get_residues_from_subset(chainB.apply(pose))
# D = np.array(pose.residue(list(lig_res)[0]).atom(2).xyz())
# E = np.array((pose.residue(list(lig_res)[len(lig_res) - 1]).atom(2).xyz()))
# de_normalized = (E - D) / np.linalg.norm(D - E)

# dg.orient(np.array([1, 0, 0]), de_normalized) # Orient X with DE

# score_function = get_fa_scorefxn()
# job_man = PyJobDistributor("pos", len(dg.points), score_function)
# p = Pose()
# inner_jobs = {}
# while not job_man.job_complete:
#     inner_jobs[job_man.current_id] = PyJobDistributor("sop", 2, score_function)
#     print(inner_jobs)
#     p.assign(pose)
#     dg.translate_selector_to_point(job_man.current_id, ChainSelector("C"), pose)
#     while not inner_jobs[job_man.current_id].job_complete:
#         inner_jobs[job_man.current_id].output_decoy(p)
#     job_man.output_decoy(p)

# dg.as_pdb("dg.pdb")
# with open("teste2.pdb", "w") as file_out:
#     for point in dg.points:
#         file_out.write(ze_utils.common.Atom(chain_name="D", x = point[0], y = point[1], z = point[2]).as_pdb())

#-----------------------

# Set up the fold tree so that only the chain C is movable
fold_tree = FoldTree()                      # Start with an empty fold tree
fold_tree.add_edge(1, 105, -1)              # Fist edge is both chain A and B
fold_tree.add_edge(105, 106, 1)             # Jump 1 is between chain B and C
fold_tree.add_edge(106, 195, -1)            # Second edge is chain C
if not fold_tree.check_fold_tree(): exit(0) # Check if fold tree is viable
pose.fold_tree(fold_tree)                   # Apply the new fold tree to pose

# switch_to_centroid = SwitchResidueTypeSetMover("centroid")
# centroid_pose      = Pose(pose)
# switch_to_centroid.apply(centroid_pose)

# ----------------------


# pre_filter = PreFilter()

# cutoff = 9.0
# interface = NeighborhoodResidueSelector(chainB, cutoff,
#     include_focus_in_subset = False)
# designable = AndResidueSelector(interface, chainC)
# repackable = chainB

# MOVE
# p = Pose(pose)
# rotation    = 3.0
# translation = 0.5
# kT = 1.0
# perturbator = RigidBodyPerturbMover(1, rotation, translation)
# # slider      = DockingSlideIntoContact(1)
# score_function = get_fa_scorefxn()

# dock = SequenceMover()
# dock.add_mover(perturbator)
# # dock.add_mover(slider)
# design_mover = get_design_mover(p, score_function, designable, repackable)
# mc = MonteCarlo(pose, score_function, kT)
# design = TrialMover(design_mover, mc)

# approved_count = 0
# n_steps        = 5000

# pmm            = PyMOLMover()
# pmm.keep_history(True)
# pmm.apply(p)

# initial_score = score_function(p)
# iem = InteractionEnergyMetric(repackable, designable)



# pose.pdb_info().name("XXXX")
# teste = PASSO(50)
# teste.apply(pose)

# print("DONE!")
# print("Acceptance ratio %12.4f" % (approved_count / n_steps))
# exit(1) # !!!!!!!!!!!!!!!!

# from matplotlib import pyplot as plt
# # fig, ax = plt.subplots(nrows=1, ncols=4)
# fig, ax = plt.subplots(nrows=1, ncols=3)

# def scatter(ax, key):
#     a_y = [y[key][0] for y in run_data if y[key][1] == True]
#     r_y_r = [y[key][0] for y in run_data if (y[key][1] == False and y["step"][1] == False)]
#     r_y_nr = [y[key][0] for y in run_data if (y[key][1] == True and y["step"][1] == False)]
#     a_x = [x["step"][0] for x in run_data if x[key][1] == True]
#     r_x_r = [x["step"][0] for x in run_data if (x[key][1] == False and x["step"][1] == False)]
#     r_x_nr = [x["step"][0] for x in run_data if (x[key][1] == True and x["step"][1] == False)]
    
#     ax.scatter(a_x, a_y, c="tab:green")
#     ax.scatter(r_x_nr, r_y_nr, c="tab:orange")
#     ax.scatter(r_x_r, r_y_r, c="tab:red")
#     ax.title.set_text(key)

# scatter(ax[0], "contacts")
# scatter(ax[1], "distance")
# scatter(ax[2], "clashes")
# # ax[3].bar([x[0] for x in energy_data], [x[1] for x in energy_data])
# ax[0].axhline(pre_filter.contact_min_count)
# ax[1].axhline(pre_filter.anchors_cutoff)
# ax[2].axhline(pre_filter.clashes_max_count)
# plt.show()


# ---------------------------

# n_atoms = 2000
# for i in range(1, n_atoms):
#     x, y, z = 50 * get_rand_vector_in_sphere() # Get translation vector!
#     z2 = z + 10
#     a = np.array([0, 0, 1])
#     b = np.array([x, y, z])
#     c = np.array([x, y, z2])
#     axis = np.cross(a, b)
#     angle = get_angle_from_two_vectors(a, b)
#     rot_matrix = get_rotation_matrix_from_axis_angle(axis, angle)
#     x2, y2, z2 = rotate_coords_from_rotation_matrix(c, rot_matrix, b)
#     print("ATOM %6d  C   UNK A   1     %7.3f %7.3f %7.3f  1.00  0.00           C" % (i, x, y, z))
#     print("ATOM %6d  C   UNK A   1     %7.3f %7.3f %7.3f  1.00  0.00           C" % (i + n_atoms, x2, y2, z2))

# for i in range(1, n_atoms):
#     print("CONECT %4d %4d" % (i, i+n_atoms))