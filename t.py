from ze_utils.molecule_manipulation import *
m = Molecule("3ch8.pdb")
m.define_chains_from_connections("ACB")
m.sort_residues_by_chain()
m.remove_residues_in_range(90, 108)
m.print_structure("3ch8_chains.pdb")