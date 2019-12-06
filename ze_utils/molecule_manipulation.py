# TODO: Documentation

# LOAD from PDB, XYZ, GRO
# ATOM SELECTION/MASKING FUNCTIONS (...)
# SET CHAIN FROM RANGE
# REMOVE RESIDUES IN RANGE
# RENUMBER RESIDUES FROM START = INT 
# SORT RESIDUES BY CHAIN
# DEFINE CHAINS FROM CONNECTIONS
# AS MATRIX
# APPLY COORDINATES
# RMSD
# ALIGN
# PRINT as PDB, XYZ, GRO

import numpy as np
import sys

# d321 dictionary converts 3 letter aminoacid representation to 1 letter
d321 = {"CYS": "C", "ASP": "D", "SER": "S", "GLN": "Q", "LYS": "K",
        "ILE": "I", "PRO": "P", "THR": "T", "PHE": "F", "ASN": "N",
        "GLY": "G", "HIS": "H", "LEU": "L", "ARG": "R", "TRP": "W",
        "ALA": "A", "VAL": "V", "GLU": "E", "TYR": "Y", "MET": "M"}


class Molecule:
    """
    Holds all information and functions regarding a molecule and
    comparisons between molecules.
    """
    def __init__(self, filename):
        self.name = filename[:-4]
        self.atoms = []
        self.conects = []
        self.sequence = ""

        self.load(filename)
        

    def __str__(self):  
        self.as_pdb(sys.stdout)
        return ""
        
    # --- LOAD ---
    def load(self, filename):
        """
        Infer the reading function depending on the input
        filename extension.
        """

        if filename[-3:] == "pdb":
            self.read_pdb(filename)
        elif filename[-3:] == "gro":
            self.read_gro(filename)
        elif filename[-3:] == "xyz":
            self.read_xyz(filename)
        else:
            print("File format unkown")
            exit(0)
            

    def read_pdb(self, filename):
        """
        Fill self.atoms and self.conects information by parsing
        a PDB format file.
        """

        current_residue = "X"
        with open(filename, "r") as pdb:
            for line in pdb:
                
                if line[0:4] == "ATOM":
                    ls = line.split()
                    elem = ls[2]
                    if elem in ['HD1', 'HD2', 'HE1', 'HE2']:
                        elem = 'HD1'
                    if elem[0].isdigit():
                        elem = elem[1:]
                    # TODO: This selection should be based on regular
                    # expressions or character position on the string
                    self.atoms.append(Atom(int(ls[1]), elem, ls[3],
                        int(ls[5]), ls[4], float(ls[6]), float(ls[7]),
                        float(ls[8]), float(ls[9]), float(ls[10])))
                    if int(ls[5]) != current_residue:
                        current_residue = int(ls[5])
                        self.sequence += d321[ls[3]]
                    
                if line[0:6] == "CONECT":
                    ls = line.split()
                    self.conects.append(Conect(int(ls[1]),
                        [int(x) for x in ls[2:]]))


    def read_gro(self, filename):
        """
        Fill self.atoms and self.conects information by parsing
        a GRO format file.
        Cannot infer atom connections.
        """

        current_residue = "X"
        with open(filename, "r") as gro:
            for line in gro:
                ls = line.split()
                
                if len(ls) == 6:
                    self.atoms.append(Atom(int(ls[2]), ls[1],
                        ls[0][-3:], int(ls[0][0:-3]), "A", float(ls[3])*10.0,
                        float(ls[4])*10.0, float(ls[5])*10.0, -1.0, -1.0))
                    if ls[0][-3:] != current_residue:
                        current_residue = ls[0][-3:]
                        self.sequence += d321[current_residue]
                        

    def read_xyz(self, filename, unit="ang"):
        """
        Fill self.atoms and self.conects information by parsing
        a XYZ format file.
        Cannot infer sequence.
        Cannot infer atom connections.
        """

        index = 1
        with open(filename, "r") as xyz:
            for line in xyz:
                ls = line.split()
                
                if len(ls) > 1:
                    self.atoms.append(Atom(index, ls[0], "NaN", -1, "A",
                    float(ls[1]), float(ls[2]), float(ls[3]), 0.0, 0.0))
                    index += 1
    

    # --- MASKING ---
    def get_mask(self, expression, parameters):
        """
        Return a Boolean list with an entry for each atom in the Molecule. Set
        the entries whose index correspond to the list of atom indexes returned
        from the given 'expression' (acting on the given 'parameters') to True.
        This function should not be used directly.
        Instead, check "Atom selection/masking functions" bellow.

        See also: get_atoms_based_on_parameter
        """

        # 1) Find all atoms to be masked as True
        atoms, atom_indexes = getattr(self, expression)(parameters)

        # 2) Create the default mask (all False)
        mask = [False] * len(self.atoms)

        # 3) Set the masked atoms to True
        for masked_atom in atom_indexes:
            mask[masked_atom] = True

        # 4) Return mask
        return mask

    
    def get_atoms_based_on_parameter(self, parameter, query):
        """
        Iterate the atoms list and extract only the atoms (and corresponding
        indexes) whose 'parameter' exists in the 'query' list of parmeters.
        This function should not be used directly.
        Instead, check "Atom selection/masking functions" bellow.

        See also: get_mask
        """

        atoms, indexes = [], []
        for entry in query:
            for index, atom in enumerate(self.atoms):
                if getattr(atom, parameter) == entry:
                    atoms.append(atom)
                    indexes.append(index)

        return atoms, indexes


    # ... Atom selection/masking functions:
    # . Based on residue index
    def get_mask_for_residue_indexes(self, res_indexes):
        return self.get_mask("get_atoms_from_residue_indexes", res_indexes)


    def get_atoms_from_residue_indexes(self, res_indexes):
        return self.get_atoms_based_on_parameter("res_index", res_indexes)


    # . Based on atom element
    def get_mask_for_elements(self, elements):
        return self.get_mask("get_atoms_from_element", elements)


    def get_atoms_from_element(self, elements):
        return self.get_atoms_based_on_parameter("elem", elements)


    # --- GENERAL MANIPULATION ---
    def set_chain_from_range(self, chain_name, range_start, range_end):
        """
        Iterate over all atoms and change the atom.chain_name to match the
        given 'chain_name' parameter, if the atom.res_index is between the given
        'range_start' and 'range_end' (both aprameter should be integers).
        """

        assert type(range_start) == int, \
            "Range start should be an integer."
        assert type(range_end) == int, \
            "Range end should be an integer."
        assert type(chain_name) == str and len(chain_name) == 1, \
            "Chain name should be a 1 character string. (Ex: 'A')"

        for atom in self.atoms:
            if atom.res_index >= range_start and atom.res_index <= range_end:
                atom.chain_name = chain_name


    def remove_residues_in_range(self, range_start, range_end,
        renumber_atoms = True):
        """
        Iterate over all atoms in the molecule. Remove atoms whose res_index is
        on the range defined between 'range_start' and 'range_end' (both 
        parameters should be integers). Additionally, remove all conect records
        involving said atoms. If 'renumber_atoms' is set to True (by default),
        renumber all atoms to compensate for the missing residues (recommended).
        """

        assert type(range_start) == int, \
            "Range start should be an integer."
        assert type(range_end) == int, \
            "Range end should be an integer."

        # Atoms cannot be deleted during iteration, it would cause problems
        # where the length of the iterated list changes during the iteration.
        # Identified atoms for removal are marked in a list of later removal,
        # as well as the corresponding conect records.
        to_remove = []

        for (index, atom) in enumerate(self.atoms):
            if atom.res_index >= range_start and atom.res_index <= range_end:
                to_remove.append(index)

                # If conect records exist, remove all instances of this atom in
                # the conect records of neighbouring atoms
                if len(self.conects) == 0:
                    for bond in self.conects[index].bonded:
                        this_atom_index = self.conects[index].index
                        self.conects[bond - 1].bonded.remove(this_atom_index)

        # Removal of marked atoms needs to be made from the end to the beginning
        # to prevent change in the index of downstream atoms.
        to_remove.reverse()
        for index in to_remove:
            del self.atoms[index]

        # If conect records exist, delete all instances of the masked atoms
        if len(self.conects) > 0:
            for index in to_remove:
                del self.conects[index]

        # Since atoms have been removed, it's recommended to renumber all atoms
        if renumber_atoms:
            self.renumber_atoms()


    def renumber_residues(self, start = 1):
        """
        Renumber residues in molecule based on the current order, starting the
        count on integer 'start'.
        """

        assert type(start) == int, "'start; value for must be an int."

        residue_index = start - 1
        current_residue = None
        for atom in self.atoms:
            if atom.res_index != current_residue:
                current_residue = atom.res_index
                residue_index += 1
            atom.res_index = residue_index


    def renumber_atoms(self, start = 1):
        """
        Renumber atoms in molecule based on the current order, starting the
        count on integer 'start'.
        """

        assert type(start) == int, "'start; value for must be an int."

        conv = {}
        for (index, atom) in enumerate(self.atoms, start = start):
            conv[atom.index] = index
            atom.index       = index

        for conect in self.conects:
            conect.index     = conv[conect.index]
            conect.bonded    = [conv[bond] for bond in conect.bonded]


    def sort_residues_by_chain(self, renumber_residues = True):
        """
        Iterate over the different chain on the molecule in alphabetical order.
        Sort residues based on chain ID name (example, all residues from chain A
        appear first than all residues of chain B). Optionally (set to True by
        default), renumber all residues so each atom res_index matches their
        position on the sequence.
        """

        assert type(renumber_residues) == bool, \
            "'renumber_residues' parameter should be a boolean (True/False)"

        # 1. Get list of unique chains
        chains = []
        for atom in self.atoms:
            if atom.chain_name not in chains: chains.append(atom.chain_name)

        # 2. Sort chain lists
        chains = sorted(chains)

        # 3. Sort atoms based on the sorted list of chains
        atoms = []
        print(chains)
        for chain in chains:
            print("Chain:", chain)
            for atom in self.atoms:
                if atom.chain_name == chain: atoms.append(atom)
        self.atoms = atoms

        # 4. (Optional) Renumber residues
        self.renumber_residues()


    def define_chains_from_connections(self, chain_id = "ABCDEFGHIJKLMOPQRSTUV",
        overwrite = True):
        """
        Follow the connection maps (defined in the molecule.conects and
        therefore requires the input file to have CONECT records) and identify
        breaks in the connection tree. Each isolated tree of connections will
        correspond to a new chain. Return the number of chains identified.
        If 'overwrite' is set to True (yes, by default), actually change each
        atom chain_name parameter to match the identified chain. The order of
        naming each chain is set by 'chain_id' parameter (alphabetical order,
        by default).
        """

        assert type(chain_id) == str, \
            "Chain ID should be a string of chain names"
        assert type(overwrite) == bool, \
            "'overwrite' parameter should be a boolean (True/False)"

        n_chains = 0
        atom_list = [atom.index for atom in self.atoms]
        while len(atom_list) > 0:
            stack = [atom_list[0]]
            chain = []
            while len(stack) > 0:
                for connection in self.conects[stack[0] - 1].bonded:
                    if connection in atom_list and connection not in stack:
                        stack.append(connection)
                if overwrite:
                    atom = self.atoms[stack[0] - 1]
                    self.atoms[stack[0] - 1].chain_name = chain_id[n_chains]
                atom_list.remove(stack[0])
                stack.pop(0)
            n_chains += 1

        return n_chains
    

    # --- BLOSUM62 ---
    def blosum62(self, reference):
        """
        Return the BLOSUM62 score between this molecule and a reference
        Molecule object.
        """

        from ze_utils.common import read_matrix_from_txt_file, blosum62

        assert len(self.sequence) > 0, \
            "The molecule doesn't seem to have a sequence: try loading from PDB"
        assert len(self.sequence) == len(reference.sequence), \
            "This molecule and the reference's sequence must be of same length"

        # Read the BLOSUM62 matrix from the default location
        matrix, entries = read_matrix_from_txt_file()

        # Calculate the BLOSUM62 score between the this molecule and reference
        # sequences
        return blosum62(matrix, entries, self.sequence, reference.sequence)[0]


    def conserved(self, reference):
        """
        Return the percentage of conserved sequence between this molecule and a
        reference Molecule object.
        """

        assert len(self.sequence) > 0, \
            "The molecule doesn't seem to have a sequence: try loading from PDB"
        assert len(self.sequence) == len(reference.sequence), \
            "This molecule and the reference's sequence must be of same length"

        conserved = 0
        for index, query_aminoacid in enumerate(self.sequence):
            if query_aminoacid == reference.sequence[index]:
                conserved += 1

        return (conserved / len(self.sequence)) * 100


    # --- ALIGN ---
    def apply_coordinates(self, new_coordinates):
        """
        Apply a matrix of coordinates to this molecule.
        """

        for atom, coord in zip(self.atoms, new_coordinates):
            atom.x = coord[0]
            atom.y = coord[1]
            atom.z = coord[2]
            

    def as_matrix(self, mask = None):
        """
        Return a matrix of coordinates based on this molecule. If a mask is
        provided, only the atoms whose index on the mask is True are returned.
        """

        if not mask == None:
            true_mask = len([x for x in mask if x])
            assert true_mask > 0, "Mask must have more than 0 atoms selected"

        coords = []
        for index, atom in enumerate(self.atoms):
            if mask == None or mask[index] == True:
                coords.append([atom.x, atom.y, atom.z])
        return np.array(coords)
        

    def rmsd(self, reference, movable_mask = None, reference_mask = None):
        """
        Calculate the RMSD between this molecule and a reference object. Also
        returns the number of atoms considered for this calculation. If movable
        and references masks are provided, only consider the atoms whose index
        in the corresponding mask is set to True.
        """

        movable_coords   = self.as_matrix(movable_mask)
        reference_coords = reference.as_matrix(reference_mask)

        assert len(movable_coords) == len(reference_coords), \
            "Movable (%d) and Reference (%d) number of atoms do not match" % \
            (len(movable_coords), len(reference_coords))

        n_atoms          = len(reference_coords)
        distance         = reference_coords - movable_coords

        return n_atoms, np.sqrt(np.sum(distance * distance) / n_atoms)
            

    def align(self, reference, movable_mask = None, reference_mask = None,
        verbose = True):
        """
        Align this molecule with a reference object. If movable and reference
        masks are provided (Default: None), only the atoms whose index in the
        corresponding mask is set to True will be considered for alignment. Both
        sets of masked atoms must be on the same size. Verbose flag determines
        if the RMSD value is computed and printed (Default: True).
        """
            
        # 1) Get coordinates as a matrix
        movable_coords    = self.as_matrix(movable_mask)
        reference_coords  = reference.as_matrix(reference_mask)

        assert len(movable_coords) == len(reference_coords), \
            "Movable (%d) and Reference (%d) number of atoms do not match" % \
            (len(movable_coords), len(reference_coords))
        
        # 2) Center on centroid
        c_movable        = np.mean(movable_coords, axis=0)
        c_reference      = np.mean(reference_coords, axis=0)
        movable_coords   = movable_coords   - c_movable
        reference_coords = reference_coords - c_reference
        
        # 3) Obtain correlation matrix
        cm = np.dot(np.transpose(movable_coords), reference_coords)
        u, d, vt = np.linalg.svd(cm)
        
        # 4) Obtain rotation
        rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
        
        # 5) Apply rotation + translation
        transformed_coords  = np.dot(self.as_matrix() - c_movable, rot)
        transformed_coords += c_reference
        
        # 6) Apply transformed coordinates
        self.apply_coordinates(transformed_coords)
        
        # 7) Calculate RMSD (Optional)
        count, rms = self.rmsd(reference, movable_mask, reference_mask)
        if verbose:
            print("RMSD: %6.3f angstrom (%3d atoms)" % (rms, count))
        
        return rms
            
            
    # --- EXPORT ---
    def export(self, filename, title = "Molecule"):
        """
        Infer the export function depending on the output filename extension and
        write the molecule contents in the appropriate format to a file.
        """

        if   filename[-3:] == "pdb":
            self.as_pdb_to_file(filename, title)
        elif filename[-3:] == "gro":
            self.as_gro(filename, title)
        elif filename[-3:] == "xyz":
            self.as_xyz(filename)
        else:
            print("File format for output unkown")
            exit(0)
                        
                         
    def as_pdb_to_file(self, filename, title = "Molecule"):
        """
        Export this molecule information in PDB format to a file.
        """

        with open(filename, "w") as pdb:
            self.as_pdb(pdb, title)

    def as_pdb(self, des, title = "Molecule"):
        """
        Export this molecule information in PDB format.
        """

        from io import TextIOWrapper

        assert type(des) == TextIOWrapper, \
            "PDB destination 'des' must be of type TextIOWrapper."

        PDB = "ATOM %6d  %-3s %3s %1s %3d" + \
                "%11.3f %7.3f %7.3f %5.2f %5.2f %10s%1s\n"
        des.write("TITLE %s\n" % (title))
        des.write("MODEL 1\n")
        for atom in self.atoms:
            des.write(atom.as_pdb())
        des.write("TER\n")
        for conect in self.conects:
            des.write("CONECT %4d" % (conect.index))
            for bond in conect.bonded:
                des.write(" %4d" % (bond))
            des.write("\n")


    def as_gro(self, filename, title = "Molecule"):
        """
        Export this molecule information in GRO format.
        """

        with open(filename, "w") as gro:
            gro.write("%s\n" % title)
            gro.write("%5d\n" % len(self.atoms))
            for atom in self.atoms:
                gro.write("%5d%3s %6s %4d %7.3f %7.3f %7.3f\n" %
                    (atom.res_index - 1, atom.res_name, atom.elem,
                    atom.index, atom.x/10.0, atom.y/10.0, atom.z/10.0))
            gro.write("%10.5f%10.5f%10.5f\n\n" % (0.0, 0.0, 0.0))
            

    def as_xyz(self, filename):
        """
        Export this molecule information in XYZ format.
        """

        with open(filename, "w") as xyz:
            xyz.write("%d\n" % len(self.atoms))
            for atom in self.atoms:
                xyz.write("%s   %f   %f   %f\n" %
                    (atom.elem, atom.x, atom.y, atom.z))
                    

class Atom:
    """
    Holds all information regarding an Atom.
    """

    def __init__(self, index=-1, elem="NaN", res_name="NaN", res_index=-1.0,
        chain_name="A", x=0.0, y=0.0, z=0.0, mass=0.0, charge=0.0):
            
        self.index      = index
        self.elem       = elem
        self.res_name   = res_name
        self.res_index  = res_index
        self.chain_name = chain_name
        self.x          = x
        self.y          = y
        self.z          = z
        self.mass       = mass
        self.charge     = charge


    def as_pdb(self):
        """
        Print this single atom in PDB format.
        """

        PDB = "ATOM %6d  %-3s %3s %1s %3d" + \
            "%11.3f %7.3f %7.3f %5.2f %5.2f %10s%1s\n"

        return PDB % \
            (self.index, self.elem, self.res_name,
            self.chain_name, self.res_index, self.x,
            self.y, self.z, self.mass, self.charge,
            " ", self.elem[0])

    def xyz(self):
        """
        Return the Atom coordinates in [X, Y, Z] format.
        """

        return [self.x, self.y, self.z]

    def __str__(self):
        return self.as_pdb()
        
    
    
class Conect:
    """
    Holds all information regarding a Connection.
    """

    def __init__(self, index=-1, bonded=[]):
        self.index  = index
        self.bonded = bonded