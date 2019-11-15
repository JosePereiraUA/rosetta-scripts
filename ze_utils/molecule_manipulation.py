# TODO: Documentation

# LOAD from PDB, XYZ, GRO
# SET CHAIN FROM RANGE
# REMOVE RESIDUES IN RANGE
# RENUMBER RESIDUES FROM START = INT 
# SORT RESIDUES BY CHAIN
# DEFINE CHAINS FROM CONNECTIONS
# COUNT ATOMS OF ELEMENT
# COMPARE ORDER
# AS MATRIX
# APPLY COORDINATES
# RMSD
# ALIGN
# PRINT as PDB, XYZ, GRO

import numpy as np

class Molecule:
    """
    Holds all information and functions regarding a molecule and
    comparisons between molecules.
    """
    def __init__(self, filename):
        self.name = filename[:-4]
        self.atoms = []
        self.conects = []
        self.load(filename)
        
        
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
                    
                if line[0:6] == "CONECT":
                    ls = line.split()
                    self.conects.append(Conect(int(ls[1]),
                        [int(x) for x in ls[2:]]))

    def read_gro(self, filename):
        """
        Fill self.atoms and self.conects information by parsing
        a GRO format file.
        """

        with open(filename, "r") as gro:
            for line in gro:
                ls = line.split()
                
                if len(ls) == 6:
                    self.atoms.append(Atom(int(ls[2]), ls[1],
                        ls[0][-3:], int(ls[0][0:-3]), "A", float(ls[3])*10.0,
                        float(ls[4])*10.0, float(ls[5])*10.0, -1.0, -1.0))
                        
    def read_xyz(self, filename, unit="ang"):
        """
        Fill self.atoms and self.conects information by parsing
        a XYZ format file.
        """

        index = 1
        with open(filename, "r") as xyz:
            for line in xyz:
                ls = line.split()
                
                if len(ls) > 1:
                    self.atoms.append(Atom(index, ls[0], "NaN", -1, "A",
                    float(ls[1]), float(ls[2]), float(ls[3]), 0.0, 0.0))
                    index += 1
    
    # --- GENERAL MANIPULATION ---
    def set_chain_from_range(self, chain_name, range_start, range_end):
        """
        """

        for atom in self.atoms:
            if atom.res_index >= range_start and atom.res_index <= range_end:
                atom.chain_name = chain_name

    def remove_residues_in_range(self, range_start, range_end):
        """
        """

        to_remove = []
        for (index, atom) in enumerate(self.atoms):
            if atom.res_index >= range_start and atom.res_index <= range_end:
                to_remove.append(index)
        to_remove.reverse()
        for index in to_remove:
            del self.atoms[index]

    def renumber_residues(self, start = 1):
        """
        """

        residue_index = start - 1
        current_residue = None
        for atom in self.atoms:
            if atom.res_index != current_residue:
                current_residue = atom.res_index
                residue_index += 1
            atom.res_index = residue_index

    def sort_residues_by_chain(self, renumber_residues = True):
        """
        """

        # 1. Get list of unique chains
        chains = []
        for atom in self.atoms:
            if atom.chain_name not in chains: chains.append(atom.chain_name)
        # 2. Sort chain lists
        chains = sorted(chains)
        # 3. Sort atoms based on the sorted list of chains
        atoms = []
        for chain in chains:
            for atom in self.atoms:
                if atom.chain_name == chain: atoms.append(atom)
        self.atoms = atoms
        # 4. (Optional) Renumber residues
        self.renumber_residues()

    def define_chains_from_connections(self, chain_id = "ABCDEFGHIJKLMOPQRSTUVWXYZ",
        overwrite = True):
        """
        """

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
                    self.atoms[stack[0] - 1].chain_name = chain_id[n_chains]
                atom_list.remove(stack[0])
                stack.pop(0)
            n_chains += 1
        return n_chains
    
    # --- ALIGN ---
    def count(self, elem = None):
        """
        Count the number of atoms in this molecule of the provided
        element `elem`. By default, counts total number of atoms.
        """

        count = 0
        for atom in self.atoms:
            if elem == None or atom.elem in elem:
                count += 1
        return count
        
    def compare_order(self, reference, elem = None):
        """
        Compare the order of the atoms between this molecule and a
        reference object. If an array of elements is provided, the
        comparison is only performed based on the subset of atoms
        of those elements.
        """

        for a1, a2 in zip(self.atoms, reference.atoms):
            if elem == None or (a1.elem in elem and a2.elem in elem):
                if not a1.elem == a2.elem:
                    return False
        return True
        
    def verify_input(self, reference, elem = None):
        """
        Verify the input for proper alignment:
        1) Number of atoms should match;
        2) Atoms should be in the same order;
        3) All elements provided in the filter must exist in both
        structures
        """

        assert self.count(elem) == reference.count(elem),\
        "Movable and Refence structures don't have the same number of atoms"

        assert self.compare_order(reference, elem),\
        "Movable and Reference structures do not have the same element order"

        if not elem == None:
            for i in elem:
                assert self.count([i]) > 0,\
                "Atom element %s not found in Movable structure" % (i)
                assert reference.count([i]) > 0,\
                "Atom element %s not found in Reference structure" % (i)
                
    def apply_coordinates(self, new_coordinates):
        """
        Apply a matrix of coordinates to this molecule.
        """

        for atom, coord in zip(self.atoms, new_coordinates):
            atom.x = coord[0]
            atom.y = coord[1]
            atom.z = coord[2]
            
    def as_matrix(self, elem = None):
        """
        Return a matrix of coordinates based on this molecule.
        """

        n = len(self.atoms)
        coords = []
        for index, atom in enumerate(self.atoms):
            if elem == None or atom.elem in elem:
                coords.append([atom.x, atom.y, atom.z])
        return np.array(coords)
        
    def rmsd(self, reference, elem = None):
        """
        Calculate the RMSD between this molecule and a reference
        object. Also returns the number of atoms considered for
        this calculation.
        """

        m = self.as_matrix(elem)
        r = reference.as_matrix(elem)
        n = len(r)
        d = r - m
        return n, np.sqrt(np.sum(d * d) / n)
            
    def align(self, reference, elem = None, verbose = True):
        """
        Align this molecule with a reference object. If an array
        of elements is provided, the alignment is only performed based
        on the subset of atoms of those elements. Verbose flag
        determines if the RMSD value is computed and printed (Default:
        True).
        """

        # 0) Input verifications
        self.verify_input(reference, elem)
            
        # 1) Get coordinates as a matrix
        movable_coords    = self.as_matrix(elem)
        reference_coords  = reference.as_matrix(elem)
        
        # 2) Center on centroid
        n                = len(self.atoms)
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
        count, rms = self.rmsd(reference, elem)
        if verbose:
            print("RMSD: %6.3f angstrom (%3d atoms)" % (rms, count))
        return rms
            
            
    # --- EXPORT ---
    def print_structure(self, filename, title="Aligned molecule"):
        """
        Infer the export function depending on the output filename
        extension.
        """

        if   filename[-3:] == "pdb":
            self.print_as_pdb(filename, title)
        elif filename[-3:] == "gro":
            self.print_as_gro(filename, title)
        elif filename[-3:] == "xyz":
            self.print_as_xyz(filename)
        else:
            print("File format for output unkown")
            exit(0)
                        
                         
    def print_as_pdb(self, filename, title="Aligned molecule"):
        """
        Export this molecule information in PDB format.
        """

        with open(filename, "w") as pdb:
            PDB = "ATOM %6d  %-3s %3s %1s %3d" + \
                  "%11.3f %7.3f %7.3f %5.2f %5.2f %10s%1s\n"
            pdb.write("TITLE %s\n" % (title))
            pdb.write("MODEL 1\n")
            for atom in self.atoms:
                pdb.write(atom.as_pdb())
            pdb.write("TER\n")
            for conect in self.conects:
                pdb.write("CONECT %4d" % (conect.index))
                for bond in conect.bonded:
                    pdb.write(" %4d" % (bond))
                pdb.write("\n")
                
    def print_as_gro(self, filename, title="Aligned molecule"):
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
            
    def print_as_xyz(self, filename):
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
        """

        PDB = "ATOM %6d  %-3s %3s %1s %3d" + \
            "%11.3f %7.3f %7.3f %5.2f %5.2f %10s%1s\n"

        return PDB % \
            (self.index, self.elem, self.res_name,
            self.chain_name, self.res_index, self.x,
            self.y, self.z, self.mass, self.charge,
            " ", self.elem[0])

    def __str__(self):
        return self.as_pdb()
        
    
    
class Conect:
    """
    Holds all information regarding a Connection.
    """

    def __init__(self, index=-1, bonded=[]):
        self.index  = index
        self.bonded = bonded