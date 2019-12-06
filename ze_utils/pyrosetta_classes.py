from pyrosetta import *
init()

class PreFilter:
    """
    Simple filter that can be applied to a pose. Returns true if all conditions
    present are valid. The conditions are the following:

     . Number of clashes (defined the number of residue pairs between 'sel_B'
     and 'sel_C' whose distance is inferior to the defined 'clash_cutoff'. The
     number of clashes must be inferior to 'clashes_max_count'.)

     . Number of contacts (defined the number of residue pairs between 'sel_B'
     and 'sel_C' whose distance is inferior to the defined 'contact_cutoff' but
     superior to the defined 'clash_cutff'. This number should, therefore, be
     higher than the 'clash_cutoff'. The number of contacts must be superior to
     the 'contact_min_count'.)

     . Distance of the anchors (defined as the euclidean distance between two
     atoms in the pose. This atoms correspond to the possible attachment points
     for a loop unityng the two regions defined in 'sel_A' and 'sel_B'. This
     distance must be inferior to 'anchors_cutoff'.

     . If set to True, both the N and C terminals of the ligand/sel_C (chain C
     in the ABC Model, by default) can have their blocking prevented. This means
     that a 'max_terminal_interaction' value can be set, which will be surpassed
     in the case of that terminal being blocked by the moving part/sel_B (chain
     B in the ABC model, by default), in which case the filter condition will
     return False. The considered interaction is the absolute value of the
     InteractionEnergyMetric (since it should consider both favourable and 
     disruptive interactions). If preventing block in one or more terminals, a
     score_function needs to be provided. If set to "auto", score function will
     be set using the get_fa_scrofxn() function.

    All cutoff distances in the PreFilter are in Angstrom.
    By default (if set to "auto"), the selections 'sel_A', 'sel_B' and 'sel_C'
    refer to the Chains A, B and C, respectively (in accordance to the ABC
    model - recommended).
    """

    def __init__(self, contact_cutoff = 9.0, contact_min_count = 10,
        anchors_cutoff = 25.0, clash_cutoff = 4.0, clashes_max_count = 2,
        sel_A = "auto", sel_B = "auto", sel_C = "auto",
        prevent_block_C_terminal = False, prevent_block_N_terminal = True,
        max_c_terminal_interaction = 0.0, max_n_terminal_interaction = 0.05,
        score_function = "auto"):

        from pyrosetta.rosetta.core.scoring import ScoreFunction
        from pyrosetta.rosetta.core.select.residue_selector import \
            ResidueSelector, ChainSelector

        assert isinstance(sel_A, ResidueSelector) or sel_A == "auto", \
            "Selection A must be a ResidueSelector or set to 'auto'"
        assert isinstance(sel_B, ResidueSelector) or sel_B == "auto", \
            "Selection B must be a ResidueSelector or set to 'auto'"
        assert isinstance(sel_C, ResidueSelector) or sel_C == "auto", \
            "Selection C must be a ResidueSelector or set to 'auto'"

        self.clashes_max_count          = clashes_max_count
        self.contact_min_count          = contact_min_count
        self.contact_cutoff             = contact_cutoff
        self.anchors_cutoff             = anchors_cutoff
        self.clash_cutoff               = clash_cutoff
        self.prevent_block_C_terminal   = prevent_block_C_terminal
        self.prevent_block_N_terminal   = prevent_block_N_terminal
        self.max_c_terminal_interaction = max_c_terminal_interaction
        self.max_n_terminal_interaction = max_n_terminal_interaction

        self.sel_A = ChainSelector("A") if sel_A == "auto" else sel_A
        self.sel_B = ChainSelector("B") if sel_B == "auto" else sel_B
        self.sel_C = ChainSelector("C") if sel_C == "auto" else sel_C

        assert type(score_function) == ScoreFunction or \
            score_function == "auto",\
            "Score function must be of type ScoreFunction or set to 'auto'"
        if score_function == "auto":
            self.score_function = get_fa_scorefxn()
        else:
            self.score_function = score_function

        
    def apply(self, pose):
        """
        Applies the filter to the given 'pose'.
        Returns a Boolean stating if the pose passes the filter (True) or not
        (False); and a data structure. The data structure is a dictionary with
        an entry for each parameter assessed. Each of this entries is a Tuple
        with two values: the first is the calculated parameter (ex: clash count)
        and the second is a Boolean stating wether this parameter in specific
        passed the filter (True) or not (False).
        """

        import numpy as np
        from ze_utils.pyrosetta_tools import \
            calc_d_table_between_selectors, get_residues_from_selector
        from pyrosetta.rosetta.core.select.residue_selector import \
            ResidueIndexSelector
        from pyrosetta.rosetta.core.simple_metrics.metrics import \
	        InteractionEnergyMetric

        to_return = True
        data = {"contacts"  : (-1,  False),
                "distance"  : (-1,  False),
                "clashes"   : (-1,  False)}

        if self.prevent_block_C_terminal:
            data["c_terminal"] = (0.0, False)
        if self.prevent_block_N_terminal:
            data["n_terminal"] = (0.0, False)

        
        # Create distance table for pose
        d_table = calc_d_table_between_selectors(self.sel_C, self.sel_B, pose)
        
        # Count the number of clashes
        clashes = len(d_table[d_table < self.clash_cutoff])
        if clashes > self.clashes_max_count:
            data["clashes"]  = (clashes, False)
            to_return        = False
            print("[PreFilter] Clash count (%d) > Max permitted clashes (%d)" \
                % (clashes, self.clashes_max_count))
        else:
            data["clashes"]  = (clashes, True)


        # Count number of contacts in pose
        contacts = len(d_table[d_table < self.contact_cutoff]) - clashes
        if contacts < self.contact_min_count:
            data["contacts"] = (contacts, False)
            to_return        = False
            print("[PreFilter] Contacts (%d) < Minimum contacts (%d)" \
                % (contacts, self.contact_min_count))
        else:
            data["contacts"] = (contacts, True)

        # Get the distance from the closing loop anchors
        upstream_residues    = get_residues_from_selector(self.sel_A, pose)
        upstream_res         = upstream_residues[-1]
        anchor_A             = upstream_res.atom("N")

        downstream_residues  = get_residues_from_selector(self.sel_B, pose)
        downstream_res       = downstream_residues[0]
        anchor_B             = downstream_res.atom("C")

        xyz_A                = np.array(anchor_A.xyz())
        xyz_B                = np.array(anchor_B.xyz())
        distance             = np.linalg.norm(xyz_B - xyz_A)

        if distance > self.anchors_cutoff:
            data["distance"] = (distance, False)
            to_return        = False
            print("[PreFilter] Anchors distance (%f) > Cutoff (%f)" \
                % (distance, self.anchors_cutoff))
        else:
            data["distance"] = (distance, True)

        # Verify if terminals are blocked. This function assumes the ABC Model
        # and that the sequence of the ligand (chain C) goes from N terminal
        # first to C terminal last.
        lig = get_residues_from_selector(self.sel_C, pose)

        # N Terminal
        if self.prevent_block_N_terminal:
            self.score_function(pose)
            n_terminal = ResidueIndexSelector(lig[0].seqpos())
            ie = InteractionEnergyMetric(n_terminal, self.sel_B)
            n_terminal_interaction = abs(ie.calculate(pose))
            if n_terminal_interaction > self.max_n_terminal_interaction:
                data["n_terminal"] = (n_terminal_interaction, False)
                to_return          = False
                print("[PreFilter] N Terminal Interaction: (%f) > Max (%f)" \
                    % (n_terminal_interaction, self.max_n_terminal_interaction))
            else:
                data["n_terminal"] = (n_terminal_interaction, True)

        # C Terminal
        if self.prevent_block_C_terminal:
            self.score_function(pose)
            c_terminal = ResidueIndexSelector(lig[-1].seqpos())
            ie = InteractionEnergyMetric(c_terminal, self.sel_B)
            c_terminal_interaction = abs(ie.calculate(pose))
            if c_terminal_interaction > self.max_c_terminal_interaction:
                data["c_terminal"] = (c_terminal_interaction, False)
                to_return          = False
                print("[PreFilter] C Terminal Interaction: (%f) > Max (%f)" \
                    % (c_terminal_interaction, self.max_c_terminal_interaction))
            else:
                data["c_terminal"] = (c_terminal_interaction, True)

        return (to_return, data)


class PASSO:
    """
                Position And Sequence Simultaneous Optimizer

    Performs simultaneous optimization of the sequence and the position
    placement of (part of) a protein. The simulation runs for 'n_steps'. Each
    step has a maximum of 3 stages:
     1 - Random placement of the movable block;
     2 - Evaluation of a PreFilter;
     3 - (If step 2 is OK) Design of the interface residues (repeated 'n_cycles'
     times - 2 by default).
    The movable region is defined in the 'dock_mover', responsible for changing
    the position of this movable region. By default (by setting this parameter
    to "auto"), it's a RigidBodyPerturbMover that randomly adds a rotation
    (maximum: 3 degrees) and a translation (maximum: 0.5 angstrom) to the pose
    edge immediatly after Jump 1 in the pose Fold Tree.
    The designable region is defined in the 'design_mover', responsible for
    changing the nature of the aminoacids in this designable region. In most 
    cases, it should be the interface between the movable region and the fixed
    region. By default (by setting this parameter to "auto"), it's a
    SequenceMover with 2 inner movers: a PackRotamersMover and a MinMover, in
    this order. For the PackRotamersMover, extra rotamers on chi1 and chi2 are
    enabled. During this step, residues on the 'designable' ResidueSelector will
    be subject to design efforts, while residues on the 'repackable'
    ResidueSelector will only change conformation to rotamers of the same
    aminoacid. During the minimization step only the sidechains are allowed to
    relax.
    Even if a custom 'design_mover' is provided, the 'designable' (interface in
    the movable region) and 'repackable' (target ligand) ResidueSelector must
    still be provided, as they are used in calculating the interaction energy
    between this two sets. By default (by setting this parameters to "auto"),
    the 'repackable' region will be the chain C and the 'designable' region will
    be the set of residues on chain B closer than a 9 ansgtrom cutoff from the
    chain C.
    By default (by setting this parameter to "auto"), the 'score_function' is
    obtained via get_fa_scorefxn().

    PASSO will create 3 reporter files, whose name is based on the input pose
    when applying the protocol:
     1 - Run data file (.JSON): holds information regarding each of the steps
    performed by the protocol. Check self.update_run_data docstring for more
    detailed information.
     2 - Energy file (.DAT): holds energy information regarding the accepted
    steps whose design stage produced mutants with energy values lower than the
    initial conformation. Check self.update_energy_data docstring for more
    detailed information.
     3 - Status file (.TXT): a single line file that gets overwrited every step.
    Simply holds information regarding the current step and the maximum number
    of steps of this simulation. Is used by status.py script to measure the
    completion of the full simulation.
    """

    def __init__(self, n_steps = -1, n_cycles = 2,
        dock_mover = "auto", design_mover = "auto", pre_filter     = "auto",
        designable = "auto", repackable   = "auto", score_function = "auto"):

        from design import get_design_mover
        from pyrosetta import get_fa_scorefxn
        from pyrosetta.rosetta.core.scoring import ScoreFunction
        from pyrosetta.rosetta.protocols.rigid import RigidBodyPerturbMover
        from pyrosetta.rosetta.core.simple_metrics.metrics import \
	        InteractionEnergyMetric
        from pyrosetta.rosetta.core.select.residue_selector import \
            AndResidueSelector, NeighborhoodResidueSelector, ChainSelector, \
            ResidueSelector

        self.n_steps      = n_steps
        self.n_cycles     = n_cycles

        # --- DESIGNABLE REGION
        assert isinstance(designable, ResidueSelector) or designable == "auto",\
            "Designable selection must be a ResidueSelector or set to 'auto'"
        if designable == "auto":
            self.designable = AndResidueSelector(
                NeighborhoodResidueSelector(
                    ChainSelector("C"),
                    9.0,
                    include_focus_in_subset = False),
                ChainSelector("B"))
        else:
            self.designable = designable

        # --- REPACKABLE REGION
        assert isinstance(repackable, ResidueSelector) or repackable == "auto",\
            "Repackable selection must be a ResidueSelector or set to 'auto'"
        if repackable == "auto":
            self.repackable = ChainSelector("C")
        else:
            self.repackable = repackable

        # --- SCORE FUNCTION
        assert type(score_function) == ScoreFunction or \
            score_function == "auto",\
            "Score function must be of type ScoreFunction or set to 'auto'"
        if score_function == "auto":
            self.score_function = get_fa_scorefxn()
        else:
            self.score_function = score_function

        # --- DOCK MOVER
        # In the RigidBodyPerturbMover, the 3 parameters are the following:
        #  First  = Jump number. The perturbation will be applied downstream.
        #  Second = Max rotation in degrees. Based on the movable part centroid.
        #  Third  = Max translation in angstrom.
        if dock_mover == "auto":
            self.dock_mover = RigidBodyPerturbMover(1, 3.0, 0.5)
        else:
            self.dock_mover = dock_mover

        # --- PRE FILTER
        if pre_filter == "auto":
            self.pre_filter = PreFilter()
        else:
            self.pre_filter = pre_filter

        # --- DESIGN MOVER
        if design_mover == "auto":
            self.design_mover = get_design_mover(
                self.score_function, self.designable, self.repackable)
        else:
            self.design_mover = design_mover

        self.ie_metric    = InteractionEnergyMetric(
            self.repackable, self.designable)
        self.data_file    = None # Is defined when applying to a pose
        self.energy_file  = None # Is defined when applying to a pose
        self.status_file  = None # Is defined when applying to a pose


    def update_run_data(self, new_data):
        """
        Updates the self.data_file in JSON format by:
         1- Reading the contents;
         2- Appending the new data;
         3- Writting the updated dictionary in JSON format.
        Each entry on the array corresponds to a single step and holds the
        following data:
         contacts: [val::int,   parameter_passes_PreFilter::boolean]
         distance: [val::float, parameter_passes_PreFilter::boolean]
         clashes : [val::int,   parameter_passes_PreFilter::boolean]
         step    : [val::int,   step_accepted::boolean, acceptance_ratio::float]
        This occurs every step after random rotation + translation.
        """
        
        import json

        with open(self.data_file, "r") as file_in:
            data = json.load(file_in)
        data.append(new_data)
        with open(self.data_file, "w") as file_out:
            json.dump(data, file_out)
            file_out.flush()


    def update_energy_data(self, step, score, interaction):
        """
        Updates the self.energy_file. This occurs when the accepted step design
        stage produces mutants with energy values lower than the initial 
        conformation. Both the 'score' and the 'interaction' energy metrics need
        to be equal or lower than the saved initial values. Both metrics are
        outputed to the energy file.
        """

        with open(self.energy_file, "a") as file_out:
            file_out.write("%5d %12.4f %12.4f\n" % (step, score, interaction))
            file_out.flush()

    def update_status_data(self, step):
        """
        Updates the self.energy_file. This occurs in the beginning of every
        step.
        """

        with open(self.status_file, "w") as file_out:
            file_out.write("%10d %10d" % (step, self.n_steps))
            file_out.flush()


    def apply(self, pose):
        """
        Apply a PASSO protocol to a pose.
        """

        from pyrosetta import Pose
        from ze_utils.pyrosetta_tools import get_pymol_selection_from_selector

        p = Pose(pose)

        # Generate the additional information files
        # 1. Generate the JSON data file.
        self.data_file = p.pdb_info().name()[:-4] + "_data.json"
        with open(self.data_file, "w") as file_out:
            file_out.write("[]")

        # 2. Generate the energy results file.
        self.energy_file = p.pdb_info().name()[:-4] + "_energy.dat"  
        with open(self.energy_file, "w") as file_out:
            file_out.write("%5s %12s %12s\n" % ("Step", "Score", "Interaction"))

        # 3. Generate the current status file.
        self.status_file = p.pdb_info().name()[:-4] + "_status.txt"

        # The initial score and interaction energy metric will be compared to
        # after design efforts: only better structures will be saved for further
        # analysis.
        init_score = self.score_function(p)
        init_interaction = self.ie_metric.calculate(p)
        self.update_energy_data(0, init_score, init_interaction)

        approved_count = 0

        for i in range(1, self.n_steps + 1):
            
            # Print simulation status to a file
            self.update_status_data(i)
            
            # Save pre movement state. The pose will be returned to this
            # conformation in the event that the new structure (after movement)
            # does not pass the pre filter. In such case, design efforts will
            # not be attempted.
            pre_movement_state = Pose(p)

            # Apply the docking movement
            self.dock_mover.apply(p)

            # Verify if the new position is valid for design efforts
            (approved, data) = self.pre_filter.apply(p)
            
            if approved == False:
                p = Pose(pre_movement_state)
                data["step"] = (i, approved, approved_count / i)
                self.update_run_data(data)
                continue
            else:
                approved_count += 1
                data["step"] = (i, approved, approved_count / i)
                self.update_run_data(data)
            
            # ------------------------------------------------------------------
            # Save the pre design state. Regardless of the outcome to the design
            # efforts, the pose will always be returned to this state, saving
            # the original sequence for the next design round.
            pre_design_state = Pose(p)

            # Apply the design mover
            for j in range(self.n_cycles):
                self.design_mover.apply(p)

            # This structure will only be saved if the overall score or the
            # interaction energy metric are equal or better than the orginal
            # state
            score = self.score_function(p)
            interaction = self.ie_metric.calculate(p)
            if (score <= init_score) or (interaction <= init_interaction):

                # Save energy results to file (including the interaction energy)
                self.update_energy_data(i, score, interaction)

                # Print the residue selection that was used for design in this
                # new position
                with open(p.pdb_info().name()[:-4] + "_%d.txt" % (i), "w") as f:
                    f.write(
                        get_pymol_selection_from_selector(self.designable, p))

                # Save PDB file of this sequence + position
                p.dump_pdb(p.pdb_info().name()[:-4] + "_%d.pdb" % (i))

            # Reset for the next simulation step
            p = Pose(pre_design_state) # Reset the structure/sequence
            self.score_function(p)     # Re-score the old structure/sequence


class DockingGrid:
    """
    Creates a simple grid, allowing for PyRosetta poses to the translated to the
    created positions. 
    The grid is centered on 'center', and extends in all 6 directions. The 'd'
    parameters ('xd', 'yd' and 'zd') control the extend in each direction: as an
    example, 'xd' parameter is a tuple with two numbers, the left and right
    extent in the x axis, respectively, starting from the 'center'. The number
    of points between this minimum and maximum positions are defined in the
    'repeats' parameters ('x_repeats', 'y_repeats'and 'z_repeats'), one number
    for each axis. Additionally, if any 'additional_points' are provided, those
    will be added in the beggining of the points list. Using 'orient' function
    will not mess with this addicional points.
    """

    def __init__(self, center = [0, 0, 0], xd=(0, 0), yd=(0, 0), zd=(0, 0),
        x_repeats = 0, y_repeats = 0, z_repeats = 0, additional_points = []):

        self.center    = center
        self.xt        = (xd[0] + xd[1]) / x_repeats
        self.yt        = (yd[0] + yd[1]) / y_repeats
        self.zt        = (zd[0] + zd[1]) / z_repeats
        self.xd        = xd
        self.yd        = yd
        self.zd        = zd
        self.x_repeats = x_repeats
        self.y_repeats = y_repeats
        self.z_repeats = z_repeats
        self.n_ap      = len(additional_points)

        if len(additional_points) > 0:
            self.points = additional_points
        else:
            self.points = []
            
        self.generate_grid_positions()


    def generate_grid_positions(self):
        """
        Generates the grid positions based on the minim/maximum distances and
        number of repeats in each direction.
        """

        import numpy as np

        self.p0 = self.center - np.array([self.xd[0], self.yd[0], self.zd[0]])
        for x_step in range(self.x_repeats):
            for y_step in range(self.y_repeats):
                for z_step in range(self.z_repeats):
                    x = x_step * self.xt
                    y = y_step * self.yt
                    z = z_step * self.zt
                    point = self.p0 + np.array([x, y, z])
                    self.points.append(point)


    def orient(self, comp, vector):
        """
        Rotate the grid so that one of the directions 'comp' of the grid matches
        the given 'vector'. For example, the 'z' component would be [0, 0, 1].
        """

        import numpy as np
        from ze_utils.algebra import get_angle_from_two_vectors, \
            get_rotation_matrix_from_axis_angle, \
            rotate_coords_from_rotation_matrix

        axis        = np.cross(comp, vector)
        angle       = get_angle_from_two_vectors(comp, vector)
        rot_matrix  = get_rotation_matrix_from_axis_angle(axis, angle)
        self.points[self.n_ap:] = rotate_coords_from_rotation_matrix(
            self.points[self.n_ap:],
            rot_matrix,
            self.center)

    
    def translate_selector_to_point(self, point_id, selector, pose):
        """
        Translate (part of) a protein so that it's centroid matches the given
        point in the grid. The selected point is given by its ID on the
        self.points array.
        """

        import numpy as np
        from pyrosetta.rosetta.core.select import \
            get_residues_from_subset
        from ze_utils.pyrosetta_tools import \
            get_centroid_coordinates_from_selector
        from pyrosetta.rosetta.numeric import xyzVector_double_t

        point = self.points[point_id]
        centroid = get_centroid_coordinates_from_selector(selector, pose)
        x, y, z = np.array(point - centroid)
        residue_list = get_residues_from_subset(selector.apply(pose))
        for residue_index in residue_list:
            residue = pose.residue(residue_index)
            for atom_index in range(1, len(residue.atoms()) + 1):
                x1, y1, z1 = list(residue.xyz(atom_index))
                trans = xyzVector_double_t(x1 + x, y1 + y, z1 + z)
                residue.set_xyz(atom_index, trans)
        return pose


    def as_pdb(self, filename):
        """
        Print the current grid in PDB format. Each point is an atom with bare
        minimum information for visualization purposes.
        """

        from ze_utils.molecule_manipulation import Atom

        with open(filename, "w") as file_out:
            for p in self.points:
                a = Atom(chain_name="D", x = p[0], y = p[1], z = p[2])
                file_out.write(a.as_pdb())


class Fragment:
    """
    Holds a continuous fragment of a 'pose', defined as the the residues whose
    ID on the corresponding PDB file are in the given 'ids' list (residue list
    is automatically extracted). The lsit of 'ids' provided does not need to be
    sequential (i.e: [1, 2, 3]), just as long as the the downstream connection
    of residue i is the residue i + 1. Currently, this Fragment class should
    be exclusively employed in loop appending using the self.append_to function.
    """

    def __init__(self, pose, ids):

        from ze_utils.pyrosetta_tools import get_residues_from_selector, ID2Rank
        from pyrosetta.rosetta.core.select.residue_selector import \
            ResidueIndexSelector

        self.pose = pose
        self.ids = ids
        
        # The Fragment to be appended receives a list of IDs, as defined in the
        # corresponding PDB file. However, ResidueIndexSelector works with rank,
        # not with ID. Therefore, a dictionary converting the given ID's to the
        # respective rank position in the pose must be created.
        residues_selector = ResidueIndexSelector(ID2Rank(pose, ids))
        self.residues = get_residues_from_selector(residues_selector, pose)

        # Residues selected will be appended in a polymer (one after the other),
        # in the exact given order. Therefore, in order to prevent unexpected
        # results, all residues provided should be in a constinuous
        # sequence/chain in the reference pose, where the downstream connection
        # of residue i is always the residue i + 1.
        for i in range(len(self.residues) - 1):
            downstream_residue = self.residues[i].connected_residue_at_upper()
            downstream_id = ID2Rank(pose, [ids[i + 1]], as_string = False)[0]
            assert downstream_residue == downstream_id, \
                "Residue %d is not the downstream residue of residue %d." % \
                    (ids[i + 1], ids[i])


    def append_to(self, pose, pos, backbone = "auto"):
        """
        Appends this fragment to a given 'pose', by attaching the sequence of
        residues at the residue in 'pos' (this value should be the rank in the
        pose, not the ID on the corresponding PDB file). A mode of backbone
        construction can be provided. The following modes are available:
         - stretched : All phi and psi angles are set to 180.0 degrees.
         - alpha     : Phi and psi angles are set to ideal alpha helix values
         - beta      : Phi and psi angles are set to ideal beta sheet values
         - auto      : Phi, psi and omega angles are set to the orginal values
                       in the fragment pose
        """

        assert type(pos) == int and pos > 0 and pos <= len(pose.sequence()), \
            "'pos' must be an int between 0 and the number of residues in pose"

        assert backbone in ["stretched", "alpha", "beta", "auto"], \
            "backbone option must be set to stretched, alpha, beta or auto."

        # Default values for non "auto" instructions.
        angles = {
            #              Phi     Psi     Omega
            "stretched": [ 180.0,  180.0, 180.0],
            "alpha":     [ -60.0,  -50.0, 180.0],
            "beta":      [-140.0,  130.0, 180.0]
        }

        # Save anchor residue phi/psi/omega as in the reference pose. When
        # appeding polymer residues, the anchor residues psi and omega angles\
        # are modified (for some reason). Therefore, by saving the reference
        # values, it's possible to fix this.
        anchor_residue_index = self.residues[0].connected_residue_at_lower()
        phi                  = self.pose.phi(  anchor_residue_index)
        psi                  = self.pose.psi(  anchor_residue_index)
        omega                = self.pose.omega(anchor_residue_index)
        saved_angles         = [[phi, psi, omega]]

        # To make iteration easier over the residues who require post-addition
        # rotation of the phi, psi and omega angles, save all positions added to
        # the pose
        residues_to_modify = [pos]

        # Add, one by one, the reference residues. Appended residues have, by
        # default, weird phi, psi and omega values. Save the reference values or
        # the new alpha, beta or stretched confromation values from 'angles'
        # dictionary for post-addition rotation.
        for residue in self.residues:

            # If backbone is set to 'auto', save the phi, psi and omega angles
            # to recover. However, the final conformation of the added residues
            # will, sometimes, not perfectly match the given reference. This is
            # because the appeding method for adding new residues to the pose
            # sets 3 atom angles to be the default idealized value, and not the
            # pre-existing value.
            if backbone == "auto":
                phi   = self.pose.phi(  residue.seqpos())
                psi   = self.pose.psi(  residue.seqpos())
                omega = self.pose.omega(residue.seqpos())
                saved_angles.append([phi, psi, omega])

            # If backbone is set to 'alpha', 'beta' or 'stretched', save the
            # corresponding value from the 'angles' dictionary.
            else:
                saved_angles.append(angles[backbone])
            
            # Add this residue to the anchor position
            pose.append_polymer_residue_after_seqpos(residue, pos, True)

            # Update the anchor position, while marking this new residue for
            # post-addition rotation of the phi, psi and omega angles
            pos += 1
            residues_to_modify.append(pos)

        # Post-addition rotation. For some reason, when appeding a residue the
        # previous residue on the chain has its psi and omega angles modified
        # (to 150.0 and 0.0 degrees, respectively). Therefore, the rotation to
        # the correct angles (either the reference of new alpha/beta/stretched
        # values) is only performed after all residues have been appended.
        for ang, r in zip(saved_angles, residues_to_modify):
            pose.set_phi(  r, ang[0])
            pose.set_psi(  r, ang[1])
            pose.set_omega(r, ang[2])


class ResidueToAtomMapping:
    """
    This simple class holds information regarding residue to atom mapping,
    regarding the sidechains of each residue:
     - 'reg' entry refers to the last atom of the sidechain of each aminoacid;
     - 'cen' entry refers to the CEN atom in centroid mode for each aminoacid;
     - 'cal' entry refers to the alpha carbon in each aminoacid;
    """

    reg = { "A":  "CB", "R": "NH2", "N":  "CG", "D":  "CG",
            "C":  "SG", "Q":  "CD", "E":  "CD", "G":  "CA",
            "H": "CE1", "I": "CD1", "L": "CD1", "K":  "NZ",
            "M":  "CE", "F":  "CZ", "P":  "CA", "S":  "OG",
            "T": "CG2", "W": "CZ3", "Y":  "OH", "V": "CG1"}
    cen = { "A": "CEN", "R": "CEN", "N": "CEN", "D": "CEN",
            "C": "CEN", "Q": "CEN", "E": "CEN", "G": "CEN",
            "H": "CEN", "I": "CEN", "L": "CEN", "K": "CEN",
            "M": "CEN", "F": "CEN", "P": "CEN", "S": "CEN",
            "T": "CEN", "W": "CEN", "Y": "CEN", "V": "CEN"}
    cal = { "A":  "CA", "R":  "CA", "N":  "CA", "D":  "CA",
            "C":  "CA", "Q":  "CA", "E":  "CA", "G":  "CA",
            "H":  "CA", "I":  "CA", "L":  "CA", "K":  "CA",
            "M":  "CA", "F":  "CA", "P":  "CA", "S":  "CA",
            "T":  "CA", "W":  "CA", "Y":  "CA", "V":  "CA"}