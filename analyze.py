import argparse
from pyrosetta import *
from ze_utils.pyrosetta_tools import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.core.select.residue_selector import *
from pyrosetta.rosetta.core.simple_metrics.metrics import \
	InteractionEnergyMetric, SasaMetric


#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\


#                            Analyze Script:
# ______________________________________________________________________________
#  Analizes a single PDB file and outputs revelant intra-pose metrics:
#    . Total energy
#    . Individual residues energy per selection
#    . Interaction energy metric (If a ligand is present)
#    . Solvent Acessible Surface Area (SASA) Metric
#    . Hydrogen bonding network energy
#
#  Necessary parts:
# - Starting Pose;
# - Score function;
# - SimpleMetrics (Are used to calculate a certain characteristic of the system.
#   They can also be used during the simulation and exported to the FASC file)
#   For more information, read the RosettaCommons page:
#   /scripting_documentation/RosettaScripts/SimpleMetrics/SimpleMetrics;
#   And the following webpages:
#   http://www.programmersought.com/article/77091668895/
# - Selections JSON *Script Specific* (This script should contain the number of 
#   the residues included in each selection in a string, divided by a comma and 
#   without blank spaces. Since the format is in JSON, the contents of the file
#   should be a dictionary where each entry should be "sel_name":"1,2,3,4,5,6".
#   Altough an arbitrary number or selections can be introduced, the interaction
#   energy metric will be calculated between two specific selection_names if
#   they are found: "interface" and "lig");

#                  Rosetta Common Latest Documentation:
# https://www.rosettacommons.org/docs/latest (*add the topic here*)

class Data():
    """
    Save relevant Data from the PDB Pose, gathered from SimpleMetrics.
    """
    def __init__(self, filename):
        self.filename           = filename
        self.total_energy       = 0.0
        self.hbne               = 0.0
        self.interaction        = None
        self.res_e              = {}
        self.sasa               = {}

    def __str__(self):
        return self.as_string()

    def as_string(self):
        """
        Print the Data content as a string.
        """
        s = ""
        s += "%40s\n" % (self.filename)
        s += "%40s: %12.4f\n" % ("Total energy", self.total_energy)
        if not self.interaction == None:
            s += "%40s: %12.4f\n" % ("Interaction energy", self.interaction)
        s += "%40s: %12.4f\n" % ("Hydrogen bonding network energy", self.hbne)
        s += "%40s:\n" % ("Solvent Accessible Surface Area (SASA)")
        for selection in self.sasa:
            s += "%39s:  %12.4f\n" % (selection, self.sasa[selection])
        s += "%40s:\n" % ("Individual residue total energies")
        for sel in self.res_e:
            s += "%39s:\n" % (sel)
            for (_, i, e) in self.res_e[sel]:
                s += "%38d    %12.4f\n" % (i, e)
        return s


class DEFAULT:
    """
    Define defaults for the script. Values can be modified using arguments.
    """
    selections  = "selections.json"
    output_file = "analyze.dat"


def validate_arguments(args):
    """
    Validate the script arguments to prevent known input errors.
    """
    if not args.input_file[-4:] == ".pdb":
        exit("ERROR: Input file should be in PDB format")
    if not args.selections[-5:] == ".json":
        exit("ERROR: Selections file should be in JSON format")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Analyze a single PDB file,
        outputing intra-structure metrics.""")
    parser.add_argument('input_file', metavar='INPUT', type=str,
        help='The input PDB file')
    parser.add_argument('-s', '--selections', metavar='', type=str,
        help='Selections JSON file (Default: %s)' % (DEFAULT.selections),
        default = DEFAULT.selections)
    parser.add_argument('-o', '--output_file', metavar='', type=str,
        help='Output file (Default: %s)' % (DEFAULT.output_file),
        default = DEFAULT.output_file)

    args = parser.parse_args()
    validate_arguments(args)
    init()

    pose           = pose_from_pdb(args.input_file)
    score_function = get_fa_scorefxn()
    selections     = load_selections_from_json_file(args.selections, pose)
    data           = Data(args.input_file)

    # 1. Total Energy
    data.total_energy = score_function(pose)

    # 2. Individual Residues Energy per Selection
    for selection in selections:
        if not selections[selection] == None:
            data.res_e[selection] = get_residue_energies_from_selector(
                selections[selection],
                pose,
                verbose = False)

    # 3. Interaction Energy Metric
    # This should only be calculated if the selections has a ligand.
    if not (selections["lig"] == None or selections["interface"] == None):
        ie = InteractionEnergyMetric(selections["lig"], selections["interface"])
        data.interaction = ie.calculate(pose)

    # 4. Solvent Acessible Surface Area (SASA) Metric
    sasa_metric = SasaMetric()
    for selection in selections:
        if not selections[selection] == None:
            sasa_metric.set_residue_selector(selections[selection])
            data.sasa[selection] = sasa_metric.calculate(pose)

    # 5. Hydrogen Bonding Network Energy
    score_terms = [hbond_sr_bb, hbond_lr_bb, hbond_bb_sc, hbond_sc]
    for score_term in score_terms:
        data.hbne += pose.energies().total_energies()[score_term]

    with open(args.output_file, "w") as file_out:
        file_out.write(data.as_string())
    print("\n\nDone!")