import argparse
from pyrosetta import *
from ze_utils.pyrosetta_tools import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.core.select.residue_selector import \
    ChainSelector, NeighborhoodResidueSelector
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


class DEFAULT:
    """
    Define defaults for the script. Values can be modified using arguments.
    """
    cut_off = 9.0
    chain   = 'C'


def validate_arguments(args):
    """
    Validate the script arguments to prevent known input errors.
    """
    if not args.input_file[-4:] == ".pdb":
        exit("ERROR: Input file should be in PDB format")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Analyze a single PDB file,
        outputing intra-structure metrics.""")
    parser.add_argument('input_file', metavar='INPUT', type=str,
        help='The input PDB file')
    parser.add_argument('-co', '--cutoff', metavar='', type=float,
        help='Interface identification cut-off in Å (Default: %5.1f)' \
            % (DEFAULT.cut_off), default = DEFAULT.cut_off)
    parser.add_argument('-c', '--chain', metavar='', type=str,
        help='Chain defined as index (Default: %s)' \
            % (DEFAULT.chain), default = DEFAULT.chain)

    args = parser.parse_args()
    validate_arguments(args)
    init()

    pose           = pose_from_pdb(args.input_file)
    score_function = get_fa_scorefxn()

    # 1. Total Energy
    # This needs to be the first step, as the score needs to be scored, for
    # example, for the InteractionEnergyMetric calculation.
    total_energy = score_function(pose)


    # 3. Interaction Energy Metric
    # This should only be calculated if a chain C is present on the pose.
    ligand_residues, interface_residues = [], []
    interaction_energy = None
    ligand = ChainSelector(args.chain)
    interface = NeighborhoodResidueSelector(ligand, args.cutoff,
        include_focus_in_subset = False)
    if len(get_residues_from_selector(ligand, pose)) > 0:
        interaction_energy_metric = InteractionEnergyMetric(ligand, interface)
        interaction_energy        = interaction_energy_metric.calculate(pose)

        # 2. Individual Residues Energy per Selection
        ligand_residues    = get_residue_energies_from_selector(ligand, pose)
        interface_residues = get_residue_energies_from_selector(interface, pose)

    # 4. Hydrogen Bonding Network Energy
    hydrogen_bonding_energy = 0.0
    score_terms = [hbond_sr_bb, hbond_lr_bb, hbond_bb_sc, hbond_sc]
    for score_term in score_terms:
        hydrogen_bonding_energy += pose.energies().total_energies()[score_term]


    print("\n\n%30s %10s\n %s" % ("Metric", "Value", "-"*40))
    print("%30s %10.3f" % ("Total energy", total_energy))
    print("%30s %10.3f" % ("Hydrogen bonding energy", hydrogen_bonding_energy))
    if not interaction_energy == None:
        print("%30s %10.3f" % ("Interaction energy", interaction_energy))
    if len(ligand_residues) > 0:
        print("\n\n%30s %10s\n %s" % ("Ligand residue index", "Energy", "-"*40))
        for (_, index, energy) in ligand_residues:
            print("%30s %10.3f" % (index, energy))
    if len(interface_residues) > 0:
        print("\n\n%30s %10s\n %s" % ("Interface residue index", "Energy", "-"*40))
        for (_, index, energy) in interface_residues:
            print("%30s %10.3f" % (index, energy))
