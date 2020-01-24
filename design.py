# Questions: jose.manuel.pereira@ua.pt

import argparse
from pyrosetta import *
from design_decoy import get_designer_mover, get_designable_region
from ze_utils.common import \
    get_number_of_jobs_in_slurm_queue, verify_launch_failed_requeued_held, \
    overwrite_dir
from ze_utils.pyrosetta_tools import get_pymol_selection_from_selector


#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\

#                                 Design Script:
# ______________________________________________________________________________
#  Performs a sequence design. This script calls design_decoy.py script to
# perform each decoy of the simulation. Read its documentation and docstrings
# for more detailed information.


class DEFAULT:
    """
    Define defaults for the script. Values can be modified using arguments.
    """
    ligand_chain   = "C"
    blocked_region = None
    cutoff         = 9.0
    n_cycles       = 4
    output_prefix  = "design"
    score_function = "auto"
    designable     = "auto"
    repackable     = "auto"
    sele_file    = None

    n_decoys       = 100
    max_slurm_jobs = 499
    keep_bash      = False
    fasc_file      = "auto"


def dump_sbatch_script(input_file, output_file = DEFAULT.output_prefix,
    cutoff = DEFAULT.cutoff, ligand_chain = DEFAULT.ligand_chain,
    blocked_region = DEFAULT.blocked_region, n_cycles = DEFAULT.n_cycles,
    keep_bash = DEFAULT.keep_bash):
    """
    Creates a new bash file at 'output_name.sh' that instructs the SLURM
    framework to spawn a new design_recoy. The location of this script is, by
    default, '~/scripts/design_decoy.py', but can be changed bellow. If
    'keep_bash' is set to False, the created bash file will be deleted after
    being run with the sbatch or srun command.
    """

    bash_filename = "%s.sh" % (output_file)
    with open(bash_filename, "w") as bash:
        bash.write("#!/bin/bash\n")
        bash.write("#SBATCH --partition=main\n")
        bash.write("#SBATCH --ntasks=1\n")
        bash.write("#SBATCH --cpus-per-task=1\n")
        bash.write("#SBATCH --mem=4GB\n")
        bash.write("#SBATCH --time=72:00:00\n")
        bash.write("#SBATCH --requeue\n")
        bash.write("#SBATCH --job-name=%s\n" % (output_file))
        bash.write("#SBATCH --output=%s.out\n" % (output_file))
        bash.write("#SBATCH --error=%s.err\n\n" % (output_file))
        bash.write("python ~/scripts/design_decoy.py %s" % (input_file))
        bash.write(" -o %s" % (output_file))

        bash.write(" -lc %s -b %s" % (ligand_chain, blocked_region))
        bash.write(" -co %.3f -nc %d" % (cutoff, n_cycles))
        
        # Auto-deletes the current bash script
        if not keep_bash:
            bash.write("\nrm -rf %s" % (bash_filename))


def deploy_decoys_on_slurm(input_file, output_prefix = DEFAULT.output_prefix,
    cutoff = DEFAULT.cutoff, ligand_chain = DEFAULT.ligand_chain,
    blocked_region = DEFAULT.blocked_region, n_cycles = DEFAULT.n_cycles,
    n_decoys = DEFAULT.n_decoys, fasc_file = DEFAULT.fasc_file,
    sele_file = DEFAULT.sele_file):
    """
    Launch the design protocol decoys in parallel mode, using SLURM on AMAREL
    framework. This function will remain ad infinitum waiting for processors to
    become available for the user (this value is set to 499 on
    DEFAULT.max_slurm_jobs, but can be modified if using a different framework
    than AMAREL). This function will output a new bash script for each decoy
    and then automatically run it when the resources become available before
    skipping to the next decoy. Output and Error files will be automatically
    stored in /out/ and /err/ directories, respectively. After all decoys are
    complete, the resulting PDB files will be loaded and their energies
    measured, saving the results to a custom .FASC file.

    Note: With this method for decoy deployment, no custom designable/repackable
    regions can be passed. Default 'auto' regions will be employed instead. Read
    the design_decoy.py documentation and docstrings for more information.
    Future upgrades to the script would allow this regions to be printed to a
    file then loaded by the design_decoy.py main script, such as the PreFilter
    in JSON format, on single_dock.py > single_dock_decoy.py scripts.
    """

    import os
    import json
    import shutil

    # Create the fasc filename if set to auto, else verify extension
    if fasc_file == "auto":
        fasc_file = "%s_custom.fasc" % (output_prefix)
    else:
        assert fasc_file[-5:] == ".fasc", \
            "Parameter 'fasc_filename' must have '.fasc' extension (%s)" % \
                (fasc_file[-5:])

    # The output and error files will be saved in the corresponding directory
    overwrite_dir("out")
    overwrite_dir("err")

    # Print a template of the executed bash scripts, for archive purposes
    output_name = "%s_%s" % (output_prefix, "template")
    dump_sbatch_script(input_file, output_name, cutoff, ligand_chain,
        blocked_region, n_cycles, True)

    # Print the designed region to a file, if sele_file is defined and the file
    # doesn't exist
    if sele_file != None and not os.path.exists(sele_file):
        with open(sele_file, "w") as file_out:
            design_reg = get_designable_region("auto", cutoff, ligand_chain,
                blocked_region)
            pose = pose_from_pdb(input_file)
            file_out.write(get_pymol_selection_from_selector(design_reg, pose))

    user = os.environ['USER']
    for decoy_index in range(n_decoys):
        
        # 1) Create the single_relax.py sbatch script
        output_name = "%s_%d" % (output_prefix, decoy_index)
        dump_sbatch_script(input_file, output_name, cutoff, ligand_chain,
            blocked_region, n_cycles)

        # 2) Ping the user slurm queue for a vacant spot to deploy the job
        while True:
            n_jobs_on_slurm_queue = get_number_of_jobs_in_slurm_queue(user)

            if n_jobs_on_slurm_queue <= DEFAULT.max_slurm_jobs:
                os.system("sbatch %s" % ("%s.sh" % (output_name)))
                break

    # 3) Yield until all decoys are completed
    while True:
        to_break = True
        for decoy_index in range(n_decoys):
            output_name = "%s_%d" % (output_prefix, decoy_index)
            if os.path.exists("%s.sh" % (output_name)):
                verify_launch_failed_requeued_held(user)
                to_break = False
                break
        if to_break:
            break

    # 4) Create a custom .fasc file a posteriori
    score_function = get_fa_scorefxn()
    with open(fasc_file, "w") as fasc:
        for decoy_index in range(n_decoys):
            inputname = "%s_%d.pdb" % (output_prefix, decoy_index)
            pose      = pose_from_pdb(inputname)
            score     = score_function(pose)
            entry     = json.dumps(
                {"pdb_name"  : inputname[:-4],
                "decoy"      : inputname,
                "total_score": score})
            fasc.write("%s\n" % (entry))
            
            # 5) Store the output and error files in the corresponding folders
            shutil.move("%s_%d.err" % (output_prefix, decoy_index), "err")
            shutil.move("%s_%d.out" % (output_prefix, decoy_index), "out")


def deploy_decoys_on_pyjobdistributor(input_file,
    output_prefix = DEFAULT.output_prefix, designable = DEFAULT.designable,
    repackable = DEFAULT.repackable, cutoff = DEFAULT.cutoff,
    ligand_chain = DEFAULT.ligand_chain, blocked_region = DEFAULT.blocked_region,
    n_cycles = DEFAULT.n_cycles, n_decoys = DEFAULT.n_decoys,
    sele_file = DEFAULT.sele_file):
    """
    Launch the design protocol decoys in serial mode, using the
    PyJobDistributor. Each decoy will occupy one processor on the machine, and a
    new decoy will be started when a processor becomes available.

    Note: As of December 2019, a bug seems to be present where the decoys are
    being simulated by a single processor only, in serial mode.
    """

    # Import pyrosetta if the calling script doesn't have it imported already
    try:
        score_function = get_fa_scorefxn()
    except:
        import pyrosetta
        init()
        score_function = get_fa_scorefxn()

    job_man  = PyJobDistributor(output_prefix, n_decoys, score_function)
    pose     = pose_from_pdb(input_file)
    # Optional but recommended for looking for NeighborhoodResidueSelector
    pose.update_residue_neighbors() 

    designer = get_designer_mover(score_function, designable, repackable,
        cutoff, ligand_chain, blocked_region)

    # Print the designed region to a file, if sele_file is defined and the file
    # doesn't exist
    if sele_file != None and not os.path.exists(sele_file):
        with open(sele_file, "w") as file_out:
            design_reg = get_designable_region(designable, cutoff, ligand_chain,
                blocked_region)
            file_out.write(get_pymol_selection_from_selector(design_reg, pose))
    
    p        = Pose()
    while not job_man.job_complete:
        p.assign(pose)
        for i in range(n_cycles):
            designer.apply(p)
        job_man.output_decoy(p)


def validate_arguments(args):
    """
    Validate the script arguments to prevent known input errors.
    """
    if not args.input_file[-4:] == ".pdb":
        exit("ERROR: Input file should be in PDB format")
    if args.cutoff < 0:
        exit("ERROR: Cut-off must be a non-negative value")
    if args.n_decoys < 0:
        exit("ERROR: Number of decoys must be a non-negative value")
    if args.n_cycles < 0:
        exit("ERROR: Number of cycles must be a non-negative value")
    if args.blocked_region != None:
        assert len(args.blocked_region.split("-")) == 2, \
            "Blocked region should be in format 'int-int'"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Design the active-site of a
        protein surrousing a ligand defined in a chain (within the provided
        cutoff), on the input PDB, using multiple n_decoys. If a blocked region
        is provided, no design nor repack will be performed on those residues. A
        blocked region is defined with the following syntax: "start-end", where
        start and end are residue index numbers. For example, blocking the
        region from residue 1 to 56, one would use "-bl 1-56". A total of
        n_cycles will be attempted back to back to increase convergence. If the
        slurm flag is provided, decoys will be deployed in the SLURM framework.
        If a sele_file is provided, the region determined for design will be
        printed to that file, if it doesn't exist.""")
    parser.add_argument('input_file', metavar='INPUT', type=str,
        help='The input PDB file')
    parser.add_argument('-o', '--output', metavar='', type=str,
        help='Output prefix (Default: %s)' % (DEFAULT.output_prefix),
        default = DEFAULT.output_prefix)
    parser.add_argument('-lc', '--ligand_chain', metavar='', type=str,
        help='The ligand chain identifier (Default: %s)' % \
        (DEFAULT.ligand_chain), default = DEFAULT.ligand_chain)
    parser.add_argument('-b', '--blocked_region', metavar='', type=str,
        help='The blocked (no design) region (Default: %s)' % \
        (DEFAULT.blocked_region), default = DEFAULT.blocked_region)
    parser.add_argument('-co', '--cutoff', metavar='', type=float,
        help='Interface identification cut-off in Å (Default: %5.1f)' \
            % (DEFAULT.cutoff), default = DEFAULT.cutoff)
    parser.add_argument('-nc', '--n_cycles', metavar='', type=int,
        help='Number of pack/min cycles (Default: %d)' % (DEFAULT.n_cycles),
        default = DEFAULT.n_cycles)
    parser.add_argument('-nd', '--n_decoys', metavar='', type=int,
        help='Number of decoys to deploy (Default: %d)' % (DEFAULT.n_decoys),
        default = DEFAULT.n_decoys)
    parser.add_argument('-s', '--slurm', action = 'store_true',
        help = "Use SLURM to launch parallel decoys (Default: False)")
    parser.add_argument('-sele', '--sele_file', metavar='', type=str,
        help='Selection output file name (Default: %s)' % (DEFAULT.sele_file),
        default = DEFAULT.sele_file)

    args = parser.parse_args()
    validate_arguments(args)
    init()

    if args.slurm:
        deploy_decoys_on_slurm(
            args.input_file,
            output_prefix  = args.output,
            cutoff         = args.cutoff,
            ligand_chain   = args.ligand_chain,
            blocked_region = args.blocked_region,
            n_cycles       = args.n_cycles,
            n_decoys       = args.n_decoys,
            sele_file      = args.sele_file)
    else:
        deploy_decoys_on_pyjobdistributor(
            args.input_file,
            args.output,
            DEFAULT.designable,
            DEFAULT.repackable,
            args.cutoff,
            args.ligand_chain,
            args.blocked_region,
            args.n_cycles,
            args.n_decoys,
            args.sele_file)