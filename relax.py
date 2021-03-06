# Questions: jose.manuel.pereira@ua.pt

import argparse
from pyrosetta import *
from relax_decoy import get_relaxer_mover
from ze_utils.common import \
    get_number_of_jobs_in_slurm_queue, verify_launch_failed_requeued_held, \
    overwrite_dir


#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\

#                                 Relax Script:
# ______________________________________________________________________________
#  Performs an energy minimization. This script calls relax_decoy.py script to
# perform each decoy of the simulation. Read its documentation and docstrings
# for more detailed information.


class DEFAULT:
    """
    Define defaults for the script. Values can be modified using arguments.
    """
    n_decoys       = 100
    c_weight       = 1.0
    max_slurm_jobs = 499
    no_constraints = False
    keep_bash      = False
    output_prefix  = "relax"
    partition      = "main"


def dump_sbatch_script(input_file, output_file = DEFAULT.output_prefix,
    no_constrains = DEFAULT.no_constraints, c_weight = DEFAULT.c_weight,
    partition = DEFAULT.partition, keep_bash = DEFAULT.keep_bash):
    """
    Creates a new bash file at 'output_name.sh' that instructs the SLURM
    framework to spawn a new relax_recoy. The location of this script is, by
    default, '~/scripts/relax_decoy.py', but can be changed bellow. If
    'keep_bash' is set to False, the created bash file will be deleted after
    being run with the sbatch or srun command.
    """

    bash_filename = "%s.sh" % (output_file)
    with open(bash_filename, "w") as bash:
        bash.write("#!/bin/bash\n")
        bash.write("#SBATCH --partition=%s\n" % (partition))
        bash.write("#SBATCH --ntasks=1\n")
        bash.write("#SBATCH --cpus-per-task=1\n")
        bash.write("#SBATCH --mem=2GB\n")
        bash.write("#SBATCH --time=72:00:00\n")
        bash.write("#SBATCH --requeue\n")
        bash.write("#SBATCH --job-name=%s\n" % (output_file))
        bash.write("#SBATCH --output=%s.out\n" % (output_file))
        bash.write("#SBATCH --error=%s.err\n\n" % (output_file))
        bash.write("python ~/scripts/relax_decoy.py %s" % (input_file))
        bash.write(" -o %s" % (output_file))
        
        if no_constrains:
            bash.write(" -nc -w %.3f" % (c_weight))
        
        # Auto-deletes the current bash script
        if not keep_bash:
            bash.write("\nrm -rf %s" % (bash_filename))


def deploy_decoys_on_slurm(input_file, output_prefix = DEFAULT.output_prefix,
    n_decoys = DEFAULT.n_decoys, no_constraints = DEFAULT.no_constraints,
    c_weight = DEFAULT.c_weight, fasc_filename = "auto",
    partition = DEFAULT.partition):
    """
    Launch the relax protocol decoys in parallel mode, using SLURM on AMAREL
    framework. This function will remain ad infinitum waiting for processors to
    become available for the user (this value is set to 499 on
    DEFAULT.max_slurm_jobs, but can be modified if using a different framework
    than AMAREL). This function will output a new bash script for each decoy
    and then automatically run it when the resources become available before
    skipping to the next decoy. Output and Error files will be automatically
    stored in /out/ and /err/ directories, respectively. After all decoys are
    complete, the resulting PDB files will be loaded and their energies
    measured, saving the results to a custom .FASC file.
    """

    import os
    import json
    import shutil

    # Create the fasc filename if set to auto, else verify extension
    if fasc_filename == "auto":
        fasc_filename = "%s_custom.fasc" % (output_prefix)
    else:
        assert fasc_filename[-5:] == ".fasc", \
            "Parameter 'fasc_filename' must have '.fasc' extension (%s)" % \
                (fasc_filename[-5:])

    # The output and error files will be saved in the corresponding directory
    overwrite_dir("out")
    overwrite_dir("err")

    # Print a template of the executed bash scripts, for archive purposes
    output_name = "%s_%s" % (output_prefix, "template")
    dump_sbatch_script(input_file, output_name, no_constraints, c_weight,
        partition, True)

    user = os.environ['USER']
    for decoy_index in range(n_decoys):
        
        # 1) Create the single_relax.py sbatch script
        output_name = "%s_%d" % (output_prefix, decoy_index)
        dump_sbatch_script(input_file, output_name, no_constraints, c_weight,
            partition)

        # 2) Ping the user slurm queue for a vacant spot to deploy the job
        while True:
            n_jobs_on_slurm_queue = get_number_of_jobs_in_slurm_queue(user)

            if n_jobs_on_slurm_queue <= DEFAULT.max_slurm_jobs:
                os.system("sbatch %s" % ("%s.sh" % (output_name)))
                print("-> Decoy %d/%-d" % (decoy_index + 1, n_decoys))
                break

    print(" All decoys spawned. Waiting for completion ...")
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
            print(" All decoys completed! Creating custom fasc file ...")
            break

    # 4) Create a custom .fasc file a posteriori
    score_function = get_fa_scorefxn()
    with open(fasc_filename, "w") as fasc:
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


def deploy_decoys_on_pyjobdistributor(input_file, output_prefix, n_decoys,
    no_constrains = False, c_weight = 1.0):
    """
    Launch the relax protocol decoys in serial mode, using the PyJobDistributor.
    Each decoy will occupy one processor on the machine, and a new decoy will be
    started when a processor becomes available.

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

    job_man = PyJobDistributor(output_prefix, n_decoys, score_function)
    pose = pose_from_pdb(input_file)
    p = Pose()
    relaxer = get_relaxer_mover(pose, score_function, no_constrains, c_weight)

    while not job_man.job_complete:
        p.assign(pose)
        relaxer.apply(p)
        job_man.output_decoy(p)


def validate_arguments(args):
    """
    Validate the script arguments to prevent known input errors.
    """
    if not args.input_file[-4:] == ".pdb":
        exit("ERROR: Input file should be in PDB format")
    if args.n_decoys < 0:
        exit("ERROR: Number of decoys must be a non-negative value")
    if args.c_weight < 0 or args.c_weight > 1.0:
        exit("ERROR: Constraints weight must be a non-negative value (< 1.0)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Relax a PDB structure using
        multiple decoys.""")
    parser.add_argument('input_file', metavar='INPUT', type=str,
        help='The input PDB file')
    parser.add_argument('-nd', '--n_decoys', metavar='', type=int,
        help='Number of decoys (Default: %d)' % (DEFAULT.n_decoys),
        default = DEFAULT.n_decoys)
    parser.add_argument('-nc', '--no_constraints',
        action='store_true', help='Turn off constraints (Default: False)')
    parser.add_argument('-w', '--c_weight', metavar='', type=int,
        help='Constraints weight (Default: %5.2f)' % (DEFAULT.c_weight),
        default = DEFAULT.c_weight)
    parser.add_argument('-o', '--output', metavar='', type=str,
        help='Output prefix (Default: %s)' % (DEFAULT.output_prefix),
        default = DEFAULT.output_prefix)
    parser.add_argument('-s', '--slurm', action = 'store_true',
        help = "Use SLURM to launch parallel decoys (Default: False)")
    parser.add_argument('-p', '--partition', metavar='', type=str,
        help='Partition for SLURM (Default: %s)' % (DEFAULT.partition),
        default = DEFAULT.partition)

    args = parser.parse_args()
    validate_arguments(args)
    init()

    if args.slurm:
        deploy_decoys_on_slurm(
            args.input_file,
            output_prefix  = args.output,
            n_decoys       = args.n_decoys,
            no_constraints = args.no_constraints,
            c_weight       = args.c_weight,
            partition      = args.partition)
    else:
        deploy_decoys_on_pyjobdistributor(
            args.input_file,
            args.output,
            args.n_decoys,
            args.no_constraints,
            args.c_weight)