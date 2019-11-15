import argparse
from pyrosetta import *
from ze_utils.pyrosetta_classes import PASSO
from single_dock_decoy import single_dock_decoy
from ze_utils.pyrosetta_tools import set_ABC_model_fold_tree
from ze_utils.common import get_number_of_jobs_in_slurm_queue

#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\

#                          Dock-Design Script:
# ______________________________________________________________________________
#  Performs docking + design algorithm. Please read the single_dock_decoy.py
# script and the PASSO class (at ze_utils.pyrosetta_classes) docstrings for more
# detailed information.
#
#  > Check multi_dock.py script for multiple starting positions PASSO simulation

class DEFAULT:
    """
    Define defaults for the script. Values can be modified using arguments.
    """
    
    n_decoys       = 100
    n_steps        = 2000
    max_slurm_jobs = 499


def validate_arguments(args):
    """
    Validate the script arguments to prevent known input errors.
    """
    if not args.input_file[-4:] == ".pdb":
        exit("ERROR: Input file should be in PDB format")
    if args.n_decoys < 0:
        exit("ERROR: Number of decoys must be a non-negative value")
    if args.n_steps < 0:
        exit("ERROR: Number of PASSO steps must be a non-negative value")


def dump_sbatch_script(output_name, input_file, n_steps):
    """
    Creates a new bash file at 'output_name'.sh path that instructs the SLURM
    framework to spawn a new single_dock_decoy. The location of this script is,
    by default, '~/scripts/single_dock_decoy.py', but can be changed bellow. In
    the context of the multi_dock.py script, the starting Pose will be obtained
    from the 'start_X.pdb' file, where X is the current 'point_id' in the
    docking grid, but a custom 'input_file' can be provided.
    Each decoy simulation will perform 'n_steps' of the PASSO protocol.
    """

    with open("%s.sh" % (output_name), "w") as bash:
        bash.write("#!/bin/bash\n")
        bash.write("#SBATCH --partition=main\n")
        bash.write("#SBATCH --ntasks=1\n")
        bash.write("#SBATCH --cpus-per-task=1\n")
        bash.write("#SBATCH --mem=4GB\n")
        bash.write("#SBATCH --time=72:00:00\n")
        bash.write("#SBATCH --requeue\n")
        bash.write("#SBATCH --job-name=%s\n" % (output_name))
        bash.write("#SBATCH --output=%s.out\n" % (output_name))
        bash.write("#SBATCH --error=%s.err\n\n" % (output_name))
        bash.write("python ~/scripts/single_dock_decoy.py %s" % (input_file))
        bash.write(" -ns %d -o %s" % (n_steps, output_name))


def deploy_decoys_on_slurm(input_file, output_prefix, n_decoys, n_steps):
    """
    Launch the PASSO protocol decoys in parallel mode, using SLURM on AMAREL
    framework. This function will remain ad infinitum waiting for processors to
    become available for the user (this value is set to 499 on
    DEFAULT.max_slurm_jobs, but can be modified if using a different framework
    than AMAREL). This function will output a new bash script for each decoy
    and then automatically run it when the resources become available before
    skipping to the next decoy.
    """

    user = os.environ['USER']
    for decoy_index in range(n_decoys):
        output_name = "%s_%d" % (output_prefix, decoy_index)
        dump_sbatch_script(output_name, input_file, n_steps)

        while True:
            n_jobs_on_slurm_queue = get_number_of_jobs_in_slurm_queue(user)

            if n_jobs_on_slurm_queue <= DEFAULT.max_slurm_jobs:
                os.system("sbatch %s" % ("%s.sh" % (output_name)))
                break

def deploy_decoys_on_pyjobdistributor(input_file, output_prefix, n_decoys, n_steps, score_function):
    """
    Launch the PASSO protocol decoys in serial mode, using the PyJobDistributor.
    Each decoy will occupy one processor on the machine, and a new decoy will be
    started when a processor becomes available.
    """

    job_man = PyJobDistributor(output_prefix, n_decoys, score_function)
    while not job_man.job_complete:
        output_name = "%s_%d" % (output_prefix, job_man.current_id)
        p = single_dock_decoy(input_file, output_name, n_steps)

        # this is the last frame of the simulation. Since the PASSO protocol is
        # based on MonteCarlo, it might not correspond to the lowest energy
        # structure. It is, nonetheless, outputted in order for the
        # PyJobDistributor to finish and skip to the next decoy.
        job_man.output_decoy(p)


# ______________________________________________________________________________
#                                 M A I N
# ______________________________________________________________________________

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Perform a docking+design
        algorithm from multiple starting positions defined in a custom grid.""")
    parser.add_argument('input_file', metavar='INPUT', type=str,
        help='The input PDB file')
    parser.add_argument('-nd', '--n_decoys', metavar='', type=int,
        help='Number of decoys (Default: %d)' % (DEFAULT.n_decoys),
        default = DEFAULT.n_decoys)
    parser.add_argument('-ns', '--n_steps', metavar='', type=int,
        help='Number of PASSO protocol steps (Default: %d)' % (DEFAULT.n_steps),
        default = DEFAULT.n_steps)
    parser.add_argument('-s', '--slurm', action = 'store_true',
        help = "Use SLURM to launch parallel decoys (Default: False)")

    args = parser.parse_args()
    validate_arguments(args)
    score_function = get_fa_scorefxn()

    if args.slurm:
        deploy_decoys_on_slurm(
            args.input_file,
            args.output,
            args.n_decoys,
            args.n_steps)
    else:
        deploy_decoys_on_pyjobdistributor(
            args.input_file,
            args.output,
            args.n_decoys,
            args.n_steps,
            score_function)