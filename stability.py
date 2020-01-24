# Questions: jose.manuel.pereira@ua.pt

import shutil
import argparse
from pyrosetta import *
from align_pyrosetta import get_mask_from_c_alphas, align
from relax import deploy_decoys_on_slurm, deploy_decoys_on_pyjobdistributor
from get_best import load_data_from_fasc_file, extract_n_decoys_by_parameter
from ze_utils.common import overwrite_dir


#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\


#                       Stability Assessment Script:
# ______________________________________________________________________________
#  The objective of this script is to calculate the RMSD variance between a
# given 'input_file' and itself after relaxation without constraints. If any
# clashes or misspositioning of residues is present, the RMSD value will be
# large.


class DEFAULT:
    """
    Define defaults for the script. Values can be modified using arguments.
    """
    n_decoys  = 250
    partition = "main"


def validate_arguments(args):
    """
    Validate the script arguments to prevent known input errors.
    """
    if not args.input_file[-4:] == ".pdb":
        exit("ERROR: Input file should be in PDB format")
    if args.n_decoys < 0:
        exit("ERROR: Number of decoys must be a non-negative value")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""The objective of this script
    is to calculate the RMSD variance between a given 'input_file' and itself
    after relaxation without constraints. the relaxation will be performed on a
    given number of 'n_decoys'. If '-s' flag is present, decoys will be deployed
    using the inteligent slurm approach.""")
    parser.add_argument('input_file', metavar='INPUT', type=str,
        help='The input PDB file')
    parser.add_argument('-nd', '--n_decoys', metavar='', type=int,
        help='Number of decoys (Default: %d)' % (DEFAULT.n_decoys),
        default = DEFAULT.n_decoys)
    parser.add_argument('-s', '--slurm', action = 'store_true',
        help = "Use SLURM to launch parallel decoys (Default: False)")
    parser.add_argument('-p', '--partition', metavar='', type=str,
        help='Partition for SLURM (Default: %s)' % (DEFAULT.partition),
        default = DEFAULT.partition)

    args = parser.parse_args()
    validate_arguments(args)
    init()

    # 1) Deploy the decoys without constraints
    if args.slurm:
        deploy_decoys_on_slurm(
            args.input_file,
            args.input_file[:-4] + "_relaxed",
            args.n_decoys,
            True,
            1.0,
            args.input_file[:-4] + "_relaxed.fasc",
            args.partition)
    else:
        deploy_decoys_on_pyjobdistributor(
            args.input_file,
            args.input_file[:-4] + "_relaxed",
            args.n_decoys,
            True)

    # 2) Get best structure based on total score
    data = load_data_from_fasc_file(args.input_file[:-4] + "_relaxed.fasc")
    best = extract_n_decoys_by_parameter(data, "total_score", 1)[0][0]

    # 3) Measure RMSD between two structures
    pose = pose_from_pdb(args.input_file)
    mask = get_mask_from_c_alphas(pose)

    reference_pose = pose_from_pdb(best)
    reference_mask = get_mask_from_c_alphas(reference_pose)

    print("All atoms:")
    all_atoms_rmsd = align(pose, reference_pose)
    print("C alphas only:")
    c_alpha_rmsd   = align(pose, reference_pose, mask, reference_mask)

    # 4) Store all results in a folder
    overwrite_dir("relaxed_structures")
    s = [f for f in os.listdir('.') if (os.path.isfile(f) and f[-4:] == ".pdb")]
    for filename in s:
        shutil.move(filename, "relaxed_structures")
    shutil.move(args.input_file[:-4] + "_relaxed.fasc", "relaxed_structures")

    # 5) Save results to a file
    with open("stability_assessment.dat", "w") as file_out:
        file_out.write("all_atoms %12.3f\n" % (all_atoms_rmsd))
        file_out.write("c_alpha   %12.3f" % (  c_alpha_rmsd))