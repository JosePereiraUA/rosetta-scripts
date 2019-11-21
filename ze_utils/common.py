from os.path import expanduser as expusr

#
#                            C O M M O N               
#
#           \\ INITIALLY CREATED BY JOSE PEREIRA, 2019 \\
#
# Functions included:
# 1) get_number_of_jobs_in_slurm_queue
#    Returns the numbers of SLURM jobs in the current user's queue
#
# 2) progress_bar
#    Returns a progress bar string with a set length and set to the given value
#
# 3) read_matrix_from_txt_file
#    Returns a BLOSUM62 matrix and its entries
#
# 4) blosum62
#    Use a BLOSUM62 matrix and entries to calculate the score between 2 strings

def get_number_of_jobs_in_slurm_queue(user):
    """
    Probe the SLURM framework for the number of jobs currently being performed
    on the 'user' profile. This function accounts for array jobs pending
    spawning. 
    """

    import re
    import subprocess

    # Number of jobs is, initially, the number of line in the squeue (minus the
    # header)
    squeue_call = subprocess.check_output(["squeue", "-u", user]).splitlines()
    number_of_jobs_in_slurm_queue = len(squeue_call) - 1

    # Scan each job to check if they are pending array jobs (should have the
    # number of elements in the array in the job description)
    for line in squeue_call:
        # elem1 corresponds to the job discription. Number of array jobs appear
        # between '[' and ']'.
        elem1 = re.split(']|\[|-', str(line.split()[0])[2:-1])
        if len(elem1) > 1:
            # If the job is defined as an array, add the number of array jobs
            number_of_jobs_in_slurm_queue += int(elem1[2]) - int(elem1[1])

    return number_of_jobs_in_slurm_queue


def progress_bar(current, maximum, length, fill = "#", empty = "-"):
    """
    Prints a progress bar of 'length' characters, which corresponds to the
    'maximum' value, with the 'fill' character repeated until the 'current'
    value. the remaining empty part of the bar is occupied by repeated 'empty'
    characters.
    """

    import math

    x = math.floor((current * length) / maximum)
    y = length - x

    return "[" + fill*x + empty*y + "]"


def read_matrix_from_txt_file(f = expusr("~/Desktop/scripts/static/b62.txt")):
    """
    Read a BLOSUM62 matrix from a TXT file. Since this file requires a specific
    format, the usage of the matrix on the default location is encourgaed. 
    Returns both the loaded matrix and the list of its entries.

    See also: blosum62
    """

    matrix = []
    with open(f, "r") as b62:
        entries = b62.readline().split()
        for line in b62:
            matrix.append([int(value) for value in line.split()[1:]])
    return matrix, entries


def blosum62(matrix, entries, query, reference):
    """
    Use a matrix (and its entries) to calculate the BLOSUM62 score between a
    query and reference structure. Returns the total score and the score per
    residue on a list.

    See also: read_matrix_from_txt_file
    """

    per_residue_score = []
    for q, r in zip(query, reference):
        residue_score = matrix[entries.index(q)][entries.index(r)]
        per_residue_score.append(residue_score)
        
    return sum(per_residue_score), per_residue_score