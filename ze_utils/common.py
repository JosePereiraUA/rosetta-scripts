# Questions: jose.manuel.pereira@ua.pt

from os.path import expanduser as expusr

#
#                              C O M M O N               
#
#            \\ INITIALLY CREATED BY JOSE PEREIRA, 2019 \\
#
#    Functions included:
#
# 1) get_number_of_jobs_in_slurm_queue
#    Returns the numbers of SLURM jobs in the current user's queue
#
# 2) verify_launch_failed_requeued_held
#    Automatically recovers jobs from the '(launch failed requeued held)' error
#
# 3) progress_bar
#    Returns a progress bar string with a set length and set to the given value
#
# 4) read_matrix_from_txt_file
#    Returns a BLOSUM62 matrix and its entries
#
# 5) blosum62
#    Use a BLOSUM62 matrix and entries to calculate the score between 2 strings
#
# 6) overwrite_dir
#    Create a new directory, overwritting any existing data
#
# 7) load_data_from_fasc_file
#    Returns an array with all Dictionary entries found in the .fasc file.


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


def verify_launch_failed_requeued_held(user):
    """
    Verifies the 'user's squeue, searching for '(launch failed requeued held)'
    errors, releasing the jobs.
    """

    import os
    import subprocess

    squeue = subprocess.check_output(["squeue", "-u", user]).splitlines()[1:]
    failed_jobs = []
    for line in squeue:
        elem = [str(e)[2:-1] for e in line.split()]
        if " ".join(elem[7:]) == '(launch failed requeued held)':
            failed_jobs.append(elem[0])

    if len(failed_jobs) > 0:
        print("Releasing held jobs: %s" % (" ".join(failed_jobs)))
        os.system("scontrol release %s" % (" ".join(failed_jobs)))


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


def overwrite_dir(path):
    """
    Create a new directory at 'path', overwritting any pre existing data.
    """

    import os
    import shutil

    if os.path.exists(path):
        print(" >> [ WARNING ] Overwritting directory %s" % (path))
        shutil.rmtree(path)
    os.mkdir(path)


def load_data_from_fasc_file(file_path):
    """
    Returns an array with all Dictionary entries found in the .fasc file.
    """
    json_path = '['
    with open(file_path, "r") as file_in:
        for line in file_in:
            json_path += line[:-1]
            json_path += ','
    json_path = json_path[:-1] + ']'
    
    try:
        return json.loads(json_path)
    except:
        import json
        return json.loads(json_path)