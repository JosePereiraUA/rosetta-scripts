# TODO: Documentation

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


def progress_bar(current, maximum, length):

    import math

    x = math.floor((current * length) / maximum)
    y = length - x
    return "[" + "#"*x + "-"*y + "]"
