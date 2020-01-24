# Questions: jose.manuel.pereira@ua.pt

import os
import re
import collections
import argparse
from ze_utils.common import get_number_of_jobs_in_slurm_queue, progress_bar

#           \\ SCRIPT INITIALLY CREATED BY JOSE PEREIRA, 2019 \\

#                        Multi Dock Status Analyzer:
# ______________________________________________________________________________
#  Prints relevant information by scanning a directory created by the
# multi-dock.py script.

class DEFAULT:
    """
    Define defaults for the script. Values can be modified using arguments.
    """
    n_best         = 5
    sort_by        = "score"
    no_comp        = False
    reverse        = False
    max_slurm_jobs = 500


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = """Analyze the evolution of a
    multi_dock process. Shows the completion status of each docking start
    position (showcasing the slowest decoy) and the best 'n' structures found so
    far.""")

    parser.add_argument('-n', '--n_best', metavar = '', type = int,
        help = 'Find the n best structures (Default: %d)' % (DEFAULT.n_best),
        default = DEFAULT.n_best)
    parser.add_argument('-s', '--sort_by', metavar = '', type = str,
        help = 'Sort by score/interface (Default: %s)' % (DEFAULT.sort_by),
        default = DEFAULT.sort_by)
    parser.add_argument('-nc', '--no_comp', action = 'store_true',
        help = 'Skip completion check (Default: %s)' % (DEFAULT.no_comp),
        default = DEFAULT.no_comp)
    parser.add_argument('-r', '--reverse', action = 'store_true',
        help = 'Show the structures with HIGHEST energy (Default: %s)' % \
            (DEFAULT.reverse), default = DEFAULT.reverse)
    
    args = parser.parse_args()

    # If in a SLURM environment, check the number of jobs running and print
    # the current resources usage percentage.
    try:
        n_slurm_jobs = get_number_of_jobs_in_slurm_queue(os.environ["USER"])
        print("\n%sSLURM resources occupied" % (" "*24))
        print(" %s %6.2f%%" % \
            (progress_bar(n_slurm_jobs, DEFAULT.max_slurm_jobs, 55),
            (n_slurm_jobs / DEFAULT.max_slurm_jobs) * 100))
    except:
        None

    energies = []
    if not 'init.conf' in os.listdir(os.getcwd()):
        exit("\n > No init.conf file found in the current working directoy." + \
        "\n > This file is created when running the multi_dock.py script.")
    else:
        # init.conf is a simple file with 3 ints: the number of docks, decoys
        # and steps in each decoy, respectively. this information allows for a
        # more efficient extraction of data from the simulation directories.
        with open("init.conf", "r") as init_conf:
            data        = [int(value) for value in init_conf.readline().split()]
            docks       = list(range(data[0]))
            decoys      = list(range(data[1]))
            max_steps   = data[2]
            total_steps = data[2] * data[1] # Total # of steps among all decoys
            init_energy = float(init_conf.readline())

        if not args.no_comp:
            # Print COMPLETION STATUS header
            print("\n%sCompletion Status" % (" "*24))
            print("\n%9s | %-9s | %-9s | %-9s %s" % \
                ("Dock #", "Finished", "Launched", "Total", "% Completion"))
            print(" %s" % ("-"*65))

        for dock in docks:
            if not args.no_comp:
                launched = len(decoys) * max_steps # How many steps launched
                status   = 0                       # How many steps finished
            dock_name = "dock_%d" % (dock)
            if dock_name in os.listdir(os.getcwd()):
                original_dir = os.getcwd()
                os.chdir(dock_name)
            else:
                if not args.no_comp:
                    # If the dock is expected (from the init.conf file) but
                    # isn't found, print as a simulation pending entry
                    print("%9s | %-9d | %-9s | %-9s %s %6.2f%% %s" % \
                        (dock, 0, 0, total_steps, 
                        progress_bar(0, total_steps, 10), 0.0,
                        "| Simulation pending"))
                continue

            for decoy in decoys:
                if not args.no_comp:
                    # Extract current status in the simulation. Stauts files are
                    # simple files with two numbers: number of steps completed
                    # and number of max steps expected, respectively.
                    status_file = "%d_%d_status.txt" % (dock, decoy)
                    try:
                        with open(status_file, "r") as s_file:
                            status += int(s_file.readline().split()[0])
                    except:
                        # If the file does not exist, it's probably because this
                        # decoy is yet to be launched
                        launched -= max_steps
    
                # Extract energy findings of the simulation. Only certain energy
                # information is printed to the energy file. Check the
                # simulation scripts (ex: PASSO class in
                # ze_utils.pyrosetta_classes) for more information
                energy_file = "%d_%d_energy.dat" % (dock, decoy)
                try:
                    with open(energy_file, "r") as e_file:
                        e_file.readline() # Skip the header
                        for line in e_file:
                            energy                = {}
                            elem                  = line.split()
                            energy["step"]        = int(elem[0])
                            energy["score"]       = float(elem[1])
                            energy["interface"]   = float(elem[2])
                            energy["dock"]        = dock
                            energy["decoy"]       = decoy
                            energies.append(energy)
                except:
                    None
            
            # Print the status information. If the dock exists but no steps were
            # completed, print as simulation pending. This information is
            # printed one dock at a time
            if not args.no_comp:
                if status == 0:
                    print("%9s | %-9d | %-9s | %-9s %s %6.2f%% %s" % \
                        (dock, 0, launched, total_steps,
                        progress_bar(status, total_steps, 10), 0.0,
                        "| Simulation pending"))
                else:
                    print("%9s | %-9d | %-9s | %-9s %s %6.2f%%" % \
                        (dock, status, launched, total_steps, 
                        progress_bar(status, total_steps, 10),
                        (status / total_steps) * 100))

            os.chdir(original_dir)

        # When information regarding the whole simulation is gathered, print the
        # energy results. Print the energy results header.
        print("\n%sBest %d structures found" % (" "*24, args.n_best))
        print("\n Initial conformation total score: %16.3f" % (init_energy))
        print("\n%10s | %10s | %10s | %12s | %12s" % \
            ("Dock #", "Decoy #", "Step #", "Total", "Interface"))
        print(" %s" % ("-"*65))

        # Sort the energies based on the defined args.sort_by. This argument
        # should be either "score" or "interface"
        energies = sorted(energies, key=lambda entry: entry[args.sort_by],
            reverse = args.reverse)
        # Print only the args.n_best entries.
        for entry in energies[0:args.n_best]:
            print("%10d | %10d | %10d | %12.3f | %12.3f" % \
                (entry["dock"],
                entry["decoy"],
                entry["step"],
                entry["score"],
                entry["interface"]))