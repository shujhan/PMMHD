#!/usr/bin/env python3
import os
import json
import subprocess
import time

def ensure_dirs(base_dir):
    """Create all output directories."""
    base = os.path.join(base_dir, "simulation_output")
    subdirs = ["xs", "ys", "w0s", "j0s", "uweights", "u1s", "u2s", "b1s", "b2s","panels"]
    for sd in subdirs:
        path = os.path.join(base, sd)
        os.makedirs(path, exist_ok=True)

def parse_species(deck_path):
    """Read species names from deck.json."""
    with open(deck_path, "r") as f:
        deck = json.load(f)
    species = []
    if "initial_list" in deck:
        for sp in deck["initial_list"]:
            if "name" in sp:
                species.append(sp["name"])
    return species

def main():
    # Get paths
    # script_dir = os.path.dirname(os.path.abspath(__file__))
    script_dir = os.path.dirname('./')
    deck_path = os.path.join(script_dir, "deck.json")

    # Parse species and make directories
    species_names = parse_species(deck_path)
    if not species_names:
        raise ValueError("No vorticity or current density found in deck.json under 'initial_list'.")

    ensure_dirs(script_dir)

    # Run simulation
    print(f"\n[RUN] Executing: PMMHD_gpu {deck_path}")
    use_gpu = True
    if use_gpu:
        PMMHD_exe = './PMMHD_gpu'
    else:
        PMMHD_exe = './PMMHD_cpu'


    PMMHD_args = [PMMHD_exe]
    out_path = os.path.join(script_dir, "sim_out.txt")
    err_path = os.path.join(script_dir, "sim_err.txt")

    t1 = time.time()
    with open(out_path,'w') as log:
        with open(err_path,'w') as err_log:
            proc = subprocess.run(PMMHD_args,stdout=log,stderr=err_log)
    t2 = time.time()

    if proc.returncode == 0:
        print('done!')
        simulation_has_run = True
    else:
        simulation_has_run = False
        print('simulation stopped with error code %i'%proc.returncode)

    with open(out_path,'a') as log:
        log.write(f'compute time {t2-t1:.03f}s')
    print(f'simulation compute time {t2-t1:.03f}s')

if __name__ == "__main__":
    main()