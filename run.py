#!/usr/bin/env python3
import os
import json
import subprocess

def ensure_dirs(base_dir):
    """Create all output directories."""
    base = os.path.join(base_dir, "simulation_output")
    subdirs = ["xs", "ys", "w0s", "j0s", "weights", "u1s", "u2s", "b1s", "b2s"]
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
    exe_path = os.path.join(script_dir, "PMMHD_cpu")
    deck_path = os.path.join(script_dir, "deck.json")
    # Verify files
    if not os.path.isfile(exe_path):
        raise FileNotFoundError(f"Executable not found: {exe_path}")
    if not os.path.isfile(deck_path):
        raise FileNotFoundError(f"Deck file not found: {deck_path}")

    # Parse species and make directories
    species_names = parse_species(deck_path)
    if not species_names:
        raise ValueError("No vorticity or current density found in deck.json under 'initial_list'.")

    ensure_dirs(script_dir)

    # Run simulation
    print(f"\n[RUN] Executing: {exe_path} {deck_path}")
    use_gpu = False
    if use_gpu:
        PMMHD_exe = './PMMHD_gpu'
    else:
        PMMHD_exe = './PMMHD_cpu'
    PMMHD_args = [PMMHD_exe]
    
    result = subprocess.run(PMMHD_args,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True)

    # # Save and show logs
    log_path = os.path.join(script_dir, "run_log.txt")
    with open(log_path, "w") as f:
        f.write("=== STDOUT ===\n")
        f.write(result.stdout)
        f.write("\n=== STDERR ===\n")
        f.write(result.stderr)

    # print(f"[DONE] Simulation finished. Log saved to {log_path}")
    if result.returncode != 0:
        print(f"[WARN] PMMHD_cpu exited with code {result.returncode}")
        print(result.stderr)
    else:
        print("[OK] PMMHD_cpu completed successfully.")

if __name__ == "__main__":
    main()