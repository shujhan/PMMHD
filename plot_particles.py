#!/usr/bin/env python3
import os
import re
import json
import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt

def load_deck_species(deck_path):
    with open(deck_path, "r") as f:
        d = json.load(f)
    return [sp["name"] for sp in d.get("initial_list", []) if isinstance(sp, dict) and "name" in sp]

def read_bin_doubles(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Missing file: {path}")
    arr = np.fromfile(path, dtype=np.float64)
    return arr

def pick_latest_tag(species_dir):
    """
    Look under simulation_output/<species>/{xs,ys,fs}/ and pick a tag that exists in all three.
    We find the newest xs_* file by mtime and use its suffix as the tag.
    """
    xs_dir = os.path.join(species_dir, "xs")
    ys_dir = os.path.join(species_dir, "ys")
    fs_dir = os.path.join(species_dir, "fs")
    xs_files = sorted(glob.glob(os.path.join(xs_dir, "xs_*")), key=os.path.getmtime)
    if not xs_files:
        raise FileNotFoundError(f"No xs_* files found in {xs_dir}")
    # Try newest first and ensure corresponding ys_/fs_ exist
    for candidate in reversed(xs_files):
        tag = os.path.basename(candidate)[len("xs_"):]  # suffix after "xs_"
        if (os.path.isfile(os.path.join(ys_dir, "ys_" + tag)) and
            os.path.isfile(os.path.join(fs_dir, "fs_" + tag))):
            return tag
    # Fallback: just use newest xs tag and hope others exist
    return os.path.basename(xs_files[-1])[len("xs_"):]

def build_paths(base_dir, time):
    base = os.path.join(base_dir, "simulation_output")
    return (
        os.path.join(base, "xs", f"xs_{time}"),
        os.path.join(base, "ys", f"ys_{time}"),
        os.path.join(base, "w0s", f"w0s_{time}"),
        os.path.join(base, "j0s", f"j0s_{time}")
    )

def plot_species(ax, x, y, f, species, time):
    n = min(len(x), len(y), len(f))
    if n == 0:
        raise ValueError(f"No particles to plot for {species} (time {time}).")
    x = x[:n]; y = y[:n]; f = f[:n]
    sc = ax.scatter(x, y, c=f, s=3, alpha=0.8)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(f"{species}   [{n} particles]   time: {time}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    cbar = plt.colorbar(sc, ax=ax)
    # cbar.set_label("f")

def main():
    p = argparse.ArgumentParser(description="Read PMMHD particle dumps and plot x–y colored by f.")
    p.add_argument("--base-dir", default=".", help="Project root (default: current dir)")
    p.add_argument("--deck", default="deck.json", help="Deck path (default: deck.json in base dir)")
    p.add_argument("--time", default=0,
                   help="File time after xs_/ys... If omitted, auto-pick latest per species.")
    p.add_argument("--out", default="particles.png",
                   help="Output image filename (default: particles.png)")
    args = p.parse_args()

    base_dir = os.path.abspath(args.base_dir)
    deck_path = os.path.join(base_dir, args.deck)
    if not os.path.isfile(deck_path):
        raise FileNotFoundError(f"Deck not found: {deck_path}")

    species_list = load_deck_species(deck_path)
    if not species_list:
        raise ValueError("No species found in deck.json: initial_list is empty.")



    fig, axes = plt.subplots(1, 2, figsize=(10,4))
    ax_left, ax_right = axes
    time = args.time
    xs_path, ys_path, w0s_path, j0s_path = build_paths(base_dir, time)

    x = read_bin_doubles(xs_path)
    y = read_bin_doubles(ys_path)
    w0 = read_bin_doubles(w0s_path)
    j0 = read_bin_doubles(j0s_path)

    # === Write xs and ys to text files ===
    # out_x_file = os.path.join(base_dir, f"xs_{time}.txt")
    # out_y_file = os.path.join(base_dir, f"ys_{time}.txt")
    # np.savetxt(out_x_file, x, fmt="%.6e")
    # np.savetxt(out_y_file, y, fmt="%.6e")

    plot_species(ax_left, x, y, w0, species_list[0], time)
    plot_species(ax_right, x, y, j0, species_list[1], time)

    fig.tight_layout()
    fig.savefig(os.path.join(base_dir, args.out), dpi=200)
    print(f"[OK] Saved figure to {os.path.join(base_dir, args.out)}")

if __name__ == "__main__":
    main()