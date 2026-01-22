#!/usr/bin/env python3
import os
import re
import json
import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Polygon

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
        os.path.join(base, "j0s", f"j0s_{time}"),
        os.path.join(base, "uweights", f"uweights_{time}"),
        os.path.join(base, "u1s", f"u1s_{time}"),
        os.path.join(base, "u2s", f"u2s_{time}"),
        os.path.join(base, "panels", f"leaf_point_inds_{time}")
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



# def plot_phase_space(sim_dir, simulation_dictionary, step_ii,flim, simulation_has_run = True, do_save = False):
def plot_phase_space(ax, xs, ys, fs, panels, species, step_ii, symmetric=False, clim=None):

    xs = np.asarray(xs)
    ys = np.asarray(ys)
    fs = np.asarray(fs)
    panels = np.asarray(panels)

    # simtime = step_ii * sd["dt"]
    simtime = int(step_ii) * 0.25
    
    num_panels = int(panels.size/9)
    panels = np.reshape(panels, (num_panels,9))
    panels_fs = np.zeros(4*num_panels)

    patches = []
    for ii, panel in enumerate(panels):
        panel_xs = xs[panel]
        panel_ys = ys[panel]
        panel_fs = fs[panel]

        p0 = [0,1,4,3]
        panels_fs[4*ii] = .25*sum(panel_fs[p0])
        rect_pts = np.vstack([panel_xs[p0],panel_ys[p0]]).T
        patches.append(Polygon(rect_pts))

        p1 = [1,2,5,4]
        panels_fs[4*ii+1] = .25*sum(panel_fs[p1])
        rect_pts = np.vstack([panel_xs[p1],panel_ys[p1]]).T
        patches.append(Polygon(rect_pts))

        p2 = [3,4,7,6]
        panels_fs[4*ii+2] = .25*sum(panel_fs[p2])
        rect_pts = np.vstack([panel_xs[p2],panel_ys[p2]]).T
        patches.append(Polygon(rect_pts))

        p3 = [4,5,8,7]
        panels_fs[4*ii+3] = .25*sum(panel_fs[p3])
        rect_pts = np.vstack([panel_xs[p3],panel_ys[p3]]).T
        patches.append(Polygon(rect_pts))

    # --- choose clim ---
    if clim is None:
        if symmetric:
            m = np.percentile(np.abs(panels_fs), 95)
            if m == 0:
                m = np.max(np.abs(panels_fs)) if np.max(np.abs(panels_fs)) > 0 else 1.0
            clim = (-m, m)
        else:
            # vmin = np.nanmin(panels_fs)
            # vmax = np.nanmax(panels_fs)
            vmin = np.nanmin(fs)
            vmax = np.nanmax(fs)
            if vmin == vmax:
                vmin, vmax = float(np.min(panels_fs)), float(np.max(panels_fs))
                if vmin == vmax:
                    vmin, vmax = vmin - 1.0, vmax + 1.0
            clim = (vmin, vmax)

    p = PatchCollection(patches, cmap=matplotlib.cm.jet)
#     if do_show_panels:
#         p.set_ec('black')

    p.set_array(panels_fs)
    p.set_clim(clim)
    ax.add_collection(p)
    # cb = fig.colorbar(p, ax=ax)

    ax.set_xlim(np.min(xs), np.max(xs))
    ax.set_ylim(np.min(ys), np.max(ys))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f"{species}  time: {simtime}")

    # colorbar on the same figure as ax
    fig = ax.figure
    fig.colorbar(p, ax=ax)

    return p
# end plot_phase_space







def main():
    p = argparse.ArgumentParser(description="Read PMMHD particle dumps and plot x–y colored by f.")
    p.add_argument("--base-dir", default=".", help="Project root (default: current dir)")
    p.add_argument("--deck", default="deck.json", help="Deck path (default: deck.json in base dir)")
    p.add_argument("--time", default=0,
                   help="File time after xs_/ys... If omitted, auto-pick latest per species.")
    p.add_argument("--plot_points", '-pts',action="store_true", help="plot points, no panels")
    p.add_argument("--plot_phase_space", '-ps', action="store_true")
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



    if args.plot_points:
        time = args.time
        xs_path, ys_path, w0s_path, j0s_path, uweights_path, u1s_path,u2s_path, panels_path = build_paths(base_dir, time)

        x = read_bin_doubles(xs_path)
        y = read_bin_doubles(ys_path)
        w0s = read_bin_doubles(w0s_path)
        j0s = read_bin_doubles(j0s_path)
        uweights = read_bin_doubles(uweights_path)
        u1s = read_bin_doubles(u1s_path)
        u2s = read_bin_doubles(u2s_path)
        panels = np.fromfile(panels_path,dtype='int32')

        # === Write xs and ys to text files ===
        out_x_file = os.path.join(base_dir, f"xs_{time}.txt")
        out_y_file = os.path.join(base_dir, f"ys_{time}.txt")
        np.savetxt(out_x_file, x, fmt="%.6e")
        np.savetxt(out_y_file, y, fmt="%.6e")

        out_u1s_file = os.path.join(base_dir, f"u1s_{time}.txt")
        out_u2s_file = os.path.join(base_dir, f"u2s_{time}.txt")
        out_w0s_file = os.path.join(base_dir, f"w0s_{time}.txt")
        out_j0s_file = os.path.join(base_dir, f"j0s_{time}.txt")
        np.savetxt(out_u1s_file, u1s, fmt="%.6e")
        np.savetxt(out_u2s_file, u2s, fmt="%.6e")
        np.savetxt(out_w0s_file, w0s, fmt="%.6e")
        np.savetxt(out_j0s_file, j0s, fmt="%.6e")



        fig, axes = plt.subplots(1, 2, figsize=(10,4))
        ax_left, ax_right = axes
        plot_species(ax_left, x, y, w0s, species_list[0], time)
        plot_species(ax_right, x, y, j0s, species_list[1], time)
        fig.tight_layout()
        out_name = os.path.join(base_dir, f"particles_{time}.png")
        fig.savefig(os.path.join(base_dir, f"particles_{time}.png"), dpi=200)
        print(f"[OK] Saved figure to {out_name}")



   
    if args.plot_phase_space:

        time = args.time
        xs_path, ys_path, w0s_path, j0s_path, uweights_path, u1s_path,u2s_path, panels_path = build_paths(base_dir, time)

        xs = read_bin_doubles(xs_path)
        ys = read_bin_doubles(ys_path)
        w0s = read_bin_doubles(w0s_path)
        j0s = read_bin_doubles(j0s_path)
        uweights = read_bin_doubles(uweights_path)
        u1s = read_bin_doubles(u1s_path)
        u2s = read_bin_doubles(u2s_path)
        panels = np.fromfile(panels_path,dtype='int32')

        fig, axes = plt.subplots(1, 2, figsize=(10,4))
        ax_left, ax_right = axes

        # plot_phase_space(sim_dir, simulation_dictionary, timestep, flim=flim)
        plot_phase_space(ax_left, xs, ys, w0s, panels, species_list[0], time)
        plot_phase_space(ax_right, xs, ys, j0s, panels, species_list[1], time)
        fig.tight_layout()
        out_name = os.path.join(base_dir, f"phase_space_{time}.png")
        fig.savefig(out_name, dpi=200)
        print(f"[OK] Saved figure to {out_name}")

if __name__ == "__main__":
    main()