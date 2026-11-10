#!/usr/bin/env python3
"""Plot radial profiles and 2D poloidal maps from GTS-S snapshot files.

Usage:
    python plot_snapshot.py SNAP.h5 [SNAP2.h5 ...] [--outdir DIR] [--instant] [--min-density FRAC] [--show]

With one file, the plots show that snapshot. With several files, every quantity
is averaged over the files (for example, ten 1-us snapshots give a 10-us average).
Temperatures are hidden (gray) where the density is below FRAC x its maximum.

Output (in --outdir, default: current directory):
    <name>_profiles.png   flux-surface-averaged n_e, n_i, T_e, T_i versus a = sqrt(psi_t)
    <name>_2D.png         Phi, n_e, T_e, T_i on one poloidal plane

Requires numpy, h5py and matplotlib.
"""
import argparse
import os

import h5py
import numpy as np

ELECTRON_COLOR = "#2a78d6"
ION_COLOR = "#eb6834"


def read_mean(files, key):
    """Average dataset `key` over the given snapshot files."""
    total = None
    for fn in files:
        with h5py.File(fn, "r") as f:
            x = np.asarray(f[key][()], dtype=float)
        total = x if total is None else total + x
    return total / len(files)


def read_first(files, key):
    with h5py.File(files[0], "r") as f:
        return np.asarray(f[key][()], dtype=float)


def times_us(files):
    t = []
    for fn in files:
        with h5py.File(fn, "r") as f:
            t.append(float(np.atleast_1d(f["time_s"][()])[0]) * 1e6)
    return t


def to_theta(gq_1d, data):
    """Drop the periodic copy at gq = 2*pi, map gq to theta in [-pi, pi), and sort."""
    if np.isclose(gq_1d[-1] - gq_1d[0], 2 * np.pi):
        gq_1d, data = gq_1d[:-1], data[:, :-1]
    theta = np.mod(gq_1d + np.pi, 2 * np.pi) - np.pi
    order = np.argsort(theta)
    return theta[order], data[:, order]


def mask_low_density(temperature, density, frac):
    """Hide T where the density is below frac * max(density): T is not meaningful there."""
    t = temperature.copy()
    t[~(density >= frac * np.nanmax(density))] = np.nan
    return t


def plot_profiles(files, label, outdir, instant, frac):
    import matplotlib.pyplot as plt

    prefix = "DG/ZON/" if instant else "DG/ZON/mean_"
    a = read_first(files, "DG/a_1D")
    ne = read_mean(files, prefix + "ne_tot_by_fh")
    ni = read_mean(files, prefix + "ni_by_ff")
    te = mask_low_density(read_mean(files, prefix + "Te_by_fh"), ne, frac)
    ti = mask_low_density(read_mean(files, prefix + "Ti_by_ff"), ni, frac)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    panels = [(axes[0], ne, ni, r"density [m$^{-3}$]"),
              (axes[1], te, ti, "temperature [eV]")]
    for ax, ye, yi, ylabel in panels:
        ax.plot(a, ye, color=ELECTRON_COLOR, lw=2, label="electrons")
        ax.plot(a, yi, color=ION_COLOR, lw=2, label="ions")
        ax.set_xlabel(r"$a=\sqrt{\psi_t}$")
        ax.set_ylabel(ylabel)
        ax.grid(color="0.9", lw=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.legend(frameon=False)
    kind = "at the snapshot step" if instant else "averaged since the previous snapshot"
    fig.suptitle(f"Flux-surface averages, {label} ({kind})")
    out = os.path.join(outdir, f"{base_name(files)}_profiles.png")
    fig.savefig(out, dpi=150)
    return out, fig


def plot_maps(files, label, outdir, frac):
    import matplotlib.pyplot as plt

    sg_a, sg_gq = read_first(files, "SG/a_1D"), read_first(files, "SG/gq_1D")
    dg_a, dg_gq = read_first(files, "DG/a_1D"), read_first(files, "DG/gq_1D")
    ne = read_mean(files, "DG/ne_tot_by_fh")
    ni = read_mean(files, "DG/ni")
    maps = [
        (r"$\Phi$ [V]", sg_a, sg_gq, read_mean(files, "SG/phi"), "magma", "potential"),
        (r"$n_e$ [m$^{-3}$]", dg_a, dg_gq, ne, "viridis", "density"),
        (r"$T_e$ [eV]", dg_a, dg_gq, mask_low_density(read_mean(files, "DG/Te"), ne, frac), "viridis", "temperature"),
        (r"$T_i$ [eV]", dg_a, dg_gq, mask_low_density(read_mean(files, "DG/Ti"), ni, frac), "viridis", "temperature"),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(16, 4), sharey=True, constrained_layout=True)
    for ax, (title, a, gq, data, cmap_name, kind) in zip(axes, maps):
        theta, z = to_theta(gq, data)
        ax.set_title(title)
        ax.set_xlim(a[0], a[-1])
        ax.set_ylim(-np.pi, np.pi)
        if kind == "potential" and np.nanmax(np.abs(z)) == 0.0:
            ax.text(0.5, 0.5, "Φ = 0\n(no Poisson solve)", transform=ax.transAxes,
                    ha="center", va="center", color="0.4")
        else:
            cmap = plt.get_cmap(cmap_name).copy()
            cmap.set_bad("0.85")  # gray: masked (density too low for a temperature)
            vmin, vmax = np.nanmin(z), np.nanmax(z)
            if kind == "density":
                vmin = 0.0  # small negative values near the wall are noise
            elif kind == "temperature":
                vmin, vmax = np.nanpercentile(z, [1, 99])
            mesh = ax.pcolormesh(a, theta, np.ma.masked_invalid(z).T, cmap=cmap,
                                 vmin=vmin, vmax=vmax, shading="nearest")
            fig.colorbar(mesh, ax=ax)
        ax.set_xlabel(r"$a=\sqrt{\psi_t}$")
        ax.set_yticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
        ax.set_yticklabels([r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    axes[0].set_ylabel(r"$\theta$")
    fig.suptitle(f"One poloidal plane, {label}")
    out = os.path.join(outdir, f"{base_name(files)}_2D.png")
    fig.savefig(out, dpi=150)
    return out, fig


def base_name(files):
    first = os.path.splitext(os.path.basename(files[0]))[0]
    if len(files) == 1:
        return first
    return f"{first}-{os.path.splitext(os.path.basename(files[-1]))[0]}"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="+", help="snapshot files, e.g. snapshots_hdf5/snap00010.h5")
    p.add_argument("--outdir", default=".", help="directory for the PNG files (default: .)")
    p.add_argument("--instant", action="store_true",
                   help="radial profiles at the snapshot step instead of the mean_* interval averages")
    p.add_argument("--min-density", type=float, default=0.05, metavar="FRAC",
                   help="hide temperatures where the density is below FRAC x its maximum (default: 0.05)")
    p.add_argument("--show", action="store_true", help="also open the figures in a window")
    args = p.parse_args()

    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    files = sorted(args.files)
    t = times_us(files)
    label = f"t = {t[0]:.1f} us" if len(t) == 1 else f"average of {len(t)} snapshots, t = {t[0]:.1f}-{t[-1]:.1f} us"
    os.makedirs(args.outdir, exist_ok=True)
    for out, _ in (plot_profiles(files, label, args.outdir, args.instant, args.min_density),
                   plot_maps(files, label, args.outdir, args.min_density)):
        print("wrote", out)
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
