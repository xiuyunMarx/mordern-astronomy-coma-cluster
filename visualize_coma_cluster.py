#!/usr/bin/env python3
"""Visualize the selected Coma cluster CSV."""

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def finite_float(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def read_cluster(csv_path):
    ra_values = []
    dec_values = []
    z_values = []
    skipped = 0

    with csv_path.open(newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            ra = finite_float(row.get("mean_fiber_ra"))
            dec = finite_float(row.get("mean_fiber_dec"))
            z = finite_float(row.get("z"))
            if ra is None or dec is None or z is None:
                skipped += 1
                continue
            ra_values.append(ra)
            dec_values.append(dec)
            z_values.append(z)

    if not ra_values:
        raise SystemExit("No valid mean_fiber_ra/mean_fiber_dec/z rows found.")

    return ra_values, dec_values, z_values, skipped


def mean(values):
    return sum(values) / len(values)


def median(values):
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return 0.5 * (ordered[mid - 1] + ordered[mid])


def plot_cluster(ra_values, dec_values, z_values, output_path, bins, point_size):
    center_ra = mean(ra_values)
    center_dec = mean(dec_values)
    mean_z = mean(z_values)
    median_z = median(z_values)

    fig, (scatter_ax, hist_ax) = plt.subplots(
        1,
        2,
        figsize=(12, 5.8),
        dpi=220,
        gridspec_kw={"width_ratios": [1.25, 1.0]},
    )

    points = scatter_ax.scatter(
        ra_values,
        dec_values,
        c=z_values,
        s=point_size,
        alpha=0.78,
        linewidths=0,
        cmap="viridis",
    )
    scatter_ax.scatter(
        [center_ra],
        [center_dec],
        s=95,
        marker="x",
        linewidths=2.0,
        color="black",
        label=f"center: RA={center_ra:.3f}, Dec={center_dec:.3f}",
    )
    scatter_ax.set_xlabel("mean_fiber_ra (deg)")
    scatter_ax.set_ylabel("mean_fiber_dec (deg)")
    scatter_ax.set_title("Selected Coma Cluster in RA-Dec")
    scatter_ax.grid(True, alpha=0.25, linewidth=0.6)
    scatter_ax.legend(loc="best", fontsize=7.5, frameon=True, framealpha=0.92)
    cbar = fig.colorbar(points, ax=scatter_ax, fraction=0.046, pad=0.04)
    cbar.set_label("z")

    hist_ax.hist(z_values, bins=bins, color="#2ca25f", alpha=0.85, edgecolor="white")
    hist_ax.axvline(
        median_z,
        color="black",
        linewidth=1.4,
        label=f"median z={median_z:.5f}",
    )
    hist_ax.set_xlabel("z")
    hist_ax.set_ylabel("Count")
    hist_ax.set_title("Redshift Distribution")
    hist_ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
    hist_ax.legend(loc="best", fontsize=7.5, frameon=True, framealpha=0.92)

    fig.suptitle(
        f"DESI Coma Cluster: n={len(ra_values)}, mean z={mean_z:.5f}",
        y=1.02,
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualize a selected Coma cluster CSV in RA/Dec and z."
    )
    parser.add_argument(
        "--input",
        default="desi-coma-gmm-center-cluster.csv",
        type=Path,
        help="Input cluster CSV file. Default: desi-coma-gmm-center-cluster.csv",
    )
    parser.add_argument(
        "--output",
        default="desi-coma-gmm-center-cluster-visualization.png",
        type=Path,
        help="Output PNG file.",
    )
    parser.add_argument(
        "--bins",
        default=50,
        type=int,
        help="Number of bins for the z histogram. Default: 50",
    )
    parser.add_argument(
        "--point-size",
        default=8.0,
        type=float,
        help="Scatter point size. Default: 8.0",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    ra_values, dec_values, z_values, skipped = read_cluster(args.input)
    plot_cluster(ra_values, dec_values, z_values, args.output, args.bins, args.point_size)

    print(f"Saved: {args.output.resolve()}")
    print(f"Plotted: {len(ra_values)} points; skipped invalid: {skipped}")
    print(f"RA range: {min(ra_values):.6f} to {max(ra_values):.6f}")
    print(f"Dec range: {min(dec_values):.6f} to {max(dec_values):.6f}")
    print(f"z range: {min(z_values):.6f} to {max(z_values):.6f}")


if __name__ == "__main__":
    main()
