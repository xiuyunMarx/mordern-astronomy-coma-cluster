#!/usr/bin/env python3
"""Plot mean_fiber_ra vs mean_fiber_dec from desi-coma.csv."""

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


def read_ra_dec(csv_path, z_low, z_high):
    ra_values = []
    dec_values = []
    skipped = 0
    below_z_low = 0
    above_z_high = 0

    with csv_path.open(newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            ra = finite_float(row.get("mean_fiber_ra"))
            dec = finite_float(row.get("mean_fiber_dec"))
            z = finite_float(row.get("z"))
            if ra is None or dec is None or z is None:
                skipped += 1
                continue
            if z_low is not None and z < z_low:
                below_z_low += 1
                continue
            if z_high is not None and z > z_high:
                above_z_high += 1
                continue
            ra_values.append(ra)
            dec_values.append(dec)

    return ra_values, dec_values, skipped, below_z_low, above_z_high


def make_title(z_low, z_high):
    title = "DESI Coma: mean_fiber_ra vs mean_fiber_dec"
    if z_low is None and z_high is None:
        return title

    low_label = "-inf" if z_low is None else f"{z_low:g}"
    high_label = "inf" if z_high is None else f"{z_high:g}"
    return f"{title} ({low_label} <= z <= {high_label})"


def plot_scatter(ra_values, dec_values, output_path, z_low, z_high):
    fig, ax = plt.subplots(figsize=(8, 7), dpi=220)
    ax.scatter(
        ra_values,
        dec_values,
        s=1.0,
        alpha=0.35,
        linewidths=0,
        color="#1f77b4",
    )
    ax.set_xlabel("mean_fiber_ra (deg)")
    ax.set_ylabel("mean_fiber_dec (deg)")
    ax.set_title(make_title(z_low, z_high))
    ax.grid(True, alpha=0.25, linewidth=0.6)
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a scatter plot of mean_fiber_ra and mean_fiber_dec."
    )
    parser.add_argument(
        "--input",
        default="desi-coma.csv",
        type=Path,
        help="Input CSV file. Default: desi-coma.csv",
    )
    parser.add_argument(
        "--output",
        default="desi-coma-mean-fiber-ra-dec-scatter.png",
        type=Path,
        help="Output PNG file.",
    )
    parser.add_argument(
        "--z-low",
        default=None,
        type=float,
        help="Minimum z to include. Default: no lower limit.",
    )
    parser.add_argument(
        "--z-high",
        default=None,
        type=float,
        help="Maximum z to include. Default: no upper limit.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.z_low is not None and args.z_high is not None and args.z_low > args.z_high:
        raise SystemExit("--z-low must be less than or equal to --z-high.")

    ra_values, dec_values, skipped, below_z_low, above_z_high = read_ra_dec(
        args.input, args.z_low, args.z_high
    )
    if not ra_values:
        raise SystemExit("No valid mean_fiber_ra/mean_fiber_dec values found after filtering.")

    plot_scatter(ra_values, dec_values, args.output, args.z_low, args.z_high)
    print(f"Saved: {args.output.resolve()}")
    print(
        f"Plotted: {len(ra_values)} points; "
        f"skipped invalid: {skipped}; "
        f"below z-low: {below_z_low}; "
        f"above z-high: {above_z_high}"
    )
    print(f"RA range: {min(ra_values):.6f} to {max(ra_values):.6f}")
    print(f"Dec range: {min(dec_values):.6f} to {max(dec_values):.6f}")


if __name__ == "__main__":
    main()
