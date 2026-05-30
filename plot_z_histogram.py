#!/usr/bin/env python3
"""Plot a histogram of z from desi-coma.csv."""

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


def read_z(csv_path, max_z):
    z_values = []
    skipped = 0
    excluded_above_max = 0

    with csv_path.open(newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            z = finite_float(row.get("z"))
            if z is None:
                skipped += 1
                continue
            if max_z is not None and z > max_z:
                excluded_above_max += 1
                continue
            z_values.append(z)

    return z_values, skipped, excluded_above_max


def plot_histogram(z_values, output_path, bins):
    fig, ax = plt.subplots(figsize=(8, 5.5), dpi=220)
    ax.hist(z_values, bins=bins, color="#2ca25f", alpha=0.85, edgecolor="white")
    ax.set_xlabel("z")
    ax.set_ylabel("Count")
    ax.set_title("DESI Coma: histogram of z (z <= 0.08)")
    ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(description="Create a histogram of z.")
    parser.add_argument(
        "--input",
        default="desi-coma.csv",
        type=Path,
        help="Input CSV file. Default: desi-coma.csv",
    )
    parser.add_argument(
        "--output",
        default="desi-coma-z-histogram.png",
        type=Path,
        help="Output PNG file.",
    )
    parser.add_argument(
        "--bins",
        default=80,
        type=int,
        help="Number of histogram bins. Default: 80",
    )
    parser.add_argument(
        "--max-z",
        default=0.08,
        type=float,
        help="Keep only z values no greater than this value. Default: 0.08",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    z_values, skipped, excluded_above_max = read_z(args.input, args.max_z)
    if not z_values:
        raise SystemExit("No valid z values found.")

    plot_histogram(z_values, args.output, args.bins)
    print(f"Saved: {args.output.resolve()}")
    print(
        f"Plotted: {len(z_values)} z values; "
        f"skipped invalid: {skipped}; "
        f"excluded above max-z: {excluded_above_max}"
    )
    print(f"z range: {min(z_values):.6f} to {max(z_values):.6f}")


if __name__ == "__main__":
    main()
