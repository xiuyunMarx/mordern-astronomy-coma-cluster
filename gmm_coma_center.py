#!/usr/bin/env python3
"""Run GMM/EM clustering with a user-provided initial Coma center."""

import argparse
import csv
import math
import warnings
from pathlib import Path

import numpy as np
from sklearn.cluster import KMeans
from sklearn.exceptions import ConvergenceWarning
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler


FEATURE_COLUMNS = ("mean_fiber_ra", "mean_fiber_dec", "z")


def finite_float(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def read_points(csv_path):
    rows = []
    features = []
    skipped = 0

    with csv_path.open(newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        if reader.fieldnames is None:
            raise SystemExit("CSV has no header.")
        missing = [col for col in FEATURE_COLUMNS if col not in reader.fieldnames]
        if missing:
            raise SystemExit(f"Missing required columns: {missing}")

        for row in reader:
            point = [finite_float(row.get(col)) for col in FEATURE_COLUMNS]
            if any(value is None for value in point):
                skipped += 1
                continue
            rows.append(row)
            features.append(point)

    if not features:
        raise SystemExit("No valid rows found for clustering.")

    return reader.fieldnames, rows, np.asarray(features, dtype=float), skipped


def make_initial_means(x_scaled, initial_center_scaled, n_components, random_state):
    if n_components == 1:
        return np.asarray([initial_center_scaled], dtype=float)

    # Use KMeans to initialize the non-target components, but force one component
    # to start exactly at the user-supplied Coma center.
    kmeans = KMeans(n_clusters=n_components, random_state=random_state, n_init=20)
    kmeans.fit(x_scaled)
    kmeans_centers = kmeans.cluster_centers_
    nearest_center = int(
        np.argmin(np.linalg.norm(kmeans_centers - initial_center_scaled, axis=1))
    )

    means = [initial_center_scaled]
    means.extend(
        kmeans_centers[index]
        for index in range(n_components)
        if index != nearest_center
    )
    return np.asarray(means, dtype=float)


def original_center(scaler, center_scaled):
    return scaler.inverse_transform(np.asarray(center_scaled, dtype=float).reshape(1, -1))[0]


def fit_gmm_with_center_trace(x_scaled, scaler, initial_center, args):
    initial_center_scaled = scaler.transform(np.asarray(initial_center).reshape(1, -1))[0]
    means_init = make_initial_means(
        x_scaled, initial_center_scaled, args.components, args.random_state
    )

    gmm = GaussianMixture(
        n_components=args.components,
        covariance_type=args.covariance_type,
        random_state=args.random_state,
        reg_covar=args.reg_covar,
        max_iter=1,
        n_init=1,
        tol=0.0,
        warm_start=True,
        means_init=means_init,
    )

    trace = []
    target_label = 0
    previous_center_scaled = initial_center_scaled

    for iteration in range(1, args.max_iter + 1):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            gmm.fit(x_scaled)

        distances = np.linalg.norm(gmm.means_ - previous_center_scaled, axis=1)
        target_label = int(np.argmin(distances))
        center_scaled = gmm.means_[target_label].copy()
        center = original_center(scaler, center_scaled)
        center_shift_scaled = float(np.linalg.norm(center_scaled - previous_center_scaled))
        labels = gmm.predict(x_scaled)
        probabilities = gmm.predict_proba(x_scaled)
        target_mask = labels == target_label
        mean_probability = float(probabilities[target_mask, target_label].mean())

        trace.append(
            {
                "iteration": iteration,
                "target_label": target_label,
                "center_ra": center[0],
                "center_dec": center[1],
                "center_z": center[2],
                "center_shift_scaled": center_shift_scaled,
                "selected_count": int(target_mask.sum()),
                "mean_probability": mean_probability,
                "lower_bound": float(gmm.lower_bound_),
            }
        )

        previous_center_scaled = center_scaled
        if center_shift_scaled <= args.center_tol:
            break

    return gmm, target_label, trace


def write_selected_rows(output_path, fieldnames, rows, labels, probabilities, target_label, min_prob):
    selected_count = 0
    with output_path.open("w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for row, label, probability in zip(rows, labels, probabilities[:, target_label]):
            if label == target_label and probability >= min_prob:
                writer.writerow(row)
                selected_count += 1
    return selected_count


def write_trace(trace_path, trace):
    fieldnames = (
        "iteration",
        "target_label",
        "center_ra",
        "center_dec",
        "center_z",
        "center_shift_scaled",
        "selected_count",
        "mean_probability",
        "lower_bound",
    )
    with trace_path.open("w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(trace)


def component_summary(x, labels, probabilities, scaler, gmm):
    summaries = []
    centers = scaler.inverse_transform(gmm.means_)
    for label in range(gmm.n_components):
        mask = labels == label
        if not mask.any():
            continue
        component = x[mask]
        summaries.append(
            {
                "label": label,
                "count": int(mask.sum()),
                "center_ra": float(centers[label, 0]),
                "center_dec": float(centers[label, 1]),
                "center_z": float(centers[label, 2]),
                "hard_mean_ra": float(component[:, 0].mean()),
                "hard_mean_dec": float(component[:, 1].mean()),
                "hard_mean_z": float(component[:, 2].mean()),
                "median_z": float(np.median(component[:, 2])),
                "mean_probability": float(probabilities[mask, label].mean()),
            }
        )
    return summaries


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Cluster desi-coma.csv with GMM/EM using a user-provided initial "
            "Coma center. The tracked Gaussian center is optimized during EM."
        )
    )
    parser.add_argument("--input", default="desi-coma.csv", type=Path)
    parser.add_argument("--output", default="desi-coma-gmm-center-cluster.csv", type=Path)
    parser.add_argument(
        "--trace-output",
        default="desi-coma-gmm-center-trace.csv",
        type=Path,
        help="CSV recording the optimized Coma center after each EM iteration.",
    )
    parser.add_argument("--center-ra", required=True, type=float)
    parser.add_argument("--center-dec", required=True, type=float)
    parser.add_argument("--center-z", required=True, type=float)
    parser.add_argument("--components", default=4, type=int)
    parser.add_argument("--max-iter", default=200, type=int)
    parser.add_argument(
        "--center-tol",
        default=1e-4,
        type=float,
        help="Stop when the tracked center shift is below this value in standardized units.",
    )
    parser.add_argument(
        "--min-prob",
        default=0.0,
        type=float,
        help="Minimum posterior probability required when writing selected rows.",
    )
    parser.add_argument(
        "--covariance-type",
        default="full",
        choices=("full", "tied", "diag", "spherical"),
    )
    parser.add_argument("--reg-covar", default=1e-6, type=float)
    parser.add_argument("--random-state", default=42, type=int)
    args = parser.parse_args()

    if args.components < 1:
        raise SystemExit("--components must be at least 1.")
    if args.max_iter < 1:
        raise SystemExit("--max-iter must be at least 1.")
    if not 0.0 <= args.min_prob <= 1.0:
        raise SystemExit("--min-prob must be between 0 and 1.")

    return args


def main():
    args = parse_args()
    fieldnames, rows, x, skipped = read_points(args.input)

    scaler = StandardScaler()
    x_scaled = scaler.fit_transform(x)
    initial_center = [args.center_ra, args.center_dec, args.center_z]

    gmm, target_label, trace = fit_gmm_with_center_trace(
        x_scaled, scaler, initial_center, args
    )
    labels = gmm.predict(x_scaled)
    probabilities = gmm.predict_proba(x_scaled)
    selected_count = write_selected_rows(
        args.output,
        fieldnames,
        rows,
        labels,
        probabilities,
        target_label,
        args.min_prob,
    )
    write_trace(args.trace_output, trace)

    final_center = original_center(scaler, gmm.means_[target_label])
    print(f"Input: {args.input.resolve()}")
    print(f"Output: {args.output.resolve()}")
    print(f"Trace output: {args.trace_output.resolve()}")
    print(f"Valid rows clustered: {len(rows)}")
    print(f"Skipped invalid rows: {skipped}")
    print(
        "Initial center: "
        f"RA={args.center_ra:.6f}, Dec={args.center_dec:.6f}, z={args.center_z:.6f}"
    )
    print(
        "Optimized center: "
        f"RA={final_center[0]:.6f}, Dec={final_center[1]:.6f}, z={final_center[2]:.6f}"
    )
    print(f"Tracked label: {target_label}")
    print(f"EM iterations: {trace[-1]['iteration']}")
    print(f"Selected rows written: {selected_count}")
    print("Components:")
    for item in component_summary(x, labels, probabilities, scaler, gmm):
        marker = "  <== tracked Coma component" if item["label"] == target_label else ""
        print(
            f"  label={item['label']} count={item['count']} "
            f"center_ra={item['center_ra']:.6f} "
            f"center_dec={item['center_dec']:.6f} "
            f"center_z={item['center_z']:.6f} "
            f"hard_mean_z={item['hard_mean_z']:.6f} "
            f"median_z={item['median_z']:.6f} "
            f"mean_prob={item['mean_probability']:.4f}{marker}"
        )


if __name__ == "__main__":
    main()
