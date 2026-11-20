"""Run T1-T7 of the email-Eu-core network-motif project.

Examples:
  python main.py --rolls 25i-0149 i240892 25i-1032
  python main.py --R 2073
  python main.py --rolls 25i-0149 i240892 25i-1032 --check-networkx
"""
import argparse
import csv
import gzip
import os
import random
import shutil
import urllib.request
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from motifs import (BIT_EDGES, CANONICAL_CODES, assert_degrees_preserved,
                    randomize, triad_census, z_scores)

BASE_URL = "https://snap.stanford.edu/data/"
FILES = ["email-Eu-core.txt", "email-Eu-core-department-labels.txt"]


def ensure_data(data_dir):
    """Download and extract the two official SNAP files if not present."""
    os.makedirs(data_dir, exist_ok=True)
    for filename in FILES:
        plain_path = os.path.join(data_dir, filename)
        if os.path.exists(plain_path):
            continue
        gz_path = plain_path + ".gz"
        if not os.path.exists(gz_path):
            print(f"Downloading {filename}.gz")
            urllib.request.urlretrieve(BASE_URL + filename + ".gz", gz_path)
        with gzip.open(gz_path, "rb") as source, open(plain_path, "wb") as target:
            shutil.copyfileobj(source, target)


def load_clean_graph(data_dir):
    """Load edges and labels; remove self-loops and duplicate directed edges."""
    edges = []
    seen = set()
    with open(os.path.join(data_dir, FILES[0]), encoding="utf-8") as file:
        for line in file:
            parts = line.split()
            if len(parts) < 2:
                continue
            u, v = map(int, parts[:2])
            if u == v or (u, v) in seen:
                continue
            seen.add((u, v))
            edges.append((u, v))

    departments = {}
    with open(os.path.join(data_dir, FILES[1]), encoding="utf-8") as file:
        for line in file:
            parts = line.split()
            if len(parts) >= 2:
                departments[int(parts[0])] = int(parts[1])
    return edges, departments


def top20_departments(edges, departments):
    """Rank by internal directed-edge count, descending; tie -> smaller ID."""
    internal = Counter()
    for u, v in edges:
        if u in departments and v in departments and departments[u] == departments[v]:
            internal[departments[u]] += 1
    ranking = sorted(set(departments.values()), key=lambda d: (-internal[d], d))
    if len(ranking) < 20:
        raise ValueError("Dataset must contain at least 20 departments")
    return ranking[:20], internal


def parameters_from_R(R, top20):
    """Calculate all group-specific settings using lexicographic pair order."""
    pairs = [(i, j) for i in range(20) for j in range(i + 1, 20)]
    pair_number = R % 190
    i, j = pairs[pair_number]
    return {
        "R": R,
        "pair_number": pair_number,
        "pair": (i, j),
        "dept_A": top20[i],
        "dept_B": top20[j],
        "n_random": 20 + (R % 31),
        "swap_factor": 5 + (R % 6),
        "seed": R,
    }


def R_from_rolls(rolls):
    """Sum the final four digits of each supplied roll number."""
    import re
    total = 0
    for roll in rolls:
        match = re.search(r"(\d{4})\s*$", roll)
        if not match:
            raise ValueError(f"Cannot read final four digits from {roll!r}")
        total += int(match.group(1))
    return total


def induced_subgraph(edges, departments, department_id):
    """Return nodes and only edges with both endpoints in this department."""
    nodes = sorted(node for node, dept in departments.items() if dept == department_id)
    node_set = set(nodes)
    sub_edges = [(u, v) for u, v in edges if u in node_set and v in node_set]
    return nodes, sub_edges


def analyze_department(sub_edges, params, label):
    """Compute real census, randomized censuses, and Z-scores for one dept."""
    real_counts = triad_census(sub_edges)
    attempts = params["swap_factor"] * len(sub_edges)
    master_rng = random.Random(params["seed"])
    null_counts = []

    for _ in range(params["n_random"]):
        network_seed = master_rng.randrange(2 ** 32)
        randomized_edges = randomize(sub_edges, attempts, network_seed)
        assert_degrees_preserved(sub_edges, randomized_edges)
        null_counts.append(triad_census(randomized_edges))

    z, mean, std = z_scores(real_counts, null_counts)
    print(f"[{label}] {params['n_random']} randomized networks; "
          f"{attempts} swap attempts per network; degree checks passed")
    return real_counts, z, mean, std


def cross_check_networkx(edges, our_counts):
    """Optional verification only; NetworkX is not used for core computation."""
    import networkx as nx
    graph = nx.DiGraph()
    graph.add_edges_from(edges)
    nx_counts = nx.triadic_census(graph)
    all_match = True
    for class_number, code in enumerate(CANONICAL_CODES, start=1):
        representative = nx.DiGraph()
        representative.add_nodes_from(range(3))
        representative.add_edges_from(
            BIT_EDGES[bit] for bit in range(6) if (code >> bit) & 1
        )
        nx_name = nx.triad_type(representative)
        matches = nx_counts[nx_name] == our_counts[class_number]
        all_match = all_match and matches
        print(f"  Class {class_number:2}: {nx_name}: "
              f"ours={our_counts[class_number]}, NetworkX={nx_counts[nx_name]} "
              f"{'OK' if matches else 'MISMATCH'}")
    if not all_match:
        raise AssertionError("NetworkX cross-check found a mismatch")
    print("NetworkX cross-check: all 13 classes match")


def write_csv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Email network triad motifs")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--rolls", nargs="+", help="roll numbers of all group members")
    group.add_argument("--R", type=int, help="provide the group R value directly")
    parser.add_argument("--data-dir", default="data")
    parser.add_argument("--out-dir", default="results")
    parser.add_argument("--check-networkx", action="store_true")
    args = parser.parse_args()

    R = args.R if args.R is not None else R_from_rolls(args.rolls)
    if R < 0:
        parser.error("R must be non-negative")
    os.makedirs(args.out_dir, exist_ok=True)

    # T1: data loading, cleaning, top-20 ranking, and group parameters.
    ensure_data(args.data_dir)
    edges, departments = load_clean_graph(args.data_dir)
    top20, internal = top20_departments(edges, departments)
    params = parameters_from_R(R, top20)

    parameter_lines = [
        f"R = {R}",
        f"Pair number = {params['pair_number']}",
        f"Pair indices = {params['pair']}",
        f"Department A = {params['dept_A']} ({internal[params['dept_A']]} internal edges)",
        f"Department B = {params['dept_B']} ({internal[params['dept_B']]} internal edges)",
        f"Random networks = {params['n_random']}",
        f"Swap multiplier = {params['swap_factor']}",
        f"Random seed = {params['seed']}",
        "",
        "Top-20 ranking (index, department ID, internal edges):",
    ]
    parameter_lines += [f"{i}: {dept_id}: {internal[dept_id]}"
                        for i, dept_id in enumerate(top20)]
    print("\n".join(parameter_lines))
    with open(os.path.join(args.out_dir, "T1_parameters.txt"), "w", encoding="utf-8") as file:
        file.write("\n".join(parameter_lines) + "\n")

    # T2: induced subgraphs.
    subgraphs = {}
    for label, department_id in (("A", params["dept_A"]), ("B", params["dept_B"])):
        nodes, sub_edges = induced_subgraph(edges, departments, department_id)
        subgraphs[label] = (department_id, nodes, sub_edges)
        print(f"Department {label} (ID {department_id}): {len(nodes)} nodes, "
              f"{len(sub_edges)} internal edges")
        write_csv(os.path.join(args.out_dir, f"T2_dept{label}_nodes.csv"),
                  ["node"], [[node] for node in nodes])
        write_csv(os.path.join(args.out_dir, f"T2_dept{label}_edges.csv"),
                  ["source", "target"], sub_edges)

    # T4-T6: triad census, null model, and Z-scores.
    analyzed = {}
    for label in ("A", "B"):
        department_id, _, sub_edges = subgraphs[label]
        analyzed[label] = analyze_department(sub_edges, params, label)
        if args.check_networkx:
            print(f"Cross-check Department {label} (ID {department_id}):")
            cross_check_networkx(sub_edges, analyzed[label][0])

    rows_counts = [
        [k, CANONICAL_CODES[k - 1], analyzed["A"][0][k], analyzed["B"][0][k]]
        for k in range(1, 14)
    ]
    rows_z = [[k, f"{analyzed['A'][1][k]:.6f}", f"{analyzed['B'][1][k]:.6f}"]
              for k in range(1, 14)]
    write_csv(os.path.join(args.out_dir, "T4_counts.csv"),
              ["class", "canonical_code", f"department_{params['dept_A']}",
               f"department_{params['dept_B']}"], rows_counts)
    write_csv(os.path.join(args.out_dir, "T6_zscores.csv"),
              ["class", f"Z_department_{params['dept_A']}",
               f"Z_department_{params['dept_B']}"], rows_z)
    rows_null = [[k, analyzed["A"][2][k], analyzed["A"][3][k],
                  analyzed["B"][2][k], analyzed["B"][3][k]] for k in range(1, 14)]
    write_csv(os.path.join(args.out_dir, "T6_null_stats.csv"),
              ["class", "A_mean", "A_population_std", "B_mean", "B_population_std"],
              rows_null)

    print("\nT4 triad counts (class, code, Department A, Department B):")
    for row in rows_counts:
        print(row)
    print("\nT6 Z-scores (class, Department A, Department B):")
    for row in rows_z:
        print(row)

    # T7: grouped bar chart of Z-scores.
    classes = list(range(1, 14))
    width = 0.38
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar([k - width / 2 for k in classes],
           [analyzed["A"][1][k] for k in classes], width,
           label=f"Department {params['dept_A']} (A)")
    ax.bar([k + width / 2 for k in classes],
           [analyzed["B"][1][k] for k in classes], width,
           label=f"Department {params['dept_B']} (B)")
    ax.axhline(2, linestyle="--", color="black", linewidth=1, label="Z = +2 / -2 threshold")
    ax.axhline(-2, linestyle="--", color="black", linewidth=1)
    ax.axhline(0, color="gray", linewidth=0.7)
    ax.set_xticks(classes)
    ax.set_xlabel("Triad class (1-13)")
    ax.set_ylabel("Z-score")
    ax.set_title(f"Network motif Z-scores: Departments {params['dept_A']} vs "
                 f"{params['dept_B']} (R = {R})")
    ax.legend()
    fig.tight_layout()
    chart_path = os.path.join(args.out_dir, "T7_zscores.png")
    fig.savefig(chart_path, dpi=180)
    plt.close(fig)
    print(f"\nT7 chart saved to {chart_path}")
    print(f"All outputs are in: {os.path.abspath(args.out_dir)}")


if __name__ == "__main__":
    main()
