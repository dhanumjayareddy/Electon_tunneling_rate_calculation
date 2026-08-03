#!/usr/bin/env python3
"""
Combined workflow for:
1. Computing chromophore-to-chromophore distance matrices from a PDB file.
2. Calculating Page-Moser-Dutton ET rates from a distance matrix.
"""

import sys
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import PDBParser
    

RHO = 0.76
BETA = 1.2 - 0.8 * RHO

CHROMOPHORE_ATOMS = {
    "TRP": "ALL",
    "TYR": "ALL",
    "HIS": "ALL",
    "MET": "ALL",
    "CYS": "ALL"
}

RESIDUE_ORDER = {
    "TRP": 1,
    "TYR": 2,
    "HIS": 3,
    "MET": 4,
    "CYS": 5
}

DELTA_G_VALUES = {
    ('TRP', 'TRP'): 0.0,
    ('TRP', 'TYR'): 0.35265974,
    ('TRP', 'CYS'): 0.45497461,
    ('TRP', 'MET'): 0.34585689,
    ('TRP', 'HIS'): 0.42667475,
    ('TYR', 'TRP'): -0.35265974,
    ('TYR', 'TYR'): 0.0,
    ('TYR', 'CYS'): 0.10231486,
    ('TYR', 'MET'): -0.00680285,
    ('TYR', 'HIS'): 0.07401501,
    ('CYS', 'TRP'): -0.45497461,
    ('CYS', 'TYR'): -0.10231486,
    ('CYS', 'CYS'): 0.0,
    ('CYS', 'MET'): -0.10911771,
    ('CYS', 'HIS'): -0.02829986,
    ('MET', 'TRP'): -0.34585689,
    ('MET', 'TYR'): 0.00680285,
    ('MET', 'CYS'): 0.10911771,
    ('MET', 'MET'): 0.0,
    ('MET', 'HIS'): 0.08081786,
    ('HIS', 'TRP'): -0.42667475,
    ('HIS', 'TYR'): -0.07401501,
    ('HIS', 'CYS'): 0.02829986,
    ('HIS', 'MET'): -0.08081786,
    ('HIS', 'HIS'): 0.0
}

LAMBDA_VALUES = {
    ('TRP', 'TRP'): None,
    ('TRP', 'TYR'): 2.29024748099,
    ('TRP', 'CYS'): 2.2440969466,
    ('TRP', 'MET'): 2.3611739951,
    ('TRP', 'HIS'): 2.56104172809,
    ('TYR', 'TRP'): 2.29024748099,
    ('TYR', 'TYR'): None,
    ('TYR', 'CYS'): 3.9768916872,
    ('TYR', 'MET'): 4.0939687357,
    ('TYR', 'HIS'): 4.29383646869,
    ('CYS', 'TRP'): 2.2440969466,
    ('CYS', 'TYR'): 3.9768916872,
    ('CYS', 'CYS'): None,
    ('CYS', 'MET'): 4.0478182013,
    ('CYS', 'HIS'): 4.24768593429,
    ('MET', 'TRP'): 2.3611739951,
    ('MET', 'TYR'): 4.0939687357,
    ('MET', 'CYS'): 4.0478182013,
    ('MET', 'MET'): None,
    ('MET', 'HIS'): 4.3647629828,
    ('HIS', 'TRP'): 2.56104172809,
    ('HIS', 'TYR'): 4.29383646869,
    ('HIS', 'CYS'): 4.24768593429,
    ('HIS', 'MET'): 4.3647629828,
    ('HIS', 'HIS'): None
}


def load_structure(pdb_file):
    if PDBParser is None:
        raise ImportError("BioPython is required to process PDB files")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    model = structure[0]
    residues = []

    for chain in model:
        for residue in chain:
            resname = residue.get_resname()

            if resname not in CHROMOPHORE_ATOMS:
                continue

            atoms = []

            for atom in residue:
                if CHROMOPHORE_ATOMS[resname] == "ALL":
                    atoms.append({
                        "name": atom.get_name(),
                        "coord": atom.get_coord()
                    })
                elif atom.get_name() in CHROMOPHORE_ATOMS[resname]:
                    atoms.append({
                        "name": atom.get_name(),
                        "coord": atom.get_coord()
                    })

            if len(atoms) == 0:
                continue

            residues.append({
                "resname": resname,
                "chain": chain.id,
                "resid": residue.get_id()[1],
                "atoms": atoms
            })

    return residues


def sort_residues(residues):
    return sorted(
        residues,
        key=lambda r: (
            RESIDUE_ORDER[r["resname"]],
            r["chain"],
            r["resid"]
        )
    )


def residue_label(residue):
    return f"{residue['resname']}_{residue['chain']}_{residue['resid']}"


def calculate_r(donor, acceptor):
    min_distance = float("inf")

    donor_atom = None
    acceptor_atom = None

    for atom_a, atom_b in product(donor["atoms"], acceptor["atoms"]):
        d = np.linalg.norm(atom_a["coord"] - atom_b["coord"])

        if d < min_distance:
            min_distance = d
            donor_atom = atom_a["name"]
            acceptor_atom = atom_b["name"]

    return round(float(min_distance), 3), donor_atom, acceptor_atom


def build_atom_list(residues):
    rows = []

    for residue in residues:
        atom_names = [atom["name"] for atom in residue["atoms"]]
        atom_list = ", ".join(atom_names)

        rows.append({
            "residue": residue_label(residue),
            "resname": residue["resname"],
            "chain": residue["chain"],
            "resid": residue["resid"],
            "num_atoms": len(residue["atoms"]),
            "atom_names": atom_list
        })

    return pd.DataFrame(rows)


def build_pairwise_table(residues):
    rows = []

    for donor in residues:
        for acceptor in residues:
            if donor == acceptor:
                continue

            r_value, donor_atom, acceptor_atom = calculate_r(donor, acceptor)

            rows.append({
                "donor": residue_label(donor),
                "acceptor": residue_label(acceptor),
                "donor_atom": donor_atom,
                "acceptor_atom": acceptor_atom,
                "R": r_value
            })

    return pd.DataFrame(rows)


def build_r_matrix(residues):
    labels = [residue_label(r) for r in residues]

    matrix = pd.DataFrame(index=labels, columns=labels)

    for donor in residues:
        donor_label = residue_label(donor)

        for acceptor in residues:
            acceptor_label = residue_label(acceptor)

            if donor_label == acceptor_label:
                matrix.loc[acceptor_label, donor_label] = 0.0
                continue

            r_value, _, _ = calculate_r(donor, acceptor)

            matrix.loc[acceptor_label, donor_label] = r_value

    return matrix


def extract_residue_info(residue_label_value):
    """Extract residue type, chain, and number from label like 'TRP_A_290'"""
    parts = residue_label_value.split('_')
    if len(parts) >= 3:
        return parts[0], parts[1], parts[2]
    return None, None, None


def calculate_pmd_rate(R, donor_type, acceptor_type):
    """
    Calculate PMD rate using equation with R-3.6 term:
    log10(k_ET) = 13 - beta*(R-3.6) - 3.1*(DeltaG + lambda)^2/lambda
    where beta = 1.2 - 0.8*rho and rho = 0.76
    """
    dG = DELTA_G_VALUES.get((donor_type, acceptor_type))
    lam = LAMBDA_VALUES.get((donor_type, acceptor_type))

    if dG is None or lam is None:
        return None, 'values_not_found'

    log_k = 13 - BETA * (R - 3.6) - (3.1 * (dG + lam) ** 2 / lam)
    return log_k, 'ok'


def build_rate_summary(rate_matrix, df, valid_pairs, valid_rates):
    """Create summary statistics and top candidate pairs for the rate matrix."""
    summary_rows = []

    if len(valid_rates) > 0:
        summary_rows.extend([
            {"metric": "total_pairs_analyzed", "value": valid_pairs, "unit": "pairs"},
            {"metric": "rate_constant_min", "value": float(np.min(valid_rates)), "unit": "s^-1"},
            {"metric": "rate_constant_max", "value": float(np.max(valid_rates)), "unit": "s^-1"},
            {"metric": "rate_constant_mean", "value": float(np.mean(valid_rates)), "unit": "s^-1"},
            {"metric": "rate_constant_median", "value": float(np.median(valid_rates)), "unit": "s^-1"},
        ])
    else:
        summary_rows.extend([
            {"metric": "total_pairs_analyzed", "value": valid_pairs, "unit": "pairs"},
            {"metric": "rate_constant_min", "value": np.nan, "unit": "s^-1"},
            {"metric": "rate_constant_max", "value": np.nan, "unit": "s^-1"},
            {"metric": "rate_constant_mean", "value": np.nan, "unit": "s^-1"},
            {"metric": "rate_constant_median", "value": np.nan, "unit": "s^-1"},
        ])

    log_rates = np.log10(valid_rates) if len(valid_rates) > 0 else []
    if len(log_rates) > 0:
        summary_rows.extend([
            {"metric": "ultrafast_contact_count", "value": int(np.sum(log_rates >= 12)), "unit": "pairs"},
            {"metric": "fast_functional_count", "value": int(np.sum((log_rates >= 10) & (log_rates < 12))), "unit": "pairs"},
            {"metric": "moderate_count", "value": int(np.sum((log_rates >= 7) & (log_rates < 10))), "unit": "pairs"},
            {"metric": "slow_count", "value": int(np.sum((log_rates >= 4) & (log_rates < 7))), "unit": "pairs"},
            {"metric": "negligible_count", "value": int(np.sum(log_rates < 4)), "unit": "pairs"},
        ])

    summary_df = pd.DataFrame(summary_rows)

    top_pairs = []
    if len(valid_rates) > 0:
        mask = rate_matrix.values != 'NA'
        rows, cols = np.where(mask)
        valid_rates_array = np.array(valid_rates, dtype=float)
        ranked_indices = np.argsort(-valid_rates_array)[:10]
        for idx in ranked_indices:
            donor = rate_matrix.index[rows[idx]]
            acceptor = rate_matrix.columns[cols[idx]]
            distance = df.loc[donor, acceptor]
            k_et = valid_rates_array[idx]
            log10_rate = np.log10(k_et)
            if log10_rate >= 12:
                category = "ultrafast"
            elif log10_rate >= 10:
                category = "fast"
            elif log10_rate >= 7:
                category = "moderate"
            elif log10_rate >= 4:
                category = "slow"
            else:
                category = "negligible"
            top_pairs.append({
                "rank": len(top_pairs) + 1,
                "donor": donor,
                "acceptor": acceptor,
                "distance_A": float(distance),
                "rate_s^-1": float(k_et),
                "log10_rate": float(log10_rate),
                "category": category,
            })

    top_pairs_df = pd.DataFrame(top_pairs)
    return summary_df, top_pairs_df


def write_rate_output_with_summary(rate_matrix, output_file, df, valid_pairs, valid_rates):
    """Write the rate matrix plus summary statistics and top candidate pairs to one CSV file."""
    summary_df, top_pairs_df = build_rate_summary(rate_matrix, df, valid_pairs, valid_rates)

    with open(output_file, "w", encoding="utf-8") as handle:
        handle.write("# Units: electron transfer rate constants in s^-1\n")
        handle.write("# Values are raw rate constants; larger values indicate faster ET\n")
        rate_matrix.to_csv(handle)

        handle.write("\n# Summary statistics\n")
        summary_df.to_csv(handle, index=False)

        handle.write("\n# Top candidate pairs\n")
        top_pairs_df.to_csv(handle, index=False)


def process_distance_matrix(matrix_file, output_file):
    """Process distance matrix and calculate PMD rates for all pairs"""
    print(f"Loading distance matrix from: {matrix_file}")
    df = pd.read_csv(matrix_file, index_col=0)

    rate_matrix = pd.DataFrame(index=df.index, columns=df.columns, dtype=object)

    print(f"Calculating PMD rates for {len(df)} chromophores...")

    valid_pairs = 0

    for donor in df.index:
        for acceptor in df.columns:
            if donor == acceptor:
                rate_matrix.loc[donor, acceptor] = 'NA'
                continue

            R = df.loc[donor, acceptor]
            if pd.isna(R) or float(R) <= 0:
                rate_matrix.loc[donor, acceptor] = 'NA'
                continue

            donor_type, donor_chain, donor_num = extract_residue_info(donor)
            acceptor_type, acceptor_chain, acceptor_num = extract_residue_info(acceptor)

            if not all([donor_type, acceptor_type]):
                rate_matrix.loc[donor, acceptor] = 'NA'
                continue

            log_k, status = calculate_pmd_rate(R, donor_type, acceptor_type)

            if log_k is None:
                rate_matrix.loc[donor, acceptor] = 'NA'
                continue

            k_ET = 10 ** log_k
            rate_matrix.loc[donor, acceptor] = k_ET
            valid_pairs += 1

    valid_rates = [x for x in rate_matrix.values.flatten() if x != 'NA']

    print(f"Saving rate matrix to: {output_file}")
    write_rate_output_with_summary(rate_matrix, output_file, df, valid_pairs, valid_rates)

    print(f"\n{'=' * 50}")
    print("PMD RATE CALCULATION SUMMARY")
    print(f"{'=' * 50}")
    print(f"Total pairs analyzed: {valid_pairs}")
    if len(valid_rates) > 0:
        valid_rates = np.array(valid_rates, dtype=float)
        print(f"Rate constant range: {valid_rates.min():.2e} to {valid_rates.max():.2e} s⁻¹")

        log_rates = np.log10(valid_rates)
        ultrafast = np.sum(log_rates >= 12)
        fast = np.sum((log_rates >= 10) & (log_rates < 12))
        moderate = np.sum((log_rates >= 7) & (log_rates < 10))
        slow = np.sum((log_rates >= 4) & (log_rates < 7))
        negligible = np.sum(log_rates < 4)

        print(f"\nET Category Distribution:")
        print(f"  Ultrafast (contact): {ultrafast}")
        print(f"  Fast (functional): {fast}")
        print(f"  Moderate: {moderate}")
        print(f"  Slow: {slow}")
        print(f"  Negligible: {negligible}")

        flat_indices = np.argsort(-valid_rates)[:10]
        mask = rate_matrix.values != 'NA'
        rows, cols = np.where(mask)
        print(f"\nTop 10 Fastest ET Pairs:")
        for i, idx in enumerate(flat_indices):
            donor = rate_matrix.index[rows[idx]]
            acceptor = rate_matrix.columns[cols[idx]]
            R = df.loc[donor, acceptor]
            k_ET = valid_rates[idx]
            print(f"  {i+1}. {donor} -> {acceptor}: R={R:.2f}Å, k={k_ET:.2e} s⁻¹")

    print(f"\n{'=' * 50}")


def generate_heatmap(rate_matrix, output_file):
    """Generate and save heatmap visualization of the rate matrix"""
    if plt is None or sns is None:
        print("matplotlib/seaborn not available; skipping heatmap generation")
        return

    print("Generating heatmap...")

    viz_matrix = rate_matrix.copy()
    viz_matrix = viz_matrix.replace('NA', np.nan)
    viz_matrix = viz_matrix.astype(float)
    viz_matrix = np.log10(viz_matrix)

    plt.figure(figsize=(12, 10))
    sns.heatmap(
        viz_matrix,
        cmap='viridis',
        cbar_kws={'label': 'log₁₀(k_ET) [s⁻¹]'},
        square=True,
        annot=False,
        fmt='.1f'
    )

    plt.title('PMD Electron Transfer Rate Matrix', fontsize=14, fontweight='bold')
    plt.xlabel('Acceptor', fontsize=12)
    plt.ylabel('Donor', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()

    heatmap_file = output_file.replace('.csv', '_heatmap.png')
    plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Heatmap saved to: {heatmap_file}")


def process_pdb_file(pdb_file):
    residues = load_structure(pdb_file)
    residues = sort_residues(residues)

    atom_list_df = build_atom_list(residues)
    pairwise_df = build_pairwise_table(residues)
    matrix_df = build_r_matrix(residues)

    pdb_basename = pdb_file.replace('.pdb', '')
    pairwise_file = f"{pdb_basename}_allatoms_edge_to_edge_pairwise_table.csv"
    matrix_file = f"{pdb_basename}_allatoms_edge_to_edge_r_matrix.csv"
    atom_list_file = f"{pdb_basename}_chromophore_atom_list.csv"
    rate_output_file = f"{matrix_file.replace('.csv', '_pmd_rates.csv')}"

    atom_list_df.to_csv(atom_list_file, index=False)
    pairwise_df.to_csv(pairwise_file, index=False)
    matrix_df.to_csv(matrix_file)

    print("Done")
    print("Saved:")
    print(f"  {atom_list_file}")
    print(f"  {pairwise_file}")
    print(f"  {matrix_file}")

    process_distance_matrix(matrix_file, rate_output_file)
    print(f"\nResults saved to: {rate_output_file}")


def main():
    """Main function"""
    if len(sys.argv) != 2:
        print("Usage:")
        print("  python pmd_rates_simple.py protein.pdb")
        print("  python pmd_rates_simple.py distance_matrix.csv")
        sys.exit(1)

    input_path = sys.argv[1]

    if not Path(input_path).exists():
        print(f"Error: File '{input_path}' not found")
        sys.exit(1)

    suffix = Path(input_path).suffix.lower()

    try:
        if suffix == '.pdb':
            process_pdb_file(input_path)
        else:
            output_file = input_path.replace('.csv', '_pmd_rates.csv')
            process_distance_matrix(input_path, output_file)
            print(f"\nResults saved to: {output_file}")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

