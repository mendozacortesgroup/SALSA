#!/usr/bin/env python3
"""
Generate substitutions_textual.txt file from ICSD compound database.
Updated to work with SALSA workflow manager.
"""

import pandas as pd
import numpy as np
import re
import sys
import os
import argparse
from datetime import datetime


def name_ion(element, charge):
    """Format ion name as Element:Charge+/- (e.g., Na:1+, O:2-)"""
    sign_string = "+" if charge >= 0 else "-"
    ion_name = f"{element}:{abs(charge)}{sign_string}"
    return ion_name


def load_element_charge_dictionary(filepath):
    """Load element charge dictionary and create necessary mappings"""
    atomic_symbol_dict = {}
    ionic_charge_dict = {}

    with open(filepath, "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            if line:
                atomic_symbol = line[0]
                atomic_number = int(line[1])
                ionic_charges = line[2:]
                atomic_symbol_dict[atomic_number] = atomic_symbol
                ionic_charge_dict[atomic_symbol] = tuple(int(q) for q in ionic_charges)

    # Create list of all possible ions
    possible_ions = []
    for element, charges in ionic_charge_dict.items():
        possible_ions += [name_ion(element, charge) for charge in np.sort(charges)]

    possible_ion_index_dict = {ion: ion_index for ion_index, ion in enumerate(possible_ions)}

    return atomic_symbol_dict, ionic_charge_dict, possible_ions, possible_ion_index_dict


def cost_function(elements, element_counts, charge_indices, ionic_charge_dict):
    """Calculate cost for a given charge assignment"""
    charge_indices = np.array(charge_indices)

    ionic_charges = [ionic_charge_dict[elements[i]][charge_indices[i]] for i in range(len(elements))]
    net_abs_charge = abs(np.dot(element_counts, ionic_charges))
    cost_per_electron = 1000
    priority_cost = charge_indices.sum()
    cost = net_abs_charge * cost_per_electron + priority_cost
    return net_abs_charge, cost


def find_best_charge_indices(elements, element_counts, ionic_charge_dict):
    """Find the best charge assignment for a compound"""
    charge_options_counts = [len(ionic_charge_dict[element]) for element in elements]

    charge_indices_candidates = [list(np.zeros(len(elements), dtype=np.int64))]
    charge_indices_candidates_evaluated = []
    best_charge_indices_candidates = []
    current_cost = np.inf
    best_cost = np.inf
    best_charge = np.inf
    cost_threshold = 1
    max_cost = 10000
    max_charge = 0.1
    suitable_candidate_found = False

    while len(charge_indices_candidates) and current_cost:
        charge_indices_candidate = charge_indices_candidates.pop(0)
        current_charge, current_cost = cost_function(elements, element_counts, charge_indices_candidate, ionic_charge_dict)

        if current_cost <= cost_threshold:
            suitable_candidate_found = True

        if current_cost < best_cost:
            best_cost = current_cost
            best_charge_indices_candidates = [charge_indices_candidate]
            best_charge = current_charge
        elif current_cost == best_cost:
            best_charge_indices_candidates.append(charge_indices_candidate)
            best_charge = min(best_charge, current_charge)

        charge_indices_candidates_evaluated.append(list(charge_indices_candidate))

        if not suitable_candidate_found:
            for element_index in range(len(elements)):
                potential_charge_indices_candidate = charge_indices_candidate.copy()
                potential_charge_indices_candidate[element_index] += 1
                if potential_charge_indices_candidate[element_index] < charge_options_counts[element_index]:
                    if list(potential_charge_indices_candidate) not in charge_indices_candidates + charge_indices_candidates_evaluated:
                        if np.count_nonzero(potential_charge_indices_candidate) <= 1:
                            charge_indices_candidates.append(potential_charge_indices_candidate)

    if best_cost < max_cost and best_charge < max_charge:
        return best_charge_indices_candidates
    else:
        return []


def process_icsd_compounds(icsd_filepath, ionic_charge_dict, possible_ions, possible_ion_index_dict):
    """Process ICSD compounds and translate to ionic formulas"""
    print("Loading ICSD compounds...")
    icsd_db = pd.read_csv(icsd_filepath, header=None)
    icsd_db.columns = ["Chemical_Formula"]
    n_icsd_compounds = len(icsd_db)
    print(f"Loaded {n_icsd_compounds} compounds")

    # Initialize columns
    icsd_db["N_Unique_Elements"] = ""
    icsd_db["Unique_Elements"] = ""
    icsd_db["Unique_Elements"] = icsd_db["Unique_Elements"].astype('object')
    icsd_db["Element_Counts"] = ""
    icsd_db["Element_Counts"] = icsd_db["Element_Counts"].astype('object')
    icsd_db["Ionic_Charges"] = ""
    icsd_db["Ionic_Charges"] = icsd_db["Ionic_Charges"].astype('object')
    icsd_db["Unique_Ions"] = ""
    icsd_db["Unique_Ions"] = icsd_db["Unique_Ions"].astype('object')
    icsd_db["Ion_Vector"] = ""
    icsd_db["Ion_Vector"] = icsd_db["Ion_Vector"].astype('object')

    ion_frequencies = [0] * len(possible_ions)

    # Process each compound
    print("Processing compounds...")
    deuterium_count = 0
    multiplicity_count = 0
    no_arrangements_count = 0
    too_many_combinations_count = 0
    max_element_count = 100

    for index, series in icsd_db.iterrows():
        if index % 10000 == 0:
            print(f"Processing compound {index}/{n_icsd_compounds}")

        original_compound_string = series["Chemical_Formula"]
        parsed_compound_array = np.array(re.findall(r'([A-z]+)(\d*\.?\d+)', original_compound_string)).transpose()

        n_unique_elements = len(parsed_compound_array[0])
        unique_elements = list(parsed_compound_array[0])
        element_counts = list(float(freq) for freq in parsed_compound_array[1])

        icsd_db.at[index, "N_Unique_Elements"] = n_unique_elements
        icsd_db.at[index, "Unique_Elements"] = unique_elements
        icsd_db.at[index, "Element_Counts"] = element_counts

        # Skip problematic compounds
        if len(unique_elements) > max_element_count:
            too_many_combinations_count += 1
            continue

        if "D" in unique_elements or "T" in unique_elements:
            deuterium_count += 1
            continue

        best_charge_indices = find_best_charge_indices(unique_elements, element_counts, ionic_charge_dict)

        if len(best_charge_indices) > 1:
            multiplicity_count += 1
            continue
        elif len(best_charge_indices) == 0:
            no_arrangements_count += 1
            continue

        my_charge_indices = best_charge_indices[0]

        ion_charges = []
        ion_names = []
        ion_vector = np.zeros(len(possible_ions))

        for i in range(n_unique_elements):
            element = unique_elements[i]
            ion_charge = ionic_charge_dict[element][my_charge_indices[i]]
            ion_charges.append(ion_charge)
            ion_name = name_ion(element, ion_charge)
            ion_names.append(ion_name)
            ion_index = possible_ion_index_dict[ion_name]
            ion_frequencies[ion_index] += 1
            ion_vector[ion_index] = element_counts[i]

        icsd_db.at[index, "Ionic_Charges"] = ion_charges
        icsd_db.at[index, "Unique_Ions"] = ion_names
        icsd_db.at[index, "Ion_Vector"] = ion_vector

    print(f"\nProcessing complete:")
    print(f"Deuterium/Tritium compounds: {deuterium_count}")
    print(f"Multiple charge arrangements: {multiplicity_count}")
    print(f"No valid arrangements: {no_arrangements_count}")
    print(f"Too many elements: {too_many_combinations_count}")

    # Filter to only translated compounds
    icsd_db_translated = icsd_db[icsd_db.Ionic_Charges != ""]
    print(f"Successfully translated {len(icsd_db_translated)} compounds")

    # Get observed ions
    observed_ions = np.array(possible_ions)[np.array(ion_frequencies) != 0]
    observed_ion_index_dict = {ion: ion_index for ion_index, ion in enumerate(observed_ions)}

    return icsd_db_translated, ion_frequencies, observed_ions, observed_ion_index_dict, possible_ion_index_dict


def build_substitution_matrix(icsd_db_translated, observed_ions, observed_ion_index_dict, ion_frequencies, possible_ion_index_dict):
    """Build the substitution matrix from compound data following MATLAB implementation"""
    print("\nBuilding substitution matrix...")
    n_observed_ions = len(observed_ions)
    n_possible_ions = len(possible_ion_index_dict)

    # Build matrix of compounds
    print("Building compound matrix...")
    compound_matrix = []
    for index, series in icsd_db_translated.iterrows():
        ion_vec = np.zeros(n_possible_ions)
        if isinstance(series.Ion_Vector, np.ndarray):
            ion_vec = series.Ion_Vector
        compound_matrix.append((ion_vec > 0).astype(int))

    compound_matrix = np.array(compound_matrix)

    # Remove duplicates and track repetitions
    print("Removing duplicate compounds...")
    unique_compounds, unique_indices, inverse_indices = np.unique(
        compound_matrix, axis=0, return_index=True, return_inverse=True
    )
    repetition_counts = np.bincount(inverse_indices)

    print(f"Reduced from {len(compound_matrix)} to {len(unique_compounds)} unique compounds")

    # Work with unique compounds
    unique_matrix = unique_compounds
    occurences = np.sum(unique_matrix, axis=0)

    # Count ions per compound
    ion_counts = np.sum(unique_matrix, axis=1)
    min_ions = 2
    max_ions = int(max(ion_counts))

    # Find all compound pairs with similar environments
    num_matches = np.zeros((n_possible_ions, n_possible_ions))

    print("Finding ion substitution pairs...")
    for n_ions in range(min_ions, max_ions + 1):
        indices_with_n = np.where(ion_counts == n_ions)[0]
        if len(indices_with_n) < 2:
            continue

        print(f"  Processing {len(indices_with_n)} compounds with {n_ions} ions")

        # Get the submatrix for compounds with n ions
        cur_matrix = unique_matrix[indices_with_n]

        # Find all pairs that differ by exactly one ion substitution
        for i in range(len(cur_matrix)):
            for j in range(i + 1, len(cur_matrix)):
                difference = cur_matrix[i] + cur_matrix[j] - 2 * cur_matrix[i] * cur_matrix[j]
                different_positions = np.where(difference == 1)[0]

                if len(different_positions) == 2:
                    # Weight by the minimum repetition count of the two compounds
                    weight = min(repetition_counts[indices_with_n[i]], repetition_counts[indices_with_n[j]])
                    num_matches[different_positions[0], different_positions[1]] += weight
                    num_matches[different_positions[1], different_positions[0]] += weight

    # Calculate fractional matches
    print("Calculating fractional matches...")
    fractional_matches = np.zeros((n_possible_ions, n_possible_ions))

    for i in range(n_possible_ions):
        for j in range(n_possible_ions):
            if occurences[i] > 0 and occurences[j] > 0:
                fractional_matches[i, j] = num_matches[i, j] / max(1, np.sqrt(occurences[i] * occurences[j]))

            # Filter out rare ions
            if occurences[i] <= 10 or occurences[j] <= 10:
                fractional_matches[i, j] = np.nan

        fractional_matches[i, i] = 1  # Set diagonal to 1

    # Calculate mean fraction for normalization
    frac_vectorized = fractional_matches.flatten()
    good_idx = ~np.isnan(frac_vectorized)
    mean_fraction = np.mean(frac_vectorized[good_idx])
    print(f"Mean fraction: {mean_fraction}")

    # Extract only observed ions
    observed_indices = [possible_ion_index_dict[ion] for ion in observed_ions]
    fractional_matches_observed = fractional_matches[np.ix_(observed_indices, observed_indices)]

    return fractional_matches_observed, mean_fraction


def write_substitutions_textual(fractional_matches, mean_fraction, observed_ions, output_filepath):
    """Write substitutions to textual format"""
    print(f"\nWriting substitutions to {output_filepath}...")

    with open(output_filepath, 'w') as f:
        f.write("ion_A ion_B alpha_AB\n")

        count = 0
        for i in range(len(observed_ions)):
            for j in range(len(observed_ions)):
                if i != j and not np.isnan(fractional_matches[i, j]) and fractional_matches[i, j] > 0:
                    # Normalize by mean fraction as in MATLAB code
                    alpha_value = fractional_matches[i, j] / mean_fraction
                    f.write(f"{observed_ions[i]} {observed_ions[j]} {alpha_value:.6f}\n")
                    count += 1

    print(f"Wrote {count} substitution pairs")


def main():
    """Main function to generate substitutions_textual.txt"""
    parser = argparse.ArgumentParser(description="Generate substitution matrix from ICSD data")
    parser.add_argument('--dataset', required=True, help="Path to dataset directory")
    parser.add_argument('--output', required=True, help="Path to project directory for outputs")

    args = parser.parse_args()

    # Get SALSA_DIR environment variable
    salsa_dir = os.getenv("SALSA_DIR")
    if not salsa_dir:
        print("Error: SALSA_DIR environment variable not set")
        return

    # Set up paths using universal data directory
    dataset_dir = args.dataset
    project_dir = args.output

    # Use existing substitution_result directory (created during project setup)
    substitution_result_dir = os.path.join(project_dir, "substitution_result")

    # Input files from dataset directory
    element_charge_dict_file = os.path.join(dataset_dir, "Element_Charge_Dictionary.txt")
    icsd_file = os.path.join(dataset_dir, "ICSD-compositions.csv")

    # Output file in substitution_result directory
    output_file = os.path.join(substitution_result_dir, "substitutions_textual.txt")

    # Check if input files exist
    if not os.path.exists(element_charge_dict_file):
        print(f"Error: Element charge dictionary not found at {element_charge_dict_file}")
        return

    if not os.path.exists(icsd_file):
        print(f"Error: ICSD compositions file not found at {icsd_file}")
        return

    print(f"Starting substitution matrix generation at {datetime.now()}")

    # Load element charge dictionary
    atomic_symbol_dict, ionic_charge_dict, possible_ions, possible_ion_index_dict = \
        load_element_charge_dictionary(element_charge_dict_file)
    print(f"Loaded {len(ionic_charge_dict)} elements with {len(possible_ions)} possible ions")

    # Process ICSD compounds
    icsd_db_translated, ion_frequencies, observed_ions, observed_ion_index_dict, possible_ion_index_dict = \
        process_icsd_compounds(icsd_file, ionic_charge_dict, possible_ions, possible_ion_index_dict)

    # Build substitution matrix
    fractional_matches, mean_fraction = build_substitution_matrix(
        icsd_db_translated, observed_ions, observed_ion_index_dict, ion_frequencies, possible_ion_index_dict
    )

    # Write output
    write_substitutions_textual(fractional_matches, mean_fraction, observed_ions, output_file)

    print(f"\nCompleted at {datetime.now()}")
    print(f"Output written to: {output_file}")


if __name__ == "__main__":
    main()
