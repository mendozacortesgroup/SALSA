#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re
import os
import argparse
import matplotlib.pyplot as plt
from IPython.display import display

def Semiconductors_Workflow(dataset_dir, project_dir):
    """
    Main semiconductors workflow function
    """

    # Use existing substitution_result directory (created during project setup)
    substitution_result_dir = os.path.join(project_dir, "substitution_result")

    SMALL_SIZE = 11
    MEDIUM_SIZE = 13
    BIG_SIZE = 15

    plt.rc('font', size=SMALL_SIZE)
    plt.rc('axes', titlesize=BIG_SIZE)
    plt.rc('axes', labelsize=MEDIUM_SIZE)
    plt.rc('xtick', labelsize=SMALL_SIZE)
    plt.rc('ytick', labelsize=SMALL_SIZE)
    plt.rc('legend', fontsize=SMALL_SIZE)
    plt.rc('figure', titlesize=BIG_SIZE)

    def refresh_icsd_db():
        icsd_filename = os.path.join(dataset_dir, "ICSD-compositions.csv")
        icsd_db = pd.read_csv(icsd_filename, header=None)
        icsd_db.columns = ["Chemical_Formula"]
        return icsd_db

    icsd_db = refresh_icsd_db()
    n_icsd_compounds = len(icsd_db)
    print("{} compounds imported into icsd_db, such as:".format(n_icsd_compounds))
    display(icsd_db.head())

    def name_ion(element, charge):
        sign_string = "+" if charge >= 0 else "-"
        ion_name = element + str(abs(charge)) + sign_string
        return ion_name

    element_charge_dict_filename = os.path.join(dataset_dir, "Element_Charge_Dictionary.txt")
    atomic_symbol_dict ={}
    ionic_charge_dict = {}
    with open(element_charge_dict_filename, "r") as element_charge_dict_file:
        lines = element_charge_dict_file.readlines()
        for line in lines:
            line = line.split()
            if line:
                atomic_symbol = line[0]
                atomic_number = int(line[1])
                ionic_charges = line[2:]
                atomic_symbol_dict[ atomic_number ] = atomic_symbol
                ionic_charge_dict[ atomic_symbol ]  = tuple(int(q) for q in ionic_charges)
    n_elements = len(atomic_symbol_dict.keys())
    n_element_charge_pairs = sum(len(val) for val in ionic_charge_dict.values())

    possible_ions = []
    for element, charges in ionic_charge_dict.items():
        possible_ions += [name_ion(element, charge) for charge in np.sort(charges)]
    n_possible_ions = len(possible_ions)
    possible_ion_index_dict = {ion:ion_index for ion_index, ion in enumerate(possible_ions)}

    print("atomic_symbol_dict maps atomic numbers to atomic symbols for {} elements.".format(n_elements))
    print("ionic_charge_dict maps atomic symbols to common ionic charges with a total of {} element-charge pairs.".format(n_element_charge_pairs))
    print("possible_ions lists all {} possible ions".format(n_possible_ions))

    def cost_function(elements, element_counts, charge_indices):
        charge_indices = np.array(charge_indices)

        ionic_charges = [ionic_charge_dict[elements[i]][charge_indices[i]] for i in range(len(elements))]
        net_abs_charge = abs(np.dot(element_counts, ionic_charges))
        cost_per_electron = 1000
        priority_cost = charge_indices.sum()
        cost = net_abs_charge * cost_per_electron + priority_cost
        return net_abs_charge, cost

    def find_best_charge_indices(elements, element_counts):
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

        counter = 0
        while len(charge_indices_candidates) and current_cost:
            counter += 1

            charge_indices_candidate = charge_indices_candidates.pop(0)
            current_charge, current_cost  = cost_function(elements, element_counts, charge_indices_candidate)
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

    # Process compounds and build matrices
    icsd_db = refresh_icsd_db()
    n_icsd_compounds = len(icsd_db)

    icsd_db["N_Unique_Elements"] = ""
    icsd_db["Unique_Elements"] = ""
    icsd_db["Unique_Elements"] = icsd_db["Unique_Elements"].astype('object')
    icsd_db["Element_Counts"] = ""
    icsd_db["Element_Counts"] = icsd_db["Element_Counts"].astype('object')
    icsd_db["Element_Count_Dictionary"] = ""
    icsd_db["Element_Count_Dictionary"] = icsd_db["Element_Count_Dictionary"].astype('object')
    icsd_db["Element_Charge_Dictionary"] = ""
    icsd_db["Element_Charge_Dictionary"] = icsd_db["Element_Charge_Dictionary"].astype('object')
    icsd_db["Ionic_Charges"] = ""
    icsd_db["Ionic_Charges"] = icsd_db["Ionic_Charges"].astype('object')
    icsd_db["Unique_Ions"] = ""
    icsd_db["Unique_Ions"] = icsd_db["Unique_Ions"].astype('object')
    icsd_db["Ionic_Formula"] = ""

    ion_frequencies = [0] * n_possible_ions

    #debugging counters
    deuterium_count = 0
    multiplicity_count = 0
    no_arrangements_count = 0
    too_many_combinations_count = 0
    max_element_count = 100

    print_frequency = 10000
    counter = 0

    for index, series in icsd_db.iterrows():
        counter += 1
        if not counter % print_frequency:
            print("Compound # {}\nCounts so far:\ndeuterium:{}, multiplicity:{}, no_arrangements:{}, too_many_combinations:{}".format(counter, deuterium_count, multiplicity_count, no_arrangements_count, too_many_combinations_count ))

        original_compound_string = series["Chemical_Formula"]
        parsed_compound_array = np.array(re.findall(r'([A-z]+)(\d*\.?\d+)', original_compound_string)).transpose()

        n_unique_elements = len(parsed_compound_array[0])
        unique_elements = list(parsed_compound_array[0])
        element_counts = list(float(freq) for freq in parsed_compound_array[1])
        icsd_db.at[index, "N_Unique_Elements"] = n_unique_elements
        icsd_db.at[index, "Unique_Elements"] = unique_elements
        icsd_db.at[index, "Element_Counts"] = element_counts
        icsd_db.at[index, "Element_Count_Dictionary"] = dict(zip(unique_elements, element_counts))

        # Skip problematic compounds
        if len(icsd_db.at[index, "Unique_Elements"]) > max_element_count:
            too_many_combinations_count += 1
            continue

        if "D" in icsd_db.at[index, "Unique_Elements"] or "T" in icsd_db.at[index, "Unique_Elements"]:
            deuterium_count +=1
            continue

        best_charge_indices = find_best_charge_indices(unique_elements, element_counts)
        if len(best_charge_indices) > 1:
            multiplicity_count += 1
            continue
        elif len(best_charge_indices) == 0:
            no_arrangements_count +=1
            continue

        my_charge_indices = best_charge_indices[0]

        ion_charges = []
        element_charge_dictionary = {}
        ion_names = []
        ion_formula_components = []
        for i in range(n_unique_elements):
            element = unique_elements[i]
            element_count = element_counts[i]
            ion_charge = ionic_charge_dict[element][my_charge_indices[i]]
            ion_charges.append(ion_charge)
            element_charge_dictionary[element] = ion_charge
            ion_name = name_ion(element, ion_charge)
            ion_frequencies[possible_ion_index_dict[ion_name]] += 1
            ion_names.append(ion_name)
            ion_formula_components.append(str(element_count) + " " + ion_name)
        ionic_formula = " ".join(ion_formula_components)
        icsd_db.at[index, "Ionic_Charges"] = ion_charges
        icsd_db.at[index, "Element_Charge_Dictionary"] = element_charge_dictionary
        icsd_db.at[index, "Unique_Ions"] = ion_names
        icsd_db.at[index, "Ionic_Formula"] = ionic_formula

    print("Total ionic compounds with identified ion make-ups = {}".format(len(icsd_db[icsd_db.Ionic_Charges != ""])))

    observed_ions = np.array(possible_ions)[np.array(ion_frequencies) != 0]
    unobserved_ions = np.array(possible_ions)[np.array(ion_frequencies) == 0]

    print("Of all {} possible ions, {} were observed in the database while {} were not observed.".format(n_possible_ions, len(observed_ions), len(unobserved_ions)))

    observed_ion_index_dict = {ion:ion_index for ion_index, ion in enumerate(observed_ions)}
    icsd_db_translated = icsd_db[icsd_db.Ionic_Charges != ""]
    n_compounds_translated = len(icsd_db_translated)
    n_observed_ions = len(observed_ions)

    min_unique_elements = max(2, min(icsd_db_translated.N_Unique_Elements))
    max_unique_elements = max(icsd_db_translated.N_Unique_Elements)

    similar_env_count_matrix = np.zeros((n_observed_ions, n_observed_ions))

    # Build substitution matrix
    for n in range(min_unique_elements, max_unique_elements + 1):
        observed_ions_by_compound_and_n = icsd_db_translated[icsd_db_translated.N_Unique_Elements == n].Unique_Ions
        sorted_icsd_indices_n = np.sort(icsd_db_translated[icsd_db_translated.N_Unique_Elements == n].index)
        print("Starting n = {} which includes {} compounds".format(n, len(observed_ions_by_compound_and_n)))

        ion_bool_matrix_n = np.zeros((n_icsd_compounds, n_observed_ions))

        for i in sorted_icsd_indices_n:
            ions_i = observed_ions_by_compound_and_n[i]

            for ion in ions_i:
                ion_bool_matrix_n[i, observed_ion_index_dict[ion]] += 1

            for j in sorted_icsd_indices_n[sorted_icsd_indices_n < i]:
                ions_j = observed_ions_by_compound_and_n[j]

                diff_vector_i_j = np.abs(ion_bool_matrix_n[i] - ion_bool_matrix_n[j])
                if np.sum(diff_vector_i_j)/2 == 1:
                    similar_ion_pair_indices = tuple(diff_vector_i_j.nonzero()[0])
                    if len(similar_ion_pair_indices) != 2:
                        print("Something went wrong. More than two similar ions in compound pair")
                    similar_env_count_matrix[similar_ion_pair_indices] +=1
                    similar_env_count_matrix[similar_ion_pair_indices[::-1]] +=1

    # Save matrix to substitution_result directory
    matrix_save_path = os.path.join(substitution_result_dir, "similar_env_count_matrix")
    np.save(matrix_save_path, similar_env_count_matrix)

    # Load matrix and process
    similar_env_count_matrix_loaded = np.load(matrix_save_path + ".npy")
    nonzero_ion_frequencies = np.array(ion_frequencies)[np.nonzero(ion_frequencies)[0]]

    substitution_matrix = similar_env_count_matrix_loaded.copy()
    substitution_matrix /= nonzero_ion_frequencies
    substitution_matrix = substitution_matrix.T
    substitution_matrix /= nonzero_ion_frequencies
    substitution_matrix += np.eye(n_observed_ions)
    substitution_matrix /= np.mean(substitution_matrix)
    nonzero_min_substitution = np.min(np.ma.masked_equal(substitution_matrix, 0.0, copy=True))
    substitution_matrix += nonzero_min_substitution
    substitution_matrix = np.log10(substitution_matrix)

    # Save plot to substitution_result directory
    fig, ax = plt.subplots(figsize=(20,20))
    plt.colorbar(ax.imshow(substitution_matrix))
    print_spacing = 3
    ax.set_xticks(range(0,250, print_spacing))
    ax.set_xticklabels(observed_ions[::print_spacing], rotation=90)
    ax.set_yticks(range(0, 250, print_spacing))
    ax.set_yticklabels(observed_ions[::print_spacing])
    plot_save_path = os.path.join(substitution_result_dir, "substitution_matrix_alpha.png")
    plt.savefig(plot_save_path)
    plt.close()

    print(f"Semiconductors workflow completed. Results saved to {substitution_result_dir}")

def main():
    parser = argparse.ArgumentParser(description="Run Semiconductors Workflow")
    parser.add_argument('--dataset', required=True, help="Path to dataset directory")
    parser.add_argument('--output', required=True, help="Path to project directory for outputs")

    args = parser.parse_args()

    Semiconductors_Workflow(args.dataset, args.output)

if __name__ == "__main__":
    main()
