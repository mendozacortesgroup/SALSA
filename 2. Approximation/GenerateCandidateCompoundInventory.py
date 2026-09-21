#!/usr/bin/env python3

import pandas as pd
import chemparse
import numpy as np
import re
import itertools
import collections
import os
import argparse
from fractions import Fraction
from datetime import datetime as dt
from IPython.display import display

def GenerateCandidateCompoundInventory(dataset_dir, project_dir):
    """
    Generate candidate compound inventory for SALSA workflow
    """

    # Use existing Approximation_result directory (created during project setup)
    approximation_result_dir = os.path.join(project_dir, "Approximation_result")

    # From stackoverflow.com/a/2158532
    # Useful if you have iterables of different shapes and types (otherwise you can use np.concatenate)
    def flatten_to_generator( nested_list ):
        for element in nested_list:
            if isinstance(element, collections.abc.Iterable) and not isinstance(element, (str, bytes)):
                yield from flatten_to_generator(element)
            else:
                yield element

    # Useful if you have iterables of different shapes and types (otherwise you can use np.concatenate)
    def flatten_to_list( *nested_list ):
        return(list(flatten_to_generator(nested_list)))


    # Assumes string format of "key1:value1::key2:value2::key3:value3"
    def dict_to_str(a_dict):
        dict_string = ""
        for key, value in a_dict.items():
            dict_string += "{}:{}::".format(key, value)
        return dict_string[:-2]

    # Assumes string format of "key1:value1::key2:value2::key3:value3"
    def str_to_dict(dict_string, value_type="int"):
        new_list = dict_string.split("::")
        new_dict = {}
        for pair in new_list:
            key, value = pair.split(":")
            if value_type == "int":
                value = int(value)
            new_dict[key] = value
        return new_dict

    def parse_compound_names(compound_names):
        elements_used = []
        compound_dicts = []
        for comp in compound_names:
            comp_dict=chemparse.parse_formula(comp)
            compound_dicts.append(comp_dict)
            for el, count in comp_dict.items():
                if el not in elements_used:
                    elements_used.append(el)
        return elements_used, compound_dicts

    def compile_frequent_charges_by_element(reference_path, H_index = 1):
        with open(reference_path, "r") as file:
            lines = file.readlines()

        atomic_symbol_dict = {}
        element_charge_dict  = {}
        for line in lines:
            line = line.split()
            if line:
                atomic_symbol = line[0]
                atomic_number = int(line[1])
                ionic_charges = line[2:]
                atomic_symbol_dict[ atomic_number ] = atomic_symbol
                element_charge_dict[ atomic_symbol ]  = tuple(int(q) for q in ionic_charges)
        n_elements = len(atomic_symbol_dict.keys())
        n_element_charge_pairs = sum(len(val) for val in element_charge_dict.values())

        all_reasonable_ions = []
        for element, charges in element_charge_dict.items():
            all_reasonable_ions += [name_ion(element, charge) for charge in np.sort(charges)]
        n_reasonable_ions = len(all_reasonable_ions)
        reasonable_ion_index_dict         = {ion                 :ion_index + H_index for ion_index, ion in enumerate(all_reasonable_ions)}
        reasonable_ion_index_dict_reverse = {ion_index + H_index :ion                 for ion_index, ion in enumerate(all_reasonable_ions)}
        return n_reasonable_ions, element_charge_dict, reasonable_ion_index_dict, reasonable_ion_index_dict_reverse

    def name_ion(element, charge):
        sign_string = "+" if charge >= 0 else "-"
        ion_name = element + str(abs(charge)) + sign_string
        return ion_name

    def ionic_charge_assignment_cost_function(elements, element_counts, charge_indices):
        charge_indices = np.array(charge_indices)

        ionic_charges = [element_charge_dict[elements[i]][charge_indices[i]] for i in range(len(elements))]
        net_abs_charge = abs(np.dot(element_counts, ionic_charges))
        cost_per_electron = 1000
        priority_cost = charge_indices.sum() #+ 0*(charge_indices.sum()**2 - charge_indices.dot(charge_indices)) / 2
        cost = net_abs_charge * cost_per_electron + priority_cost
        return net_abs_charge, cost


    def find_best_charge_indices(elements, element_counts, leniency=1):

        charge_options_counts = [len(element_charge_dict[element]) for element in elements]

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
            current_charge, current_cost  = ionic_charge_assignment_cost_function(elements, element_counts, charge_indices_candidate)
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
                            if np.count_nonzero(potential_charge_indices_candidate) <= leniency:
                                charge_indices_candidates.append(potential_charge_indices_candidate)

        if best_cost < max_cost and best_charge < max_charge:
            return best_charge_indices_candidates
        else:
            return []

    def find_best_charges(unique_elements, element_counts, max_charge_options, leniency=0):
        # Find reasonable charge index assignments, listed in order of likelihood
        charge_index_possibilities = find_best_charge_indices(unique_elements, element_counts, leniency=leniency)

        if len (charge_index_possibilities) == 0:
            return []

        charge_possibilities = []
        for charge_index_possibility in charge_index_possibilities:
            charge_possibility = [element_charge_dict[unique_elements[i]][charge_index_possibility[i]] for i in range(len(unique_elements))]
            charge_possibilities.append(charge_possibility)
        return charge_possibilities[:max_charge_options]

    def parse_ion_list(ion_list):
        return np.array(re.findall(r'([A-z]+)(\d*\.?\d+)(\+|\-)', " ".join(ion_list))).T

    # vulnerable to "0+"
    def ion_charges_not_same_sign(ion_list):
        parsed_ion_list = parse_ion_list(ion_list)
        if "0" in parsed_ion_list[1]:
            print("At least one ion has charge 0 ... this function is not designed to account for that.")
        if "+" in parsed_ion_list[2] and "-" in parsed_ion_list[2]:
            return True
        else:
            return False

    # Using a systematic naming process should make searching easier
    def argsort_parsed_ionic_array(parsed_ionic_array):
        a = parsed_ionic_array
        b = np.array([a[0], [-1*int(chg) for chg in a[1]], a[2]])
        c = np.array([tuple([*d]) for d in np.array(b).T],  dtype=[('el', 'U10'), ('charge', 'i8'), ('charge_sign', 'U10')])
        ordering = np.argsort(c, order=['charge_sign', 'charge','el'], axis=0)
        return ordering

    def standardize_compound_names(compound_names, verbosity=1):
        if type(compound_names) == str:
            compound_names = [compound_names]
        _, parsed_name_dicts = parse_compound_names(compound_names)
        standardized_names = []
        for parsed_name_dict in parsed_name_dicts:
            cant_standardize=False
            elements = list(parsed_name_dict.keys())
            counts = list(parsed_name_dict.values())
            best_charge_indices = find_best_charge_indices( elements, counts, leniency=2)
            if len(best_charge_indices) > 1:
                if verbosity >= 1:
                    print("Found multiple equally good options. Using the first one listed, arbitrarily.")
            elif len(best_charge_indices) == 0:
                if verbosity >= 1:
                    print("Found no acceptable charge arrangements. Can't standardize this name.")
                cant_standardize=True
            if cant_standardize:
                standardized_name = "???"
            else:
                relative_charge_indices = best_charge_indices[0]
                ion_charges =  [element_charge_dict[elements[i]][relative_charge_indices[i]] for i in range(len(elements))]
                ion_names = [name_ion(elements[i], ion_charges[i]) for i in range(len(elements))]
                parsed_ionic_array = np.array(re.findall(r'([A-z]+)(\d*\.?\d+)(\+|\-)', " ".join(ion_names))).T
                naming_order = argsort_parsed_ionic_array(parsed_ionic_array)
                ordered_elements    = np.array(elements)[naming_order]
                ordered_composition = np.array(np.array(counts, dtype=int),dtype=str)[naming_order]
                filtered_ordered_composition = np.where(ordered_composition == "1", "",ordered_composition)
                standardized_name = "".join(np.array([ordered_elements, filtered_ordered_composition]).T.flatten())
            standardized_names.append(standardized_name)
        return standardized_names

    def ion_list_to_index_list(ion_list):
        return [reasonable_ion_index_dict[ion] for ion in ion_list]

    def ion_index_list_from_df(df):
        all_unique_ions = df.Unique_Ions
        all_unique_indices = []
        for ions in all_unique_ions:
            all_unique_indices.append(ion_list_to_index_list(ions))
        return all_unique_indices


    def index_list_to_bool_vector(index_list, true_bools = False):
        if true_bools:
            ion_bool_vector = np.zeros(len(reasonable_ions)+1, dtype=bool)
        else:
            ion_bool_vector = np.zeros(len(reasonable_ions)+1, dtype=int)
        for ion_index in index_list:
            ion_bool_vector[ion_index] += 1
        return ion_bool_vector

    def ion_list_to_bool_vector(ion_list, true_bools = False):
        return(index_list_to_bool_vector(ion_list_to_index_list(ion_list), true_bools))

    def index_list_to_ion_list(index_list):
        ion_list = [reasonable_ion_index_dict_reverse[index] for index in index_list]
        return ion_list

    def bool_vector_to_ion_list(bool_vec):
        bool_vec= np.array(bool_vec)
        ion_indices = np.where(bool_vec>0)[0]
        return index_list_to_ion_list(ion_indices)


    def ion_bool_vectors_from_df(df, true_bools = False):
        init_ion_index_pairs = ion_index_list_from_df(df)
        ion_bool_vectors = []
        for ion_index_pair in init_ion_index_pairs:
            if len(ion_index_pair):
                ion_bool_vectors.append(index_list_to_bool_vector(ion_index_pair, true_bools))
        return np.array(ion_bool_vectors)

    def create_descriptive_compound_dataframe(df, additional_columns=()):
        if type(df) == list or type(df) == np.ndarray:
            df = pd.DataFrame(df, columns = ["Compound"])

        columns = ["N_Unique_Elements", "Chemical_Formula", "Unique_Elements", "Element_Counts", "N_Atoms", "Composition_Dictionary",
                "Element_Charge_Dictionary", "Ionic_Charges", "Unique_Ions", "Ionic_Formula"]
        # These columns are of type list or dictionary
        collection_columns = ["Unique_Elements", "Element_Counts", "Composition_Dictionary", "Element_Charge_Dictionary", "Ionic_Charges",
                            "Unique_Ions"]

        for col in additional_columns:
            columns.append(col)

        for col in columns:
            if col not in df.columns:
                df[col] = ""
            if col in collection_columns:
                df[col] = df[col].astype('object')

        for index, series in df.iterrows():
            # Non-ionic compositional information
            unique_elements, [composition_dictionary] = parse_compound_names([series["Compound"]])
            n_unique_elements = len(unique_elements)
            element_counts = list(composition_dictionary.values())
            df.at[index, "N_Unique_Elements"]         = n_unique_elements
            df.at[index, "Unique_Elements"]           = unique_elements
            df.at[index, "Element_Counts"]            = element_counts
            df.at[index, "N_Atoms"]                   = int(sum(element_counts))
            df.at[index, "Composition_Dictionary"]    = composition_dictionary
            df.at[index, "Chemical_Formula"]          = " ".join([key + str(int(value)) for key, value in composition_dictionary.items()])

            # Ionic compositional information
            most_likely_charges = find_best_charges(unique_elements, element_counts, max_charge_options=1, leniency=2)

            if len(most_likely_charges) == 0:
                print("Could not determine ionic charges for {}".format(series["Compound"]))
                continue
            else:
                most_likely_charges = most_likely_charges[0]

            element_charge_dictionary = {element:charge for element, charge in zip(unique_elements, most_likely_charges)}
            unique_ions = [ name_ion(element, charge) for element, charge in zip(unique_elements, most_likely_charges)]
            df.at[index, "Ionic_Charges"] = most_likely_charges
            df.at[index, "Element_Charge_Dictionary"] = element_charge_dictionary
            df.at[index, "Unique_Ions"] = unique_ions
            df.at[index, "Ionic_Formula"] = " ".join([str(count) + " " + ion for count, ion in zip(element_counts, unique_ions)])

        return df

    def compile_substitution_correlations(reference_path):
        with open(reference_path, 'r') as file:
            lines = file.readlines()

        correlation_dictionary = {}
        for line in lines[1:]:
            correlation = float(re.findall(r'\d+\.?\d*', line.split(" ")[-1])[0])
            correlation_label = ' '.join(line.split(" ")[:-1])
            correlation_dictionary[correlation_label] = correlation
        return correlation_dictionary

    def filter_pairings_by_substitutability(pairs, df, substitution_dictionary, correlation_threshold = 1, verbose=True):
        n_pairings_allowed_by_sub_matrix = 0
        allowed_pairs = []
        for comps in pairs:
            ion_lists = []
            for i in range(2):
                element_charge_dictionary = df[df["Compound"] == comps[i]]["Element_Charge_Dictionary"].to_numpy()[0]
                if type(element_charge_dictionary) == dict:
                    ion_list = []
                    for element, charge in element_charge_dictionary.items():
                        charge_sign = "+" if charge > 0 else "-"
                        ion_string = element + ":" + str(np.abs(charge)) + charge_sign
                        ion_list.append(ion_string)
                    ion_lists.append(ion_list)
                else:
                    ion_lists.append([])

            comp_partners_forbidden_by_sub_matrix = False
            for i in range(2):
                ion_list_1 = ion_lists[i]
                ion_list_2 = ion_lists[i-1]

                for ion_1 in ion_list_1:
                    sub_partner_absent_in_other_compound = True
                    if ion_1 in ion_list_2:
                        sub_partner_absent_in_other_compound = False
                    for ion_2 in ion_list_2:
                        if sub_partner_absent_in_other_compound:
                            corr_1_2 = substitution_dictionary["{} {}".format(ion_1, ion_2)]
                            if corr_1_2 >= correlation_threshold:
                                sub_partner_absent_in_other_compound = False
                    if sub_partner_absent_in_other_compound:
                        comp_partners_forbidden_by_sub_matrix = True
            if comp_partners_forbidden_by_sub_matrix and verbose:
                print("No interpolation between {} and {}".format(*comps))
            else:
                n_pairings_allowed_by_sub_matrix += 1
                allowed_pairs.append(comps)
        if verbose:
            print( "{} / {} interpolation pairings are allowed by substitution probability matrix".format(n_pairings_allowed_by_sub_matrix, len(pairs)))
        return allowed_pairs

    def count_quaternary_constrained_max_atoms(n_atoms_a, n_atoms_b, N_atoms_max, verbosity=2):
        if verbosity > 1:
            print("Counting the number of unique quaternary compounds that can be formed with {}-atom and {}-atom compounds".format(n_atoms_a, n_atoms_b), end=" ")
            print("under the constraint of {} or fewer total atoms in a unit cell.".format(N_atoms_max))
        simplified_pairs = []
        for i, j in itertools.product(range(1,N_atoms_max+1), repeat=2):
            CD_ij = np.gcd(i, j)
            i_simplified = int(i/CD_ij)
            j_simplified = int(j/CD_ij)
            n_atoms_i = i_simplified* n_atoms_a
            n_atoms_j = j_simplified* n_atoms_b
            n_atoms_total = n_atoms_i + n_atoms_j
            if n_atoms_total <= N_atoms_max:
                if (i_simplified, j_simplified) not in simplified_pairs:
                    simplified_pairs.append((i_simplified, j_simplified))
                    frac_i = n_atoms_i/ n_atoms_total
                    frac_j = n_atoms_j/ n_atoms_total
                    if verbosity > 2:
                        pair_ij_printout = "{} A + {} B => A{}B{} which has {} atoms total".format(i, j, i_simplified, j_simplified, n_atoms_total)
                        if verbosity > 3:
                            pair_ij_printout += "   ({} A   {} B)".format(Fraction(frac_i).limit_denominator(), Fraction(frac_j).limit_denominator())
                        print(pair_ij_printout)
        if verbosity > 1:
            print("{} total unique compounds formable.".format(len(simplified_pairs)))
        if verbosity == 1:
            print("({},{})| {}  => {}".format(n_atoms_a, n_atoms_b, N_atoms_max,len(simplified_pairs)))
        return simplified_pairs

    def count_ternary_constrained_max_atoms(composition_a, composition_b, N_atoms_max, verbosity=2, return_fracs=False):
        if len(composition_a) != 2 or len(composition_b) != 2:
            print("Cannot verify that both compositions were binary. Cannot move forward with count.")
            return 1
        n_atoms_a = np.sum(composition_a)
        n_atoms_b = np.sum(composition_b)
        CD_a = np.gcd(*np.array(composition_a, dtype=int))
        CD_b = np.gcd(*np.array(composition_b, dtype=int))
        composition_a = np.array(composition_a, dtype=int) // CD_a
        composition_b = np.array(composition_b, dtype=int) // CD_b
        simplified_triads = []
        fracs = []
        if verbosity > 1:
            print("Counting the number of unique ternary compounds that can be formed with compounds of the form C{}A{} and C{}B{}".format(*composition_a, *composition_b), end=" ")
            print("under the constraint of {} or fewer total atoms in a unit cell.".format(N_atoms_max))

        for i, j in itertools.product(range(1,N_atoms_max+1), repeat=2):
            CD_ij = np.gcd(i, j) # Common Divisor
            i_simplified = int(i/CD_ij)
            j_simplified = int(j/CD_ij)
            new_composition = np.array([ i_simplified * composition_a[0] + j_simplified * composition_b[0], i_simplified * composition_a[1], j_simplified * composition_b[1]], dtype=int)
            CD_ternary = np.gcd(new_composition[0], np.gcd(*new_composition[1:]))
            new_composition //= CD_ternary
            n_atoms_total = np.sum(new_composition)
            if n_atoms_total <= N_atoms_max:
                if tuple(new_composition) not in simplified_triads:
                    simplified_triads.append(tuple(new_composition))
                    frac_a = n_atoms_a*i_simplified/CD_ternary/ n_atoms_total
                    fracs.append(frac_a)
                    frac_b = n_atoms_b*j_simplified/CD_ternary/ n_atoms_total
                    Frac_a = Fraction( frac_a ).limit_denominator()
                    numer_a = Frac_a.numerator
                    denom_a = Frac_a.denominator
                    Frac_b = Fraction(frac_b ).limit_denominator()
                    numer_b = Frac_b.numerator
                    denom_b = Frac_b.denominator
                    if verbosity > 2:
                        new_comp_string = "C{}A{}B{}".format(*new_composition)
                        ternary_comp_printout = "{:2d} C{}A{} + {:2d} C{}B{} => {:9s} which has {:2d} atoms total".format(i, *composition_a, j, *composition_b, new_comp_string, n_atoms_total)
                        if verbosity > 3:
                            ternary_comp_printout += "   ({:2d}/{:2d} C{}A{}  {:2d}/{:2d} C{}B{})".format(numer_a, denom_a, *composition_a, numer_b, denom_b, *composition_b)
                        print(ternary_comp_printout)
        if verbosity > 1:
            print("{} total unique compounds formable.".format(len(simplified_triads)))
        if verbosity == 1:
            print("(({},{}),({},{}))| {}  => {}".format(*composition_a, *composition_b, N_atoms_max, len(simplified_triads)))
        if return_fracs:
            return fracs, simplified_triads
        else:
            return simplified_triads

    def populate_interp_comp_df(indices, composition, parent_compound_df, parent_a_frac):
        interp_comp_info = {}
        ions = [reasonable_ion_index_dict_reverse[index] for index in indices]
        parsed_ionic_array = np.array(re.findall(r'([A-z]+)(\d*\.?\d+)(\+|\-)', " ".join(ions))).T
        name_ordering = argsort_parsed_ionic_array(parsed_ionic_array)
        ions = list(np.array(ions)[name_ordering])
        composition = list(np.array(composition)[name_ordering])
        elements = list(parsed_ionic_array[0, name_ordering])
        element_info = np.array([elements, composition], dtype=str)
        charges = [int(charge+sign) for sign, charge in parsed_ionic_array[1:,name_ordering].T]
        parent_b_frac = 1 - parent_a_frac

        compound_name = "".join(np.where(element_info == "1", "",element_info).T.flatten())
        standardized_name = standardize_compound_names(compound_name, verbosity=0)[0]
        interp_comp_info["Compound"] = standardized_name if standardized_name != "???" else compound_name
        interp_comp_info["Chemical_Formula"] = " ".join(["".join(element_i_info) for element_i_info in element_info.T])
        interp_comp_info["N_Unique_Elements"] = len(elements)
        interp_comp_info["Unique_Elements"] = elements
        interp_comp_info["Element_Counts"] = list(composition)
        interp_comp_info["Composition_Dictionary"] = { element_i_info[0]:int(element_i_info[1]) for element_i_info in element_info.T}
        interp_comp_info["Element_Charge_Dictionary"] = {elements[i]:charges[i] for i in range(len(elements))}
        interp_comp_info["Ionic_Charges"] = charges
        interp_comp_info["Unique_Ions"] = ions
        interp_comp_info["N_Atoms"] = np.sum(composition)
        interp_comp_info["Parents"] = list(parent_compound_df["Compound"].to_numpy())
        interp_comp_info["Parental_Makeup_Dictionary"] = {interp_comp_info["Parents"][0]:parent_a_frac, interp_comp_info["Parents"][1]:parent_b_frac}
        interp_comp_info["Ionic_Formula"] = " ".join([str(count) + " " + ion for count, ion in zip(list(composition), ions)])
        return interp_comp_info

    def find_possible_interpolated_candidate_compounds(interpolation_pairs, parent_df, max_size, return_duplicates=False):
        hybrid_df = create_descriptive_compound_dataframe([], ["Parents", "Parental_Makeup_Dictionary"])
        all_interpolated_compositions = []
        duplicate_compositions = []
        for comp_names in interpolation_pairs:
            compounds = parent_df[parent_df["Compound"].isin(comp_names)].copy()
            all_unique_indices = ion_index_list_from_df(compounds)
            n_unique_ions = len(np.unique(flatten_to_list(all_unique_indices)))
            all_species_counts = compounds.Element_Counts.to_numpy()
            two_binary_compounds = np.all([len(indices)==2 for indices in all_unique_indices])
            if two_binary_compounds:

                if n_unique_ions == 3:  # <---------------- ternary
                    index_in_common = np.intersect1d(*all_unique_indices)[0]
                    index_in_common_index1 = np.where(all_unique_indices[0] == index_in_common)[0][0]
                    index_in_common_index2 = np.where(all_unique_indices[1] == index_in_common)[0][0]
                    compositions = [(all_species_counts[0][index_in_common_index1], all_species_counts[0][index_in_common_index1-1]), (all_species_counts[1][index_in_common_index2], all_species_counts[1][index_in_common_index2-1])]
                    fracs_for_interpolation, interpolated_compositions = count_ternary_constrained_max_atoms(*compositions, max_size, verbosity=0, return_fracs=True)
                    indices_for_interpolated_compositions = [index_in_common, all_unique_indices[0][index_in_common_index1-1], all_unique_indices[1][index_in_common_index2-1]]
                    for j in range(len(interpolated_compositions)):
                        interpolated_composition = interpolated_compositions[j]
                        interpolated_ion_vector = np.zeros(n_reasonable_ions+1, dtype=int)
                        for i in range(3):
                            interpolated_ion_vector[indices_for_interpolated_compositions[i]] += interpolated_composition[i]
                        if tuple(interpolated_ion_vector) not in all_interpolated_compositions:
                            all_interpolated_compositions.append(tuple(interpolated_ion_vector))
                            interp_comp_dict = populate_interp_comp_df(indices_for_interpolated_compositions, interpolated_composition, compounds, fracs_for_interpolation[j])
                            hybrid_df = pd.concat([hybrid_df, pd.DataFrame(np.array([interp_comp_dict[key] for key in hybrid_df.columns], dtype=object).reshape(1,-1), columns=hybrid_df.columns)], ignore_index=True)
                            if not(len(all_interpolated_compositions) % 1000):
                                print("Status Update: {:5d} hybrid compounds found so far".format(len(all_interpolated_compositions)))
                        else:
                            if return_duplicates:
                                duplicate_compositions.append(tuple(interpolated_ion_vector))
                                print("Repeat found", indices_for_interpolated_compositions, interpolated_composition)

                elif n_unique_ions == 4:  # <------------------- quaternary
                    all_species_counts_np=np.array([[*counts] for counts in all_species_counts], dtype=int).reshape([2,2])
                    n_atoms_per_binary_comp = np.sum([[*counts] for counts in all_species_counts], axis=1, dtype=int)
                    composition_ratios = count_quaternary_constrained_max_atoms(*n_atoms_per_binary_comp,  max_size, verbosity=0)
                    indices_for_interpolated_compositions = np.array(all_unique_indices).flatten()
                    for composition_ratio in composition_ratios:
                        interpolated_composition = (all_species_counts_np*np.array(composition_ratio).reshape(2,1)).flatten()
                        interpolated_ion_vector = np.zeros(n_reasonable_ions+1, dtype=int)
                        for i in range(4):
                            interpolated_ion_vector[indices_for_interpolated_compositions[i]] += interpolated_composition[i]
                        if tuple(interpolated_ion_vector) not in all_interpolated_compositions:
                            all_interpolated_compositions.append(tuple(interpolated_ion_vector))
                            interp_comp_dict = populate_interp_comp_df(indices_for_interpolated_compositions, interpolated_composition, compounds, parent_a_frac=composition_ratio[0]/sum(composition_ratio))
                            hybrid_df = pd.concat([hybrid_df, pd.DataFrame(np.array([interp_comp_dict[key] for key in hybrid_df.columns], dtype=object).reshape(1,-1), columns=hybrid_df.columns)], ignore_index=True)
                            if not(len(all_interpolated_compositions) % 1000):
                                print("Status Update: {:5d} hybrid compounds found so far".format(len(all_interpolated_compositions)))
                        else:
                            if return_duplicates:
                                duplicate_compositions.append(tuple(interpolated_ion_vector))
                                print("Repeat found", indices_for_interpolated_compositions, interpolated_composition)
                else:
                    print("Did not expect interpolation pairs with {} unique ions, but this pair seems to have that many: {} --- {}".format(n_unique_ions, *comp_names))
        if return_duplicates:
            return hybrid_df, all_interpolated_compositions, duplicate_compositions
        else:
            return hybrid_df, all_interpolated_compositions

    def assign_hybrid_properties(parent_df, hybrid_df, props):
        for prop in props:
            if prop not in hybrid_df.columns:
                hybrid_df[prop] = ""

        for index, series in hybrid_df.iterrows():
            parent_1, parent_2 = series["Parents"]
            parent_1_frac, parent_2_frac = series["Parental_Makeup_Dictionary"][parent_1], series["Parental_Makeup_Dictionary"][parent_2]
            parent_prop_values = parent_df[parent_df["Compound"].isin(series["Parents"])][props].to_numpy()
            prop_values = parent_prop_values[0] * parent_1_frac + parent_prop_values[1] * parent_2_frac
            for i in range(len(props)):
                hybrid_df.at[index, props[i]] = prop_values[i]

        return hybrid_df

    def compile_property_space_restrictions(df, include_tolerance=False):
        restriction_list = []
        for index, series in df.iterrows():
            property_name = series["Property"]
            if series["Minimum"] != "":
                minimum = series["Minimum"]
                if include_tolerance and series["Tolerance"] != "":
                    minimum = float(minimum) - float(series["Tolerance"])
                restriction_list.append("{} >= {}".format(property_name, minimum))
            if series["Maximum"] != "":
                maximum = series["Maximum"]
                if include_tolerance and series["Tolerance"] != "":
                    maximum = float(maximum) + float(series["Tolerance"])
                restriction_list.append("{} <= {}".format(property_name, maximum))
        restriction_string = " & ".join(restriction_list)
        return restriction_string

    def export_compounds_for_SALSA_on_hpcc(df, example_inventory_path, prop_df):
        inventory_cols = pd.read_csv(example_inventory_path).columns.to_numpy()
        out_df = pd.DataFrame(df.Compound.to_numpy(), columns = ["compound_name"])
        for col in inventory_cols:
            if col not in out_df.columns:
                out_df[col] = ""
        out_in_index_dict = {out_i:in_i for in_i, out_i in zip(df.index.to_list(), out_df.index.to_list())}

        props = prop_df["Property"].to_list()
        ideal_conditions = compile_property_space_restrictions(prop_df, include_tolerance=False)
        ideal_df = df.query(ideal_conditions)
        target_conditions = compile_property_space_restrictions(prop_df, include_tolerance=True)
        target_df = df.query(target_conditions)
        for index, series in out_df.iterrows():
            in_index = out_in_index_dict[index]
            out_df.at[index, "composition_dictionary"] = dict_to_str(df.at[in_index, "Composition_Dictionary"])
            out_df.at[index, "source"] = "new_candidate_compound"
            out_df.at[index, "interpolated_property_dictionary"] = dict_to_str({prop:"{:.3f}".format(df.at[in_index, prop]) for prop in props})
            parent_1, parent_2 = df.at[in_index, "Parents"]
            out_df.at[index, "parent_1_compound_name"] = parent_1
            out_df.at[index, "parent_2_compound_name"] = parent_2
            out_df.at[index, "parent_1_fraction"] = df.at[in_index, "Parental_Makeup_Dictionary"][parent_1]
            out_df.at[index, "parent_2_fraction"] = df.at[in_index, "Parental_Makeup_Dictionary"][parent_2]
            out_df.at[index, "interpolated_properties_in_ideal_region"] = 1 if (in_index in ideal_df.index.to_list()) else 0
            out_df.at[index, "interpolated_properties_in_target_region"] = 1 if (in_index in target_df.index.to_list()) else 0

        return out_df

    # Main processing starts here - using paths from arguments
    initial_compounds_path = os.path.join(dataset_dir, "known_compositions_and_properties.csv")
    frequent_charges_by_element_path = os.path.join(dataset_dir, "Element_Charge_Dictionary.txt")
    substitution_correlations_path = os.path.join(project_dir, "substitution_result", "substitutions_textual.txt")
    property_space_parameters_path = os.path.join(dataset_dir, "property_space_selection_criteria.csv")
    inventory_template_path = os.path.join(project_dir, "template_inventory.csv")
    hpcc_inventory_output_path = os.path.join(approximation_result_dir, "new_candidate_compound_inventory.csv")

    #Import names and properties of known (initial) compounds
    initial_compounds_df = pd.read_csv(initial_compounds_path)
    initial_compound_name_list = initial_compounds_df["Compound"].to_numpy()
    print("Loaded {} initial compounds".format(len(initial_compounds_df)))

    #Generate elemental and ionic composition information for initial compounds
    n_reasonable_ions, element_charge_dict, reasonable_ion_index_dict, reasonable_ion_index_dict_reverse = compile_frequent_charges_by_element(frequent_charges_by_element_path)
    reasonable_ions = [reasonable_ion_index_dict_reverse[i] for i in range(1, n_reasonable_ions+1)]
    initial_compounds_df = create_descriptive_compound_dataframe(initial_compounds_df)

    print("Total ionic compounds with identified ion make-ups = {}".format(len(initial_compounds_df[initial_compounds_df.Ionic_Charges != ""])))

    # Determine a list of compound pairings allowed by a substitution probability matrix
    substitution_correlation_dictionary = compile_substitution_correlations(substitution_correlations_path)
    all_possible_pairings = [*itertools.combinations(initial_compound_name_list, r=2)]
    interpolation_allowed_pairings = filter_pairings_by_substitutability(all_possible_pairings, initial_compounds_df, substitution_correlation_dictionary, correlation_threshold = 1)

    max_unitcell_size=20
    candidate_compounds_df,_=find_possible_interpolated_candidate_compounds(interpolation_allowed_pairings, initial_compounds_df, max_unitcell_size, return_duplicates=False)
    print("{} candidate compounds interpolated with the constraint of {} atoms max in a unit cell".format(len(candidate_compounds_df), max_unitcell_size))

    secondary_size_limit = 4
    candidate_compounds_limited = candidate_compounds_df[candidate_compounds_df.N_Atoms <= secondary_size_limit]
    n_trinary = len(candidate_compounds_limited[candidate_compounds_limited.N_Unique_Elements == 3])
    n_quaternary = len(candidate_compounds_limited[candidate_compounds_limited.N_Unique_Elements == 4])
    print("{} candidate compounds can be interpolated in total under the constraint of {} atoms max in a unit cell".format(len(candidate_compounds_limited), secondary_size_limit))
    print("This includes {} trinary and {} quaternary compounds.".format(n_trinary, n_quaternary))
    print("The average interpolated compounds per interpolation pair was {:.0f} among the {} substitution-allowed interpolations.".format(len(candidate_compounds_limited)/len(interpolation_allowed_pairings), len(interpolation_allowed_pairings),))

    #Property space
    property_space_parameters_df = pd.read_csv(property_space_parameters_path, keep_default_na=False)
    property_list = property_space_parameters_df["Property"].to_list()
    candidate_compounds_df = assign_hybrid_properties(initial_compounds_df, candidate_compounds_df, property_list)

    # Property Restrictions
    ideal_region_conditions = compile_property_space_restrictions(property_space_parameters_df, include_tolerance=False)
    target_region_conditions = compile_property_space_restrictions(property_space_parameters_df, include_tolerance=True)
    candidate_compounds_df.query(ideal_region_conditions)
    target_candidate_compounds_df = candidate_compounds_df.query(target_region_conditions)

    #Export Candidate compounds
    # To save all the data from the dataframe used in this notebook in the same way it is represented here
    today = dt.now().strftime("%Y_%m_%d")
    saved_output_path = os.path.join(approximation_result_dir, "candidate_compounds_upto_{}_atom_generated_on_{}".format(max_unitcell_size, today))
    candidate_compounds_df.to_csv(saved_output_path+".csv", index=False)
    candidate_compounds_df.to_pickle(saved_output_path+".pickle")

    # For use with the hpc portion of SALSA (the LSA of SALSA)
    hpcc_inventory_df = export_compounds_for_SALSA_on_hpcc(target_candidate_compounds_df, inventory_template_path, property_space_parameters_df)
    hpcc_inventory_df.to_csv(hpcc_inventory_output_path, index=False)

def main():
    parser = argparse.ArgumentParser(description="Generate candidate compound inventory")
    parser.add_argument('--dataset', required=True, help="Path to dataset directory")
    parser.add_argument('--output', required=True, help="Path to project directory for outputs")

    args = parser.parse_args()

    GenerateCandidateCompoundInventory(args.dataset, args.output)

if __name__ == "__main__":
    main()
