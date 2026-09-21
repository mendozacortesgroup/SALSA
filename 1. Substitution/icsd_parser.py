#!/usr/bin/env python3

import sys
import re
import os
import argparse
from dataclasses import dataclass
from typing import List

@dataclass
class Compound:
    ions: List[int]
    band_gap: float
    oxid_potential: float
    red_potential: float

def divide_by_common_divisor(a, b):
    """Find common divisor for charges"""
    while (a % 2 == 0) and (b % 2 == 0):
        a = a //2
        b = b // 2

    while (a % 3 == 0) and (b % 3 == 0):
        a = a // 3
        b = b // 3

    while (a % 5 == 0) and (b % 5 == 0):
        a = a // 5
        b = b // 5

    while (a % 7 == 0) and (b % 7 == 0):
        a = a // 7
        b = b // 7

    return a, b

def main():
    parser = argparse.ArgumentParser(description="Parse ICSD substitutions and semiconductors data")
    parser.add_argument('--dataset', required=True, help="Path to dataset directory")
    parser.add_argument('--output', required=True, help="Path to project directory for outputs")

    args = parser.parse_args()

    # Set up paths
    dataset_dir = args.dataset
    project_dir = args.output

    # Use existing substitution_result directory (created during project setup)
    substitution_result_dir = os.path.join(project_dir, "substitution_result")

    # Input file from substitution_result (created by previous step)
    substitutions_file = os.path.join(substitution_result_dir, "substitutions_textual.txt")

    # Output files in substitution_result directory
    ions_code_file_path = os.path.join(substitution_result_dir, "ions_codes")
    ions_values_file_path = os.path.join(substitution_result_dir, "ions_values")

    # Check if substitutions file exists
    if not os.path.exists(substitutions_file):
        print(f"Error: substitutions_textual.txt not found at {substitutions_file}")
        print("Make sure to run the substitution generation script first.")
        return

    # Open input file
    try:
        with open(substitutions_file, "r", encoding="utf-8") as input_file:
            lines = input_file.readlines()
    except FileNotFoundError:
        print("Error: substitutions_textual.txt not found")
        return

    # Open output files
    ions_code_file = open(ions_code_file_path, "w", encoding="utf-8")
    ions_values_file = open(ions_values_file_path, "w", encoding="utf-8")

    ions = []
    letters_pattern = re.compile(r'[a-zA-Z]')

    # Skip first line (header)
    cur = 0

    # Process substitutions file
    for line in lines[1:]:  # Skip first line
        line = line.strip()
        if not line:
            continue

        str_parsed = line.split()
        if len(str_parsed) < 3:
            continue

        ion_first = str_parsed[0]
        ion_second = str_parsed[1]
        likelihood = float(str_parsed[2])

        # Extract charge from first ion
        ion_first_charge = ion_first
        ion_first_charge = ion_first_charge.replace("+", "")
        ion_first_charge = ion_first_charge.replace("-", "")
        ion_first_charge = ion_first_charge.replace(":", "")
        ion_first_charge = letters_pattern.sub("", ion_first_charge)

        # Extract charge from second ion
        ion_second_charge = ion_second
        ion_second_charge = ion_second_charge.replace("+", "")
        ion_second_charge = ion_second_charge.replace("-", "")
        ion_second_charge = ion_second_charge.replace(":", "")
        ion_second_charge = letters_pattern.sub("", ion_second_charge)

        try:
            first_ion_num = int(ion_first_charge) if ion_first_charge else 1
            second_ion_num = int(ion_second_charge) if ion_second_charge else 1
        except ValueError:
            first_ion_num = 1
            second_ion_num = 1

        first_ion_num, second_ion_num = divide_by_common_divisor(first_ion_num, second_ion_num)

        # Handle first ion
        if ion_first not in ions:
            cur += 1
            first_num = cur
            ions.append(ion_first)
            ions_code_file.write(f"{cur} {ion_first}\n")
        else:
            first_num = ions.index(ion_first) + 1

        # Handle second ion
        if ion_second not in ions:
            cur += 1
            second_num = cur
            ions.append(ion_second)
            ions_code_file.write(f"{cur} {ion_second}\n")
        else:
            second_num = ions.index(ion_second) + 1

        ions_values_file.write(f"{first_num} {second_num} {likelihood} {first_ion_num} {second_ion_num}\n")

    # Process semiconductors file from dataset
    semiconductors = []
    sem_file_name = os.path.join(dataset_dir, "semiconductors_existing.txt")

    # Check if there are additional command line arguments for semiconductors file variant
    if len(sys.argv) > 5:  # account for script name + 4 required args (--dataset, path, --output, path)
        variant = sys.argv[5]
        sem_file_name = os.path.join(dataset_dir, f"semiconductors_existing{variant}.txt")

    try:
        with open(sem_file_name, "r", encoding="utf-8") as semiconductors_file:
            for line in semiconductors_file:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                str_parsed = line.split()
                length = len(str_parsed)

                if length < 3:
                    continue

                try:
                    band_gap = float(str_parsed[length - 3])
                    oxid_potential = float(str_parsed[length - 2])
                    red_potential = float(str_parsed[length - 1])
                except (ValueError, IndexError):
                    continue

                help_compound = Compound([], band_gap, oxid_potential, red_potential)

                for i in range(length - 3):
                    ion_name = str_parsed[i]
                    if ion_name in ions:
                        help_compound.ions.append(ions.index(ion_name) + 1)
                    else:
                        cur += 1
                        help_compound.ions.append(cur)
                        ions.append(ion_name)
                        ions_code_file.write(f"{cur} {ion_name}\n")

                semiconductors.append(help_compound)

    except FileNotFoundError:
        print(f"Warning: {sem_file_name} not found, skipping semiconductors processing")

    # Write semiconductors values to substitution_result directory
    sem_output_file = os.path.join(substitution_result_dir, "semiconductors_existing_values.csv")
    if len(sys.argv) > 5:
        variant = sys.argv[5]
        sem_output_file = os.path.join(substitution_result_dir, f"semiconductors_existing{variant}_values.csv")

    with open(sem_output_file, "w", encoding="utf-8") as semiconductors_values_file:
        for compound in semiconductors:
            # Write ion indices
            for ion_idx in compound.ions:
                semiconductors_values_file.write(f"{ion_idx},")

            # Write properties
            semiconductors_values_file.write(f"{compound.band_gap},")
            semiconductors_values_file.write(f"{compound.oxid_potential},")
            semiconductors_values_file.write(f"{compound.red_potential}\n")

    # Close files
    ions_code_file.close()
    ions_values_file.close()

    print(f"ICSD parser completed successfully")
    print(f"Outputs written to: {substitution_result_dir}")

if __name__ == "__main__":
    main()
