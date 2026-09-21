#!/usr/bin/env python3
# This script finalized by Sean
# The primary function in this script was inspired by and initially adapted from CIF2D12 function, written by Danny and Marcus


import sys 
import argparse
import ase
from ase.io import read
import ase.spacegroup as spacegroup
import math
import os
import re

#print("ase version = {}".format(ase.__version__))

D12_template_keywords = {"Title":"CALCULATION_NAME", "Reference Geometry":"GEOMETRY_INPUT_BLOCK", "Geometry Editing":"GEOMETRY_OPTIMIZATION_BLOCK", "Basis Set":"BASIS_SET_BLOCK", "k-Points":"K_POINTS_BLOCK"}

# Inspired by and initially adapted from CIF2D12 function, written by Danny and Marcus
def extract_BULK_reference_geometry_with_sym(structure, input_file=''):
    reference_geometry_text = ""

    # GET SPACE GROUP
    if 'spacegroup' in structure.info:
        sg = structure.info['spacegroup'].no
    else:
        no_sg_warning = "Warning: When extracting reference geometry "
        if input_file:
            no_sg_warning += "from {} ".format(input_file)
        no_sg_warning += "could not identify spacegroup so 'P 1' will be used."
        print(no_sg_warning)
        sg = 1
    reference_geometry_text += str(sg) + "\n"

    # GET LATTICE PARAMETERS
    cells_params = structure.cell.get_bravais_lattice().vars()
    min_params_printout = " ".join(["{:.6f}   ".format(param) for param in cells_params.values()])
    reference_geometry_text += min_params_printout + "\n"
    
    # GET ATOM POSTIONS
    pos = spacegroup.get_basis(structure)
    n_atoms = len(pos)
    reference_geometry_text += str(n_atoms) + "\n"
    
    atomic_symbols = structure.get_chemical_symbols() # structure.info['_atom_site_type_symbol']
    atomic_nums    = structure.get_atomic_numbers()   # [atomic_numbers[symb] for symb in atomic_symbols]
    reference_geometry_text += "\n".join([ "{:<3d} {:10.6f}  {:10.6f}  {:10.6f} # {:<2s}".format(atomic_nums[i], *pos[i], atomic_symbols[i]) for i in range(n_atoms)])

    return reference_geometry_text

# Inspired by and initially adapted from CIF2D12 function, written by Danny and Marcus
def compile_basis_set(structure, basis_set_directory, input_file=''):
    basis_set_text = ''
    
    atomic_nums = structure.get_atomic_numbers()   # [atomic_numbers[symb] for symb in atomic_symbols]
    
    for atomic_num in atomic_nums:
        basis_set_file_matches = []
        for basis_set_file in os.listdir(basis_set_directory):
            basis_set_path = basis_set_directory + "/" + basis_set_file
            file_atomic_nums = re.findall(r'\d+',basis_set_file)
            if str(atomic_num) in file_atomic_nums:
                if len(file_atomic_nums) > 1:
                    unparsable_error = "Basis set file {} cannot be parsed. Cannot convert ".format(basis_set_path)
                    if input_file:
                        unparsable_error += "{} ".format(input_file)
                    unparsable_error += "to d12."
                    print(unparsable_error)
                    return
                else:
                    basis_set_file_matches.append(basis_set_path)
        if len(basis_set_file_matches) == 1:
            with open(basis_set_file_matches[0], 'r') as file:
                lines = file.readlines()
            basis_set_text += "".join(lines) 
        elif len(basis_set_file_matches) == 0:
            missing_basis_set_error = "Basis set for element {} is missing from '{}'. Cannot convert ".format(atomic_num, basis_set_directory)
            if input_file:
                missing_basis_set_error += "{} ".format(input_file)
            missing_basis_set_error += "to d12."
            print(missing_basis_set_error)
            return
        elif len(basis_set_file_matches) > 1:
            multiple_basis_set_error = "Multiple basis set files {} seem to exist for element {}. Cannot convert ".format(basis_set_file_matches, atomic_num)
            if input_file:
                multiple_basis_set_error += "{} ".format(input_file)
            multiple_basis_set_error += "to d12."
            print(multiple_basis_set_error)
            return
    
    basis_set_text += "99 0\nEND"
    
    return basis_set_text
                    
# Inspired by and initially adapted from CIF2D12 function, written by Danny and Marcus
def determine_k_point_grid(structure, k_point_criterion):
    shrink_text = ''
    
    a, b, c = structure.cell.cellpar()[:3]
    
    ka = math.ceil(k_point_criterion/a)
    kb = math.ceil(k_point_criterion/b)
    kc = math.ceil(k_point_criterion/c)
    k_ISP = 2* max(ka, kb, kc)
    
    shrink_text += 'SHRINK\n'
    shrink_text += "{:2d} {:2d}\n".format(0, k_ISP)
    shrink_text += "{:2d} {:2d} {:2d}".format(ka, kb, kc)
    
    return shrink_text
    
# Inspired by and initially adapted from CIF2D12 function, written by Danny and Marcus
# Currently using settings consistent with those used for CRYSTAL09 SALSA calculations
def instantiate_D12(input_geometry_file, output_D12, template_D12, calculation_name, basis_set_directory, D12_type="BULK", k_point_criterion=40):
    if D12_type not in ["BULK", "SLAB"]:
        print("Warning: Do not recognize D12 type input '{}'. Using 'BULK' by default.".format(D12_type))
        D12_type = "BULK"
    elif D12_type == "SLAB":
        SLAB_Warning = "Warning: You request a 'SLAB' D12, but this script is not functionalized for that.\n"
        SLAB_Warning += "However, the script helpscripts/code/create_d12_from_cif.py, located within the mendozacortesgroup GitHub, can produce SLABs.\n"
        SLAB_Warning += "I suggest adapting that script into a function with the same inputs/outputs as this function (instantiate_D12).\n"
        SLAB_Warning += "Alternatively you could adapt extract_BULK_reference_geometry_with_sym for SLABs.\n"
        SLAB_Warning += "Cannot convert {} to d12.".format(input_geometry_file)
        print(SLAB_Warning)
        
    
    structure = read(input_geometry_file)
    
    with open(template_D12, 'r') as file:
        D12_lines = file.readlines()
    
    for i in range(len(D12_lines)):
        if D12_template_keywords["Title"] in D12_lines[i]:
            D12_lines[i] = D12_lines[i].replace(D12_template_keywords["Title"], calculation_name)
        if D12_template_keywords["Reference Geometry"] in D12_lines[i]:
            reference_geometry_text = extract_BULK_reference_geometry_with_sym(structure, input_file=input_geometry_file)
            D12_lines[i] = D12_lines[i].replace(D12_template_keywords["Reference Geometry"], reference_geometry_text)
        # if D12_template_keywords["Geometry Editing"] in D12_lines[i]:
        #     D12_lines[i] = D12_lines[i].replace(D12_template_keywords["Geometry Editing"], "END")
        if D12_template_keywords["Basis Set"] in D12_lines[i]:
            basis_set_text = compile_basis_set(structure, basis_set_directory, input_file=input_geometry_file)
            D12_lines[i] = D12_lines[i].replace(D12_template_keywords["Basis Set"], basis_set_text)
        if D12_template_keywords["k-Points"] in D12_lines[i]:
            k_points_text = determine_k_point_grid(structure, k_point_criterion)
            D12_lines[i] = D12_lines[i].replace(D12_template_keywords["k-Points"], k_points_text)
    
    #print("\nEXAMPLE OUTPUT BELOW\n"+"-"*60)
    #print("".join(D12_lines))
    with open(output_D12, 'w') as file:
        for line in D12_lines:
            file.write( line)
 
parser = argparse.ArgumentParser()

parser.add_argument('-f', '--filename', required=True, help="Reference geometry filename")
parser.add_argument('-of', '--output_filename', help="Output filename. Default is input with file extention changed to 'd12'.")
parser.add_argument('-tf', '--template_filename', required=True, help="Template filename. This file contains all d12input except calculation name, reference geometry, basis sets and k-point grid.")
parser.add_argument('-k', '--k_point', default=40, type=float, help="k-point criterion. The product of the number of k points and the lattice constant along some direction shall equal this criterion at minimum.")
parser.add_argument('-bs', '--basis_sets', required=True, help="Directory containing basis sets.")
parser.add_argument('-t', '--title', help="Calculation name. Default is derived from input geometry file.")
args = parser.parse_args()

if not os.path.isfile(args.filename): 
    print("Reference geometry file '{}' does not exist. Cannot create d12.".format(args.filename))
    sys.exit( 1 )

if args.output_filename :
    output_filename = args.output_filename
else:
    output_filename = args.filename.replace(".cif", ".d12")

if not os.path.isfile(args.template_filename): 
    print("Template d12 file '{}' does not exist. Cannot create d12.".format(args.template_filename))
    sys.exit( 2 )

if not os.path.isdir(args.basis_sets): 
    print("Basis file directory '{}' does not exist. Cannot create d12.".format(args.basis_sets))
    sys.exit( 3 )

if args.title:
    title = args.title
else:
    title = "D12 created from '{}' using instantiate_D12.py script".format(args.filename.split("/")[-1])



test_input_geometry_file = args.filename #"a_proj/a_calc/example_PbCuSeCl.cif"
test_output_D12 = output_filename #"a_proj/a_calc/example_PbCuSeCl.d12"
basis_sets_path = args.basis_sets #"basis_sets/"
template_D12_path = args.template_filename #"a_proj/template_d12.txt"
test_calculation_name = title #"example_PbCuSeCl"
test_k_point_criterion = args.k_point #40


instantiate_D12(input_geometry_file = test_input_geometry_file,
                output_D12=test_output_D12,
                template_D12=template_D12_path,
                calculation_name=test_calculation_name,
                basis_set_directory=basis_sets_path,
                k_point_criterion=test_k_point_criterion)
