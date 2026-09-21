#!/usr/bin/env python3
# This script finalized by Sean
# Adapted slightly from a function created by Marcus which was adapted from code by Danny

import argparse
import glob
import matplotlib
matplotlib.use('tkagg') # <-- THIS MAKES IT FAST!
import os

# Determines the bandgap in eV as well as the bandgap type
def get_bandgap(output_content, return_type=True, verbose=False):
    # converts Ha to eV
    HeV    = 27.2114 
    
    #initialize
    index_alpha = 0
    index_beta = 0
    index_direct = 0
    index_cond = 0 
    index_indirect = 0
    alpha_cond_band =[]
    alpha_val_band =[]
    beta_cond_band = []
    beta_val_band = []
    bandgap_type = "Undetermined"
    
    #check which keyword applies
    for line in output_content:
        if "ALPHA      ELECTRONS" in line:
            index_alpha = output_content.index(line)
        if "BETA       ELECTRONS" in line:
            index_beta = output_content.index(line)
        if line.startswith(" DIRECT ENERGY BAND GAP"):
            index_direct = output_content.index(line)
        if "POSSIBLY CONDUCTING STATE" in line:
            index_cond = output_content.index(line)
        if "INDIRECT ENERGY BAND GAP" in line:
            index_indirect = output_content.index(line)
    if verbose:
        print('alpha line #:' + str(index_alpha))
        print('beta line #:' + str(index_beta))
        print('direct line #:' + str(index_direct))
        print('cond line #:' + str(index_cond))
        print('indirect line #:' + str(index_indirect))
    # direct is last
    if index_direct > index_alpha and index_direct > index_beta and index_direct > index_cond and index_direct > index_indirect:
        for index in range(index_direct-4,index_direct+1):
            if "TOP OF VALENCE" in output_content[index]:
                unclean_alpha = [x for x in output_content[index].split(" ") if    x != "" ]
                alpha_val_band.append(unclean_alpha[10].split(';')[0])
            if "BOTTOM OF VIRTUAL" in output_content[index]:
                unclean_alpha = [x for x in output_content[index].split(" ") if    x != "" ]
                alpha_cond_band.append(unclean_alpha[10].split(';')[0])
        beta_cond_band = alpha_cond_band
        beta_val_band = alpha_val_band
        bandgap_type = "Direct"
    #cond is last
    elif index_cond > index_alpha and index_cond > index_beta and index_cond > index_direct and index_cond > index_indirect:
        unclean_alpha = [x for x in output_content[index_cond].split(" ") if    x != "" ]
        alpha_val_band.append(unclean_alpha[5].split(';')[0])
        alpha_cond_band = alpha_val_band
        beta_cond_band = alpha_cond_band
        beta_val_band = alpha_val_band
        bandgap_type = "None"
    # alpha/beta is last
    elif index_alpha > index_cond and index_alpha > index_direct and index_alpha > index_indirect:
        for index in range(index_alpha+2,index_alpha+7):
            if "TOP OF VALENCE" in output_content[index]:
                unclean_alpha = [x for x in output_content[index].split(" ") if    x != "" ]
                alpha_val_band.append(unclean_alpha[10].split(';')[0])
            if "BOTTOM OF VIRTUAL" in output_content[index]:
                unclean_alpha = [x for x in output_content[index].split(" ") if    x != "" ]
                alpha_cond_band.append(unclean_alpha[10].split(';')[0])
        for index in range(index_beta+2,index_beta+7):
            if "TOP OF VALENCE" in output_content[index]:
                unclean_beta = [x for x in output_content[index].split(" ") if    x != "" ]
                beta_val_band.append(unclean_beta[10].split(';')[0])
            if "BOTTOM OF VIRTUAL" in output_content[index]:
                unclean_beta = [x for x in output_content[index].split(" ") if    x != "" ]
                beta_cond_band.append(unclean_beta[10].split(';')[0])
    #indirect is last
    elif index_indirect > index_direct and index_indirect > index_cond and index_indirect > index_alpha:
        for index in range(index_indirect-4,index_indirect+1):
            if "TOP OF VALENCE" in output_content[index]:
                unclean_alpha = [x for x in output_content[index].split(" ") if    x != "" ]
                alpha_val_band.append(unclean_alpha[10].split(';')[0])
            if "BOTTOM OF VIRTUAL" in output_content[index]:
                unclean_alpha = [x for x in output_content[index].split(" ") if    x != "" ]
                alpha_cond_band.append(unclean_alpha[10].split(';')[0])
        beta_cond_band = alpha_cond_band
        beta_val_band = alpha_val_band
        bandgap_type = "Indirect"
    else:
            print('error, no keyword')
    # placeholder
    min_alpha_c = alpha_cond_band[0]
    max_alpha_v = alpha_val_band[0]
    min_beta_c = beta_cond_band[0]
    max_beta_v = beta_val_band[0]
    # check
    if len(alpha_cond_band) == 1:
        min_alpha_c = alpha_cond_band[0]
    else:
        for shell in alpha_cond_band:
            if float(shell) < float(alpha_cond_band[0]):
                min_alpha_c = shell
                
    if len(alpha_val_band) == 1:
        max_alpha_v = alpha_val_band[0]
    else:
        for shell in alpha_val_band:
            if float(shell) > float(alpha_val_band[0]):
                max_alpha_v = shell

    if len(beta_cond_band) == 1:
        min_beta_c = beta_cond_band[0]
    else:
        for shell in beta_cond_band:
            if float(shell) < float(beta_cond_band[0]):
                min_beta_c = shell

    if len(beta_val_band) == 1:
        max_beta_v = beta_val_band[0]
    else:
        for shell in beta_val_band:
            if float(shell) < float(beta_val_band[0]):
                max_beta_v = shell
    
    if float(min_alpha_c) > float(min_beta_c):
        min_c = min_beta_c
    else:
        min_c = min_alpha_c
    if float(max_alpha_v) > float(max_beta_v):
        max_v = max_alpha_v
    else:
        max_v = max_beta_v
        
    conduction_min, valence_max = float(min_c)*HeV, float(max_v)*HeV
    
    if return_type:
        return conduction_min - valence_max, bandgap_type
    else:
        return conduction_min - valence_max

#################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filename', required=True, help="Crystal output filename")
parser.add_argument('-nt', '--no_type', action='store_true', help="Don't print bandgap type. Otherwise it is printed by default.")
args = parser.parse_args()

#################################################################################
output_content = []

with open(args.filename, 'r') as f:
    for line in f:
        if "SCF ENDED" in line:
            break
        else:
            output_content.append(line)

if len(output_content):
    bg_val, bg_type = get_bandgap(output_content)
    bg_print = "{:.5f}".format(bg_val)
    if (not args.no_type):
        bg_print += "     {:15s}".format(bg_type)
    print(bg_print)
else:
    print("No output detected.")
