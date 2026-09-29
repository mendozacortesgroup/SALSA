# SALSA, user setup guide

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23045693.svg)](https://doi.org/10.5281/zenodo.23045693)

## SALSA directory tree:

```sh
├── 1. Substitution
│   ├── generate_substitutions_textual.py
│   ├── icsd_parser.py
│   └── Semiconductors_Workflow.py
├── 2. Approximation
│   └── GenerateCandidateCompoundInventory.py
├── 3. USPEX_setup
│   ├── DefaultUSPEXinput
│   ├── POTCARs                # you supply the POTCAR files; none ship here
│   └── USPEX_scripts
├── 4. CRYSTAL_setup
│   ├── CRYSTAL_Scripts
│   ├── DefaultBasisSetsCRYSTAL
│   └── DefaultCRYSTALinput
├── data
│   ├── known_compositions_and_properties.csv
│   ├── property_space_selection_criteria.csv
│   └── template_inventory.csv
├── original_datasets
│   └── semiconductors_genome_dataset
├── Projects                  # created by setup_project.sh; not tracked
├── banner.txt
├── propagate.py
├── SALSA_WorkFlow_Manager.py
├── setup_env.sh              # writes SALSA_DIR / VASP_PP_PATH to ~/.bashrc
├── setup_project.sh
├── requirements.txt
├── CITATION.cff
├── CHANGELOG.md
├── LICENSE
└── README.md
```


## initial setup

### Install USPEX
Here is a link from Mendoza-Cortes Wiki on how to install USPEX and have it set up. [Click Here](https://sites.google.com/view/mendozagroup/codestutorials/uspex)


### Python Installation
SALSA requires **Python 3.10 or newer**. Install its Python dependencies with:

```sh
pip install --user -r requirements.txt
```

### **Important Setup**
These lines must be added to your .bashrc

These PATHS must be updated to where USPEX was installed
```sh
#### ------------- USPEX v.10.5.0 ------------- ####
export MCRROOT=/.../bin/USPEX_v10.5
export PATH=/.../application/archive/:$PATH
export USPEXPATH=/...bin/USPEX_v10.5/application/archive/src
export PYTHONPATH=/...bin/USPEX_v10.5/application/archive/src/FunctionFolder:$PYTHONPATH
###----------------------------------------------
```


The scripts locate this repository through the `SALSA_DIR` environment
variable. The easiest way to set it is to run `./setup_env.sh` from the repo
root (see [Setup](#setup)), which writes the block to your `~/.bashrc` for you.

To do it by hand instead, add this to `~/.bashrc` - note there are no spaces
around the `=`:
```sh
export SALSA_DIR="/path/to/SALSA"
```

### Projects Directory Explained
The project directory is dedicated to your specific project and is not intended for general use. All calculations and data related to your specific project will be stored here. An example project would be the artifical photosynthesis project for the SALSA paper. If your project is named *artificial_photosynthesis*, a dedicated directory will be created, and all related data and calculations will be stored there. Each project contains results from each of the SALSA steps, it also stores the USPEX and CRYSTAL inputs for the specific project. Any data specific to the project is store in its directory created bythe workflow manager script.

### 3. USPEX_setup and 4. CRYSTAL_setup Directories Explained:
Holds basis sets and reusable USPEX/CRYSTAL input templates that are not specific to any one project. VASP POTCAR files are **not** distributed with SALSA - see [Pseudopotentials](#pseudopotentials-stage-3-only) for how to supply your own.


## SALSA_WorkFlow_Manager.py
This is the main script, and it runs the other scripts that make up the SALSA process.

### Project Management:

* Creates and manages SALSA project directories with proper structure
* Handles project naming, validation, and conflict resolution
* Provides options for full or partial project overwrites
* Maintains detailed project logs with timestamps

### Dataset Integration:

* Lists and allows selection from available datasets in the original_datasets directory
* Supports multiple material property datasets for analysis
* Automatically links selected datasets to project workflows

### Multi-Step Workflow Execution:

#### Substitution Step - Create ionic substiution matrix

* Generates substitution matrices from crystallographic databases
* Identifies which ions can substitute for each other in compounds
* Creates probability matrices for viable substitutions

Utlizes scripts in the ``1. Substitution`` directory which contains generate_substitutions_textual.py, icsd_parser.py and Semiconductors_Workflow.py scripts.


#### Approximation Step - Filter for candidate compounds

Uses substitution data to interpolate new compound compositions
Predicts properties of candidate materials based on parent compounds
Creates inventories of promising materials for further study

Utlizes scripts in the ``2. Approximation`` directory which contains GenerateCandidateCompoundInventory.py script.

#### Propagation Step - Initiates computational calculations

* Sets up USPEX inputs
* Prepares CRYSTAL calculation templates
* Launches high-throughput computational screening

Utlizes propagate.py script

## generate_substitutions_textual.py

### Summary:
The program takes a large database of chemical compounds and figures out which ions (charged atoms) can replace each other in crystal structures. For example, it might discover that Na+ (sodium ion) frequently substitutes for K+ (potassium ion) because they appear in similar chemical environments across many compounds

### input:
Program expects two input files located in the ```../data/``` directory

1.) ```Element_Charge_Dictionary.txt```


A reference file that lists each chemical element and it's possible ionic charges. Format:
    
```sh
ElementSymbol AtomicNumber PossibleCharge1 PossibleCharge2 ...
```

Example:
```sh
Na 11 1
O 8 -2
Fe 26 2 3
```

2.) ```ICSD-compositions.csv```

A CSV file containing chemical formulas from the ICSD database, one per line. Format:

```sh
Na2O
CaF2
Fe2O3
Al2SiO5
```

### Process
Parse Chemical Formulas: Breaks down each compound (like "Na2O") into elements and their quantities. Assign Ionic Charges: For each compound, determines the most likely charge state for each element that makes the compound electrically neutral. Find Substitution Patterns: Compares compounds that are nearly identical except for one ion substitution (e.g., "NaCl" vs "KCl"). Calculate Substitution Probabilities: Counts how often each ion pair substitutes and normalizes by frequency


### Output
A text file ```substitutions_textual.txt``` with three columns:

```sh
ion_A ion_B alpha_AB
Na:1+ K:1+ 2.456789
O:2- S:2- 1.234567
```
Where:

* ```ion_A``` and ```ion_B``` are ions in the format "Element:Charge+/-"
* ```alpha_AB``` is the substitution probability (higher = more likely to substitute)

alphaAB is a substitution score that quantifies how likely ion A can be substituted by ion B in crystal

structures. It's based on empirical data from the ICSD (Inorganic Crystal Structure Database).

How it's calculated:
 1. Count Similar Environments: The algorithm finds pairs of compounds that differ by exactly one ion
 substitution. For example:
   - Compound 1: Na₂SO₄ (contains Na⁺, S⁶⁺, O²⁻)
   - Compound 2: K₂SO₄ (contains K⁺, S⁶⁺, O²⁻)
   - This counts as evidence that Na⁺ :left_right_arrow: K⁺ substitution is possible
 2. Build Substitution Matrix: For each ion pair (A,B), it counts how many times they appear in similar
    chemical environments (line 264-267):
    num_matches[A, B] += weight
 
 3. Normalize by Occurrence: The raw count is normalized by the geometric mean of how often each ion
    appears (line 276):
    fractional_matches[i, j] = num_matches[i, j] / sqrt(occurences[i] * occurences[j])
 
 4. Final alphaAB: The fractional match is divided by the mean fraction across all ion pairs (line 309):
    
    alpha_value = fractional_matches[i, j] / mean_fraction
    Interpretation:
    - alphaAB > 1: Ions A and B substitute more frequently than average
    - alphaAB = 1: Average substitution frequency
    - alphaAB < 1: Below-average substitution frequency
    - alphaAB = 0: No observed substitutions
    
    This metric helps predict which ions can replace each other in new materials, useful for materials   
    discovery and design.

### Breakdown of the substitution matrix:
This applies to Semiconductors_Workflow.py as well.

**Substitution Matrix Equations**
 Based on the code analysis, here are the key mathematical equations:
 1. Ion-Ion Match Counting
 For compounds that differ by exactly one ion substitution:
 num_matches[i,j] = Σ weight(compound_pair)
 where weight = min(repetition_count₁, repetition_count₂)
 2. Fractional Matches (Normalization)
 fractional_matches[i,j] = num_matches[i,j] / √(occurrences[i] × occurrences[j])
 This normalizes by the geometric mean of ion occurrences.
 3. Mean Fraction Calculation
 mean_fraction = mean(fractional_matches[all valid pairs])
 where valid pairs exclude rare ions (occurrences ≤ 10).
 4. Final Substitution Score (αAB)
 αAB = fractional_matches[A,B] / mean_fraction
 5. Compound Similarity Criterion
 Two compounds are considered similar if:
 difference = compound₁ + compound₂ - 2×(compound₁ · compound₂)
 sum(difference) = 2 (exactly 2 positions differ)
 6. Cost Function for Charge Assignment
 cost = |Σ(element_count × ionic_charge)| × 1000 + Σ(charge_indices)
 This ensures charge neutrality when assigning oxidation states.
 
 Mathematical Interpretation:
 - The matrix element MAB represents how frequently ions A and B appear in similar chemical environments
 - The normalization by √(nA × nB) accounts for the relative abundance of each ion
 - Division by mean_fraction scales the values so that αAB = 1 represents average substitutability
 - The approach is based on empirical co-occurrence patterns in crystal structures


## icsd_parser.py

### Summary
This script converts chemical data from text format to numerical format for algorithmic processing, handling two input types: ion substitution data (converting ion names to codes with substitution probabilities) and semiconductor data (converting compound names and properties to numerical values).

### inputs:
Program expects two input input files located ```./data/```

1.) ```substitutions_textual.txt``` (Required)

Generated by the previous ICSD substitution script. Contains ion substitution data:

```sh
ion_A ion_B alpha_AB
Na:1+ K:1+ 2.456789
O:2- S:2- 1.234567
```

2.) ```semiconductors_existing.txt``` (Optional)

Contains semiconductor compounds with their properties. Format:

```sh
Ga:3+ N:3- 3.4 -0.0624 -0.355
Cd:2+ S:2- 2.58 0.353 -0.667
```

Where the last three numbers are: band_gap, oxidation_potential, reduction_potential

### Process
The script performs four key operations: creates numerical codes for unique ions, extracts and simplifies ion charges to their lowest terms, encodes ion substitution pairs with probability values, and converts semiconductor formulas into lists of ion codes with associated properties.

### Output:

1.) ```ion_codes``` 

Maps numerical codes to ion names:

```sh
1 Na:1+
2 K:1+
3 O:2-
4 Cl:1-
```

2.) ```ions_values``` 

Contains substitution data in numerical format:

```sh
first_ion_code second_ion_code likelihood simplified_charge1 simplified_charge2
1 2 2.456789 1 1
3 5 1.234567 2 2
```
simplified charges are the charges after being reduced to their simplest ratio by dividing out common factors of their original charges.

3.) ```semiconductors_existing_values.csv``` (if semiconductor data exists)

Semiconductor compounds and properties in CSV format:

```sh
1,4,5.7,-3.2,2.1
2,3,3,12.1,-4.8,3.9
```
Where each row contains: ion code1, ion code2, ..., band gap, oxididation potential, reduction potential

## Semiconductors_Workflow.py
This script analyzes a crystal structure database to generate ion substitution probabilities by parsing chemical formulas into charge-balanced ionic compositions, identifying compound pairs that differ by single ionic substitutions, calculating substitution frequencies across similar chemical environments, and producing a normalized probability matrix that quantifies the likelihood of each possible ion exchange while filtering out improbable substitution solutions.

### Summary:

### Input:
The program expects these files in the directory specified by ```local_directory```:

1.) ```ICSD-compositions.csv``` (Required)

Contains chemical formulas from the ICSD database, one per line:

```sh
NaCl
CaF2
Al2O3
SiO2
```

2.) ```Element_Charge_Dictionary.txt``` (Required)

Lists elements and their common ionic charges:

```sh
Na 11 1
Ca 20 2
Al 13 3
O 8 -2
F 9 -1
Cl 17 -1
```

### Configuration
Important: You must set the correct path at the top of the script:

```sh
local_directory = "/path/to/your/data/"
```

### Process
This script generates ion substitution probability matrices from crystal structure databases by parsing chemical formulas, assigning optimal ionic charges through cost minimization, detecting single-ion substitution patterns via pairwise compound comparisons, and producing normalized probability matrices through frequency adjustment, symmetrization, and logarithmic scaling to quantify ion exchange likelihoods.

### Output:
1.) ```../results/similar_env_count_matrix.npy```

Raw count matrix showing how many times each ion pair was observed in similar environments.

2.) ```../results/substitution_matrix_alpha.png```

Heatmap visualization of the final substitution probability matrix.

3.) Terminal Output
Detailed statistics including:

Number of compounds processed successfully
Ion frequency distributions
Processing progress updates
Analysis of charge distribution patterns

## Generate_Candidate_Compound_Inventory.py
This script generates new semiconductor compounds by taking a database of known materials and their properties, then creating hybrid candidates through ionic substitution and interpolation. It identifies valid ion exchanges using substitution probability data, generates charge-neutral compositions with size constraints, and estimates properties through weighted per-atom linear interpolation of parent compounds. The system filters candidates against target property ranges with tolerance margins and exports qualifying materials as CSV files, effectively creating new compound possibilities like generating "NaBr" or "KCl" from parent compounds "NaCl" and "KBr" when substitution data shows Na↔K and Cl↔Br exchanges are favorable.

### Input:
1.) ```known_compositions_and_properties.csv``` (Required)

Contains known semiconductor compounds and their measured properties:

```sh
Compound,band_gap,oxidation_potential,reduction_potential
NaCl,5.7,-3.2,2.1
CaF2,12.1,-4.8,3.9
```


2.) ```frequent_charges_by_element.txt``` (Required)

Same format as the first script - lists elements and their possible ionic charges:

```sh
Na 11 1
Ca 20 2
F 9 -1
Cl 17 -1
```

3.) ```substitutions_textual.txt``` (Required)

Generated by the first script - contains ion substitution probabilities:

```sh
ion_A ion_B alpha_AB
Na:1+ K:1+ 2.456789
Cl:1- Br:1- 1.234567
```

4.) ```property_space_selection_criteria.csv``` (Required)

Defines the desired property ranges for candidate selection:

```sh
Property,Minimum,Maximum,Tolerance
band_gap,1.0,3.0,0.2
oxidation_potential,-5.0,-2.0,0.5
reduction_potential,1.0,4.0,0.3
```
5.) ```Example_Project/inventory.csv``` (Required)

Template file that defines the output format for high-performance computing applications.

##### Key Column Descriptions

| Section | Column | Purpose | Example Value |
|---------|--------|---------|---------------|
| **Basic** | `compound_name` | Chemical formula | `Na2KCl3` |
| **Basic** | `composition_dictionary` | Element ratios | `Na:2::K:1::Cl:3` |
| **Basic** | `source` | Compound origin | `new_candidate_compound` |
| **Basic** | `status` | Processing status | `pending`, `complete` |
| **Properties** | `interpolated_property_dictionary` | Predicted properties | `band_gap:2.34::energy:-1.23` |
| **Properties** | `interpolated_properties_in_target_region` | Meets target criteria | `1` (yes), `0` (no) |
| **Properties** | `interpolated_properties_in_ideal_region` | Meets ideal criteria | `1` (yes), `0` (no) |
| **Parents** | `parent_1_compound_name` | First parent compound | `NaCl` |
| **Parents** | `parent_1_fraction` | Parent 1 contribution | `0.67` |
| **Parents** | `parent_2_compound_name` | Second parent compound | `KCl` |
| **Parents** | `parent_2_fraction` | Parent 2 contribution | `0.33` |
| **USPEX** | `USPEX_status` | Structure prediction status | `running`, `completed` |
| **USPEX** | `USPEX_current_generation` | Algorithm generation | `15` |
| **USPEX** | `USPEX_best_enthalpy` | Lowest energy found | `-45.67` |
| **USPEX** | `USPEX_done` | Calculation complete | `True`, `False` |
| **CRYSTAL** | `CRYSTAL_status` | Property calculation status | `running`, `done` |
| **CRYSTAL** | `CRYSTAL_output_energy` | Formation energy | `-2.34` |
| **CRYSTAL** | `CRYSTAL_properties_in_target_region` | Meets computed targets | `1`, `0` |
| **CRYSTAL** | `CRYSTAL_done` | Calculation complete | `True`, `False` |
| **Final** | `champion_compound` | Best in batch | `True`, `False` |
| **Final** | `row_last_updated` | Last modified | `2024-01-16 20:00:00` |

### Process;
This script generates hybrid semiconductor compounds by parsing known chemical formulas into balanced ionic compositions, filtering compound pairs using substitution probability data to identify valid interpolation candidates, creating all possible hybrid compositions within size constraints, predicting properties through weighted compositional averages, and applying selection criteria to export qualifying candidates that meet specified property requirements.

### Output:
1.) ```candidate_compounds_upto_[size]_atom_generated_on_[date].csv/.pickle``` 

Complete database of all generated candidates with full compositional and property information.

2.) ```new_candidate_compound_inventory.csv``` 

Filtered candidates formatted for high-performance computing workflows, containing only compounds meeting the property criteria.

## substitution_result directory:
This directory stores the results of the scripts located in the ``1. Substitution`` directory, and each project has this directory

## Approximation_result directory:
This directory stores the results of the scripts located in the ``2. Approximation`` directory, and each project has this directory

## data directory:
This directory contains data used as inputs for the scripts located in Substitution and Approximation directories.
This directory also contains the outputs of the scripts, most of these outputs will be used as inputs for other scripts.
The data that the user will have to provide will be sotred here, and its unviversally used in all projects.


## Setup

Prerequisites: **Python 3.10 or newer**, and [`jq`](https://jqlang.github.io/jq/),
which `setup_project.sh` uses to read `defaults.json`.

```bash
./setup_env.sh                          # sets SALSA_DIR
./setup_env.sh /path/to/vasp/potentials # also sets VASP_PP_PATH (stage 3 only)
./setup_env.sh --print                  # preview without touching ~/.bashrc
```

Re-running replaces the block it previously wrote rather than appending a
duplicate, so it is safe to run as often as you like. Set `RC_FILE` to write
to a file other than `~/.bashrc`.

## Pseudopotentials (stage 3 only)

VASP PAW pseudopotentials (POTCAR files) are proprietary and are **not**
distributed with SALSA. Supply your own licensed copies, then populate
`3. USPEX_setup/POTCARs/` with:

```bash
# using $VASP_PP_PATH, as set by setup_env.sh:
"3. USPEX_setup/USPEX_scripts/compile_POTCARs_locally.sh"

# or naming the directory explicitly:
"3. USPEX_setup/USPEX_scripts/compile_POTCARs_locally.sh" /path/to/potpaw_PBE
```

The source must be laid out one directory per element, each holding a file
named `POTCAR` - the layout of the VASP `potpaw_PBE` distribution:

```
potpaw_PBE/Ag/POTCAR
potpaw_PBE/Si/POTCAR
```

This produces one file per element named `POTCAR_<Symbol>` in
`3. USPEX_setup/POTCARs/`; set `SALSA_POTCAR_DIR` to send them elsewhere. The
script reports how many it copied and exits non-zero if it found none.

**Only stage 3 (USPEX structure prediction) needs them** - stages 1, 2 and 4 do
not require VASP at all.

## Citing SALSA

If you use this software, please cite both the software and the paper it implements.

Software (this release):

> M. Djokic, G. Martinez, A. Aduenko, J. L. Mendoza-Cortes,
> *SALSA: Substitution Approximation evoLutionary Search and Ab-initio*,
> version 1.0.0, Zenodo (2026). https://doi.org/10.5281/zenodo.23045693

Paper:

> S. M. Stafford, A. Aduenko, M. Djokic, Y.-H. Lin, J. L. Mendoza-Cortes,
> *Transforming Materials Discovery for Artificial Photosynthesis:
> High-Throughput Screening of Earth-Abundant Semiconductors*,
> arXiv:2310.00118 (2023). https://doi.org/10.48550/arXiv.2310.00118

Machine-readable metadata is in `CITATION.cff`.

Software authorship and paper authorship are not the same list, and neither
contains the other. `CITATION.cff` lists the people who wrote this code;
`preferred-citation` within it gives the paper's byline, unchanged. Gabriel
Martinez adapted the group's separate working scripts into the single
generalizable pipeline released here, under the guidance of Marcus Djokic;
Alexander Aduenko wrote the similarity-matrix work that the substitution stage
was adapted from.

## Versions

| Version | Contents |
|---|---|
| **v1.0.0** | SALSA as published in 2023: substitution, approximation, USPEX structure prediction, CRYSTAL setup (stages 1-4). |
| **v2.0.0** | Adds stage 5, hybrid-DFT convex-hull phase stability analysis. |

## License

MIT - see `LICENSE`. Note that VASP pseudopotentials are separately licensed
and are not covered by it.
