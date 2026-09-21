#!/bin/bash

ml GCC/9.3.0  OpenMPI/4.0.3  ASE/3.21.1-Python-3.8.2 > /dev/null 2>&1

if [[ -z "$SALSA_DIR" ]]; then
	echo "ERROR: SALSA_DIR is not set. Run setup_env.sh from the repository root." >&2
	exit 1
fi

if [ -d "$2" ]; then

	cd "$2" || exit 1
	calcname=$1
	infile="../reference_geometry.cif"
	outfile="${1}.d12"
	basis_sets="$SALSA_DIR/4. CRYSTAL_setup/DefaultBasisSetsCRYSTAL/"
	template_file=template_d12.txt
	python3 "$SALSA_DIR/4. CRYSTAL_setup/CRYSTAL_Scripts/instantiate_D12.py" -f "$infile" -of "$outfile" -k 40 -t "$calcname" -tf "$template_file" -bs "$basis_sets"
	mv template_submission.slurm submission.slurm
	sed -i "s/CALCULATION_NAME/$calcname/g" submission.slurm
fi
