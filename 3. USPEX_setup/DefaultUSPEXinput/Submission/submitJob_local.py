from __future__ import with_statement
from __future__ import absolute_import
from subprocess import check_output
import re
import sys
from io import open

# Adapted from script made by 'etikhonov'

def submitJob_local(index, executable, script_path): 
	
	# Step 1: Create Submission Script
	with open(script_path, 'r') as file:
		script_lines = file.readlines()
	for i in range(len(script_lines)):
		script_lines[i] = script_lines[i].replace("JOB_NAME", "USPEX_CALCULATION_NAME_"+str(index))
		script_lines[i] = script_lines[i].replace("COMMAND_EXECUTABLE", executable)
	with open("myrun", "w") as file:
		for line in script_lines:
			file.write(line)

	# Step 2:
	output = str(check_output('sbatch myrun', shell=True))
	
	# Step 3
	# Here we parse job ID from the output of previous command
	print(output)
	jobNumber = int(re.findall(r'\d+', output)[0])
	
	return jobNumber


import os
# Resolved from the SALSA_DIR environment variable (see setup_env.sh).
# Falls back to a path relative to this file so the repo works when
# SALSA_DIR is unset.
_salsa_dir = os.environ.get("SALSA_DIR")
if _salsa_dir:
    local_submit_script_path = os.path.join(
        _salsa_dir, "3. USPEX_setup", "DefaultUSPEXinput", "Submission",
        "local_submit_script")
else:
    local_submit_script_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "local_submit_script")

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', dest='index', type=int)
	parser.add_argument('-c', dest='executable')
	args = parser.parse_args()

	jobNumber = submitJob_local(args.index, args.executable, local_submit_script_path)
	print('<CALLRESULT>')
	print(int(jobNumber))
