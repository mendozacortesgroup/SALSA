#!/usr/bin/env python3
"""
SALSA Workflow Manager
A comprehensive script to setup and run the complete SALSA computational materials workflow.
"""

import os
import sys
import subprocess
import shutil
import json
import pandas as pd
from datetime import datetime
import argparse
from pathlib import Path

class SALSAWorkflowManager:
    def __init__(self):
        self.salsa_dir = os.getenv("SALSA_DIR")
        if not self.salsa_dir:
            raise EnvironmentError("SALSA_DIR environment variable not set")

        self.project_name = ""
        self.project_dir = ""
        self.dataset_path = ""
        self.selected_dataset = ""

        # Directory paths
        self.substitution_dir = os.path.join(self.salsa_dir, "1. Substitution")
        self.approximation_dir = os.path.join(self.salsa_dir, "2. Approximation")
        self.uspex_setup_dir = os.path.join(self.salsa_dir, "3. USPEX_setup")
        self.crystal_setup_dir = os.path.join(self.salsa_dir, "4. CRYSTAL_setup")
        self.projects_dir = os.path.join(self.salsa_dir, "Projects")
        self.datasets_dir = os.path.join(self.salsa_dir, "original_datasets")
        self.data_dir = os.path.join(self.salsa_dir, "data")  # Universal data directory

        # Log setup
        self.timestamp_format = "%H:%M:%S %Y-%m-%d"

        # Initialize workflow flags
        self.run_substitution = False
        self.run_approximation = False
        self.run_propagate = False
        self.substitution_start_script = 0

    def get_timestamp(self):
        return datetime.now().strftime(self.timestamp_format)

    def log_message(self, message, print_to_console=True):
        """Log message to project log file and optionally print to console"""
        timestamp = self.get_timestamp()
        log_entry = f"[{timestamp}] {message}"

        if print_to_console:
            print(log_entry)

        if hasattr(self, 'log_file') and self.log_file:
            with open(self.log_file, 'a') as f:
                f.write(log_entry + "\n")

    def prompt_project_name(self):
        """Prompt user for project name and validate it"""
        while True:
            self.project_name = input("Enter project name: ").strip()

            if not self.project_name:
                print("Project name cannot be empty. Please try again.")
                continue

            # Check for invalid characters
            invalid_chars = ['/', '\\', ':', '*', '?', '"', '<', '>', '|', ' ']
            if any(char in self.project_name for char in invalid_chars):
                print(f"Project name contains invalid characters: {invalid_chars}")
                continue

            self.project_dir = os.path.join(self.projects_dir, self.project_name)

            # Check if project already exists
            if os.path.exists(self.project_dir):
                choice = self.handle_existing_project()
                if choice == "retry":
                    continue
                elif choice == "full_overwrite":
                    shutil.rmtree(self.project_dir)
                    self._is_partial_overwrite = False
                    break
                elif choice == "partial_overwrite":
                    # Project exists, user wants partial overwrite
                    self._is_partial_overwrite = True
                    break
                elif choice == "cancel":
                    print("Project setup canceled.")
                    return False
            else:
                self._is_partial_overwrite = False
                break

        print(f"Project name set to: {self.project_name}")
        return True

    def handle_existing_project(self):
        """Handle existing project with multiple options"""
        print(f"\nProject '{self.project_name}' already exists.")
        print("Choose an option:")
        print("1. Full overwrite (delete everything and start fresh)")
        print("2. Partial overwrite (keep some results, re-run selected steps)")
        print("3. Choose different project name")
        print("4. Cancel")

        while True:
            choice = input("Enter your choice (1-4): ").strip()

            if choice == "1":
                confirm = input(f"This will DELETE ALL data in '{self.project_name}'. Are you sure? (y/N): ").lower()
                if confirm == 'y':
                    return "full_overwrite"
                else:
                    continue
            elif choice == "2":
                return "partial_overwrite"
            elif choice == "3":
                return "retry"
            elif choice == "4":
                return "cancel"
            else:
                print("Please enter 1, 2, 3, or 4")

    def handle_partial_overwrite(self):
        """Handle partial overwrite by selecting which steps to re-run"""
        print(f"\nPartial overwrite for project '{self.project_name}'")

        # Set up log file for partial overwrite
        self.log_file = os.path.join(self.project_dir, "log")

        # Mark this as partial overwrite
        self._is_partial_overwrite = True

        print("Select which step to start from:")
        print("1. Substitution step (re-run all substitution scripts)")
        print("2. Approximation step (keep substitution results)")
        print("3. Run propagate.py only (keep all analysis results)")

        # Check for existing substitution scripts if user wants to start from substitution
        substitution_scripts = [
            "generate_substitutions_textual.py",
            "icsd_parser.py",
            "Semiconductors_Workflow.py"
        ]

        while True:
            step_choice = input("Enter step to start from (1-3): ").strip()

            if step_choice == "1":
                # Starting from substitution - let user choose which script
                print("\nSubstitution step has multiple scripts:")
                for i, script in enumerate(substitution_scripts, 1):
                    print(f"  {i}. {script}")

                while True:
                    script_choice = input(f"Start from which script? (1-{len(substitution_scripts)}): ").strip()
                    try:
                        script_idx = int(script_choice) - 1
                        if 0 <= script_idx < len(substitution_scripts):
                            # Clear results from this script onwards
                            self.clear_partial_results("substitution", script_idx)
                            self.run_substitution = True
                            self.run_approximation = True
                            self.run_propagate = True
                            self.substitution_start_script = script_idx
                            return
                        else:
                            print(f"Please enter a number between 1 and {len(substitution_scripts)}")
                    except ValueError:
                        print("Please enter a valid number")

            elif step_choice == "2":
                # Starting from approximation
                self.clear_partial_results("approximation", 0)
                self.run_substitution = False
                self.run_approximation = True
                self.run_propagate = True
                return

            elif step_choice == "3":
                # Starting from propagate only
                self.run_substitution = False
                self.run_approximation = False
                self.run_propagate = True
                return

            else:
                print("Please enter 1, 2, or 3")

    def clear_partial_results(self, start_step, script_index=0):
        """Clear results from specified step onwards"""
        print(f"Clearing results from {start_step} step onwards...")

        if start_step == "substitution":
            substitution_result_dir = os.path.join(self.project_dir, "substitution_result")
            if os.path.exists(substitution_result_dir):
                # Clear specific files based on script index
                files_to_clear = []
                if script_index <= 0:  # Starting from generate_substitutions_textual
                    files_to_clear.extend(["substitutions_textual.txt"])
                if script_index <= 1:  # Starting from icsd_parser or earlier
                    files_to_clear.extend(["ions_codes", "ions_values", "semiconductors_existing_values.csv"])
                if script_index <= 2:  # Starting from Semiconductors_Workflow or earlier
                    files_to_clear.extend(["similar_env_count_matrix.npy", "substitution_matrix_alpha.png"])

                for file in files_to_clear:
                    file_path = os.path.join(substitution_result_dir, file)
                    if os.path.exists(file_path):
                        os.remove(file_path)
                        print(f"  Removed: {file}")

            # Also clear approximation results since they depend on substitution
            self.clear_partial_results("approximation", 0)

        elif start_step == "approximation":
            approximation_result_dir = os.path.join(self.project_dir, "Approximation_result")
            if os.path.exists(approximation_result_dir):
                # Clear all approximation results
                for file in os.listdir(approximation_result_dir):
                    file_path = os.path.join(approximation_result_dir, file)
                    if os.path.isfile(file_path):
                        os.remove(file_path)
                        print(f"  Removed: {file}")

        print("Partial cleanup completed.")

    def list_available_datasets(self):
        """List all available datasets in original_datasets directory"""
        if not os.path.exists(self.datasets_dir):
            print(f"Error: Datasets directory not found at {self.datasets_dir}")
            return []

        datasets = []
        for item in os.listdir(self.datasets_dir):
            dataset_path = os.path.join(self.datasets_dir, item)
            if os.path.isdir(dataset_path):
                datasets.append(item)

        return sorted(datasets)

    def prompt_dataset_selection(self):
        """Prompt user to select a dataset"""
        datasets = self.list_available_datasets()

        if not datasets:
            print("No datasets found in original_datasets directory")
            return False

        print("\nAvailable datasets:")
        for i, dataset in enumerate(datasets, 1):
            print(f"  {i}. {dataset}")

        while True:
            try:
                choice = input(f"\nSelect dataset (1-{len(datasets)}): ").strip()
                if not choice:
                    continue

                choice_idx = int(choice) - 1
                if 0 <= choice_idx < len(datasets):
                    self.selected_dataset = datasets[choice_idx]
                    self.dataset_path = os.path.join(self.datasets_dir, self.selected_dataset)
                    break
                else:
                    print(f"Please enter a number between 1 and {len(datasets)}")
            except ValueError:
                print("Please enter a valid number")

        print(f"Selected dataset: {self.selected_dataset}")
        return True

    def setup_project_structure(self):
        """Create the basic project directory structure"""
        print(f"\nSetting up project structure for '{self.project_name}'...")

        # Create main project directory
        os.makedirs(self.project_dir, exist_ok=True)

        # Setup log file and initial log entries (like bash script)
        self.log_file = os.path.join(self.project_dir, "log")

        with open(self.log_file, 'w') as f:
            f.write(f"SALSA project log for {self.project_name}\n")

            # Add banner (like bash script)
            banner_file = os.path.join(self.salsa_dir, "banner.txt")
            if os.path.exists(banner_file):
                with open(banner_file, 'r') as banner:
                    f.write(banner.read())

            f.write(f"Project set up at {self.get_timestamp()}\n")

        # Copy empty inventory template (like bash script)
        empty_inventory = os.path.join(self.data_dir, "template_inventory.csv")

        if os.path.exists(empty_inventory):
            shutil.copy2(empty_inventory, os.path.join(self.project_dir, "inventory.csv"))
            # Also copy as blank_inventory.csv (matching your structure)
            shutil.copy2(empty_inventory, os.path.join(self.project_dir, "template_inventory.csv"))
        else:
            # Create minimal inventory.csv if template doesn't exist
            inventory_columns = ["index", "compound_name", "compound_ID", "composition_dictionary",
                               "source", "status", "interpolated_property_dictionary"]
            df = pd.DataFrame(columns=inventory_columns)
            df.to_csv(os.path.join(self.project_dir, "inventory.csv"), index=False)
            df.to_csv(os.path.join(self.project_dir, "blank_inventory.csv"), index=False)

        # Create project directories (just basic structure, let propagate.py handle computational setup)
        directories_to_create = [
            "substitution_result",
            "Approximation_result"
        ]

        for directory in directories_to_create:
            os.makedirs(os.path.join(self.project_dir, directory), exist_ok=True)

        self.log_message("Created project directory structure with result folders")

        print(f"Project structure created at: {self.project_dir}")

    def prompt_workflow_steps(self):
        """Prompt user which workflow steps to run"""
        print("\nSALSA Workflow Steps:")
        print("1. Substitution - Generate substituted compounds")
        print("2. Approximation - Generate candidate compound inventory")
        print("3. Run propagate.py - Start USPEX and CRYSTAL calculations")

        self.run_substitution = self.prompt_yes_no("\nRun substitution step? (y/N): ")
        self.run_approximation = self.prompt_yes_no("Run approximation step? (y/N): ")
        self.run_propagate = self.prompt_yes_no("Run propagate.py? (y/N): ")

        # Initialize script starting point
        self.substitution_start_script = 0

    def prompt_yes_no(self, message, default=False):
        """Prompt user for yes/no response"""
        while True:
            response = input(message).strip().lower()
            if response in ['y', 'yes']:
                return True
            elif response in ['n', 'no', '']:
                return default
            else:
                print("Please enter 'y' for yes or 'n' for no")

    def run_substitution_step(self):
        """Run the substitution workflow in the correct order"""
        if not self.run_substitution:
            return

        print("\n" + "="*50)
        print("RUNNING SUBSTITUTION STEP")
        print("="*50)

        self.log_message("Starting substitution step")

        # Define the order of substitution scripts
        substitution_script_order = [
            "generate_substitutions_textual.py",
            "icsd_parser.py",
            "Semiconductors_Workflow.py"
        ]

        # Run substitution scripts starting from specified script
        for i, script in enumerate(substitution_script_order):
            if i < self.substitution_start_script:
                print(f"Skipping: {script} (starting from script {self.substitution_start_script + 1})")
                continue

            script_path = os.path.join(self.substitution_dir, script)

            if not os.path.exists(script_path):
                self.log_message(f"Warning: {script} not found in substitution directory")
                continue

            self.log_message(f"Running substitution script: {script}")
            print(f"Currently running: {script}")

            try:
                # Change to project directory for script execution
                os.chdir(self.project_dir)

                # Run the script with dataset and output paths
                result = subprocess.run([
                    sys.executable, script_path,
                    '--dataset', self.dataset_path,
                    '--output', self.project_dir
                ], capture_output=True, text=True, timeout=3600)

                if result.returncode == 0:
                    self.log_message(f"Successfully completed {script}")
                    print(f"✓ Completed: {script}")
                    if result.stdout:
                        self.log_message(f"Output: {result.stdout.strip()}")
                else:
                    self.log_message(f"Error in {script}: {result.stderr}")
                    print(f"✗ Failed: {script}")
                    # Don't continue if a script fails
                    self.log_message(f"Stopping substitution workflow due to error in {script}")
                    return

            except subprocess.TimeoutExpired:
                self.log_message(f"Timeout running {script}")
                print(f"⏱ Timeout: {script}")
                return
            except Exception as e:
                self.log_message(f"Exception running {script}: {str(e)}")
                print(f"✗ Exception in {script}: {str(e)}")
                return

        self.log_message("Substitution step completed successfully")

    def run_approximation_step(self):
        """Run the approximation workflow"""
        if not self.run_approximation:
            return

        print("\n" + "="*50)
        print("RUNNING APPROXIMATION STEP")
        print("="*50)

        self.log_message("Starting approximation step")

        # Find approximation scripts
        approximation_scripts = []
        for script in os.listdir(self.approximation_dir):
            if script.endswith('.py') and not script.startswith('__'):
                approximation_scripts.append(script)

        if not approximation_scripts:
            self.log_message("No approximation scripts found")
            return

        # Run each approximation script
        for script in sorted(approximation_scripts):
            script_path = os.path.join(self.approximation_dir, script)
            self.log_message(f"Running approximation script: {script}")
            print(f"Currently running: {script}")

            try:
                os.chdir(self.project_dir)

                result = subprocess.run([
                    sys.executable, script_path,
                    '--dataset', self.dataset_path,
                    '--output', self.project_dir
                ], capture_output=True, text=True, timeout=3600)

                if result.returncode == 0:
                    self.log_message(f"Successfully completed {script}")
                    print(f"✓ Completed: {script}")
                    if result.stdout:
                        self.log_message(f"Output: {result.stdout.strip()}")
                else:
                    self.log_message(f"Error in {script}: {result.stderr}")
                    print(f"✗ Failed: {script}")

            except subprocess.TimeoutExpired:
                self.log_message(f"Timeout running {script}")
                print(f"⏱ Timeout: {script}")
            except Exception as e:
                self.log_message(f"Exception running {script}: {str(e)}")
                print(f"✗ Exception in {script}: {str(e)}")

        self.log_message("Approximation step completed")

        # Replace inventory.csv with new candidate compounds for propagate.py
        self.update_inventory_with_candidates()

    def update_inventory_with_candidates(self):
        """Replace inventory.csv with new candidate compounds from approximation step"""
        approximation_result_dir = os.path.join(self.project_dir, "Approximation_result")
        candidate_inventory_file = os.path.join(approximation_result_dir, "new_candidate_compound_inventory.csv")
        main_inventory_file = os.path.join(self.project_dir, "inventory.csv")

        if os.path.exists(candidate_inventory_file):
            try:
                # Read the new candidate compounds
                candidate_df = pd.read_csv(candidate_inventory_file)

                # Save current inventory as backup
                if os.path.exists(main_inventory_file):
                    backup_file = os.path.join(self.project_dir, "inventory_backup.csv")
                    shutil.copy2(main_inventory_file, backup_file)
                    self.log_message(f"Backed up existing inventory to inventory_backup.csv")

                # Replace main inventory with candidate compounds
                candidate_df.to_csv(main_inventory_file, index=False)

                self.log_message(f"Updated inventory.csv with {len(candidate_df)} candidate compounds from approximation step")
                print(f"✓ Updated inventory.csv with {len(candidate_df)} candidate compounds")

            except Exception as e:
                self.log_message(f"Error updating inventory: {str(e)}")
                print(f"Warning: Could not update inventory.csv - {str(e)}")
        else:
            self.log_message("No new_candidate_compound_inventory.csv found - inventory.csv unchanged")
            print("Warning: No candidate compounds found to update inventory")

    def run_propagate_step(self):
        """Run propagate.py to start calculations"""
        if not self.run_propagate:
            return

        print("\n" + "="*50)
        print("RUNNING PROPAGATE.PY")
        print("="*50)

        self.log_message("Starting propagate.py")

        # Set up USPEX and CRYSTAL directories first (propagate.py expects these)
        self.setup_calculation_templates()

        # Path to propagate.py script
        propagate_script = os.path.join(self.salsa_dir, "propagate.py")

        if not os.path.exists(propagate_script):
            self.log_message(f"Error: propagate.py not found at {propagate_script}")
            print(f"Error: propagate.py not found at {propagate_script}")
            return

        self.log_message(f"Running propagate.py for project: {self.project_name}")
        print(f"Currently running: propagate.py -d {self.project_name}")

        try:
            # Change to SALSA directory for propagate.py execution
            os.chdir(self.salsa_dir)

            # Run propagate.py with project directory argument
            result = subprocess.run([
                sys.executable, propagate_script,
                '-d', self.project_name
            ], capture_output=True, text=True, timeout=7200)  # 2 hour timeout

            if result.returncode == 0:
                self.log_message("Successfully completed propagate.py")
                print("✓ Completed: propagate.py")
                if result.stdout:
                    self.log_message(f"Output: {result.stdout.strip()}")
            else:
                self.log_message(f"Error in propagate.py: {result.stderr}")
                print("✗ Failed: propagate.py")
                print(f"Error details: {result.stderr}")

        except subprocess.TimeoutExpired:
            self.log_message("Timeout running propagate.py (2 hours)")
            print("⏱ Timeout: propagate.py (2 hours)")
        except Exception as e:
            self.log_message(f"Exception running propagate.py: {str(e)}")
            print(f"✗ Exception in propagate.py: {str(e)}")

        self.log_message("Propagate step completed")

    def setup_calculation_templates(self):
        """Setup USPEX and CRYSTAL input templates for propagate.py"""
        print("Setting up calculation input templates...")

        # Create USPEX directories and copy templates
        uspex_inputs_dir = os.path.join(self.project_dir, "USPEX", "inputs")
        os.makedirs(uspex_inputs_dir, exist_ok=True)

        # Copy USPEX templates if directory is empty
        if not os.listdir(uspex_inputs_dir):
            default_uspex_input = os.path.join(self.salsa_dir, "3. USPEX_setup", "DefaultUSPEXinput")
            if os.path.exists(default_uspex_input):
                for item in os.listdir(default_uspex_input):
                    src = os.path.join(default_uspex_input, item)
                    dst = os.path.join(uspex_inputs_dir, item)
                    if os.path.isdir(src):
                        shutil.copytree(src, dst, dirs_exist_ok=True)
                    else:
                        shutil.copy2(src, dst)
                self.log_message("Copied USPEX input templates")
                print("✓ Copied USPEX input templates")

        # Create CRYSTAL directories and copy templates
        crystal_inputs_dir = os.path.join(self.project_dir, "CRYSTAL", "inputs")
        os.makedirs(crystal_inputs_dir, exist_ok=True)

        # Copy CRYSTAL templates if directory is empty
        if not os.listdir(crystal_inputs_dir):
            default_crystal_input = os.path.join(self.salsa_dir, "4. CRYSTAL_setup", "DefaultCRYSTALinput")
            if os.path.exists(default_crystal_input):
                for item in os.listdir(default_crystal_input):
                    src = os.path.join(default_crystal_input, item)
                    dst = os.path.join(crystal_inputs_dir, item)
                    if os.path.isdir(src):
                        shutil.copytree(src, dst, dirs_exist_ok=True)
                    else:
                        shutil.copy2(src, dst)
                self.log_message("Copied CRYSTAL input templates")
                print("✓ Copied CRYSTAL input templates")

    def setup_calculation_inputs(self):
        """Legacy method - now handled by propagate.py"""
        pass

    def run_workflow(self):
        """Run the complete SALSA workflow"""
        print("SALSA Workflow Manager")
        print("="*50)

        try:
            # Step 1: Get project name
            if not self.prompt_project_name():
                return

            # Step 2: Handle workflow based on project type
            if self._is_partial_overwrite:
                # For partial overwrite, still need dataset selection
                if not self.prompt_dataset_selection():
                    return
                # Handle partial overwrite - this sets the workflow flags
                self.handle_partial_overwrite()
            else:
                # Normal workflow for new projects or full overwrite
                # Step 2: Select dataset
                if not self.prompt_dataset_selection():
                    return

                # Step 3: Setup project structure
                self.setup_project_structure()

                # Step 4: Prompt for workflow steps
                self.prompt_workflow_steps()

            # Step 5: Run selected workflow steps
            self.run_substitution_step()
            self.run_approximation_step()
            self.run_propagate_step()

            # Step 7: Summary
            print("\n" + "="*50)
            print("WORKFLOW COMPLETED")
            print("="*50)

            self.log_message("SALSA workflow completed successfully")

            print("\nProject log printout:")
            with open(self.log_file, 'r') as f:
                print(f.read())

        except KeyboardInterrupt:
            print("\n\nWorkflow interrupted by user")
            if hasattr(self, 'log_file'):
                self.log_message("Workflow interrupted by user")
        except Exception as e:
            print(f"\nError: {str(e)}")
            if hasattr(self, 'log_file'):
                self.log_message(f"Workflow failed with error: {str(e)}")

def main():
    workflow = SALSAWorkflowManager()
    workflow.run_workflow()

if __name__ == "__main__":
    main()
