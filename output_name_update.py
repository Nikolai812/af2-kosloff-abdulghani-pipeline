import argparse
import os
import json
import re
import shutil
from typing import List

def get_sequences_names(job_folder: str, path: str = ".") -> List[str]:
    """
    Returns the list of first-level subfolders (sequence names) inside job_folder.
    Prints a warning if files are found at the first level.

    :param job_folder: Name of the job folder to check.
    :param path: Path to the directory containing job_folder (default: current directory).
    :return: List of subfolder names.
    """
    full_path = os.path.join(path, job_folder)

    if not os.path.isdir(full_path):
        raise NotADirectoryError(f"{full_path} is not a valid directory.")

    entries = os.listdir(full_path)
    subfolders = []
    files = []

    for entry in entries:
        entry_path = os.path.join(full_path, entry)
        if os.path.isdir(entry_path):
            subfolders.append(entry)
        elif os.path.isfile(entry_path):
            files.append(entry)

    if files:
        print(f"Warning: The following files were found in {full_path}: {', '.join(files)}")

    return subfolders


def include_sequence_name_into_output_filenames(job_folder: str, path: str = "."):
    """
    Renames files in each sequence subfolder by prefixing them with the subfolder (sequence) name.

    Example:
        old_name: result.pdb
        new_name: Seq1_result.pdb

    :param job_folder: Name of the main job folder containing sequence subfolders.
    :param path: Path to the directory containing job_folder (default: current directory).
    """
    # Import the helper defined earlier
    sequences = get_sequences_names(job_folder, path)

    full_job_path = os.path.join(path, job_folder)

    for seq_name in sequences:
        seq_path = os.path.join(full_job_path, seq_name)

        for entry in os.listdir(seq_path):
            entry_path = os.path.join(seq_path, entry)

            # Skip subfolders, only process files
            if os.path.isfile(entry_path):
                new_name = f"{seq_name}_{entry}"
                new_path = os.path.join(seq_path, new_name)

                # Only rename if not already renamed
                if not entry.startswith(f"{seq_name}_"):
                    os.rename(entry_path, new_path)
                    print(f"Renamed: {entry} → {new_name}")
                else:
                    print(f"Skipped (already renamed): {entry}")


def process_ranking_and_write_summary(job_folder: str, path: str = "."):
    job_path = os.path.join(path, job_folder)
    processed_folder = f"{job_folder}_processed"
    processed_path = os.path.join(path, processed_folder)

    # Create the processed folder if it doesn't exist
    os.makedirs(processed_path, exist_ok=True)

    # Dictionary to store best models and their values
    best_models_data = {}

    # Iterate over each subfolder in the job_folder
    for subfolder_name in os.listdir(job_path):
        subfolder_path = os.path.join(job_path, subfolder_name)

        # Check if it's a directory
        if not os.path.isdir(subfolder_path):
            continue

        # Try to read the ranking_debug.json file
        ranking_file = os.path.join(subfolder_path, f"{subfolder_name}_ranking_debug.json")
        if not os.path.exists(ranking_file):
            ranking_file = os.path.join(subfolder_path, "ranking_debug.json")
            if not os.path.exists(ranking_file):
                print(f"Ranking file not found in {subfolder_name}")
                continue

        # Read the ranking file
        with open(ranking_file, 'r') as f:
            ranking_data = json.load(f)

        # Get the first value from the "order" array
        best_model = ranking_data["order"][0]

        # Get the corresponding pLDDT or (in case of multimer) iptm+ptm value
        scores = ranking_data.get("plddts") or ranking_data.get("iptm+ptm")
        if scores is None:
            print("")
            print (" \N{ballot x} !!!KEY ERROR!!! Neither 'plddts' nor 'iptm+ptm' found in ranking file.")
            print("")

        best_value = scores[best_model]

        # Add to the best_models_data dictionary
        best_models_data[f"{subfolder_name}_{best_model}"] = best_value

        # Extract the best number from the best_model string
        best_number = extract_best_number(best_model)

        # Find the corresponding PDB file
        pdb_file = os.path.join(subfolder_path, f"unrelaxed_{best_model}.pdb")
        if not os.path.exists(pdb_file):
            pdb_file = os.path.join(subfolder_path, f"{subfolder_name}_unrelaxed_{best_model}.pdb")
            if not os.path.exists(pdb_file):
                print(f"PDB file not found for {best_model} in {subfolder_name}")
                continue

        # Copy the PDB file to the processed directory with the new name
        new_pdb_name = f"{subfolder_name}_{best_number}.pdb"
        new_pdb_path = os.path.join(processed_path, new_pdb_name)
        shutil.copy(pdb_file, new_pdb_path)
        print(f"Copied {pdb_file} to {new_pdb_path}")

    # Write the best_models_data to best_models.json in the processed folder
    best_models_file = os.path.join(processed_path, "best_models.json")
    with open(best_models_file, 'w') as f:
        json.dump(best_models_data, f, indent=4)
    print(f"Best models data written to {best_models_file}")


def extract_best_number(best_model: str) -> int:
    # Extract the number between 'model_' and '_pred'
    match = re.search(r'model_(\d+)_(pred|multimer)', best_model)
    if match:
        return int(match.group(1))
    else:
        raise ValueError(f"Could not extract best number from {best_model}")


def select_unrelaxed_to_or_based_output(job_folder: str, path: str, or_output_path: str) -> None:
    job_path = os.path.join(path, job_folder)

    # Create the root OR-output directory if needed
    os.makedirs(or_output_path, exist_ok=True)

    # Iterate over the 1st-level subfolders inside the job_folder
    for or_subfolder_name in os.listdir(job_path):
        or_subfolder_path = os.path.join(job_path, or_subfolder_name)

        # Skip non-directories
        if not os.path.isdir(or_subfolder_path):
            continue

        print(f"\nProcessing subfolder: {or_subfolder_name}")
        print(f"Source path      : {or_subfolder_path}")
        print(f"Output root      : {or_output_path}")

        # Create destination subfolder under OR-output path
        dest_subfolder = os.path.join(or_output_path, or_subfolder_name)
        os.makedirs(dest_subfolder, exist_ok=True)
        print(f"Destination path : {dest_subfolder}")

        # List all files in the current OR subfolder
        for filename in os.listdir(or_subfolder_path):

            # Files to copy:
            # 1. those starting with 'unrelaxed_model'
            # 2. or exactly 'ranking_debug.json'
            if (
                filename.startswith("unrelaxed_model")
                or filename == "ranking_debug.json"
            ):
                src_file = os.path.join(or_subfolder_path, filename)

                # ---- Construct new filename ----
                # Prefix with OR subfolder name
                new_filename = f"{or_subfolder_name}_{filename}"

                # For unrelaxed_model files → remove "_pred_0" from the base name (if present),
                # but keep the original file extension (e.g. .pdb)
                if filename.startswith("unrelaxed_model"):
                    name, ext = os.path.splitext(filename)  # split off extension (".pdb", ".txt", etc.)
                    if name.endswith("_pred_0"):
                        name = name[:-len("_pred_0")]  # remove the suffix from the base name
                    new_filename = f"{or_subfolder_name}_{name}{ext}"

                dest_file = os.path.join(dest_subfolder, new_filename)

                # ---- Copy with overwrite ----
                if os.path.exists(dest_file):
                    print(f"⚠ File exists, will be overwritten: {dest_file}")
                else:
                    print(f"Copying file to {dest_file}")

                shutil.copy2(src_file, dest_file)  # overwrite enabled

    print("\n select_unrelaxed_to_or_based_output: \N{CHECK MARK} Done.\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input_folder", dest="input_folder",
        help="input jobnumber folder",
        default="000000"
    )

    parser.add_argument(
        "-p", "--path", dest="folder_path",
        help="jobnumber folder path",
        default="."
    )

    parser.add_argument(
        "-o", "--or_output_path", dest="or_output_path",
        help="output path based on OR, not job number",
        default="."
    )

    args = parser.parse_args()
    job_folder = args.input_folder
    folder_path = args.folder_path
    or_output_path = args.or_output_path

    # Selecting unrelaxed models and adding or_name in front,
    # copying to or_name based folder (not recreating if folder exists),
    # rewriting existing files
    select_unrelaxed_to_or_based_output(job_folder, folder_path, or_output_path)

    # Here all output of the script is written to the newly created folder: old_name+processed
    process_ranking_and_write_summary(job_folder, folder_path)


if __name__ == '__main__':
    main()
