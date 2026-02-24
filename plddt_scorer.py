import os
import csv
import argparse


class PLDDTScorer:

    @classmethod
    def extract_plddt_by_residue(cls, pdb_file) -> list[tuple[str, int, float]]:
        rows = []
        seen = set()

        with open(pdb_file) as f:
            for line in f:
                if line.startswith("ATOM"):
                    fields = line.split()
                    chain = fields[4]
                    seq_id = int(fields[5])
                    chain_seq_id = (chain, seq_id)

                    if chain_seq_id not in seen:
                        seen.add(chain_seq_id)
                        plddt = float(fields[10])
                        rows.append((chain, seq_id, plddt))

        return rows

    @classmethod
    def write_to_csv(cls, pdb_name: str, rows: list[tuple[str, int, float]], output_dir: str) -> None:
        output_path = os.path.join(output_dir, f"{pdb_name}_plddt.csv")

        with open(output_path, mode="w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["chain", "residue_number", "plddt"])
            writer.writerows(rows)

    @classmethod
    def process_directory(cls, directory: str) -> int:
        """
        Process a single directory (non-recursive).
        Returns number of processed PDB files.
        """
        count = 0

        for filename in os.listdir(directory):
            if filename.endswith(".pdb"):
                pdb_path = os.path.join(directory, filename)
                pdb_name = os.path.splitext(filename)[0]

                rows = cls.extract_plddt_by_residue(pdb_path)
                cls.write_to_csv(pdb_name, rows, directory)

                count += 1
                print(f"Processed: {filename}")

        return count

    @classmethod
    def get_plddts_for_directory(cls, input_directory: str, recursive: bool = False) -> None:
        """
        Process PDB files in a directory.
        If recursive=True, traverse subdirectories.
        """
        total_processed = 0

        if recursive:
            for root, _, _ in os.walk(input_directory):
                count = cls.process_directory(root)
                print(f"{root}: {count} file(s) processed")
                if count > 0:
                    total_processed += count
        else:
            count = cls.process_directory(input_directory)
            print(f"{input_directory}: {count} file(s) processed")
            total_processed += count

        print(f"\nTotal processed PDB files: {total_processed}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract per-residue PLDDT scores from PDB files."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input directory containing PDB files"
    )

    parser.add_argument(
        "-r",
        action="store_true",
        help="Process subdirectories recursively"
    )

    args = parser.parse_args()

    PLDDTScorer.get_plddts_for_directory(
        input_directory=args.input,
        recursive=args.r
    )
