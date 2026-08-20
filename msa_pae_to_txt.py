import argparse
import csv
import logging
from pathlib import Path


logger = logging.getLogger(__name__)


#
# Expected AlphaFold2 MSA files.
#
MSA_FILES = [
    "bfd_uniclust_hits.a3m",
    "mgnify_hits.sto",
    "pdb_hits.hhr",
    "uniref90_hits.sto",
]


def parse_a3m(file_path):
    """
    Parse an A3M alignment file.

    Returns:
        list of tuples:
            (sequence_id, sequence)
    """

    records = []

    sequence_id = None
    sequence_parts = []

    with open(file_path, "r") as f:
        for line in f:
            line = line.rstrip()

            if not line:
                continue

            if line.startswith(">"):
                #
                # Save previous sequence.
                #
                if sequence_id is not None:
                    records.append(
                        (
                            sequence_id,
                            "".join(sequence_parts),
                        )
                    )

                #
                # Everything after '>' up to the first whitespace
                # is used as the sequence identifier.
                #
                sequence_id = line[1:].split()[0]
                sequence_parts = []

            else:
                sequence_parts.append(line)

    #
    # Save the last sequence.
    #
    if sequence_id is not None:
        records.append(
            (
                sequence_id,
                "".join(sequence_parts),
            )
        )

    return records


def parse_stockholm(file_path):
    """
    Parse a Stockholm (.sto) alignment file.

    Returns:
        list of tuples:
            (sequence_id, sequence)
    """

    sequences = {}

    with open(file_path, "r") as f:
        for line in f:
            line = line.rstrip()

            #
            # Ignore empty lines and Stockholm comments/metadata.
            #
            if not line:
                continue

            if line.startswith("#"):
                continue

            if line == "//":
                break

            fields = line.split()

            #
            # Normal Stockholm sequence line:
            #
            # sequence_name    sequence
            #
            if len(fields) >= 2:
                sequence_id = fields[0]
                sequence = fields[1]

                if sequence_id not in sequences:
                    sequences[sequence_id] = []

                sequences[sequence_id].append(sequence)

    #
    # A sequence can be split over several lines in Stockholm.
    #
    records = []

    for sequence_id, sequence_parts in sequences.items():
        records.append(
            (
                sequence_id,
                "".join(sequence_parts),
            )
        )

    return records


def read_hhr(file_path):
    """
    Read an HHsearch .hhr file.

    HHR is not an MSA file. It contains template-search results.

    For this first implementation we preserve every non-empty line
    as a separate CSV row so that the complete HHR information is
    available in the fourth section of the CSV.
    """

    lines = []

    with open(file_path, "r") as f:
        for line in f:
            line = line.rstrip()

            if line:
                lines.append(line)

    return lines


def write_section(writer, database_name, records, section_type="sequence"):
    """
    Write one database section into the CSV file.

    The section is preceded by a section marker.
    """

    writer.writerow(
        [
            f"### {database_name} ###",
            "",
            "",
            "",
        ]
    )

    if section_type == "sequence":

        writer.writerow(
            [
                "database",
                "record_type",
                "sequence_id",
                "sequence",
            ]
        )

        for sequence_id, sequence in records:
            writer.writerow(
                [
                    database_name,
                    "MSA_SEQUENCE",
                    sequence_id,
                    sequence,
                ]
            )

    elif section_type == "hhr":

        writer.writerow(
            [
                "database",
                "record_type",
                "line_number",
                "content",
            ]
        )

        for line_number, line in enumerate(records, start=1):
            writer.writerow(
                [
                    database_name,
                    "HHR_LINE",
                    line_number,
                    line,
                ]
            )

    #
    # Empty row separating sections.
    #
    writer.writerow([])


def create_msa_csv(or_name_dir, output_dir):
    """
    Create {OR_NAME}_msa.csv for one AlphaFold2 prediction.

    Parameters:
        or_name_dir: Path to the OR_NAME AlphaFold2 output directory.
        output_dir: Directory where the CSV should be written.

    Returns:
        Path to the created CSV file.
    """

    or_name = or_name_dir.name

    msas_dir = or_name_dir / "msas"

    if not msas_dir.is_dir():
        logger.warning(
            "MSA directory does not exist for '%s': %s",
            or_name,
            msas_dir,
        )
        return None

    #
    # Verify all expected files.
    #
    missing_files = []

    for msa_file in MSA_FILES:
        file_path = msas_dir / msa_file

        if not file_path.is_file():
            missing_files.append(msa_file)

    if missing_files:
        logger.warning(
            "Missing MSA files for '%s': %s",
            or_name,
            ", ".join(missing_files),
        )

    #
    # We can still create the CSV from the files that exist.
    #
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_file = output_dir / f"{or_name}_msa.csv"

    logger.info(
        "Creating MSA CSV for '%s': %s",
        or_name,
        output_file,
    )

    with open(
        output_file,
        "w",
        newline="",
        encoding="utf-8",
    ) as f:

        writer = csv.writer(f)

        #
        # General information.
        #
        writer.writerow(
            [
                "AlphaFold2 MSA data",
                "",
                "",
                "",
            ]
        )

        writer.writerow(
            [
                "OR_NAME",
                or_name,
                "",
                "",
            ]
        )

        writer.writerow([])

        #
        # ---------------------------------------------------------
        # 1. BFD / UniClust
        # ---------------------------------------------------------
        #
        file_path = msas_dir / "bfd_uniclust_hits.a3m"

        if file_path.is_file():
            logger.info(
                "Reading %s",
                file_path,
            )

            records = parse_a3m(file_path)

            write_section(
                writer,
                "bfd_uniclust_hits.a3m",
                records,
                section_type="sequence",
            )

        #
        # ---------------------------------------------------------
        # 2. Mgnify
        # ---------------------------------------------------------
        #
        file_path = msas_dir / "mgnify_hits.sto"

        if file_path.is_file():
            logger.info(
                "Reading %s",
                file_path,
            )

            records = parse_stockholm(file_path)

            write_section(
                writer,
                "mgnify_hits.sto",
                records,
                section_type="sequence",
            )

        #
        # ---------------------------------------------------------
        # 3. PDB HHR
        # ---------------------------------------------------------
        #
        file_path = msas_dir / "pdb_hits.hhr"

        if file_path.is_file():
            logger.info(
                "Reading %s",
                file_path,
            )

            records = read_hhr(file_path)

            write_section(
                writer,
                "pdb_hits.hhr",
                records,
                section_type="hhr",
            )

        #
        # ---------------------------------------------------------
        # 4. UniRef90
        # ---------------------------------------------------------
        #
        file_path = msas_dir / "uniref90_hits.sto"

        if file_path.is_file():
            logger.info(
                "Reading %s",
                file_path,
            )

            records = parse_stockholm(file_path)

            write_section(
                writer,
                "uniref90_hits.sto",
                records,
                section_type="sequence",
            )

    logger.info(
        "Created MSA CSV: %s",
        output_file,
    )

    return output_file


def main():
    """
    Main entry point.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Extract AlphaFold2 MSA data from woutputs/<jobnumber> "
            "and create CSV files."
        )
    )

    parser.add_argument(
        "jobnumber",
        help="Previously completed AlphaFold2 job number.",
    )

    args = parser.parse_args()

    #
    # Directory containing this Python script.
    #
    script_dir = Path(__file__).resolve().parent

    #
    # AlphaFold2 input/output directory.
    #
    woutputs_dir = script_dir / "woutputs"

    #
    # Input directory for this job.
    #
    job_dir = woutputs_dir / str(args.jobnumber)

    #
    # Output directory.
    #
    processed_dir = woutputs_dir / f"{args.jobnumber}_processed"

    logger.info(
        "Script directory: %s",
        script_dir,
    )

    logger.info(
        "AlphaFold2 job directory: %s",
        job_dir,
    )

    logger.info(
        "Processed output directory: %s",
        processed_dir,
    )

    if not job_dir.is_dir():
        logger.error(
            "AlphaFold2 job directory does not exist: %s",
            job_dir,
        )
        return 1

    #
    # Find all OR_NAME directories.
    #
    or_name_dirs = sorted(
        path
        for path in job_dir.iterdir()
        if path.is_dir()
    )

    if not or_name_dirs:
        logger.warning(
            "No OR_NAME directories found in: %s",
            job_dir,
        )
        return 0

    logger.info(
        "Found %d OR_NAME directories.",
        len(or_name_dirs),
    )

    #
    # Process every OR_NAME.
    #
    for or_name_dir in or_name_dirs:

        logger.info(
            "Processing OR_NAME: %s",
            or_name_dir.name,
        )

        create_msa_csv(
            or_name_dir,
            processed_dir,
        )

    #
    # This is where the next processing step will be added:
    #
    # create_pae(...)
    #

    logger.info(
        "MSA processing completed for job %s.",
        args.jobnumber,
    )

    return 0


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    raise SystemExit(main())