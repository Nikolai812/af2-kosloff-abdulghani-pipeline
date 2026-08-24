import argparse
import csv
import logging
import pickle
import sys
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

    For this implementation we preserve every non-empty line
    as a separate TSV/CSV row.
    """

    lines = []

    with open(file_path, "r") as f:
        for line in f:
            line = line.rstrip()

            if line:
                lines.append(line)

    return lines


def write_section(
    writer,
    database_name,
    records,
    section_type="sequence",
):
    """
    Write one database section into the output file.
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


def create_msa_csv(
    or_name_dir,
    output_dir,
    separator="\t",
    extension=".tsv",
):
    """
    Create {OR_NAME}_msa.tsv or {OR_NAME}_msa.csv.

    Parameters:
        or_name_dir: AlphaFold2 OR_NAME output directory.
        output_dir: Directory where the output file is written.
        separator: Field delimiter. Tab by default.
        extension: Output file extension.
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
    # Verify expected files.
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
    # Create output directory if necessary.
    #
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_file = output_dir / f"{or_name}_msa{extension}"

    logger.info(
        "Creating MSA %s for '%s': %s",
        extension[1:].upper(),
        or_name,
        output_file,
    )

    with open(
        output_file,
        "w",
        newline="",
        encoding="utf-8",
    ) as f:

        writer = csv.writer(
            f,
            delimiter=separator,
        )

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
        "Created MSA %s: %s",
        extension[1:].upper(),
        output_file,
    )

    return output_file


def create_pae_csv(
    or_name_dir,
    output_dir,
    separator="\t",
    extension=".tsv",
):
    """
    Extract the PAE matrix from the result pickle corresponding
    to the best-ranked AlphaFold2 model and write it to:

        {OR_NAME}_pae.tsv

    or:

        {OR_NAME}_pae.csv

    Parameters:
        or_name_dir: AlphaFold2 OR_NAME output directory.
        output_dir: Directory where the output file is written.
        separator: Field delimiter. Tab by default.
        extension: Output file extension.

    Returns:
        Path to the created file, or None if PAE could not be created.
    """

    or_name = or_name_dir.name

    ranking_file = or_name_dir / "ranking_debug.json"

    if not ranking_file.is_file():
        logger.warning(
            "ranking_debug.json does not exist for '%s': %s",
            or_name,
            ranking_file,
        )
        return None

    #
    # Read ranking information.
    #
    import json

    with open(
        ranking_file,
        "r",
        encoding="utf-8",
    ) as f:
        ranking = json.load(f)

    #
    # AlphaFold2 stores models in ranking order in the
    # "order" list. The first model is the best-ranked model.
    #
    model_order = ranking.get("order")

    if not model_order:
        logger.warning(
            "No model order found in ranking_debug.json for '%s'.",
            or_name,
        )
        return None

    best_model = model_order[0]

    logger.info(
        "Best-ranked model for '%s': %s",
        or_name,
        best_model,
    )

    #
    # Corresponding result pickle.
    #
    result_file = (
        or_name_dir
        / f"result_{best_model}.pkl"
        #f"result_{best_model}_pred_0.pkl"
    )

    if not result_file.is_file():
        logger.warning(
            "Result file for best model does not exist for '%s': %s",
            or_name,
            result_file,
        )
        return None

    logger.info(
        "Reading PAE from %s",
        result_file,
    )

    #
    # Load AlphaFold2 result.
    #
    with open(
        result_file,
        "rb",
    ) as f:
        result = pickle.load(f)

    #
    # Extract PAE.
    #
    if "predicted_aligned_error" not in result:
        logger.warning(
            "predicted_aligned_error is not present in %s",
            result_file,
        )
        return None

    pae = result["predicted_aligned_error"]

    logger.info(
        "PAE matrix shape for '%s': %s",
        or_name,
        pae.shape,
    )

    #
    # Create output directory if necessary.
    #
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_file = output_dir / f"{or_name}_pae{extension}"

    logger.info(
        "Creating PAE %s for '%s': %s",
        extension[1:].upper(),
        or_name,
        output_file,
    )

    #
    # Write PAE matrix.
    #
    with open(
        output_file,
        "w",
        newline="",
        encoding="utf-8",
    ) as f:

        writer = csv.writer(
            f,
            delimiter=separator,
        )

        #
        # First row:
        #
        #        Residue 1  Residue 2  ...
        #
        header = ["residue"]

        for residue_number in range(
            1,
            pae.shape[1] + 1,
        ):
            header.append(residue_number)

        writer.writerow(header)

        #
        # Each subsequent row contains:
        #
        # residue_number + PAE values
        #
        for residue_number, pae_row in enumerate(
            pae,
            start=1,
        ):
            writer.writerow(
                [
                    residue_number,
                    *pae_row,
                ]
            )

    logger.info(
        "Created PAE %s: %s",
        extension[1:].upper(),
        output_file,
    )

    return output_file


def main():
    """
    Main entry point.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Extract AlphaFold2 MSA and PAE data from "
            "woutputs/<jobnumber>."
        )
    )

    parser.add_argument(
        "-j",
        "--jobnumber",
        required=True,
        help="Previously completed AlphaFold2 job number.",
    )

    parser.add_argument(
        "-s",
        "--separator",
        default="\t",
        help=(
            "Output delimiter. Use '\\t' for tab "
            "or ',' for comma. Default: '\\t'."
        ),
    )

    args = parser.parse_args()

    #
    # argparse receives the literal characters '\\t' when the
    # user specifies --separator "\\t".
    #
    # Convert them to an actual tab character.
    #
    separator = args.separator

    if separator == r"\t":
        separator = "\t"

    #
    # Determine output extension.
    #
    if separator == "\t":
        extension = ".tsv"
    elif separator == ",":
        extension = ".csv"
    else:
        logger.error(
            "Unsupported separator: %r. "
            "Only '\\t' and ',' are supported.",
            separator,
        )
        return 1

    #
    # Directory containing this Python script.
    #
    script_dir = Path(__file__).resolve().parent

    #
    # AlphaFold2 output directory.
    #
    woutputs_dir = script_dir / "woutputs"

    #
    # Input directory for this job.
    #
    job_dir = woutputs_dir / str(args.jobnumber)

    #
    # Output directory.
    #
    processed_dir = (
        woutputs_dir
        / f"{args.jobnumber}_processed"
    )

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

    logger.info(
        "Output separator: %r",
        separator,
    )

    logger.info(
        "Output extension: %s",
        extension,
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

        #
        # Create MSA output.
        #
        try:
            create_msa_csv(
                or_name_dir,
                processed_dir,
                separator,
                extension,
            )
        except Exception as ex:
            logger.exception(
                "Failed to create MSA for '%s'. "
                f"exeption: {ex}"
                "Continuing with PAE processing.",
                or_name_dir.name,
            )

        #
        # Create PAE output.
        #
        try:
            create_pae_csv(
                or_name_dir,
                processed_dir,
                separator,
                extension,
            )
        except Exception as ex:
            logger.exception(
                "Failed to create PAE for '%s'. "
                f"exeption: {ex}"
                "Continuing with the next OR_NAME.",
                or_name_dir.name,
            )
    logger.info(
        "MSA and PAE processing completed for job %s.",
        args.jobnumber,
    )

    return 0


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    sys.exit(main())