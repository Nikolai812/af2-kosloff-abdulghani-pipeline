import os
import argparse

def split_to_monosequence(input_path, output_subdir):
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_subdir):
        os.makedirs(output_subdir)

    # Extract the pure file name from the input path
    pure_file_name = os.path.splitext(os.path.basename(input_path))[0]

    output_files = []
    with open(input_path, 'r') as file:
        sequences = file.read().split('>')[1:]  # Split sequences by '>' and ignore the first empty element

        for i, sequence in enumerate(sequences):
            lines = sequence.split('\n', 1)  # Split each sequence into header and body
            header = lines[0].strip()
            seq_file_name = header.replace('|', '-').replace(',', '-').replace(' ', '-')
            body = lines[1].replace('\n', '') if len(lines) > 1 else ''  # Remove newlines from the body

            # Create a filename for the output file
            output_filename = f"{seq_file_name}.fasta"
            output_filepath = os.path.join(output_subdir, output_filename)

            # Write the sequence to the output file
            with open(output_filepath, 'w') as output_file:
                output_file.write(f">{header}\n{body}\n")

            output_files.append(output_filename)

    return output_files

# Example usage:
# input_path = 'path/to/input.fasta'
# output_subdir = 'path/to/output_subdir'
# split_to_monosequence(input_path, output_subdir)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", dest="input_path",
                        help="input .fasta file",
                        default="./aaa")
    parser.add_argument("-o", "--output_subdirectory", dest="subdir",
                        help="output subdirectoy for splitted files",
                        default="SUBMONO")

    args = parser.parse_args()
#   print(f"args=\n{args}")
    input_path = args.input_path
    output_subdir = args.subdir

    output_files = split_to_monosequence(input_path, output_subdir)
    print(" ".join(output_files))
