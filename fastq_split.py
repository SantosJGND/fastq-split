import dnaio
import xopen
from typing import Iterable
import os

from fastq_filter import file_to_fastq_records

DEFAULT_COMPRESSION_LEVEL = 2


def estimate_fastq_records(filepath: str) -> int:
    """
    estimate number of fastq records in a file
    """

    count = 0
    for record in file_to_fastq_records(filepath):
        count += 1

    return count


def determine_extention(filepath: str) -> str:
    """
    determine file extension
    """

    if filepath.endswith("fastq.gz"):
        return "fastq.gz"
    elif filepath.endswith("fastq"):
        return "fastq"
    elif filepath.endswith("fq.gz"):
        return "fq.gz"
    elif filepath.endswith("fq"):
        return "fq"
    else:
        raise ValueError("File extension not recognized")


def subset_program(filepath: str, output_dir: str, nfiles: int, total_records: int):
    """
    return list of filepaths for subsetted fastq files"""

    filepath_ext = determine_extention(filepath)
    filepath_prefix = filepath[: -len(filepath_ext) - 1]

    records_per_file = total_records // nfiles

    subset_filepaths = [
        f"{output_dir}/{filepath_prefix}_{i+1}.{filepath_ext}" for i in range(nfiles)
    ]

    return subset_filepaths, records_per_file


def fastq_records_to_file_w_max(
    records: Iterable[dnaio.Sequence],
    filepath: str,
    subset_filepaths: list,
    records_per_file: int,
    compression_level: int = DEFAULT_COMPRESSION_LEVEL,
):
    """
    write fastq records to file, stop when max size is reached
    """

    written_records_current_file = 0
    current_file = subset_filepaths[0]

    with xopen.xopen(
        current_file, mode="wb", threads=0, compresslevel=compression_level
    ) as output_h:

        for record in records:

            record = record.fastq_bytes()

            if written_records_current_file == records_per_file:

                if len(subset_filepaths) > 1:
                    output_h.close()
                    subset_filepaths.pop(0)
                    current_file = subset_filepaths[0]
                    output_h = xopen.xopen(
                        current_file,
                        mode="wb",
                        threads=0,
                        compresslevel=compression_level,
                    )
                    written_records_current_file = 0

            output_h.write(record)
            written_records_current_file += 1

        output_h.close()


def fastq_split(
    filepath: str,
    output_dir: str,
    nfiles: int,
    compression_level: int = DEFAULT_COMPRESSION_LEVEL,
):
    """
    split fastq file into nfiles
    """

    total_records = estimate_fastq_records(filepath)
    subset_filepaths, records_per_file = subset_program(
        filepath, output_dir, nfiles, total_records
    )

    records = file_to_fastq_records(filepath)

    fastq_records_to_file_w_max(
        records,
        filepath,
        subset_filepaths,
        records_per_file,
        compression_level=compression_level,
    )

    return subset_filepaths


def main():

    import argparse

    parser = argparse.ArgumentParser(
        description="Split a fastq file into multiple files"
    )

    parser.add_argument("filepath", type=str, help="Path to the fastq file to split")

    parser.add_argument("output_dir", type=str, help="Path to the output directory")

    parser.add_argument(
        "nfiles", type=int, help="Number of files to split the fastq file into"
    )

    parser.add_argument(
        "--compression_level",
        type=int,
        default=DEFAULT_COMPRESSION_LEVEL,
        help="Compression level for output files",
    )

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    fastq_split(
        args.filepath,
        args.output_dir,
        args.nfiles,
        compression_level=args.compression_level,
    )


if __name__ == "__main__":
    main()
