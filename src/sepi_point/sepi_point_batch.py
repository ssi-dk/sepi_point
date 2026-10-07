#!/usr/bin/env python3
import logging
logger = logging.getLogger(__name__)

import sys
from pathlib import Path
from sepi_point.mutation_finder import MutationFinder
from sepi_point.seqtools import WgsData
import argparse
import sepi_point.sepi_point as sp
from importlib import resources


def parse_args(argv):
    parser = argparse.ArgumentParser(description='Run SepiPOINT on a batch of isolates.')
    parser.add_argument("-r", "--read_dir",
                        help = "Folder with paired end read files (fastq or fastq.gz)",
                        type=Path,
                        required = False)
    parser.add_argument("-a", "--assembly_dir",
                        help = "Folder with genome assemblies (fasta)",
                        type=Path,
                        required = False)
    parser.add_argument("-o", "--output",
                        help = "Output directory.",
                        type=Path,
                        required = True)
    parser.add_argument("-n", "--no_clean",
                        help = "Do not clean up sam and bam files after variant calling. Default False.",
                        action= "store_true",
                        default=False)
    parser.add_argument("-l", "--log_level",
                        help = "Logging depth. Default: INFO",
                        type=str,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        default="INFO",
                        required = False)
    args = parser.parse_args()
    return args


def main_cli():
    args = parse_args(argv=sys.argv)
    
    ### Check if necessary inputs are provided and that folders exist
    if args.read_dir is not None:
        if args.read_dir.exists():
            print(f"Checking for read files in {args.read_dir}")
        else:
            print(f"Provided read input folder {args.read_dir} does not exist. Exitting.")
            sys.exit()
    if args.assembly_dir is not None:
        if args.assembly_dir.exists():
            print(f"Checking for assembly files in {args.assembly_dir}")
        else:
            print(f"Provided assembly input folder {args.assembly_dir} does not exist. Exitting.")
            sys.exit()
    elif args.read_dir is None:
        print("No inputs provided. Please provide an input folder with paired end reads (-r) and/or input folder with assemblies (-a). Exitting.")
        sys.exit()
    
    ### Set up output directory
    try:
        args.output.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        print(f"Failed to setup output directory at {args.output}: {e}. Exitting.")
        sys.exit()

    log_file = args.output.joinpath("SepiPOINT.log")
    logger = sp.setup_logger(log_file=log_file, log_level=args.log_level)
    mutation_db_tsv = resources.files("sepi_point").joinpath("db").joinpath("mutations.tsv")
    mutation_db_fasta = resources.files("sepi_point").joinpath("db").joinpath("sequences.fasta")
    logger.info(f"#### RUNNING SepiPOINT in batch mode ####")
    logger.info(f"Checking for mutations found in {mutation_db_tsv}")
    logger.info(f"Reference fasta file: {mutation_db_fasta}")

    ###
    wgs_data = WgsData.from_folders(assembly_data_folder=args.assembly_dir, paired_end_read_data_folder=args.read_dir)
    print(f"Loaded data for {len(wgs_data)} samples with read and/or assembly data.")
    logger.debug("Loaded data for %s samples with read and/or assembly data.", len(wgs_data))
    mf = MutationFinder()
    mf.load_and_check_db(mutation_db_tsv=mutation_db_tsv,sequence_db_fasta=mutation_db_fasta)

    all_sample_mutations, all_sample_putative_mutations = sp.run_sepi_point(
        mf = mf,
        wgs_data=wgs_data,
        output_dir=args.output,
        no_clean=args.no_clean,
    )
    logger.debug("Done running on %s samples. Now summarising results...", len(wgs_data))

    if args.output:
        summary_output_file = args.output.joinpath("results.tsv")
        summary_putative_output_file = args.output.joinpath("results.putative.tsv")
        matrix_output_file = args.output.joinpath("results.matrix.tsv")
        if not args.output.exists():
            args.output.mkdir()
    else:
        summary_output_file = args.results_dir.joinpath("results.tsv")
        summary_putative_output_file = args.results_dir.joinpath("results.putative.tsv")
        matrix_output_file = args.results_dir.joinpath("results.matrix.tsv")
    mf.print_sample_mutations_batch(
        mutation_summaries=all_sample_mutations,
        putative_mutation_summaries=all_sample_putative_mutations,
        summary_output_file=summary_output_file,
        summary_putative_output_file=summary_putative_output_file,
        matrix_output_file=matrix_output_file,
        )

    logger.info("#### Done running SepiPOINT on %s samples ####", len(all_sample_mutations))

if __name__ == "__main__":
    main_cli()