#!/usr/bin/env python3

import os
import sys
from pathlib import Path
import logging
import subprocess
from sepi_point.mutation_finder import MutationFinder
from sepi_point.seqtools import WgsData
import argparse
import re
import sepi_point.sepi_point as sp
from importlib import resources


def parse_args(argv):
    parser = argparse.ArgumentParser(description='Run sepi_point on a single isolate')
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

    log_file = args.output.joinpath("sepi_point.log")
    logger = sp.setup_logger(log_file=log_file)
    mutation_db_tsv = resources.files("sepi_point").joinpath("db").joinpath("mutations.tsv")
    mutation_db_fasta = resources.files("sepi_point").joinpath("db").joinpath("sequences.fasta")

    ###
    wgs_data = WgsData.from_folders(assembly_data_folder=args.assembly_dir, paired_end_read_data_folder=args.read_dir)
    print(f"Loaded data for {len(wgs_data)} samples with read and/or assembly data.")
    mf = MutationFinder()
    mf.load_and_check_db(mutation_db_tsv=mutation_db_tsv,sequence_db_fasta=mutation_db_fasta)

    all_sample_mutations= {}
    for sample_name, file_paths, info in wgs_data:
        if file_paths["r1_file"] is not None and file_paths["r2_file"] is not None:
            if file_paths["assembly_file"] is not None:
                output_name_fasta = f"{sample_name}_fasta"
                output_name_fastq = f"{sample_name}_fastq"
                output_dir_fasta = args.output.joinpath(output_name_fasta)
                output_dir_fasta.mkdir(exist_ok=True)
                logger.info(f"#####################################################################")
                logger.info(f"Running sepi_point on {sample_name}")
                logger.info(f"Using input {file_paths['assembly_file']}")
                logger.info(f"Printing output to {output_dir_fasta}")
                print(f"Running sepi_point on {sample_name}")
                snp_file = sp.run_nucmer_and_showsnps(assembly_file=file_paths["assembly_file"],
                                                   output_dir=output_dir_fasta, output_prefix=sample_name,
                                                   reference_fasta=mutation_db_fasta, logger=logger)
                sample_mutations = mf.get_mutations_from_nucmer_snps(nucmer_snp_file=snp_file)
                sample_mutation_summary = mf.summarize_sample_mutations(sample_mutations=sample_mutations).copy()
                mf.print_sample_mutations(mutation_summary=sample_mutation_summary,summary_output_file=output_dir_fasta.joinpath("results.tsv"))
                all_sample_mutations[output_name_fasta] = sample_mutation_summary
            else:
                output_name_fastq = sample_name
            output_dir_fastq = args.output.joinpath(output_name_fastq)
            output_dir_fastq.mkdir(exist_ok=True)
            logger.info(f"#####################################################################")
            logger.info(f"Running sepi_point on {sample_name}")
            logger.info(f"Using inputs {file_paths['r1_file']}, {file_paths['r2_file']}")
            logger.info(f"Printing output to {output_dir_fastq}")
            print(f"Running sepi_point on {sample_name}")
            vcf = sp.run_mapping_and_variant_calling(r1_file=file_paths["r1_file"], r2_file=file_paths["r2_file"],
                                                     output_dir=output_dir_fastq, output_prefix=sample_name,
                                                     reference_fasta=mutation_db_fasta,
                                                     no_clean=args.no_clean, logger=logger)
            sample_mutations = mf.get_mutations_from_vcf(vcf_file=vcf)
            sample_mutation_summary = mf.summarize_sample_mutations(sample_mutations=sample_mutations).copy()
            mf.print_sample_mutations(mutation_summary=sample_mutation_summary,summary_output_file=output_dir_fastq.joinpath("results.tsv"))
            all_sample_mutations[output_name_fastq] = sample_mutation_summary
        elif file_paths["assembly_file"] is not None:
            output_name_fasta = sample_name
            output_dir_fasta = args.output.joinpath(output_name_fasta)
            output_dir_fasta.mkdir(exist_ok=True)
            logger.info(f"#####################################################################")
            logger.info(f"Running sepi_point on {sample_name}")
            logger.info(f"Using input {file_paths['assembly_file']}")
            logger.info(f"Printing output to {output_dir_fasta}")
            print(f"Running sepi_point on {sample_name}")
            snp_file = sp.run_nucmer_and_showsnps(assembly_file=file_paths["assembly_file"],
                                                  output_dir=output_dir_fasta, output_prefix=sample_name,
                                                  reference_fasta=mutation_db_fasta, logger=logger)
            sample_mutations = mf.get_mutations_from_nucmer_snps(nucmer_snp_file=snp_file)
            sample_mutation_summary = mf.summarize_sample_mutations(sample_mutations=sample_mutations).copy()
            mf.print_sample_mutations(mutation_summary=sample_mutation_summary,summary_output_file=output_dir_fasta.joinpath("results.tsv"))
            all_sample_mutations[output_name_fasta] = sample_mutation_summary
            
    
    mf.print_sample_mutations_batch(mutation_summaries=all_sample_mutations,
                                    summary_output_file=args.output.joinpath("results.tsv"), matrix_output_file=args.output.joinpath("results.matrix.tsv"))
    

if __name__ == "__main__":
    main_cli()