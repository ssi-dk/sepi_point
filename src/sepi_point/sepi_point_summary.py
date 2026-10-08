#!/usr/bin/env python3
import logging
logger = logging.getLogger(__name__)

import sys
from pathlib import Path
from sepi_point.mutation_finder import MutationFinder
from sepi_point.sepi_point import setup_logger
from importlib import resources
import argparse
import logging

def parse_args(argv):
    parser = argparse.ArgumentParser(description='Summarize sepi_point results')
    parser.add_argument("-r", "--results_dir",
                        help = "Folder with output folders from sepi_point runs.",
                        type=Path,
                        required = True)
    parser.add_argument("-o", "--output",
                        help = "Output Folder. Default <results_dir>",
                        type=Path,
                        required = False)
    parser.add_argument("-m", "--mutation_tsv",
                        help = "A custom .tsv file with mutations to be called. Must agree with the reference fasta.",
                        type=Path,
                        required = False)
    parser.add_argument("-f", "--mutation_fasta",
                        help = "A custom reference fasta file. Must agree with the mutation db .tsv.",
                        type=Path,
                        required = False)
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
    logger = setup_logger(log_file='EepiPOINTsummary.log', log_level=args.log_level)
    mutation_db_tsv = args.mutation_tsv if args.mutation_tsv else resources.files("sepi_point").joinpath("db").joinpath("mutations.tsv")
    mutation_db_fasta = args.mutation_fasta if args.mutation_fasta else resources.files("sepi_point").joinpath("db").joinpath("sequences.fasta")
    ###
    mf = MutationFinder()
    mf.load_and_check_db(mutation_db_tsv=mutation_db_tsv,sequence_db_fasta=mutation_db_fasta)
    all_sample_mutations = {}
    all_sample_putative_mutations = {}
    logger.debug("SepiPOINT summary lookiung for samples in '%s' ...", args.results_dir)
    for folder in args.results_dir.iterdir():
        if folder.is_dir():
            sample_name = folder.name
            mf.sample_name = sample_name
            snps_file = folder.joinpath(f"{sample_name}.snps")
            vcf_file = folder.joinpath(f"{sample_name}.vcf")
            if vcf_file.is_file():
                sample_mutations = mf.get_mutations_from_vcf(vcf_file=vcf_file)
            elif snps_file.is_file():
                sample_mutations = mf.get_mutations_from_nucmer_snps(nucmer_snp_file=snps_file)
            else:
                logger.debug("No .vcf or .snps file found for folder '%s'", folder.absolute())
                continue
            sample_mutation_summary, sample_putative_mutation_summary = mf.summarize_sample_mutations(sample_mutations=sample_mutations)
            all_sample_mutations[sample_name] = sample_mutation_summary
            all_sample_putative_mutations[sample_name] = sample_putative_mutation_summary
    logger.debug("Summarized SepiPOINT results from %s samples. Now writing these to file or stdout...", len(all_sample_mutations))
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
    logger.info("SepiPOINT summary done for %s samples.", len(all_sample_mutations))


if __name__ == "__main__":
    main_cli()