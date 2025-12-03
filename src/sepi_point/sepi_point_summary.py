#!/usr/bin/env python3

import sys
from pathlib import Path
from sepi_point.mutation_finder import MutationFinder
from sepi_point.seqtools import WgsData
import argparse
from importlib import resources


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
    args = parser.parse_args()
    return args


def main_cli():
    args = parse_args(argv=sys.argv)
    mutation_db_tsv = resources.files("sepi_point").joinpath("db").joinpath("mutations.tsv")
    mutation_db_fasta = resources.files("sepi_point").joinpath("db").joinpath("sequences.fasta")

    ###
    mf = MutationFinder()
    mf.load_and_check_db(mutation_db_tsv=mutation_db_tsv,sequence_db_fasta=mutation_db_fasta)
    all_sample_mutations = {}
    for folder in args.results_dir.iterdir():
        if folder.is_dir():
            sample_name = folder.name
            snps_file = folder.joinpath(f"{sample_name}.snps")
            vcf_file = folder.joinpath(f"{sample_name}.vcf")
            if vcf_file.exists():
                sample_mutations = mf.get_mutations_from_vcf(vcf_file=vcf_file)
                sample_mutation_summary = mf.summarize_sample_mutations(sample_mutations=sample_mutations).copy()
                all_sample_mutations[sample_name] = sample_mutation_summary
            elif snps_file.exists():
                sample_mutations = mf.get_mutations_from_nucmer_snps(nucmer_snp_file=snps_file)
                sample_mutation_summary = mf.summarize_sample_mutations(sample_mutations=sample_mutations).copy()
                all_sample_mutations[sample_name] = sample_mutation_summary
    print(f"Summarizing sepi_point results from {len(all_sample_mutations)} samples")
    if args.output:
        summary_output_file = args.output.joinpath("results.tsv")
        matrix_output_file = args.output.joinpath("results.matrix.tsv")
        if not args.output.exists():
            args.output.mkdir()
    else:
        summary_output_file = args.results_dir.joinpath("results.tsv")
        matrix_output_file = args.results_dir.joinpath("results.matrix.tsv")
    mf.print_sample_mutations_batch(mutation_summaries=all_sample_mutations,
                                    summary_output_file=summary_output_file, matrix_output_file=matrix_output_file)
    


if __name__ == "__main__":
    main_cli()