#!/usr/bin/env python3

import os
import sys
from pathlib import Path
import logging
import subprocess
from sepi_point.mutation_finder import MutationFinder
import argparse
import re
from importlib import resources


def parse_args(argv):
    parser = argparse.ArgumentParser(description='Run sepi_point on a single isolate')
    parser.add_argument("-1", "--r1_file",
                        help = "Forward Illumina read file (fastq or fastq.gz)",
                        type=Path,
                        required = False)
    parser.add_argument("-2", "--r2_file",
                        help = "Reverse Illumina read file (fastq or fastq.gz)",
                        type=Path,
                        required = False)
    parser.add_argument("-a", "--assembly",
                        help = "Genome assembly (fasta)",
                        type=Path,
                        required = False)
    parser.add_argument("-o", "--output",
                        help = "Output directory.",
                        type=Path,
                        required = True)
    parser.add_argument("-s", "--sample_name",
                        help = "Sample name (will be auto-detected from input files if not supplied)",
                        type=str,
                        required = False)
    parser.add_argument("-n", "--no_clean",
                        help = "Do not clean up sam and bam files after variant calling. Default False.",
                        action= "store_true",
                        default=False)
    args = parser.parse_args()
    return args



def execute_cmd_and_log(cmd, logger, log_stdout=True, log_stderr=True) -> tuple[str, str]:
    logger.info(f"Running command: {cmd}")
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        encoding="utf-8",
    )
    stdout, stderr = process.communicate()
    if log_stdout and stdout and stdout is not None:
        logger.info(f"Shell command STDOUT: {stdout}")
    if log_stderr and stderr and stderr is not None:
        logger.error(f"Shell command STDERR: {stderr}")
    return stdout, stderr

def setup_logger(log_file, log_level="INFO") -> logging.RootLogger:
    logger = logging.getLogger()
    logging.basicConfig(
        level=log_level,
        filename=str(log_file),
        encoding="utf-8",
        filemode="w",
        format="{asctime} - {levelname} - {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logger


def check_fastq_inputs(r1_file: Path, r2_file: Path, mutation_db_tsv: Path, mutation_db_fasta: Path, logger):
    """"
    Check if all inputs and database files are present for paired end read input
    """
    all_files_found = True
    if not r1_file.exists():
        all_files_found = False
        warning = f"No R1 input file found at {r1_file}"
        print(warning)
        logger.error(warning)
    if not r2_file.exists():
        all_files_found = False
        warning = f"No R2 input file found at {r2_file}"
        print(warning)
        logger.error(warning)
    if not mutation_db_fasta.exists():
        all_files_found = False
        warning = f"No mutation db tsv found at {mutation_db_tsv}"
        print(warning)
        logger.error(warning)
    if not mutation_db_fasta.exists():
        all_files_found = False
        warning = f"No mutation db fasta file found at {mutation_db_fasta}"
        print(warning)
        logger.error(warning)
    return(all_files_found)


def check_fasta_inputs(assembly_file: Path, mutation_db_tsv: Path, mutation_db_fasta: Path, logger):
    """"
    Check if all inputs and database files are present for assembly input
    """
    all_files_found = True
    if not assembly_file.exists():
        all_files_found = False
        warning = f"No assembly file found at {assembly_file}"
        print(warning)
        logger.error(warning)
    if not mutation_db_fasta.exists():
        all_files_found = False
        warning = f"No mutation db tsv found at {mutation_db_tsv}"
        print(warning)
        logger.error(warning)
    if not mutation_db_fasta.exists():
        all_files_found = False
        warning = f"No mutation db fasta file found at {mutation_db_fasta}"
        print(warning)
        logger.error(warning)
    return(all_files_found)


def run_mapping_and_variant_calling(r1_file: Path, r2_file: Path,
                                    output_dir: Path, output_prefix: str, 
                                    reference_fasta: Path,
                                    no_clean: bool, logger) -> Path:
    """"
    Run read-mapping and variant calling on paired end read files.
    Return path to vcf file
    """

    prefix = output_dir.joinpath(output_prefix)
    sam = Path(f"{prefix}.sam")                                   # sam file produced from bwa mem mapping to 
    bam = Path(f"{prefix}.bam")                                   # bam file filtered on q30 and only including mapped reads from primary alignments
    sorted_bam = Path(f"{prefix}.sorted.bam")                     # sorted bam file
    sorted_bam_idx = Path(f"{prefix}.sorted.bam.bai")             # sorted bam index file
    sorted_sam = Path(f"{prefix}.sorted.sam")                     # sorted sam file for parsing
    vcf = Path(f"{prefix}.vcf")                                   # variant calls

    if not vcf.exists() and not sorted_sam.exists() and not sorted_bam.exists() and not bam.exists() and  not sam.exists():
        cmd = f"bwa mem -v 1 -o {sam} {reference_fasta} {r1_file} {r2_file} 2> /dev/null"
        stdout, stderr = execute_cmd_and_log(cmd=cmd, logger=logger)
    else:
        logger.info(f"Sam file found at {sam}. Skipping bwa mem read mapping.")

    # Run samtools view to filter unmapped reads and convert to bam
    # -F 260 to only include mapped and exclude secondary alignments,
    if not vcf.exists() and not sorted_sam.exists() and not sorted_bam.exists() and not bam.exists():
        cmd = f"samtools view -q 30 -h -F 4 -O BAM -o {bam} {sam}"
        stdout, stderr = execute_cmd_and_log(cmd=cmd, logger=logger)
    else:
        logger.info(f"Bam file found at {bam}.")

    # Sort and index bam with samtools
    if not vcf.exists() and not sorted_sam.exists() and (not sorted_bam.exists() or not sorted_bam_idx.exists()):
        cmd = f"samtools sort {bam} -o {sorted_bam}; samtools index -o {sorted_bam_idx} {sorted_bam}"
        stdout, stderr = execute_cmd_and_log(cmd=cmd, logger=logger)
    else:
        logger.info(f"Sorted bam file found at {sorted_bam}.")


    # Run samtools view to convert bam to sam
    if not vcf.exists() and not sorted_sam.exists():
        cmd = f"samtools view -h -O SAM -o {sorted_sam} {sorted_bam}"
        stdout, stderr = execute_cmd_and_log(cmd=cmd, logger=logger)
    else:
        logger.info(f"Sorted sam file found at{sorted_sam}.")
    
    # run bcftools call to generate vcf
    if not vcf.exists():
        cmd = f"bcftools mpileup -A -f {reference_fasta} {sorted_sam} | bcftools call -p 0.5 --ploidy 2 -mv -Ov -o {vcf}"
        stdout, stderr = execute_cmd_and_log(cmd=cmd, logger=logger,log_stdout=False, log_stderr=False)
    else:
        logger.info(f"Vcf file found at {vcf}.")
    
    if not no_clean:
        if sam.exists():
            sam.unlink()
        if bam.exists():
            bam.unlink()
        if sorted_bam.exists():
            sorted_bam.unlink()
        if sorted_bam_idx.exists():
            sorted_bam_idx.unlink()
        if sorted_sam.exists():
            sorted_sam.unlink()
        logger.info(f"Cleaned up bam and sam files from output_folder")

    return(vcf)


def run_nucmer_and_showsnps(assembly_file: Path, 
                            output_dir: Path, output_prefix: str,
                            reference_fasta: Path, logger) -> Path:
    """"
    Run nucmer and show-snps on assembled genome
    Return path to .snps file.
    """
    prefix = output_dir.joinpath(output_prefix)
    delta_file = Path(f"{prefix}.delta")
    snps_file = Path(f"{prefix}.snps")
    cmd = f"nucmer --minmatch 15 --prefix={prefix} {reference_fasta} {assembly_file}; show-snps -HrT {delta_file} > {snps_file}"
    if not snps_file.exists():
        stdout, stderr = execute_cmd_and_log(cmd=cmd, logger=logger,log_stdout=False, log_stderr=False)
    else:
        logger.info(f"Nucmer snps file found at {snps_file}.")
    return(snps_file)

def run_blast(assembly_file: Path,
              output_dir: Path, output_prefix: str,
              reference_fasta: Path, logger) -> Path:
    """"
    Blast genes in reference fasta file
    """
    prefix = output_dir.joinpath(output_prefix)
    blast_output_tsv = Path(f"{prefix}.blast.tsv")
    cmd = f"blastn -query {reference_fasta} -subject {assembly_file} -qcov_hsp_perc 60 -perc_identity 80 -outfmt '6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send qseq sseq evalue bitscore' -out {blast_output_tsv}"
    if not blast_output_tsv.exists():
        stdout, stderr = execute_cmd_and_log(cmd=cmd, logger=logger,log_stdout=False, log_stderr=False)
    else:
        logger.info(f"Blast output file found at {blast_output_tsv}.")
    return(blast_output_tsv)


def run_on_reads(r1_file: Path, r2_file: Path, output_dir: Path, output_prefix: str,
                 mutation_db_tsv: Path, mutation_db_fasta: Path, no_clean:bool, logger):
    all_files_found = check_fastq_inputs(r1_file=r1_file, r2_file=r2_file, mutation_db_tsv=mutation_db_tsv, mutation_db_fasta=mutation_db_fasta, logger=logger)
    if all_files_found:
        vcf = run_mapping_and_variant_calling(r1_file=r1_file, r2_file=r2_file,
                                        output_dir=output_dir, output_prefix=output_prefix,
                                        reference_fasta=mutation_db_fasta,
                                        no_clean=no_clean, logger=logger)
        mf = MutationFinder()
        mf.load_and_check_db(mutation_db_tsv=mutation_db_tsv,sequence_db_fasta=mutation_db_fasta)
        sample_mutations = mf.get_mutations_from_vcf(vcf_file=vcf)
    else:
        print(f"Input or reference DB files not found. Exiting.")
        sys.exit()
    return(mf, sample_mutations)

def run_on_assembly(assembly_file: Path, output_dir: Path, output_prefix: str,
                 mutation_db_tsv: Path, mutation_db_fasta: Path, logger):
    all_files_found = check_fasta_inputs(assembly_file=assembly_file, mutation_db_tsv=mutation_db_tsv, mutation_db_fasta=mutation_db_fasta,logger=logger)
    if all_files_found:
        snp_file = run_nucmer_and_showsnps(assembly_file=assembly_file,
                                            output_dir=output_dir, output_prefix=output_prefix,
                                            reference_fasta=mutation_db_fasta, logger=logger)
        mf = MutationFinder()
        mf.load_and_check_db(mutation_db_tsv=mutation_db_tsv,sequence_db_fasta=mutation_db_fasta)
        sample_mutations = mf.get_mutations_from_nucmer_snps(nucmer_snp_file=snp_file)
    else:
        sys.exit()
    return(mf, sample_mutations)

def main_cli():
    args = parse_args(argv=sys.argv)

    ### Set up output directory
    try:
        args.output.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        print(f"Failed to setup output directory at {args.output}: {e}")
        sys.exit()

    ### Check if input files are provided and set sample name if not provided
    if args.r1_file and args.r2_file:
        if not args.sample_name:
            re_match = re.match(r"(?P<sample_name>.+?)(?P<sample_number>(_S[0-9]+)?)(?P<lane>(_L[0-9]+)?)[\._]R?(?P<paired_read_number>[1|2])(?P<set_number>(_[0-9]+)?)(?P<file_extension>\.fastq\.gz)", args.r1_file.name)
            args.sample_name = re_match.group("sample_name")
    elif args.assembly:
        if not args.sample_name:
            args.sample_name = args.assembly.stem
    else:
        print("Error: Input files must be provided. Either an assembled genome with -a option or paired end reads with -1 and -2 options. See epi_point.py -h for more info.")
        sys.exit()
    log_file = Path(args.output).joinpath(f"{args.sample_name}.log")
    results_file = Path(args.output).joinpath(f"{args.sample_name}.results.tsv")

    logger = setup_logger(log_file=log_file)
    mutation_db_tsv = resources.files("sepi_point").joinpath("db").joinpath("mutations.tsv")
    mutation_db_fasta = resources.files("sepi_point").joinpath("db").joinpath("sequences.fasta")

    logger.info(f"#### RUNNING SEPI_POINT ####")
    logger.info(f"Checking for mutations found in {mutation_db_tsv}")
    logger.info(f"Reference fasta file: {mutation_db_fasta}")
    if args.r1_file and args.r2_file:
        logger.info(f"Using inputs {args.r1_file}, {args.r2_file}")
        mf, sample_mutations = run_on_reads(r1_file=args.r1_file, r2_file=args.r2_file,
                                            output_dir=args.output, output_prefix=args.sample_name,
                                            mutation_db_tsv=mutation_db_tsv, mutation_db_fasta=mutation_db_fasta,
                                            no_clean=args.no_clean, logger=logger)
            
    elif args.assembly:
        logger.info(f"Using input {args.assembly}")
        mf, sample_mutations = run_on_assembly(assembly_file=args.assembly,
                                               output_dir=args.output, output_prefix=args.sample_name, 
                                               mutation_db_tsv=mutation_db_tsv, mutation_db_fasta=mutation_db_fasta,
                                               logger=logger)
    else:
        logger.error("Input files must be provided. Either an assembled genome with -a option or paired end reads with -1 and -2 options.")
        print("Error: Input files must be provided. Either an assembled genome with -a option or paired end reads with -1 and -2 options. See epi_point.py -h for more info.")
    
    sample_mutation_summary = mf.summarize_sample_mutations(sample_mutations=sample_mutations)
    mf.print_sample_mutations(mutation_summary=sample_mutation_summary,summary_output_file=results_file)


if __name__ == "__main__":
    main_cli()
    
