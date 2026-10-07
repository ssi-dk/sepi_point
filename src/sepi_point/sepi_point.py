#!/usr/bin/env python3

import logging
logger = logging.getLogger(__name__)

import sys
from pathlib import Path
import subprocess
from sepi_point.mutation_finder import MutationFinder
from sepi_point.seqtools import WgsData
import argparse
import re
from importlib import resources

def parse_args(argv):
    parser = argparse.ArgumentParser(description='Run SepiPOINT on a single isolate')
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
    parser.add_argument("-l", "--log_level",
                        help = "Logging depth. Default: INFO",
                        type=str,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        default="INFO",
                        required = False)
    args = parser.parse_args()
    return args



def execute_cmd_and_log(cmd, log_stdout=True, log_stderr=True) -> tuple[str, str]:
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
        format="{asctime} {name} {levelname}: {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logger


def check_fastq_inputs(r1_file: Path, r2_file: Path, mutation_db_tsv: Path, mutation_db_fasta: Path):
    """"
    Check if all inputs and database files are present for paired end read input
    """
    all_files_found = True
    if not r1_file.is_file():
        all_files_found = False
        warning = f"No R1 input file found at {r1_file}"
        print(warning)
        logger.error(warning)
    if not r2_file.is_file():
        all_files_found = False
        warning = f"No R2 input file found at {r2_file}"
        print(warning)
        logger.error(warning)
    if not mutation_db_tsv.exists():
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


def check_fasta_inputs(assembly_file: Path, mutation_db_tsv: Path, mutation_db_fasta: Path):
    """"
    Check if all inputs and database files are present for assembly input
    """
    all_files_found = True
    if not assembly_file.is_file():
        all_files_found = False
        warning = f"No assembly file found at {assembly_file}"
        print(warning)
        logger.error(warning)
    if not mutation_db_tsv.exists():
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


def run_mapping_and_variant_calling(
        r1_file: Path, 
        r2_file: Path,
        output_dir: Path, 
        output_prefix: str, 
        reference_fasta: Path,
        no_clean: bool, 
        threads: int = 1,
        ) -> Path:
    """"
    Run read-mapping and variant calling on paired end read files.
    Return path to vcf file
    """

    prefix = output_dir.joinpath(output_prefix)
    sorted_bam = Path(f"{prefix}.sorted.bam")                     # sorted bam file filtered on q30 and only including mapped reads from primary alignments
    sorted_bam_idx = Path(f"{prefix}.sorted.bam.bai")             # sorted bam index file
    vcf = Path(f"{prefix}.vcf")                                   # variant calls

    # Run read mapping using bwa mem, discard unmapped reads (-F 4) and low quality 
    # mappings (-q 30) using samtools view and sort using samtools sort
    if not vcf.exists() and not sorted_bam_idx.exists() and not sorted_bam.exists():
        cmd = f"bwa mem -v 1 -t {threads} {reference_fasta} {r1_file} {r2_file} | \
            samtools view -b -F 4 -q 30 - | \
            samtools sort -o {sorted_bam} 2> /dev/null"
        stdout, stderr = execute_cmd_and_log(cmd=cmd)
    else:
        logger.info(f"Bam file found at {sorted_bam}. Skipping bwa mem read mapping and samtools view and sort.")

    # Index bam with samtools index
    if not vcf.exists() and not sorted_bam_idx.exists():
        cmd = f"samtools index -o {sorted_bam_idx} {sorted_bam}"
        stdout, stderr = execute_cmd_and_log(cmd=cmd)
    else:
        logger.info(f"Sorted, indexed bam file found at {sorted_bam_idx}.")

    # run bcftools mpileup including anomalous read pairs (-A) and max 10000 reads (-d), 
    # then call variants to generate vcf file using alternative model for multiallelic and 
    # rare-variant calling (-m) and output variant sites only (-v)
    if not vcf.exists():
        cmd = f"bcftools mpileup -A -Ou -d 10000 -f {reference_fasta} {sorted_bam} | \
            bcftools call -mv -Ov -o {vcf}"
        stdout, stderr = execute_cmd_and_log(cmd=cmd, log_stdout=False, log_stderr=False)
    else:
        logger.info(f"Vcf file found at {vcf}.")
    
    if not no_clean:
        if sorted_bam.exists():
            sorted_bam.unlink()
        if sorted_bam_idx.exists():
            sorted_bam_idx.unlink()
        logger.info(f"Cleaned up bam and sam files from output_folder")

    return(vcf)


def run_nucmer_and_showsnps(
        assembly_file: Path, 
        output_dir: Path,
        output_prefix: str,
        reference_fasta: Path,
        ) -> Path:
    """"
    Run nucmer and show-snps on assembled genome
    Return path to .snps file.
    """
    prefix = output_dir.joinpath(output_prefix)
    delta_file = Path(f"{prefix}.delta")
    snps_file = Path(f"{prefix}.snps")
    cmd = f"nucmer --minmatch 15 --prefix={prefix} {reference_fasta} {assembly_file}; show-snps -HrT {delta_file} > {snps_file}"
    if not snps_file.exists():
        stdout, stderr = execute_cmd_and_log(cmd=cmd, log_stdout=False, log_stderr=False)
    else:
        logger.info(f"Nucmer snps file found at {snps_file}.")
    return(snps_file)

def run_blast(assembly_file: Path,
              output_dir: Path, output_prefix: str,
              reference_fasta: Path) -> Path:
    """"
    Blast genes in reference fasta file
    """
    prefix = output_dir.joinpath(output_prefix)
    blast_output_tsv = Path(f"{prefix}.blast.tsv")
    cmd = f"blastn -query {reference_fasta} -subject {assembly_file} -qcov_hsp_perc 60 -perc_identity 80 -outfmt '6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send qseq sseq evalue bitscore' -out {blast_output_tsv}"
    if not blast_output_tsv.exists():
        stdout, stderr = execute_cmd_and_log(cmd=cmd, log_stdout=False, log_stderr=False)
    else:
        logger.info(f"Blast output file found at {blast_output_tsv}.")
    return(blast_output_tsv)


def run_on_reads(
        mf: MutationFinder,
        sample_name: str,
        r1_file: Path,
        r2_file: Path,
        output_dir: Path,
        no_clean:bool,
    ) -> dict:
    logger.info(f"#####################################################################")
    logger.info(f"Running SepiPOINT on {sample_name}")
    logger.info(f"Using inputs {r1_file}, {r2_file}")
    logger.info(f"Printing output to {output_dir}")
    print(f"Running SepiPOINT on {sample_name}")

    vcf = run_mapping_and_variant_calling(
        r1_file=r1_file,
        r2_file=r2_file,
        output_dir=output_dir,
        output_prefix=sample_name,
        reference_fasta=mf.fasta_path,
        no_clean=no_clean,
        )
    sample_mutations = mf.get_mutations_from_vcf(vcf_file=vcf)
    return(sample_mutations)

def run_on_assembly(
        mf: MutationFinder,
        sample_name: str,
        assembly_file: Path,
        output_dir: Path,
    ) -> dict:
    # all_files_found = check_fasta_inputs(assembly_file=assembly_file, mutation_db_tsv=mutation_db_tsv, mutation_db_fasta=mutation_db_fasta,logger=logger)
    # if all_files_found:
    #     snp_file = run_nucmer_and_showsnps(assembly_file=assembly_file,
    #                                         output_dir=output_dir, output_prefix=output_prefix,
    #                                         reference_fasta=mutation_db_fasta, logger=logger)
    #     mf = MutationFinder()
    #     mf.load_and_check_db(mutation_db_tsv=mutation_db_tsv,sequence_db_fasta=mutation_db_fasta)
    #     sample_mutations = mf.get_mutations_from_nucmer_snps(nucmer_snp_file=snp_file)
    # else:
    #     sys.exit()
    # return(mf, sample_mutations)
    logger.info("#####################################################################")
    logger.info("Running SepiPOINT on %s", sample_name)
    logger.info("Using input %s", assembly_file)
    logger.info("Printing output to %s", output_dir)
    print(f"Running SepiPOINT on {sample_name}")
    snp_file = run_nucmer_and_showsnps(
        assembly_file=assembly_file,
        output_dir=output_dir,
        output_prefix=sample_name,
        reference_fasta=mf.fasta_path,
        )
    sample_mutations = mf.get_mutations_from_nucmer_snps(nucmer_snp_file=snp_file)
    return(sample_mutations)

def run_sepi_point(
        mf: MutationFinder,
        wgs_data: WgsData,
        output_dir: Path,
        no_clean: bool = False,
) -> tuple[dict, dict]:
    """
    This function runs SepiPOINT on a WgsData object. Note that it assumes paths have been checked before.
    """
    all_sample_mutations= {}
    all_sample_putative_mutations = {}
    for sample_name, file_paths, info in wgs_data:
        if file_paths.get("r1_file") is not None and file_paths.get("r2_file") is not None and file_paths.get("assembly_file") is not None:
            output_name_fasta = f"{sample_name}_fasta"
            output_name_fastq = f"{sample_name}_fastq"
        else:
            output_name_fasta = sample_name
            output_name_fastq = sample_name

        results_file = f"{sample_name}.results.tsv"
        if file_paths.get("r1_file") is not None and file_paths.get("r2_file") is not None:
            output_dir_fastq = output_dir.joinpath(output_name_fastq)
            output_dir_fastq.mkdir(exist_ok=True)
            sample_mutations = run_on_reads(
                mf=mf,
                sample_name=sample_name,
                r1_file=file_paths["r1_file"],
                r2_file=file_paths["r2_file"],
                output_dir=output_dir_fastq,
                no_clean=no_clean,
            )
            sample_mutation_summary, sample_putative_mutation_summary = mf.summarize_sample_mutations(sample_mutations=sample_mutations)
            mf.print_sample_mutations(mutation_summary=sample_mutation_summary,summary_output_file=output_dir_fastq.joinpath(results_file))
            all_sample_mutations[output_name_fastq] = sample_mutation_summary
            all_sample_putative_mutations[sample_name] = sample_putative_mutation_summary
        if file_paths.get("assembly_file") is not None:
            output_dir_fasta = output_dir.joinpath(output_name_fasta)
            output_dir_fasta.mkdir(exist_ok=True)
            sample_mutations = run_on_assembly(
                mf=mf,
                sample_name=sample_name,
                assembly_file=file_paths["assembly_file"],
                output_dir=output_dir_fasta,
                )
            sample_mutation_summary, sample_putative_mutation_summary = mf.summarize_sample_mutations(sample_mutations=sample_mutations)
            mf.print_sample_mutations(mutation_summary=sample_mutation_summary,summary_output_file=output_dir_fasta.joinpath(results_file))
            all_sample_mutations[output_name_fasta] = sample_mutation_summary
            all_sample_putative_mutations[sample_name] = sample_putative_mutation_summary
        else:
            logger.info("No files found to run on for sample %s with file_paths: %s", sample_name, file_paths)

    return(all_sample_mutations, all_sample_putative_mutations)

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

    logger = setup_logger(log_file=log_file, log_level=args.log_level)
    mutation_db_tsv = resources.files("sepi_point").joinpath("db").joinpath("mutations.tsv")
    mutation_db_fasta = resources.files("sepi_point").joinpath("db").joinpath("sequences.fasta")

    logger.info("#### RUNNING SepiPOINT ####")
    logger.info("Checking for mutations found in %s", mutation_db_tsv)
    logger.info("Reference fasta file: %s", mutation_db_fasta)
    file_paths = {}
    if args.r1_file and args.r2_file:
        if check_fastq_inputs(
            r1_file=args.r1_file,
            r2_file=args.r2_file,
            mutation_db_tsv=mutation_db_tsv,
            mutation_db_fasta=mutation_db_fasta,
        ):
            logger.info("Using inputs %s, %s", args.r1_file, args.r2_file)
            file_paths["r1_file"] = str(args.r1_file)
            file_paths["r2_file"] = str(args.r2_file)
        else:
            print("Some files were not found. Check log at %s", log_file)
            sys.exit()
     
    elif args.assembly:
        if check_fasta_inputs(
            assembly_file=args.assembly,
            mutation_db_tsv=mutation_db_tsv,
            mutation_db_fasta=mutation_db_fasta,
        ):
            logger.info(f"Using input {args.assembly}")
            file_paths["assembly_file"] = str(args.assembly)
        else:
            print("Some files were not found. Check log at %s", log_file)
            sys.exit()
    else:
        logger.error("Input files must be provided. Either an assembled genome with -a option or paired end reads with -1 and -2 options.")
        print("Error: Input files must be provided. Either an assembled genome with -a option or paired end reads with -1 and -2 options. See epi_point.py -h for more info.")
        sys.exit()

    wgs_data = WgsData(file_paths={args.sample_name: file_paths})
    mf = MutationFinder()
    mf.load_and_check_db(mutation_db_tsv=mutation_db_tsv,sequence_db_fasta=mutation_db_fasta)

    run_sepi_point(
        mf = mf,
        wgs_data=wgs_data,
        output_dir=args.output,
        no_clean=args.no_clean,
    )

    print("Done.")
    logger.info("#### DONE RUNNING SepiPOINT. ####")


if __name__ == "__main__":
    main_cli()
    
