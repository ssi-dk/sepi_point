import os
import sys
from sepi_point.seqtools import DnaSeq, ProteinSeq, NucleotideFasta, ProteinFasta, translate_dna
from pathlib import Path
import logging
import subprocess


class MutationFinder:

    def __init__(self):
        return(None)

    def load_and_check_db(self, mutation_db_tsv: Path, sequence_db_fasta: Path) -> None:
        self.tsv_path = os.path.abspath(mutation_db_tsv)
        self.fasta_path = os.path.abspath(sequence_db_fasta)
        self.sequences = NucleotideFasta.from_file(sequence_db_fasta)
        self.protein_sequences = self.sequences.translate()
        aa_mutation_dict = {}
        aa_mutation_list = []
        codon_mutation_dict = {}
        nt_mutation_dict = {}
        aa_to_codon = setup_aa_to_codon_table()
        with open(self.tsv_path) as f:
            firstline = True
            for line in f:
                line = line.rstrip('\n').split('\t')
                if firstline:
                    firstline = False
                    gene_idx = line.index("gene")
                    type_idx = line.index("type")
                    mutation_idx = line.index("mutation")
                    category_idx = line.index("category")
                else:
                    gene = line[gene_idx]
                    mutation = line[mutation_idx]
                    position = int(mutation[1:-1])
                    if line[type_idx] == "protein":
                        ref_aa = mutation[0]
                        if not self.protein_sequences[gene][0].sequence[position-1] == ref_aa:
                            print(f"Warning. Mutation tsv file contains mutation {gene}::{mutation}, but reference aa does not match fasta file")
                        alt_aa = mutation[-1]
                        alt_codons = aa_to_codon[alt_aa]
                        aa_mutation_list.append(gene+"::"+mutation)
                        if not gene in aa_mutation_dict:
                            aa_mutation_dict[gene] = {}
                            aa_mutation_dict[gene][str(position)] = [mutation]
                            codon_mutation_dict[gene] = {}
                            codon_mutation_dict[gene][str(position)] = alt_codons
                        elif not str(position) in aa_mutation_dict[gene]:
                            aa_mutation_dict[gene][str(position)] = [mutation]
                            codon_mutation_dict[gene][str(position)] = alt_codons
                        else:
                            aa_mutation_dict[gene][str(position)].append(mutation)
                            codon_mutation_dict[gene][str(position)] += alt_codons
                    else:
                        ref_nt = mutation[0]
                        if not self.sequences[gene][0].sequence[position-1] == ref_nt:
                            print(f"Warning. Mutation tsv file contains mutation {gene}::{mutation}, but reference nt does not match fasta file")
                        alt_nt = mutation[-1]
                        if not gene in nt_mutation_dict:
                            nt_mutation_dict[gene] = {}
                            nt_mutation_dict[gene][str(position)] = [alt_nt]
                        elif not str(position) in nt_mutation_dict[gene]:
                            nt_mutation_dict[gene][str(position)] = [alt_nt]
                        else:
                            nt_mutation_dict[gene][str(position)].append(alt_nt)

        self.aa_mutation_list = aa_mutation_list
        self.aa_mutation_dict = aa_mutation_dict
        self.codon_mutation_dict = codon_mutation_dict
        self.nt_mutation_dict = nt_mutation_dict
        return(None)

    def get_mutations_from_vcf(self, vcf_file: Path) -> dict:
        sample_mutation_dict = {}
        for gene in self.aa_mutation_dict:
            sample_mutation_dict[gene] = {}
        for gene in self.nt_mutation_dict:
            sample_mutation_dict[gene] = {}
        with open(vcf_file) as f:
            for line in f:
                if not line.startswith("##"):
                    if line.startswith("#"):
                        line = line.rstrip("\n").lstrip("#").split("\t")
                        CHROM_index = line.index("CHROM")
                        POS_index = line.index("POS")
                        REF_index = line.index("REF")
                        ALT_index = line.index("ALT")
                        FORMAT_index = line.index("FORMAT")
                        INFO_index = FORMAT_index + 1
                    else:
                        line = line.rstrip("\n").split("\t")
                        format = line[FORMAT_index].split(":")
                        info = line[INFO_index].split(":")
                        GT_index = format.index("GT")
                        if info[GT_index] == "1":
                            gene = line[CHROM_index]
                            pos = line[POS_index]
                            alt = line[ALT_index]
                            sample_mutation_dict[gene][pos] = alt
        return(sample_mutation_dict)

    def get_mutations_from_nucmer_snps(self, nucmer_snp_file: Path):
        sample_mutation_dict = {}
        for gene in self.nt_mutation_dict:
            sample_mutation_dict[gene] = {}
        with open(nucmer_snp_file) as f:
            for line in f:
                line = line.rstrip('\n').split('\t')
                gene = line[10]
                pos = line[0]
                alt = line[2]
                if not gene in sample_mutation_dict:
                    sample_mutation_dict[gene] = {}
                sample_mutation_dict[gene][pos] = alt
        return(sample_mutation_dict)

    def summarize_sample_mutations(self, sample_mutations: dict) -> dict:
        mutation_summary = {}
        for gene, alt_dict in self.nt_mutation_dict.items():
            for nt_position, alt_nts in alt_dict.items():
                if nt_position in sample_mutations[gene]:
                    for alt_nt in alt_nts:
                        if alt_nt in sample_mutations[gene][nt_position]:
                            ref_nt = self.sequences[gene][0].sequence[int(nt_position)-1]
                            nt_mut = ref_nt+nt_position+alt_nt
                            mutation_summary[gene+"::"+nt_mut] = [gene,nt_position,ref_nt,alt_nt,"",""]
        for gene, codon_dict in self.codon_mutation_dict.items():
            for aa_position, codons in codon_dict.items():
                start_position = int(aa_position)*3-2
                sample_codon = ""
                for position in range(start_position, start_position+3):
                    if gene in sample_mutations and str(position) in sample_mutations[gene]:
                        sample_codon += sample_mutations[gene][str(position)]
                    else:
                        sample_codon += self.sequences[gene][0].sequence[position-1]
                if sample_codon in codons:
                    ref_codon = self.sequences[gene][0].sequence[start_position-1:start_position+2]
                    ref_aa = translate_dna(ref_codon)
                    alt_aa = translate_dna(sample_codon)
                    aa_mut = ref_aa+aa_position+alt_aa
                    mutation_summary[gene+"::"+aa_mut] = [gene,aa_position,translate_dna(ref_codon),translate_dna(sample_codon),ref_codon,sample_codon]

        return(mutation_summary)

    @staticmethod
    def print_sample_mutations(mutation_summary: dict, summary_output_file: Path = None) -> None:
        print_header = ["Mutation","Gene","Position","Ref","Alt","Ref_codon","Alt_codon"]
        if summary_output_file is None:
            print("\t".join(print_header))
            for mutation, details in mutation_summary.items():
                printlist = [mutation]+details
                print("\t".join(printlist))
        else:
            o = open(summary_output_file,'w')
            o.write("\t".join(print_header)+"\n")
            for mutation, details in mutation_summary.items():
                printlist = [mutation]+details
                o.write("\t".join(printlist)+"\n")
            o.close()
        return(None)
    

    def print_sample_mutations_batch(self, mutation_summaries: dict[dict], summary_output_file: Path = None, matrix_output_file: Path = None) -> None:
        print_header = ["Sample","Mutation","Gene","Position","Ref","Alt","Ref_codon","Alt_codon"]
        matrix_header = ["Sample"]+self.aa_mutation_list
        o = open(summary_output_file,'w')
        o.write("\t".join(print_header)+"\n")
        om = open(matrix_output_file, 'w')
        om.write("\t".join(matrix_header)+"\n")
        for sample_name, mutation_summary in mutation_summaries.items():
            # write to long format output
            for mutation, details in mutation_summary.items():
                printlist = [sample_name, mutation]+details
                o.write("\t".join(printlist)+"\n")
            # write to presence/absence matrix output
            matrix_printlist = [sample_name]
            for mutation in self.aa_mutation_list:
                if mutation in mutation_summary:
                    matrix_printlist.append("1")
                else:
                    matrix_printlist.append("0")
            om.write("\t".join(matrix_printlist)+"\n")
        o.close()
        om.close()
        return(None)

    def iter_db_codons(self):
        for gene, codon_dict in self.codon_mutation_dict.items():
            for aa_position, codons in codon_dict.items():
                yield(gene, aa_position, codons)


def setup_aa_to_codon_table() -> dict:
    table = { 
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', 
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 
            '---': '-'
        }
    aa_to_codon_table = {}
    for codon, aa in table.items():
        if aa in aa_to_codon_table:
            aa_to_codon_table[aa].append(codon)
        else:
            aa_to_codon_table[aa] = [codon]
    return(aa_to_codon_table)

