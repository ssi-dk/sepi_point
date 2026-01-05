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
        mutation_list = []
        codon_mutation_dict = {}
        nt_mutation_dict = {}
        indel_dict = {}
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
                    req_frequency_idx = line.index("req_frequency")
                else:
                    gene = line[gene_idx]
                    mutation = line[mutation_idx]
                    category = line[category_idx]
                    req_frequency = line[req_frequency_idx]
                    if line[type_idx] == "protein":
                        if mutation.endswith("del") or mutation.startswith("ins"):
                            if mutation.endswith("del"):
                                position = int(mutation[1:-3])
                            else:
                                position = int(mutation[3:])
                            if not gene in indel_dict:
                                indel_dict[gene] = {position: {"mutation": mutation, "category": category, "req_frequency": req_frequency}}
                            else:
                                indel_dict[gene][str(position)] = {"mutation": mutation, "category": category, "req_frequency": req_frequency}
                        else:
                            position = int(mutation[1:-1])
                            ref_aa = mutation[0]
                            if not self.protein_sequences[gene][0].sequence[position-1] == ref_aa:
                                print(f"Warning. Mutation tsv file contains mutation {gene}::{mutation}, but reference aa does not match fasta file")
                            alt_aa = mutation[-1]
                            alt_codons = aa_to_codon[alt_aa]
                            if not gene in codon_mutation_dict:
                                codon_mutation_dict[gene] = {}
                            if not str(position) in codon_mutation_dict[gene]:
                                codon_mutation_dict[gene][str(position)] = {}
                            for codon in alt_codons:
                                codon_mutation_dict[gene][str(position)][codon] = {"mutation": mutation, "ref": ref_aa, "alt": alt_aa, "category": category, "req_frequency": req_frequency}.copy()
                    else:
                        position = int(mutation[1:-1])
                        ref_nt = mutation[0]
                        if not self.sequences[gene][0].sequence[position-1] == ref_nt:
                            print(f"Warning. Mutation tsv file contains mutation {gene}::{mutation}, but reference nt does not match fasta file")
                        alt_nt = mutation[-1]
                        if not gene in nt_mutation_dict:
                            nt_mutation_dict[gene] = {str(position): {"mutation": mutation, "ref": ref_nt, "alt": alt_nt, "category": category, "req_frequency": req_frequency}}
                        else:
                            nt_mutation_dict[gene][str(position)] = {"mutation": mutation, "ref": ref_nt, "alt": alt_nt, "category": category, "req_frequency": req_frequency}
                    mutation_list.append(gene+"::"+mutation)
        self.mutation_list = mutation_list
        self.codon_mutation_dict = codon_mutation_dict
        self.nt_mutation_dict = nt_mutation_dict
        self.indel_dict = indel_dict
        return(None)

    def get_mutations_from_vcf(self, vcf_file: Path) -> dict:
        sample_mutation_dict = {}
        for gene in self.sequences.seq_names:
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
                        gene = line[CHROM_index]
                        if gene in sample_mutation_dict:
                            format = line[FORMAT_index].split(":")
                            info = line[INFO_index].split(":")
                            GT_index = format.index("GT")
                            if not info[GT_index] == "0/0":
                                pos = line[POS_index]
                                ref = line[REF_index]
                                alt = line[ALT_index]
                                AD_index = format.index("AD")
                                AD = info[AD_index].split(",")
                                sample_mutation_dict[gene][pos] = {alt: {"mutation": ref+pos+alt, "ref": ref, "alt_depth": int(AD[1]), "total_depth": int(AD[0])+int(AD[1])}}
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
                ref = line[1]
                alt = line[2]
                if not gene in sample_mutation_dict:
                    sample_mutation_dict[gene] = {}
                sample_mutation_dict[gene][pos] = {alt: {"mutation": ref+pos+alt, "ref": ref, "alt_depth": 1, "total_depth": 1}}
        return(sample_mutation_dict)

    def get_mutations_from_blast_tsv(self, blast_output_file: Path):
        blast_hit_dict = {}
        with open(blast_output_file) as f:
            for line in f:
                line = line.rstrip('\n').split('\t')
                gene = line[0]
                bitscore = int(line[14])
                if gene not in blast_hit_dict or bitscore > blast_hit_dict[gene][0]:
                    qseq = line[11]
                    sseq = line[12]
                    blast_hit_dict[gene] = line
        f.close()

        sample_mutation_dict = {}
        for gene in self.nt_mutation_dict:
            sample_mutation_dict[gene] = {}

        for gene, stats in blast_hit_dict.items():      
            qseqid = stats[0]
            qstart = stats[7]
            qseq = stats[11]
            sseq = stats[12]
            bitscore = stats[14]

            qstart = int(qstart)
            ref_pos = qstart

            for ref_base, alt_base in zip(qseq, sseq):

                if ref_base == "-" or alt_base == "-":
                    # skip indels — only substitutions desired
                    if ref_base != "-":  # reference advances only if ref_base is a nucleotide
                        ref_pos += 1
                    continue

                if ref_base.upper() != alt_base.upper():
                    sample_mutation_dict[qseqid].append(f"{ref_base}{ref_pos}{alt_base}")

                ref_pos += 1

        return sample_mutation_dict


    def summarize_sample_mutations(self, sample_mutations: dict) -> dict:
        mutation_summary = {}
        for gene, position_dict in self.nt_mutation_dict.items():
            for nt_position, nt_dict in position_dict.items():
                if gene in sample_mutations and nt_position in sample_mutations[gene]:
                    alt_nt = nt_dict["alt"]
                    if alt_nt in sample_mutations[gene][nt_position]:
                        info_dict = sample_mutations[gene][str(nt_position)]
                        nt = list(info_dict.keys())[0]
                        alt_depth = info_dict[nt]["alt_depth"]
                        total_depth = info_dict[nt]["total_depth"]
                        try:
                            alt_freq_req = float(nt_dict["req_frequency"])
                        except ValueError:
                            alt_freq_req = 0
                        if alt_depth/total_depth >= alt_freq_req:
                            ref_nt = nt_dict["ref"]
                            nt_mut = ref_nt+nt_position+alt_nt
                            mut_string = gene+"::"+nt_mut
                            category = nt_dict["category"]
                            mutation_summary[mut_string] = [gene,nt_position,ref_nt,alt_nt,"","",f"{alt_depth}/{total_depth}",category]
        for gene, position_dict in self.codon_mutation_dict.items():
            for aa_position, codon_dict in position_dict.items():
                start_position = int(aa_position)*3-2
                sample_codon = ""
                for position in range(start_position, start_position+3):
                    if gene in sample_mutations and str(position) in sample_mutations[gene]:
                        nt_dict = sample_mutations[gene][str(position)]
                        nt = list(nt_dict.keys())[0]
                        sample_codon += nt
                        alt_depth = nt_dict[nt]["alt_depth"]
                        total_depth = nt_dict[nt]["total_depth"]
                    else:
                        sample_codon += self.sequences[gene][0].sequence[position-1]
                if sample_codon in codon_dict:
                    try:
                        alt_freq_req = float(codon_dict[sample_codon]["req_frequency"])
                    except ValueError:
                        alt_freq_req = 0
                    if alt_depth/total_depth >= alt_freq_req:
                        ref_codon = self.sequences[gene][0].sequence[start_position-1:start_position+2]
                        ref_aa = codon_dict[sample_codon]["ref"]
                        alt_aa = codon_dict[sample_codon]["alt"]
                        aa_mut = codon_dict[sample_codon]["mutation"]
                        category = codon_dict[sample_codon]["category"]
                        mutation_summary[gene+"::"+aa_mut] = [gene,aa_position,ref_aa,alt_aa,ref_codon,sample_codon,f"{alt_depth}/{total_depth}",category]

        for gene, position_dict in self.indel_dict.items():
            for aa_position, info_dict in position_dict.items():
                nt_position = int(aa_position)*3-2
                check_start = nt_position - 6
                for nt_pos in range(check_start, nt_position):
                    if str(nt_pos) in sample_mutations[gene]:
                        nt_dict = sample_mutations[gene][str(nt_pos)]
                        nt = list(nt_dict.keys())[0]
                        alt_depth = nt_dict[nt]["alt_depth"]
                        total_depth = nt_dict[nt]["total_depth"]
                        try:
                            alt_freq_req = float(info_dict["req_frequency"])
                        except ValueError:
                            alt_freq_req = 0
                        ref = nt_dict[nt]["ref"]
                        if (not len(ref) == len(nt) or ref == "." or nt == ".") and alt_depth/total_depth >= alt_freq_req:
                            category = info_dict["category"]
                            aa_mut = info_dict["mutation"]
                            mutation_summary[gene+"::"+aa_mut] = [gene,str(aa_position),ref,nt,"","",f"{alt_depth}/{total_depth}",category]
        return(mutation_summary)

    @staticmethod
    def print_sample_mutations(mutation_summary: dict, summary_output_file: Path = None) -> None:
        print_header = ["Mutation","Gene","Position","Ref","Alt","Ref_codon","Alt_codon","Alt_frequency","Category"]
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
        print_header = ["Sample","Mutation","Gene","Position","Ref","Alt","Ref_codon","Alt_codon","Alt_frequency","Category"]
        matrix_header = ["Sample"]+self.mutation_list
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
            for mutation in self.mutation_list:
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

