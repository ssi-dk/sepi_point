import os
import sys
from pathlib import Path
import pandas
import re
import numpy

class BioSeq:
    def __init__(self, name: str,sequence: str) -> None:
        self.name = name
        self.sequence = sequence.upper()

    def __len__(self):
        return(len(self.sequence))

    def __iter__(self):
        for nt in self.sequence:
            yield nt
        
    
    def __repr__(self):
        return f"< Biological sequence object {self.name} of length {len(self)} bp >"
    
    def __add__(self,other):
        return(BioSeq(self.name,self.sequence+other.sequence))

    @classmethod
    def from_fasta(cls, input_file: Path):
        startflag = True
        with open(input_file) as f:
            for line in f:
                line = line.rstrip('\n')
                if line[0] == ">":
                    if startflag:
                        seq = ""
                        startflag = False
                    else:
                        print(f"Fasta file {input_file} contains more than one entry. Only first entry is loaded into DnaSeq object.\nTo load all entries into Fasta object use Fasta.from_file({input_file})\n", file = sys.stderr)
                        f.close()
                        return(cls(sequence_name, seq))
                    sequence_name = line[1:]

                else:
                    seq += line
        f.close()
        return cls(sequence_name, seq)
    

    
    def print_fasta(self,file: Path,linelength: int = None) -> None:
        printlines = []
        if linelength is None:
            printlines += [">"+self.name,self.sequence]
        else:
            printlines += [">"+self.name] +  [(self.sequence[i:i+linelength]) for i in range(0, len(self.sequence), linelength)]
        
        o = open(file,"w")
        o.write("\n".join(printlines)+"\n")
        o.close
        return None
    


class DnaSeq(BioSeq):
    def __init__(self, name: str,sequence: str, phred_scores: str = None) -> None:
        self.name = name
        self.sequence = sequence.upper()
        self.phred_scores = phred_scores
    
    
    def __repr__(self):
        return f"< DNA sequence object {self.name} of length {len(self)} bp >"
    
    def reverse_complement(self):
        complement = {'A': 'T', 'T': 'A',
                      'C': 'G', 'G': 'C',
                      'R': 'Y', 'Y': 'R',
                      'S': 'W', 'W': 'S',
                      'K': 'M', 'M': 'K',
                      'B': 'V', 'V': 'B',
                      'D': 'H', 'H': 'D',
                      'N': 'N', 'X': 'X',
                      '-': '-', '.': '.'}
        reverse_complement = "".join(complement.get(base, base) for base in self.sequence[::-1])
        illegal_characters = set(self.sequence)-set(complement.keys())
        if len(illegal_characters) > 0:
            print("Non standard UIPAC characters found in sequence. Replacing with N when reverse complementing.", file=sys.stderr)
            for c in illegal_characters:
                reverse_complement = reverse_complement.replace(c,"N")
        if self.phred_scores is None:
            return(DnaSeq(self.name,reverse_complement))
        else:
            return(DnaSeq(self.name,reverse_complement,self.phred_scores[::-1]))
    
    def translate(self):
        seq = self.sequence.replace("-","")
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
        } 
        protein_seq = "" 
        mod = len(seq)%3
        #if not mod == 0:
        #    print(f"Warning: number of nucleotides in sequence is not divisible by 3. Ignoring {mod} nt at the end of sequence when translating to protein.", file=sys.stderr)
        illegal_codons = 0
        for i in range(0, len(seq)-mod, 3):
            codon = seq[i:i + 3]
            if len(set(codon)-{"A","C","G","T"}) > 0:
                illegal_codons += 1
                protein_seq += "X"
            else:
                protein_seq+= table[codon]
        if illegal_codons > 0:
            print(f"WARNING: {illegal_codons} codons containing non A/C/G/T found in {self.name}. X has been inserted on those positions.", file=sys.stderr)
        return(ProteinSeq(self.name,protein_seq))
    
    
    def __add__(self,other):
        if self.phred_scores is None or other.phred_scores is None:
            return(DnaSeq(self.name,self.sequence+other.sequence))
        else:
            return(DnaSeq(self.name,self.sequence+other.sequence,self.phred_scores,other.phred_scores))

    @property
    def iupac_counts(self):
        base_types = ["A","C","G","T"]
        amb_types = ["R","Y","S","W","K","M","B","D","H","M"]
        gap_types = [".","-"]
        nt_counts = {"A": 0,"C": 0,"G": 0,"T": 0, "N": 0, 
                     "R": 0,"Y": 0,"S": 0,"W": 0,"K": 0,"M": 0,"B": 0,"D": 0,"H": 0,"M": 0,
                     ".": 0,"-": 0,
                     "ACGT": 0, "amb": 0, "N": 0, "gap": 0, "nonIUPAC": 0}
        for nt in self.sequence:
            try:
                nt_counts[nt] += 1
            except KeyError:
                nt_counts["nonIUPAC"] += 1
            if nt in base_types:
                nt_counts["ACGT"] += 1
            elif nt in amb_types:
                nt_counts["amb"] += 1
            elif nt in gap_types:
                nt_counts["gap"] += 1
            
                

        return(nt_counts)
        


    @staticmethod
    def phred_score_to_int(phred_scores: str) -> list:
        int_scores = []
        for score in phred_scores:
            int_scores.append(ord(score)-33)
        return(int_scores)
        


class ProteinSeq(BioSeq):
    def __init__(self, sequence_name: str, sequence: str) -> None:
        self.name = sequence_name
        self.sequence = sequence.upper()
    

    def __repr__(self):
        return f"< Protein sequence object {self.name} of length {len(self)} AA >"
    
    def __add__(self, other):
        return(ProteinSeq(self.name, self.sequence+other.sequence))


    @classmethod
    def from_fasta(cls, input_file: Path):
        startflag = True
        with open(input_file) as f:
            for line in f:
                line = line.rstrip('\n')
                if line[0] == ">":
                    if startflag:
                        seq = ""
                        startflag = False
                    else:
                        f.close()
                        print(f"Fasta file {input_file} contains more than one entry. Only first entry is loaded into DnaSeq object. To load all entries use Fasta.from_file()\n", file = sys.stderr)
                        return(cls(sequence_name, seq))
                    sequence_name = line[1:]

                else:
                    seq += line
        f.close()
        return cls(sequence_name, seq)
    

    def print_fasta(self, output_file: str, linelength: int = None) -> None:
        printlines = []
        if linelength is None:
            printlines += [">"+self.name,self.sequence]
        else:
            printlines += [">"+self.name] +  [(self.sequence[i:i+linelength]) for i in range(0, len(self.sequence), linelength)]
        
        o = open(output_file,"w")
        o.write("\n".join(printlines)+"\n")
        o.close
        return None


class Fasta:

    def __init__(self, bio_seqs: list[BioSeq]) -> None:
        self.entries = bio_seqs

        ### Indexing for entries to allow for easy retrieval by sequence name, i.e. Fasta["sequence_name"] = [DnaSeq object]
        self.name_index = {}
        for i, bio_seq in enumerate(self.entries):
            if bio_seq.name in self.name_index:
                self.name_index[bio_seq.name].append(i)
            else:
                self.name_index[bio_seq.name] = [i]
    
    # Alternative constructor that initiates from a fasta file
    @classmethod
    def from_file(cls, input_file: str):
        entries = "Start"
        skipped_lines_count = 0
        with open(input_file) as f:
            for line in f:
                line = line.rstrip('\n')
                if not len(line) == 0:
                    if line[0] == ">":
                        if entries == "Start":
                            entries = []
                        else:
                            entries.append(BioSeq(sequence_name,new_seq))
                        new_seq = ""
                        sequence_name = line[1:]

                    else:
                        new_seq += line
                else:
                    skipped_lines_count += 1
        entries.append(BioSeq(sequence_name,new_seq))
        f.close()
        if skipped_lines_count > 0:
            print(f"Warning: {skipped_lines_count} blank lines were skipped when loading fasta file")
        return cls(entries)
    
    def check_sequence_type(self):
        
        return

    def __iter__(self):
        for entry in self.entries:
            yield entry
    
    def items(self) -> tuple[str, str]:
        for entry in self.entries:
            yield entry.name, entry.sequence

    def iterseqs(self) -> str:
        for entry in self.entries:
            yield entry.sequence
    
    def __getitem__(self,name) -> BioSeq:
        return [ self.entries[i] for i in self.name_index[name] ]
    
    # len returns total number of nucleotides in fasta file
    def __len__(self):
        length = 0
        for sequence in self.iterseqs():
            length += len(sequence)
        return(length)
    
    def __repr__(self):
        return f"< Fasta object containing {len(self.entries)} sequences, total length {len(self)} bp >"
    
    @property
    def seq_names(self):
        names = []
        for entry in self:
            names.append(entry.name)
        return(names)

            
    def write(self, output_file: str = None, linelength: int = None) -> None:
        printlines = []
        if linelength is None:
            for name, sequence in self.items():
                printlines += [">"+name,sequence]
        else:
            for name, sequence in self.items():
                printlines += [">"+name] +  [(sequence[i:i+linelength]) for i in range(0, len(sequence), linelength)]
        if output_file is not None:
            o = open(output_file,"w")
            o.write("\n".join(printlines)+"\n")
            o.close
        else:
            for line in printlines:
                print(line)
        return None

    # Append a new BioSeq object to exisisting Fasta object
    def append(self, entry: BioSeq) -> None:
        self.entries.append(entry)
        if entry.name in self.name_index:
            self.name_index[entry.name].append(len(self.entries)-1)
        else:
            self.name_index[entry.name] = [len(self.entries)-1]
    
    ### Add all entries from another Fasta object to entries in self
    def add(self, other: BioSeq):
        for entry in other:
            self.append(entry)
        return None

    def concat_seq(self):
        combined_seq = ""
        for name, sequence in self.items():
            combined_seq += sequence
        return(combined_seq)


    @property
    def sequence_names(self):
        seq_names = []
        for entry in self.entries:
            seq_names.append(entry.name)
        return seq_names
    
    @property
    def sequences(self) -> BioSeq:
        seqs = []
        for entry in self.entries:
            seqs.append(entry.sequence)
        return seqs

    def uniquify(self) -> BioSeq:
        seq_to_names = {}
        for name, sequence in self.items():
            if sequence in seq_to_names:
                seq_to_names[sequence].append(name)
            else:
                seq_to_names[sequence] = [name]
        ordered_seqs = [k for k, v in sorted(seq_to_names.items(), key=lambda item: len(item[1]))]
        count_str_length = len(str(len(ordered_seqs)))
        seq_count = 0
        seq_to_unique_name = {}
        unique_fasta = Fasta([]) ### initialize empty Fasta object
        for sequence in ordered_seqs:
            seq_count += 1
            seq_count_str = "0"*(count_str_length-len(str(seq_count)))+str(seq_count)
            seq_unique_name = f"unique_sequence_{seq_count_str}"
            unique_fasta.append(BioSeq(name=seq_unique_name,sequence=sequence))
            seq_to_unique_name[sequence] = seq_unique_name
        
        seq_info = []
        for name, sequence in self.items():
            seq_info.append( (name, seq_to_unique_name[sequence], len(seq_to_names[sequence]), sequence) )
        
        unique_fasta.uniquify_info = seq_info
        return unique_fasta


class NucleotideFasta(Fasta):


    def __init__(self, dna_seqs: list[DnaSeq]) -> None:
        self.entries = dna_seqs

        ### Indexing for entries to allow for easy retrieval by sequence name, i.e. Fasta["sequence_name"] = [DnaSeq object]
        self.name_index = {}
        for i, dna_seq in enumerate(self.entries):
            if dna_seq.name in self.name_index:
                self.name_index[dna_seq.name].append(i)
            else:
                self.name_index[dna_seq.name] = [i]
    
    # Alternative constructor that initiates from a fasta file
    @classmethod
    def from_file(cls, input_file: str):
        entries = "Start"
        skipped_lines_count = 0
        with open(input_file) as f:
            for line in f:
                line = line.rstrip('\n')
                if not len(line) == 0:
                    if line[0] == ">":
                        if entries == "Start":
                            entries = []
                        else:
                            entries.append(DnaSeq(sequence_name,new_seq))
                        new_seq = ""
                        sequence_name = line[1:]

                    else:
                        new_seq += line
                else:
                    skipped_lines_count += 1
        entries.append(DnaSeq(sequence_name,new_seq))
        f.close()
        if skipped_lines_count > 0:
            print(f"Warning: {skipped_lines_count} blank lines were skipped when loading fasta file")
        return cls(entries)
    
    def __repr__(self):
        return f"< NucleotideFasta object containing {len(self.entries)} sequences, total length {len(self)} bp >"
        
    def uniquify(self) -> DnaSeq:
        seq_to_names = {}
        for name, sequence in self.items():
            if sequence in seq_to_names:
                seq_to_names[sequence].append(name)
            else:
                seq_to_names[sequence] = [name]
        ordered_seqs = [k for k, v in sorted(seq_to_names.items(), key=lambda item: len(item[1]))]
        count_str_length = len(str(len(ordered_seqs)))
        seq_count = 0
        seq_to_unique_name = {}
        unique_fasta = Fasta([]) ### initialize empty Fasta object
        for sequence in ordered_seqs:
            seq_count += 1
            seq_count_str = "0"*(count_str_length-len(str(seq_count)))+str(seq_count)
            seq_unique_name = f"unique_sequence_{seq_count_str}"
            unique_fasta.append(DnaSeq(name=seq_unique_name,sequence=sequence))
            seq_to_unique_name[sequence] = seq_unique_name
        
        seq_info = []
        for name, sequence in self.items():
            seq_info.append( (name, seq_to_unique_name[sequence], len(seq_to_names[sequence]), sequence) )
        
        unique_fasta.uniquify_info = seq_info
        return unique_fasta

    def translate(self):
        protein_fasta = ProteinFasta([])
        for entry in self:
            protein_fasta.append(entry.translate())
        return protein_fasta


    @property
    def iupac_counts(self):
        combined_seq = self.concat_seq()
        base_types = ["A","C","G","T"]
        amb_types = ["R","Y","S","W","K","M","B","D","H","M"]
        gap_types = [".","-"]
        nt_counts = {"A": 0,"C": 0,"G": 0,"T": 0, "N": 0, 
                     "R": 0,"Y": 0,"S": 0,"W": 0,"K": 0,"M": 0,"B": 0,"D": 0,"H": 0,"M": 0,
                     ".": 0,"-": 0,
                     "ACGT": 0, "amb": 0, "N": 0, "gap": 0, "nonIUPAC": 0}
        for nt in combined_seq:
            try:
                nt_counts[nt] += 1
            except KeyError:
                nt_counts["nonIUPAC"] += 1
            if nt in base_types:
                nt_counts["ACGT"] += 1
            elif nt in amb_types:
                nt_counts["amb"] += 1
            elif nt in gap_types:
                nt_counts["gap"] += 1

        return(nt_counts)
    
    def batch_iupac_counts(input_directory):
        input_directory = Path(input_directory)
        iupac_counts_dict = {}
        for p in input_directory.iterdir():
            if p.name.endswith( (".fna", ".fa", ".fasta") ):
                fasta = NucleotideFasta.from_file(p)
                iupac_counts_dict[p.name] = fasta.iupac_counts
        return(iupac_counts_dict)



class ProteinFasta(Fasta):

    def __init__(self, protein_seqs: list[ProteinSeq]) -> None:
        self.entries = protein_seqs

        ### Indexing for entries to allow for easy retrieval by sequence name, i.e. ProteinFasta["sequence_name"] = ProteinSeq object
        self.name_index = {}
        for i, seq in enumerate(self.entries):
            self.name_index[seq.name] = i
    
    # Alternative constructor that initiates from a fasta file
    @classmethod
    def from_file(cls, input_file: str):
        entries = "Start"
        skipped_lines_count = 0
        with open(input_file) as f:
            for line in f:
                line = line.rstrip('\n')
                if not len(line) == 0:
                    if line[0] == ">":
                        if entries == "Start":
                            entries = []
                        else:
                            entries.append(ProteinSeq(sequence_name,new_seq))
                        new_seq = ""
                        sequence_name = line[1:]

                    else:
                        new_seq += line
                else:
                    skipped_lines_count += 1
        entries.append(ProteinSeq(sequence_name,new_seq))
        f.close()
        if skipped_lines_count > 0:
            print(f"Warning: {skipped_lines_count} blank lines were skipped when loading fasta file")
        return cls(entries)

    
    # Alternative constructor that initiates from a fasta file
    @classmethod
    def from_file(cls, input_file: str):
        entries = "Start"
        with open(input_file) as f:
            for line in f:
                line = line.rstrip('\n')
                if line[0] == ">":
                    if entries == "Start":
                        entries = []
                    else:
                        entries.append(ProteinSeq(sequence_name,new_seq))
                    new_seq = ""
                    sequence_name = line[1:]

                else:
                    new_seq += line
        entries.append(ProteinSeq(sequence_name,new_seq))
        f.close()
        return cls(entries)
    
    def __repr__(self):
        return f"< ProteinFasta object containing {len(self.entries)} sequences, total length {len(self)} bp >"

    
    def uniquify(self) -> ProteinSeq:
        seq_to_names = {}
        for name, sequence in self.items():
            if sequence in seq_to_names:
                seq_to_names[sequence].append(name)
            else:
                seq_to_names[sequence] = [name]
        ordered_seqs = [k for k, v in sorted(seq_to_names.items(), key=lambda item: len(item[1]))]
        count_str_length = len(str(len(ordered_seqs)))
        seq_count = 0
        seq_to_unique_name = {}
        unique_fasta = ProteinFasta([]) ### initialize empty Fasta object
        for sequence in ordered_seqs:
            seq_count += 1
            seq_count_str = "0"*(count_str_length-len(str(seq_count)))+str(seq_count)
            seq_unique_name = f"unique_sequence_{seq_count_str}"
            unique_fasta.append(ProteinSeq(name=seq_unique_name,sequence=sequence))
            seq_to_unique_name[sequence] = seq_unique_name
        
        seq_info = []
        for name, sequence in self.items():
            seq_info.append( (name, seq_to_unique_name[sequence], len(seq_to_names[sequence]), sequence) )
        
        unique_fasta.uniquify_info = seq_info
        return unique_fasta


    def __getitem__(self,name) -> ProteinSeq:
        return [ self.entries[i] for i in self.name_index[name] ]
    
def translate_dna(dna_string: str):
    seq = dna_string.replace("-","")
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
    } 
    protein_seq = "" 
    mod = len(seq)%3
    if not mod == 0:
        print(f"Warning: number of nucleotides in sequence is not divisible by 3. Ignoring {mod} nt at the end of sequence when translating to protein.", file=sys.stderr)
    illegal_codons = 0
    for i in range(0, len(seq)-mod, 3):
        codon = seq[i:i + 3]
        if len(set(codon)-{"A","C","G","T"}) > 0:
            illegal_codons += 1
            protein_seq += "X"
        else:
            protein_seq+= table[codon]
    if illegal_codons > 0:
        print(f"WARNING: {illegal_codons} codons containing non A/C/G/T found. X has been inserted on those positions.", file=sys.stderr)
    return(protein_seq)



class WgsData:

    def __init__(self, file_paths: dict[dict], sample_data: dict[dict] = None):
        """
        Default initialization is from file paths provided in a dictionary with sample name as keys, like this:
        file_paths =
        {<sample_name>: {
            "assembly_file": <assembly_file_path>,
            "r1_file": <r1_file_path>,
            "r2_file": <r2_file_path>
            }
        }
        and metadata in a similar format with sample names as keys and a dictionary of attribute: value pairs as values.
        {<sample_name>: {
            "serotype": "a",
            "MLST": "15",
            "source": blood,
            ....
            }
        }



        To initialize from folders without using a metadata sheet, use WgsData.from_folders(<assembly_data_folder>, <read_data_folder>)
        Both assembly_data_folder and read_data_folder can be set to None if you have only reads or only assemblies.
        """
        self.sample_names = file_paths.keys()
        self.file_paths = file_paths
        self.paired_end_read_files = []
        self.assembly_files = []
        for files in file_paths.values():
            if "r1_file" in files and "r2_file" in files and files["r1_file"] is not None and files["r2_file"]:
                self.paired_end_read_files.append((files["r1_file"], files["r2_file"]))
            else:
                self.paired_end_read_files.append((None, None))
            if "assembly_file" in files and files["assembly_file"] is not None:
                self.assembly_files.append(files["assembly_file"])
            else:
                self.assembly_files.append(None)
        self.sample_data = sample_data
        if sample_data:
            self.sample_data_df = pandas.DataFrame.from_dict(
                sample_data, orient="index"
            )
            self.sample_data_df.index = self.sample_names
            self.sample_data_columns = list(self.sample_data_df.columns.values)
        else:
            self.sample_data_columns = None
            self.sample_data_df = None

    @classmethod
    def from_samplesheet(cls, sample_sheet: Path, base_folder: Path = None):
        """
        Load data from samplesheet into dictionary.
        Data sheet must contain column names "sample_name", "assembly_file", "r1_file", "r2_file" (if files exist)
        If there are duplicates in "sample_name" column, entries after the first instance will have added _2, _3, _4 etc to their name
        """

        df = pandas.read_csv(sample_sheet, sep="\t")
        file_paths_df = df[["sample_name", "assembly_file", "r1_file", "r2_file"]]
        file_paths = {}
        duplicate_sample_count = 0
        ignored_assembly_count = 0
        ignored_r1_count = 0
        ignored_r2_count = 0
        uniq_sample_names = []
        for index, row in file_paths_df.iterrows():
            sample_name = row["sample_name"]
            if not sample_name in uniq_sample_names:
                file_paths[sample_name] = row.to_dict()
                uniq_sample_names.append(sample_name)
            else:
                name_counter = 1
                sample_name_uniq = f"{sample_name}_{name_counter}"
                while sample_name_uniq in uniq_sample_names:
                    name_counter += 1
                    sample_name_uniq = f"{sample_name}_{name_counter}"
                file_paths[sample_name_uniq] = row.to_dict()
                uniq_sample_names.append(sample_name_uniq)
                duplicate_sample_count += 1
        if base_folder:
            base_folder = Path(base_folder).absolute()
            for files in file_paths.values():
                files["assembly_file"] = base_folder.joinpath(files["assembly_file"])
                files["r1_file"] = base_folder.joinpath(files["r1_file"])
                files["r2_file"] = base_folder.joinpath(files["r2_file"])
        else:
            base_folder = Path(sample_sheet).parent
            for sample, files in file_paths.items():
                if files["assembly_file"] is not None:
                    try:
                        files["assembly_file"] = base_folder.joinpath(
                            files["assembly_file"]
                        ).resolve()
                    except TypeError as e:
                        ignored_assembly_count += 1
                        files["assembly_file"] = None
                if files["r1_file"] is not None:
                    try:
                        files["r1_file"] = base_folder.joinpath(
                            files["r1_file"]
                        ).resolve()
                    except TypeError as e:
                        ignored_r1_count += 1
                        files["r1_file"] = None
                if files["r2_file"] is not None:
                    try:
                        files["r2_file"] = base_folder.joinpath(
                            files["r2_file"]
                        ).resolve()
                    except TypeError as e:
                        ignored_r2_count += 1
                        files["r2_file"] = None
        if ignored_assembly_count > 0:
            print(f"{ignored_assembly_count} assembly file paths ignored")
        if ignored_r1_count > 0:
            print(f"{ignored_r1_count} assembly file paths ignored")
        if ignored_assembly_count > 0:
            print(f"{ignored_r2_count} assembly file paths ignored")
        df["sample_name_uniq"] = uniq_sample_names
        df.set_index("sample_name_uniq", inplace=True)
        sample_data = df.to_dict("index")
        return cls(file_paths=file_paths, sample_data=sample_data)

    @classmethod
    def from_folders(
        cls,
        assembly_data_folder: Path = None,
        paired_end_read_data_folder: Path = None,
        fasta_file_pattern=r"(?P<sample_name>.+?)(\.fa|\.fna|\.fasta)$",
        paired_end_reads_pattern=r"(?P<sample_name>.+?)(?P<sample_number>(_S[0-9]+)?)(?P<lane>(_L[0-9]+)?)[\._]R?(?P<paired_read_number>[1|2])(?P<set_number>(_[0-9]+)?)(?P<file_extension>\.fastq\.gz)",
    ):
        if assembly_data_folder is not None:
            assembly_data_dict = cls.parse_folder_for_fasta_files(
                data_folder=assembly_data_folder, fasta_file_pattern=fasta_file_pattern
            )
        if assembly_data_folder is not None and len(assembly_data_dict) > 0:
            assembly_df = pandas.DataFrame.from_dict(assembly_data_dict, orient="index")
            assembly_df.columns = ["assembly_file"]
        else:
            assembly_df = pandas.DataFrame(
                {
                    "assembly_file": [],
                }
            )
        if paired_end_read_data_folder is not None:
            paired_end_read_data_dict, unmatched_read_files = (
                cls.parse_folder_for_paired_end_reads(
                    data_folder=paired_end_read_data_folder,
                    paired_end_reads_pattern=paired_end_reads_pattern,
                )
            )
        if paired_end_read_data_folder is not None and len(paired_end_read_data_dict) > 0:
            paired_end_read_df = pandas.DataFrame(paired_end_read_data_dict).transpose()
        else:
            paired_end_read_df = pandas.DataFrame(
                {"sample_name": [], "r1_file": [], "r2_file": []}
            )

        combined_df = assembly_df.merge(
            paired_end_read_df, how="outer", left_index=True, right_on="sample_name"
        )
        combined_df.replace(numpy.nan, None, inplace=True)
        combined_df = combined_df[
            ["sample_name", "assembly_file", "r1_file", "r2_file"]
        ]
        combined_dict = combined_df.to_dict("index")

        return cls(file_paths=combined_dict)

    def __iter__(self):
        if self.sample_data is None:
            for sample_name, files in self.file_paths.items():
                yield sample_name, files, None
        else:
            for sample_name, files in self.file_paths.items():
                yield sample_name, files, self.sample_data[sample_name]
    
    def __len__(self):
        return len(self.sample_names)

    def get_missing_files(self):
        missing_assembly_files = {}
        missing_read_files = {}
        for sample_name, assembly_file, r1_file, r2_file in self.filepaths():
            if assembly_file is None or not assembly_file.exists():
                missing_assembly_files[sample_name] = assembly_file
            if (
                r1_file is None
                or r2_file is None
                or not r1_file.exists()
                or not r2_file.exists()
            ):
                missing_read_files[sample_name] = (r1_file, r2_file)
        return missing_assembly_files, missing_read_files

    def wgs_files_exist(self):
        missing_assembly_files, missing_read_files = self.get_missing_files()
        if len(missing_assembly_files) == 0 and len(missing_read_files) == 0:
            return True
        else:
            return False

    @classmethod
    def from_dataframe(cls, pandas_dataframe: pandas.DataFrame):
        sample_data = pandas_dataframe.to_dict("index")
        sample_data.replace(numpy.nan, None, inplace=True)
        return cls(file_paths=sample_data)

    def filepaths(self):
        for i, sample_name in enumerate(self.sample_names):
            yield sample_name, self.assembly_files[i], self.paired_end_read_files[i][
                0
            ], self.paired_end_read_files[i][1]

    def readpaths(self):
        for i, sample_name in enumerate(self.sample_names):
            yield sample_name, self.paired_end_read_files[0], self.paired_end_read_files[1]

    def assemblypaths(self):
        for i, sample_name in enumerate(self.sample_names):
            yield sample_name, self.assembly_files[i]

    def __repr__(self):
        if self.sample_data:
            return f"< WgsData object with {len(self.sample_names)} samples and {len(self.sample_data_columns)} metadata columns >"
        else:
            return f"< WgsData object with data from {len(self.sample_names)} samples >"

    @staticmethod
    def parse_folder_for_paired_end_reads(
        data_folder: Path,
        paired_end_reads_pattern: str = r"(?P<sample_name>.+?)(?P<sample_number>(_S[0-9]+)?)(?P<lane>(_L[0-9]+)?)[\._]R?(?P<paired_read_number>[1|2])(?P<set_number>(_[0-9]+)?)(?P<file_extension>\.fastq\.gz)",
    ) -> tuple[dict[dict], set]:
        """
        Parses a provided folder for paired end reads matching the provided regex pattern
        Returns a dictionary with sample names as keys and a dict with keys "r1_file" and "r2_file" and filenames as values
        In case of duplicate sample names, the <sample_name>_<sample_number> combination will be used instead
        If <sample_name>_<sample_number> is also not unique, the key used for sample name will simply be the name of the r1_file
        """
        sample_list = []
        unmatched_files = set()
        sample_names = set()
        duplicate_sample_names = set()
        duplicate_sample_names_uniq = set()
        data_folder = Path(data_folder).resolve()
        for f in data_folder.iterdir():
            fname = f.name
            fname_match = re.match(paired_end_reads_pattern, fname)
            if fname_match:
                sample_name = fname_match.group("sample_name")
                if fname_match.group("paired_read_number") is not None:
                    read_number = fname_match.group("paired_read_number")
                    if fname_match.group("sample_number") is not None:
                        sample_name_uniq = (
                            f"{sample_name}{fname_match.group('sample_number')}"
                        )
                    else:
                        sample_name_uniq = f"{sample_name}"
                    if read_number == "1":
                        r1_file = f.absolute()
                        r2_file = data_folder.joinpath(
                            fname[: fname_match.start("paired_read_number")]
                            + "2"
                            + fname[fname_match.end("paired_read_number") :]
                        )
                        if r2_file.exists():
                            sample_list.append(
                                (sample_name, sample_name_uniq, r1_file, r2_file)
                            )
                            if sample_name in sample_names:
                                duplicate_sample_names.add(sample_name)
                                if sample_name_uniq in duplicate_sample_names:
                                    duplicate_sample_names_uniq.add(sample_name_uniq)
                            else:
                                sample_names.add(sample_name)
                        else:
                            unmatched_files.add(fname)
                    if read_number == "2":
                        r1_file = data_folder.joinpath(
                            fname[: fname_match.start("paired_read_number")]
                            + "1"
                            + fname[fname_match.end("paired_read_number") :]
                        )
                        if not r1_file.exists():
                            unmatched_files.add(fname)
        sample_file_dict = {}
        for sample_name, sample_name_uniq, r1_file, r2_file in sample_list:
            if sample_name in duplicate_sample_names:
                if sample_name in duplicate_sample_names_uniq:
                    sample_file_dict[r1_file] = {"r1_file": r1_file, "r2_file": r2_file}
                else:
                    sample_file_dict[sample_name_uniq] = {
                        "sample_name": sample_name,
                        "r1_file": r1_file,
                        "r2_file": r2_file,
                    }
            else:
                sample_file_dict[sample_name] = {
                    "sample_name": sample_name,
                    "r1_file": r1_file,
                    "r2_file": r2_file,
                }

        return sample_file_dict, unmatched_files

    @staticmethod
    def parse_folder_for_fasta_files(
        data_folder: Path,
        fasta_file_pattern: str = r"(?P<sample_name>.+?)(\.fa|\.fna|\.fasta)$",
    ) -> tuple[dict[dict], set]:
        """
        Parses a provided folder for paired end reads matching the provided regex pattern
        Returns a dictionary with sample names as keys and a dict with keys "r1_file" and "r2_file" and filenames as values
        In case of duplicate sample names, the <sample_name>_<sample_number> combination will be used instead
        If <sample_name>_<sample_number> is also not unique, the key used for sample name will simply be the name of the r1_file
        """
        sample_file_dict = {}
        data_folder = Path(data_folder).resolve()
        for f in data_folder.iterdir():
            fname = f.name
            fname_match = re.match(fasta_file_pattern, fname)
            if fname_match:
                sample_name = fname_match.group("sample_name")
                sample_file_dict[sample_name] = f
        return sample_file_dict