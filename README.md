# SepiPOINT #

SepiPOINT is a simple tool for identifying mutations associated with antimicrobial resistance in whole genome sequencing data from *Staphylococcus epidermidis* isolates.



## Installation ##

### From pypi ###

```
pip install sepi_point
```

### From Conda ###

```
conda create -n sepi_point

conda activate sepi_point

conda install thej-ssi::sepi_point

sepi_point -h
```


### Dependencies ###

 - Python >= 3.9
 - pandas
 - numpy
 - bwa (Only for paired-end read input)
 - samtools >= 1.22.1 (Only for paired-end read input)
 - bcftools >= 1.22 (Only for paired-end read input)
 - mummer (Only genome assembly input)


Dependencies will be installed automatically if installed with conda install, but when installing via pypi, non-python dependencies will have to be added manually (bwa, samtools, bcftools, mummer)


## Usage ##

To run on paired-end read input from a single isolate:

`sepi_point -1 <R1_file.fastq.gz> -2 <R2_file.fastq.gz> -o <output_folder>`



To run on a single genome assembly:

`sepi_point -a <assembly_file.fasta> -o <output_folder>`



To run on all genome assemblies in a folder

`sepi_point_batch -a <path_to_assembly_files_folder> -o <output_folder>`



To run on all paired-end reads in a folder

`sepi_point_batch -r <path_to_read_files_folder> -o <output_folder>`



To run on both assemblies and reads

`sepi_point_batch -a <path_to_assembly_files_folder> -r <path_to_read_files_folder> -o <output_folder>`


## Inputs ##

SepiPoint expects inputs as .fastq.gz files for paired end reads and fasta-format for assembled genomes.

In batch mode the specified folder will be parsed for fastq and or fasta-files. Fasta files are identified by standard prefixes (.fasta, .fa, .fna) and fastq files must follow standard Illumna naming or simple SRA-like naming convention (*_R1.fastq.gz, *_1.fastq.gz, *.R1.fastq.gz, *.1.fasta.gz)


## Outputs ##

In single isolate mode, all resistance-associated mutations identified will be presented in a tsv-file (*.results.tsv) like this:

| Mutation |	Gene |	Position |	Ref |	Alt |	Ref_codon |	Alt_codon |	Alt_frequency |	Category |
| -------- |	---- |	-------: |	--- |	--- |	--------- |	--------- |	------------- |	-------- |
| thyA-intergenic::A286G |	thyA-intergenic |	286 |	A |	G |	 |	 |	25/35 |	Trimethoprim-sulfamethoxazole |
| gyrA::S84F |	gyrA |	84 |	S |	F |	TCT |	TTT |	144/144 |	Fluoroquinolone |
| parC::S80F |	parC |	80 |	S |	F |	TCT |	TTT |	98/99 |	Fluoroquinolone |
| qacA4::A157G |	qacA4 |	157 |	A |	G |	GCT |	GGT |	190/190 |	Chlorhexidine |
| rpoB::D471E |	rpoB |	471 |	D |	E |	GAC |	GAA |	164/164 |	Rifampicin, vancomycin |
| rpoB::I527M |	rpoB |	527 |	I |	M |	ATA |	ATG |	163/163 |	Rifampicin, vancomycin |



In batch mode each isolate in the provided input folder(s) will have their own subfolder within the output folder containing results from that isolate.

In addition, the base output folder will contain *results.tsv* with the combined results from all isolates, as well as *results.matrix.tsv*, a 0/1 filled matrix with the presence/absence of each mutation in each isolate.











