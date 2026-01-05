### SepiPOINT ###

SepiPOINT is a simple tool for identifying mutations associated with antimicrobial resistance in whole genome sequencing data from Staphylococcus epidermidis isolates.




### Dependencies ###

 - Python >= 3.9
 - pandas
 - numpy
 - bwa (Only for paired-end read input)
 - samtools >= 1.22.1 (Only for paired-end read input)
 - bcftools >= 1.22 (Only for paired-end read input)
 - mummer (Only genome assembly input)

### Installation ###

#### From pypi ####

pip install sepi_point

#### From Conda ####

conda install thej-ssi::sepi_point
