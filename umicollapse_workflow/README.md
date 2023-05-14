# [UMICollapse workflow](https://github.com/Daniel-Liu-c0deb0t/UMICollapse)
## Input parameters
* String sample_name - prefix added to every name of output files
* File fastq_read_1 - input file containing read1 in FASTQ format compressed to .gz
* File? fastq_read_2 - (optional) input file containing read2 in FASTQ format compressed to .gz
* Int umi_len - UMIs length, format according to [description](https://github.com/OpenGene/fastp#all-options)
* Int umi_skip - bases skip, format according to [description](https://github.com/OpenGene/fastp#all-options)
* String umi_loc - UMIs location, format according to [description](https://github.com/OpenGene/fastp#all-options)
* String read_structure - UMIs pattern for both read files, format according to [description](https://github.com/OpenGene/fastp#all-options)
* File bwa_index - bwa index for the reference genome compressed to .tar.gz
* File ref_file - reference genome file in FASTA format compressed to .gz
* String docker_image - file id of the docker image on DNAnexus platform