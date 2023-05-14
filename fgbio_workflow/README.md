# [fgbio workflow](https://github.com/fulcrumgenomics/fgbio)
## Input parameters
* String sample_name - prefix added to every name of output files
* File fastq_r1 - input file containing read1 in FASTQ format compressed to .gz
* File? fastq_r2 - (optional) input file containing read2 in FASTQ format compressed to .gz
* String read_structure - UMIs pattern for both read files, format according to [documentation](http://fulcrumgenomics.github.io/fgbio/tools/latest/FastqToBam.html)
* String library - info parameter, can be changed accordingly
* File ref_bwa_index - bwa index for the reference genome compressed to .tar.gz
* File ref_fa - reference genome file in FASTA format compressed to .gz
* File ref_dict - reference dictionary file in .dict format, can be obtained with [create_ref_dict.wdl](https://github.com/hannahoprchal/Unique-molecular-identifiers-handling-tools-benchmarking/blob/main/additional_tasks/create_ref_dict.wdl)
* String docker_image - file id of the docker image on DNAnexus platform
