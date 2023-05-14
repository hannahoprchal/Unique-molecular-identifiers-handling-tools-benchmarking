# [UMI-tools workflow](https://github.com/CGATOxford/UMI-tools)
## Input parameters
* String sample_name - prefix added to every name of output files
* File input_fastq_r1 - input file containing read1 in FASTQ format compressed to .gz
* File? input_fastq_r2 - (optional) input file containing read2 in FASTQ format compressed to .gz
* String bc_pattern - UMIs pattern in read1 file, format according to [documentation](https://umi-tools.readthedocs.io/en/latest/reference/extract.html#barcode-extraction)
* String? bc_pattern2 - (optional) UMIs pattern in read2 file, format according to [documentation](https://umi-tools.readthedocs.io/en/latest/reference/extract.html#barcode-extraction)
* String output_read1 - name of the output FASTQ when UMIs extracted
* String? output_read2 - (optional) name of the output FASTQ when UMIs extracted
* File ref_file - reference genome file in FASTA format compressed to .gz
* File ref_index_bwa - bwa index for the reference genome compressed to .tar.gz
* String docker_image - file id of the docker image on DNAnexus platform
* The rest of the input parameters are only informative, can be changed accordingly


