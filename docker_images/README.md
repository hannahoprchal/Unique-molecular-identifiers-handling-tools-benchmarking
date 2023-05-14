# Docker images
fgbio - image to run [fgbio workflow](https://github.com/hannahoprchal/Unique-molecular-identifiers-handling-tools-benchmarking/blob/main/fgbio_workflow/workflow_fgbio.wdl)
picard - image to run [create_ref_dict task](https://github.com/hannahoprchal/Unique-molecular-identifiers-handling-tools-benchmarking/blob/main/additional_tasks/create_ref_dict.wdl)
samtools - image to run [samtools_stats task](https://github.com/hannahoprchal/Unique-molecular-identifiers-handling-tools-benchmarking/blob/main/additional_tasks/samtools_stats.wdl)
umicollapse - image to run [UMICollapse workflow](https://github.com/hannahoprchal/Unique-molecular-identifiers-handling-tools-benchmarking/blob/main/umicollapse_workflow/workflow_umicollapse.wdl)
umi_tools - image to run [UMI-tools workflow](https://github.com/hannahoprchal/Unique-molecular-identifiers-handling-tools-benchmarking/blob/main/umi_tools_workflow/workflow_umi_tools.wdl)



# Build and upload to DNAnexus platform
## Build
* Build docker image locally with command: `docker build . -t tool_name`
* Save docker image to tar.gz with command: `docker save tool_name | gzip -c > tool_docker_image.tar.gz`
## Upload
* Upload tar.gz to DNAnexus platform with command: `dx upload tool_docker_image.tar.gz`

Or
* Upload tar.gz to DNAnexus platform manually via GUI

