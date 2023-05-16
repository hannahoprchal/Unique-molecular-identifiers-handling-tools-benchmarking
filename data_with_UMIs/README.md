# Datasets containing unique molecular identifiers
Datasets can be downloaded all at once with the command: `make datasets_all`

Optionally, you can download a specific dataset by defining the dataset name, an example  command: `make dataset1`
## Data origin
* dataset1, dataset2, dataset3 comes from [UMIc tool by BiodataAnalysisGroup](https://github.com/BiodataAnalysisGroup/UMIc/tree/master/data) repository
* dataset4 comes from [UMI-tools by CGATOxford](https://github.com/CGATOxford/UMI-tools/releases/tag/v0.2.3) repository

## UMIs location
* dataset1
    * read1: UMIs location = first 12 bases
    * read2: UMIs location = None

* dataset2
    * read1: UMIs location = first 10 bases
    * read2: UMIs location = first 10 bases

* dataset3
    * read1: UMIs location = first 12 bases

* dataset4
    * read1: UMIs location = first 9 bases