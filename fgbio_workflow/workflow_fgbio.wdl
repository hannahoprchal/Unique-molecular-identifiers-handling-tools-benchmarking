
version 1.0

workflow fgbio_WF {
    input {
        String sample_name = "fgbio"
        File fastq_r1
        File? fastq_r2
        String read_structure
        String library = "ILLUMINA"
        File ref_bwa_index
        File ref_fa
        File ref_dict
        String docker_image = "dx://file-GPFjqk80Pf0zjzQB73xvB0qy"
    }

    call fastq_to_ubam {
        input:
            fastq_r1 = fastq_r1
            fastq_r2 = fastq_r2
            read_structure = read_structure
            sample_name = sample_name
            library = library
            docker_image = docker_image
    }

    call ubam_to_bam {
        input:
            ubam = fastq_to_ubam.unmapped_bam
            ref_index = ref_bwa_index
            ref_fa = ref_fa
            ref_dict = ref_dict
            docker_image = docker_image
    }

    call bam_to_grouped {
        input:
            mapped_bam = ubam_to_bam.mapped_bam
            docker_image = docker_image
    }

    call grouped_to_consensus {
        input:
           grouped_bam = bam_to_grouped.grouped_bam
           ref_index = ref_bwa_index
           ref_fa = ref_fa
           ref_dict = ref_dict
           docker_image = docker_image
    }

    output {
        File cons_bam = grouped_to_consensus.cons_bam
        File cons_mapped_bam = grouped_to_consensus.cons_mapped_bam
        File cons_filtered_bam = grouped_to_consensus.cons_filtered_bam
    }  
}


task fastq_to_ubam {
    input{
        File fastq_r1
        File? fastq_r2
        String read_structure
        String sample_name 
        String library
        String docker_image
    }
    command <<<
        set -exo pipefail
        
        java -jar $FGBIO FastqToBam \
            --input ~{fastq_r1} ~{if defined(fastq_r2) then "~{fastq_r2} " else ""} \ 
            --read-structures ~{read_structure} \
            --sample ~{sample_name} \
            --platform-unit ABCDEAAXX.1 \
            --library ~{library} \
            --output  ~{sample_name}_unmapped.bam
    >>>
    runtime {
        docker: docker_image
    }

    meta {
         title: "Fgbio Fastq to ubam"
        description: "Extracts the UMI bases and stores them in standard SAM tags."
    }

    output {
        File unmapped_bam = "~{sample_name}_unmapped.bam"
    }
}

task ubam_to_bam {
    input{
        File ubam
        File ref_index
        File ref_fa
        File ref_dict

        String docker_image
    }

    String sample_name = basename(ubam,"_unmapped.bam")

    command <<<
        set -exo pipefail

        mkdir genome

        tar zxvf ~{ref_index} -C genome
        ref_file=$(ls genome/*.bwt)
        ref_file="${ref_file%.bwt}"

        mv ~{ref_fa} .
		mv ~{ref_dict} .
		reference_genome=$(basename ~{ref_fa})

        samtools fastq ~{ubam} \
            | bwa mem -t $(nproc) \
                -p \
                -K 100000000 \
                -Y ${ref_file} \
                - \
            | java -jar $FGBIO --compression 1 --async-io ZipperBams \
                --unmapped ~{ubam} \
                --ref ${reference_genome} \
                --output ~{sample_name}_mapped.bam
    >>>
    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x8"
    }

    meta {
         title: "Fgbio ubam to mapped bam"
        description: "Converts unmapped bam to fastq, maps to ref genome."
    }


    output {
        File mapped_bam = "~{sample_name}_mapped.bam"
    }
}

task bam_to_grouped {
    input{
        File mapped_bam

        String docker_image
    }

    String sample_name = basename(mapped_bam,"_mapped.bam")

    command <<<
        set -exo pipefail

        java -jar $FGBIO --compression 1 --async-io GroupReadsByUmi \
            --input ~{mapped_bam} \
            --strategy Adjacency \
            --edits 1 \
            --output ~{sample_name}_grouped.bam \
            --family-size-histogram  ~{sample_name}.tag-family-sizes.txt
        
    >>>
    runtime {
        docker: docker_image
    }

    meta {
         title: "Fgbio mapped to grouped"
        description: "Identifies reads that originate from the same source molecule based on genomic positions and UMI."
    }


    output {
        File grouped_bam = "~{sample_name}_grouped.bam"
        File histogram = "~{sample_name}.tag-family-sizes.txt"
    }
}

task grouped_to_consensus {
    input{
        File grouped_bam
        File ref_index
        File ref_fa
        File ref_dict
        String docker_image
    }

    String sample_name = basename(grouped_bam,"_grouped.bam")

    command <<<
        set -exo pipefail

        mkdir genome

        tar zxvf ~{ref_index} -C genome
        ref_file=$(ls genome/*.bwt)
        ref_file="${ref_file%.bwt}"

        mv ~{ref_fa} .
		mv ~{ref_dict} .
		reference_genome=$(basename ~{ref_fa})

        java -jar $FGBIO --compression 1 --async-io CallMolecularConsensusReads \
            --input ~{grouped_bam} \
            --output ~{sample_name}_cons.unmapped.bam \
            --min-reads 1 \
            --min-input-base-quality 5 \
            --threads $(nproc)

        samtools fastq ~{sample_name}_cons.unmapped.bam \
        | bwa mem -t $(nproc) -p -K 100000000 -Y ${ref_file} - \
        |  java -jar $FGBIO --compression 1 --async-io ZipperBams \
            --unmapped ~{sample_name}_cons.unmapped.bam \
            --ref ${reference_genome} \
            --tags-to-reverse Consensus \
            --tags-to-revcomp Consensus \
            #| samtools sort --threads $(nproc) -o ~{sample_name}_cons.mapped.bam --write-index
            --output ~{sample_name}_cons.mapped.bam

        java -jar $FGBIO --compression 0 FilterConsensusReads \
            --input ~{sample_name}_cons.mapped.bam \
            --output /dev/stdout \
            --ref ${reference_genome} \
            --min-reads 3 \
            --min-base-quality 45 \
            --max-base-error-rate 0.2 \
            | samtools sort --threads $(nproc) -o ~{sample_name}_cons.filtered.bam --write-index
        
    >>>
    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x8"
    }

    meta {
         title: "Fgbio grouped to consensus"
        description: "Generates unmapped consensus reads from the grouped reads, remaps the consensus reads, filters them and sorts them"
    }


    output {
        File cons_bam = "~{sample_name}_cons.unmapped.bam"
        File cons_mapped_bam = "~{sample_name}_cons.mapped.bam"
        File cons_filtered_bam = "~{sample_name}_cons.filtered.bam"
    }
}

