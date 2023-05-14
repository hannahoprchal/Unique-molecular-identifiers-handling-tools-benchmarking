version 1.0



workflow umitools_WF {
  input {
        String sample_name = "umi_tools"
        File input_fastq_r1
        File? input_fastq_r2
        String bc_pattern       
        String? bc_pattern2
        String output_read1 = "extracted_R1.fastq.gz"
        String? output_read2
        
        String seq_platform = "ILLUMINA"
        Array[String]id = ["0","umi_tools"]
        Array[Int] pi = [0]
        String stats_prefix = "deduplicated_stats"

        File ref_file
        File ref_index_bwa
        String docker_image = "dx://file-GPqXBQQ0Pf0Xb3fz27f5JFPz"
  }

    call umi_tools_extract {
    input: 
        input_fastq_r1 = input_fastq_r1,
        input_fastq_r2 = input_fastq_r2,
        bc_pattern = bc_pattern,
        bc_pattern2 = bc_pattern2,
        output_read1 = output_read1,
        output_read2 = output_read2,
        sample_name = sample_name,
        docker_image = docker_image
  }

    call umi_tools_bwa_mem {
    input: 
        extracted_fastq1 = umi_tools_extract.extracted_read1,
        extracted_fastq2 = umi_tools_extract.extracted_read2,
        ref_file = ref_file,
        ref_index = ref_index_bwa,
        sample_name = sample_name,
        seq_platform = seq_platform,
        id = id,
        pi = pi,
        docker_image = docker_image
  }


    call umi_tools_dedup {
        input: 
            mapped_input_bam = umi_tools_bwa_mem.out_bam,
            mapped_input_bam_index = umi_tools_bwa_mem.out_bai,
            sample_name = sample_name,
            stats_prefix = stats_prefix,
            is_paired = defined(input_fastq_r2),
            docker_image = docker_image
    }

    call umi_tools_group {
        input: 
            mapped_input_bam = umi_tools_bwa_mem.out_bam,
            mapped_input_bam_index = umi_tools_bwa_mem.out_bai,
            sample_name = sample_name,
            is_paired = defined(input_fastq_r2),
            docker_image = docker_image
    }


  output {
        File dedupped_bam = umi_tools_dedup.dedupped_bam
        File edit_distance = umi_tools_dedup.edit_distance
        File umi_stats = umi_tools_dedup.umi_stats
        File position_stats =  umi_tools_dedup.position_stats

        File groupped_bam = umi_tools_group.groupped_bam
        File groupped_info_bam = umi_tools_group.groupped_info_bam
  }
}


task umi_tools_extract {
    input {
        File input_fastq_r1    
        File? input_fastq_r2

        String bc_pattern
        String? bc_pattern2
        String output_read1
        String? output_read2
        String sample_name
        String docker_image 
    }

    command {
        umi_tools extract \
        --stdin ~{input_fastq_r1} \
        ~{"--read2-in " + input_fastq_r2} \
        --bc-pattern ~{bc_pattern} \
        ~{"bc-pattern2 " + bc_pattern2} \
        ~{"--stdout " + sample_name + "_" + output_read1} \
        ~{if defined(input_fastq_r2) then "--read2-out " + "~{sample_name}_" + output_read2 else ""}
    }

    runtime {
        docker: docker_image
    }

    meta {
         title: "Umi-tools extract"
        description: "Extracts UMIs from reads in .fq.gz files."
    }

    output {
        File extracted_read1 = "~{sample_name}_~{output_read1}"
        File? extracted_read2 = "~{sample_name}_~{output_read2}"
    }


}


task umi_tools_bwa_mem {
    input {
        File extracted_fastq1
        File? extracted_fastq2
        File ref_file
        File ref_index
        String sample_name
        String seq_platform
        Array[String]id
        Array[Int] pi
        String docker_image

    }

    command <<<
        set -exo pipefail
        mkdir genome

        tar zxvf ~{ref_index} -C genome
        ref_file=$(ls genome/*.bwt)
        ref_file="${ref_file%.bwt}"

        bwa mem \
            -M \
            -t $(nproc) \
            -R $(echo "@RG\tID:~{id[0]}\tSM:~{sample_name}\tPL:~{seq_platform}\tPI:~{id[0]}") \
            ${ref_file} \
            ~{extracted_fastq1}  ~{if defined(extracted_fastq2) then "~{extracted_fastq2} " else ""}\ 
            -K 100000000 \
            -Y \
        | samtools sort \
            -O bam \
            -@ $(nproc) \
            -T "genome/" \
            -o ~{sample_name}_alligned.bam -

            samtools index ~{sample_name}_alligned.bam -o ~{sample_name}_alligned.bai


    >>>

    runtime {
        docker: docker_image
        dx_instance: "mem1_ssd1_v2_x8"
    }

    meta {
         title: "BWA mem, Samtools sort, Samtools index"
        description: "Maps on reference genome, sorts the output and creates index for mapped file."
    }

    output {
        File out_bam = "~{sample_name}_alligned.bam"
        File out_bai = "~{sample_name}_alligned.bai"

    }
}

task umi_tools_dedup {
    input {
        File mapped_input_bam
        File mapped_input_bam_index
        String sample_name
        Boolean is_paired
        String stats_prefix
        String docker_image
    }


    command <<<
        umi_tools dedup \
            -I ~{mapped_input_bam} \
            ~{if is_paired then "--paired " else ""} \ 
            --output-stats "~{sample_name}_~{stats_prefix}" \
            -S ~{sample_name}_deduplicated.bam

    >>>

    runtime {
        docker: docker_image
    }

    meta {
         title: "Umi-tools dedup"
        description: "Dedups (collapses) by UMI, outputs stats about data."
    }

    output {
        File dedupped_bam = "~{sample_name}_deduplicated.bam"
        
        File edit_distance = "~{sample_name}_~{stats_prefix}_edit_distance.tsv"
        File umi_stats = "~{sample_name}_~{stats_prefix}_per_umi.tsv"
        File position_stats =  "~{sample_name}_~{stats_prefix}_per_umi_per_position.tsv"
    }
}

task umi_tools_group {
    input {
        File mapped_input_bam
        File mapped_input_bam_index
        String sample_name
        Boolean is_paired

        String docker_image
    }


    command <<<
        umi_tools group \
         -I ~{mapped_input_bam} \
         ~{if is_paired then "--paired " else ""} \ 
         --group-out ~{sample_name}_groups.tsv \
         --output-bam \
         -S ~{sample_name}_mapped_grouped.bam

    >>>

    runtime {
        docker: docker_image
    }

    meta {
         title: "Umi-tools group"
        description: "Groups by UMI, outputs info about groups."
    }

    output {
        File groupped_bam = "~{sample_name}_mapped_grouped.bam"
        File groupped_info_bam = "~{sample_name}_groups.tsv"
    }
}
