version 1.0



workflow umi_collapse_WF {
  input {
        String sample_name = "UMIcollapse"
        Int umi_len
        Int umi_skip
        String umi_loc      #location of umi - read1, read2 or per_read
        File fastq_read_1
        File? fastq_read_2
        String read_structure

        File bwa_index
        File ref_file

        String docker_image = "dx://file-GPQZF980Pf0b6YxGxq9ZKBQ7"
  }
    call umicollapse_fastp_umi_extract{
        input: sample_name = sample_name,
            umi_len = umi_len,
            umi_skip = umi_skip,
            umi_loc = umi_loc,
            fastq_read_1 = fastq_read_1, 
            fastq_read_2 = fastq_read_2,
            read_structure = read_structure,
            docker_image = docker_image
    }
  

    call umicollapse_bwa_mem {
        input: trimmed_fastq_R1 = umicollapse_fastp_umi_extract.trimmed_fastq1,
            trimmed_fastq_R2 = umicollapse_fastp_umi_extract.trimmed_fastq2,
            bwa_index = bwa_index,
            ref_file = ref_file,
            sample_name = sample_name,
            docker_image = docker_image
  }


    call umi_collapse_final {
        input: precollapsed_bam = umicollapse_bwa_mem.out_bam,
                sample_name = sample_name,
                docker_image = docker_image,
                precollapsed_bai = umicollapse_bwa_mem.out_bai,
                is_paired = defined(fastq_read_2)
    }



  output {
    File final_output_bam = umi_collapse_final.out_bam
    File final_output_bai = umi_collapse_final.out_bai
  }
}


task umicollapse_fastp_umi_extract {
    input {
        String sample_name
        Int umi_len
        Int umi_skip
        String umi_loc      #location of umi - read1, read2 or per_read
        File fastq_read_1
        File? fastq_read_2
        String read_structure

        String docker_image
    
    }

    command {
        fastp \
        -i ~{fastq_read_1} \
        ~{if defined(fastq_read_2) then "-I ~{fastq_read_2} " else ""} \
        -o ~{sample_name}.trimmed.R1.fastq.gz \
        ~{if defined(fastq_read_2) then "-O ~{sample_name}.trimmed.R2.fastq.gz " else ""} \ 
        --json Metrics/~{sample_name}.metrics.fastp.json \
        --html Metrics/~{sample_name}.metrics.fastp.html \
        --report_title "~{sample_name} trimming report" \
        --umi \
        --umi_loc ~{umi_loc} \
        --umi_len ~{umi_len} \
        --umi_skip ~{umi_skip}

        

    }

    runtime {
        docker: docker_image
    }

    meta {
         title: "Fastp preprocessing"
        description: "Preprocesses, trimms UMIs from reads in.fq.gz files."
    }


    output {
        File trimmed_fastq1 = "~{sample_name}.trimmed.R1.fastq.gz"
        File? trimmed_fastq2 = "~{sample_name}.trimmed.R2.fastq.gz"


    }
}



task umicollapse_bwa_mem{
    input {
        File trimmed_fastq_R1
        File? trimmed_fastq_R2
        File bwa_index
        File ref_file
        String sample_name
        
        String docker_image

    }

    command <<<
        set -exo pipefail
        mkdir genome

        tar zxvf ~{bwa_index} -C genome
        ref_file=$(ls genome/*.bwt)
        ref_file="${ref_file%.bwt}"

        bwa mem \
        -M \
        -R $(echo "@RG\tID:~{sample_name}\tSM:~{sample_name}\tPL:ILLUMINA\tLB:~{sample_name}") \
        -t 4 \
        ${ref_file} \
        ~{trimmed_fastq_R1} ~{if defined(trimmed_fastq_R2) then "~{trimmed_fastq_R2} " else ""} \ 
        | \
        samtools sort \
        -O bam \
        --threads 3 \
        -T "genome/" \
        -o ~{sample_name}.precollapsed.bam -
        
        samtools index ~{sample_name}.precollapsed.bam -o ~{sample_name}.precollapsed.bai


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
        File out_bam = "~{sample_name}.precollapsed.bam"
        File out_bai = "~{sample_name}.precollapsed.bai"

    }
}


task umi_collapse_final {
    input {
        File precollapsed_bam
        File precollapsed_bai
        String sample_name
        String docker_image
        Boolean is_paired

    }

    command <<<
            java -jar /UMICollapse/umicollapse.jar bam \
            --two-pass \
            ~{if is_paired then "--paired " else ""} \ 
            --umi-sep : \
            -i ~{precollapsed_bam} \
            -o ~{sample_name}.final.bam &&
            samtools index ~{sample_name}.final.bam -o ~{sample_name}_final_index.bai
    >>>

    runtime {
        docker: docker_image
    }

    meta {
         title: "UmiCollapse bam"
        description: "Deduplicates (collapses) UMIs"
    }

    output {
        File out_bam = "~{sample_name}.final.bam"
        File out_bai = "~{sample_name}_final_index.bai"

    }
}
