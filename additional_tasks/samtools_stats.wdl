version 1.0

task samtools_stats {
    input {
        File input_bam
        String sample_name
        String docker_image = "dx://file-GQXY47Q0Pf0v43v88yY1FX3Y"

    }
    command <<<
        samtools stats ~{input_bam} | grep ^SN | cut -f 2- > ~{sample_name}_samtools_stats.txt
    >>>

    runtime{
        docker: docker_image
    }

    output {
        File output_stats = "~{sample_name}_samtools_stats.txt"
    }
}