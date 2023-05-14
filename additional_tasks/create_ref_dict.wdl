version 1.0

task fgbio_create_ref_dict {
    input {
        File ref_file
        String docker_image = "dx://file-GPjVKVQ0Pf0bGpVF5bQfVGP3"
    
    }

    command <<<
        java -jar $PICARD CreateSequenceDictionary \ 
            -R ~{ref_file} \ 
            -O reference.dict

    >>>

    runtime {
        docker: docker_image
    }


    output {
        File ref_dict = "reference.dict"
    }

}
            