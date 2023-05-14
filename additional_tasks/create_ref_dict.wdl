version 1.0

task create_ref_dict {
    input {
        File ref_file
	String dict_name
        String docker_image
    
    }

    command <<<
        java -jar $PICARD CreateSequenceDictionary \ 
            -R ~{ref_file} \ 
            -O ~{dict_name}.dict

    >>>

    runtime {
        docker: docker_image
    }


    output {
        File ref_dict = "~{dict_name}.dict"
    }

}
            
