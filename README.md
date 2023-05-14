# Unique-molecular-identifiers-handling-tools-benchmarking
Implementation for the bachelor thesis Unique molecular identifiers handling tools benchmarking

## WDL task and workflow compilation
* To compile WDL task or workflow to the DNAnexus platform you have to download [dxCompiler](https://github.com/dnanexus/dxCompiler#setup)
* Make sure you have also installed [dx-toolkit](https://documentation.dnanexus.com/downloads)
* When you have WDL task(or workflow) ready to compile run following command: `java -jar /path/to/dxcompiler_jar/ compile task.wdl `

## Task and workflow run on DNAnexus platform
When WDL task or workflow is successfully compiled into platform into runnable applet (application), you can use following options to run your applet (application):
* Run the task (workflow) via command: `dx run applet_id -iinput_variable -ianother_input_variable`

Or
* Run the task (workflow) on the platform [GUI](https://platform.dnanexus.com/) and fill manually input variables

