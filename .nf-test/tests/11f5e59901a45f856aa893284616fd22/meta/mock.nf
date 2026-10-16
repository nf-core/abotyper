// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies


// include test workflow
include { MINIMAP2_ALIGN_READS } from '/mnt/storage_c2/bioinformatics/abotyper/subworkflows/local/minimap2_align_reads/tests/../main.nf'

workflow {

    // define custom rules for JSON that will be generated.
    def jsonOutput = createJsonOutput()
    def jsonWorkflowOutput = createJsonWorkflowOutput()

    def input = []

    // run dependencies
    

    // workflow mapping
    input = []
    
                input[0] = channel.of([ [ id:'ABO_AA_barcode01', single_end:true ],
                                        [ file("/mnt/storage_c2/bioinformatics/abotyper/assets/samplesheet.csv", checkIfExists: true) ] ])
                input[1] = channel.value([ [ id:'NG_006669v2' ], file("/mnt/storage_c2/bioinformatics/abotyper/assets/refs/NG_006669v2.fasta", checkIfExists: true) ])
                input[2] = channel.value([ [ id:'NG_006669v2' ], file("/mnt/storage_c2/bioinformatics/abotyper/assets/refs/NG_006669v2.fasta.fai", checkIfExists: true) ])
                
    //----

    //run workflow
    MINIMAP2_ALIGN_READS.run(input.toArray())
    
    if (MINIMAP2_ALIGN_READS.output){

        // consumes all named output channels and stores items in a json file
        MINIMAP2_ALIGN_READS.out.getNames().each { name ->
            serializeChannel(name, MINIMAP2_ALIGN_READS.out.getProperty(name), jsonOutput, params.nf_test_output)
        }	  
    
        // consumes all unnamed output channels and stores items in a json file
        def array = MINIMAP2_ALIGN_READS.out as List<Object>
        def i = 0
        array.each { output ->
            serializeChannel(i, output, jsonOutput, params.nf_test_output)
            i += 1
        }    	

    }

    workflow.onComplete = {

        def result = [
            success: workflow.success,
            exitStatus: workflow.exitStatus,
            errorMessage: workflow.errorMessage,
            errorReport: workflow.errorReport
        ]
        new File("${params.nf_test_output}/workflow.json").text = jsonWorkflowOutput.toJson(result)
        
    }
}


def serializeChannel(name, channel, jsonOutput, outputDir) {
    def _name = name
    def list = [ ]
    channel.subscribe(
        onNext: { entry ->
            list.add(entry)
        },
        onComplete: {
            def map = new HashMap()
            map[_name] = list
            def filename = "${outputDir}/output_${_name}.json"
            new File(filename).text = jsonOutput.toJson(map)		  		
        } 
    )
}

def createJsonOutput(_input = null) {
    // _input is needed because a closure is provided to all functions called in the process
    return [
        toJson: { obj ->
            def converted = convertPathsToStrings(obj)
            return groovy.json.JsonOutput.toJson(converted)
        }
    ]
}

def convertPathsToStrings(obj) {
    if (obj instanceof java.nio.file.Path) {
        return obj.toAbsolutePath().toString()
    } else if (obj instanceof Map) {
        return obj.collectEntries { k, v -> [k, convertPathsToStrings(v)] }
    } else if (obj instanceof Collection) {
        return obj.collect { it -> convertPathsToStrings(it) }
    } else {
        return obj
    }
}

def createJsonWorkflowOutput(_input = null) {
    // _input is needed because a closure is provided to all functions called in the workflow
    return [
        toJson: { obj ->
            def filtered = removeNullValues(obj)
            return groovy.json.JsonOutput.toJson(filtered)
        }
    ]
}

def removeNullValues(obj) {
    if (obj instanceof Map) {
        return obj.findAll { _k, v -> v != null }.collectEntries { k, v -> [k, removeNullValues(v)] }
    } else if (obj instanceof Collection) {
        return obj.findAll { it -> it != null }.collect { it -> removeNullValues(it) }
    } else {
        return obj
    }
}