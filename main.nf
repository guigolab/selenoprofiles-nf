#!/usr/bin/env nextflow

/*
 * @authors
 * Ferriol Calvet <ferriol.calvet@crg.eu>
 * Emilio Righi <emilio.righi@crg.eu>
 */

/*
 * Input parameters: genome, protein evidences, parameter file,
 * additional values for the generation of the parameter file.
 * Params are stored in the params.config file
 *
 * The configuration is in nextflow.config file
 */

nextflow.enable.dsl=2


process PARSE_FASTA_HEADER {

    tag "${id}"

    input:
    tuple val(id), val(taxid), val(scientific_name),path(fasta)

    output:
    tuple val(id), val(taxid),val(scientific_name), path(out_fasta)

    script:
    out_fasta = "parsed_${fasta.BaseName}.fa"
    
    """
    #!/usr/bin/env python3

    with open("${fasta}", 'r') as infile, open("${out_fasta}", 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                header = line.strip()
                if 'ENA' in header:
                    parts = header.split('|')
                    if len(parts) > 2:
                        new_header = f'>{parts[-1].split()[0]}\\n'
                        outfile.write(new_header)
                    else:
                        outfile.write(header + '\\n')
                else:
                    outfile.write(header + '\\n')
            else:
                outfile.write(line)
    """
}

process UNZIP_FASTA {

    tag "${id}"

    input:
    tuple val(id), val(taxid), val(scientific_name), path(genome)

    output:
    tuple val(id), val(taxid), val(scientific_name), path(unzipped_genome)

    script:
    unzipped_genome = "${id}_unzipped.fa"
    """
    gunzip -c ${genome} > ${unzipped_genome};
    """
}

process RUN_SELENOPROFILES {

    container "sp:latest"

    publishDir "${params.outdir}", mode: 'copy'

    tag "${id}"

    input:
    tuple val(id),val(taxid),val(scientific_name),path(genome)

    output:
    path("${taxid}/*", type:'dir')

    script:
    genome_name = genome.BaseName

    """
    selenoprofiles -o ./${taxid} -t ${genome} -s "${scientific_name}" -p ${params.selenoprofile}
    """
}

workflow {

    if (params.tsv) { tsv_input = file(params.tsv) } else {
        exit 1, 'TSV samplesheet not specified!'
    }

    tsv = channel.fromPath(params.tsv).splitCsv( sep: '\t', header:true )
    .collectFile(storeDir:params.assemblies_dir){ row -> ["${row[params.column_id_value]}-${row[params.column_taxid_value]}-${row[params.column_scientific_name]}-.fa.gz", file(row[params.column_path_value])]}
    .map { it -> 
        def elements = it.baseName.tokenize('-')
        tuple(elements[0],elements[1],elements[2], it)
         }
    .branch { 
        validFasta: it[3].countFasta() > 0
        invalidFasta: true 
        } 

    invalid_fasta = tsv.invalidFasta
    if(invalid_fasta.count()){
        invalid_fasta.count().view { it -> "A total of ${it} invalid assemblies have been found"}
        invalid_fasta.collectFile(storeDir:params.outdir){ item ->
            [ "INVALID_ASSEMBLIES.txt", "${item[0]}\t${item[1]}\t${item[2]}" + '\n' ]
        }
        .subscribe {
            println "Invalid assemblies have been saved in: $it"
        }
    } 

    assemblies = tsv.validFasta.map{it-> tuple(it[0], it[1], it[2], it[3])}

    results = UNZIP_FASTA(assemblies) | PARSE_FASTA_HEADER | RUN_SELENOPROFILES

}


workflow.onComplete {
	println ( "\nDone!\n" )
}

