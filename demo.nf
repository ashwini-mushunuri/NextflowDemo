/* 
 * Proof of concept Nextflow pipeline
 * 
 * Authors:
 * Ashwini Mushunuri 
 */ 

 
/*
 * Defines some parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.reads = "$baseDir/data/*_{1,2}.fastq"
params.outdir = 'results'

params.randNums = 10

log.info """\
         NEXTFLOW    PIPELINE
         =============================
         reads : ${params.reads}
         outdir: ${params.outdir}
         """
         .stripIndent()

/*
 * Create the `read_pairs_ch` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file 
 */
Channel
    .fromPath( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch }

/*
 * Step 1. Builds the genome index required by the mapping process
 */
/*
 * A Python script task which parses the output of the previous script
 */

process t1{
    output:
    stdout into randNums

    """
    echo 10
    """
}

process pyTask { 
    echo true
 
    input:
    stdin from randNums

    '''
    #!/usr/bin/env python
    import sys
    import random
 
    x = 0
    y = 0
    lines = list()
    for line in sys.stdin:
        print(line)
    #for _ in range(sys.stdin):
    #    lines.append(random.randint(1, 100))
    #print(f"sum: {sum(lines)}, avg: {sum(lines)/len(lines)}")
    '''
}