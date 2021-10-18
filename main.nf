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
params.outdir = "$baseDir/results/"
params.genome = "$baseDir/data/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"

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

process qualityControl{
    output:
    stdout into randNums1

    """
    fastqc ${params.reads}
    """
}

/*
 * Step 1. Builds the genome index required by the mapping process
 */
process trimmingReads {
    output:
    stdout into randNums2

    """
    cutadapt -q 33 -o ${params.outdir}/out1.fastq -p ${params.outdir}/out2.fastq ${params.reads}
    """
}

/*
 * Step 1. Builds the genome index required by the mapping process
 */
process genomeIndexing {
    input:
    path genome from params.genome

    output:
    stdout into randNums3

    """
    bowtie2-build --threads ${task.cpus} ${genome} ${params.outdir}/genome.index
    """
}

/*
 * The below print the output from the above steps; I use them for debugging.
 */
process printstep1 { 
    echo true
 
    input:
    stdin from randNums1

    '''
    #!/usr/bin/env python
    import sys
    
    for line in sys.stdin:
        #print(line)
        pass
    '''
}

process printstep2 { 
    echo true
 
    input:
    stdin from randNums2

    """
    echo $randNums2
    """
}

process printstep3 { 
    echo true
 
    input:
    stdin from randNums3

    '''
    #!/usr/bin/env python
    import sys
    
    for line in sys.stdin:
        #print(line)
        pass
    '''
}