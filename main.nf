def helpMessage() {
    log.info"""
    =========================================
     Nextflow Pipeline Demo v${params.author}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf --fastqs '/project/*_{1,2}*.fastq' --outdir '/project/'

    Required arguments:
         --fastqs                      Directory pattern for fastq files: /project/*{R1,R2}*.fastq (Required if --sras not specified)
         --sras                        Directory pattern for SRA files: /project/*.sras (Required if --fastqs not specified)

    Input File options:
        --singleEnd                    Specifies that the input files are not paired reads (default is paired-end).

    Save options:
        --outdir                       Specifies where to save the output from the nextflow run.
        --savefq                       Compresses and saves raw fastq reads.
        --saveTrim                     Compresses and saves trimmed fastq reads.
        --skipBAM                      Skip saving BAM files. Only CRAM files will be saved with this option.
        --saveAll                      Compresses and saves all fastq reads.
        
    QC Options:
        --skipRSeQC                    Skip running RSeQC.
        
    Analysis Options:
        --noTrim                       Skip trimming and map only. Will also skip flip/flipR2 (any BBMap) steps.    
        --count                        Run RSeQC FPKM count over RefSeq annotated genes.

    """.stripIndent()
}

// Configurable variables
params.fastqs = "$baseDir/data/fastq/*_{1,2}.fastq"
params.genome = "$baseDir/data/genome/*.fa"
params.genome_refseq = "$baseDir/data/genome/*.bed.gz"
params.outdir = "$baseDir/results"
params.hisat2_indices = "$baseDir/data/hisat_indices/"
params.author = "Ashwini"
// params.singleEnd = false

// Validate inputs
// TODO: If the input directories are not present as in the above structure, write a script to create them.
if ( params.genome ){
    genome = file(params.genome)
    if( genome.size() == 0 ) exit 1, "Genome directory not found: ${params.genome}"
}

if ( params.hisat2_indices ){
    hisat2_indices = file("${params.hisat2_indices}")
}

if ( params.genome_refseq ){
    genome_refseq = file("${params.genome_refseq}")
}

/*
 * Create a channel for input read files
 */
if (params.fastqs) {
    if (params.singleEnd) {
        println("SingleEnd True")
        fastq_reads_qc = Channel
                            .fromPath(params.fastqs)
                            .map { file -> tuple(file.simpleName, file) }
        fastq_reads_trim = Channel
                            .fromPath(params.fastqs)
                            .map { file -> tuple(file.simpleName, file) }
        fastq_reads_subsample = Channel
                            .fromPath(params.fastqs)
                            .map { file -> tuple(file.simpleName, file) }        
    } else {
        println("SingleEnd False")
        println(params.fastqs)
        Channel
            .fromFilePairs( params.fastqs, size: params.singleEnd ? 1 : 2 )
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.fastqs}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
            .into { fastq_reads_qc; fastq_reads_trim; fastq_reads_subsample }
    }
}

else {
    Channel
        .empty()
        .into { fastq_reads_qc; fastq_reads_trim }
}

if (params.sras) {
    println("pattern for SRAs provided")
    read_files_sra = Channel
                        .fromPath(params.sras)
                        .map { file -> tuple(file.baseName, file) }
}

else {
    read_files_sra = Channel.empty()
}

/*
 * Parse software version numbers
 */
process get_software_versions {
    validExitStatus 0,1,127
    publishDir "${params.outdir}/software_versions/", mode: 'copy', pattern: '*.txt'

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file '*.txt' into software_versions_text

    script:
    """
    echo $params.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    hisat2 --version > v_hisat2.txt
    samtools --version > v_samtools.txt
    fastq-dump --version > v_fastq-dump.txt
    igv version > v_igv-tools.txt

    for X in `ls *.txt`; do
        cat \$X >> all_versions.txt;
    done
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

/*
 * Step 1 -- get fastq files from downloaded sras
 */

process sra_dump {
    tag "$prefix"
    if (params.threadfqdump) {
        cpus 4 }
    else {
        cpus 1
    }
    if (params.savefq || params.saveAllfq) {
        publishDir "${params.outdir}/fastq", mode: 'copy'
    }
    
    input:
    set val(prefix), file(reads) from read_files_sra

    output:
    set val(prefix), file("*.fastq.gz") into fastq_reads_qc_sra, fastq_reads_trim_sra, fastq_reads_hisat2_notrim_sra

    script:
    prefix = reads.baseName
    println("#154 prefix")
    println(prefix)
    println("#156 params.threadfqdump")
    println(params.threadfqdump)
    println("#158 params.singleEnd")
    println(params.singleEnd)
    if (!params.threadfqdump) {
        """
        echo ${prefix}
        fastq-dump ${reads} --gzip
        """
    } else if (!params.singleEnd) {
        println("Running this")
         """
        export PATH=~/.local/bin:$PATH
        parallel-fastq-dump \
            --threads 8 \
            --gzip \
            --split-3 \
            --outdir ${params.outdir}/sras \
            --sra-id ${reads}
        """
    } else if (!params.threadfqdump && !params.singleEnd) {
        """
        echo ${prefix}
        fastq-dump --split-3 ${reads} --gzip
        """
    } else {
        """
        export PATH=~/.local/bin:$PATH
        parallel-fastq-dump \
            --threads 8 \
            --gzip \
            --sra-id ${reads}
        """
    }
}

/*
 * Step 2 -- Perform qualityControl on the fastq files
 */

process qualityControl {
    output:
    stdout into outputString1

    """
    fastqc ${params.fastqs}
    """
}

/*
 * Step 3 -- Trim the fastq files
 */

process trimmingReads {
    tag "$prefix"

    if (params.sras) {
        println("#211 fastq_reads_trim.mix(fastq_reads_qc_sra)")
        println("#212 fastq_reads_trim.mix(fastq_reads_trim_sra)")
        println("#213 fastq_reads_trim.mix(fastq_reads_hisat2_notrim_sra)")
    } else {
        println("#215 params.fastqs")
        println(params.fastqs)
    }

    input:
    set val(name), file(reads) from fastq_reads_trim.mix(fastq_reads_trim_sra)

    output:
    stdout into outputString2
    
    """
    echo name${name}
    echo reads${reads}
    cutadapt -q 33 \
    -o ${params.outdir}/out1.fastq \
    -p ${params.outdir}/out2.fastq ${params.fastqs}
    
    """
}

/*
 * Step 4 -- Perform post trim qualityControl on the fastq files
 */

process postTrimQualityControl{
    output:
    stdout into outputString3

    """
    fastqc ${params.outdir}/*.fastq
    """
}

/*
 * Step 5 -- Perform Read mapping
 */

 process readMapping{

    input:
    file(indices) from hisat2_indices
    val(indices_path) from hisat2_indices   

    output:
    set val(name), file("*.sam") into hisat2_sam
    
    script:
    prefix_pe = reads[0].toString() - ~/(_1\.)?(_R1)?(\.fq)?(fq)?(\.fastq)?(fastq)?(\.gz)?$/
    prefix_se = reads[0].toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/

    if (!params.singleEnd){
     """             
        hisat2  -p 6 \
        -x ${hisat2_indices} \
        -1 ${prefix_pe}_1.fastq.gz \
        -2 ${prefix_pe}_2.fastq.gz 
        -S ${${prefix_pe}}aligned.sam
    
    """
    }else if(params.singleEnd){
    """             
        hisat2 -p 6 \
        -x ${hisat2_indices} \
        -U ${prefix_se}.fastq.gz  
        -S ${${prefix_pe}}aligned.sam
    
    """
    }
   
}
// Convert sam to bam

process samToBam {
    input:
    set val(name), file(mapped_sam) from hisat2_sam

    output:
    set val(name), file("${name}.sorted.bam") into sorted_bam_ch
    set val(name), file("${name}.sorted.bam.bai") into sorted_bam_indices_ch

    """
    samtools view -@ 16 -bS -o ${name}.bam ${mapped_sam}
    samtools sort -@ 16 ${name}.bam > ${name}.sorted.bam
    samtools flagstat ${name}.sorted.bam > ${name}.flagstat
    samtools view -@ 16 -F 0x40 ${name}.sorted.bam | cut -f1 | sort | uniq | wc -l > ${name}.millionsmapped
    samtools index ${name}.sorted.bam ${name}.sorted.bam.bai
    samtools view -@ 16 -C -T ${genome} -o ${name}.cram ${name}.sorted.bam
    samtools sort -@ 16 -O cram ${name}.cram > ${name}.sorted.cram
    samtools index -c ${name}.sorted.cram ${name}.sorted.cram.crai
    """
}

//  Analyze read distributions using RSeQC

process rseqc {
    when:
    !params.skipRSeQC

    input:
    set val(name), file(bam_file) from sorted_bam_ch
    file(bam_indices) from sorted_bam_indices_ch

    output:
    file "*.{txt,pdf,r,xls}" into rseqc_results

    script:
    """
    
    export PATH=~/.local/bin:$PATH
    
    read_distribution.py -i ${bam_file} \
                         -r ${genome_refseq} \
                         > ${name}.read_distribution.txt
    read_duplication.py -i ${bam_file} \
                        -o ${name}.read_duplication
    infer_experiment.py -i ${bam_file} \
                        -r ${genome_refseq} \
                        > ${name}.infer_experiment.txt
                        
    bam_stat.py -i ${bam_file} \
                        > ${name}.bam_stat.txt
                        
    junction_annotation.py -i ${bam_file} \
                           -o ${name}.junction_annotation \
                           -r ${genome_refseq}
                           
    junction_saturation.py -i ${bam_file} \
                           -o ${name}.junction_saturation \
                           -r ${genome_refseq}                       
    """
 }


 process rseqc_count {
    publishDir "${params.outdir}/rseqc_counts" , mode: 'copy', pattern: "*FPKM.xls"
    
    when:
    params.count

    input:
    set val(name), file(bam_file) from sorted_bam_ch
    file(bam_indices) from sorted_bam_indices_ch

    output:
    file "*FPKM.xls" into rseqc_counts

    script:
    """
    export PATH=~/.local/bin:$PATH
    
    FPKM_count.py -i ${bam_file} \
                           -o ${name}.counts \
                           -r ${genome_refseq} \
                           -e
    """
 }

// Assign aligned reads to genes

// process htseq{
// """
// htseq-count -f bam -r name -i gene_id -t exon -m intersection-nonempty ${algnmentfile} ${gtffile} > ${textfile}
// """
// }