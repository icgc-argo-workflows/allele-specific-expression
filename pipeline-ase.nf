nextflow.enable.dsl=2

scripts_folder = "${projectDir}/scripts/"
ase_cleanup_path = scripts_folder + "ase_cleanup.py"
gene_annotation_path = scripts_folder + "gene_annotation.py"
hap_table_path = scripts_folder + "hap_table.py"
plot_VAF_path = scripts_folder + "plot_VAF.py"

params.min_depth = 16
params.min_mapping_quality = 20
params.min_base_quality = 10
params.min_mappability = 0.05
params.out_dir ="./out"
params.debug_out = ""

clean_suffix = params.containsKey("gtf") ? "clean" : "tsv"

process bgzip 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""
    
    input:
    path(fa)

    output:
    path("*.gz")
    
    """
    bgzip $fa 
    """  
}

process samIndex 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    tuple(val(name), path(bam))

    output:
    tuple(val(name), path(bam), path("*.bai"))
    
    """
    samtools index $bam 
    """
}

process gatkDict 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    path(fa)

    output:
    path("*.dict")
    
    """
    gatk CreateSequenceDictionary -R $fa 
    """
}

process samFaidx 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    path(fa)

    output:
    path("*.fai"), emit: fai
    path("*.gzi"), emit: gzi

    """
    samtools faidx $fa
    """
}

// For each pair (chr, pos) keeps only the first occurence if there are multiple. Only keep the bi-allelic positions
process filterVcf {
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input: 
    tuple val(name), path(vcf_file)

    output:
    tuple val(name), path("${vcf_file.baseName}.filtered.vcf")
    
    script:
    cat_cmd = "$vcf_file".endsWith(".gz") ? "zcat" : "cat"
    """
    $cat_cmd $vcf_file |
    cut -f 1,2 | 
    sort -k1,1V -k2,2n |
    uniq -u > chr_file      
    bcftools view $vcf_file --min-alleles 2 --max-alleles 2 -T chr_file -O v -o ${vcf_file.baseName}.filtered.vcf
    """
}   

process gatkFeatureIndex 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    tuple val(name), path(vcf_file)

    output:
    tuple val(name), path(vcf_file), path("*.idx")
    
    """
    gatk IndexFeatureFile -I $vcf_file
    """
}

process aseReadCounterPair
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input: 
    tuple val(name), path(bam), path(bai), path(vcf), path(idx)
    path(fa)
    path(fai)
    path(gzi)
    path(dict)

    output:
    tuple(val(name), path("${name}.read"))

    """
    gatk ASEReadCounter -R $fa -I $bam -V $vcf -O ${name}.read \
        --min-depth-of-non-filtered-base $params.min_depth \
        --min-mapping-quality $params.min_mapping_quality \
        --min-base-quality $params.min_base_quality
    """
}

process faGenome 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""
    
    input:
    path(fai)

    output:
    path("${fai.baseName}.genome")

    """
    awk -v OFS=\"\\t\" \"{print \\\$1,\\\$2}\" $fai | sort -V -s > ${fai.baseName}.genome
    """
}

process bedFile 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    tuple(val(name), path(ase))

    output:    
    tuple(val(name), path("${name}.bed"))

    """   
    awk -v OFS=\"\\t\" \"{print \\\$1,\\\$2-1,\\\$2 }\" $ase | tail -n+2 | sort -V -k1,1 -s > ${name}.bed
    """
}

process bedtoolsMap
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    path(genome)
    path(mapp_path)
    tuple(val(name), path(bed))
    
    output:
    tuple(val(name), path("${name}.mapp"))
 
    """    
    echo -e \"contig\tpos\tmappability\" > ${name}.mapp
    bedtools map -a $bed -b $mapp_path -o min -c 4 -g $genome | cut -f1,3,4 >> ${name}.mapp
    """
}

process aseCleanUpMapp 
{
    publishDir params.out_dir, mode: "copy"

    input:
    tuple(val(name), path(ase), path(mapp))
    path(ase_cleanup_py)

    output:
    tuple val(name), path("*." + clean_suffix), emit: ase_file
    path("*.log"), emit: ase_log

    script:  
    """
    python $ase_cleanup_py --ase $ase --min_SNP_depth $params.min_depth  --output ${ase.baseName}.${clean_suffix} --mappability $mapp --filter_mapp $params.min_mappability
    mv ase_cleanup.log ${ase.baseName}.ase.log
    """
}

process aseCleanUp 
{
    publishDir params.out_dir, mode: "copy"

    input:
    tuple(val(name), path(ase))
    path(ase_cleanup_py)

    output:
    tuple val(name), path("*." + clean_suffix), emit: ase_file
    path("*.log"), emit: ase_log

    script:  
    """
    python $ase_cleanup_py --ase $ase --min_SNP_depth $params.min_depth  --output ${ase.baseName}.${clean_suffix}
    mv ase_cleanup.log ${ase.baseName}.ase.log
    """
}

process aseGeneAnnotation 
{
    publishDir params.out_dir, mode: "copy"

    input: 
    path(gene_annotation_py)
    tuple val(name), path(ase_clean)
    path(gtf)
    path(gtf_db)
    val assembly

    output: 
    tuple val(name), path("${name}.tsv"), emit: gene_file
    path("*.log"), emit: gene_log
    
    """
    python $gene_annotation_py -I $ase_clean -O ${name}.hap --gtf $gtf --ref $assembly
    mv ${name}.hap ${name}.tsv
    mv gene_annotation.log ${name}.hap.log
    """
}

process aseVAFPlot 
{
    publishDir params.out_dir, mode: "copy"

    input: 
    tuple val(name), path(ase_clean)

    output: 
    path("${name}.vaf.png")

    """
    $plot_VAF_path -I $ase_clean -O ${name}.vaf.png
    """
}

process computeHapTable
{
    publishDir params.out_dir, mode: "copy"

    input:
    path(hap_table_py)
    tuple val(name), path(ase_gene_file), path(vcf_file)
    
    output:
    path("${name}.hap.tsv")

    """    
    python $hap_table_py -I $ase_gene_file -V $vcf_file -O ${name}.hap.tsv
    """
}

process createGtfDB
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    path(gtf_path)
    val(assembly)
    val(filename)

    output:
    path("*.db")

    """
    #!/usr/bin/env python
    
    import pyensembl
    data = pyensembl.Genome(reference_name='$assembly', annotation_name='genome_annotation', gtf_path_or_url='$gtf_path')
    data.index()
    """
}

def getVcfFiles(bam_names) 
{
    if (params.containsKey("vcf")) 
    {
        vcf_files = Channel.fromPath(params.vcf)
        vcf_files.ifEmpty{error "No VCF files found matching the naming of the bam files."}
        vcf = bam_names.combine(vcf_files)
        vcf.ifEmpty{error "Failed to match VCF files onto BAM files."}
    }
    else 
    {
        vcf_filename = params.bam - ".bam" + ".vcf"
        vcf = Channel.fromFilePairs(vcf_filename, size: 1)
    }
    return vcf
}

workflow 
{
    bam = Channel.fromFilePairs(params.bam, size: 1)
    bam_names = bam.map{name, file -> name}
    vcf_raw = getVcfFiles(bam_names)
    ref_fa = file(params.fa)
    fa_is_zipped = params.fa.endsWith(".gz")
    fa = fa_is_zipped ? ref_fa : bgzip(ref_fa)
    ase_clean_script = file(ase_cleanup_path)

    // Create indices
    bam_bai = samIndex(bam)
    faidx = samFaidx(fa)
    dict = gatkDict(fa)
    vcf_filtered = filterVcf(vcf_raw)
    idx = gatkFeatureIndex(vcf_filtered)
    
    // merge data
    joint_files = bam_bai.join(idx)

    // Call read counter
    ase_count = aseReadCounterPair(joint_files, fa, faidx.fai, faidx.gzi, dict)
    
    // Call cleaner (with optinal mapp input)
    if (params.containsKey("mapp")) // Apply mappability
    {
        bedgraph_mapp = file(params.mapp)
        genome = faGenome(faidx.fai)
        bed = bedFile(ase_count)
        bed_mapp = bedtoolsMap(genome, bedgraph_mapp, bed)
        mapped_files = ase_count.join(bed_mapp)
        ase_clean = aseCleanUpMapp(mapped_files, ase_clean_script)
    } 
    else 
    {
        ase_clean = aseCleanUp(ase_count, ase_clean_script)
    }
    
    vaf_plot = aseVAFPlot(ase_clean.ase_file)

    if (params.containsKey("gtf")) 
    {        
        if (params.gtf.endsWith(".gz")) 
        {
            gtf_db_path = params.gtf - ~"\\.gz\$" + ".db"
        }
        else 
        {
            gtf_db_path = params.gtf - ~"\\.gtf\$" + ".db"
        }
        gene_table_script = file(gene_annotation_path)
        hap_table_script = file(hap_table_path)
        gtf_file_ch = file(params.gtf)

        if (file(gtf_db_path).exists()) 
        {
            gtf_db_ch = file(gtf_db_path)
        }
        else 
        {
            gtf_db_ch = createGtfDB(gtf_file_ch, params.assembly, file(gtf_db_path).name)
            gtf_db_ch.subscribe 
            {
                file(gtf_db_path).bytes = it.bytes  
            }
        }
        ase_gene_annotation = aseGeneAnnotation(gene_table_script, ase_clean.ase_file, gtf_file_ch, gtf_db_ch, params.assembly)          
        gene_merge = ase_gene_annotation.gene_file.join(vcf_filtered)
        computeHapTable(hap_table_script, gene_merge)
    } 
}
