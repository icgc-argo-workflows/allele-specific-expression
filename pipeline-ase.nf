nextflow.enable.dsl=2

params.min_depth = 16
params.min_mapping_quality = 20
params.min_base_quality = 10
params.min_mappability = 0.05
params.out_dir ="./out"
params.debug_out = ""
params.vcf_sample = null

process collectVersions
{
    publishDir params.out_dir, mode: "copy"

    input:
    val(pipeline_version)
    val(nextflow_version)

    output:
    path("software_versions.yml")

    script:
    """
    python_version=\$(python --version 2>&1)
    samtools_version=\$(samtools --version | head -n 1)
    bcftools_version=\$(bcftools --version | head -n 1)
    bedtools_version=\$(bedtools --version)
    gatk_version=\$(gatk --version 2>&1 | tail -n 1)

    printf '%s\n' \
        "pipeline_ase: '${pipeline_version}'" \
        "nextflow: '${nextflow_version}'" \
        "python: '\${python_version}'" \
        "samtools: '\${samtools_version}'" \
        "bcftools: '\${bcftools_version}'" \
        "bedtools: '\${bedtools_version}'" \
        "gatk: '\${gatk_version}'" \
        > software_versions.yml
    """
}

process bgzip 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""
    
    input:
    path(fa)

    output:
    path("*.gz")

    script:
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

    script:
    """
    samtools index $bam 
    """
}

process samToBam
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    tuple val(name), path(sam)

    output:
    tuple val(name), path("${name}.bam")

    script:
    """
    samtools sort -o ${name}.bam $sam
    """
}

process gatkDict 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    path(fa)

    output:
    path("*.dict")

    script:
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

    script:
    """
    samtools faidx $fa
    """
}

process validateInputs
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    tuple val(name), path(bam), path(bai), path(vcf)
    path(fa)
    path(fai)
    path(gzi)
    path(validate_inputs_py)

    output:
    tuple val(name), path(vcf), emit: vcf
    tuple val(name), path("${name}.vcf_sample.txt"), emit: sample

    script:
    def sample_arg = params.vcf_sample ? "--sample '${params.vcf_sample}'" : ""
    """
    samtools quickcheck -v $bam
    samtools view -H $bam > ${name}.header.sam
    python $validate_inputs_py \
        --alignment-header ${name}.header.sam \
        --variants $vcf \
        --reference-index $fai \
        --bam-name '${name}' \
        --sample-output ${name}.vcf_sample.txt \
        $sample_arg
    bcftools norm \
        --check-ref e \
        --fasta-ref $fa \
        --output-type u \
        --output /dev/null \
        $vcf
    """
}

// Remove duplicated coordinates and retain only biallelic SNV records.
process filterVcf {
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input: 
    tuple val(name), path(vcf_file)

    output:
    tuple val(name), path("${vcf_file.baseName}.filtered.vcf")
    
    script:
    def cat_cmd = "$vcf_file".endsWith(".gz") ? "zcat" : "cat"
    """
    $cat_cmd $vcf_file |
    cut -f 1,2 | 
    sort -k1,1V -k2,2n |
    uniq -u > chr_file      
    bcftools view $vcf_file --types snps --min-alleles 2 --max-alleles 2 -T chr_file -O v -o ${vcf_file.baseName}.filtered.vcf
    """
}   

process gatkFeatureIndex 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    tuple val(name), path(vcf_file)

    output:
    tuple val(name), path(vcf_file), path("*.idx")

    script:
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

    script:
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

    script:
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

    script:
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

    script:
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
    val(output_suffix)

    output:
    tuple val(name), path("*.${output_suffix}"), emit: ase_file
    path("*.log"), emit: ase_log

    script:  
    """
    python $ase_cleanup_py --ase $ase --min_SNP_depth $params.min_depth  --output ${ase.baseName}.${output_suffix} --mappability $mapp --filter_mapp $params.min_mappability
    mv ase_cleanup.log ${ase.baseName}.ase.log
    """
}

process aseCleanUp 
{
    publishDir params.out_dir, mode: "copy"

    input:
    tuple(val(name), path(ase))
    path(ase_cleanup_py)
    val(output_suffix)

    output:
    tuple val(name), path("*.${output_suffix}"), emit: ase_file
    path("*.log"), emit: ase_log

    script:  
    """
    python $ase_cleanup_py --ase $ase --min_SNP_depth $params.min_depth  --output ${ase.baseName}.${output_suffix}
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
    path("${name}.gene_annotation.log"), emit: gene_log

    script:
    """
    python $gene_annotation_py -I $ase_clean -O ${name}.hap --gtf $gtf --ref $assembly
    mv ${name}.hap ${name}.tsv
    mv gene_annotation.log ${name}.gene_annotation.log
    """
}

process aseVAFPlot 
{
    publishDir params.out_dir, mode: "copy"

    input: 
    tuple val(name), path(ase_clean)
    path(plot_vaf_py)

    output: 
    path("${name}.vaf.png")

    script:
    """
    python $plot_vaf_py -I $ase_clean -O ${name}.vaf.png
    """
}

process computeHapTable
{
    publishDir params.out_dir, mode: "copy"

    input:
    path(hap_table_py)
    tuple val(name), path(ase_gene_file), path(vcf_file), path(sample_file)
    
    output:
    tuple val(name), path("${name}.hap.tsv"), emit: hap_file
    path("${name}.hap.log"), emit: hap_log

    script:
    """
    python $hap_table_py -I $ase_gene_file -V $vcf_file -O ${name}.hap.tsv --bam-name '${name}' --sample-file $sample_file
    mv hap_table.log ${name}.hap.log
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
    path(filename)

    script:
    """
    #!/usr/bin/env python
    
    import os
    import shutil
    import pyensembl

    output_path = '$filename'
    data = pyensembl.Genome(
        reference_name='$assembly',
        annotation_name='genome_annotation',
        gtf_path_or_url='$gtf_path',
        cache_directory_path='.')
    data.index()

    database_path = data.db.local_db_path
    if os.path.abspath(database_path) != os.path.abspath(output_path):
        shutil.copy2(database_path, output_path)

    if not os.path.isfile(output_path):
        raise RuntimeError(f'PyEnsembl database was not created: {output_path}')
    """
}

workflow
{
    if (!params.bam)
        error "Missing required parameter: --bam"
    if (!params.fa)
        error "Missing required parameter: --fa"
    if (params.containsKey("gtf") && !params.assembly)
        error "--assembly is required when --gtf is provided"

    collectVersions(workflow.manifest.version ?: "dev", workflow.nextflow.version)

    clean_suffix = params.containsKey("gtf") ? "clean" : "tsv"
    alignment = channel.fromFilePairs(params.bam, size: 1, checkIfExists: true)
        .map { name, files -> tuple(name, files[0]) }
    bam_names = alignment.map{name, _file -> name}
    if (params.containsKey("vcf"))
    {
        vcf_files = channel.fromPath(params.vcf, checkIfExists: true)
        vcf_files.ifEmpty{error "No VCF files found matching the naming of the bam files."}
        vcf_raw = bam_names.combine(vcf_files.collect()).map { combined ->
            def name = combined[0]
            def candidates = combined[1..-1]

            if (candidates.size() == 1)
            {
                return tuple(name, candidates[0])
            }

            def matches = candidates.findAll { candidate ->
                candidate.name.replaceFirst(/(?i)\.vcf(\.gz)?$/, '') == name
            }
            if (matches.size() != 1)
            {
                error "Expected exactly one VCF matching BAM '${name}', found ${matches.size()}."
            }
            tuple(name, matches[0])
        }
    }
    else
    {
        vcf_filename = params.bam.replaceFirst(/(?i)\.(bam|sam)$/, '.vcf')
        vcf_raw = channel.fromFilePairs(vcf_filename, size: 1, checkIfExists: true)
            .map { name, files -> tuple(name, files[0]) }
    }
    ref_fa = file(params.fa, checkIfExists: true)
    fa_is_zipped = params.fa.endsWith(".gz")
    fa = fa_is_zipped ? ref_fa : bgzip(ref_fa)
    ase_clean_script = file("${projectDir}/scripts/ase_cleanup.py", checkIfExists: true)
    validate_inputs_script = file("${projectDir}/scripts/validate_inputs.py", checkIfExists: true)

    bam_input = alignment.filter { _name, alignment_file ->
        alignment_file.name.toLowerCase().endsWith(".bam")
    }
    sam_input = alignment.filter { _name, alignment_file ->
        alignment_file.name.toLowerCase().endsWith(".sam")
    }
    converted_sam = samToBam(sam_input)
    bam = bam_input.mix(converted_sam)

    // Create indices
    bam_bai = samIndex(bam)
    faidx = samFaidx(fa)
    dict = gatkDict(fa)
    validation_input = bam_bai.join(vcf_raw)
    validation_result = validateInputs(
        validation_input,
        fa,
        faidx.fai,
        faidx.gzi,
        validate_inputs_script)
    vcf_filtered = filterVcf(validation_result.vcf)
    idx = gatkFeatureIndex(vcf_filtered)
    
    // merge data
    joint_files = bam_bai.join(idx)

    // Call read counter
    ase_count = aseReadCounterPair(joint_files, fa, faidx.fai, faidx.gzi, dict)
    
    // Call cleaner (with optional mappability input)
    if (params.containsKey("mapp")) // Apply mappability
    {
        bedgraph_mapp = file(params.mapp)
        genome = faGenome(faidx.fai)
        bed = bedFile(ase_count)
        bed_mapp = bedtoolsMap(genome, bedgraph_mapp, bed)
        mapped_files = ase_count.join(bed_mapp)
        ase_clean = aseCleanUpMapp(mapped_files, ase_clean_script, clean_suffix)
    } 
    else 
    {
        ase_clean = aseCleanUp(ase_count, ase_clean_script, clean_suffix)
    }
    
    aseVAFPlot(ase_clean.ase_file, file("${projectDir}/scripts/plot_VAF.py", checkIfExists: true))

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
        gene_table_script = file("${projectDir}/scripts/gene_annotation.py", checkIfExists: true)
        hap_table_script = file("${projectDir}/scripts/hap_table.py", checkIfExists: true)
        gtf_file_ch = file(params.gtf, checkIfExists: true)

        if (file(gtf_db_path).exists()) 
        {
            gtf_db_ch = file(gtf_db_path)
        }
        else 
        {
            gtf_db_ch = createGtfDB(gtf_file_ch, params.assembly, file(gtf_db_path).name)
        }
        ase_gene_annotation = aseGeneAnnotation(gene_table_script, ase_clean.ase_file, gtf_file_ch, gtf_db_ch, params.assembly)          
        gene_merge = ase_gene_annotation.gene_file.join(vcf_filtered).join(validation_result.sample)
        computeHapTable(hap_table_script, gene_merge)
    } 
}
