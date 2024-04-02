//
nextflow.enable.dsl=2

// Define the default parameters
params.dragen_path = ''
params.target_bed = ''
params.snpeff_db = 'GRCh38.99'
params.reference_path = '/home/ubuntu/annotation/fasta/hg38.fa'
params.outdir = './'
params.samtools_path = '/home/ubuntu/bin/samtools'
params.rmsk_path = '/home/ubuntu/annotation/hg38.rmsk.bed'

process SNV_VARIANT_ANNOTATE {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    val db
    path target_bed

    output:
    tuple val(meta), path("*.ann.vcf"), emit: vcf
    tuple val(meta), path("*.ann.high_risk.csv"), emit: table

    script:
    """
    snpeff -canon -noInteraction -motif $db -filterInterval $target_bed -stats ${meta}.summary.html $vcf > ${meta}.ann.vcf

    # only keep the variants with PASS and high risk
    bcftools view ${meta}.ann.vcf -i 'FILTER="PASS" && ANN ~ "HIGH"' \\
        | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t[%GT]\\t%QD\\t%ANN\\n' \\
        | awk '{
            OFS="\\t";
            split(\$9,ann,",");
            annotation="";
            for(i in ann)
            {
                split(ann[i],f,"|");
                if(f[3]=="HIGH")
                {
                    annotation=annotation"\\t"f[4]"\\t"f[1]"|"f[2]"|"f[3]
                }
            }
            print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8annotation
            }' > ${meta}.ann.high_risk.csv
    """

}

process CNV_VARIANT_ANNOTATE {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    val db
    path target_bed

    output:
    tuple val(meta), path("*.ann.vcf"), emit: vcf
    tuple val(meta), path("*.ann.high_risk.csv"), emit: table

    script:
    """
    snpeff -canon -noInteraction -motif $db -filterInterval $target_bed -stats ${meta}.summary.html $vcf > ${meta}.cnv.ann.vcf

    # only keep the variants with PASS and high risk
    bcftools view ${meta}.cnv.ann.vcf -i 'FILTER="PASS" && SVTYPE="CNV" && ANN ~ "HIGH"' \\
        | bcftools query -f '%ID\\t%CHROM\\t%POS\\t%END\\t%ALT\\t%QUAL\\t%FILTER\\t[%GT]\\t%ANN\\n' \\
        | awk '{
            OFS="\\t";
            split(\$9,ann,",");
            annotation="";
            for(i in ann)
            {
                split(ann[i],f,"|");
                if(f[3]=="HIGH")
                {
                    annotation=annotation"\\t"f[4]"\\t"f[1]"|"f[2]"|"f[3]
                }
            }
            print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8annotation
            }' > ${meta}.cnv.ann.high_risk.csv
    """

}

process SV_VARIANT_ANNOTATE {
    publishDir "${params.outdir}", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    val db
    path target_bed

    output:
    tuple val(meta), path("*.sv.ann.std.vcf.gz"), emit: vcf
    tuple val(meta), path("*.ann.high_risk.csv"), emit: table

    script:
    """
    # convert BND to inversion
    convertInversion.py ${params.samtools_path} ${params.reference_path} $vcf | bcftools view - -o tmp.vcf.gz -O z

    # snpEff annotation
    snpeff -canon -noInteraction -motif $db -filterInterval $target_bed -stats ${meta}.summary.html tmp.vcf.gz > ${meta}.sv.ann.vcf

    # standardize vcf
    convertBND.py ${meta}.sv.ann.vcf | bcftools view - -o ${meta}.sv.ann.std.vcf.gz -O z

    # only keep the variants with PASS and high risk
    bcftools view ${meta}.sv.ann.std.vcf.gz -i 'FILTER="PASS" && ANN ~ "HIGH"' \\
        | bcftools query -f '%ID\\t%CHROM\\t%POS\\t%CHR2\\t%POS2\\t%QUAL\\t%SVLEN\\t%CHREND\\t%FILTER\\t[%GT]\\t%ANN\\n' \\
        | awk '{
            OFS="\\t";
            split(\$11,ann,",");
            annotation="";
            for(i in ann)
            {
                split(ann[i],f,"|");
                if(f[3]=="HIGH")
                {
                    annotation=annotation"\\t"f[4]"\\t"f[1]"|"f[2]"|"f[3]
                }
            }
            print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10annotation
            }' > tmp.high_risk.csv

    # add rmsk annotation
    awk -F "\\t" '{OFS="\\t"; print \$2,\$3-1,\$3}' tmp.high_risk.csv | bedtools intersect -a stdin -b ${params.rmsk_path} -wao | cut -f7 > tmp1
    awk -F "\\t" '{OFS="\\t"; print \$4,\$5-1,\$5}' tmp.high_risk.csv | bedtools intersect -a stdin -b ${params.rmsk_path} -wao | cut -f7 > tmp2
    paste tmp.high_risk.csv tmp1 tmp2 > ${meta}.sv.ann.high_risk.csv
    """
}


workflow {

    // load SNV vcf
    snv_vcf = Channel.fromPath("${params.dragen_path}/*.hard-filtered.vcf.{gz,gz.tbi}")
        .collect()
        .map { if (it[0].extension == "tbi")
                { [it[1].name.minus(".hard-filtered.vcf.gz"), it[1], it[0]] }
               else
                { [it[0].name.minus(".hard-filtered.vcf.gz"), it[0], it[1]] }
             }

    // load CNV vcf
    cnv_vcf = Channel.fromPath("${params.dragen_path}/*.cnv.vcf.{gz,gz.tbi}")
        .collect()
        .map { if (it[0].extension == "tbi")
                { [it[1].name.minus(".cnv.vcf.gz"), it[1], it[0]] }
               else
                { [it[0].name.minus(".cnv.vcf.gz"), it[0], it[1]] }
             }

    // load SV vcf
    sv_vcf = Channel.fromPath("${params.dragen_path}/*.sv.vcf.{gz,gz.tbi}")
        .collect()
        .map { if (it[0].extension == "tbi")
                { [it[1].name.minus(".sv.vcf.gz"), it[1], it[0]] }
               else
                { [it[0].name.minus(".sv.vcf.gz"), it[0], it[1]] }
             }

    // small variants annotation
    SNV_VARIANT_ANNOTATE(snv_vcf, params.snpeff_db, params.target_bed)
    // CNV variants annotation
    CNV_VARIANT_ANNOTATE(cnv_vcf, params.snpeff_db, params.target_bed)
    // SV variants annotation
    SV_VARIANT_ANNOTATE(sv_vcf, params.snpeff_db, params.target_bed)

}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}