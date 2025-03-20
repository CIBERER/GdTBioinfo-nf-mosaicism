//
// Run GATK Mutect2 in tumor-only mode + BCFTOOLS_REHEADER & TABIX_TABIX
//

include { GATK4_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/gatk4/createsequencedictionary'
include { GATK4_MUTECT2 as MUTECT2 } from '../../modules/nf-core/gatk4/mutect2/main'
include { BCFTOOLS_REHEADER } from '../../modules/nf-core/bcftools/reheader/main'
include { TABIX_TABIX as TABIX_TABIX } from '../../modules/nf-core/tabix/tabix/main'

workflow BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK {
    take:
    ch_input
    ch_fasta
    ch_fai
    ch_germline_resource
    ch_germline_resource_tbi
    ch_interval_file

    main:
    ch_versions = Channel.empty()

    // Convert FASTA & FAI into Nextflow-compatible format
    ch_fasta.map{ create_fasta_input(it) }.set{ ch_fasta_meta }
    ch_fai.map{ create_fai_input(it) }.set{ ch_fai_meta }

    // Generate sequence dictionary for FASTA
    GATK4_CREATESEQUENCEDICTIONARY (ch_fasta_meta)
    ch_dict_meta = GATK4_CREATESEQUENCEDICTIONARY.out.dict
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)

    // Run Mutect2
    ch_input.combine(ch_interval_file, by: 0).set{ ch_input_mutect2 }

    MUTECT2 (
        ch_input_mutect2,
        ch_fasta_meta,
        ch_fai_meta,
        ch_dict_meta,
        ch_germline_resource,
        ch_germline_resource_tbi
    )

    ch_versions = ch_versions.mix(MUTECT2.out.versions)

    // Prepare input for BCFTOOLS_REHEADER (empty header, sample "Mutect")
    MUTECT2.out.vcf
        .map { meta, vcf -> tuple(meta, vcf, [], []) }  // Empty header, fixed sample name
        .set { ch_reheader_input }

    // Run BCFTOOLS_REHEADER (renaming sample)
    BCFTOOLS_REHEADER (ch_reheader_input, tuple([], []))

    ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions)

    // Run TABIX on the renamed VCF
    TABIX_TABIX (BCFTOOLS_REHEADER.out.vcf)

    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    // Define output channels
    ch_vcf = BCFTOOLS_REHEADER.out.vcf
    ch_tabix = TABIX_TABIX.out.tbi

    emit:
    mutect2_vcf    = ch_vcf   // Store ONLY the final VCF
    mutect2_index  = ch_tabix // Store ONLY the final .tbi index
    versions       = ch_versions  // Version tracking
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUPPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Convert FASTA to metadata format
def create_fasta_input (fasta) {
    def meta = [:]
    meta.id = fasta.baseName
    if (!file(fasta).exists()) exit 1, "ERROR: FASTA file does not exist! -> ${fasta}"
    return [ meta, file(fasta) ]
}

// Convert FASTA index to metadata format
def create_fai_input (fai) {
    def meta = [:]
    meta.id = fai.baseName
    if (!file(fai).exists()) exit 1, "ERROR: FASTA index file does not exist! -> ${fai}"
    return [ meta, file(fai) ]
}
