include { BCFTOOLS_MERGE } from '../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_REHEADER } from '../../modules/nf-core/bcftools/reheader/main'
include { RENAME } from '../../modules/local/rename/main'

workflow MERGE_WORKFLOW {
    take:
    ch_merged_input  // Receives grouped input

    main:
    ch_versions = Channel.empty()

    // Run merge
    BCFTOOLS_MERGE(
        ch_merged_input,
        tuple([], []), 
        tuple([], []),
        []
    )

    ch_merge = BCFTOOLS_MERGE.out.merged_variants
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

    
    emit:
    merge_vcf    = ch_merge
    versions = ch_versions
}
