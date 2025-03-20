include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_MPILEUP } from '../../modules/nf-core/samtools/mpileup/main'
include { VARSCAN } from '../../modules/local/varscan/main'
include { TABIX_BGZIP} from "../../modules/nf-core/tabix/bgzip/main"

workflow VARSCAN_WF {

take:
    input_check_bam_out // file: [ [ meta ], bam ] file, fasta
    input_check_bed_out
    fasta

    main:

    ch_versions = Channel.empty()

    SAMTOOLS_SORT (
        input_check_bam_out
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())


    //
    // MODULE: Run Samtools mpileup
    //

    ch_mpileup = SAMTOOLS_SORT.out.bam.join(input_check_bed_out)


    SAMTOOLS_MPILEUP (
        ch_mpileup, fasta
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP.out.versions.first())


    //
    // MODULE: Run TABIX_BGZIP para descomprimir el vcf
    //

    TABIX_BGZIP (
        SAMTOOLS_MPILEUP.out.mpileup
    )
    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions.first())

    


    //
    // MODULE: Run VarScan
    //


    VARSCAN (
        TABIX_BGZIP.out.output
    )

    varscan_out = VARSCAN.out.vcf_varscan
    ch_versions = ch_versions.mix(VARSCAN.out.versions.first())

    emit:

    varscan_out
    
    ch_versions

}