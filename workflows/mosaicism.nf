/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMosaicism.initialise(params, log)

// Validate paths
def checkPathParamList = [params.input, params.multiqc_config, params.fasta, params.fai]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters

if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Fasta file not specified!' }
if (params.fai) { ch_fasta_fai = file(params.fai) } else { exit 1, 'Fai file not specified!' }
if (params.chrom_sizes) { ch_chrom_sizes = file(params.chrom_sizes) } else { exit 1, 'chrom.sizes file not specified!' }
if (params.germline_resource) { ch_germline_resource = file(params.germline_resource) } else { exit 1, 'germline_resource VCF file not specified!' }
if (params.germline_resource_tbi) { ch_germline_resource_tbi = file(params.germline_resource_tbi) } else { exit 1, 'germline_resource TBI file not specified!' }


// Parse tools selected
def tool_list = params.tool ? params.tool.split(",") as List : []
def selected_tools = tool_list.findAll { it in ["vardict", "varscan", "mutect"] }


// Define empty channels first
ch_vardictjava = Channel.empty()
ch_varscan = Channel.empty()
ch_mutect2 = Channel.empty()



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { VARSCAN_WF } from '../subworkflows/local/varscan_workflow'
include { PROCESSING_VARSCAN } from '../subworkflows/local/processing_varscan'
include { BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK } from '../subworkflows/local/bam_tumor_only_somatic_variant_calling_gatk'
include { MERGE_WORKFLOW } from '../subworkflows/local/merge_workflow'

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { VARDICTJAVA } from '../modules/nf-core/vardictjava/main'
include { TABIX_BGZIP as TABIX_BGZIP_VARDICTJAVA } from "../modules/nf-core/tabix/bgzip/main"
include { TABIX_BGZIP as TABIX_BGZIP_JOIN } from "../modules/nf-core/tabix/bgzip/main"
include { TABIX_TABIX as TABIX_TABIX } from '../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_MERGE } from '../modules/nf-core/bcftools/merge/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MOSAICISM {
    ch_versions = Channel.empty()

    // Run input validation
    INPUT_CHECK ( ch_input )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Run Vardict
    if (selected_tools.contains("vardict")) {
        VARDICTJAVA (
            INPUT_CHECK.out.reads_bam_bai_bed, 
            tuple([], ch_fasta), 
            tuple([], ch_fasta_fai)
        )
        ch_versions = ch_versions.mix(VARDICTJAVA.out.versions.first())

        TABIX_BGZIP_VARDICTJAVA ( VARDICTJAVA.out.vcf )
        ch_versions = ch_versions.mix(TABIX_BGZIP_VARDICTJAVA.out.versions.first())

        TABIX_TABIX ( TABIX_BGZIP_VARDICTJAVA.out.output )
        ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())
        ch_vardictjava = TABIX_BGZIP_VARDICTJAVA.out.output.join(TABIX_TABIX.out.tbi)

    }

    // Run Varscan
    if (selected_tools.contains("varscan")) {
        VARSCAN_WF ( INPUT_CHECK.out.reads_bam, INPUT_CHECK.out.reads_bed, ch_fasta )
        ch_versions = ch_versions.mix(VARSCAN_WF.out.ch_versions)

        PROCESSING_VARSCAN ( VARSCAN_WF.out.varscan_out )
        ch_versions = ch_versions.mix(PROCESSING_VARSCAN.out.ch_versions)

        ch_varscan = PROCESSING_VARSCAN.out.annotation_out.join(PROCESSING_VARSCAN.out.annotation_tabix)
    }

    // Run Mutect
    if (selected_tools.contains("mutect")) {
        BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK(
            INPUT_CHECK.out.reads_bam_bai, 
            Channel.fromPath(ch_fasta).first(), 
            Channel.fromPath(ch_fasta_fai).first(), 
            Channel.fromPath(ch_germline_resource).first(), 
            Channel.fromPath(ch_germline_resource_tbi).first(), 
            INPUT_CHECK.out.reads_bed
        )
        ch_versions = ch_versions.mix(BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK.out.versions.first())

        ch_mutect2 = BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK.out.mutect2_vcf.join(BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK.out.mutect2_index)
    }

    
// Merge channels using `.groupTuple()`

ch_merged_input = ch_vardictjava.mix(ch_varscan).mix(ch_mutect2).groupTuple()

// Print for debugging


// Run MERGE_WORKFLOW if multiple tools were used
if (selected_tools.size() > 1) {
    MERGE_WORKFLOW(ch_merged_input)
    ch_versions = ch_versions.mix(MERGE_WORKFLOW.out.versions)
    
    TABIX_BGZIP_JOIN ( MERGE_WORKFLOW.out.merge_vcf )
    ch_versions = ch_versions.mix(TABIX_BGZIP_JOIN.out.versions.first())
}

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
