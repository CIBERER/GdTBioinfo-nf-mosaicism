/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMosaicism.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta , params.fai]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters

if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Fasta file not specified!' }
if (params.fai) { ch_fasta_fai = file(params.fai) } else { exit 1, 'Fai file not specified!' }
if (params.chrom_sizes) { ch_chrom_sizes = file(params.chrom_sizes) } else { exit 1, 'chrom.sizes file not specified!' }
if (params.germline_resource) { ch_germline_resource = file(params.germline_resource) } else { exit 1, 'germline_resource VCF file not specified!' }
if (params.germline_resource_tbi) { ch_germline_resource_tbi = file(params.germline_resource_tbi) } else { exit 1, 'germline_resource TBI file not specified!' }
ch_known_sites = params.known_snps            ? Channel.fromPath(params.known_snps).collect()              : Channel.value([])
ch_known_sites_tbi = params.known_snps_tbi ? Channel.fromPath(params.known_snps_tbi) : Channel.empty()

//ch_intervals = params.intervals ? Channel.fromPath(params.intervals).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.value("")

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
//ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { MAPPING } from '../subworkflows/local/mapping'
include { VARSCAN_WF } from '../subworkflows/local/varscan_workflow'
include { PROCESSING_VARSCAN } from '../subworkflows/local/processing_varscan'
include { BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK } from '../subworkflows/local/bam_tumor_only_somatic_variant_calling_gatk'
include { MERGE_WORKFLOW } from '../subworkflows/local/merge_workflow'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { BWA_INDEX } from '../modules/nf-core/bwa/index/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from '../modules/nf-core/picard/createsequencedictionary/main'
include { VARDICTJAVA } from '../modules/nf-core/vardictjava/main'
include { TABIX_BGZIP } from "../modules/nf-core/tabix/bgzip/main"
include { TABIX_TABIX } from '../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_MERGE } from '../modules/nf-core/bcftools/merge/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MOSAICISM {

  ch_versions = Channel.empty()
 

  //
  // SUBWORKFLOW: Read in samplesheet, validate and stage input files
  //
  INPUT_CHECK (
      ch_input
  )
  ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

  if (params.index) { 
        ch_index = Channel.fromPath(params.index).map{ it -> [ [id:it.baseName], it ] }.collect() 
    } else {
    BWA_INDEX(
        ch_fasta
    )
    ch_index = BWA_INDEX.out.index
   }

   if (params.refdict) { 
        ch_refdict = Channel.fromPath(params.refdict).map{ it -> [ [id:it.baseName], it ] }.collect() 
    } else { 
    PICARD_CREATESEQUENCEDICTIONARY (
        tuple([], ch_fasta)
    )
    ch_refdict = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict
   }

    //
    // SUBWORKFLOW: Align raw reads
    //

    ch_intervals = params.intervals ? INPUT_CHECK.out.reads.map{ meta, fastqs -> tuple(meta, file(params.intervals)) } : INPUT_CHECK.out.reads.map{ meta, fastqs -> tuple(meta, []) }

    MAPPING (
        INPUT_CHECK.out.reads,
        ch_index,
        ch_fasta,
        ch_fai,
        ch_refdict,
        ch_intervals,
        ch_known_sites,
        ch_known_sites_tbi
    )

  //
  // SUBWORKFLOW: Run Samtools_sort, Samtools_mpileup and Varscan
  //

  VARSCAN_WF(

    INPUT_CHECK.out.reads_bam, INPUT_CHECK.out.reads_bed, ch_fasta

    )

ch_versions = ch_versions.mix(VARSCAN_WF.out.ch_versions)

  PROCESSING_VARSCAN(

  VARSCAN_WF.out.varscan_out
  )

  
  ch_versions = ch_versions.mix(PROCESSING_VARSCAN.out.ch_versions)

  
//
  // SUBWORKFLOW: Run Mutect2 and related software
  //
  BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK(
    INPUT_CHECK.out.reads_bam_bai, Channel.fromPath(ch_fasta).first(), Channel.fromPath(ch_fasta_fai).first(), Channel.fromPath(ch_germline_resource).first(), Channel.fromPath(ch_germline_resource_tbi).first(), INPUT_CHECK.out.reads_bed
  )
  
  
 



  ch_versions = ch_versions.mix(BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK.out.versions.first())




  //
  // MODULE: Run Vardictjava
  //




  VARDICTJAVA (
    INPUT_CHECK.out.reads_bam_bai_bed, tuple([],ch_fasta), tuple([],ch_fasta_fai)
  )



  ch_versions = ch_versions.mix(VARDICTJAVA.out.versions.first())

  //ch_versions.mix(VARDICTJAVA.out.versions.first().view())
  //ch_versions.mix(VARDICTJAVA.out.versions.view())

  TABIX_BGZIP (
        VARDICTJAVA.out.vcf
    )

   
   ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions.first())



   // MODULE: Run tabix

   TABIX_TABIX (
       TABIX_BGZIP.out.output
   )

   

ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())
  

ch_vardictjava=TABIX_BGZIP.out.output.join(TABIX_TABIX.out.tbi)
ch_varscan=PROCESSING_VARSCAN.out.annotation_out.join(PROCESSING_VARSCAN.out.annotation_tabix)
ch_mutect2=BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK.out.mutect2_vcf.join(BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK.out.mutect2_index)

  // Include BCFTOOLS_MERGE process



 MERGE_WORKFLOW (
 ch_vardictjava,
 ch_varscan,
ch_mutect2

)


ch_versions = ch_versions.mix(MERGE_WORKFLOW.out.versions)



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
