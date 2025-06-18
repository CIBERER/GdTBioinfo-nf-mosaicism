//
// Check input samplesheet and get read, sample, and case channels
//

include { samplesheetToList } from 'plugin/nf-schema'

workflow INPUT_CHECK {
    take:
        samplesheet // file: /path/to/samplesheet.csv

    main:

    canal_sp = Channel.fromList(samplesheetToList(samplesheet, "assets/schema_input.json"))

    canal_sp
        .map {create_bam_bai_bed_channel(it)}
        .set { reads_bam_bai_bed }

    canal_sp
        .map {create_bam_channel(it)}
        .set { reads_bam }
    
    canal_sp
        .map {create_bam_bai_channel(it)}
        .set { reads_bam_bai }

    canal_sp
        .map {create_bed_channel(it)}
        .set { reads_bed }

    //reads_bam = reads.flatten().first().concat(reads.flatten().filter(~/.*.bam/).toList()).toList()

    emit:
    reads_bam_bai_bed // channel: [ meta, [ bam, bai, bed] ]
    reads_bam // channel: [ meta, [ bam ] ]
    reads_bam_bai // channel: [ meta, [ bam, bai ] ]
    reads_bed // channel: [ meta, [ bed ] ]
    versions = Channel.empty() // SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]

}

//Function to get list of [ meta, [ bam, bai, bed] ]

def create_bam_bai_bed_channel(ArrayList row) {
    // gather meta
    meta = row.get(0)

    // add path(s) of the bam/bai/bed files to the meta map
    bam_meta = [ meta, file(row.get(1)), file(row.get(2)), file(row.get(3)) ]

    return bam_meta
}

def create_bam_channel(ArrayList row) {
    // gather meta
    meta = row.get(0)

    // add path(s) of the bam files to the meta map
    bam_meta = [ meta, file(row.get(1)) ]

    return bam_meta
}

def create_bam_bai_channel(ArrayList row) {
    // gather meta
    meta = row.get(0)

    // add path(s) of the bam/bai files to the meta map
    bam_bai_meta = [ meta, file(row.get(1)), file(row.get(2)) ]

    return bam_bai_meta
}

def create_bed_channel(ArrayList row) {
    // gather meta
    meta = row.get(0)

    // add path(s) of the bed files to the meta map
    bed_meta = [ meta, file(row.get(3)) ]

    return bed_meta
}
