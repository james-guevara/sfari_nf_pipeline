/* 
 * Workflow to run VEP on VCF files
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 */

include { VEP } from 'modules/local/vep.nf'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_DROP_GENOTYPES } from 'modules/nf-core/bcftools/view/main'
include { TABIX_TABIX as TABIX } from 'modules/nf-core/tabix/tabix/main'
include { PLINK2_VCF } from 'modules/nf-core/plink2/vcf/main'

if (params.help) {
  log.info ''
  log.info 'Pipeline to run VEP and PLINK2 on VCFs'
  log.info '-------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info '  nextflow run main.nf --ped <ped_filepath> --fasta <fasta_filepath> --vcfs <vcf_folder> --vep_config <vep_ini_filepath>'
  log.info ''
  log.info 'Options:'
  log.info '  --cpus INT                       Number of CPUs to use. Default 1.'
  log.info '  --ped PATH                       Path to PED file.'
  log.info '  --fasta PATH                     Path to reference FASTA file.'
  log.info '  --vcfs PATH                      Path to VCF files.'
  log.info '  --vep_config FILENAME            VEP config file. Default: vep.ini'
  exit 1
}

log.info 'Starting workflow.....'

if (!params.fasta) { exit 1, "Undefined --fasta parameter. Please provide the FASTA reference filepath." }
if (!params.ped == "") { exit 1, "Undefined --ped parameter. Please provide the PED filepath." }
if (!params.vcfs) { exit 1, "Undefined --vcfs parameter. Please provide the VCF folder." }
if (!params.vep_config) { exit 1, "Undefined --vep_config parameter. Please provide the VEP configuration filepath." }

Channel.fromPath(params.vcfs + "/" + "*.vcf.gz").map { it -> tuple([id: it.simpleName], it, file(it + ".tbi", checkIfExists: true)) }.set { vcf_channel }

workflow RUN_PLINK {
    PLINK2_VCF( vcf_channel )
}

workflow DROP_GENOTYPES {
    BCFTOOLS_VIEW_DROP_GENOTYPES( vcf_channel, [], [], [] )
    TABIX( BCFTOOLS_VIEW_DROP_GENOTYPES.out.vcf )
}

workflow RUN_VEP {
    VEP( vcf_channel, file(params.vep_config, checkIfExists: true), file(params.fasta, checkIfExists: true) )
}

workflow {
    vcf_channel.view()
}
