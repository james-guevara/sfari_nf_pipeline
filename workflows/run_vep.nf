/* 
 * Workflow to run VEP on VCF files
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 */

params.help = false
params.cpus = 1
params.vep_config = ""

include { VEP } from '../modules/local/vep.nf'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_DROP_GENOTYPES } from '../modules/nf-core/bcftools/view/main'
include { TABIX_TABIX as TABIX } from '../modules/nf-core/tabix/tabix/main'
include { PLINK2_VCF } from '../modules/nf-core/plink2/vcf/main'

if (params.help) {
  log.info ''
  log.info 'Pipeline to run VEP and PLINK2 on VCFs'
  log.info '-------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info '  nextflow -C nf_config/nextflow.config run workflows/run_vep.nf --vep_config <path_to_VEP_configuration_file>'
  log.info ''
  log.info 'Options:'
  log.info '  --vcf_folder PATH                Path to VCF files.'
  log.info '  --vep_config FILENAME            VEP config file. Default: nf_config/vep.ini'
  log.info '  --cpus INT	               Number of CPUs to use. Default 1.'
  exit 1
}

log.info 'Starting workflow.....'

workflow RUN_PLINK {
    if (!params.vcf_folder) { exit 1, "Undefined --vcf_folder parameter. Please provide the folder containing the VCFs." }
    vcf_folder = params.vcf_folder
    Channel.fromPath(vcf_folder + "/" + "*.vcf.gz").map { it -> tuple([id: it.simpleName], it }.set { vcf_channel }
    ped_file = file(params.ped_file, checkIfExists: true)
    PLINK2_VCF(vcf_channel)
}

workflow DROP_GENOTYPES {
    if (!params.vcf_folder) { exit 1, "Undefined --vcf_folder parameter. Please provide the folder containing the VCFs." }
    vcf_folder = params.vcf_folder
    Channel.fromPath(vcf_folder + "/" + "*.vcf.gz").map { it -> tuple([id: it.simpleName], it, file(it + ".tbi", checkIfExists: true)) }.set { vcf_channel }
    Channel.fromPath(vcf_folder + "/" + "*.vcf.gz").map { it -> tuple([id: it.simpleName], it, file(it + ".tbi", checkIfExists: true)) }.set { vcf_channel }
    BCFTOOLS_VIEW_DROP_GENOTYPES( vcf_channel, [], [], [] )
    TABIX( BCFTOOLS_VIEW_DROP_GENOTYPES.out.vcf )
}

workflow RUN_VEP {
    vep_config = file(params.vep_config, checkIfExists: true)
    fasta = file(params.fasta, checkIfExists: true)
    Channel.fromPath("/expanse/lustre/scratch/j3guevar/temp_project/nf_output_dir/*.vcf.gz").map { it -> tuple([id: it.simpleName], it, file(it + ".tbi", checkIfExists: true)) }.set { vcf_channel }
    VEP( vcf_channel, vep_config, fasta )
}

workflow {
  Channel.of( tuple( [id: "test"], vcf_file, vcf_index ) ).set { my_channel }
}
