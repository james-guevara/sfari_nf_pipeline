/* 
 * Workflow to run VEP on VCF files
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 */

params.help = false
params.cpus = 1
params.outdir = "outdir"
params.singularity_dir = ""
params.vep_config = ""
params.chroms = ""
params.chroms_file = ""

include { VEP } from '../modules/vep.nf'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_DROP_GENOTYPES } from '../../modules/nf-core/bcftools/view/main'
include { TABIX_TABIX as TABIX } from '../../modules/nf-core/tabix/tabix/main'
include { PLINK2_VCF } from '../../modules/nf-core/plink2/vcf/main'

if (params.help) {
  log.info ''
  log.info 'Pipeline to run VEP on VCFs'
  log.info '-------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info '  nextflow -C nf_config/nextflow.config run workflows/run_vep.nf --vep_config <path_to_VEP_configuration_file>'
  log.info ''
  log.info 'Options:'
  log.info '  --vep_config FILENAME            VEP config file. Default: nf_config/vep.ini'
  log.info '  --cpus INT	               Number of CPUs to use. Default 1.'
  exit 1
}

if (!params.vcf) { exit 1, "Undefined --vcf parameter. Please provide the path to a VCF file" }
vcf_file = file(params.vcf, checkIfExists: true)
vcf_index = file("${params.vcf}.tbi", checkIfExists: true)

ped_file = file(params.ped_file, checkIfExists: true)

if (!params.vep_config) { exit 1, "Undefined --vep_config parameter. Please provide a VEP config file" }
vep_config = file(params.vep_config, checkIfExists: true)

log.info 'Starting workflow.....'

workflow DROP_GENOTYPES {
  // vcf_channel = Channel.fromPath("/expanse/lustre/scratch/j3guevar/temp_project/bcf_test/*.vcf.gz")
  // vcf_channel.map { it -> tuple([id: it.simpleName], it, file(it + ".tbi", checkIfExists: true)) }.set { my_channel2 }

  Channel.fromPath("/expanse/lustre/scratch/j3guevar/temp_project/bcf_test/*.vcf.gz").map { it -> tuple([id: it.simpleName], it, file(it + ".tbi", checkIfExists: true)) }.set { vcf_channel }
  BCFTOOLS_VIEW_DROP_GENOTYPES( vcf_channel, [], [], [] )
  TABIX( BCFTOOLS_VIEW_DROP_GENOTYPES.out.vcf )

  // my_channel2.view()
  // Channel.of( tuple( [id: "test"], vcf_file, vcf_index ) ).set { my_channel }
  // TABIX_TABIX( BCFTOOLS_VIEW.out.vcf )
  // vcf_ch = BCFTOOLS_VIEW.out.vcf.join(TABIX_TABIX.out.tbi)
  // BCFTOOLS_VIEW.out.vcf.view()
}

workflow RUN_VEP {
  fasta = file(params.fasta, checkIfExists: true)
  Channel.fromPath("/expanse/lustre/scratch/j3guevar/temp_project/nf_output_dir/*.vcf.gz").map { it -> tuple([id: it.simpleName], it, file(it + ".tbi", checkIfExists: true)) }.set { vcf_channel }
  VEP( vcf_channel, vep_config, fasta )

}

workflow RUN_PLINK {

}

workflow {
  Channel.of( tuple( [id: "test"], vcf_file, vcf_index ) ).set { my_channel }
  // my_channel.view()
  // BCFTOOLS_VIEW( my_channel, [], [], [] )
  // VEP(vcf_file, vcf_index, vep_config)
  PLINK2_VCF( my_channel )
}
