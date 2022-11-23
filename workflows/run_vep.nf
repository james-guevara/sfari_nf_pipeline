/* 
 * Workflow to run VEP on VCF files
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 */

params.vep_config = "vep.ini"

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
  log.info '  nextflow run workflows/run_vep.nf --ped <ped_filepath> --fasta <fasta_filepath> --vcfs <vcf_folder> --vep_config <vep_ini_filepath>'
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

Channel.fromPath(params.vcfs + "/" + "*.vcf.gz").map { it -> tuple([id: it.simpleName, ped: file(params.ped, checkIfExists: true)], it, file(it + ".tbi", checkIfExists: true) ) }.set { vcf_channel }
// Channel.fromPath(params.vcfs + "/" + "*.vcf.gz").map { it -> tuple([id: it.simpleName, ped: params.ped], it, file(it + ".tbi", checkIfExists: true) ) }.set { vcf_channel }

workflow RUN_PLINK {
    PLINK2_VCF( vcf_channel )
}

workflow DROP_GENOTYPES {
    Channel.fromPath(params.vcfs + "/" + "*.vcf.gz").map { it -> tuple([id: it.simpleName], it, file(it + ".tbi", checkIfExists: true)) }.set { vcf_channel }
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
    print(params.fasta)
    print(params.ped)
    vcf_channel.view()
}


