/* 
 * Script to run VEP on VCF files
 */

process VEP {
  /*
  Function to run VEP on VCF files

  Returns
  -------
  Returns 2 files per chromosome:
      1) VEP output file for each VCF
      2) A tabix index for that VCF output file
  */
  container "${params.singularity_dir}/vep.sif"

  input:
  tuple val(meta), path(vcf), path(index)
  path(vep_config)
  path(fasta)
  
  output:
  tuple val(meta), path("*.vep.vcf.gz"), path("*.vep.vcf.gz.tbi") ,  emit: vcf_with_index
  path("*.vep.vcf.gz_summary.*") ,                                   emit: report

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def fasta_file = fasta ? "--fasta ${fasta}" : ""
  """
  vep --input_file ${vcf} \\
      --format vcf \\
      --output_file ${prefix}.vep.vcf.gz \\
      --vcf \\
      --compress_output bgzip \\
      --config ${vep_config} \\
      ${fasta_file} 
  tabix ${prefix}.vep.vcf.gz
  """
}
