nextflow run main.nf --fasta /expanse/projects/sebat1/vep_resources/SFARI/iWGS_v1/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa --ped /expanse/lustre/projects/ddp195/eiovino/ita_cohort/SNVs/ita_fam_nov_sorted.ped  --vcfs test_folder --vep_config vep.ini -entry RUN_VEP -profile slurm

nextflow run main.nf --fasta /expanse/projects/sebat1/vep_resources/SFARI/iWGS_v1/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa --ped /expanse/lustre/projects/ddp195/eiovino/ita_cohort/SNVs/ita_fam_nov_sorted.ped  --vcfs test_folder --vep_config vep.ini -entry  -profile slurm
