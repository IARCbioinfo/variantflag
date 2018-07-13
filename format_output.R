stk_file <- paste0( data_path, "raw_VCF_strelka/", l_tag[1], ".bam_vs_E210_S6.bam.somatic.snvs.vcf.gz"  )
system( paste( "Rscript read_strelka.R", stk_file ) )

mtc_file <- paste0( data_path, "raw_VCF_mutect/", l_tag[1], "_calls.vcf"  )
system( paste( "Rscript read_mutect.R", mtc_file ) )
## not working ??

