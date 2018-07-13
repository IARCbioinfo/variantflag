library( data.table )
library( ggplot2 )


############
### variantFlag implementation with MNU data

data_path <- "~/data/DATA-1/MMB_MEF_WGS_signatureRatio/" 
genome <- "mm10"
annot_path <- paste0( "~/data/DATA-1/annotations/", genome, "/" )


###### 0. filter out PASS mutation + format 
l_tag <- c( "MNU-1", "MNU-2", "0-5mM-OTA-S9-1_S7", "0-5mM-OTA-S9-2_S8", "MK01-OTA-Spont-S9-1_S1", "MK02-OTA-Spont-S9-2_S2" )
l_new_tag <- c( "MNU-1", "MNU-2", "OTA-1", "OTA-2", "Spont_OTA-1", "Spont_OTA-2" )
l_norm <- c( "E210_S6", "E210_S6", "E210_S6", "E210_S6", "MK06-E210", "MK06-E210" )

for (i in 1:6 ){
  tag <- l_tag[i]
  tag_norm <- l_norm[i]
  new_tag <- l_new_tag[i]
  
  stk_name <- paste0( tag, ".bam_vs_", tag_norm, ".bam.somatic.snvs.vcf"  )
  #nf <- paste0( data_path, "raw_VCF_strelka/",  tag, ".bam_vs_", tag_norm, ".bam.somatic.snvs.vcf.gz"  )
  nf <- paste0( data_path, "raw_VCF_strelka/",  stk_name, ".gz"  )
  system( paste( "Rscript read_strelka.R", nf ) )
  system( paste( "mv", gsub( ".gz", "", nf ), paste0( data_path, new_tag, "_pass_strelka.vcf" ) ) )
  
  nf <- paste0( data_path, "raw_VCF_mutect/",  tag, "_calls.vcf"  )
  mtc_name <- paste0( tag, "_calls.vcf" )
  system( paste( "Rscript read_mutect.R", nf ) )
  system( paste( "mv", paste0( data_path, "raw_VCF_mutect/", tag, "_pass_mutect.vcf" ), paste0( data_path, new_tag, "_pass_mutect.vcf" ) ) )
}


l_tag <- c( "MNU-1", "MNU-2", "0-5mM-OTA-S9-1_S7", "0-5mM-OTA-S9-2_S8", "MK01-OTA-Spont-S9-1_S1", "MK02-OTA-Spont-S9-2_S2" )
l_new_tag <- c( "MNU-1", "MNU-2", "OTA-1", "OTA-2", "Spont_OTA-1", "Spont_OTA-2" )
l_norm <- c( "E210_S6", "E210_S6", "E210_S6", "E210_S6", "MK06-E210", "MK06-E210" )




###### 1. Get mutations and intersect / annotate with DB

#source( "variant_parsing.R" )
l_tag <- c( "MNU-1", "MNU-2", "OTA-1", "OTA-2", "Spont_OTA-1", "Spont_OTA-2" )
system( "bash get_header.sh" )
for ( i in c( 1:4, 6 ) ){
  system( paste( "bash parse_file.sh", data_path, l_tag[i], annot_path, genome ) )
}
## run annovar

###### 2. Parse outputs and merge all infos
dt_all_mut <- NULL

l_tag <- c( "MNU-1", "MNU-2", "OTA-1", "OTA-2", "Spont_OTA-1", "Spont_OTA-2" )

for ( i in c( 1:4, 6 ) ){
  source( "load_merge_data.R" )
  ### Usual filters
  th_Cov <- 8
  dt_mut_annot[ , keep_mutation := 0 ]
  dt_mut_annot[ is.na( snp142 ) & is.na( is_str ) & VAF_N == 0.0 & Cov_T >= th_Cov & Cov_N >= th_Cov, keep_mutation := 1 ]
  
  dt_all_mut <- rbind( dt_all_mut, dt_mut_annot )
  save( dt_mut_annot, file = paste0( tag, "_dt_mut_annot.Rdata" ) )
  
  rmarkdown::render( "report_output.Rmd", params = list( tag = "MNU-1" ), output_file = paste0( "report_output_", tag, ".html" ) )
}



###### 3. Process mutations and propose filtering flags




