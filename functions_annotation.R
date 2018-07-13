annotate_SNP <- function( tab ){
  tab[ , is_SNP := 0 ]
  tab[ ( X1000g2015aug_all >= 0.001 
         | ExAC_nontcga_ALL  >= 0.001 
         | gnomAD_genome_ALL >= 0.001 
         | esp6500siv2_all >= 0.00001 ), 
       is_SNP := 1 ]
}

annotate_cosmic <- function( tab ){
  tab[ , is_cosmic := 0 ]
  tab[ ( ! is.na( cosmic83 ) ), 
       is_cosmic := 1 ]
}

annotate_dup <- function( tab ){
  tab[ , is_dup := 0 ]
  tab[ ( ! is.na( genomicSuperDups ) ), 
       is_dup := 1 ]
}

