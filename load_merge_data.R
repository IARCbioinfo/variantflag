###############################"
#### Read VCF

dt_stk <- fread( paste0( data_path, tag, "_pass_strelka.vcf" ), stringsAsFactors = F )[ grep( "chrUn|alt|random|M", X.CHROM, invert = T),  ]
dt_stk[ , ID := NULL ]
setnames( dt_stk, c( "AF_N", "AF_T" ), c( "VAF_N", "VAF_T" ) )
dt_stk[ , is_strelka := 1 ]

dt_mtc <- fread( paste0( data_path, tag, "_pass_mutect.vcf" ),
                             stringsAsFactors = F )[ grep( "chrUn|alt|random|M", X.CHROM, invert = T ), ]
dt_mtc[ , ID := NULL ]
setnames( dt_mtc, c( "AF_N", "AF_T" ), c( "VAF_N", "VAF_T" ) )
dt_mtc <- dt_mtc[ nchar( ALT ) == 1 & nchar( REF ) == 1 ]
dt_mtc[ , is_mutect := 1 ]

dt_mut_all <-  merge( dt_mtc, dt_stk,
                      by = c( "X.CHROM", "POS", "REF", "ALT" ),
                      suffix = c( ".M", ".S" ),
                      all = T )

dt_mut_all[ , nb_caller := apply( .SD[, .( is_mutect, is_strelka ) ], 1, sum, na.rm = T ) ]


################################
####### Get annotations
dt_annot <- fread( paste0( data_path, tag, "_allMutations_annot.txt" ), 
                               stringsAsFactors = F ) 

source('~/local_analyses/functions_mutations/functions_annotation.R')
source('~/local_analyses/functions_mutations/function_process.R')

add_mutation_type2( dt_annot )
dt_annot <- add_new_context( dt_annot )
add_new_strand( dt_annot )


annotate_dup( dt_annot )

dt_repeat <- setDT( read.table( paste0( data_path, tag, "_allPositions_intersectRepeat.bed" ), stringsAsFactors = F ) )[ , 1:4 ]
names( dt_repeat ) <- c( "Chr", "Start", "End", "Ovlap" )
dt_repeat <- merge( dt_repeat, dt_repeat[ Ovlap == "rpt", .N, by = c( "Chr", "Start" ) ][ ,.( Chr, Start, is_rpt = N ) ], 
                    by = c( "Chr", "Start" ), all.x = T )
dt_repeat <- merge( dt_repeat, dt_repeat[ Ovlap == "str", .N, by = c( "Chr", "Start" ) ][ ,.( Chr, Start, is_str = N ) ], 
                    by = c( "Chr", "Start" ), all.x = T )

dt_annot <- merge( dt_annot, unique( dt_repeat[ , .( Chr, Start, is_rpt, is_str ) ] ),
                   by = c( "Chr", "Start" ), all.x = T )


########## Multi allelic positions
dt_annot <- merge( dt_annot, dt_annot[ , .N, by = c( "Chr", "Start", "Ref" ) ][ , .( Chr, Start, Ref, N_Alt = N ) ],
                   by = c( "Chr", "Start", "Ref" ), all.x = T )


############## select annotations

# dt_annot[, .( Chr, Start, Ref, Alt, N_Alt,
#               mut_class, trinucleotide_context, new_context, 
#               Func.knownGene, ExonicFunc.knownGene, Strand, new_tx_status,
#               snp142, is_dup, is_rpt, is_str ) ]


dt_mut_annot <- merge( dt_annot[, .( Chr, Start, Ref, Alt, N_Alt,
                                     mut_class, trinucleotide_context, new_context, 
                                     Func.knownGene, ExonicFunc.knownGene, Strand, new_tx_status,
                                     snp142, is_dup, is_rpt, is_str ) ],
                       dt_mut_all,
                       by.x = c( "Chr", "Start", "Ref", "Alt" ), by.y = c( "X.CHROM", "POS", "REF", "ALT" ), all.y = T )



######## Set global Cov and VAF
dt_mut_annot[ nb_caller == 2, Cov_N := apply( .SD[ , .( Cov_N.M, Cov_N.S ) ], 1, max ) ]
dt_mut_annot[ nb_caller == 1 & is_strelka, Cov_N := Cov_N.S ]
dt_mut_annot[ nb_caller == 1 & is_mutect, Cov_N := Cov_N.M ]

dt_mut_annot[ nb_caller == 2, Cov_T := apply( .SD[ , .( Cov_T.M, Cov_T.S ) ], 1, max ) ]
dt_mut_annot[ nb_caller == 1 & is_strelka, Cov_T := Cov_T.S ]
dt_mut_annot[ nb_caller == 1 & is_mutect, Cov_T := Cov_T.M ]

dt_mut_annot[ nb_caller == 2, VAF_N := apply( .SD[ , .( VAF_N.M, VAF_N.S ) ], 1, min ) ]
dt_mut_annot[ nb_caller == 1 & is_strelka, VAF_N := VAF_N.S ]
dt_mut_annot[ nb_caller == 1 & is_mutect, VAF_N := VAF_N.M ]

dt_mut_annot[ nb_caller == 2, VAF_T := apply( .SD[ , .( VAF_T.M, VAF_T.S ) ], 1, min ) ]
dt_mut_annot[ nb_caller == 1 & is_strelka, VAF_T := VAF_T.S ]
dt_mut_annot[ nb_caller == 1 & is_mutect, VAF_T := VAF_T.M ]


