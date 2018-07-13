###### Define mutation and context

add_mutation_type <- function( tab ){
  tab[ , mut_class := "C>A" ]
  tab[ ( Ref == "C" & Alt == "G" ) | ( Ref == "G" & Alt == "C" ), mut_class := "C>G" ]
  tab[ ( Ref == "C" & Alt == "T" ) | ( Ref == "G" & Alt == "A" ), mut_class := "C>T" ]
  tab[ ( Ref == "T" & Alt == "A" ) | ( Ref == "A" & Alt == "T" ), mut_class := "T>A" ]
  tab[ ( Ref == "T" & Alt == "C" ) | ( Ref == "A" & Alt == "G" ), mut_class := "T>C" ]
  tab[ ( Ref == "T" & Alt == "G" ) | ( Ref == "A" & Alt == "C" ), mut_class := "T>G" ]
  return( tab )
}

add_mutation_type2 <- function( tab ){
  tab[ , mut_class := "-" ]
  tab[ ( Ref == "C" & Alt == "A" ) | ( Ref == "G" & Alt == "T" ), mut_class := "C>A" ]
  tab[ ( Ref == "C" & Alt == "G" ) | ( Ref == "G" & Alt == "C" ), mut_class := "C>G" ]
  tab[ ( Ref == "C" & Alt == "T" ) | ( Ref == "G" & Alt == "A" ), mut_class := "C>T" ]
  tab[ ( Ref == "T" & Alt == "A" ) | ( Ref == "A" & Alt == "T" ), mut_class := "T>A" ]
  tab[ ( Ref == "T" & Alt == "C" ) | ( Ref == "A" & Alt == "G" ), mut_class := "T>C" ]
  tab[ ( Ref == "T" & Alt == "G" ) | ( Ref == "A" & Alt == "C" ), mut_class := "T>G" ]
  tab[ , change_strand := 0 ]
  tab[ Ref %in% c( "A", "G" ), change_strand := 1 ]
  
  return( tab )
}



get_complementary <- function( tab, inCol, outCol ){
  tab[ , c( outCol ) := "N" ]
  tab[ get( inCol ) == "A", c( outCol ) := "T" ]
  tab[ get( inCol ) == "C", c( outCol ) := "G" ]
  tab[ get( inCol ) == "G", c( outCol ) := "C" ]
  tab[ get( inCol ) == "T", c( outCol ) := "A" ]
  return( tab )
}

# get_complementary_3nt <- function( v ){
#   return( paste0( 
#     get_complementary( substr( v, 1, 1 ) ),
#     substr( v, 2, 2 ),
#     get_complementary( substr( v, 3, 3 ) )
#   ) )
# }

add_new_context <- function( tab ){
  tab[ , c( "L_ctxt", "R_ctxt" ) := tstrsplit( trinucleotide_context, "", keep = c( 1,3 ) ) ]
  get_complementary( tab, "L_ctxt", "L_comp" )
  get_complementary( tab, "R_ctxt", "R_comp" )
  
  tab[ , new_context := trinucleotide_context ]
  tab[ change_strand == 1, new_context := paste0( R_comp, "x", L_comp ) ]
  
  tab[ , c( "L_ctxt", "R_ctxt", "L_comp", "R_comp", "change_strand" ) := NULL ]
}

add_new_context2 <- function( tab ){
  tab[ , c( "L_ctxt", "R_ctxt" ) := tstrsplit( trinucleotide_context, "", keep = c( 1,3 ) ) ]
  get_complementary( tab, "L_ctxt", "L_comp" )
  get_complementary( tab, "R_ctxt", "R_comp" )
  
  tab[ , new_context := trinucleotide_context ]
  tab[ change_strand == 1, new_context := paste0( R_comp, "x", L_comp ) ]
  
  tab[ , new_tx_status := tx_status ]
  tab[ change_strand == 1 & Strand == "+", new_tx_status := "Untranscribed" ]
  tab[ change_strand == 1 & Strand == "-", new_tx_status := "Transcribed" ]
  
  tab[ , c( "L_ctxt", "R_ctxt", "L_comp", "R_comp", "change_strand" ) := NULL ]
}


add_new_strand <- function( tab ){
  tab[ Strand == "+", tx_status := "Transcribed" ]
  tab[ Strand == "-", tx_status := "Untranscribed" ]
  tab[ , new_tx_status := tx_status ]
  tab[ Ref %in% c( "A", "G" ) & Strand == "+", new_tx_status := "Unstranscribed" ]
  tab[ Ref %in% c( "A", "G" ) & Strand == "-", new_tx_status := "Transcribed" ]
}




############ count matrix

totMut <- data.table( mut_class = c( rep( "C>A", 16 ), rep( "C>G", 16 ), rep( "C>T", 16 ), 
                                 rep( "T>A", 16 ), rep( "T>C", 16 ), rep( "T>G", 16 ) ), 
                  context = c( "AxA", "AxC", "AxG", "AxT", "CxA", "CxC", "CxG", "CxT", 
                               "GxA", "GxC", "GxG", "GxT", "TxA", "TxC", "TxG", "TxT" ) )


get_tab96 <- function( tab, ctxt = "new_context" ){
  t <- tab[ grep( "N", get( ctxt ), invert = T ), .N, by = c( "mut_class", ctxt ) ]
  setnames( t, c( ctxt, "N" ), c( "context", "Count" ) ) 
  t <- merge( totMut, t, by = c( "mut_class", "context" ), all.x = T)
  t[ is.na( Count ), Count := 0 ]
  t[ , Percent := Count / nrow( tab ) ]
  return( t[ order( mut_class, context ) ] )
}


get_tab96_v1 <- function( tab, ctxt = "new_context", compare.with = NULL, correct.count = NULL ){
  t <- tab[ grep( "N", get( ctxt ), invert = T ), .N, by = c( "mut_class", ctxt ) ]
  #setnames( t, c( ctxt, "N", "frac" ), c( "context", "Count", "Percent" ) ) 
  setnames( t, c( ctxt, "N" ), c( "context", "Count" ) ) 
  if( ! is.null( compare.with ) ){
    t <- merge( compare.with[ , .( mut_class, context ) ], t, 
                by = c( "mut_class", "context" ),
                all.x = T)
    t[ is.na( Count ), Count := 0 ]
  }
  t[ , Percent := Count / nrow( tab ) ]
  return( t[ order( mut_class, context ) ] )
}



### normalization
get_tab96_opportunity <- function( tab, ctxt = "new_context", norm_tab ){
  t <- tab[ grep( "N", get( ctxt ), invert = T ), .N, by = c( "mut_class", ctxt ) ]
  setnames( t, c( ctxt, "N" ), c( "context", "Count" ) ) 
  t <- merge( totMut, t, by = c( "mut_class", "context" ), all.x = T)
  t[ is.na( Count ), Count := 0 ]
  ####
  t[ , Percent := Count / nrow( tab ) ]
  return( t[ order( mut_class, context ) ] )
}




######### Specific signatures

get_count_sign17 <- function( tab ){
  return( tab[ mut_class == "T>A" & context %in% c( "CxC", "CxG", "CxT" )
               | mut_class == "T>C" & context %in% c( "CxC", "CxG", "CxT" )
               | mut_class == "T>G" & context %in% c( "AxT", "CxC", "CxG", "CxT", "GxT", "TxT" ), 
               sum( Count ) ] )
}

flag_mut_sig17 <- function( tab ){
  tab[ , is_mut_sig17 := 0 ]
  tab[( ( mut_class == "T>A" & new_context %in% c( "CxC", "CxG", "CxT" ) ) 
        | ( mut_class == "T>C" & new_context %in% c( "CxC", "CxG", "CxT" ) )
        | ( mut_class == "T>G" & new_context %in% c( "AxT", "CxC", "CxG", "CxT", "GxT", "TxT" ) ) ), 
      is_mut_sig17 := 1 ]
}


flag_mut_sig18 <- function( tab ){
  tab[ , is_mut_sig18 := 0 ]
  tab[( mut_class == "C>A"  ), is_mut_sig18 := 1 ]
}


flag_mut_sigCpG <- function( tab ){
  tab[ , is_mut_sigCpG := 0 ]
  tab[( mut_class == "C>T" & new_context %in% c( "AxG", "CxG", "GxG", "TxG" ) ), is_mut_sigCpG := 1 ]
}


