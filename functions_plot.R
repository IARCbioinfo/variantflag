v_cols <- c( "skyblue3", "black" , "red2", "grey", "green3", "pink" )
v_cols <- c( "skyblue2", "black" , "red2", "grey", "palegreen3", "pink" )



#############################
#### Plot 6 mutation types

plot_count6_v1 <- function( tab, tit = "" ){
  return( ggplot( tab ) 
          + geom_col( mapping = aes( x = mut_class, y = N, fill = mut_class ) ) 
          + scale_fill_manual( values = v_cols ) 
          + ylab( label = "Count" )
          #+ scale_y_continuous( name = "Count", size = 12 )  
          + xlab( label = "" )
          + ggtitle( label = tit )
          + theme_minimal()
          + theme( legend.position =  "none", 
                   #axis.line.y = element_line( colour = "black", size = 0.4 ), 
                   #axis.ticks.y = element_line( colour = "black" ), 
                   axis.text = element_text( size = 12 ),
                   axis.title.y = element_text( size = 15 ) ) 
          )
}


plot_count6_facets <- function( tab, fact, tit = "" ){
  return( ggplot( tab ) 
          + geom_col( mapping = aes( x = mut_class, y = N, fill = mut_class ) ) 
          + facet_wrap(  ~ as.factor( get( fact ) ) )
          + scale_fill_manual( values = v_cols ) 
          + ylab( label = "Count" )
          + xlab( label = "" )
          + ggtitle( label = tit )
          + theme_minimal()
          + theme( legend.position =  "none", 
                   axis.line.y = element_line( colour = "black", size = 0.4 ), 
                   axis.ticks.y = element_line( colour = "black" ), 
                   axis.text = element_text( size = 10 ) ) 
    )
  
}


plot_count6 <- function( tab, fact = NULL, tit = "" ){
  p <- ( ggplot( tab ) 
         + geom_col( mapping = aes( x = mut_class, y = N, fill = mut_class ) ) 
         + scale_fill_manual( values = v_cols ) 
         + ylab( label = "Count" )
         + xlab( label = "" )
         + ggtitle( label = tit )
         + theme_minimal()
         + theme( legend.position =  "none", 
                  axis.text = element_text( size = 12 ),
                  axis.title.y = element_text( size = 15 ) ) 
  )
  if( ! is.null( fact ) ){
    p <- p + facet_wrap(  ~ as.factor( get( fact ) ) ) + theme( strip.text = element_text( size = 14 ) )
  }
  return( p )
}


#########################
### Strand bias

plot_strand <- function( tab, nm_col = "N" ){
  return( ggplot( tab ) 
          #+ geom_col( mapping = aes( x = mut_class, y = get( nm_col ), fill = Strand ), position = "dodge" ) 
          + geom_col( mapping = aes( x = mut_class, y = get( nm_col ), fill = new_strand ), position = "dodge" ) 
          + scale_fill_manual( values = c( "+" = "blue3", "-" = "red3" ), labels = c( "-" = "Untranscribed", "+" = "Transcribed" ) )
          + theme(legend.position = "bottom")
          + xlab( label = "" )
          + ylab( label = "Count" )
  )
  
}

plot_strand2 <- function( tab, nm_col = "N" ){
  return( ggplot( tab ) 
          + geom_col( mapping = aes( x = mut_class, y = get( nm_col ), fill = new_tx_status ), position = "dodge" ) 
          + scale_fill_manual( values = c( "Transcribed" = "blue3", "Untranscribed" = "red3" ) )
          + theme(legend.position = "bottom")
          + xlab( label = "" )
          + ylab( label = "Count" )
  )
  
}


##########################
#### Mutation spectrum (96 mut, Count or percent)

plot_spectra_easy <- function( tab, val, tit = "", add_number = F ){
  ###val can be "Count" or "Percent"
  vCount <- tab[ , .N, by = mut_class ]$N
  mx <- max( tab[ , get( val ) ] )
  tsum <- tab[ , sum( Count ), by = mut_class ]
  if( add_number == T){ lbl <- tsum[ , paste( mut_class, V1, sep = "\n" ) ] 
  }else{ lbl <- tsum$mut_class }
  #names( lbl ) <- zsum$mut_class
  names( lbl ) <- tsum$mut_class
  
  p <- ( ggplot( tab, 
                 aes( x= context, y= get( val ), fill= mut_class ) ) 
         + geom_col( width=0.5) 
         #+ facet_grid( . ~ mut_class, scales="free_y")
         + facet_grid( . ~ mut_class, scales="free_y", labeller = labeller( mut_class = lbl  ) )
         + scale_fill_manual( values = v_cols ) 
         + ylab( label = val )
         + xlab( label = "" )
         + ggtitle( label = tit )
         + theme_minimal()
         + theme( legend.position =  "none", 
                  panel.grid = element_blank(), panel.grid.major.y = element_line( colour = "grey90" ),
                  #axis.line.y = element_line( colour = "black", size = 0.5 ), 
                  #axis.ticks.y = element_line( colour = "black" ), 
                  axis.title.y = element_text( size = 15 ),
                  axis.text.x = element_text( size = 12, face = "bold" ),
                  axis.text.y = element_text( size = 12, face = "bold" ),
                  strip.text = element_text( size = 15, face = "bold", colour = "black" ) ) 
         + scale_x_discrete(breaks = c( "AxA","AxC","AxG","AxT", "CxA","CxC","CxG","CxT", 
                                        "GxA","GxC","GxG","GxT", "TxA","TxC","TxG","TxT" ),
                            labels =c( 'A\nA',"\nC","\nG","\nT", 'C\nA',"\nC","\nG","\nT",
                                       'G\nA',"\nC","\nG","\nT", 'T\nA',"\nC","\nG","\nT" ) )
         + geom_rect( aes( xmin = "AxA", xmax = "TxT", ymin = mx + 0.01*mx, ymax = mx + 0.05*mx ) )
  )
  return(p)
  
}






