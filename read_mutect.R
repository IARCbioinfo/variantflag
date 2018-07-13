library(data.table)
args <-  commandArgs(trailingOnly=TRUE)
print(args[1])
vcf<-fread(args[1], skip = "CHROM")
vcf <- vcf[ FILTER == "PASS" ]
names( vcf )[1] <- "X.CHROM"
vcf <- vcf[ grep( "chrUn|alt|random|M", X.CHROM, invert = T ) ]

vcf <- vcf[ nchar( REF ) == 1 & nchar( ALT ) == 1 ]

Cov_N<-vcf[ , tstrsplit( NORMAL, ":", keep = 2 ) ][ , tstrsplit( V1, "," ) ]
Cov_N<-Cov_N[ , as.numeric(V1)+as.numeric(V2) ]

AF_N<- vcf[ , c("AF_N") := tstrsplit( NORMAL, ":", keep = 3 ) ][,AF_N]

Cov_T<-vcf[ , tstrsplit( TUMOR, ":", keep = 2 ) ][ , tstrsplit( V1, "," ) ]
Cov_T<-Cov_T[ , as.numeric(V1)+as.numeric(V2) ]

AF_T<- vcf[ , c("AF_T") := tstrsplit( TUMOR, ":", keep = 3 ) ][,AF_T]

# foo<-cbind(vcf[,1:5], Cov_N, Cov_T, AF_N, AF_T, vcf[,6:(ncol(vcf)-2)])
foo<-cbind( vcf[,1:5], Cov_N, Cov_T, AF_N, AF_T )
# write.table(foo, file=paste0( substr( args[1], 1, nchar(args[1])-14), "_pass_vcf"), sep="\t", row.names=F, quote=F)
write.table(foo, file=paste0( gsub( "_calls.vcf", "_pass_mutect.vcf", args[1] ) ), sep="\t", row.names=F, quote=F ) 