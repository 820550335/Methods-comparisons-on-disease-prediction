library(EBPRS)
library(optparse)
library(data.table)

option_list = list(
	make_option(c("--sst"),action="store",default=FALSE,type="character",help="Input path and filename of summary statistics file"),
	make_option(c("--out"),action="store",default=FALSE,type="character",help="Input path and filename of output"),
	make_option(c("--cc"),action="store",default=FALSE,type="character",help="Input path and filename of .frq.cc file (from plink)"))

opt = parse_args(OptionParser(option_list = option_list))
df1 <- fread(opt$sst,sep=' ')
ss <-data.table(df1)
ss <- ss[,c(-1,-3,-6,-8,-9,-(11:14))]
setcolorder(ss,c("A1","A2","OR","P","SNP"))
cc <- fread(opt$cc,sep=' ')
cc <- data.table(cc)
set1<-c()
for (i in seq_along(cc$SNP)){
	if (!(cc$SNP[i] %in% ss$SNP)){
		set1<-c(set1,i)
	}
}
cc <- cc[-set1]

out <- EBPRS(train=ss,N1=cc$NCHROBS_A/2,N0=cc$NCHROBS_U/2)
eff <- data.table(out$result$SNP,out$result$A1,out$result$effectsize)
write.table(eff,file=opt$out,sep=' ',row.names=FALSE,col.names=FALSE,quote=FALSE)
