# pass in argument infile, outfile and plotTitle

d <- read.table(infile, col.names=c('BarcodeID','BarcodeSeq','Count'))
d = d[order(d$Count),]

bc_tot <- sum(d$Count)
bc_assigned <- sum(d[d$BarcodeSeq!="NONE",]$Count)

p_title <- sprintf("Barcode distribution (%s)\n%d out of %d assigned (%.2f%%)", plotTitle, bc_assigned, bc_tot, bc_assigned/bc_tot*100)

postscript(file=sprintf("%s.ps", outfile))
par(mar=c(10,4,4,2))
with(d, barplot(d$Count, names=d$BarcodeSeq, las=3, main=p_title))
dev.off()
