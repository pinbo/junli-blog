dd = read.delim("QTL7AS.txt")

source("linkmap.R")

svg(file="QTl7AS-4.svg", width=3, height=4, bg=NA)
linkmap(dd,m.cex=1, ruler=F, interval = T)
dev.off()
