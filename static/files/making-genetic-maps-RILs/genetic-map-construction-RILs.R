# 2021-05-22
# Map construction

library(ASMap)
library(lattice)

##### functions

## check marker distribution on physical maps after linkage group separation
## help assign LG to chromosomes
# marker name should be like this: IWB65373_1A_3.8 = marker name, physical chromosome and Mb, separated by "_"
check_chrom = function(maps){
  aa = strsplit(rownames(maps), "_")
  cc = sapply(aa, function(x) x[2]) # chrom names
  maps$chrom = cc
  # xtabs(~chr+chrom,maps)
  write.table(xtabs(~chr+chrom,maps), "chromosome_check.txt", sep="\t")
}

## check genetic orientation compared to physical order
# marker name should be like this: IWB65373_1A_3.8 = marker name, physical chromosome and Mb, separated by "_"
check_orientation = function(maps){
  aa = strsplit(rownames(maps), "_")
  bb = sapply(aa, function(x) as.numeric(x[3]))
  cc = sapply(aa, function(x) x[2]) # chrom names
  maps$Mb = bb
  print(xyplot(Mb~pos|chr, data=maps[maps$chr==cc,],as.table=TRUE)) # check orientation and make sure to use markers on the same chrom
  write.table(maps, "orientation_check.txt", col.names=NA, quote=F)
}

## recombination frequency to Kosambi distance
kosambi = function(y) 0.25*log((1+2*y)/(1-2*y))

## pull imputed marker data
pull.imputed.geno <- function (cross, chr) 
{
  if (!inherits(cross, "cross")) 
    stop("Input should have class \"cross\".")
  if (!missing(chr)) 
    cross <- subset(cross, chr = chr)
  X <- cross$imputed.geno[[1]]$data
  if (nchr(cross) > 1) 
    for (i in 2:nchr(cross)) X <- cbind(X, cross$imputed.geno[[i]]$data)
  X
}

## function to reorder a map manually
reorder.marker = function(cross, chr, new.marker.order) {# new.marker.order can be either numbers or marker names
  genodata = cross$geno[[chr]]$data
  mapdata = cross$geno[[chr]]$map
  new.geno = genodata[,new.marker.order]
  new.map = mapdata[new.marker.order]
  cross$geno[[chr]]$data = new.geno
  cross$geno[[chr]]$map = new.map
  cross2 = quickEst(cross, chr)
  print(pull.map(cross2, chr, as.table = T))
  #new.pos = est.map(cross, chr, map.function="kosambi")[[chr]]
  #print(new.pos)
  #ll = length(new.pos) # n marker
  #cross$geno[[chr]]$map[1:ll] = new.pos[1:ll]
  return(cross2)
}

## replace too big interval (>100) with 100
decrease.interval = function(cross, chr, threshold = 100) {
  mm = pull.map(cross, chr)
  mm1 = lapply(mm, function(x) {
    ll = length(x) # map length
    ii = x[2:ll] - x[1:(ll-1)] # intervals
    ii[ii>threshold] = threshold
    y = cumsum(c(x[1], ii))
    x[1:ll] = y[1:ll]
    return(x)
  })
  return(mm1)
}


#########################
list.files(".", "csv")

# example input is for durum wheat (14 pairs of chromosomes) biparental F7 90K SNPs
snp <- read.cross("csvr", ".", "rqtl-RILs-example.csv", map.function="kosambi", na.strings="-", crosstype="riself")

##### now check and remove bad data
plotMissing(snp)

# identify the genotypes with a certain number of missing values
sg <- statGen(snp, bychr = FALSE, stat.type = "miss")
hist(sg$miss)
# hist(sg$miss[sg$miss<1000])
sg$miss[sg$miss>500]
sum(sg$miss>1000)
map1 <- subsetCross(snp, ind = sg$miss < 1000) # only keep individuals with missing data < 1000
snp$pheno$Genotype[sg$miss >= 1000] # names of individuals with missing data > 1000
# [1] RIL-154 RIL-156
nind(map1)
# another way to check individuals with too many missing data
plot(ntyped(map1), ylab="No. typed markers", main="No. genotypes by individual")
# check markers with too many missing data
plot(ntyped(map1, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")
# map1 <- subset(snp, ind=(ntyped(snp)>6000))

## check highly related individuals
gc <- genClones(map1, tol = 0.9) # get line pairs with similarity > 0.9 or 0.95
dim(gc$cgd)
gc$cgd # similar pairs
write.table(gc$cgd, "very_similar_lines_0.9.txt",sep="\t")
# similarity histogram, help determine which level to use, most of time 0.95 is good
hist(gc$cgm[lower.tri(gc$cgm)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes",main="Genotype Similarity")
rug(gc$cgm[lower.tri(gc$cgm)])


# based on the histogram and table above
# I will remove similar lines
cgd <- gc$cgd
map1.1 <- fixClones(map1, cgd, consensus = TRUE)
nind(map1.1) #[1] 143
nind(map1)
## check segregation distortion
profileMark(map1.1, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l", cex = 0.5)

gt <- geno.table(map1.1)
nrow(gt[gt$P.value < 0.05/totmar(map1),]) # number of markers that have significant seg distoration at bonf level 0.05
nrow(gt[gt$P.value < 0.05,]) # no bonf correction

totmar(map1.1)
## pull out 3 types of markers before making genetic maps
# 1. markers with too many missing data, usually use 10% as the threshold
map2 <- pullCross(map1.1, type = "missing", pars = list(miss.thresh = 0.1)) # using high stringency here
totmar(map2)
# 2. markers that are colocated
map2 <- pullCross(map2, type = "co.located")
totmar(map2)
# 3. markers that are highly distorted
map2 <- pullCross(map2, type = "seg.distortion", pars = list(seg.thresh = "bonf"))
# map2 <- pullCross(map2, type = "seg.distortion") # no bonf correction
totmar(map2) # 6220
names(map2)

# check individual genotype frequencies
g <- pull.geno(map2)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:2)))
gfreq <- t(t(gfreq) / colSums(gfreq))

par(mfrow=c(1,2), las=1)
for(i in 1:2) plot(gfreq[i,], ylab="Genotype frequency", main=c("A", "B")[i],ylim=c(0,1))
par(mfrow=c(1,1))

# a better plot
dimnames(gfreq) <- list(c("A","B"), 1:ncol(gfreq))
barplot(gfreq,main="Individuals’ genotype frequencies", xlab="Individuals", ylab="Genotype frequencies", col=c("skyblue","magenta","cyan"),legend=T, args.legend = list(x=92))

#########################
## start map construction
# may try different p values to separate the markers properly, 1e-10 is good start
map3 <- mstmap(map2, bychr = FALSE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-8)
nchr(map3)
nmar(map3)
sum(nmar(map3)>2)
plotMap(map3, chr=(nmar(map3)>4))

# check lines with too many crossovers
# pg <- profileGen(map3.1, bychr = FALSE, stat.type = c("xo", "dxo","miss"), id = "Genotype", xo.lambda = 168, layout = c(1, 3), lty = 2, cex = 0.7) # 168=14x2x6 for F7 tetraploid wheat
# map3.2 <- subsetCross(map3.1, ind = !pg$xo.lambda)
# nind(map3)


# check ind with too many crossovers
plot(countXO(map3), ylab="Number of crossovers")
abline(h=400,col="red")
map3$pheno$Genotype[countXO(map3) > 350] #bad lines
# [1] RIL-11  RIL-116 RIL-140 RIL-155 RIL-158 RIL-157 RIL-159
map3.1 <- subsetCross(map3, ind=(countXO(map3) < 350))
map3.1 = subsetCross(map3.1, chr=(nmar(map3.1)>4))
nchr(map3.1)
# after removing some lines, more markers will be colocated
# but we need to push back the old co.located first, then pull again
map3.2 <- pushCross(map3.1, type = "co.located")
names(map3.2)
totmar(map3.2)
map3.2 <- pullCross(map3.2, type = "co.located")
totmar(map3.2) # 4579
totmar(map3.1)
map3.2 <- mstmap(map3.2, bychr = T, trace = TRUE, dist.fun = "kosambi", p.value = 1e-10)

totmar(map3.2)
nmar(map3.2)
nchr(map3.2)
maps = pull.map(map3.2,as.table=T, chr=(nmar(map3.2)>1))
check_chrom(maps)
write.table(maps,"maps_map3.2_for_checking.txt",sep="\t", col.names=NA,quote=F)
# check poteintal switched alleles
dum <- est.rf(map3.2) # most of them belong to groups with 1 marker
# In est.rf(map3.1) : Alleles potentially switched at markers 
# IWB41059_5A_270.83

## break and merge chromosomes
# 1A 1B: L.10 remove
# heatMap(dum, chr = c("L.3.1","L.3.2", "L.1.4"),lmax=4)


# merge
map4 <- mergeCross(map3.2, merge = list(
  "3B"=c("L.11.2","L.9.1"),
  "5A"=c("L.15.1","L.15.2","L.15.3"),
  "7A"=c("L.8.1","L.8.2","L.8.3"),
  "6B"=c("L.5.1","L.5.2"),
  "7B"=c("L.13.1","L.13.3"),
  "1A"=c("L.2.2","L.2.5"),
  "2A"=c("L.11.1","L.7.1","L.7.2"),
  "3A"=c("L.14.1","L.70")
))
nmar(map3.2)
nmar(map4)
map4 = subsetCross(map4, chr=(nmar(map4)>20)) # only keeps chromosomes with >20 markers
nchr(map4)
# reorder markers by chromosome, set p.value >= 1
map4 <- mstmap(map4, bychr = T, trace = TRUE, dist.fun = "kosambi", p.value = 2)

nchr(map4)
nmar(map4)
names(map4$geno)
names(map4$geno)[9:14] = c("1B", "4B", "4A", "6A", "5B", "2B")

maps=pull.map(map4,as.table=T)
write.table(maps,"maps_map4_for_checking.txt",sep="\t", col.names=NA,quote=F)
check_orientation(maps)

dum <- est.rf(map4)
pdf(file="heatmap-map4-homeolog-chrom.pdf")
plotRF(map4, chr=c("1A","1B"), alternate.chrid=TRUE, zmax=6)
plotRF(map4, chr=c("2A","2B"), alternate.chrid=TRUE, zmax=6)
plotRF(map4, chr=c("3A","3B"), alternate.chrid=TRUE, zmax=6)
plotRF(map4, chr=c("4A","4B"), alternate.chrid=TRUE, zmax=6)
plotRF(map4, chr=c("5A","5B"), alternate.chrid=TRUE, zmax=6)
plotRF(map4, chr=c("6A","6B"), alternate.chrid=TRUE, zmax=6)
plotRF(map4, chr=c("7A","7B"), alternate.chrid=TRUE, zmax=6)
dev.off()

## check
pg <- profileGen(map4, bychr = FALSE, stat.type = c("xo", "dxo","miss"), id = "Genotype", xo.lambda = 168, layout = c(1, 3), lty = 2, cex = 0.7)

pdf(file="profileMark.pdf", width = 32)
profileMark(map4, stat.type = c("seg.dist","prop","dxo", "recomb","miss"), layout = c(1, 6), type = "l")
dev.off()

plot(countXO(map4), ylab="Number of crossovers")
# map5$pheno$Genotype[countXO(map5) > 400] #bad lines
# map6 <- subsetCross(map5, ind=(countXO(map5) < 400))


## check big gaps and remove some bad markers
# I also checked the original genome studio file to make sure these markers have bad clusters
badmar = c("IWB181_7A_675.1","IWB847_7B_5.51","IWB1385_2A_751.89","IWB1471_7D_622","IWB2881_5B_692.68","IWB3472_3A_12.38","IWB4084_3D_7.59","IWB4311_5A_6.2","IWB4477_2A_712.86","IWB4669_7D_5.34","IWB4864_2B_761.99","IWB5822_3D_612.27","IWB6561_Un_172.24","IWB6942_7A_691.26","IWB7396_1D_16.79","IWB8950_5B_387.9","IWB9132_2D_309.58","IWB9672_1B_250.8","IWB9688_1A_570.62","IWB9790_7D_603.66","IWB10550_3B_556.94","IWB11213_1B_16.7","IWB11284_1B_484.38","IWB12290_2B_793.13","IWB13518_7D_587.1","IWB13543_5D_258.7","IWB15732_4D_124.53","IWB15736_1D_486.9","IWB15832_4D_73.37","IWB16198_6D_472.3","IWB16957_6D_58.69","IWB21023_1D_66.2","IWB21342_5A_693.3","IWB22096_6B_673.77","IWB22300_7A_100.76","IWB22524_2B_716.8","IWB23769_2D_638.6","IWB24059_1A_155.08","IWB25399_4A_722.9","IWB27227_1A_9.11","IWB27349_7A_251.8","IWB28576_1A_8.3","IWB28781_3B_760.1","IWB29606_2D_638.73","IWB30246_1B_603.3","IWB30456_1B_22.71","IWB32219_2A_564.5","IWB32358_3B_583.8","IWB33553_2B_773.35","IWB33855_6D_4.78","IWB34073_4A_645.4","IWB34301_7A_678.16","IWB34743_unknown_0","IWB35970_7A_62.5","IWB36187_1B_28.56","IWB36283_1D_448.66","IWB36367_1B_2.1","IWB38040_2B_761.96","IWB38866_2A_293.89","IWB39072_4B_544.45","IWB40149_2D_250.3","IWB40409_5B_406.17","IWB40939_2A_306.2","IWB41436_1A_584.4","IWB41851_6D_3.73","IWB43627_1B_269.22","IWB43707_2A_775.68","IWB44409_1A_405.78","IWB44632_2B_9.73","IWB44773_2A_762.3","IWB45086_5A_155.22","IWB45291_4B_347.6","IWB45528_2A_171.8","IWB46179_7B_1.28","IWB47794_1D_11.1","IWB47980_1B_22.71","IWB50287_2B_701.76","IWB50552_6A_194.72","IWB51034_2D_2.05","IWB51355_7A_49.1","IWB52231_6D_256.6","IWB54379_3D_597.5","IWB54430_6B_704.19","IWB54670_5B_531.41","IWB55453_7A_8.7","IWB56191_6A_25.63","IWB57562_3B_660.79","IWB58613_5A_617.98","IWB59330_3B_31.81","IWB59680_2B_769.91","IWB60043_7A_26.88","IWB60250_2B_66.19","IWB60694_3A_725.52","IWB62045_2D_161.7","IWB62515_4A_79.05","IWB63836_4A_170.12","IWB65017_Un_88.2","IWB65312_3B_491.8","IWB65331_5A_547.17","IWB66309_5B_580.51","IWB68225_3A_372.46","IWB68566_7B_680.21","IWB68994_7B_713.86","IWB69115_1A_362.41","IWB69200_2A_695.52","IWB69780_7B_254.12","IWB70040_2B_523.8","IWB70644_5B_248.7","IWB70831_3A_455.51","IWB71657_1D_487.84","IWB71931_1D_165.9","IWB72787_1B_28.56","IWB73626_1A_558.59","IWB74817_7B_743.62","IWA1312_3D_600.2","IWA3562_7A_675.11","IWA4767_5A_19.23","IWA5203_3A_10.3","IWA5222_3D_606.88","IWA7849_Un_82.52","IWB48008_7A_134.63")

map5 <- drop.markers(map4, badmar)
map5 = mstmap(map5, bychr = T, trace = TRUE, dist.fun = "kosambi", p.value = 2)
plotMap(map5)
maps = pull.map(map5,as.table=T)
check_orientation(maps)
check_orientation(maps, F)
write.table(maps,"map5_for_checking.txt",sep="\t", col.names=NA,quote=F)

# check potential switched alleles
dum = est.rf(map5)
rf <- pull.rf(dum)
lod <- pull.rf(dum, what="lod")
# plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score") # takes too long, I just need to know the maximum of rf
max(rf,na.rm = T)

# check individual
profileGen(map5, bychr = FALSE, stat.type = c("xo", "dxo", "miss"), id ="Genotype", xo.lambda = 120, layout = c(1, 3), lty = 2)
# check marker
profileMark(map5, stat.type = c("seg.dist", "dxo", "erf", "lod"), id ="Genotype", layout = c(1, 4), type = "l")

# flip some chromosomes after orientation check compared to 90K consensus
map6 = flip.order(map5, c("1B","2A", "3B", "4A", "4B", "5A", "5B", "6B","7A","7B"))
plotMap(map6)
maps = pull.map(map6,as.table=T)
check_orientation(maps)
write.cross(map6, "csvr", filestem = "map6_with_Phy_no_coolocated")

write.table(maps,"maps_final_NO_coolocated.txt",sep="\t", col.names=NA,quote=F)

map <- pushCross(map6, type = "co.located")
totmar(map)
maps = pull.map(map,as.table=T)
check_orientation(maps)
write.table(maps,"maps_final_with_coolocated_pushed.txt",sep="\t", col.names=NA,quote=F)

save(map,file="map.RData")
save(map6,file="map6_no_colocated.RData") # no identical markers
write.cross(map, "csvr", filestem = "map_final_all_markers_with_Phy")

write.table(map6$pheno$Genotype, "Lines-used-to-construct-maps.txt", sep="\t")


# use detecterr
map7 = mstmap(map6, bychr = T, trace = TRUE, dist.fun = "kosambi", p.value = 2, detectBadData = T)
maps = pull.map(map7,as.table=T)
check_orientation(maps)

map7 = flip.order(map7, c("2A", "3B", "4A", "6A","7B"))
plotMap(map6)
maps = pull.map(map7,as.table=T)
check_orientation(maps)

plotMap(map6, map7)
plotMap(map6)
plotMap(map7)
write.table(maps,"maps_final_map7_detectErr_No_coolocated.txt",sep="\t", col.names=NA,quote=F)
