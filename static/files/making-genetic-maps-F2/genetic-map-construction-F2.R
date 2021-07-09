# 2021-05-18
# Map construciton

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
## Using ASMap to construct maps
list.files(".", "csv")
## directly read crossr data with a phsudo phenotype
# genotypes=c("A","H","B"): has to be c("A","H","B") NOT c("A","B","H"), seems corresponding AA, AB, BB. Maybe just leave them untouched

# snp <- read.cross("csvr", ".", "rqtl-F2-example.csv", genotypes=c("A","H","B"),alleles=c("A", "B"), map.function="kosambi", na.strings="N", F.gen = 2, BC.gen = 0)

# example input is for durum wheat (14 pairs of chromosomes) biparental F2 90K SNPs
snp <- read.cross("csvr", ".", "rqtl-F2-example.csv", map.function="kosambi", na.strings="N", F.gen = 2, BC.gen = 0) # to use ASMap, crosstype should not be F2, but BC0F2

#summary(snp)

##### now check and remove bad data
plotMissing(snp)

# identify the genotypes with a certain number of missing values
sg <- statGen(snp, bychr = FALSE, stat.type = "miss")
hist(sg$miss)
hist(sg$miss[sg$miss<1000])
sg$miss[sg$miss>1000]
sum(sg$miss<1000)
# only keep individuals with missing data < 1000
map1 <- subsetCross(snp, ind = sg$miss < 1000)

nind(map1) # number of individuals after filtering

# another way to check individuals with too many missing data for F2 pop
# par(mfrow=c(1,2), las=1)
plot(ntyped(snp), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(snp, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")
# map1 <- subset(snp, ind=(ntyped(snp)>6000))


## check highly related individuals
gc <- genClones(map1, tol = 0.95) # no clones
dim(gc$cgd) # rows are pairs of highly relatved genotypes
gc$cgd

# similarity histogram
hist(gc$cgm[lower.tri(gc$cgm)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes",main="Genotype Similarity")
rug(gc$cgm[lower.tri(gc$cgm)])

# based on the hisogram and table above
# if there are highly similar genotypes, remove them with fixClones function
# cgd <- gc$cgd
# map1 <- fixClones(map1, cgd, consensus = TRUE)
# nind(map1) #[1] 74

## check segregation distortion
# we can use bonferroni criteria = chi squre pvalue / number of markers < 0.05 to remove distorted markers
# profileMark(map1, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l", cex = 0.5)

gt <- geno.table(map1)
nrow(gt[gt$P.value < 0.05/totmar(map1),]) # number of distorted markers with bonferroni criteria

totmar(map1)

## pull out 3 types of markers before making genetic maps
# 1. colocated markers: markers with identical genotying data
map2 <- pullCross(map1, type = "co.located")
totmar(map2)
# 2. markers with too many missing data, I only keep markers with missing data <10%, we can push back markers with more missing data to the map after making the genetic maps
map2 <- pullCross(map2, type = "missing", pars = list(miss.thresh = 0.1)) # using high stringency here
totmar(map2)
# 3. markers that are highly distorted.
map2 <- pullCross(map2, type = "seg.distortion") #, pars = list(seg.thresh = "bonf"))
totmar(map2)

names(map2) # colocated, missing, seg.distrotion were added to store the pulled markers

# Check individuals' genotype frequencies
g <- pull.geno(map2)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))

par(mfrow=c(1,3), las=1)
for(i in 1:3) plot(gfreq[i,], ylab="Genotype frequency", main=c("A","H", "B")[i],ylim=c(0,1))
par(mfrow=c(1,1))

# Anothere way to see the distribution
dimnames(gfreq) <- list(c("A","H", "B"), 1:ncol(gfreq))
barplot(gfreq,main="Individuals’ genotype frequencies", xlab="Individuals", ylab="Genotype frequencies", col=c("skyblue","magenta","cyan"),legend=T, args.legend = list(x=92))

#########################
## sart map construction
# check different pvalues to see which one get the best grouping
# you can start from 1e-20 if you have >2000 markers
map3 <- mstmap(map2, bychr = FALSE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-15)
nchr(map3)
nmar(map3)
plotMap(map3, chr=(nmar(map3)>1))

maps = pull.map(map3, chr=(nmar(map3)>1), as.table=T)
write.table(maps,"maps_map3_for_checking.txt",sep="\t", col.names=NA,quote=F)
check_chrom(maps) # check output file "chromosome_check.txt" to see whether all chromosomes have been sepated and which LG is which chromosome, so we can regroup them together.

# check ind with too many crossovers
# looks okay
dum = map3
class(dum) = c("f2", "cross")
plot(countXO(dum), ylab="Number of crossovers")
abline(h=100,col="red")

## break and merge chromosomes
## for the example input, 2A, 2B, 3B, 5A and 7A need to merge LG
# check whether there are linages from heat map
# 2A: not very strong
heatMap(map3, chr = c("L.10","L.15"),lmax=4)
# 2B: good
heatMap(map3, chr = c("L.11","L.12","L.14"),lmax=4)
# 3B: not very good
heatMap(map3, chr = c("L.7","L.20"), lmax=4)
# 5A: good
heatMap(map3, chr = c("L.13","L.6"),lmax=4)
# 7A: not strong
heatMap(map3, chr = c("L.2","L.21"),lmax=4)


# if everything looks good, then merge
map4 <- mergeCross(map3, merge = list(
  "2A" = c("L.10","L.15"),
  "2B" = c("L.11","L.12","L.14"),
  "3B" = c("L.7","L.20"),
  "5A" = c("L.13","L.6"),
  "7A" = c("L.2","L.21")
))

nchr(map4)
nmar(map4)
names(map4$geno)
names(map4$geno) = c("2A","2B","3B","5A","7A","1A","4A","3A","L.18","6B","L.22","1B","4B","6A","5B","7B") # rename them to chromosome names


## We then need to break L.1 (= 1A + 1B)
map5 <- mstmap(map4, bychr = T, trace = TRUE, dist.fun = "kosambi", p.value = 1e-20, chr = "1A")
# map6 <- mstmap(map4, bychr = FALSE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-30, chr = "L.1")
nchr(map5)
nmar(map5)
pull.map(map5,as.table=T,chr = "1A.3")
pull.map(map5,as.table=T,chr = "1A.2")
write.table(maps,"maps_map5_for_checking.txt",sep="\t", col.names=NA,quote=F)

## merge chromosomes
# 1B
heatMap(map5, chr = c("1B","1A.2","1A.3"),lmax=4)

# merge
map6 <- mergeCross(map5, merge = list(
  "1B" = c("1B","1A.2","1A.3")
))
nmar(map6)
names(map6$geno)
names(map6$geno)[1] = "1A"
# reorder markers
map7 <- mstmap(map6, bychr = T, trace = TRUE, dist.fun = "kosambi", p.value = 2)
map7 = subsetCross(map7, chr=(nmar(map7)>1))
plotMap(map7)

maps = pull.map(map7,as.table=T)
check_orientation(maps)
write.table(maps,"maps_map7_for_checking.txt",sep="\t", col.names=NA,quote=F)

# check lines with too many crossovers
dum = map7
class(dum) = c("f2", "cross")
plot(countXO(dum), ylab="Number of crossovers") # no bad lines

# check whether 4A and 7A have close linked markers
heatMap(map7, chr = c("4A","7A"),lmax=6)

## check
pg <- profileGen(map7, bychr = FALSE, stat.type = c("xo", "dxo","miss"), id = "Genotype", xo.lambda = 42, layout = c(1, 3), lty = 2, cex = 0.7)
profileMark(map7, stat.type = c("seg.dist","dxo", "recomb"), layout = c(1, 3), type = "l")

# check potential switched alleles
dum = est.rf(map7) # est.rf can check markers with potential switched alleles
rf <- pull.rf(dum)
lod <- pull.rf(dum, what="lod")
max(rf,na.rm = T)
# ploting takes too long
# plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")
# reorder 7A
# map8 = reorder.marker(map7, "7A", c(50:1, 117:51))

## remove some markers causing big gaps
# usually at ends of near the border of big gaps
totmar(map7)
badmar = c("IWB24549_6D_7.1", "IWB3832_3D_559.1", "IWB31740_1D_474")
map9 <- drop.markers(map7, badmar)
map9 = mstmap(map9, bychr = T, trace = TRUE, dist.fun = "kosambi", p.value = 2)
maps = pull.map(map9,as.table=T)
check_orientation(maps)

# flip some chromosomes after orientation check compared to 90K consensus
map10 = flip.order(map9, c("1A", "3A", "3B", "4B", "5A", "6B", "7A", "7B"))
maps = pull.map(map10,as.table=T)
check_orientation(maps)
plotMap(map10)
write.table(maps,"maps_final_without_coolocated.txt",sep="\t", col.names=NA,quote=F)

map <- pushCross(map10, type = "co.located")
totmar(map)
maps = pull.map(map,as.table=T)
check_orientation(maps)
write.table(maps,"maps_map_final_with_coolocated_pushed.txt",sep="\t", col.names=NA,quote=F)

save(map10,file="map10.RData") # no identical markers
write.cross(map10, "csvr", filestem = "maps_no_coolocated")
write.cross(map, "csvr", filestem = "maps_all_markers")

######################## check potential genotyping errors ############################
### use detectBadData
# not work well for F2
# map11 =  mstmap(map10, bychr = T, trace = TRUE, dist.fun = "kosambi", p.value = 2, detectBadData = T, return.imputed = T)

# use r/qtl to correct potential genotyping errors
# does not work well for F2
# mapthis = map10
# mapthis <- calc.errorlod(mapthis, error.prob=0.005, map.function="kosambi")
# toperr <- top.errorlod(mapthis, cutoff=2.5)
# dim(toperr)
# head(toperr)

#### use my own function to check genotyping error for F2 population
# okay for F2: hhhh-a-hhhh or hhhh-b-hhhh or aaaa-h-bbbb or bbbb-h-aaaa (generated by close CO in the two crhomosomes coming from male and female meiosis)
#  bad for F2: aaaa-h-aaaa or bbbb-h-bbbb or bbbb-a-bbbb or aaaa-b-aaaa 
# If they are at the end of chromosome there are not a problem, because they are not double crossovers (or in the border of a big gap[ (e.g. >20 cM)

# 1 = A, 2 = H, 3 = B in the geno table
check.geno.err = function(x) {
  # first convert NA to 9 for easy comparison
  x[is.na(x)] = 9
  dum = sapply(2:(length(x)-1), function(i){
    a = x[i-1]
    b = x[i]
    c = x[i+1]
    if (b==9 | a==b | b==c) return(0)
    else if (b == 1) {# b is A
      if (a == 3 | c == 3) return (1)
      else return(0)
    } else if (b == 2) { # b is H
      if (a == c & (a == 1 | a == 3)) return(1)
      else return(0)
    } else {# b is B
      if (a == 1 | c == 1) return (1)
      else return(0)
    }
  })
  return(c(0,dum,0))
}

tt = map10$geno$`1A`$data
cc = apply(tt, 1, function(x) check.geno.err(x))
which(cc==1, arr.ind=TRUE)

dum = subsetCross(map10, chr = nmar(map10)>1)
# get positions of genotyping error for all chromosomes
ll = lapply(dum$geno, function(x){
  tt = x$data
  cc = apply(tt, 1, function(x) check.geno.err(x))
  which(cc==1, arr.ind=TRUE)
})

sum(sapply(ll, nrow)) # 310 total geno err

# get the absolute positions for converting them to NA
map11 = subsetCross(map10, chr = nmar(map10)>1)
ap = lapply(dum$geno, function(x){
  tt = x$data
  cc = apply(tt, 1, function(x) check.geno.err(x))
  t(cc)
})

cn = chrnames(map11) # chromosome names
for (x in cn){
  map11$geno[[x]]$data[which(ap[[x]]==1)] <- NA
}
# re-estimate maps
map12 = quickEst(map11)

plotMap(map11, map12) # better now
plotMap(map12)
write.cross(map12, "csvr", filestem = "map12_with_Phy_no_coolocated_geno_corrected")
# use shell command to change to ABG coding
# sed  -i 's/AA/A/g; s/AB/H/g; s/BB/B/g' map12_with_Phy_no_coolocated_geno_corrected.csv

maps = pull.map(map12,as.table=T)
check_orientation(maps)
write.table(maps,"maps_final_without_coolocated_pushed_geno_corrected.txt",sep="\t", col.names=NA,quote=F)
save(map12,file="map12-geno-corrected.RData") # no identical markers

# get the map for all markers
map <- pushCross(map12, type = "co.located")
totmar(map)
maps = pull.map(map,as.table=T)
check_orientation(maps)
write.table(maps,"maps_map_final_with_coolocated_pushed_geno_corrected.txt",sep="\t", col.names=NA,quote=F)



