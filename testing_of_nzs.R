library(regioneReloaded)
library(ggplot2)
library(patchwork)

bestScore<-function(GR){
  max(GR$score)
  GR<-GR[GR$score>=quantile(GR$score)[4],]
  return(GR)
}


sets <- readRDS("/Users/ilmalli/Downloads/permResENCODE.RDS")

setX<-lapply(sets,bestScore)
setX<-lapply(setX,filterChromosomes) # setX is the databases "cleaned" only with sequence wit > 75% of score
 newSets<-lapply(setX,sample,size=min(sapply(setX,length)))  # newSets is a rarefacted sample in which all the GRanges has the same number of regions (the one of the smallest)

gR<-GRanges() # gR is the genome calculated for the regions of newSets
for(i in 1:length(newSets)){ 
  gR<-mergeRegions(gR,newSets[[i]])
}

expN<-round(mean(width(gR1))/15)

gR1<-GRanges() # gR1 is the genome calculated for the regions of estX more 15% more
for(i in 1:length(setX)){
  gR1<-mergeRegions(gR1,setX[[i]])
  mean(width(gR1))
  gR1<-extendRegions(gR1,extend.end = expN, extend.start = expN)
}

A<-setX$JUN_ENCFF079UOJ # test for one region set

vecLength<- round(seq(10,length(A),length.out=10))
listA<-list(A1 = A[sample(vecLength[1])],
            A2 = A[sample(vecLength[2])],
            A3 = A[sample(vecLength[3])],
            A4 = A[sample(vecLength[4])],
            A5 = A[sample(vecLength[5])],
            A6 = A[sample(vecLength[6])],
            A7 = A[sample(vecLength[7])],
            A8 = A[sample(vecLength[8])],
            A9 = A[sample(vecLength[9])],
            A10 = A[sample(vecLength[10])])


cwPT2<-crosswisePermTest(Alist = listA,Blist = setX, ntimes=100, ranFUN = "resampleGenome", genome =gR1) # is possible to try with different ranFUN
# cwPT<-crosswisePermTest(Alist = listA,Blist = setX, ntimes=100, ranFUN = "resampleGenome", genome ="hg19") # activate if you want to use all vs hg19

cwPT2<-makeCrosswiseMatrix(cwPT2,zs.type = "z_score")
cwPT2n<-makeCrosswiseMatrix(cwPT2,zs.type = "norm_zscore")

matX<-t(cwPT@matrix$GMat)[c(paste0("A",1:10)),] matrix calculated for 
matX<-as.data.frame(matX[order(rownames(matX)),])

matXn<-t(cwPTn@matrix$GMat)
matXn<-as.data.frame(matXn[order(rownames(matXn)),])

nameX<-colnames(matX)[1] # name of Granges that we want to compare with 

vecLength<-round(seq(10,length(A),length.out=10))
dt_test<-data.frame(frac=paste0("V",1:10),zs=matX[,nameX],nzs=matXn[,nameX],nreg=vecLength)

p1<-ggplot(data=dt_test,aes(x=nreg,y = zs),) +
  geom_line( color="red")+
  geom_point(color="red") + ggtitle(paste0("CREB1_ENCFF314KET vs ",nameX) )+
  xlab("num. of regions") + ylab("z-score") + theme(plot.title = element_text(size=8))

p2<-ggplot(data=dt_test,aes(x=nreg,y = nzs),) +
  geom_line( color="green")+
  geom_point(color="green") + ggtitle(paste0("CREB1_ENCFF314KET vs ",nameX) )+
  xlab("num. of regions") + ylab("norm z-score") +  theme(plot.title = element_text(size=8)) + ylim(0,20)

p1|p2 # plot results
