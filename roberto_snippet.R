# Code snippet from Roberto to subsample a list of granges to the set with the 
# minimum number of regions in the list (making all the sets the same size)

bestScore<-function(GR){
  max(GR$score)
  GR<-GR[GR$score>=quantile(GR$score)[4],]
  return(GR)
}

sets<-readRDS("/Users/ilmalli/Downloads/ENCODE_sets.RDS")
setX<-lapply(sets,bestScore)
setX<-lapply(setX,filterChromosomes)
newSets<-lapply(setX,sample,size=min(sapply(setX,length)))

gR<-GRanges()
for(i in 1:length(newSets)){
  gR<-mergeRegions(gR,newSets[[i]])
}