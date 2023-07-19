#Run this to test the effect of normal ZScore vs Zscore


library(patchwork)


distT<-runif(100000,min = 5,max = 30)
distT2<-runif(100000,min = 10,max = 50)
distT3<-runif(100000,min = 1,max = 20)
distT4<-runif(100000,min = 7,max = 35)

calMean<-function(distT,n=10){
  res<-sum(sample(distT,n,replace = TRUE) > 18)
  return(res)
}

calculateZscore<-function(N=1000,d1,d2){
  print(N)
  ab<-calMean(d2,N)
  rdist<-sapply(rep(N,1000),FUN =calMean,distT=d1)
  zscore<-(ab-mean(rdist))/sd(rdist)
  nzs<-zscore/sqrt(N)
  ls<-c(zs=zscore,nzs=nzs,N=N)
  return(ls)
}

step=2000

lp<-lapply(seq(10,50000,by=step),calculateZscore,d1=distT,d2=distT)
lp2<-lapply(seq(10,50000,by=step),calculateZscore,d1=distT,d2=distT2)
lp3<-lapply(seq(10,50000,by=step),calculateZscore,d1=distT,d2=distT3)
lp4<-lapply(seq(10,50000,by=step),calculateZscore,d1=distT,d2=distT4)

df<-do.call("rbind",lp) %>% as.data.frame %>% mutate(sample="df",) 
df2<-do.call("rbind",lp2) %>% as.data.frame %>%mutate(sample="df2")
df3<-do.call("rbind",lp3) %>% as.data.frame %>%mutate(sample="df3")
df4<-do.call("rbind",lp4) %>% as.data.frame %>%mutate(sample="df4")

df<-rbind(df,df2,df3,df4)

pn <-df %>% ggplot(aes(x=N,y=zs,group.by=sample,color=sample)) +
  geom_line() +ggtitle("zscore") 

p1n <- df %>% ggplot(aes(x=N,y=nzs,group.by=sample,color=sample)) +
  geom_line() +ggtitle("normalized zscore") 

p <-df %>% ggplot(aes(x=N,y=zs,group.by=sample,color=sample)) +
  geom_line() +ggtitle("zscore") + scale_x_log10()

p1 <- df %>% ggplot(aes(x=N,y=nzs,group.by=sample,color=sample)) +
  geom_line() +ggtitle("normalized zscore") + scale_x_log10()
(pn| p1n) /(p | p1)
