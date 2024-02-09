ll<-c(1,1,1)

KK<-combn(100,2)

Delta<-rep(c(1,rep(0,24)),4)
Stage<-cumsum(Delta)
Stage_P<-(Stage[KK[1,]]==Stage[KK[2,]])

JJ1<-NULL
JJ2<-NULL
JJ3<-NULL
for(wa in 1:3)
{
  aaa = c("I","II","III")[wa]
  Xuan<-letters[c(ll[wa]:10)]
  Ya<-length(Xuan)
  Result<-NULL
  for(a in Xuan)
  {
    Result<-rbind(Result,read.csv(paste(c("Day100",aaa,"a_",a,".csv"),collapse = ""),header = T)[,-1])
  }
  Stage_dist<-t(apply(Result[1+0:(Ya*10-1)*6,],1,cumsum))
  ARI<-NULL
  MI<-NULL
  for(i in 0:(Ya*10-1)+1)
  {
    Stage_P_i<-(Stage_dist[i,KK[1,]]==Stage_dist[i,KK[2,]])
    
    
    TP<-mean(Stage_P*Stage_P_i)
    FP<-mean((1-Stage_P)*Stage_P_i)
    FN<-mean(Stage_P*(1-Stage_P_i))
    TN<-mean((1-Stage_P)*(1-Stage_P_i))
    Enumerator<-TP+TN-(TP+FP)*(TP+FN)-(TN+FP)*(TN+FN)
    Denominator<-1-(TP+FP)*(TP+FN)-(TN+FP)*(TN+FN)
    ARI<-c(ARI,Enumerator/Denominator)
    
    TT<-table(Stage,Stage_dist[i,])
    row_TT<-rowSums(TT)
    col_TT<-colSums(TT)
    MI<-c(MI,sum(TT/100*log((TT+(TT==0))*100/(row_TT%*%t(col_TT)))))
  }
  JJ1<-c(JJ1,as.character(round(mean(ARI),3)),"(",as.character(round(sd(ARI),3)),")")
  JJ2<-c(JJ2,as.character(round(mean(MI),3)),"(",as.character(round(sd(MI),3)),")")
  # JJ3<-c(JJ3,as.character(round(mean(Stage_dist[,100]),3)),"(",as.character(round(sd(Stage_dist[,100]),3)),")")
}

write.table(t(c(JJ1,JJ2)),"Quan_large.txt",sep=" ")