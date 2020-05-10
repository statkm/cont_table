#simulation (not using parallel)
#exact logrank simulation

#setwd("/Users/mizumakoutarou/Desktop/codes")
#getwd()
#source(file="cas/MCMC/MetHalg.R") #提案MCMC
#source(file="cas/cmhtest/mehta1985.R") #exact
#source(file="cas/logrank/cont_table_surv.R") #survival to contingency table

seednum <- (runif(1) * 1000) %>% ceiling
set.seed(seednum)

N1<-50; N2<-50

odds_true <- 1 ; odds_null <- 1

name_dist <- "exponential"
# "geometric" "exponential"

lambda<-.9 ; 
lambda1<-lambda*odds_true; lambda2<-lambda

Nsim <- 10000

pval_mh <- matrix(0,nrow=Nsim,ncol=2)

pb <- txtProgressBar(min = 1, max = Nsim, style = 3)

for(i in 1:Nsim){
  #  cat(i,"-th simulation\n")
  if(name_dist=="exponential"){
    x1test<-rexp(n = N1,rate=lambda1)
    x2test<-rexp(n = N2,rate=lambda2)
  }else{
    if(name_dist=="geometric"){
      x1test<-rgeom(n = N1,prob=lambda1)+1
      x2test<-rgeom(n = N2,prob=lambda2)+1
    }else{
      cat("error")
    }
  }


  
  x1delta <- rbinom(N1, 1, 1)
  x2delta <- rbinom(N2, 1, 1)
  
#  cont_logrank<-cont_surv(cbind(x1test,x1delta),cbind(x2test,x2delta))
#  res_MCMC <- Methalg(Nmcmc,cont_logrank,odds=odds_null,Nburnin,method="each", midp = T, alternative = "two.sided")
  #  res_mh2 <- Methalg(Nmcmc,cont_logrank,odds=1,Nburnin,method="each", midp = T)$p.value
  #  res_lr <- mantelhaen.test(cont_logrank)$p.value
  #  res_mh  <- mantelhaen.test(cont_logrank, correct = T)$p.value
  #  res_mh2 <- mantelhaen.test(cont_logrank, correct = F)$p.value
  res_lr <- survdiff(Surv(c(x1test,x2test), c(x1delta,x2delta)) ~ c(rep(1,N1),rep(0,N2)))
  #  res_exact <- mantelhaen.test(cont_logrank,exact = T)$p.value
  # res_exact2 <- mantelhaen.test(cont_logrank,correct  = F)$p.value
  
  pval_mh[i,]<-c(res_lr$chisq, pchisq( res_lr$chisq, df=1, lower.tail=FALSE) )  
  #  pval_mh[i,]<-c( res_lr)
  
  setTxtProgressBar(pb, i) 
}


hist(pval_mh[,1], breaks=20, probability = T, ylim=c(0,1), xlim=c(0,max(pval_mh[,1])))
par(new=T)
curve(dchisq(x, df=1), col="red", ylim=c(0,1), xlim=c(0,max(pval_mh[,1])))


mean(pval_mh[,2]<0.05)

mean(pval_mh[,2]<0.05) - 2*sqrt(0.05*0.95/Nsim)
mean(pval_mh[,2]<0.05) + 2*sqrt(0.05*0.95/Nsim)

hist(pval_mh[,2], breaks=20, probability = T)
abline(v=mean(pval_mh[,2]<0.05), col="orange")
abline(v=0.05, col="blue")
abline(v=c(
  mean(pval_mh[,2]<0.05) - 2*sqrt(0.05*0.95/Nsim),
  mean(pval_mh[,2]<0.05) + 2*sqrt(0.05*0.95/Nsim)
), col="red"
)






pval_mh <- data.frame(pval_mh)
names(pval_mh) <- c("chisq", "p.value")

name2save <- paste("result_", seednum, sep = "")

dir.create(name2save)
write.csv(pval_mh, file = paste(name2save, "/pval.csv", sep = ""))

png(paste(name2save, "/p_value.png", sep = ""), width = 800, height = 600)
hist(pval_mh[,2], breaks=20, probability = T)
abline(v=mean(pval_mh[,2]<0.05), col="orange")
abline(v=0.05, col="blue")
abline(v=c(
  mean(pval_mh[,2]<0.05) - 2*sqrt(0.05*0.95/Nsim),
  mean(pval_mh[,2]<0.05) + 2*sqrt(0.05*0.95/Nsim)
), col="red"
)
dev.off()

#write(paste("simulation result: seed is ", seednum,"\n baseline distribution is geometric"), 
#      file=paste(name2save, "/detail.txt", sep = ""), append = F)
write(paste("simulation result: seed is ", seednum,"\n" ,
            "baseline distribution is", name_dist), 
      file=paste(name2save, "/detail.txt", sep = ""), append = F)


write(paste("N1, N2 = ", N1, ", ", N2), 
      file=paste(name2save, "/detail.txt", sep = ""), append = T)
write(paste("true, null = ", odds_true, ", ", odds_null), 
      file=paste(name2save, "/detail.txt", sep = ""), append = T)
write(paste("lambda1, lambda2 = ", lambda1, ", ", lambda2), 
      file=paste(name2save, "/detail.txt", sep = ""), append = T)
write(paste("Nsim = ", Nsim), 
      file=paste(name2save, "/detail.txt", sep = ""), append = T)
write(paste("p.value = ", mean(pval_mh[,2]<0.05)), 
      file=paste(name2save, "/detail.txt", sep = ""), append = T)
write(paste("C.I. of p.value = (", 
            mean(pval_mh[,2]<0.05) - 2*sqrt(0.05*0.95/Nsim) %>% round(.,digits = 3),", ",
            mean(pval_mh[,2]<0.05) + 2*sqrt(0.05*0.95/Nsim) %>% round(.,digits = 3), ")"),
      file=paste(name2save, "/detail.txt", sep = ""), append = T)



cat("random seed was ", seednum, "\n")