library(runjags)
library(dplyr)
library(coda)
library(bayesplot)

options(scipen = 999)
options(max.print = 10000000)
runjags.options(force.summary = TRUE)
options(mc.cores = parallel::detectCores())


###########################################################
# input data  for mmse @ 25 and @25 thresholds 
# e.g. for burkatt: 
# 20 s.t. < 24, 1 at 24, 2 >=25, for diseased pop (TP+FN) and 2,5,226 resp. for non-diseased (FP+TN)
# for now just use the 7 studies with full data:

setwd("C:/Users/eceru/Documents/latent variable modelling/latent trait")

data <- read.csv("MMSE_data.csv")
head(data,50)
data1 <- list()  ; data2 <- list()
#View(data)
for (i in 1:length(unique(data$Threshold))){
data1[[i]] <- dplyr::filter(data, Threshold == i+13)
              }

for (i in 1:length(unique(data$Threshold))){
r1 <- (paste("r1", i, sep = "_"))
r2 <- (paste("r2", i, sep = "_"))
r3 <- (paste("r3", i, sep = "_"))
r4 <- (paste("r4", i, sep = "_"))
data2[[i]] <-  dplyr::rename(data1[[i]], r1 = TP, r2 = FP, r3 = FN, r4=TN)
}
data2[[11]]

for (i in 1:length(unique(data$Threshold))) {
r1[i]!!(paste("r1", i, sep = "_")) <-  data2[[i]]$r1 }

# for now just use the 7 studies with full data:
r1 <- c(20,14 , 31, 132, 39, 25) ; length(r1)
r2 <- c(3,   1    ,6,23,4,44)  ; length(r2)
###
r3 <- c(2,  44,      3,54,536, 3) ; length(r3)
r4 <- c(231, 285,   248, 226, 227,  170)  ; length(r4)
sens24 = r1/(r1+r2)
spec24 = r4/(r3+r4)
# for mmse @ 25
r5 <- c(21,15,34,141,39,36) ; length(r5)
r6 <- c(2,0,3,14,4,33) ; length(r6)
r7 <- c(7,53,10,76,567,7) ; length(r7)
r8 <- c(226,276,240,204,196,166) ; length(r8)

#convert dataset to format suitable for the multinomial likelihood: 
# diseased
c1 <- r1
c2 <- r5-r1
c3 <- (r1+r2) - (c1+c2) # r1+r2-r5
# non-diseased
c4 <- (r3+r4) - (r7 - r3+r8) #= r4-r7-r8
c5 <- r7 - r3
c6 <- r8


num <- 6
c <- matrix(ncol = 6, nrow = num, c(c1,c2,c3,c4,c5,c6)) ; c

ns <- c()
for (i in 1:num) {ns[i] <- c1[i] + c2[i] + c3[i] + c4[i] + c5[i] + c6[i]} #construct the data 
ns

data <- list()
data <- list( r = c, n = ns, m=5, NS=num)
data


res1<-list() ; res2 <- list()
n1 <-c(); n2<-c() ; n3 <- c(); n4<- c(); n5 <- c() ; n6 <- c() ; n7 <- c() ; n8 <- c()

# data for studies with thr. 1 
# thr 24, aevarsson
r1 = 97 ; r2 = 20 ; r3 = 14 ; r4 = 297
s=1
n1[s] = r1 ; n2[s] = r2 ; n3[s] = r3 ; n4[s] = r4
res1[[s]] <-  as.numeric(c(rep(2, n1[s]), rep(1,n2[s]), rep(2, n3[s]),  rep(1,n4[s])))  ; res1[s] #index test, 2 cats 
res2[[s]] <-  as.numeric(c(rep(2, n1[s]), rep(2,n2[s]), rep(1, n3[s]),  rep(1,n4[s])))  ; res2[s] # ref test, 2 cats 

# data for studies with thr. 2
# # thr 25, Brodaty
r1 = 66; r2 = 16; r3 = 23; r4 = 71
s=2
n1[s] = r1 ; n2[s] = r2 ; n3[s] = r3 ; n4[s] = r4
res1[[s]] <-  as.numeric(c(rep(2, n1[s]), rep(1,n2[s]), rep(2, n3[s]),  rep(1,n4[s])))  ; res1[s] #index test, 2 cats 
res2[[s]] <-  as.numeric(c(rep(2, n1[s]), rep(2,n2[s]), rep(1, n3[s]),  rep(1,n4[s])))  ; res2[s] # ref test, 2 cats 


# data for studies with 2 thresholds 
a <- c(1:6)
for (s in a) {
  n1[s] = data$r[s,1] ; n2[s] = data$r[s,2] ; n3[s] = data$r[s,3] ; n4[s] = data$r[s,4]
  n5[s] = data$r[s,5] ; n6[s] = data$r[s,6] ;
  res1[[s+2]] <-  as.numeric(c(rep(1, n1[s]), rep(2,n2[s]), rep(3, n3[s]),  rep(1,n4[s]), rep(2,n5[s]), rep(3, n6[s])))  ; res1[s] #index test, 3 cats 
  res2[[s+2]] <-  as.numeric(c(rep(2, n1[s]), rep(2,n2[s]), rep(2, n3[s]), rep(1,n4[s]), rep(1, n5[s]), rep(1, n6[s])))  ; res2[s] # ref test, 2 cats 
}
  
n = c(length(res1[[1]]), length(res1[[2]]), length(res1[[3]]), length(res1[[4]]), length(res1[[5]]), 
      length(res1[[6]]), length(res1[[7]]), length(res1[[8]]))
n

end.ids <-   c(n[1], sum(n[1:2]) , sum(n[1:3]) ,  sum(n[1:4]),  sum(n[1:5]),  sum(n[1:6]), sum(n[1:7]), sum(n[1:8])) ; end.ids
start.ids <- c(1,end.ids +1)[1:length(end.ids)] ; start.ids

library(plyr)
mats1 <- mapply(function(x, y) t(c(unlist(res1))[seq(x, y)]), start.ids, end.ids) 
mats2 <- mapply(function(x, y) t(c(unlist(res2))[seq(x, y)]), start.ids, end.ids) 
res1a <- t(rbind.fill.matrix(mats1))
res2a <- t(rbind.fill.matrix(mats2))


##### EEG data Xu et al 2013
res1 <- c(rep(1, 262), rep(2, 24), rep(3, 26))
res2 <- c(rep(1, 183), rep(2, 61), rep(3, 18), rep(1, 13), rep(2, 5), 
          rep(3, 6), rep(1, 3), rep(2, 4), rep(3, 19))
res3 <- c(rep(1, 181), rep(2, 1), rep(3, 1), rep(1, 56), rep(2, 5),
          rep(1,16), rep(2, 1), rep(3, 1), rep(1, 10), rep(3, 3), 
          rep(1, 5), rep(1, 2), rep(2, 2), rep(3, 2),
          rep(1,1), rep(2, 2), rep(3, 4), rep(2, 1), rep(3, 18))
length(res1)
res1

sum(res3==1)
sum(res3==2)
sum(res3 == 3)
#############
# models - use at least 200k iter and half as burn-in to get all MC%ofSD <5%. Much more for non-kappa models (>400k)
mod<-list()
files <- c("model_latenttrait_3tests_3cat_singlestudy.txt", 
           "model_latenttrait_3tests_3cat_singlestudy_cd_Qu.txt")

## single study
mod[[i]] <- run.jags(model = "model_latenttrait_3tests_3cat_singlestudy.txt", 
                     monitor = c("bb0", "bb1", "int_nd", "int_d","int_nd", "Se","Sp","SeI","SpI","b1", "c1","b0", "c0", "p",
                                 "a10", "a11", "a21", "a20",  "c11", "c10"), # , "dic", "pd"),
                    # data = list(res1=as.matrix(res1a), res2=as.matrix(res2a-1), n6=n , NS=8, N_1=1,N_2=1, N_1_2=6 ), 
                    data = list(res1=res1, res2=res2, res3=res3, n19 = 312), 
                   # data = list( n19 = 312), 
                     n.chains = 1,
                     #inits=inits,
                     burnin = 1000, adapt = 100, thin = 1 , sample = 1000,  keep.jags.files = FALSE)
mod[[i]]

         
#save(mod, file = "eeg_test.R")         
#load("eeg_test.R")         


c11 <- rnorm(1000, 0,0.444)
a11 <- rnorm(1000, 0, 0.444)
b1 <- rnorm(1000, 0, 0.444)
t <- rnorm(1000, 0, 1)
se <- 1/(1+exp(-1.7*(c11 - a11 + b1*t)))
plot(density(se))

