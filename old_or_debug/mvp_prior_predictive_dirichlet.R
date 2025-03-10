
###################################################################################################
setwd("/mnt/c/Users/Enzo/Documents/latent variable modelling/latent trait/MMSE")
setwd("/media/enzo/A05C06A35C0673F6/Users/Enzo/Documents/latent variable modelling/latent trait/MMSE")
setwd("C:/Users/Enzo/Documents/latent variable modelling/latent trait/MMSE")

require(MCMCpack)
require(ggplot2)
require(dplyr)
require(hexbin)
require(patchwork)
library(latex2exp)
require(truncnorm)


N <- 1e4
K <- 3

res <- array(NA_real_, dim = c(K, N))

kappa <- truncnorm::rtruncnorm(N, a = 0, mean =  0, sd = 50)

phi <- MCMCpack::rdirichlet(N, rep(1, K)) #Uniform simplex
for(i in 1:N) {
  res[,i] <- MCMCpack::rdirichlet(1, kappa[i]*phi[i,])
}

df_to_plot_1 <- data.frame(p1 = res[1,], p2 = res[2,], dist = 1) %>% filter(!is.na(p1), !is.na(p2))     


g1 <- ggplot(df_to_plot_1, aes(p1, p2)) + 
  geom_hex() + 
  scale_fill_continuous(trans = "log10") +
  theme_bw() + 
  # ggtitle(TeX("$ \\kappa  \\sim     N(0, 5) $")) + 
  xlab(TeX("$P_{i}$")) + 
  ylab(TeX("$P_{j}$")) + 
  xlim(0,1) + 
  ylim(0,1)# + 
#facet_wrap( ~ dist,   labeller=label_parsed)

g1



kappa <- rnorm(N, 5, 2.5)
phi <- MCMCpack::rdirichlet(N, rep(1, K)) #Uniform simplex
for(i in 1:N) {
  res[,i] <- MCMCpack::rdirichlet(1, exp(kappa[i]) * phi[i,])
}
df_to_plot_2 <- data.frame(p1 = res[1,], p2 = res[2,], dist = 2) %>% filter(!is.na(p1), !is.na(p2))   

df_to_plot = rbind(df_to_plot_1, df_to_plot_2)

df_to_plot$dist <- factor(df_to_plot$dist, labels = 
                                          c(TeX("$ \\kappa  \\sim     N(0, 5) $"),
                                            TeX("$ \\kappa  \\sim     N(5, 2.5) $")))


tiff("prior_pred_check_mmse.tif",units = "in", width = 6, height=5, res=500, compression = "lzw")

g1 

dev.off()



g1 <- ggplot(df_to_plot, aes(p1, p2)) + 
  geom_hex() + 
  scale_fill_continuous(trans = "log10") +
  theme_bw() + 
  # ggtitle(TeX("$ \\kappa  \\sim     N(0, 5) $")) + 
  xlab(TeX("$P_{i}$")) + 
  ylab(TeX("$P_{j}$")) + 
  xlim(0,1) + 
  ylim(0,1) + 
  facet_wrap( ~ dist,   labeller=label_parsed)

g1


tiff("prior_pred_check_mmse.tif",units = "in", width = 10, height=5, res=500, compression = "lzw")

g1 

dev.off()


#################################
# Wells (3 categories) 

N <- 1e6
K <- 3

res <- array(NA_real_, dim = c(K, N))

kappa <- rnorm(N, 0, 5)
phi <- MCMCpack::rdirichlet(N, rep(1, K)) #Uniform simplex
for(i in 1:N) {
  res[,i] <- MCMCpack::rdirichlet(1, exp(kappa[i]) * phi[i,])
}

df_to_plot_1 <- data.frame(p1 = res[1,], p2 = res[2,], dist = 1) %>% filter(!is.na(p1), !is.na(p2))     

kappa <- rnorm(N, 5, 2.5)
phi <- MCMCpack::rdirichlet(N, rep(1, K)) #Uniform simplex
for(i in 1:N) {
  res[,i] <- MCMCpack::rdirichlet(1, exp(kappa[i]) * phi[i,])
}
df_to_plot_2 <- data.frame(p1 = res[1,], p2 = res[2,], dist = 2) %>% filter(!is.na(p1), !is.na(p2))   

df_to_plot = rbind(df_to_plot_1, df_to_plot_2)

df_to_plot$dist <- factor(df_to_plot$dist, labels = 
                            c(TeX("$ \\kappa  \\sim     N(0, 5) $"),
                              TeX("$ \\kappa  \\sim     N(5, 2.5) $")))



g2 <- ggplot(df_to_plot_1, aes(p1, p2)) + 
  geom_hex() + 
  scale_fill_continuous(trans = "log10") +
  theme_bw() + 
  # ggtitle(TeX("$ \\kappa  \\sim     N(0, 5) $")) + 
  xlab(TeX("$P_{i}$")) + 
  ylab(TeX("$P_{j}$")) + 
  xlim(0,1) + 
  ylim(0,1)# + 
#facet_wrap( ~ dist,   labeller=label_parsed)

g2


tiff("prior_pred_check_wells.tif",units = "in", width = 6, height=5, res=500, compression = "lzw")

g2

dev.off()



g2 <- ggplot(df_to_plot, aes(p1, p2)) + 
  geom_hex() + 
  scale_fill_continuous(trans = "log10") +
  theme_bw() + 
 # ggtitle(TeX("$ \\kappa  \\sim     N(0, 5) $")) + 
  xlab(TeX("$P_{i}$")) + 
  ylab(TeX("$P_{j}$")) + 
  xlim(0,1) + 
  ylim(0,1) + 
  facet_wrap( ~ dist,   labeller=label_parsed)

g2


tiff("prior_pred_check_wells.tif",units = "in", width = 10, height=5, res=500, compression = "lzw")

g2

dev.off()









