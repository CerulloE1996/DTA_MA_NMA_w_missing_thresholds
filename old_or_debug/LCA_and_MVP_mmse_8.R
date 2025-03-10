


require(rstan)
require(bayesplot)
require(dplyr)
require(shinystan)
require(devtools)
require(posterior)
require(ggplot2)

#devtools::install_github("stan-dev/cmdstanr", force = TRUE)

require(cmdstanr)
#install_cmdstan(overwrite = TRUE)


#path_to_opencl_lib <- "/mnt/c/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.0/lib/x64"
#write(paste0("LDFLAGS= -L\"",path_to_opencl_lib,"\" -lOpenCL"), file.path(cmdstan_path(), "make", "local"), append = TRUE)


###################################################################################################
setwd("/mnt/c/Users/Enzo/Documents/latent variable modelling/latent trait/MMSE")
setwd("/media/enzo/A05C06A35C0673F6/Users/Enzo/Documents/latent variable modelling/latent trait/MMSE")
setwd("C:/Users/Enzo/Documents/latent variable modelling/latent trait/MMSE")

rstan_options(auto_write = TRUE)            # Cache compiled Stan programs
options(mc.cores = parallel::detectCores()) # Parallelize chains
#Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
#Sys.setenv(LOCAL_CPPFLAGS = '-march=native -mtune=native')

options(scipen = 999)
options(max.print = 1000000000)

#############################################################################################
## MMSE data - re-structure the data into pts. falling between each threshold
#############################################################################################

################################################
## read in data
################################################
mmse_list <- readRDS("mmse_data2.rds")
#mmse_list

## indicator for the six primary care studies (3, 6, 7, 14, 15, 23)
primary_care_ind = c(0, 0, 1, 0, 0, 
                     1, 1, 0, 0, 0, 
                     0, 0, 0, 1, 1, 
                     0, 0, 0, 0, 0, 
                     0, 0, 1, 0, 0, 0)
primary_care_ind

## indicator for which ref. test used in each of the 27 studies
# 1 = DSM-III-R (10 studies)
# 2 = DSM-IV (7 studies)
# 3 = ICD-10 + DSM-III-R (3 studies)
# 4 = CDR (2 studies)
# 5 = ICD-10+DSM-IV (1)
# 6 = DSM-IV + NINCDS-ADRDA (1)
# 7 = DSM-III (2)

ref <- c(1, 1, 2, 3, 3, #5
         2, 2, 1, 1, 2, #10
         2, 4, 1, 4, 5, #15
         6, 1, 1, 1, 1, #20
         7, 1, 3, 2, 7, 2) #26
table(ref)
num_ref = length(unique(ref))
re <-  c(1, 1, 1, 1, 1, 
         1, 1, 1, 1, 1, 
         1, 1, 1, 1, 0, 
         0, 1, 1, 1, 1, 
         1, 1, 1, 1, 1, 1)


y = array(data= NA, dim = c( max(mmse_list$ns) , 2, mmse_list$NS))

str(y)

for (s in 1:mmse_list$NS) {
    y[, 1, s] <- mmse_list$y1[s,]
    y[, 2, s] <- mmse_list$y2[s,]
}

str(y)

y[1:1000, 2, 4]
y[1:1000, 1, 4]
(mmse_list$cutpoints2[studies, 1:16] -13 ) 

pa2 <- aperm(mmse_list$pa, c(2,1))
str(pa2)

twoway_thr = c(24,24,25,25,24,
               24,19,24,24, # 9
               19,24,23,24,23, # 14
               24,26,24,24,24,#19
               24,24,18,24,25,
               18, 25)

## make 2x2 table for each study - dichotomise each study at a given threshold

## mixed
r1 <- c(97, 1, 66, 21, 14, 77, 12, 93, 31, 42, 132, 9, 305, 19, 63, 34, 17, 44, 54, 39, 64, 7, 21, 36, 10, 24) # tp
r2 <- c(20, 0, 16, 2, 1, 0, 3, 30, 6, 4, 23, 5, 9, 9, 15, 0, 19, 5, 22, 4, 0, 2, 36, 33, 0, 0)  # fn
r3 <- c(14, 9, 23, 7, 44, 153, 20, 183,  3, 45, 54, 8, 98, 23, 78, 74, 31, 61, 63, 536, 612, 123, 16, 7, 23, 58)
r4 <- c(297, 45, 71, 226, 285, 130, 125, 1845, 248, 144, 226, 15, 256, 263, 147, 48, 580, 552, 980, 227, 2051, 368, 295, 166, 57, 32)

prevs <- c((r1+r2)/(r1+r2+r3+r4))

mmse_data <- read.csv("MMSE_data.csv")

mmse_data <- tibble(mmse_data)

mmse_data <- dplyr::rename(mmse_data,Study = ?..Study )
#View(mmse_data)


# label each study in order (1, ..., 26)
mmse_data_obs <-  filter(mmse_data, !(TP == 0 & FN == 0),  !(FP == 0 & TN == 0), Study != "ADAMS 2007", Study != "AMSTEL 1997 ", Study != "FillenBaum 1990", Study !="Frank 1996", 
                         Study != "Keskinoglu 2009", Study != "Lindesay 1997", Study != "Ramlall 2013" , Study !="Rummans 1996",
                         Study != 	"Scazufca 2009", Study != "Helsinki Aging Study 1994")
#View(mmse_data_obs)

mmse_data_obs2 <- mmse_data_obs %>% arrange(Study)

mmse_data_obs3 <-  mmse_data_obs2 %>% 
                  mutate(Threshold = Threshold - 13 , Study2 = as.numeric(as.factor(Study))) %>%
                  select( -N)
mmse_data_obs3

#View(mmse_data_obs3)

# make empty dataframe to merge with 

mmse_empty_frame <- tibble(Threshold = rep(seq(1, 16, 1),26), 
                           Study2 = rep(1:26, each = 16), 
                           TP = rep(999999, 26*16),
                           FP = rep(999999, 26*16),
                           TN = rep(999999, 26*16),
                           FN = rep(999999, 26*16))

#View(mmse_empty_frame)

mmse_data_obs4 <- full_join(mmse_empty_frame,mmse_data_obs3 , by = c("Threshold", "Study2")) %>% arrange(Study2, Threshold) %>% 
  mutate(TP.y = case_when(is.na(TP.y) ~ 999999,
         TRUE ~ as.numeric(as.character(.$TP.y))),
         FP.y = case_when(is.na(FP.y) ~ 999999,
         TRUE ~ as.numeric(as.character(.$FP.y ))),
         FN.y = case_when(is.na(FN.y) ~ 999999,    
         TRUE ~ as.numeric(as.character(.$FN.y ))),
         TN.y = case_when(is.na(TN.y) ~ 999999,
         TRUE ~ as.numeric(as.character(.$TN.y ))))

mmse_data_obs4 

#View(mmse_data_obs4)


r1 <- list(); r2 <- list(); 
r3 <- list(); r4 <- list(); 

for (i in 1:16) { 
  r1[[i]] <- filter(mmse_data_obs4, Threshold == i )$TP.y ; 
  r2[[i]] <- filter(mmse_data_obs4, Threshold == i )$FN.y ; 
  r3[[i]] <- filter(mmse_data_obs4, Threshold == i )$FP.y ; 
  r4[[i]] <- filter(mmse_data_obs4, Threshold == i )$TN.y ; 
}

r1

r = array(dim = c( 26,4,1,16)) ; r

for (s in 1:26) {
  for (thr in 1:16) { 
    r[s, 1, 1, thr] = r1[[thr]][s];
    r[s, 2, 1, thr] = r2[[thr]][s];
    r[s, 3, 1, thr] = r3[[thr]][s];
    r[s, 4, 1, thr] = r4[[thr]][s];
    }
}

r

#r = array(data = c(r1,r2,r3,r4), dim = c( 26,4,1)) ; r

prevs <- c((r1+r2)/(r1+r2+r3+r4))
sens <- c(r1/(r1+r2))
spec <- c(r4/(r3+r4))
studies25 <- c(0,0,0,4,5,0,0,8,9,0,11,0,0,0,0,0,0,0,19,20,0,0,0,24,0,26)
r2 = data.frame(r1,r2,r3,r4, twoway_thr,
                primary_care_ind, studies25, prevs, sens = round(sens,2) , spec = round(spec,2))
r2

xx <- rbeta(10000,10,10)*2-1 #
xx <- sort(xx)
xx[5000] ; xx[250] ; xx[9750]


xx <- rbeta(10000,1,6) # beta(1,5) has median 0.13 and CrI = (0.01, 0.52) 
                      # moderately informative for low-prev / community settings
xx <- sort(xx)
mean(xx)
xx[5000] ; xx[250] ; xx[9750]
xx[9999]
plot(density(xx))

xx <- pnorm(rnorm(10000,0.75,0.65)) #
xx <- pnorm(rnorm(10000,0,1)) #
xx <- pnorm(rnorm(10000,0.80,0.50)) #
xx <- pnorm(rnorm(10000,0.85,0.60)) #
mean(xx)
plot(density(xx))
xx <- sort(xx)
xx[5000] ; xx[250] ; xx[9750];
xx[9990];



# (-0.82, 0.82) for 2
# (-0.72, 0.72) for 3
# (-0.62, 0.62) for 4
# (-0.59, 0.58) for 5
# (-0.52, 0.52) for 6
# (-0.47, 0.47) for 8
# (-0.42, 0.42) for 10
# (-0.38, 0.38) for 12
# (-0.28, 0.28) for 24


qnorm(0.8)

pnorm(0.8416 + 1*0.3)
pnorm(0.8416 - 1*0.3)
pnorm(0.8416 + 1*0.6)
pnorm(0.8416 - 1*0.6)

pnorm(0.8416 + 2*0.3)
pnorm(0.8416 - 2*0.3)
pnorm(0.8416 + 2*0.6)
pnorm(0.8416 - 2*0.6)


mmse_list$ns

ns_cumsum <- cumsum(mmse_list$ns)
ns_cumsum

# indicator variable for thresholds

#############################################################################################################################
### all 26 studies  w/ full data
n <-   26
num <- 2727
studies <- c(1:26) 
max_cutpoint <- 16
length_cutpoint_vec <- 16


#n <-   17
#num <- 2727
#studies <- c(1,2,4,5,6,8,9,11,13,15,17,18,19,20,21,23,24)
#max_cutpoint <- 16
#length_cutpoint_vec <- 16


mmse_list$cutpoints2 - 13


data <- list( n_studies= n ,
             ns= mmse_list$ns[studies], 
             y=y[1:num, 1:2, studies], 
             n_binary_tests = 1, 
             nt = 2,
             pa = pa2[1:num, studies] ,
             n_patterns = 34, 
             ############################### threshold data
             n_thresholds_study = mmse_list$num_thresh_params[studies]-1     ,
             length_cutpoint_vec = length_cutpoint_vec,
             cutpoints = array(data = as.numeric(as.factor(mmse_list$cutpoints2[studies, 1:length_cutpoint_vec] -13 )) , 
                               dim = c(n,length_cutpoint_vec))   ,
             n_thresholds = length(unique(c(array(data = mmse_list$cutpoints2[studies, 1:length_cutpoint_vec] -13  ,
                                                      dim = c(n,length_cutpoint_vec))))) - 1     ,
             max_cutpoint = max_cutpoint ,
             twoway_thr = twoway_thr[studies]-13 ,
             ###############################
             numg = 2000, 
             p = prevs[studies],
             primary_care_ind = primary_care_ind[studies]     , 
            r = array(data = r[studies, , 1 , 1:16], dim = c(n, 4, 1, 16))   , 
            # r=r,
             num_refs = num_ref    ,
             ref =  ref[studies]  ,
             re = re[studies]   ,
             total_n = sum(mmse_list$ns[studies])  ,
             ind = c(0, rep(1, (length(studies[studies]) - 1)))   ,
             thresh = seq(-5, 5, length = 500)   ,
             ns_cumsum = c(0, cumsum(mmse_list$ns[studies])[1:25])       )



array(data = as.numeric(as.factor(mmse_list$cutpoints2[studies, 1:(length_cutpoint_vec+1)] -13 )) , 
      dim = c(n,(length_cutpoint_vec+1)))


m <- matrix(data = c(1, 0, 0, 0.80), nrow = 2, ncol = 2)
m <- matrix(data = c(1, 0, 0, 0.80), nrow = 2, ncol = 2)
m2 <- array(data = rep(m, n), dim = c(2,2,n))
m3 <- array(dim = c(n,2,2))

for (s in 1:n) {
  for (i in 1:2) { 
    for (j in 1:2) {
      m3[s,i,j] = m2[i,j,s]
    }
  }
}

m3





init = list(
            a2_m_raw = c(-1.5, 1.5), 
            a1_m_raw = array( data= c(-2,   -2,   -2,   -2,  -2,  -2,  -2, 
                                       1.5, 1.5,  1.5,  1.5, 1.5, 1.5 , 1.5), dim = c(7,2)),
            #b_primary_care = c(0.5, 0),
          # alpha = rep(1, 17),
         C_d = t(array(data = rep(c(seq(-1.7, 2.5, length = length_cutpoint_vec)), 16), dim = c(16, n))),
         C_nd = t(array(data = rep(c(seq(-1.7, 2.5, length = length_cutpoint_vec)), 16), dim = c(16, n))),
         p=prevs,
         kappa = 5,
         phi = c(0.22, 0.02, 0.03, 0.03, 0.04, 
                 0.04, 0.06, 0.05, 0.05, 0.04, 
                 0.05, 0.06, 0.06, 0.04, 0.07, 
                 0.07, 0.07),
    alpha_d = 0.1,
    alpha_nd = 0.1,
    L_Omega_global_d = m3[1,,],
    L_Omega_global_nd = m3[1,,],
    L_Omega_d = m3[ ,,],
    L_Omega_nd = m3[ ,,]
      )



#############################################################################################
# MCMC (HMC algorithm)
#############################################################################################
## indicator for which ref. test used in each of the 27 studies
# 1 = DSM-III-R (10 studies)
# 2 = DSM-IV (7 studies)
# 3 = ICD-10 + DSM-III-R (3 studies)
# 4 = CDR (2 studies)
# 5 = ICD-10+DSM-IV (1)
# 6 = DSM-IV + NINCDS-ADRDA (1)
# 7 = DSM-III (2)

ref <- c(1, 1, 2, 3, 3, #5
         2, 2, 1, 1, 2, #10
         2, 4, 1, 4, 5, #15
         6, 1, 1, 1, 1, #20
         7, 1, 3, 2, 7, 2) #

## conditional independence , perfect GS (for reference), random cutpoints
file <- file.path(file = "mvp_mmse_CI_random_cutpoints_perfectGS.stan") 
#file <- file.path(file = "mvp_mmse_CI_random_cutpoints_perfectGS_roughroc.stan") 

## conditional independence 
file <- file.path(file = "mvp_mmse_CI_random_cutpoints.stan") 
#file <- file.path(file = "mvp_mmse_CI_random_cutpoints_roughroc.stan") 
#file <- file.path(file = "mvp_mmse_CI_random_cutpoints_logit.stan") 

## conditional dependence (attempts)
#file <- file.path(file = "mvp_mmse_CD_random_cutpoints.stan") 
file <- file.path(file = "mvp_mmse_CD_random_cutpoints_logit2.stan")  # indep corr, lkj8
file <- file.path(file = "mvp_mmse_pCD_random_cutpoints_logit2.stan")  # pooled corr, lkj4

mod <- cmdstan_model(file)


meta_model2 <- mod$sample(
  data =  data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 500, 
  refresh = 10, 
  init = list(init, init,init,init),
  adapt_delta = 0.95,
  max_treedepth = 7)

#rstan_mod <- stan_model("mvp_mmse_CD_random_cutpoints_logit2.stan")



meta_model2r <- rstan::read_stan_csv(meta_model2$output_files())

print(meta_model2r, pars= c("SeR", "SpR", "SeI", "SpI", "SeIp", "SpIp", "SeI_pred","p"),  probs = c(0.025, 0.5,    0.975))
print(meta_model2r, pars= c( "sd1","sd2", "alpha", "L_Omega_bs1", "phi" , "p_dm", "kappa"),probs  = c(0.025, 0.5,    0.975))

#saveRDS(meta_model2r,  file = "mvp_mmse_pooledCD_covariate_randomthresh_kappa5_lkj4_sdr0.3_logit_td9.rds")

saveRDS(meta_model2r,  file = "mvp_mmse_nopoolingCD_covariate_randomthresh_kappa5_lkj8_sdr0.3_logit_td7.rds")

#meta_model2r <- readRDS("mvp_mmse_pooledCD_covariate_randomthresh_kappa5_lkj4_sdr0.3_logit.rds") 
#meta_model2r <- readRDS("mvp_mmse_nopoolingCD_covariate_randomthresh_kappa5_lkj8_sdr0.3_logit_td8.rds") 

stan_dens(meta_model2r, pars = c("p"), separate_chains = TRUE)
stan_dens(meta_model2r, pars = c("se"), separate_chains = TRUE)
stan_dens(meta_model2r, pars = c("sp"), separate_chains = TRUE)
stan_dens(meta_model2r, pars = c("L_Omega_d"), separate_chains = TRUE)
stan_dens(meta_model2r, pars = c("L_Omega_nd"), separate_chains = TRUE)

stan_trace(meta_model2r, pars = c("p"))
stan_trace(meta_model2r, pars = c("phi"))
stan_trace(meta_model2r, pars = c("se"))
stan_trace(meta_model2r, pars = c("sp"))
stan_trace(meta_model2r, pars = c("nu"))
stan_trace(meta_model2r, pars = c("z2"))
stan_trace(meta_model2r, pars = c("L_Omega_d"))
stan_trace(meta_model2r, pars = c("L_Omega_nd"))

#############################################################################################
# observed - expected correlation plots 
#############################################################################################
dc <- round(summary(meta_model2r, probs = c(0.025,  0.5, 0.975), pars = c("dc"))$summary[,5],3) ; dc
dc_l <- round(summary(meta_model2r, probs = c(0.025,  0.5, 0.975), pars = c("dc"))$summary[,4],3) ; dc_l
dc_u <- round(summary(meta_model2r, probs = c(0.025,  0.5, 0.975), pars = c("dc"))$summary[,6],3) ; dc_u

dc_data <- tibble(dc, dc_l, dc_u, Comparison = as.factor(c(rep("Ref vs MMSE", length(dc)))), 
                  obs = seq(1, length(dc), by = 1), 
                  Threshold =  rep(1:16, each = 26) )

#View(dc_data)
cutoff <- data.frame( x = c(-Inf, Inf), y = 0, cutoff = factor(0) )

#tiff("mvp_mmse_nopoolingCD_covariate_randomthresh_alpha5_logit.tif",units = "in", width = 30, height=16, res=500, compression = "lzw")
ggplot(data = dc_data, aes(y = dc, x=obs, colour = Comparison)) + geom_point(size = 1.5) + 
  geom_errorbar(aes(ymin=dc_l, ymax=dc_u), width= 1, position=position_dodge(.9)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  ylim(-0.4, 0.4) + 
  ylab("Correlation Residuals") + 
  #  scale_y_continuous(breaks = c(seq(-0.5, 0.5, by = 0.10))) + 
  theme(legend.position="none") +
  xlab(" ") + 
  theme(text = element_text(size=14),
        axis.text.x = element_text())  + facet_wrap( ~ Threshold, scales = "free")
dev.off()


dc_data2 <- dc_data %>%
  filter(dc != 999999) %>%
  mutate( ind = ifelse( dc_u < 0 | dc_l > 0, 1 , 0))

1 - sum(dc_data2$ind)/nrow(dc_data2) # fits 95% of data well.

dc_data2$ind

#################################################################################
#print(meta_model2r, pars= c("SeR", "SpR", "SeI", "SpI", "SeIp", "SpIp", "L_Omega_bs"))
print(meta_model2r, pars= c("SeR", "SpR", "SeI", "SpI", "SeIp", "SpIp", "b_primary_care"),  probs = c(0.025, 0.5,    0.975))
print(meta_model2r, pars= c( "dc"),probs = c(0.025, 0.5, 0.975),digits = 3)
print(meta_model2r, pars= c( "p"),probs  = c(0.025, 0.5,    0.975))
print(meta_model2r, pars= c( "se", "sp"),probs  = c(0.025, 0.5,    0.975))
print(meta_model2r, pars= c( "C_d"),probs  = c(0.025, 0.5,    0.975))
print(meta_model2r, pars= c( "sd1","sd2", "alpha"),probs  = c(0.025, 0.5,    0.975))
print(meta_model2r, pars= c( "L_Omega_global_d","L_Omega_global_nd", "rho_d"),probs  = c(0.025, 0.5,    0.975))
print(meta_model2r, pars= c( "L_Omega_d","L_Omega_nd"),probs  = c(0.025, 0.5,    0.975))
#################################################################################
###################
print(meta_model2r, pars= c( "rho_d", "rho_nd", "L_Omega_d", "L_Omega_nd"),probs  = c(0.025, 0.5,    0.975)) 
print(meta_model2r, pars= c( "se", "sp"),probs  = c(0.025, 0.5,    0.975)) 
print(meta_model2r, pars= c( "e"),probs = c(0.025, 0.5, 0.975),digits = 2)
print(meta_model2r, pars= c( "dt"),probs = c(0.025, 0.5, 0.975),digits = 2)
print(meta_model2r, pars= c( "L_Omega_d"),probs  = c(0.025, 0.5,    0.975))
print(meta_model2r, pars= c( "L_Omega_nd"),probs  = c(0.025, 0.5,    0.975))
print(meta_model2r, pars= c( "L_Omega_bs"),probs  = c(0.025, 0.5,    0.975))
print(meta_model2r, pars= c( "alpha_d"),probs  = c(0.025, 0.5, 0.975)) 

#########################
##### from specific chains
s<- summary(subset_draws(as_draws(meta_model2r),chain=c(2:4),
                         variable = c("SeR", "SpR", "SeI", "SpI","SeIp", "SpIp", "p","L_Omega_d","L_Omega_nd")), "mean",
            ~quantile(.x, probs = c(0.025, 0.975), na.rm = TRUE), "rhat")
View(s)


s2<-  summary(subset_draws(as_draws(meta_model2r),chain=c(2:4),
                     variable = c("a2_m_raw", "dc")), "mean", 
        ~quantile(.x, probs = c(0.025, 0.5,  0.975), na.rm = TRUE), "rhat")

View(s2)


#############################################################################################
# observed - expected table count plots 
#############################################################################################
dt <-   round(summary(meta_model2r, probs = c(0.025,  0.5, 0.975), pars = c("dt"))$summary[,5],3) ; dt
dt_l <- round(summary(meta_model2r, probs = c(0.025,  0.5, 0.975), pars = c("dt"))$summary[,4],3) ; dt_l
dt_u <- round(summary(meta_model2r, probs = c(0.025,  0.5, 0.975), pars = c("dt"))$summary[,6],3) ; dt_u

dt_data <- tibble(dt, dt_l, dt_u, Comparison = as.factor(c(rep("Ref vs MMSE", length(dt)))), 
                  obs = seq(1, length(dt), by = 1), 
                  Threshold = rep(c(1:16), 104),
                  Cell = factor(rep(rep(1:4, each = 16), 26), labels = c("++", "+-", "-+", "--")))

cutoff <- data.frame( x = c(-Inf, Inf), y = 0, cutoff = factor(0) )

#tiff("mvp_mmse_CI_covariate_random_thresh_phi_approx_tables.tif",units = "in", width = 8, height=6, res=500, compression = "lzw")
ggplot(data = dt_data, aes(y = dt, x=obs, colour = Cell)) + 
  geom_point(size = 1) + 
  geom_errorbar(aes(ymin=dt_l, ymax=dt_u), width= 0.75, size = 0.3, position=position_dodge(.9)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  ylim(-60, 60) + 
  ylab("Observed - expected table counts") + 
  xlab(" ") + 
  theme(text = element_text(size=14),
        #    axis.text.x = element_text(),
        axis.text.x=element_blank()) +
  # theme(legend.position="none") +
  facet_wrap( ~ Threshold, scales = "free")

#dev.off()
ind <- c()

dt_data2 <- dt_data %>%
  filter(dt != 999999) %>%
  mutate( ind = ifelse( dt_u < 0 | dt_l > 0, 1 , 0))

1 - sum(dt_data2$ind)/nrow(dt_data2) # fits 93% of data well.




#############################################################################################
# LOO
#############################################################################################
library(loo)
log_lik <- extract_log_lik(meta_model2r, parameter_name = "log_lik")
str(log_lik)

length <- length(log_lik[1,][log_lik[1,] != "NaN" & !is.na(log_lik[1,]) ])

draws <- 300
loglik2 <- array(data = log_lik[log_lik != "NaN" & !is.na(log_lik) ], dim = c(draws*4, length ))
str(loglik2)
sum(is.na(loglik2))


r_eff <- relative_eff(exp(loglik2), cores = 2, chain_id = c(rep(1, draws), rep(2, draws),
                                                            rep(3, draws), rep(4,draws)
                                                            ))

mod_loo4 <- loo(loglik2, r_eff = r_eff, cores = 2)
mod_loo4

loo_model <- loo(meta_model2r, cores = 1)
loo_model

loo3 <- loo::loo_moment_match.default(
  x = rstan_mod,
  loo = loo_model,
  post_draws = post_draws_stanfit,
  log_lik_i = log_lik_i_stanfit,
  unconstrain_pars = unconstrain_pars_stanfit,
  log_prob_upars = log_prob_upars_stanfit,
  log_lik_i_upars = log_lik_i_upars_stanfit,
  cores = 1
)
loo3

loo(meta_model2$draws("log_lik"))
loo(meta_model2$draws("log_lik"), r_eff=relative_eff(exp(meta_model2$draws("log_lik"))),chain_id = c(rep(1, draws), rep(2, draws),
                                                                                                rep(3, draws), rep(4,draws)) )

loo_compare(mod_loo4, mod_loo3)

# according to LOO, dependence model (M7) better than indep (M4)


#############################################################################################
#############################################################################################
# Dichotomous MMSE models
#############################################################################################
#############################################################################################

rr <- r

###############################################
#### 24 cut off (community settings - 14 studies)

# 5, 10 and 16 are primary care - exclude 

num <- 14
T1 <- c(1,2,2,2,2,3,2,3,2,2,2,2,1,3) ; length(T1)
T2 <- rep(4, times = 14)
T <- matrix(c(T1, T2), ncol=2, nrow=14)
r1 <- r[, 1, ,11][r[, 1, ,11] != 999999][-c(5, 10, 16)]
r2 <- r[, 2, ,11][r[, 2, ,11] != 999999][-c(5, 10, 16)]
r3 <- r[, 3, ,11][r[, 3, ,11] != 999999][-c(5, 10, 16)]
r4 <- r[, 4, ,11][r[, 4, ,11] != 999999][-c(5, 10, 16)]
ns <- c()
for (i in 1:num) {ns[i] <- r1[i] + r2[i] + r3[i] + r4[i]}
# order by test
data <- data.frame(r1,r2,r3,r4, ns, t1 = T[,1], t2= T[,2]) #%>% arrange(t1)
r1 <- data$r1 ; r2 <- data$r2 ;  r3 <- data$r3 ; r4 <- data$r4
r <- matrix(ncol = 4, nrow = num, c(r1,r2,r3,r4)) ; r
ns <- data$ns
data24 <-list()
pos <- r1+r2
neg <- r3+r4
data24 <- list( r = r, n = ns, NS= num , pos=pos, neg=neg, T=data$t1, num_ref=3, nt=2)
NS=14
sum(ns) 

prevs <- (r1+r2)/(r1+r2+r3+r4)

y_list <- list()
y1a <- list(length = max(ns))
y1b <- list(length = max(ns))
pa <-  list(length = max(ns))

max <- max(ns)

for (i in 1:NS) {
  y1a[[i]] = c(rep(1, r[i,1]), rep(1, r[i,2]), rep(0, r[i,3]), rep(0, r[i,4]), rep(100,  max - ns[i] )) # ref test
  y1b[[i]] = c(rep(1, r[i,1]), rep(0, r[i,2]), rep(1, r[i,3]), rep(0, r[i,4]), rep(100,  max - ns[i] )) # MMSE 
  pa[[i]] = c(rep(1, r[i,1]), rep(2, r[i,2]), rep(3, r[i,3]), rep(4, r[i,4]), rep(100,  max - ns[i] ) )
  y_list[[i]] =   matrix(ncol = 2, c(y1a[[i]] , y1b[[i]])) 
}

y = array(data= unlist(y_list), dim = c( max(ns), 2, NS))
pa2 = array(data = unlist(pa), dim = c( max(ns), 1, NS))



######################################################
#### 25 cut off (community settings - 9 studies)

# 1, 4 and  8 are primary care - exclude 

num <- 9
T1 <- c(1, 1, 2, 2, 3, 2, 2, 3, 3) ; length(T1)  # 2 studies ref 1, 4 studies ref 2, 3 studies ref 3. 
T2 <- rep(4, times = num)
T <- matrix(c(T1, T2), ncol=2, nrow=num)
r1 <- r[, 1, ,12][r[, 1, ,12] != 999999][-c(1, 4, 8)]
r2 <- r[, 2, ,12][r[, 2, ,12] != 999999][-c(1, 4, 8)]
r3 <- r[, 3, ,12][r[, 3, ,12] != 999999][-c(1, 4, 8)]
r4 <- r[, 4, ,12][r[, 4, ,12] != 999999][-c(1, 4, 8)]
ns <- c()
for (i in 1:num) {ns[i] <- r1[i] + r2[i] + r3[i] + r4[i]}
# order by test
data <- data.frame(r1,r2,r3,r4, ns, t1 = T[,1], t2= T[,2]) #%>% arrange(t1)
r1 <- data$r1 ; r2 <- data$r2 ;  r3 <- data$r3 ; r4 <- data$r4
r <- matrix(ncol = 4, nrow = num, c(r1,r2,r3,r4)) ; r
ns <- data$ns
data25 <-list()
pos <- r1+r2
neg <- r3+r4
data25 <- list( r = r, n = ns, NS= num , pos=pos, neg=neg, T=data$t1, num_ref=3, nt=2)
NS=num
sum(ns) 

prevs <- (r1+r2)/(r1+r2+r3+r4)

y_list <- list()
y1a <- list(length = max(ns))
y1b <- list(length = max(ns))
pa <-  list(length = max(ns))

max <- max(ns)

for (i in 1:NS) {
  y1a[[i]] = c(rep(1, r[i,1]), rep(1, r[i,2]), rep(0, r[i,3]), rep(0, r[i,4]), rep(100,  max - ns[i] )) # ref test
  y1b[[i]] = c(rep(1, r[i,1]), rep(0, r[i,2]), rep(1, r[i,3]), rep(0, r[i,4]), rep(100,  max - ns[i] )) # MMSE 
  pa[[i]] = c(rep(1, r[i,1]), rep(2, r[i,2]), rep(3, r[i,3]), rep(4, r[i,4]), rep(100,  max - ns[i] ) )
  y_list[[i]] =   matrix(ncol = 2, c(y1a[[i]] , y1b[[i]])) 
}

y = array(data= unlist(y_list), dim = c( max(ns), 2, NS))
pa2 = array(data = unlist(pa), dim = c( max(ns), 1, NS))


################################
#### run models
file <- file.path(file = "mvp_mmse_dichotomous_CI_perfectGS.stan") # dichotomous - M1
#file <- file.path(file = "mvp_mmse_dichotomous_CD.stan") # dichotomous - M2
file <- file.path(file = "mvp_mmse_dichotomous_CD_logit.stan") # dichotomous - M2

mod <- cmdstan_model(file)

#############################################################################################
# MCMC (HMC algorithm)
#############################################################################################

init = list(a2_m_raw = c(-1.5, 1.5), 
               a1_m_raw = array( data= c(-2, -2, -2,
                                         1,  1, 1), dim = c(3,2)),
            p = prevs,
            alpha_d = 0.1,
            alpha_nd = 0.1,
            L_Omega_global_d = m3[1,,],
            L_Omega_global_nd = m3[1,,],
            L_Omega_d = m3[1:n ,,],
            L_Omega_nd = m3[1:n ,,]
)

n <- NS
num <- max(ns[1:n])
data =  list( 
              NS = n, ns=ns[1:n], y=y[1:num,,1:n], num_binary_tests = 2,       
              nt = 2, r = array(r[1:n,], dim = c(n,4,1)), pa = pa2[1:num ,1,1:n] , numg = 2000, 
              n_patterns = 4,
              ns_cumsum = cumsum(ns[1:n])  ,
              total_n = sum(ns[1:n])  ,
              ind = c(0, rep(1, n-1))   , 
              prev = prevs,
              num_refs = 3,
              ref = T1,
              re = rep(1, NS))

meta_model2 <- mod$sample(
  data = data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1500, 
  refresh = 50, 
  init = list(init,init,init,init), 
  adapt_delta = 0.90, 
  max_treedepth = 7)

meta_model2r <- rstan::read_stan_csv(meta_model2$output_files())

print(meta_model2r, pars= c("SeR", "SpR", "SeI", "SpI","p","dc"),  probs = c(0.025, 0.5,    0.975))
print(meta_model2r, pars= c("L_Omega_d","L_Omega_nd"),  probs = c(0.025, 0.5,    0.975))

saveRDS(meta_model2r,  file = "mvp_mmse_dichotomous_CD_24_cutoff_logit.rds")

#meta_model2r <- readRDS("mvp_mmse_dichotomous_CI_perfectGS_24_cutoff_phi.rds") 

meta_model2r <- readRDS("mvp_mmse_dichotomous_CI_perfectGS_24_cutoff_logit.rds") 

#meta_model2r <- readRDS("mvp_mmse_dichotomous_CI_perfectGS_25_cutoff_phi.rds") 
meta_model2r <- readRDS("mvp_mmse_dichotomous_CD_24_cutoff_logit.rds") 

print(meta_model2r, pars= c("SeR", "SpR", "SeI", "SpI","dc","rho_d","rho_nd"),  probs = c(0.025, 0.5,    0.975))


#############################################################################################
# LOO
library(loo)
log_lik <- extract_log_lik(meta_model2r, parameter_name = "log_lik")
str(log_lik)

length <- length(log_lik[1,][log_lik[1,] != "NaN" & !is.na(log_lik[1,]) ])

draws <- 250
loglik2 <- array(data = log_lik[log_lik != "NaN" & !is.na(log_lik) ], dim = c(draws*4, length ))
str(loglik2)
sum(is.na(loglik2))


r_eff <- relative_eff(exp(loglik2), cores = 2, chain_id = c(rep(1, draws), rep(2, draws),
                                                            rep(3, draws), rep(4,draws)))

mod_loo <- loo(loglik2, r_eff = r_eff, cores = 2)
mod_loo


loo_compare(mod_loo, mod_loo3)

######
stan_trace(meta_model2r, pars = c("SeI", "SpI"))

s<- summary(subset_draws(as_draws(meta_model2r),chain=c(1),
                         variable = c("SeR", "SpR", "SeI", "SpI", "p","dc")), "mean",
            ~quantile(.x, probs = c(0.025, 0.50, 0.975), na.rm = TRUE), "rhat")
View(s)


########################################################################
#PLOTS
#########################################################################

m1 <- readRDS(file = "mvp_mmse_CI_covariate_randomthresh_kappa5_perfectGS_prevparams.rds")   ###  M1 - perfect GS
m2 <- readRDS(file = "mvp_mmse_CI_covariate_randomthresh_kappa5.rds")   ###  M2 - IGS, CI
m3 <- readRDS(file = "mvp_mmse_nopoolingCD_covariate_randomthresh_kappa5_lkj8_sdr0.3_logit_td7.rds")   ###  M3 - IGS, CD

print(m1, pars= c("SeR", "SpR", "SeI", "SpI", "SeIp", "SpIp"),  probs = c(0.025, 0.5,    0.975))

models <- list(m1, m2, m3)

dichot_24_GS <- readRDS(file = "mvp_mmse_dichotomous_CI_perfectGS_prevparams_24_cutoff_logit.rds")
dichot_24_CD <- readRDS(file = "mvp_mmse_dichotomous_CD_24_cutoff_logit.rds")

dichot_25_GS <- readRDS(file = "mvp_mmse_dichotomous_CI_perfectGS_prevparams_25_cutoff_logit.rds")
dichot_25_CD <- readRDS(file = "mvp_mmse_dichotomous_CD_25_cutoff_logit.rds")

models_dichot <- list(dichot_24_GS, dichot_24_CD, 
                      dichot_25_GS, dichot_25_CD)

print(dichot_24_CD, pars= c("SeR", "SpR","SeI", "SpI"),  probs = c(0.025, 0.5,    0.975))
print(dichot_25_CD, pars= c("SeR", "SpR","SeI", "SpI"),  probs = c(0.025, 0.5,    0.975))
#############################################################################################
# LOO TABLE 
#meta_model2r <- readRDS("Wells_fixed_thr_CI.rds") 

log_lik <- list()
log_lik2 <- list()
draws <- list()
length <- list()
r_eff <- list()
mod_loo <- list()

require(loo)
for (i in 1:length(models)) { 
  log_lik[[i]] <- extract_log_lik(models[[i]], parameter_name = "log_lik")
  
  draws[[i]] <- dim(log_lik[[i]])[1]/4
  length[[i]] <- length(log_lik[[i]][1,][log_lik[[i]][1,] != "NaN" & !is.na(log_lik[[i]][1,]) ])
 # log_lik2[[i]] <- array(data = log_lik[[i]][ log_lik[[i]] != "NaN" & !is.na(log_lik[[i]] ) ] , dim = c(draws[[i]]*4, length[[i]] ))
  
  
  r_eff[[i]] <- relative_eff(exp(log_lik[[i]]), cores = 1, chain_id = c(rep(1, draws[[i]]), rep(2, draws[[i]]), rep(3, draws[[i]]), rep(4,draws[[i]])))
  
  mod_loo[[i]] <- loo (log_lik[[i]], cores = 1, r_eff = r_eff[[i]])
} 

mod_loo

loo_compare(mod_loo)


loo_m <- loo(meta_model2$draws("log_lik"), r_eff=relative_eff(exp(meta_model2$draws("log_lik"))),chain_id = c(rep(1, draws), rep(2, draws),
                                                                                                     rep(3, draws), rep(4,draws)) )

loo_m

#############################################################################################
m_dichot <-c(); l_dichot  <-c(); u_dichot  <- c()
m2_dichot  <- c(); l2_dichot  <- c(); u2_dichot  <- c()
labels1_dichot  <- c(); labels2_dichot <- c()

meta_model_dichot2 <- m2
m_dichot  <- round(summary(meta_model_dichot2, probs = c(0.025,  0.5, 0.975), pars = c("SeI"))$summary[,5],2)  
m2_dichot  <- (1 -  round(summary(meta_model_dichot2, probs = c(0.025,  0.5, 0.975), pars = c("SpI"))$summary[,5],2))

l_dichot  <- round(summary(meta_model_dichot2, probs = c(0.025,  0.5, 0.975), pars = c("SeI"))$summary[,4],2)  
l2_dichot  <-  (1 -  round(summary(meta_model_dichot2, probs = c(0.025,  0.5, 0.975), pars = c("SpI"))$summary[,4],2))

u_dichot  <-  round(summary(meta_model_dichot2, probs = c(0.025,  0.5, 0.975), pars = c("SeI"))$summary[,6],2)  
u2_dichot   <- (1 -  round(summary(meta_model_dichot2, probs = c(0.025,  0.5, 0.975), pars = c("SpI"))$summary[,6],2))

l_dichot 
m_dichot 
u_dichot 

l2_dichot 
m2_dichot 
u2_dichot 

for (i in 1:1) { 
  labels1_dichot[i] <- paste0(m_dichot[i]*100, " [", l_dichot[i]*100, ",", u_dichot[i]*100, "]")
  labels2_dichot[i] <- paste0(m_dichot[i]*100, " [", round(l_dichot[i]*100, 3),  ",", u_dichot[i]*100, "]")
}

m_dichot  <- c(m_dichot , m2_dichot )
l_dichot  <- c(l_dichot , l2_dichot )
u_dichot  <- c(u_dichot , u2_dichot )

m_dichot  ; l_dichot ; m_dichot  

m_dichot2 <- rep(NA, 32)
l_dichot2 <- rep(NA, 32)
u_dichot2 <- rep(NA, 32)


round(summary(dichot_24_CD, pars= c("SeI", "SpI"),  probs = c(0.025, 0.5,    0.975))$summary[,5],2) 


## 24
m_dichot2[11] <- round(summary(dichot_24_CD, pars= c("SeI"),  probs = c(0.025, 0.5,    0.975))$summary[,5],2) 
m_dichot2[16 + 11] <- round(1 - summary(dichot_24_CD, pars= c("SpI"),  probs = c(0.025, 0.5,    0.975))$summary[,5],2) 
l_dichot2[11] <- round(summary(dichot_24_CD, pars= c("SeI"),  probs = c(0.025, 0.5,    0.975))$summary[,4],2) 
l_dichot2[16 + 11] <- round(1 - summary(dichot_24_CD, pars= c("SpI"),  probs = c(0.025, 0.5,    0.975))$summary[,4],2) 
u_dichot2[11] <- round(summary(dichot_24_CD, pars= c("SeI"),  probs = c(0.025, 0.5,    0.975))$summary[,6],2) 
u_dichot2[16 + 11] <- round(1 - summary(dichot_24_CD, pars= c("SpI"),  probs = c(0.025, 0.5,    0.975))$summary[,6],2) 

## 25
m_dichot2[12] <- round(summary(dichot_25_CD, pars= c("SeI"),  probs = c(0.025, 0.5,    0.975))$summary[,5],2) 
m_dichot2[16 + 12] <- round(1 - summary(dichot_25_CD, pars= c("SpI"),  probs = c(0.025, 0.5,    0.975))$summary[,5],2) 
l_dichot2[12] <- round(summary(dichot_25_CD, pars= c("SeI"),  probs = c(0.025, 0.5,    0.975))$summary[,4],2) 
l_dichot2[16 + 12] <- round(1 - summary(dichot_25_CD, pars= c("SpI"),  probs = c(0.025, 0.5,    0.975))$summary[,4],2) 
u_dichot2[12] <- round(summary(dichot_25_CD, pars= c("SeI"),  probs = c(0.025, 0.5,    0.975))$summary[,6],2) 
u_dichot2[16 + 12] <- round(1 - summary(dichot_25_CD, pars= c("SpI"),  probs = c(0.025, 0.5,    0.975))$summary[,6],2) 


#############################################################################################
# PLOT RESULTS
#############################################################################################
mm <- list() ; 
l<-list() ; 
u<- list() 
l_pred<-list() ; 
u_pred<- list() 
mm2 <- list() ; 
l2 <- list() ;
u2 <- list() 
l2_pred <- list() ;
u2_pred <- list() ;
mm_pc <- list() ; 
l_pc<-list() ; 
u_pc<- list() 
l_pred_pc<-list() ; 
u_pred_pc<- list() 
mm2_pc <- list() ; 
l2_pc <- list() ;
u2_pc <- list() 
l2_pred_pc <- list() ;
u2_pred_pc <- list() 

labels1 <- c(); labels2 <- c()


for (i in 1:length(models)) { 
  # community settings
  mm[[i]] <- summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("SeI"))$summary[,5]
  mm2[[i]] <- 1 -  summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("SpI"))$summary[,5]
  
  l[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("SeI"))$summary[,4],2)  
  l2[[i]] <-  (1 - round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("SpI"))$summary[,4],2))
  l_pred[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("SeI_pred2"))$summary[,4],2)  
  l2_pred[[i]] <-  (1 - round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("SpI_pred2"))$summary[,4],2))
  
  u[[i]] <-  round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("SeI"))$summary[,6],2)  
  u2[[i]]  <- (1 -  round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("SpI"))$summary[,6],2))
  u_pred[[i]] <-  round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("SeI_pred2"))$summary[,6],2)  
  u2_pred[[i]]  <- (1 -  round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("SpI_pred2"))$summary[,6],2))
  
  # primary care
  mm_pc[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("SeIp"))$summary[,5],5)
  mm2_pc[[i]] <- (1 -  round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("SpIp"))$summary[,5],5))
  
  l_pc[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("SeIp"))$summary[,4],2)  
  l2_pc[[i]] <-  (1 - round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("SpIp"))$summary[,4],2))
  l_pred_pc[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("SeIp_pred2"))$summary[,4],2)  
  l2_pred_pc[[i]] <-  (1 - round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("SpIp_pred2"))$summary[,4],2))
  
  u_pc[[i]] <-  round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("SeIp"))$summary[,6],2)  
  u2_pc[[i]]  <- (1 -  round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("SpIp"))$summary[,6],2))
  u_pred_pc[[i]] <-  round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("SeIp_pred2"))$summary[,6],2)  
  u2_pred_pc[[i]]  <- (1 -  round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("SpIp_pred2"))$summary[,6],2))
}

for (i in 1:length(models)) { 
  mm[[i]] <- c(mm[[i]], mm2[[i]])
  l[[i]] <- c(l[[i]], l2[[i]])
  u[[i]]<- c(u[[i]], u2[[i]])
  l_pred[[i]] <- c(l_pred[[i]], l2_pred[[i]])
  u_pred[[i]] <- c(u_pred[[i]], u2_pred[[i]])
  mm_pc[[i]] <- c(mm_pc[[i]], mm2_pc[[i]])
  l_pc[[i]] <- c(l_pc[[i]], l2_pc[[i]])
  u_pc[[i]] <- c(u_pc[[i]], u2_pc[[i]])
  l_pred_pc[[i]] <- c(l_pred_pc[[i]], l2_pred_pc[[i]])
  u_pred_pc[[i]] <- c(u_pred_pc[[i]], u2_pred_pc[[i]])
}


for (i in 1:16) { 
  labels1[i] <- paste0("M3: ", mm[[3]][i]*100, " [", l[[3]][i]*100, ",", u[[3]][i]*100, "]")
  labels2[i] <- paste0("M3: ", mm2[[3]][i]*100, " [", round(u2[[3]][i]*100, 3),  ",", l2[[3]][i]*100, "]")
}


require(rlist)
require(data.table)

labels <- c(labels1, labels2)
labels_2 <- rep(" ", length(labels))
labels_2[11:12] <- labels[11:12]
labels_2[27:28] <- labels[(27):(28)]


labels_dichot2 <- rep(NA, length(labels))
# Sens
labels_dichot2[11] <- paste0("Dichot: ", m_dichot2[11]*100, " [", l_dichot2[11]*100, ",", u_dichot2[11]*100, "]") # 24
labels_dichot2[12] <- paste0("Dichot: ", m_dichot2[12]*100, " [", l_dichot2[12]*100, ",", u_dichot2[12]*100, "]") # 25
# Spec
labels_dichot2[27] <- paste0("Dichot: ",m_dichot2[16 + 11]*100, " [", l_dichot2[16 + 11]*100, ",", u_dichot2[16 + 11]*100, "]") # 24
labels_dichot2[28] <- paste0("Dichot: ",m_dichot2[16 + 12]*100, " [", l_dichot2[16 + 12]*100, ",", u_dichot2[16 + 12]*100, "]") # 25

m_label <- rep(NA, length(labels))

m_label[11:12] <- unlist(mm[[3]][11:12])
m_label[27:28] <- unlist(mm[[3]][27:28])

Threshold = rep(c(14:29))
Measure <- c(rep("Sensitivity", 16), rep("1-Specificity", 16))
data_thresh <- data.frame(mm = c(mm[[1]][1:32], mm[[2]][1:32], mm[[3]][1:32],
                                 mm_pc[[1]][1:32], mm_pc[[2]][1:32], mm_pc[[3]][1:32]),
                          l  = c( l[[3]][1:32],  l[[3]][1:32],  l[[3]][1:32],
                                  l_pc[[3]][1:32],  l_pc[[3]][1:32],  l_pc[[3]][1:32]), # CrI's only for model 3
                          u  = c( u[[3]][1:32],  u[[3]][1:32],  u[[3]][1:32],
                                  u_pc[[3]][1:32],  u_pc[[3]][1:32],  u_pc[[3]][1:32]),
                          l_pred  = c( l_pred[[3]][1:32],  l_pred[[3]][1:32],  l_pred[[3]][1:32], 
                                       l_pred_pc[[3]][1:32],  l_pred_pc[[3]][1:32],  l_pred_pc[[3]][1:32]), # PrI's only for model 3
                          u_pred  = c( u_pred[[3]][1:32],  u_pred[[3]][1:32],  u_pred[[3]][1:32], 
                                       u_pred_pc[[3]][1:32],  u_pred_pc[[3]][1:32],  u_pred_pc[[3]][1:32]),
                          Model = factor(c(rep(1, 32), rep(2, 32), rep(3, 32)), 
                                         labels = c("M1: Perfect GS", "M2: Imperfect GS, CI", "M3: Imperfect GS, CD")),
                          Model_CrI = c(rep(NA, 32), rep(NA, 32), rep(3, 32)),
                          Threshold = Threshold,
                          Measure,
                          Setting = factor(c(rep(0, 96), rep(1, 96)), 
                                                labels = c("Community Settings (N=20)", "Primary Care (N=6)")),
                          m_label = c(unlist(m_label), rep(NA, 64)),
                          m_dichot2, 
                          l_dichot2, 
                          u_dichot2)
data_thresh


###############################################################################################
############
## Figure 1 : Community Settings for M1, M2, M3 + results from stratified dichot. models @ 24, 25 thresholds 
pd <- position_dodge(.2)

tiff("mmse_figure_1.tif",units = "in", width = 11, height=6, res=500, compression = "lzw")
require(ggrepel)
ggplot(data = filter(data_thresh, Setting == "Community Settings"), aes(x = Threshold, y = mm, col = Measure, linetype = Model)) + 
  geom_ribbon(aes(ymin=l, ymax=u, fill = Measure), position=pd,   alpha = 0.10, size = 0) + 
  geom_ribbon(aes(ymin=l_pred, ymax=u_pred, fill = Measure), position=pd,   alpha = 0.05, size = 0) + 
  geom_line(position=pd, alpha = 0.7) +
  geom_point(position=pd, size = 2, alpha = 0.7) + 
 geom_text_repel(aes(y=m_label, label = labels_2[1:96]), hjust = 1.0,  vjust = -0,
            size=4.00, label.size = 0.05, nudge_x = 0.5, nudge_y = 0.04, alpha = 0.6, col = "black",
            show.legend = FALSE) +
  geom_point(aes(y = m_dichot2, x = Threshold),  size = 3, shape = 1) + 
 geom_text_repel(aes(y=m_dichot2, label = labels_dichot2[1:96]), hjust = 1.0,  vjust = -0,
            size=4.00, label.size = 0.05,nudge_x = 1.5, nudge_y = -0.07,alpha = 0.6, col = "black",
            show.legend = FALSE) +
# geom_vline(xintercept=24, col = "green", linetype = 2, size = 2, alpha = 0.2) +
# geom_vline(xintercept=25, col = "green", linetype = 2, size = 2, alpha = 0.2) +
  theme_bw() + 
  xlab("Threshold") +
  ylab("Accuracy") + 
  scale_y_continuous(breaks = seq(0, 1, 0.10)) +
  scale_x_continuous(breaks = seq(14,29,1)) + 
  theme(legend.position="bottom", 
        legend.title = element_text(size = 11), 
        legend.text  = element_text(size = 11))

dev.off()

#  ggtitle("Diagnostic accuracy of MMSE in Community Settings at various cut-offs\n'Multiple Thresholds' Multivariate Probit (MVP) Model [Solid points] &\nStratified dichotomous MVP [Hollow points]")
# 95% posterior and prediction intervals for M3 

############
## Figure 2 : Primary care + Community Settings for M2 (top) and M3 (bottom)


g1<- ggplot(data = filter(data_thresh, Model == "M2: Imperfect GS, CI"), aes(x = Threshold, y = mm, col = Measure, linetype = Setting)) + 
      geom_ribbon(aes(ymin=l, ymax=u, linetype = Setting), position=pd,   alpha = 0.10, size = 0.25) + 
     # geom_ribbon(aes(ymin=l_pred, ymax=u_pred, fill = Measure), position=pd,   alpha = 0.05, size = 0) + 
      geom_line(position=pd, alpha = 0.7) +
      geom_point(position=pd, size = 2, alpha = 0.7) + 
      theme_bw() + 
      xlab("Threshold") +
      ylab("Accuracy") + 
      scale_y_continuous(breaks = seq(0, 1, 0.10)) +
      scale_x_continuous(breaks = seq(14,29,1)) + 
      theme(legend.position = "none") + 
     theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + 
  ggtitle("M2: Imperfect GS, CI")
      
  
g1

g2<- ggplot(data = filter(data_thresh, Model == "M3: Imperfect GS, CD"), aes(x = Threshold, y = mm, col = Measure, linetype = Setting)) + 
  geom_ribbon(aes(ymin=l, ymax=u, linetype = Setting), position=pd,   alpha = 0.10, size = 0.25) + 
  # geom_ribbon(aes(ymin=l_pred, ymax=u_pred, fill = Measure), position=pd,   alpha = 0.05, size = 0) + 
  geom_line(position=pd, alpha = 0.7) +
  geom_point(position=pd, size = 2, alpha = 0.7) + 
  theme_bw() + 
  xlab("Threshold") +
  ylab("Accuracy") + 
  scale_y_continuous(breaks = seq(0, 1, 0.10)) +
  scale_x_continuous(breaks = seq(14,29,1)) + 
  theme(legend.position="bottom") + 
  ggtitle("M3: Imperfect GS, CD")
g2

require(patchwork)
tiff("mmse_figure_2.tif",units = "in", width = 8, height=10, res=500, compression = "lzw")
g1+g2 + plot_layout(ncol = 1)
dev.off()
############
## Figure 3 : Reference test estimates for Imperfect GS models (M2, M3) 

# 1 = DSM-III-R (10 studies)
# 2 = DSM-IV (7 studies)
# 3 = ICD-10 + DSM-III-R (3 studies)
# 4 = CDR (2 studies)
# 5 = ICD-10+DSM-IV (1)
# 6 = DSM-IV + NINCDS-ADRDA (1)
# 7 = DSM-III (2)


m_refs <- list(); m2_refs <- list()
l_refs <- list(); l2_refs <- list()
u_refs <- list(); u2_refs <- list()
data_Se_mod <- list() 
data_Sp_mod <- list() 


for (i in 1:length(models)) {
  m_refs[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), 
                               pars = c("SeR"))$summary[,5],2)
  
  m2_refs[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), 
                                pars = c("SpR"))$summary[,5],2)
  
  
  l_refs[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), 
                               pars = c("SeR"))$summary[,4],2)  
  
  l2_refs[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), 
                                pars = c("SpR"))$summary[,4],2)
  
  
  u_refs[[i]] <-  round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), 
                                pars = c("SeR"))$summary[,6],2)  
  
  u2_refs[[i]]  <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), 
                                 pars = c("SpR"))$summary[,6],2)
  
  
  data_Se_mod[[i]] <- data.frame(m=m_refs[[i]],l=l_refs[[i]],u=u_refs[[i]], location =  c(0.9+0.5*(i-1), 
                                                                                          3.9+0.5*(i-1), 
                                                                                          6.9+0.5*(i-1), 
                                                                                          9.9+0.5*(i-1), 
                                                                                          12.9+0.5*(i-1),
                                                                                          15.9+0.5*(i-1),
                                                                                          18.9+0.5*(i-1)), 
                                 
                                 label  = factor(c("DSM-III-R (N=10)", 
                                                   "DSM-IV (N=7)",
                                                   "ICD-10 + DSM-III-R (N=3)", 
                                                   "CDR (N=2)",
                                                   "ICD-10+DSM-IV (N=1)",
                                                   "DSM-IV + NINCDS-ADRDA (N=1)",
                                                   "DSM-III (N=2)"), 
                                                 levels = c("DSM-III-R (N=10)", 
                                                            "DSM-IV (N=7)",
                                                            "ICD-10 + DSM-III-R (N=3)", 
                                                            "CDR (N=2)",
                                                            "ICD-10+DSM-IV (N=1)",
                                                            "DSM-IV + NINCDS-ADRDA (N=1)",
                                                            "DSM-III (N=2)")),
                                 Model = factor(rep(i,  7)))
  
  data_Sp_mod[[i]] <- data.frame(m=m2_refs[[i]],l=l2_refs[[i]],u=u2_refs[[i]], location =  c(0.9+0.5*(i-1), 
                                                                                             3.9+0.5*(i-1), 
                                                                                             6.9+0.5*(i-1), 
                                                                                             9.9+0.5*(i-1), 
                                                                                             12.9+0.5*(i-1),
                                                                                             15.9+0.5*(i-1),
                                                                                             18.9+0.5*(i-1)),
                                 label  = factor(c("DSM-III-R (N=10)", 
                                                   "DSM-IV (N=7)",
                                                   "ICD-10 + DSM-III-R (N=3)", 
                                                   "CDR (N=2)",
                                                   "ICD-10+DSM-IV (N=1)",
                                                   "DSM-IV + NINCDS-ADRDA (N=1)",
                                                   "DSM-III (N=2)"), 
                                                 levels = c("DSM-III-R (N=10)", 
                                                            "DSM-IV (N=7)",
                                                            "ICD-10 + DSM-III-R (N=3)", 
                                                            "CDR (N=2)",
                                                            "ICD-10+DSM-IV (N=1)",
                                                            "DSM-IV + NINCDS-ADRDA (N=1)",
                                                            "DSM-III (N=2)")),
                                 Model = factor(rep(i, 7)))
}

require(data.table)
# use rbindlist to put all models in same data frame
data_Se <- rbindlist(data_Se_mod)
data_Sp <- rbindlist(data_Sp_mod)


print(levels(data_Se$label))

data_Se$Model<- factor(data_Se$Model, labels = c("M1: Perfect GS",
                                                 "M2: Imperfect GS, CI", 
                                                 "M3: Imperfect GS, CD"))

data_Se$Mod = c(rep("M1: Perfect GS", 7),
                rep("M2: Imperfect GS, CI", 7),
                rep("M3: Imperfect GS, CD", 7))


data_Sp$Model<- factor(data_Sp$Model, labels = c("M1: Perfect GS",
                                                 "M2: Imperfect GS, CI", 
                                                 "M3: Imperfect GS, CD"))

data_Sp$Mod = c(rep("M1: Perfect GS", 7),
                rep("M2: Imperfect GS, CI", 7),
                rep("M3: Imperfect GS, CD", 7))


################# plot figure 3
se_plot <- ggplot(filter(tibble(data_Se), Model != "M1: Perfect GS"), aes(x=as.numeric(location), y = m, ymin= l,ymax= u, shape = Model)) + 
  geom_point(size=4, alpha=0.2) + 
  geom_pointrange(aes(colour=label)) +
  coord_flip() + 
  theme_bw()  + 
  xlab("") + 
  ylab("Sensitivity") + 
  theme(legend.position = "none") + 
  theme(axis.ticks = element_blank()) + 
  theme(text = element_text(size=16))+ 
  scale_y_continuous(breaks = seq(0, 1, 0.10), limits = c(0.2,1)) + 
  #scale_x_discrete( labels = rep("  ", length(data_Se$location))) + 
  geom_text(aes( label = paste0( m*100," ", "[",l*100,",", u*100 , "]") ,hjust = 1.1, vjust = -0.25), alpha = 0.35, size = 3.5) +
  scale_x_reverse() + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

se_plot


sp_plot <- ggplot(filter(tibble(data_Sp), Model != "M1: Perfect GS"), aes(x=location, y = m, ymin= l,ymax= u, shape = Model )) + 
  geom_point(size=4, alpha=0.2) + 
  geom_pointrange(aes(colour=(label))) +
  coord_flip() + 
  theme_bw()  + 
  xlab("") + 
  ylab("Specificity")  + 
  theme(legend.title = element_blank()) + 
  theme(text = element_text(size=16))+ 
  scale_y_continuous(breaks = seq(0, 1, 0.10), limits = c(0.4, 1)) + 
 # scale_x_discrete( labels = rep("  ", length(data_Se$location))) + 
  geom_text(aes( label = paste0( m*100," ", "[",l*100,",", u*100 , "]") ,hjust = 1.1, vjust = -0.25), alpha = 0.35, size = 3.5) +
  scale_x_reverse() + 
  theme(legend.title = element_text(size = 0), 
        legend.text  = element_text(size = 12),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
      #  legend.position="bottom"
        )

sp_plot

tiff("mmse_figure_3.tif",units = "in", width = 10, height=6, res=500, compression = "lzw")
se_plot + sp_plot
dev.off()
############
## Figure 4 : Posterior Predictive Check - correlation residual plot for M3 

dc <- round(summary(models[[3]], probs = c(0.025,  0.5, 0.975), pars = c("dc"))$summary[,5],3) ; dc
dc_l <- round(summary(models[[3]], probs = c(0.025,  0.5, 0.975), pars = c("dc"))$summary[,4],3) ; dc_l
dc_u <- round(summary(models[[3]], probs = c(0.025,  0.5, 0.975), pars = c("dc"))$summary[,6],3) ; dc_u

dc_data <- tibble(dc, dc_l, dc_u, Comparison = as.factor(c(rep("Ref vs MMSE", length(dc)))), 
                  obs = seq(1, length(dc), by = 1), 
                  Threshold =  rep(1:16, each = 26) )

#View(dc_data)
cutoff <- data.frame( x = c(-Inf, Inf), y = 0, cutoff = factor(0) )

tiff("mmse_figure_4_corr_resid.tif",units = "in", width = 15, height=10, res=500, compression = "lzw")

ggplot(data = dc_data, aes(y = dc, x=obs, colour = Comparison)) + geom_point(size = 3) + 
  geom_errorbar(aes(ymin=dc_l, ymax=dc_u), width= 0.75, position=position_dodge(.9)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  ylim(-0.4, 0.4) + 
  ylab("Correlation residuals") + 
  #  scale_y_continuous(breaks = c(seq(-0.5, 0.5, by = 0.10))) + 
  theme(legend.position="none") +
  xlab(" ") + 
  theme(text = element_text(size=14),
        axis.text.x = element_text())  + 
  facet_wrap( ~ Threshold, scales = "free")

dev.off()


############
## Figure 5 (appendix) : Posterior Predictive Check - 2x2 table residual plot for M3 

dt <-   round(summary(models[[3]], probs = c(0.025,  0.5, 0.975), pars = c("dt"))$summary[,5],3) ; dt
dt_l <- round(summary(models[[3]], probs = c(0.025,  0.5, 0.975), pars = c("dt"))$summary[,4],3) ; dt_l
dt_u <- round(summary(models[[3]], probs = c(0.025,  0.5, 0.975), pars = c("dt"))$summary[,6],3) ; dt_u

dt_data <- tibble(dt, dt_l, dt_u, Comparison = as.factor(c(rep("Ref vs MMSE", length(dt)))), 
                  obs = seq(1, length(dt), by = 1), 
                  Threshold = rep(c(1:16), 104),
                  Cell = factor(rep(rep(1:4, each = 16), 26), labels = c("++", "+-", "-+", "--")))

cutoff <- data.frame( x = c(-Inf, Inf), y = 0, cutoff = factor(0) )

tiff("mmse_figure_5_table_resid.tif",units = "in", width = 15, height=10, res=500, compression = "lzw")

ggplot(data = dt_data, aes(y = dt, x=obs, colour = Cell)) + 
  geom_point(size = 1) + 
  geom_errorbar(aes(ymin=dt_l, ymax=dt_u), width= 0.75, size = 0.3, position=position_dodge(.9)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  ylim(-60, 60) + 
  ylab("2x2 table residuals") + 
  xlab(" ") + 
  theme(text = element_text(size=14),
        #    axis.text.x = element_text(),
        axis.text.x=element_blank()) +
 # theme(legend.position="none") +
  facet_wrap( ~ Threshold, scales = "free")

dev.off()


##############################
####### ROC plot

########################################
# load in the model
mod <- m3


print(mod, pars= c("SeR", "SpR"),
      probs = c(0.025,0.5, 0.975))


require(bayesSurv)
require(scales)


pnorm2 <- function(x) { plogis( 1.702*x ) }

qnorm2 <- function(x) { (1/1.702)*qlogis(x) }

#####
## credible region
cred_1 <- list()
cred_1p <- list()

for (t in 1:7) { 
  cred_1[[t]] <- tibble(y = qnorm2(extract(mod, pars = "SeR")$Se[,t]), x = qnorm2(extract(mod, pars = "SpR")$Sp[,t]))
  cred_1p[[t]] <- tibble(y = (extract(mod, pars = "SeR")$Se[,t]), x =  (extract(mod, pars = "SpR")$Sp[,t]))
} 


require(data.table)
cred <- rbindlist(cred_1, idcol = TRUE)
cred_p <- rbindlist(cred_1p, idcol = TRUE)

cred2 <- mutate(cred,  Test = factor(.id,   label  = factor(c("DSM-III-R (N=10)", 
                                                              "DSM-IV (N=7)",
                                                              "ICD-10 + DSM-III-R (N=3)", 
                                                              "CDR (N=2)",
                                                              "ICD-10+DSM-IV (N=1)",
                                                              "DSM-IV + NINCDS-ADRDA (N=1)",
                                                              "DSM-III (N=2)"), 
                                                            levels = c("DSM-III-R (N=10)", 
                                                                       "DSM-IV (N=7)",
                                                                       "ICD-10 + DSM-III-R (N=3)", 
                                                                       "CDR (N=2)",
                                                                       "ICD-10+DSM-IV (N=1)",
                                                                       "DSM-IV + NINCDS-ADRDA (N=1)",
                                                                       "DSM-III (N=2)"))))


# in inv_probit space
g <- ggplot(data = cred2, aes(x = x, y = y, colour = Test))  + 
  # geom_point() + 
  stat_ellipse()  

g

# Get ellipse coordinates from plot
pb <-  ggplot_build(g)
View(pb$data[[2]])

el = pb$data[[1]][c("x","y", "group")]


credible_region <- tibble(x = pnorm(el$x), y = pnorm(el$y), Test = factor(el$group, label  = factor(c("DSM-III-R (N=10)", 
                                                                                                      "DSM-IV (N=7)",
                                                                                                      "ICD-10 + DSM-III-R (N=3)", 
                                                                                                      "CDR (N=2)",
                                                                                                      "ICD-10+DSM-IV (N=1)",
                                                                                                      "DSM-IV + NINCDS-ADRDA (N=1)",
                                                                                                      "DSM-III (N=2)"), 
                                                                                                    levels = c("DSM-III-R (N=10)", 
                                                                                                               "DSM-IV (N=7)",
                                                                                                               "ICD-10 + DSM-III-R (N=3)", 
                                                                                                               "CDR (N=2)",
                                                                                                               "ICD-10+DSM-IV (N=1)",
                                                                                                               "DSM-IV + NINCDS-ADRDA (N=1)",
                                                                                                               "DSM-III (N=2)"))))
credible_region


####
## prediction region

pred_1 <- list()
pred_1p <- list()

for (t in 1:7) { 
  pred_1[[t]] <- tibble(y = qnorm2(extract(mod, pars = "SeR_pred")$SeR_pred[,t]), x = qnorm2(extract(mod, pars = "SpR_pred")$SpR_pred[,t]))
  pred_1p[[t]] <- tibble(y = (extract(mod, pars = "SeR_pred")$SeR_pred[,t]), x =  (extract(mod, pars = "SpR_pred")$SpR_pred[,t]))
} 

require(data.table)
pred <- rbindlist(pred_1, idcol = TRUE)
pred_p <- rbindlist(pred_1p, idcol = TRUE)

pred2 <- mutate(pred,  Test = factor(.id, label  = factor(c("DSM-III-R (N=10)", 
                                                           "DSM-IV (N=7)",
                                                           "ICD-10 + DSM-III-R (N=3)", 
                                                           "CDR (N=2)",
                                                           "ICD-10+DSM-IV (N=1)",
                                                           "DSM-IV + NINCDS-ADRDA (N=1)",
                                                           "DSM-III (N=2)"), 
                                                         levels = c("DSM-III-R (N=10)", 
                                                                    "DSM-IV (N=7)",
                                                                    "ICD-10 + DSM-III-R (N=3)", 
                                                                    "CDR (N=2)",
                                                                    "ICD-10+DSM-IV (N=1)",
                                                                    "DSM-IV + NINCDS-ADRDA (N=1)",
                                                                    "DSM-III (N=2)"))))

# in inv_probit space
g <- ggplot(data = pred2, aes(x = x, y = y, colour = Test))  + 
  # geom_point() + 
  stat_ellipse()  

g

# Get ellipse coordinates from plot
pb <-  ggplot_build(g)
View(pb$data[[2]])

el = pb$data[[1]][c("x","y", "group")]


pred_region <- tibble(x = pnorm(el$x), y = pnorm(el$y), Test = factor(el$group,label  = factor(c("DSM-III-R (N=10)", 
                                                                                                 "DSM-IV (N=7)",
                                                                                                 "ICD-10 + DSM-III-R (N=3)", 
                                                                                                 "CDR (N=2)",
                                                                                                 "ICD-10+DSM-IV (N=1)",
                                                                                                 "DSM-IV + NINCDS-ADRDA (N=1)",
                                                                                                 "DSM-III (N=2)"), 
                                                                                               levels = c("DSM-III-R (N=10)", 
                                                                                                          "DSM-IV (N=7)",
                                                                                                          "ICD-10 + DSM-III-R (N=3)", 
                                                                                                          "CDR (N=2)",
                                                                                                          "ICD-10+DSM-IV (N=1)",
                                                                                                          "DSM-IV + NINCDS-ADRDA (N=1)",
                                                                                                          "DSM-III (N=2)"))))
pred_region


## medians
print(mod, pars= c("SeIp", "SpIp"),
      probs = c(0.025,0.5, 0.975))




median_sens <- c(round(summary(mod, probs = c(0.025,  0.5, 0.975), pars = c("SeR"))$summary[,5], 4))
median_spec <- c(round(summary(mod, probs = c(0.025,  0.5, 0.975), pars = c("SpR"))$summary[,5], 4))

medians <- tibble(median_sens = median_sens, median_spec = median_spec, Test = factor( c(1:7), label  = factor(c("DSM-III-R (N=10)", 
                                                                                                                 "DSM-IV (N=7)",
                                                                                                                 "ICD-10 + DSM-III-R (N=3)", 
                                                                                                                 "CDR (N=2)",
                                                                                                                 "ICD-10+DSM-IV (N=1)",
                                                                                                                 "DSM-IV + NINCDS-ADRDA (N=1)",
                                                                                                                 "DSM-III (N=2)"), 
                                                                                                               levels = c("DSM-III-R (N=10)", 
                                                                                                                          "DSM-IV (N=7)",
                                                                                                                          "ICD-10 + DSM-III-R (N=3)", 
                                                                                                                          "CDR (N=2)",
                                                                                                                          "ICD-10+DSM-IV (N=1)",
                                                                                                                          "DSM-IV + NINCDS-ADRDA (N=1)",
                                                                                                                          "DSM-III (N=2)"))))

print(mod, pars= c("SeR", "SpR"),   probs = c(0.025,0.5, 0.975))



#############################
## plot

g <- ggplot(data = medians, aes(y=median_sens, x = 1 - median_spec, colour = Test)) +    # summary points
  geom_point( size=2 ) + 
  #  geom_point(data = ss, aes(y=s_sens, x = 1 - s_spec, colour = Test))  + 
  # xlim(0, 1) + 
  #ylim(0, 1) + 
  theme_bw() + 
  scale_x_continuous(breaks = seq(0,1,0.1), limits = c(0,1))  + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1))  +
  theme(legend.title=element_blank()) + 
  xlab("1 - Specificity") + 
  ylab("Sensitivity") + 
  geom_polygon(data = credible_region, aes(x = 1  - x, y = y, colour = Test), alpha=0.05, size=0.4)  + 
  geom_path(data = pred_region, aes(x = 1  - x, y = y, colour = Test), linetype = 2, size=0.4) + 
  theme(legend.position =  "bottom")
g


tiff("mmse_sroc.tif",units = "in", width = 7, height=5, res=500, compression = "lzw")
g
dev.off()



data_thresh_m3 <- tibble(tp = filter(data_thresh, Model == "M3: Imperfect GS, CD", Measure == "Sensitivity")$mm, 
                         tp_u = filter(data_thresh, Model == "M3: Imperfect GS, CD", Measure == "Sensitivity")$u,
                         tp_l = filter(data_thresh, Model == "M3: Imperfect GS, CD", Measure == "Sensitivity")$l,
                         tp_up = filter(data_thresh, Model == "M3: Imperfect GS, CD", Measure == "Sensitivity")$u_pred,
                         tp_lp = filter(data_thresh, Model == "M3: Imperfect GS, CD", Measure == "Sensitivity")$l_pred,
                         fp = filter(data_thresh, Model == "M3: Imperfect GS, CD", Measure == "1-Specificity")$mm, 
                         fp_u = filter(data_thresh, Model == "M3: Imperfect GS, CD", Measure == "1-Specificity")$u,
                         fp_l = filter(data_thresh, Model == "M3: Imperfect GS, CD", Measure == "1-Specificity")$l,
                         fp_up = filter(data_thresh, Model == "M3: Imperfect GS, CD", Measure == "1-Specificity")$u_pred,
                         fp_lp = filter(data_thresh, Model == "M3: Imperfect GS, CD", Measure == "1-Specificity")$l_pred,
                         Setting = filter(data_thresh, Model == "M3: Imperfect GS, CD", Measure == "Sensitivity")$Setting,
                         Test = "MMSE")


tiff("mmse_sroc_2.tif",units = "in", width = 11, height= 7, res=500, compression = "lzw")

g + geom_path(data = data_thresh_m3, aes(y = tp, x = fp, linetype = Setting))+ 
  geom_point(data = data_thresh_m3, aes(y = tp, x = fp))
             
dev.off()
#geom_path(data = data_thresh_m3, aes(y = tp_u, x = fp_u) )







