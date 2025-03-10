

require(rstan)
require(bayesplot)
require(dplyr)
require(shinystan)
require(devtools)
require(posterior)
require(patchwork)
require(loo)


#devtools::install_github("stan-dev/cmdstanr", force = TRUE)

require(cmdstanr)
#install_cmdstan(overwrite = TRUE)


#path_to_opencl_lib <- "/mnt/c/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.0/lib/x64"
#write(paste0("LDFLAGS= -L\"",path_to_opencl_lib,"\" -lOpenCL"), file.path(cmdstan_path(), "make", "local"), append = TRUE)


###################################################################################################
setwd("/mnt/c/Users/Enzo/Documents/latent variable modelling/latent trait/Wells")
setwd("/media/enzo/A05C06A35C0673F6/Users/Enzo/Documents/latent variable modelling/latent trait/Wells")
setwd("C:/Users/Enzo/Documents/latent variable modelling/latent trait/Wells")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
#Sys.setenv(LOCAL_CPPFLAGS = '-march=native -mtune=native')

options(scipen = 999)
options(max.print = 1000000000)


################################################################################################
# Wells, D-Dimer and ref test(s) data (11 studies w/ complete data)
num <- 11
#T1 <- c(1,1,2,2,1,3,2,2,1,2,4,2,2) 
T1 <- rep(1, times = 11) # need to see what refs the studies used 
T2 <- rep(2, times = 11)
T <- matrix(c(T1, T2), ncol=2, nrow=11)
### GS=1
#DDIMER=1
r1 <- c(1, 3, 4, 15, 1, 4, 17, 18, 23, 6, 28) # low
r2 <- c(6, 9, 17, 31,22,6, 61, 16, 37, 15, 117) # moderate
r3 <- c(8, 30, 33, 54,35,13,79,21, 9,  10, 109) # high
# 2nd dichot (default ; Novielli et al)
r1_b <- r1 # wells = 0 
r2_b <- r2+r3 #  wells = 1 , GS = 1, DDIMER = 1
# 1st dichot 
#r1_b <- r1 + r2 # wells = 0 
#r2_b <- r3 #  wells = 1 , GS = 1, DDIMER = 1
#DDIMER=0
r4 <- c(0, 1, 1, 1, 0, 0, 3, 0, 5, 0, 1)
r5 <- c(0, 3, 7, 0, 0, 3, 15,1, 7, 3, 0)
r6 <- c(2, 0, 2, 1, 0, 2, 15, 0, 0, 1, 0)
r3_b <- r4
r4_b <- r5+r6 #  wells = 1 , GS = 1, DDIMER = 0
#r3_b <- r4 + r5 # wells = 0 
#r4_b <- r6 #  wells = 1 , GS = 1, DDIMER = 1
### GS=0
#DDIMER=1
r7 <- c(8, 8, 25, 49, 20, 17, 113, 85, 1, 42, 233)
r8 <- c(18, 12,51,51, 23, 9, 93,  83,  6, 59, 104)
r9 <- c(2, 8, 8, 36,  9,  2, 55,  30, 3, 16,  29)
r5_b <- r7
r6_b <- r8+r9 #  wells = 1 , GS = 0, DDIMER = 1
#r5_b <- r7 + r8 # wells = 0 
#r6_b <- r9 #  wells = 1 , GS = 1, DDIMER = 1
 #DDIMER=0
r10 <- c(32, 76, 176, 70, 17, 97, 313, 193, 3, 41, 243)
r11 <- c(20, 43, 113, 54, 19, 48, 23,  89, 5, 46, 16)
r12 <- c(5,  7,  6,  21,  12, 13, 50, 20,  2, 4, 3)
r7_b <- r10
r8_b <- r11+r12 #  wells = 1 , GS = 0, DDIMER = 0
#r7_b <- r10 + r11 # wells = 0 
#r8_b <- r12 #  wells = 1 , GS = 1, DDIMER = 1

#marginalise over d-dimer (ref vs wells) 
r1_t1 <- r1_b + r3_b ; r1_t1
r2_t1 <- r2_b + r4_b ; r2_t1
r3_t1 <- r5_b + r7_b ; r3_t1
r4_t1 <- r6_b + r8_b ; r4_t1

#marginalise over Wells (ref vs d-dimer) -
r1_t2 <- r1_b + r2_b
r2_t2 <- r3_b + r4_b
r3_t2 <- r5_b + r6_b
r4_t2 <- r7_b + r8_b

#marginalise over ref (wells vs d-dimer) -
r1_t3 <- r1_b + r5_b
r2_t3 <- r2_b + r6_b
r3_t3 <- r3_b + r7_b
r4_t3 <- r4_b + r8_b

prev <- (r1_t1+r2_t1)/(r2_t1+ r1_t1 + r3_t1 + r4_t1)
prevs <- round(prev,2)

ns <- c()
for (i in 1:num) {ns[i] <- r1[i] + r2[i] + r3[i] + r4[i] +
                           r5[i] + r6[i] + r7[i] + r8[i] +
                           r9[i] + r10[i] + r11[i] + r12[i] }
# order by test
data <- data.frame(r1,r2,r3,r4,r5, r6, r7, r8, r9, r10, r11, r12, ns, t1 = T[,1], t2= T[,2]) #%>% arrange(t1)

r <- matrix(ncol = 12, nrow = num, c(r1,r2,r3,r4,r5, r6, r7, r8, r9, r10, r11, r12)) ; r

t1 <- matrix(ncol = 4, nrow = num, c(r2_t1, r1_t1, r4_t1, r3_t1)) ;t1 ; sum(t1) # this is fine 

t2 <- matrix(ncol = 4, nrow = num, c(r1_t2, r2_t2, r3_t2, r4_t2)) # this is fine 

t3 <- matrix(ncol = 4, nrow = num, c(r2_t3, r1_t3, r4_t3, r3_t3)) ; sum(t3)
ns <- data$ns
data24 <-list()
data24 <- list( r = r, n = ns, NS= num , T=data$t1, num_ref=4, nt=2)
NS <- 11
sum(ns) # N= 4120

y_list <- list()
y1a <- list(length = max(ns))
y1b <- list(length = max(ns))
y1c <- list(length = max(ns))
pa <-  list(length = max(ns))

max <- max(ns)

for (i in 1:NS) {
  y1a[[i]] = c(rep(1, r[i,1]), rep(1, r[i,2]), rep(1, r[i,3]), rep(1, r[i,4]), rep(1, r[i,5]), rep(1, r[i,6]), rep(0, r[i,7]), rep(0, r[i,8]), rep(0, r[i,9]), rep(0, r[i,10]), rep(0, r[i,11]), rep(0, r[i,12]), rep(100,  max - ns[i] )) # ref test
  y1b[[i]] = c(rep(1, r[i,1]), rep(1, r[i,2]), rep(1, r[i,3]), rep(0, r[i,4]), rep(0, r[i,5]), rep(0, r[i,6]), rep(1, r[i,7]), rep(1, r[i,8]), rep(1, r[i,9]), rep(0, r[i,10]), rep(0, r[i,11]), rep(0, r[i,12]), rep(100,  max - ns[i] )) # D-Dimer
  y1c[[i]] = c(rep(1, r[i,1]), rep(2, r[i,2]), rep(3, r[i,3]), rep(1, r[i,4]), rep(2, r[i,5]), rep(3, r[i,6]), rep(1, r[i,7]), rep(2, r[i,8]), rep(3, r[i,9]), rep(1, r[i,10]), rep(2, r[i,11]), rep(3, r[i,12]), rep(100,  max - ns[i] ) ) # Wells
  pa[[i]] = c(rep(1, r[i,1]), rep(2, r[i,2]), rep(3, r[i,3]), rep(4, r[i,4]), rep(5, r[i,5]), rep(6, r[i,6]), rep(7, r[i,7]), rep(8, r[i,8]), rep(9, r[i,9]), rep(10, r[i,10]), rep(11, r[i,11]), rep(12, r[i,12]), rep(100,  max - ns[i] ) ) #non-dichotomous wells 
   y_list[[i]] =   matrix(ncol = 3, c(y1a[[i]] , y1b[[i]], y1c[[i]] )) 
}

y = array(data= unlist(y_list), dim = c( max(ns), 3, NS))
pa2 = array(data = unlist(pa), dim = c( max(ns), 1, NS))
r = array(data = c(t2, t1, t3), dim = c( NS,4,3)) ; r

r_full <- array(data= c(r1,r2,r3,r4,r5, r6, r7, r8, r9, r10, r11, r12), dim = c(11, 12))
r_full
######################################################################################################

# 2 S.D from mean of 0.6 if normal(0,0.5) prior for logistic or N(0, 0.3) prior for probiot
pnorm(0.6 + 0.6) ; pnorm(0.6 - 0.6) ; pnorm(0.6 ) ; 
plogis(1 + 1) ; plogis(1 - 1) ; plogis(1 ) ;
# 2 S.D from mean of 0.6 if normal(0,1) prior for logistic or normal(0, 0.6) prior for probit
pnorm(0.6 + 1.2) ; pnorm(0.6 - 1.2) ; pnorm(0.6 ) ; 
plogis(1 + 2) ; plogis(1 - 2) ; plogis(1 ) ; 
## piors for means on probit scale
prior <- pnorm(rnorm(10000, 0.75, 0.40)) # 49 - 94
#prior <- pnorm(rnorm(10000, 1.7, 0.40)) # 82 - 99
plot(density(prior))
sort(prior)[250] ; sort(prior)[5000] ; sort(prior)[9750]

prior <- (rbeta(10000, 1,5))
plot(density(prior))
sort(prior)[250] ; sort(prior)[5000] ; sort(prior)[9750]
######################################################################################################

n = 11
m <- matrix(data = c(1, 0, 0, 0.80), nrow = 2, ncol = 2)
m2 <- array(data = rep(m,n), dim = c(2,2,n))
m3 <- array(dim = c(n,2,2))

for (s in 1:n) {
  for (i in 1:2) { 
    for (j in 1:2) {
      m3[s,i,j] = m2[i,j,s]
    }
  }
}
init = list(
  mu = t(array(dim = c(3, 2), data = c(  0.85, 1.60, 0.80,
                                        -2, -0.35, -0.82 ) )) , 
            a1_m_raw = c(-2, 0.85), 
            a2_m_raw = c(-0.35, 1.60),
            a3_m_raw = c(-0.82, 0.80),
        #    C_dm= c(-0.5, 0.5), 
            C_d = t(array(dim = c(2,11), data = rep(c(-0.5, 0.5),22 ) ))  ,
            alpha_d = 0.49,
            alpha_nd = 0.51, 
            alpha = c(7,7,7),
            p = prevs,
            L_Omega_global_d = m3[1,,],
            L_Omega_global_nd = m3[1,,],
            mu_s = array(dim = c(2,3), data = c(1,1,1,1,1,1)))


######################################################################################################
################
# dichotomous Wells  ( imperfect GS )
#file <- file.path(file = "mvp_wells_CI_dichotomous_Wells.stan") # CI
#file <- file.path(file = "mvp_wells_CD_dichotomous_Wells.stan") # CD

#################
# polytomous Wells 

# perfect GS  
file <- file.path(file = "mvp_wells_CI_random_cutpoints_perfectGS_2.stan") # CI, random thr
#file <- file.path(file = "mvp_wells_CI_random_cutpoints_perfectGS_logit_2.stan") # CI, random thr
#file <- file.path(file = "mvp_wells_CD_random_cutpoints_perfectGS_2.stan") # CD, random thr
file <- file.path(file = "mvp_wells_CD_random_cutpoints_perfectGS_logit_2.stan") # CD, random thr

# imperfect GS
#file <- file.path(file = "mvp_wells_CI_random_cutpoints.stan") # CI, random thr
#file <- file.path(file = "mvp_wells_CI_random_cutpoints_logit.stan") # CI, random thr

#file <- file.path(file = "mvp_wells_CD_random_cutpoints2.stan") # CD, random thr - homog sd's (hetero N.I)
file <- file.path(file = "mvp_wells_CD_random_cutpoints_logit_2.stan") # CD, random thr - homog sd's (hetero N.I)


mod <- cmdstan_model(file)
#r_mod <- stan_model( file)

#############################################################################################
# MCMC (HMC algorithm)
#############################################################################################
n <- 11
num <- max(ns[1:n])
data =  list( x = 2, y3 = 0, y4 = 0, 
              n_studies=n, NS = n, ns=ns[1:n], y=y[1:num,,1:n], num_binary_tests = 2, Thr = c(1,1,2),         
              nt = 3, r = array(r[1:n, , ], dim = c(n,4,3)), pa = pa2[1:num ,1,1:n] , numg = 2000, 
              n_patterns = 12,
              ns_cumsum = c(0, cumsum(ns[1:n-1])),
              total_n = sum(ns[1:n]),
              r_full = r_full[1:n,], 
              ind = c(0, rep(1, n-1)), 
              n_thresholds = 2,
              p = prevs,
              ind2 = c(0, 0, 1))

meta_model2 <- mod$sample(
  data = data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1500, 
  refresh = 10, 
  init = list(init,init,init,init), 
  adapt_delta = 0.95, 
  max_treedepth = 7)

#.rs.restartR()

meta_model2r <- rstan::read_stan_csv(meta_model2$output_files())


print(meta_model2r, pars= c("Se", "Sp", "Wells_DDimer_BTN_Se", "Wells_DDimer_BTN_Sp",
                            "Wells_DDimer_BTP_Se", "Wells_DDimer_BTP_Sp", "sd", "cov_global_d","cov_global_nd",
                            "phi", "kappa", "Se_pred","Sp_pred",
                            "Wells_DDimer_BTN_Se_pred","Wells_DDimer_BTN_Sp_pred" ,"Wells_DDimer_BTP_Se_pred","Wells_DDimer_BTP_Sp_pred" ,
                            "p"),probs = c(0.025,0.5, 0.975))

saveRDS(meta_model2r, file = "Wells_random_thr_CD_kappa5_samesd_td7.rds")

#meta_model2r <-  readRDS("Wells_random_thr_CD_kappa5_samesd_td9.rds")


print(meta_model2r, pars= c("Se", "Sp", "L","M", "H", "Wells_DDimer_BTN_Se", "Wells_DDimer_BTN_Sp",
                            "Wells_DDimer_BTP_Se", "Wells_DDimer_BTP_Sp"),probs = c(0.025,0.5, 0.975))

##### from specific chains
s<- summary(subset_draws(as_draws(meta_model2r),chain=c(1:3),
                     variable = c("Se", "Sp", "Wells_DDimer_BTN_Se", "Wells_DDimer_BTN_Sp",
                                 "Wells_DDimer_BTP_Se", "Wells_DDimer_BTP_Sp", "sd1", "cov_global_d", "cov_global_nd",
                                 "rho_global_d", "rho_global_nd")), "mean",
        ~quantile(.x, probs = c(0.025, 0.50, 0.975), na.rm = TRUE), "rhat")
View(s)


s2 <-  summary(subset_draws(as_draws(meta_model2r),chain=c(1:3),
                     variable = c("a1_m_raw", "a2_m_raw", "a3_m_raw", "dc")), "mean",
        ~quantile(.x, probs = c(0.025, 0.50,  0.975), na.rm = TRUE), "rhat")

View(s2)
s2$`97.5%`
#############################################################################################
# observed - expected correlation plots 
#############################################################################################

dc <- round(summary(meta_model2r, probs = c(0.025,  0.5, 0.975), pars = c("dc"))$summary[,5],3) ; dc
dc_l <- round(summary(meta_model2r, probs = c(0.025,  0.5, 0.975), pars = c("dc"))$summary[,4],3) ; dc_l
dc_u <- round(summary(meta_model2r, probs = c(0.025,  0.5, 0.975), pars = c("dc"))$summary[,6],3) ; dc_u

#dc <- s2$`50%` ; dc
#dc_l <- s2$`2.5%` ; dc_l
#dc_u <-s2$`97.5%`; dc_u


dc_data <- tibble(dc, dc_l, dc_u, `test-pair` = as.factor(c(rep("Ref vs D-Dimer", n), rep("Ref vs Wells",n),rep("D-Dimer vs Wells",n))), 
                  obs = seq(1, length(dc), by = 1))

cutoff <- data.frame( x = c(-Inf, Inf), y = 0, cutoff = factor(0) )

#tiff("wells_figure_3_CD_random_cutpoints_homog.tif",units = "in", width = 10, height=4, res=500, compression = "lzw")

ggplot(data = dc_data, aes(y = dc, x=obs, colour = `test-pair`)) + geom_point(size = 3) + 
  geom_errorbar(aes(ymin=dc_l, ymax=dc_u), width= 0.75, position=position_dodge(.9)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  ylim(-0.31, 0.31) + 
  ylab("Correlation Residuals") + 
  xlab("Study / test-pair") + 
  theme(text = element_text(size=14),
        axis.text.x = element_text()) 

#dev.off()

#############################################################################################
# observed - expected table count plots 
#############################################################################################
dt <-   round(summary(meta_model2r, probs = c(0.025,  0.5, 0.975), pars = c("dt"))$summary[,5],3) ; dt
dt_l <- round(summary(meta_model2r, probs = c(0.025,  0.5, 0.975), pars = c("dt"))$summary[,4],3) ; dt_l
dt_u <- round(summary(meta_model2r, probs = c(0.025,  0.5, 0.975), pars = c("dt"))$summary[,6],3) ; dt_u

dt_data <- tibble(dt, dt_l, dt_u, `test-pair` = as.factor(c(rep("Ref vs D-Dimer", n*4), rep("Ref vs Wells",n*4),
                                                           rep("D-Dimer vs Wells",n*4))), 
                                                           obs = seq(1, length(dt), by = 1))

cutoff <- data.frame( x = c(-Inf, Inf), y = 0, cutoff = factor(0) )

#tiff("wells_figure_4_CD_random_cutpoints_tables_homog.tif",units = "in", width = 8, height=6, res=500, compression = "lzw")

    ggplot(data = dt_data, aes(y = dt, x=obs, colour = `test-pair`)) + 
      geom_point(size = 1) + 
      geom_errorbar(aes(ymin=dt_l, ymax=dt_u), width= 0.75, position=position_dodge(.9)) +
      geom_hline(yintercept = 0) +
      theme_bw() +
      ylim(-50, 50) + 
      ylab("2x2 table residuals") + 
      xlab("Study / table cell") + 
      theme(text = element_text(size=14),
           # axis.text.x = element_text(),
            axis.text.x=element_blank()) +
      theme(legend.position="none") +
      facet_wrap( ~ `test-pair`, scales = "free", dir = "v")

dev.off()
#############################################################################################
# DIAGNOSTIC PLOTS
#############################################################################################
stan_trace(meta_model2r, pars = c("sd"))
stan_trace(meta_model2r, pars = c("p"))
stan_trace(meta_model2r, pars = c("C_dm", "C_ndm"))
stan_trace(meta_model2r, pars = c("Se", "Sp"))
stan_trace(meta_model2r, pars = c("se", "sp"))

stan_trace(meta_model2r, pars = c("alpha"))
stan_trace(meta_model2r, pars = c("alpha_d", "alpha_nd"))

stan_trace(meta_model2r, pars = c("a1_m_raw", "a2_m_raw", "a3_m_raw"))


stan_trace(meta_model2r, pars = c("L_Omega_d"))
stan_trace(meta_model2r, pars = c("L_Omega_nd"))

stan_dens(meta_model2r, pars = c("Omega_d"))
stan_dens(meta_model2r, pars = c("L_Omega_nd"))
stan_dens(meta_model2r, pars =  c( "nu"))
stan_dens(meta_model2r, pars =  c( "sd", "z1"))
stan_dens(meta_model2r, pars =  c("Se", "Sp"))
stan_dens(meta_model2r, pars =  c("se", "sp"))
stan_dens(meta_model2r, pars =  c("p"), separate_chains = TRUE)

stan_dens(meta_model2r, pars =  c("L", "M", "H"))

#############################################################################################
# LOO
#############################################################################################
#meta_model2r <- readRDS("Wells_fixed_thr_CI.rds") 

log_lik <- extract_log_lik(meta_model2r, parameter_name = "log_lik")
str(log_lik)

draws <- 500
length <- length(log_lik[1,][log_lik[1,] != "NaN" & !is.na(log_lik[1,]) ])
loglik2 <- array(data = log_lik[log_lik != "NaN" & !is.na(log_lik) ], dim = c(draws*4, length ))
str(loglik2)
sum(is.na(loglik2))


r_eff <- relative_eff(exp(loglik2), cores = 4, chain_id = c(rep(1, draws), rep(2, draws), rep(3, draws), rep(4,draws)))

mod_loo <- loo:loo(loglik2, cores = 4, r_eff = r_eff, moment_match = TRUE)

loo:loo(meta_model2r, moment_match = TRUE)


print(mod_loo)

loo3 <- loo::loo_moment_match.default(
  x = meta_model2r,
  loo = mod_loo,
  post_draws = post_draws_stanfit,
  log_lik_i = log_lik_i_stanfit,
  unconstrain_pars = unconstrain_pars_stanfit,
  log_prob_upars = log_prob_upars_stanfit,
  log_lik_i_upars = log_lik_i_upars_stanfit
)


##############################################################################
## Figure for summary of results from 5 models for section 4.1.2 
##############################################################################

m1 <- readRDS(file = "Wells_random_thr_CI_perfectGS.rds")   ###  M1 - perfect GS, CI, random thresholds , diff SDs
m2 <- readRDS(file = "Wells_random_thr_CD_perfectGS_2.rds")   ###  M1 - perfect GS, CD, random thresholds , same SDs

#m3 <- readRDS(file = "Wells_fixed_thr_CI.rds")   ###  M2 - IGS, CI, fixed thresholds, diff SDs
m4 <- readRDS(file = "Wells_random_thr_CI.rds")   ###  M3 - IGS, CI, random thresholds, diff SDs

#m4 <- readRDS(file = "Wells_fixed_thr_CD_diffsd.rds")   ###  M2 - IGS, CD, fixed thresholds, diff SDs
#m5 <- readRDS(file = "Wells_fixed_thr_CD_td8_homogsd.rds")   ###  M2 - IGS, CD, fixed thresholds, diff SDs

#m5 <- readRDS(file = "Wells_random_thr_CD.rds")   ###  M3 - IGS, CD, random thresholds, diff SDs (I.D. issues)
m6 <- readRDS(file = "Wells_random_thr_CD_kappa5_samesd_td7.rds")   ###  M3 - IGS, CD, random thresholds, same SDs

models <- list(m1, m2, m4, m6)

print(m6, pars= c("Se", "Sp", "Wells_DDimer_BTN_Se", "Wells_DDimer_BTN_Sp",
                            "Wells_DDimer_BTP_Se", "Wells_DDimer_BTP_Sp", "rho_global_d","rho_global_nd",
                  "SeI_pred2", "SpI_pred2"),
      probs = c(0.025,0.5, 0.975))

#############################################################################################
# LOO TABLE 
#meta_model2r <- readRDS("Wells_fixed_thr_CI.rds") 

log_lik <- list()
loglik2 <- list()
draws <- list()
length <- list()
r_eff <- list()
mod_loo <- list()

require(loo)
for (i in 1:length(models)) { 
  log_lik[[i]] <- extract_log_lik(models[[i]], parameter_name = "log_lik")
  
  draws[[i]] <- dim(log_lik[[i]])[1]/4
  length[[i]] <- length(log_lik[[i]][1,][log_lik[[i]][1,] != "NaN" & !is.na(log_lik[[i]][1,]) ])
  loglik2[[i]] <- array(data = log_lik[[i]][ log_lik[[i]] != "NaN" & !is.na(log_lik[[i]] ) ] , dim = c(draws[[i]]*4, length[[i]] ))
  
  
  r_eff[[i]] <- relative_eff(exp(loglik2[[i]]), cores = 4, chain_id = c(rep(1, draws[[i]]), rep(2, draws[[i]]), rep(3, draws[[i]]), rep(4,draws[[i]])))
  
  mod_loo[[i]] <- loo(loglik2[[i]], cores = 4, r_eff = r_eff[[i]])
} 

mod_loo

loo_compare(mod_loo)

###############################################################

m_mod1 <- list(); m2_mod1 <- list()
l_mod1 <- list(); l2_mod1 <- list()
u_mod1 <- list(); u2_mod1 <- list()
m_mod1_cat <- list(); m2_mod1_cat <- list()
l_mod1_cat <- list(); l2_mod1_cat <- list()
u_mod1_cat <- list(); u2_mod1_cat <- list()
data_Se_mod <- list() 
data_Sp_mod <- list() 
data_d_mod1 <- list() 
data_nd_mod1 <- list() 

for (i in 1:length(models)) {
    m_mod1[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), 
                            pars = c("Se","Wells_DDimer_BTN_Se","Wells_DDimer_BTP_Se"))$summary[,5],2)
    
    m2_mod1[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), 
                            pars = c("Sp","Wells_DDimer_BTN_Sp","Wells_DDimer_BTP_Sp"))$summary[,5],2)
    
    
    l_mod1[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), 
                            pars = c("Se", "Wells_DDimer_BTN_Se","Wells_DDimer_BTP_Se"))$summary[,4],2)  
    
    l2_mod1[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), 
                             pars = c("Sp","Wells_DDimer_BTN_Sp","Wells_DDimer_BTP_Sp"))$summary[,4],2)
    
    
    u_mod1[[i]] <-  round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), 
                             pars = c("Se", "Wells_DDimer_BTN_Se","Wells_DDimer_BTP_Se"))$summary[,6],2)  
    
    u2_mod1[[i]]  <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), 
                              pars = c("Sp","Wells_DDimer_BTN_Sp","Wells_DDimer_BTP_Sp"))$summary[,6],2)
    
    
    # wells categories
    m_mod1_cat[[i]] <- c(round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("L[1]", "M[1]", "H[1]"))$summary[,5],2))
    m2_mod1_cat[[i]]  <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("L[2]", "M[2]", "H[2]"))$summary[,5],2)
    
    
    l_mod1_cat[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("L[1]", "M[1]", "H[1]"))$summary[,4],2)  
    l2_mod1_cat[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975),pars = c("L[2]", "M[2]", "H[2]"))$summary[,4],2)
    
    
    u_mod1_cat[[i]] <-  round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("L[1]", "M[1]", "H[1]"))$summary[,6],2)  
    u2_mod1_cat[[i]]  <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("L[2]", "M[2]", "H[2]"))$summary[,6],2)
    
    
    data_Se_mod[[i]] <- data.frame(m=m_mod1[[i]],l=l_mod1[[i]],u=u_mod1[[i]], location =  c(0.9+0.4*(i-1), 
                                                                                            3.9+0.4*(i-1), 
                                                                                            6.9+0.4*(i-1), 
                                                                                            9.9+0.4*(i-1), 
                                                                                            12.9+0.4*(i-1)), 
                               label  = factor(c("Ref", "D-Dimer", "Wells", 
                                                 "Wells & D-Dimer, BTN",
                                                 "Wells & D-Dimer, BTP"), 
                                               levels = c("Wells & D-Dimer, BTP",
                                                          "Wells & D-Dimer, BTN",
                                                          "Wells", 
                                                          "D-Dimer",
                                                          "Ref")),
                               Model = factor(rep(i,  5)))
    
    data_Sp_mod[[i]] <- data.frame(m=m2_mod1[[i]],l=l2_mod1[[i]],u=u2_mod1[[i]], location =  c(0.9+0.4*(i-1), 
                                                                                               3.9+0.4*(i-1), 
                                                                                               6.9+0.4*(i-1), 
                                                                                               9.9+0.4*(i-1), 
                                                                                               12.9+0.4*(i-1)),
                               label  = factor(c("Ref", "D-Dimer", "Wells", 
                                                 "Wells & D-Dimer, BTN",
                                                 "Wells & D-Dimer, BTP"), 
                                               levels = c("Wells & D-Dimer, BTP",
                                                          "Wells & D-Dimer, BTN",
                                                          "Wells", 
                                                          "D-Dimer",
                                                          "Ref")),
                               Model = factor(rep(i, 5)))
    
    ####### wells categories
    ## mod 1
    data_d_mod1[[i]] <- data.frame(m=m_mod1_cat[[i]], l=l_mod1_cat[[i]], u=u_mod1_cat[[i]], location =  c(0.9+0.4*(i-1), 
                                                                                                          3.9+0.4*(i-1), 
                                                                                                          6.9+0.4*(i-1)), 
                              Category  = factor(c("Low", "Medium", "High"),
                                                 levels = c("Low", "Medium", "High")),
                              Model = factor(rep(i, 3)))
    
    data_nd_mod1[[i]] <- data.frame(m=m2_mod1_cat[[i]],l=l2_mod1_cat[[i]],u=u2_mod1_cat[[i]], location =  c(0.9+0.4*(i-1), 
                                                                                                            3.9+0.4*(i-1), 
                                                                                                            6.9+0.4*(i-1)),
                               Category  = factor(c("Low", "Medium", "High"),
                                                  levels = c("Low", "Medium", "High")),
                               Model = factor(rep(i, 3)))
    
    
}

require(data.table)
# use rbindlist to put all models in same data frame
data_Se <- rbindlist(data_Se_mod)
data_Sp <- rbindlist(data_Sp_mod)
data_d <- rbindlist(data_d_mod1)
data_nd <- rbindlist(data_nd_mod1)

data_Se$Model<- factor(data_Se$Model, labels =c("M1: Perfect GS, CI", 
                                                "M2: Perfect GS, CD",
                                                "M3: Imperfect GS, CI", 
                                                "M4: Imperfect GS, CD"))

data_Se$Mod = c(rep("Perfect GS, CI", 5), 
                rep("Perfect GS, CD", 5),
                rep("Imperfect GS, CI", 5),
                rep("Imperfect GS, CD", 5))
                
  
data_Sp$Model<- factor(data_Sp$Model, labels =c("M1: Perfect GS, CI", 
                                                "M2: Perfect GS, CD",
                                                "M3: Imperfect GS, CI", 
                                                "M4: Imperfect GS, CD"))

data_Sp$Mod = c(rep("Perfect GS, CI", 5),
                rep("Perfect GS, CD", 5),
                rep("Imperfect GS, CI", 5),
                rep("Imperfect GS, CD", 5))


data_d$Model<- factor(data_d$Model, labels =c("M1: Perfect GS, CI", 
                                              "M2: Perfect GS, CD",
                                              "M3: Imperfect GS, CI", 
                                              "M4: Imperfect GS, CD"))


data_nd$Model<- factor(data_nd$Model, labels =c("M1: Perfect GS, CI", 
                                                "M2: Perfect GS, CD",
                                                "M3: Imperfect GS, CI", 
                                                "M4: Imperfect GS, CD"))


################# plot 

cat_d_plot <- ggplot(data_d, aes(x=location, y = m, ymin= l,ymax= u, shape = Model)) + 
  geom_point(size=4, alpha=0.2) + 
  geom_pointrange(aes(colour=Category)) +
  coord_flip() + 
  theme_bw()  + 
  xlab("") + 
  ylab("Diseased") + 
  theme(legend.position = "none") + 
  theme(axis.ticks = element_blank()) + 
  theme(text = element_text(size=16)) + 
  scale_y_continuous(breaks = seq(0, 1, 0.10)) + 
  scale_x_discrete( labels = rep("  ", length(data_Se$location))) + 
  geom_text(aes( label = paste0( m*100," ", "[",l*100,",", u*100 , "]") ,hjust = -0.2, vjust = -0.25), alpha = 0.35, size = 4.5) 

cat_d_plot

cat_nd_plot <- ggplot(data_nd, aes(x=location, y = m, ymin= l,ymax= u, shape = Model)) + 
  geom_point(size=4, alpha=0.2) + 
  geom_pointrange(aes(colour=Category)) + 
  coord_flip() + 
  theme_bw()  + 
  xlab("") + 
  ylab("Non-diseased") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank()) + 
  theme(text = element_text(size=16))+ 
  scale_y_continuous(breaks = seq(0, 1, 0.10)) + 
  scale_x_discrete( labels = rep("  ", length(data_Se$location))) + 
  geom_text(aes( label = paste0( m*100," ", "[",l*100,",", u*100 , "]") ,hjust = -0.2, vjust = -0.25), alpha = 0.35, size = 4.5) +
  theme(legend.title = element_text(size = 0), 
        legend.text  = element_text(size = 12))

cat_nd_plot

require(patchwork)
cat_d_plot + cat_nd_plot

####################
se_plot <- ggplot(tibble(data_Se), aes(x=as.numeric(location), y = m, ymin= l,ymax= u, shape = Model)) + 
  geom_point(size=4, alpha=0.2) + 
  geom_pointrange(aes(colour=label)) +
  coord_flip() + 
  theme_bw()  + 
  xlab("") + 
  ylab("Sensitivity") + 
  theme(legend.position = "none") + 
  theme(axis.ticks = element_blank()) + 
  theme(text = element_text(size=16))+ 
  scale_y_continuous(breaks = seq(0, 1, 0.10), limits = c(0.45,1)) + 
  scale_x_discrete( labels = rep("  ", length(data_Se$location))) + 
  geom_text(aes( label = paste0( m*100," ", "[",l*100,",", u*100 , "]") ,hjust = 1.1, vjust = -0.25), alpha = 0.35, size = 4.5) 

se_plot


sp_plot <- ggplot(data_Sp, aes(x=location, y = m, ymin= l,ymax= u, shape = Model )) + 
  geom_point(size=4, alpha=0.2) + 
  geom_pointrange(aes(colour=(label))) +
  coord_flip() + 
  theme_bw()  + 
  xlab("") + 
  ylab("Specificity")  + 
  theme(legend.title = element_blank()) + 
  theme(text = element_text(size=16))+ 
  scale_y_continuous(breaks = seq(0, 1, 0.10), limits = c(0.1, 1)) + 
  scale_x_discrete( labels = rep("  ", length(data_Se$location))) + 
  geom_text(aes( label = paste0( m*100," ", "[",l*100,",", u*100 , "]") ,hjust = 1.1, vjust = -0.25), alpha = 0.35, size = 4.5) +
  theme(legend.title = element_text(size = 0), 
        legend.text  = element_text(size = 12))


sp_plot

tiff("wells_figure_5.tif",units = "in", width = 11, height=8, res=500, compression = "lzw")
se_plot + sp_plot
dev.off()

tiff("wells_figure_6.tif",units = "in", width = 11, height=8, res=500, compression = "lzw")
cat_d_plot + cat_nd_plot
dev.off()

stan_plot(models[[1]], pars = c("p"), ncol = 1, show_density = T, ci_level = 0.80, outer_level = 0.95)
stan_trace(models[[1]], pars = c("p"))

stan_plot(models[[2]], pars = c("p"), ncol = 1, show_density = T, ci_level = 0.80, outer_level = 0.95)
stan_trace(models[[2]], pars = c("p"))

stan_plot(models[[3]], pars = c("p"), ncol = 1, show_density = T, ci_level = 0.80, outer_level = 0.95)
stan_trace(models[[3]], pars = c("p"))

stan_plot(models[[4]], pars = c("p"), ncol = 1, show_density = T, ci_level = 0.80, outer_level = 0.95)
stan_trace(models[[4]], pars = c("p"))

stan_plot(models[[5]], pars = c("p"), ncol = 1, show_density = T, ci_level = 0.80, outer_level = 0.95)
stan_trace(models[[5]], pars = c("p"))


##############################
####### ROC plot
########################################
# load in the model
mod <- models[[4]]


#tiff("figure_4.tif",units = "in", width = 11, height=7, res=500, compression = "lzw")

       
print(mod, pars= c("Se", "Sp", "Wells_DDimer_BTN_Se", "Wells_DDimer_BTN_Sp",
                  "Wells_DDimer_BTP_Se", "Wells_DDimer_BTP_Sp", "rho_global_d","rho_global_nd",
                  "Se_pred", "Sp_pred"),
      probs = c(0.025,0.5, 0.975))


require(bayesSurv)
require(scales)


#####
## credible region
cred_1 <- tibble(y = qnorm2(extract(mod, pars = "Se")$'Se'[,1]), x = qnorm2(extract(mod, pars = "Sp")$'Sp'[,1]))
cred_1p <- tibble(y = (extract(mod, pars = "Se")$'Se'[,1]), x =  (extract(mod, pars = "Sp")$'Sp'[,1]))

# in inv_probit space
g <- ggplot(data = cred_1, aes(x = x, y = y))  + 
  geom_point() + 
  stat_ellipse(type = "t")  

g

# Get ellipse coordinates from plot
pb = ggplot_build(g)
el = pb$data[[2]][c("x","y")]

credible_region <- tibble(x = pnorm(el$x), y = pnorm(el$y))
credible_region


# in probability space
ggplot(data = cred_1p, aes(x = 1- x, y =  y))  + 
  #geom_point(size = 0.1) + 
  geom_path(data = credible_region, aes(x = 1  - x, y = y), colour = "blue") + 
  ylim(0,1) +  
  xlim(0,1)

####
## prediction region

pred_1 <- tibble(y = qnorm2(extract(mod, pars = "SeI_pred2")$'SeI_pred2'[,1]), x = qnorm2(extract(mod, pars = "SpI_pred2")$'SpI_pred2'[,1]))
pred_1p <- tibble(y = (extract(mod, pars = "SeI_pred2")$'SeI_pred2'[,1]), x =  (extract(mod, pars = "SpI_pred2")$'SpI_pred2'[,1]))

# in inv_probit space
g <- ggplot(data = pred_1, aes(x = x, y = y))  + 
  geom_point() + 
  stat_ellipse(type = "t")  

g

# Get ellipse coordinates from plot
pb = ggplot_build(g)
el = pb$data[[2]][c("x","y")]

pred_region <- tibble(x = pnorm(el$x), y = pnorm(el$y))
pred_region


# in probability space
ggplot(data = pred_1p, aes(x = 1- x, y =  y))  + 
  #geom_point(size = 0.1) + 
  geom_path(data = credible_region, aes(x = 1  - x, y = y), colour = "blue") + 
  ylim(0,1) +  
  xlim(0,1)

#############################
## plot

g <- ggplot(data = means2, aes(y=y, x = x, colour = Test)) +    # summary points
  geom_point( size=4 ) + 
 # xlim(0, 1) + 
  #ylim(0, 1) + 
  geom_path(data = filter(pred_region2,Test == 1), aes(x=X, y=Y, colour = Test), linetype = 2, size = 0.4, inherit.aes = F) +                         # prediction region
  geom_polygon(data = filter(conf_region2,Test == 1), aes(x=X, y=Y, colour = Test), alpha=0.05, size=0.4, linetype = 2,inherit.aes = F) + # conf region
  theme_bw() + 
  scale_x_continuous(breaks = seq(0,1,0.1), limits = c(0,1))  + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1))  +
  theme(legend.title=element_blank()) + 
  xlab("1 - Specifitity") + 
  ylab("Sensitivity") + 
  geom_path(data = credible_region, aes(x = 1  - x, y = y), colour = "blue")  + 
  geom_path(data = pred_region, aes(x = 1  - x, y = y), colour = "blue") 
g


dev.off()


##############################################################################
## Figure for summary of results from 4 models for section 4.1.1 (dichotomous Wells)
##############################################################################

## first dichot 
d1_ci <- readRDS(file = "Wells_dichotomous_CI_1st_dichot.rds")   ###  1st dichot, M1 - perfect GS, CI, random thresholds , diff SDs
d2_ci <- readRDS(file = "Wells_dichotomous_CI_2nd_dichot.rds")   ###  M2 - IGS, CI, fixed thresholds, diff SDs

d1_cd <- readRDS(file = "Wells_dichotomous_CD_1st_dichot_td9.rds")   ###  M3 - IGS, CI, random thresholds, diff SDs
d2_cd <- readRDS(file = "Wells_dichotomous_CD_2nd_dichot.rds")   ###  M2 - IGS, CD, fixed thresholds, diff SDs

models2 <- list(d1_ci, d2_ci, d1_cd, d2_cd)

###############################################################

m_mod1 <- list(); m2_mod1 <- list()
l_mod1 <- list(); l2_mod1 <- list()
u_mod1 <- list(); u2_mod1 <- list()
m_mod1_cat <- list(); m2_mod1_cat <- list()
l_mod1_cat <- list(); l2_mod1_cat <- list()
u_mod1_cat <- list(); u2_mod1_cat <- list()
data_Se_mod <- list() 
data_Sp_mod <- list() 
data_d_mod1 <- list() 
data_nd_mod1 <- list() 

for (i in 1:length(models2)) {
  m_mod1[[i]] <- round(summary(models2[[i]], probs = c(0.025,  0.5, 0.975), 
                               pars = c("Se","Wells_DDimer_BTN_Se","Wells_DDimer_BTP_Se"))$summary[,5],2)
  
  m2_mod1[[i]] <- round(summary(models2[[i]], probs = c(0.025,  0.5, 0.975), 
                                pars = c("Sp","Wells_DDimer_BTN_Sp","Wells_DDimer_BTP_Sp"))$summary[,5],2)
  
  
  l_mod1[[i]] <- round(summary(models2[[i]], probs = c(0.025,  0.5, 0.975), 
                               pars = c("Se", "Wells_DDimer_BTN_Se","Wells_DDimer_BTP_Se"))$summary[,4],2)  
  
  l2_mod1[[i]] <- round(summary(models2[[i]], probs = c(0.025,  0.5, 0.975), 
                                pars = c("Sp","Wells_DDimer_BTN_Sp","Wells_DDimer_BTP_Sp"))$summary[,4],2)
  
  
  u_mod1[[i]] <-  round(summary(models2[[i]], probs = c(0.025,  0.5, 0.975), 
                                pars = c("Se", "Wells_DDimer_BTN_Se","Wells_DDimer_BTP_Se"))$summary[,6],2)  
  
  u2_mod1[[i]]  <- round(summary(models2[[i]], probs = c(0.025,  0.5, 0.975), 
                                 pars = c("Sp","Wells_DDimer_BTN_Sp","Wells_DDimer_BTP_Sp"))$summary[,6],2)
  
  data_Se_mod[[i]] <- data.frame(m=m_mod1[[i]],l=l_mod1[[i]],u=u_mod1[[i]], location =  c(0.9+0.3*(i-1), 
                                                                                          2.9+0.3*(i-1),
                                                                                          4.9+0.3*(i-1), 
                                                                                          6.9+0.3*(i-1), 
                                                                                          8.9+0.3*(i-1)), 
                                 label  = factor(c("Ref", "D-Dimer", "Wells", 
                                                   "Wells & D-Dimer, BTN",
                                                   "Wells & D-Dimer, BTP"), 
                                                 levels = c("Wells & D-Dimer, BTP",
                                                            "Wells & D-Dimer, BTN",
                                                            "Wells", 
                                                            "D-Dimer",
                                                            "Ref")),
                                 Model = factor(rep(i,  5)))
  
  data_Sp_mod[[i]] <- data.frame(m=m2_mod1[[i]],l=l2_mod1[[i]],u=u2_mod1[[i]], location =  c(0.9+0.3*(i-1), 
                                                                                             2.9+0.3*(i-1),
                                                                                             4.9+0.3*(i-1),
                                                                                             6.9+0.3*(i-1),
                                                                                             8.9+0.3*(i-1)),
                                 label  = factor(c("Ref", "D-Dimer", "Wells", 
                                                   "Wells & D-Dimer, BTN",
                                                   "Wells & D-Dimer, BTP"), 
                                                 levels = c("Wells & D-Dimer, BTP",
                                                            "Wells & D-Dimer, BTN",
                                                            "Wells", 
                                                            "D-Dimer",
                                                            "Ref")),
                                 Model = factor(rep(i, 5)))
}

require(data.table)
# use rbindlist to put all models in same data frame
data_Se2 <- rbindlist(data_Se_mod)
data_Sp2 <- rbindlist(data_Sp_mod)

data_Se2$Dichot = c(rep("1st dichot.", 5),
                    rep("2nd dichot.", 5),
                    rep("1st dichot.", 5),
                    rep("2nd dichot.", 5))

data_Se2$Mod   =  c(rep("CI", 10),
                    rep("CD", 10))

data_Se2$Model<- factor(data_Se2$Model, labels =c("1st dichot., CI", 
                                                "2nd dichot., CI", 
                                                "1st dichot., CD", 
                                                "2nd dichot., CD"))

data_Sp2$Model<- factor(data_Sp2$Model, labels =c("1st dichot., CI", 
                                                "2nd dichot., CI", 
                                                "1st dichot., CD", 
                                                "2nd dichot., CD"))

data_Sp2$Dichot = c(rep("1st dichot.", 5),
                    rep("2nd dichot.", 5),
                    rep("1st dichot.", 5),
                    rep("2nd dichot.", 5))

data_Sp2$Mod   =  c(rep("CI", 10),
                    rep("CD", 10))

################# plot 

require(ggrepel)

se_plot <- ggplot(tibble(data_Se2), aes(x=(location), y = m, ymin= l,ymax= u, shape = Mod)) + 
  geom_point(size=4, alpha=0.2) + 
  geom_pointrange(aes(colour=label, linetype = Dichot)) +
  coord_flip() + 
  theme_bw()  + 
  xlab("") + 
  ylab("Sensitivity") + 
  theme(legend.position = "none") + 
  theme(axis.ticks = element_blank()) + 
  theme(text = element_text(size=16))+ 
  scale_y_continuous(breaks = seq(0, 1, 0.10)) + 
  scale_x_discrete( labels = rep("  ", length(data_Se2$location))) + 
  geom_text(aes( label = paste0( m*100," ", "[",l*100,",", u*100 , "]") ,hjust = 1.1, vjust = -0.25),alpha = 0.35) 
#  scale_x_discrete( labels = c(" ", " ", " ", "Ref", 
#                               " ", " ", " ", "D-Dimer",
 #                              " ", " ", " ", "Wells", 
  #                             " ", " ", " ", "Wells & D-Dimer, BTN",
   #                            " ", " ", " ", "Wells & D-Dimer, BTP"))

se_plot


sp_plot <- ggplot(data_Sp2, aes(x=location, y = m, ymin= l,ymax= u, shape = Mod)) + 
  geom_point(size=4, alpha=0.2) + 
  geom_pointrange(aes(colour=label, linetype = Dichot)) +
  coord_flip() + 
  theme_bw()  + 
  xlab("") + 
  ylab("Specificity")  + 
  theme(legend.title = element_blank()) + 
  theme(text = element_text(size=16))+ 
  scale_y_continuous(breaks = seq(0, 1, 0.10))+ 
  scale_x_discrete( labels = rep("  ", length(data_Se2$location))) + 
  geom_text(aes( label = paste0( m*100," ", "[",l*100,",", u*100 , "]") ,hjust = 1.1, vjust = -0.25), alpha = 0.35) 


sp_plot


tiff("wells_figure_1.tif",units = "in", width = 11, height=7, res=500, compression = "lzw")
se_plot + sp_plot
dev.off()

### plot of disease prevalences for 1st dichot vs 2nd dichot (CI on left, CD on right)



g1 <- stan_plot(models2[[1]], pars = c("p"), ncol = 1, show_density = T, ci_level = 0.95, outer_level = 0.95)
g1

g2<- stan_plot(models2[[2]], pars = c("p"), ncol = 1, show_density = T, ci_level = 0.95, outer_level = 0.95)
g2

params <- extract(models2[[1]])



g1 + g2

g1 <- stan_plot(models2[[3]], pars = c("p"), ncol = 1, show_density = T, ci_level = 0.95, outer_level = 0.95)
g1

g2<- stan_plot(models2[[4]], pars = c("p"), ncol = 1, show_density = T, ci_level = 0.95, outer_level = 0.95)
g2

g1 + g2
#########################################################################################################################
#### standard latent class model
file <- file.path(file = "latent_class_wells.stan")
mod <- cmdstan_model(file)


data <- list(x=2, y1 = 1,
             NS= 11 ,
             ns = ns, 
             n=ns,
             nt = 3,
             r =   r_full, 
             r2 = r,
             num_ref = 1,
             T = rep(1, 11),
             kappa_model = 0,
             c2_m = 0, c1_m = 1, s2_m = 0, s1_m = 1,
             c2_sd = 1.5, c1_sd = 0.4, s1_sd = 1.5, s2_sd = 0.4, 
             c2_sdv = 1, c1_sdv = 1, s1_sdv = 1, s2_sdv = 1, 
             ci = 1, 
             c1r = 1, s1r = 1, s2r = 2, c2r = 1,
             sd_kappa_s1 = array(2),
             sd_kappa_c1 = array(2))

init <- list(mu_beta1 = array(2), mu_alpha1 = c(2), 
             mu_beta2 = c(2), mu_alpha2 = c(2),
             mu_beta3 = c(-0.5), mu_alpha3 = c(-0.5), 
             mu_beta4 = c(-0.5), mu_alpha4 = c(-0.5), 
             p = prev)

meta_model_lcm <- mod$sample(
  data =  data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 200,
  iter_sampling = 200, 
  refresh = 20, 
  init = list(init,init,init,init), 
  max_treedepth = 15,
  adapt_delta = 0.98)


meta_model_lcm2 <- rstan::read_stan_csv(meta_model_lcm$output_files())


print(meta_model_lcm2, pars = c("SeR", "SpR", "SeDD", "SpDD", "L_d", "M_d", "H_d", "L_nd", "M_nd", "H_nd"),  probs = c(0.025, 0.5,    0.975))
print(meta_model_lcm2, pars = c("s2", "c2"),  probs = c(0.025, 0.5,    0.975))
print(meta_model_lcm2, pars = c("dc"),  probs = c(0.025, 0.5,    0.975))


## extract draws from specific chains
chain <- 1
round(summary(meta_model_lcm2, 
              probs = c(0.025,  0.5, 0.975), 
              pars= c( "SeR" ) )$c_summary[, , chain][,3],2) 



round(summary(meta_model_lcm2, 
              probs = c(0.025,  0.5, 0.975), 
              pars= c( "dc") )$c_summary[, , chain][,3],2) 


stan_dens(meta_model_lcm2, pars = c("SeI", "SpI", "SeR", "SpR"))






#######################
dc <- round(summary(meta_model_lcm2, probs = c(0.025,  0.5, 0.975), pars = c("dc"))$summary[,5],3) ; dc
dc_l <- round(summary(meta_model_lcm2, probs = c(0.025,  0.5, 0.975), pars = c("dc"))$summary[,4],3) ; dc_l
dc_u <- round(summary(meta_model_lcm2, probs = c(0.025,  0.5, 0.975), pars = c("dc"))$summary[,6],3) ; dc_u

dc_data <- tibble(dc, dc_l, dc_u, Comparison = as.factor(c(rep("Ref vs D-Dimer", 11), rep("Ref vs Wells",11),rep("D-Dimer vs Wells",11))), 
                  obs = seq(1, length(dc), by = 1))

cutoff <- data.frame( x = c(-Inf, Inf), y = 0, cutoff = factor(0) )

#tiff("wells_m3.tif",units = "in", width = 15, height=7, res=500, compression = "lzw")
ggplot(data = dc_data, aes(y = dc, x=obs, colour = Comparison)) + geom_point(size = 3) + 
  geom_errorbar(aes(ymin=dc_l, ymax=dc_u), width= 0.75, position=position_dodge(.9)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  ylim(-0.3, 0.3) + 
  ylab("Correlation Residuals (95% CrI)") + 
  xlab(" ") + 
  theme(text = element_text(size=14),
        axis.text.x = element_text()) 










