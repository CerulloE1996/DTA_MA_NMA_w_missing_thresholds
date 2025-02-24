K <- 3
prior_kappa_sd <- 50


#############################
# prior predictive check plot

require(ggplot2)
# require(hexbin)
# require(latex2exp)


N <- 1e4
K <- 5
prior_kappa_sd <- 10



induced_Dirichlet_ppc_plot <- function(  use_alpha_directly, 
                                         use_log_kappa, 
                                         prior_sd, 
                                         n_cat, 
                                         N = 1e4) {
  
        p_if_uniform <- 1/n_cat
          
        message(paste("p_if_uniform = ", p_if_uniform))
        
 
        if (use_alpha_directly == TRUE) {
            if (use_log_kappa == TRUE) {
                log_kappa <- rnorm(n = N, mean = 0, sd = prior_sd)
                kappa <- exp(log_kappa)
            } else { 
                kappa <- truncnorm::rtruncnorm(N, a = 0, mean =  0, sd = prior_sd)
            }
        }
        
        phi <- MCMCpack::rdirichlet(N, rep(1, n_cat)) 
        res <- array(NA_real_, dim = c(n_cat, N))
        
        for (i in 1:N) {
          if (use_alpha_directly == TRUE) {
            log_alpha <- rnorm(n = n_cat, mean = 0, sd = prior_sd)
            alpha <- exp(log_alpha)
          } else { 
            alpha <- kappa[i]*phi[i,]
          }
          
          res[,i] <- MCMCpack::rdirichlet(1, alpha)
        }
        
        df <- data.frame(p1 = res[1,], p2 = res[2,], dist = 1) %>% filter(!is.na(p1), !is.na(p2))     
        
        g1 <- ggplot(df, aes(p1, p2)) + 
          geom_hex() + 
          scale_fill_continuous(trans = "log10") +
          theme_bw() + 
          # xlab(TeX("$P_{i}$")) + 
          # ylab(TeX("$P_{j}$")) + 
          xlim(0, 3.0*p_if_uniform) + 
          ylim(0, 3.0*p_if_uniform) + 
          geom_hline(yintercept = p_if_uniform) + 
          geom_vline(xintercept = p_if_uniform)
        
        print(g1)

}










