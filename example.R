#### load packages and functions ####
library(data.table)
library(dplyr)
library(cmdstanr)
library(posterior)

#set working directory to this repo
setwd("~/repos/parse-Stan/")

#read in functions
source("R/functions.R")

#### load in the Stan model and fitted MCMC output ####
base_dir <- "~/repos/parse-Stan/"
stan_model_dir <- paste0(base_dir, "models/")
stan_data_dir <- paste0(base_dir, "data/")
stan_output_dir <- paste0(base_dir, "output/")
stan_progress_dir <- paste0(base_dir, "progress/")

#specify model and data object name
model_name <- "beta-binomial"

#load in Stan model file
stan_code <- readLines(paste0(stan_model_dir, model_name, ".stan")) #load stan model spec

#specify input data as list
print(retrieve_inputs(stan_code))
dat <- list(n = 100)
dat$k <- sample(1:1000, size = dat$n, replace = T)
dat$x <- rbinom(dat$n, size = dat$k, prob = rbeta(dat$n, 5, 5))

# process Stan file for prior predictive simulation
r_code <- parse_Stan(stan_code, 
                     dat = dat, 
                     samps = NA, 
                     output_file = NA, 
                     sample_index = NA, 
                     post_pred_sim = F, 
                     sim = TRUE)

#evaluate code
stan_env <- new.env()
sample_index <- sample(nrow(samps), 1)
eval(parse(text = r_code), envir = stan_env)
out_sim <- stan_env$out

#write data to disk for inference in Stan
dat_sim <- out_sim[names(dat)]
data_name <- "simulated_data"
write_stan_json(dat_sim, file = paste0(stan_data_dir, data_name,".json"), 
                always_decimal = FALSE)
write_stan_json(out_sim, file = paste0(stan_output_dir, data_name,"_params.json"), 
                always_decimal = FALSE)

#fit model to simulated data


#read back in fitted model and data
load(paste0(stan_output_dir, model_name, ".cmdStanR.fit")) #load Stan output
samps <- setDT(data.frame(as_draws_df(out$draws())))
dat <- jsonlite::fromJSON(paste0(stan_data_dir, model_name, ".json"))

#check MCMC diagnostics just in case
load(paste0(stan_output_dir, model_name, "_", data_name, ".cmdStanR.summ"))
summ[order(summ$ess_bulk), c(colnames(summ)[1:2], "rhat", "ess_bulk")][1:20,]
summ[order(summ$rhat, decreasing = T),c(colnames(summ)[1:2], "rhat", "ess_bulk")][1:20,]

#process Stan model for posterior predictive simulation

#process Stan model to evaluate posterior predictive densities / massess

#check calibration
n_samp <- 500
sample_inds <- sample(1:nrow(samps_sim), n_samp)
mcprint <- function(...){
  system(sprintf('echo "%s"', paste0(..., collapse="")))
}

calib_list <- parallel::mclapply(1:n_samp, function(si){
  mcprint(paste0(si, " "))
  
  #generate and process code
  r_sim_code <- parse_Stan(stan_code, dat, samps_sim, output_file = NA, 
                           sample_index = sample_inds[si], post_pred_sim = TRUE, sim = TRUE)
  r_sim_code <- strsplit(x = r_sim_code, split = "shape_alpha")[[1]][1]
  
  #evaluate code
  stan_sim_env <- new.env()
  eval(parse(text = r_sim_code), envir = stan_sim_env)
  out_sim_sim <- stan_sim_env$out
  
  return(out_sim_sim)
}, mc.cores = 8)

# Initialize calib with the same structure as out_sim but with zeros
out_sim_calib <- out_sim[names(calib_list[[1]])]
calib <- lapply(out_sim_calib, function(item) {
  if (length(dim(item)) == 2) {
    return(matrix(0, nrow = nrow(item), ncol = ncol(item)))
  } else {
    return(numeric(length(item)))
  }
})

# count when the inferred parameter falls above or below the true value
for (i in seq_along(calib_list)) {
  cat(paste0(i, " "))
  for (param_name in names(out_sim_calib)) {
    calib[[param_name]] <- calib[[param_name]] + (out_sim_calib[[param_name]] < calib_list[[i]][[param_name]])
  }
}

#set all values in parameters fixed to 0 to have calibration NA
for (i in seq_along(calib_list)) {
  for (param_name in names(out_sim_calib)) {
    calib[[param_name]][out_sim_calib[[param_name]] == 0] <- NA
  }
}

# After accumulation, divide by the number of draws to get the proportion
calib <- lapply(calib, function(x) x / length(calib_list))
calib <- calib[!(names(calib) %in% c("log_x", names(dat)))]
hist(unlist(calib), breaks = 0:10/10, xlab = "simulated parameter quantile in posterior", 
     main = "Calibration Test of 1x Simulated Data", freq = F)

text(x = mean(par("usr")[1:2]), y = par("usr")[4], cex = 0.75, xpd = NA, pos = 3,
     labels = "(NOTE: data simulated from empirical posterior. Perfect calibration only expected when using prior predictive in aggregate)")


#try making a similar figure to before

#generate and process code
r_sim_code <- parse_Stan(stan_code, dat, samps_sim, output_file = NA, 
                         sample_index = sample_inds[si], post_pred_sim = TRUE, sim = TRUE)

#evaluate code
stan_sim_env <- new.env()
eval(parse(text = r_sim_code), envir = stan_sim_env)
out_sim_sim <- stan_sim_env$out

#process output
flip_sim <- abs(out_sim_sim$count - dat_sim$count) > abs((out_sim_sim$total - out_sim_sim$count) - dat_sim$count) & (out_sim_sim$count_mixcomp != 1)
best_count_sim <- out_sim_sim$count
best_count_sim[flip_sim] <- (out_sim_sim$total - out_sim_sim$count)[flip_sim]

#plot empirical vs simulated outcomes
plot((best_count_sim/dat_sim$total), (dat_sim$count/dat_sim$total), pch = 19, 
     col = adjustcolor(viridisLite::viridis(101), 0.1)[ceiling(log(dat_sim$total) / max(log(dat_sim$total)) * 100)], 
     xlab = "log10(mixture-flipped count)", ylab = "log10(observed count)", main = "single posterior predictive sample")
abline(0,1,col=2,lty=2,lwd=2)
cor(log10((best_count_sim+1)/dat$total), log10((dat$count+1)/dat$total), use = "pairw")

#### check retrodictive accuracy with posterior mean ####
n_samp <- 500
sample_inds <- sample(1:nrow(samps), n_samp)

post_list <- parallel::mclapply(1:n_samp, function(si){
  mcprint(paste0(si, " "))
  
  #generate and process code
  r_sim_code <- parse_Stan(stan_code, dat, samps, output_file = NA, 
                           sample_index = sample_inds[si], post_pred_sim = TRUE, sim = TRUE)
  
  #evaluate code
  stan_sim_env <- new.env()
  eval(parse(text = r_sim_code), envir = stan_sim_env)
  out_sim_sim <- stan_sim_env$out
  
  return(out_sim_sim)
}, mc.cores = 8)

sim_counts <- do.call(cbind, lapply(1:n_samp, function(i) post_list[[i]]$count))
sim_counts_side <- do.call(cbind, lapply(1:n_samp, function(i) post_list[[i]]$count_mixcomp))

#look at modes
mode_heights <- function(x){
  dens_est <- density(x)
  deriv_1 <- diff(dens_est$y) / diff(dens_est$x)
  cross_0 <- cumsum(rle(deriv_1 > 0)$lengths)
  deriv_2 <- c(0, diff(deriv_1) / diff(dens_est$x)[-1])
  return(dens_est$y[cross_0[deriv_2[cross_0] < 0]])
}

modiness <- sapply(1:nrow(sim_counts), function(i){
  mode_locs <- mode_heights(sim_counts[i,])
  length(mode_locs) / var(mode_locs)
})
modiness[is.na(modiness)] <- 0

#now compute posterior mean

#first flip the signs if the true obs is on the other side 
#(for count_mixcomp, 1 : prob = 0.5, 2 : prob > 0.5, 3 : prob < 0.5)
sim_counts_dev <- abs(sim_counts - dat_sim$total / 2)
ppd_mean <- apply(sim_counts_dev, 1, mean)
dat_sim_dev <- abs(dat_sim$count - dat_sim$total / 2)

#plot deviation from 0.5
ppd_prop_dev <- ppd_mean /dat_sim$total
sim_prop_dev <- dat_sim_dev / dat_sim$total

# ppd_prop_dev <- log10(ppd_prop_dev)
# sim_prop_dev <- log10(sim_prop_dev)

nInf_inds <- (ppd_prop_dev == -Inf | sim_prop_dev == -Inf)
ppd_prop_dev <- ppd_prop_dev[!nInf_inds]
sim_prop_dev <- sim_prop_dev[!nInf_inds]

pearson_r <- cor(ppd_prop_dev, sim_prop_dev)
col_scale <- (viridisLite::viridis(100))
pt_col_inds <- ceiling(log10(dat_sim$total) / max(log10(dat_sim$total)) * 100)
legend_vals <- seq(min(log10(dat_sim$total)), max(log10(dat_sim$total)), length.out = 10)
legend_col_inds <- ceiling(legend_vals / max(log10(dat_sim$total)) * 100)
ifelse2 <- function(test, yes, no) if(test){return(yes)}else{return(no)}
plot(x = ppd_prop_dev, pch = 19,
     y = sim_prop_dev, 
     col = adjustcolor(col_scale, 0.1)[pt_col_inds], 
     xlab = "posterior mean proportion deviation from 0.5", 
     ylab = "observed freq deviation from 0.5", 
     main = paste0(ifelse(any(sim_prop_dev < 0), "log10 scale, ", ""), "posterior predictive mean, r = ", round(pearson_r, 3)), 
     xlim = ifelse2(any(ppd_prop_dev < 0), NULL, c(0,0.5)), 
     ylim = ifelse2(any(sim_prop_dev < 0), NULL, c(0,0.5)))
legend(x = "bottomright", legend = round(10^legend_vals), pch = 19, col = col_scale[legend_col_inds], bty="n", title = "sample\nsize", title.font = 2)
abline(0,1,col=2,lty=2,lwd=2)


#can also look at width of hpdi necessary to contain true observation

