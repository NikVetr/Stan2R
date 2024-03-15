set.seed(1)

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
stan_figure_dir <- paste0(base_dir, "figures/")
stan_mcmc_dir <- paste0(stan_output_dir, "mcmc/")
stan_progress_dir <- paste0(base_dir, "progress/")

#specify model and data object name
model_name <- "beta-binomial"

#load in Stan model file
model_path <- paste0(stan_model_dir, model_name, ".stan")
stan_code <- readLines(model_path) #load stan model spec

#specify input data as list
print(retrieve_block_objects(stan_code, "data", "declaration"))
dat <- list(n = 1000)
dat$k <- sample(1:500, size = dat$n, replace = T)
dat$x <- rbinom(dat$n, size = dat$k, prob = rbeta(dat$n, 5, 5))
dat$n_groups <- 50
dat$group <- sample(1:dat$n_groups, size = dat$n, replace = T)

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
eval(parse(text = r_code), envir = stan_env)
prior_predictive_sim <- stan_env$out

#write data to disk for inference in Stan
dat_sim <- prior_predictive_sim[names(dat)]
write_stan_json(dat_sim, file = paste0(stan_data_dir, model_name, "_simulated-data.json"), 
                always_decimal = FALSE)
write_stan_json(prior_predictive_sim, file = paste0(stan_output_dir, model_name,"_simulated-params.json"), 
                always_decimal = FALSE)

#fit model to simulated data
mod <- cmdstan_model(model_path)
fit <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, 
                  data = dat_sim, parallel_chains = 4, adapt_delta = 0.9, 
                  refresh = 100, max_treedepth = 10, output_dir = stan_mcmc_dir,
                  thin = 1, init = 1)

#check MCMC diagnostics just in case
summ <- fit$summary(variables = retrieve_block_objects(stan_code, "parameters", "declaration")$var_name)
print(summ[order(summ$ess_bulk),])
print(summ[order(summ$rhat, decreasing = T),])

#save fitted model to disk
save(fit, file = paste0(stan_output_dir, model_name, ".cmdStanR.fit"))
save(summ, file = paste0(stan_output_dir, model_name, ".cmdStanR.summ"))

#read back in fitted model and data
load(paste0(stan_output_dir, model_name, ".cmdStanR.fit")) #load Stan output
samps <- setDT(data.frame(as_draws_df(fit$draws())))
dat <- jsonlite::fromJSON(paste0(stan_data_dir, model_name, "_simulated-data.json"))

#process Stan model for posterior predictive simulation
sample_index <- sample(nrow(samps), 1)

#process Stan model to evaluate posterior predictive densities / massess
r_code <- parse_Stan(stan_code, 
                     dat = dat, 
                     samps = samps, 
                     output_file = NA, 
                     sample_index = sample_index, 
                     post_pred_sim = T, 
                     sim = TRUE)

#evaluate code
stan_env <- new.env()
eval(parse(text = r_code), envir = stan_env)
posterior_predictive_sim <- stan_env$out

#check calibration of inference on simulated data
n_samp <- 500
sample_inds <- sample(1:nrow(samps), n_samp)
mcprint <- function(...){
  system(sprintf('printf "%s"', paste0(..., collapse="")))
}

post_pred_list <- parallel::mclapply(1:n_samp, function(si){
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

# Initialize calib with the same structure as out_sim but with zeros
sim_calib <- prior_predictive_sim[names(post_pred_list[[1]])]
calib <- lapply(sim_calib, function(item) {
  if (length(dim(item)) == 2) {
    return(matrix(0, nrow = nrow(item), ncol = ncol(item)))
  } else {
    return(numeric(length(item)))
  }
})

# count when the inferred parameter falls above or below the true value
for (i in seq_along(post_pred_list)) {
  cat(paste0(i, " "))
  for (param_name in names(sim_calib)) {
    calib[[param_name]] <- calib[[param_name]] + (sim_calib[[param_name]] < post_pred_list[[i]][[param_name]])
  }
}

#set all values in parameters fixed to 0 to have calibration NA
for (i in seq_along(post_pred_list)) {
  for (param_name in names(sim_calib)) {
    calib[[param_name]][sim_calib[[param_name]] == 0] <- NA
  }
}

# After accumulation, divide by the number of draws to get the proportion
calib <- lapply(calib, function(x) x / length(post_pred_list))

#remove input data and transformed data
variables_to_ignore <- setdiff(c(retrieve_block_objects(stan_code, "data", "declaration")$var_name, 
                                 retrieve_block_objects(stan_code, "transformed data", "declaration")$var_name), 
                               retrieve_block_objects(stan_code, "model", "sampling")$var_name)
calib <- calib[!(names(calib) %in% variables_to_ignore)]

#examine calibration of output
png(filename = paste0(stan_figure_dir, "outcome_calibration.png"), width = 800, height = 500, pointsize = 15)
par(mar = c(5,5,4,2))
hist(unlist(calib$x), breaks = 0:10/10, xlab = "simulated parameter quantile in posterior", 
     main = "Calibration Test of 1x Simulated Data", freq = F)
text(x = mean(par("usr")[1:2]), y = par("usr")[4], cex = 0.75, xpd = NA, pos = 3,
     labels = "(NOTE: data simulated from empirical posterior. Perfect calibration only expected when using prior predictive in aggregate)")
dev.off()

#examine calibration of parameters, ignoring scale
png(filename = paste0(stan_figure_dir, "parameter_calibration.png"), width = 800, height = 500, pointsize = 15)
par(mar = c(5,5,4,2))
hist(unlist(calib[setdiff(names(calib), "x")]), breaks = 0:10/10, xlab = "simulated parameter quantile in posterior", 
     main = "Calibration Test of 1x Simulated Data", freq = F)
text(x = mean(par("usr")[1:2]), y = par("usr")[4], cex = 0.75, xpd = NA, pos = 3,
     labels = "(NOTE: data simulated from empirical posterior. Perfect calibration only expected when using prior predictive in aggregate)")
dev.off()

#check retrodictive accuracy with posterior mean
#also assessing posterior predictive density (mass for our beta-binomial model)
unobs_vars

#first let's do the model parameters -- can get directly from samps, or by processing post_pred_list
unobs_var_names <- names(prior_predictive_sim)[!(names(prior_predictive_sim) %in% c(variables_to_ignore, "x"))]
unobs_var_names <- setNames(unobs_var_names, unobs_var_names)
unobs_vars <- prior_predictive_sim[unobs_var_names]
post_pred_dists_unobs_vars <- lapply(unobs_var_names, function(uv_name) 
  do.call(rbind, lapply(post_pred_list, function(post_pred_sim) post_pred_sim[[uv_name]])))
post_pred_mean_unobs_vars <- lapply(unobs_var_names, function(uv_name) apply(post_pred_dists_unobs_vars[[uv_name]], 2, mean))
post_pred_95q_unobs_vars <- lapply(unobs_var_names, function(uv_name) apply(post_pred_dists_unobs_vars[[uv_name]], 2, quantile, probs = c(0.05, 0.95)))

#now do some plotting
dim_param_fig <- ceiling(sqrt(length(unobs_vars)))
if((dim_param_fig^2-dim_param_fig) > length(unobs_var_names)){
  dim_param_fig <- c(dim_param_fig, dim_param_fig-1)
} else {
  dim_param_fig <- c(dim_param_fig, dim_param_fig)
}

png(filename = paste0(stan_figure_dir, "parameter_goodness-of-fit.png"), width = dim_param_fig[2] * 500, height = dim_param_fig[1] * 500, pointsize = 22)
par(mar = c(5,5,3.5,2), mfrow = dim_param_fig)
for(uv_name in unobs_var_names){
  
  varlims <- range(c(post_pred_95q_unobs_vars[[uv_name]], unobs_vars[[uv_name]]))
  plot(unobs_vars[[uv_name]], post_pred_mean_unobs_vars[[uv_name]], pch = 19, col = adjustcolor(1, 0.5), cex = 0.5, #dat$k / max(dat$k),
      xlab = paste0("observed ", uv_name, " (from prior predictive sample)"), ylab = paste0("posterior ", uv_name), ylim = varlims, xlim = varlims,
       main = paste0(uv_name))
  segments(x0 = unobs_vars[[uv_name]], x1 = unobs_vars[[uv_name]], y0 = post_pred_95q_unobs_vars[[uv_name]][1,], y1 = post_pred_95q_unobs_vars[[uv_name]][2,], col = adjustcolor(1, 0.5))
  abline(0,1,col=2,lty=2,lwd=2)
  legend(x = "topleft", legend = c("posterior mean", "95% credible interval", "1-to-1 line"), pch = c(19, NA,NA), lty = c(NA,1,2), col = c(1,1,2))
  text(x = mean(par("usr")[1:2]), y = par("usr")[4], pos = 3, xpd = NA,
       labels = paste0("Pearson's r(observed value, posterior mean) = ", round(cor(unobs_vars[[uv_name]], post_pred_mean_unobs_vars[[uv_name]]), 3)))
  
}
dev.off()

#now do the same for posterior means
post_pred_dists_outcome <- do.call(rbind, lapply(post_pred_list, function(post_pred_sim) post_pred_sim$x))
post_pred_mean_outcome <- apply(post_pred_dists_outcome, 2, mean)
post_pred_95q_outcome <- apply(post_pred_dists_outcome, 2, quantile, probs = c(0.05, 0.95))

post_dens_list <- parallel::mclapply(1:n_samp, function(si){
  mcprint(paste0(si, " "))
  
  #generate and process code
  r_sim_code <- parse_Stan(stan_code, dat, samps, output_file = NA, 
                           sample_index = sample_inds[si], post_pred_sim = T, sim = F)
  
  #evaluate code
  stan_sim_env <- new.env()
  eval(parse(text = r_sim_code), envir = stan_sim_env)
  out_sim_sim <- stan_sim_env$out
  
  return(out_sim_sim)
}, mc.cores = 8)
post_pred_dens_outcome <- do.call(rbind, lapply(post_dens_list, function(post_pred_dens) post_pred_dens$x))
nlppm <- log(apply(post_pred_dens_outcome, 2, mean))
nlppm_col_ind <- nlppm - min(nlppm) 
nlppm_col_ind <- ceiling(nlppm_col_ind / max(nlppm_col_ind) * 100) + 1
colvals <- viridisLite::cividis(101)

#and plot it
png(filename = paste0(stan_figure_dir, "outcome_goodness-of-fit.png"), width = 800, height = 500, pointsize = 15)
par(mar = c(4,4,4,8))
plot(dat$x , post_pred_mean_outcome / dat$k, pch = 19, col = colvals[nlppm_col_ind], cex = 0.5, #dat$k / max(dat$k),
     xlim = c(0,1), ylim = c(0,1), xlab = "observed proportion", ylab = "posterior proportion", 
     main = paste0("Pearson's r = ", round(cor(dat$x / dat$k, post_pred_mean_outcome / dat$k), 3)))
segments(x0 = dat$x / dat$k, x1 = dat$x / dat$k, y0 = post_pred_95q_outcome[1,] / dat$k, y1 = post_pred_95q_outcome[2,] / dat$k, col = colvals[nlppm_col_ind])
abline(0,1,col=2,lty=2,lwd=2)
legend(x = "topleft", legend = c("posterior mean", "95% credible interval"), pch = c(19, NA), lty = c(NA,1))
add_continuous_legend(colvals, positions = (nlppm_col_ind-1) / 100, labels = nlppm, x = 1.075, y = 1, xpd = NA, main = "lppm", n_rect = 50, left_below = F)
dev.off()
