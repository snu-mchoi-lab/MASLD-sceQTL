#########################################################
# power analysis for poisson regression
###################################################
# Function to simulate data and calculate power
simulate_eqtl_power_poisson <- function(sample_size, allele_freq, effect_size, alpha = 0.05) {
  # Parameters
  n <- sample_size                # Number of individuals
  f <- allele_freq                # Minor allele frequency
  beta1 <- effect_size            # Genotype effect size (log-scale)
  
  # Simulate genotype (0, 1, 2) based on allele frequency
  genotypes <- rbinom(n, 2, f)
  
  # Simulate covariates (e.g., age, sex)
  covariate1 <- rnorm(n, mean = 0, sd = 1)  # Continuous covariate
  covariate2 <- rbinom(n, 1, 0.5)           # Binary covariate
  covariate3 <- rnorm(n, mean = 0, sd = 1)  # Continuous covariate
  covariate4 <- rnorm(n, mean = 0, sd = 1)  # Continuous covariate
  covariate5 <- rnorm(n, mean = 0, sd = 1)  # Continuous covariate
  
  # Simulate gene expression as Poisson counts
  random_covariates = rnorm(6, mean=0, sd=1)
  log_lambda <- random_covariates[6] + beta1 * genotypes +random_covariates[1] * covariate1 - random_covariates[2] * covariate2 + random_covariates[3] * covariate3 + random_covariates[4] * covariate4 + random_covariates[5] *covariate5
  expression <- rpois(n, lambda = exp(log_lambda))
  
  # Fit Poisson regression model
  model <- glm(expression ~ genotypes + covariate1 + covariate2 + covariate3 + covariate4 + covariate5, family = poisson(link = "log"))
  
  # Extract test statistic and p-value for genotype
  summary_model <- summary(model)
  p_value <- coef(summary_model)["genotypes", "Pr(>|z|)"]
  
  # Return whether the effect is significant at the given alpha level
  return(p_value < alpha)
}

# Function to calculate power through simulation
calculate_power_poisson <- function(sample_size, allele_freq, effect_size, alpha = 0.05, num_sim = 100) {
  significant_results <- replicate(num_sim, simulate_eqtl_power_poisson(sample_size, allele_freq, effect_size, alpha))
  power <- mean(significant_results)
  return(power)
}


################################
#### simulation: 10000 times ###
################################
library(parallel)

# Define parameter grids
sample_sizes <- c(100, 200, 400, 800, 1600, 3200, 6400, 12800)
allele_freqs <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
effect_sizes <- c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6)
# effect_sizes <- c(0.01, 0.05, 0.1)

# Create a grid of all parameter combinations
param_grid <- expand.grid(sample_size = sample_sizes, 
                          allele_freq = allele_freqs, 
                          effect_size = effect_sizes)

# Function to calculate power for a given parameter combination
calculate_power_wrapper <- function(params) {
  sample_size_i <- params$sample_size
  allele_freq_i <- params$allele_freq
  effect_size_i <- params$effect_size
  
  # Log progress
  message(paste0(Sys.time(), "  :  ", sample_size_i, " / ", allele_freq_i, " / ", effect_size_i))
  
  # Perform power calculation
  power_result <- calculate_power_poisson(
    sample_size = sample_size_i, 
    allele_freq = allele_freq_i, 
    effect_size = effect_size_i, 
    alpha = 0.05 / 20000, 
    num_sim = 10000
  )
  
  # Return results as a data frame
  return(data.frame(
    sample_size = sample_size_i, 
    allele_freq = allele_freq_i, 
    effect_size = effect_size_i, 
    alpha = 0.05 / 20000, 
    num_sim = 10000, 
    power = power_result
  ))
}

# Use mclapply to process the parameter grid in parallel
results_list <- mclapply(seq_len(nrow(param_grid)), function(i) {
  calculate_power_wrapper(param_grid[i, ])
}, mc.cores = 20)

# Combine results into a single data frame
res <- do.call(rbind, results_list)
write.table(res, './poisson_regression_simulation_power_full_n10000.txt', sep='\t', row.names = F)


################################
## allele freq ~ Prob of 2 hets ##
################################
# Function to calculate the probability of at least 2 heterozygotes

calculate_het_probability <- function(cohort_size, allele_freq) {
  # Probability of heterozygote under Hardy-Weinberg
  prob_het <- 2 * allele_freq * (1 - allele_freq)
  # Binomial probabilities for X = 0 and X = 1
  p0 <- dbinom(0, cohort_size, prob_het)  # P(X = 0)
  p1 <- dbinom(1, cohort_size, prob_het)  # P(X = 1)
  # P(X >= 2) = 1 - P(X = 0) - P(X = 1)
  p_ge_2 <- 1 - p0 - p1
  return(p_ge_2)
}

cc=1
for(cohort_size in c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,80,90,100,120,140,160,180,200)){
  for(allele_freq in seq(0.05, 0.5, by=0.05)){
    two_het_prob <- calculate_het_probability(cohort_size, allele_freq)
    res_tmp = data.frame(cohort_size=cohort_size, allele_freq=allele_freq, probability=two_het_prob)
    if(cc==1) { prob_res = res_tmp }
    else{prob_res = rbind(prob_res, res_tmp)}
    cc = cc+1
  }
}
write.table(prob_res, './probability_of_two_hets.txt', sep='\t')
ggplot(prob_res, aes(x=cohort_size, y=probability, color=as.factor(allele_freq))) + geom_point() + geom_line() + theme_classic() +  geom_hline(yintercept=0.95, linetype='dashed') + scale_color_discrete(name="Allele frequency") + ylab("Probability")
ggsave('./probability_of_two_hets_50_0.95.pdf', width=5, height=3)

