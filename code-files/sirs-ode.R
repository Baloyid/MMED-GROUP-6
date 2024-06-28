



# Libs --------------------------------------------------------------------

library(pacman)
p_load(dplyr, ggplot2, tidyr, deSolve, stringr, purrr, broom)

# Load data from ABM ------------------------------------------------------

data_files <- list.files(path = paste0('data'), pattern = "*.csv", full.names = TRUE, recursive = T)

#' Function to read and process a single CSV file
#'
#' @param file the file path and name and extension e.g. 'data/p500/sim.csv'
#'
#' @return a single file read
process_file <- function(file) {
  # Read the file as a text string
  file_content <- read_file(file)
  file_content <- gsub(";", ",", file_content)
  # print(file)
  
  read_csv(I(file_content)) |>
    select(time, S = Susceptible, I = `Infected-infectious`, R = Recovered) |>
    mutate(P = I / (S + I + R),
           file = file) |>
    select(S, I, R, time, P, file)
}

abm_data <- map_dfr(data_files, process_file)

abm_data = abm_data |>
  mutate(
    pop = str_split(file, '/') |> lapply(FUN = \(x) x[2]) |> unlist() |> str_replace_all('p', '') |> as.numeric(),
    nb = str_split(file, '/') |> lapply(FUN = \(x) x[3]) |> unlist()
  ) |>
  # filter(pop == 2000 & nb == 'no-neighbours') |>
  select(-pop, -nb)

# Functions --------------------------------------------------------------

#' Function for computing binomial confidence intervals, and point estimates
#'
#' @param proportions the binomial proportions
#'
#' @return a dataframe with mean, median and confidence intervals
calculate_summary_stats <- function(proportions) {
  
  mean_proportion <- mean(proportions)
  median_proportion <- median(proportions)
  
  se <- sd(proportions) / sqrt(length(proportions))
  crit <- qt(0.975, df = length(proportions) - 1)
  margin_of_error <- crit * se
  
  # Calculate the confidence interval
  ci_lower <- mean_proportion - margin_of_error
  ci_upper <- mean_proportion + margin_of_error
  
  # Create a data frame with the results
  results <- data.frame(
    mean = mean_proportion,
    median = median_proportion,
    conf.low = ci_lower,
    conf.high = ci_upper
  )
  
  return(results)
}


#' Define the SIRS model with the beta parametrization
#'
#' @param time # a vector of time points at which to evaluate 
#' @param state a vector with numbr of S, I, R at a given point
#' @param parameters a vector of parameters: prob_infect, prob_recover, prob_reinfect, contact_rate
#'
#' @return a list with the dS, dI, dR
sirs_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S + I + R
    
    # Differential equations
    dS <- xi * R - beta * S * I / N
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I - xi * R
    
    list(c(dS, dI, dR))
  })
}

#' This is the SIRS model with c*P implementation
#'
#' @param time # a vector of time points at which to evaluate 
#' @param state a vector with numbr of S, I, R at a given point
#' @param parameters a vector of parameters: prob_infect, prob_recover, prob_reinfect, contact_rate
#'
#' @return a list with the dS, dI, dR
sirs_model_cp <- function(time, state, par) {
  
  with(as.list(c(state, par)), {
    N <- S + I + R
    beta = c*p
    
    # Differential equations
    dS <- xi * R - beta * S * I / N
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I - xi * R
    
    list(c(dS, dI, dR))
  })
}

#' Function for sampling from the ABM world dataset: returns either full or reduced sample
#'
#' @param simDat The data from the ABM 
#'
#' @return a dataset with time, positive cases, sample size, prevalence, lower and upper CI
sampleEpidemic <- function(simDat # Simulated "data" which we treat as real 
                           , sampleDates = seq(from = 1, to = 200, by = 1) # Sample every 3 years 
                           , numSamp = rep(1000, length(sampleDates)) # Number of individuals sampled at each time point
){
  prev_at_sample_times <- simDat |> pull(P)
  samp_size = lag(simDat$S) |> zoo::na.locf(fromLast = T)
  numSamp = samp_size
  
  numPos <- rbinom(length(numSamp), round(numSamp, 0), prev_at_sample_times)
  
  lci = qbeta(.025, numPos, round(numSamp, 0) - numPos + 1)
  uci = qbeta((1-.025), numPos + 1, round(numSamp, 0) - numPos)
  # lci <- mapply(function(x,n) binom.test(x,n)$conf.int[1], x = numPos, n = round(numSamp, 0))
  # uci <- mapply(function(x,n) binom.test(x,n)$conf.int[2], x = numPos, n = round(numSamp, 0)) 
  
  return(data.frame(time = 1:length(numPos), numPos, numSamp, sampPrev =  numPos/numSamp,
                    lci = lci, uci = uci))
}

# myDat = sampleEpidemic(abm_data)

#' Function for computing the negative log likelihood
#'
#' @param parms the parameters of the SIRS model
#' @param obsDat the observation data i.e. a sample from the ABM data
#' @param abm_data # the actual ABM data
#'
#' @return numeric value: the sum of the negative log-likelihoods
nllikelihood <- function(par = par, obsDat=myDat, abm_data = abm_data) {
  
  initial_state = (abm_data)[1, 1:3] |> unlist()
  time = seq(from = 0, to = max(abm_data$time), by = 1)
  
  simDat <- ode(y = initial_state, times = time, func = sirs_model_cp, par = par) |> # sirs_model
    as.data.frame() |>
    mutate(P = I / (S + I + R))
  
  # Ensure that simDat$P contains valid probabilities, non valid are hugely penalized
  if(any(is.na(simDat$P) | simDat$P < 0 | simDat$P > 1)) {
    return(1000000)  # Return a large value to penalize invalid parameter sets
  }
  
  nlls <- -dbinom(obsDat$numPos, obsDat$numSamp, prob = simDat$P, log = T)
  return(sum(nlls))
}


# nllikelihood(par =  c(c = 1, p = .2, gamma = .4, xi = .5),
#              obsDat=myDat, abm_data = abm_data) ## loglikelihood of the true parameters (which we usually never know)



# Fitting the model by MLE ------------------------------------------------

# data_files = abm_data |> pull(file) |> unique()
for (i in 1:length(data_files)) {
  message('Running simulation: ', i, ' out of : ', length(data_files))
  
  # using only the chosen simulation
  filtered_data = abm_data |> 
    filter(file == data_files[i]) |>
    select(-file)
  
  # sampling from the ABM: Here we choose to use all the abm data
  myDat = sampleEpidemic(filtered_data)
  
  # optimizing to get most appropriate parameters 
  optim.vals <- optim(par = c(c = 3, p = .3, gamma = .05, xi = .2)
                      , nllikelihood
                      , obsDat = myDat
                      , abm_data = filtered_data
                      , lower = c(1, rep(0, 3))
                      , upper = c(5, rep(1, 3))
                      , control = list(trace = 0, maxit = 1000)
                      , method = "L-BFGS-B")

  
  pars = optim.vals$par |> data.frame() |> t()
  
  if (i == 1) fullpar = pars
  else fullpar = rbind(fullpar, pars)
  
}

# Saving results back in memory -------------------------------------------

fullpar2 = fullpar |>
  data.frame(check.names = F)
rownames(fullpar) = 1:nrow(fullpar)
write.csv(fullpar, paste0('output/sirs-ode/', 'full-pars.csv'), row.names = F)


# Visualizations ----------------------------------------------------------

# creating the datasets showing population size and neighbourhood type
# this adds the population and neighbourhood structure to the parameters dataframe
pop_params = abm_data |>
  select(file) |>
  distinct() |>
  mutate(
    pop = str_split(file, '/') |> lapply(FUN = \(x) x[2]) |> unlist() |> str_replace_all('p', '') |> as.numeric(),
    nb = str_split(file, '/') |> lapply(FUN = \(x) x[3]) |> unlist()
  ) |>
  select(-file)
  

# Plotting the parameter values 
fullpar2 = cbind(fullpar, pop_params) |>
  data.frame() |>
  setNames(c('Contact rate', 'Probability of \ninfection', 'Rate of \nrecovery', 'Rate of \nre-infection', 'Population size', 'Neighbourhood'))
fullpar2 = fullpar2 |>
  mutate(`Probability of 
infection` = ifelse(`Probability of 
infection` < 0, 0, `Probability of 
infection`))


# boxplot of results
(p1 = fullpar2 |> ggplot() + 
  geom_boxplot(aes(x = as.factor(`Population size`),
                   y = `Probability of \ninfection`,
                   fill = Neighbourhood)) + 
  geom_hline(aes(yintercept = .4), col = 'blue') + 
  facet_wrap(~Neighbourhood, scales = 'free') +
  labs(title = 'Probability of infection over population size ranges',
       x = 'Population size', y = 'Infection probability',
       fill = 'Neighbourhood: ') +
  theme_bw(base_line_size = 0) + 
  theme(legend.position = 'bottom'))



# The confidence intervals plot ----------------------------

(p2 = fullpar2 %>%
  group_by(`Population size`, Neighbourhood) %>%
  summarise(
    calculate_summary_stats(`Probability of \ninfection`)
  ) %>%
  ungroup() %>%
  ggplot(aes(x = mean, y = factor(`Population size`), col = factor(Neighbourhood))) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) + 
   facet_wrap(~Neighbourhood, scales = 'free') +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, position = position_dodge(width = 0.7)) +
  geom_vline(xintercept = 0.4, linetype = "dashed", color = "blue") +
  theme_classic() +
  labs(
    x = "Median and 95% CI\nfor the probability of infection",
    y = "Population size",
    col = 'Neighbourhood type: ',
    title = "Coverage for the probability of infection parameter"
  ) +
  theme(legend.title = element_blank(),
        legend.position = 'bottom'))

grid.arrange(p1, p2, nrow = 1)

# true parameters used in ABM simulation
# truepars = data.frame(name = c('Contact rate', 'Probability of \ninfection', 'Rate of \nrecovery', 'Rate of \nre-infection'),
#                       value = c(NA, .4, 1/7, 1/14))
# 
# ggplot(fullpar2 |> filter(name != 'Contact rate')) +
#   geom_density(aes(x = value, col = Neighbourhood)) + 
#   facet_wrap(~name, scales = "free", nrow = 2) +
#   geom_vline(data = truepars |> filter(name != 'Contact rate'),
#              aes(xintercept = value, group = name), col = 'blue', size = 1) + 
#   labs(title = 'Parameter values for N = 4000', x = 'Parameter', y = 'Density') +
#   theme_bw(base_line_size = 0) +
#   theme(panel.spacing = unit(1, "lines"))
# 
# # a boxplot
# ggplot(fullpar2 |> filter(name != 'Contact rate')) +
#   geom_boxplot(aes(x = name, y = value, col = Neighbourhood)) + 
#   facet_wrap(~name, scales = "free", nrow=2) +
#   geom_point(data = truepars |> filter(name != 'Contact rate'),
#              aes(y = value, x = name), col = 'black', size = 3) + 
#   labs(title = 'Parameter values for N = 4000', x = 'Parameter', y = 'Density') +
#   theme_bw(base_line_size = 0) +
#   theme(panel.spacing = unit(1, "lines"))



# Plots for transmission dynamics: I(t) -----------------------------------

# a list to hold all plots
plts = list()

for (i in 1:length(data_files)) {
  
  # using only the chosen simulation
  filtered_data = abm_data |> 
    filter(file == data_files[i]) |>
    select(-file)
  
  
  # Solving the system of equations: SIRS-ODE using the fitted parameters
  # using a single set of parameter values
  out <- ode(y = (filtered_data)[1, 1:3] |> unlist(),
             times = 1:nrow(filtered_data),
             func = sirs_model_cp,
             parms = fullpar[i, ]) |>
    data.frame() |>
    mutate(time = row_number(),
           N = S + I + R)
  
  # getting binomial confidence intervals using I and N
  bci.out = with(out, DescTools::BinomCI(x = I, n = N, method = 'wilson') |>
                   data.frame() |>
                   mutate(
                     across(
                       .cols = everything(),
                       .fns = \(x) x * N
                     )
                   )
  )
  
  # using actual parameters, 
  # the c parameter distribution was obtained in a small simulation done outside this main script, as c is a nuisance parameter
  actual.out = ode(y = (filtered_data)[1, 1:3] |> unlist(),
                   times = 1:nrow(filtered_data),
                   func = sirs_model_cp,
                   parms = c(c = runif(1, min = 0.13787, max = 0.18683), p = 0.4, gamma = 1/7, xi = 1/14)) |>
    data.frame() |>
    mutate(time = row_number(),
           N = S + I + R)
  
  # confidence intervals for actual parameters
  bci.actual = with(actual.out, DescTools::BinomCI(x = I, n = N, method = 'wilson') |>
                      data.frame() |>
                      mutate(
                        across(
                          .cols = everything(),
                          .fns = \(x) x * N
                        )
                      )
  )
  
  # storing the simulated data
  if (i == 1) simdatfull = cbind(out, bci.out) |> mutate(file = data_files[i])
  else simdatfull = rbind(simdatfull, cbind(out, bci.out) |> mutate(file = data_files[i]))
  
  if (i == 1) full.actual = cbind(actual.out, bci.actual) |> mutate(file = data_files[i])
  else full.actual = rbind(full.actual, cbind(actual.out, bci.actual) |> mutate(file = data_files[i]))
  
  plts[[i]] = ggplot(filtered_data, aes(x = time)) +
    geom_point(aes(y = I), cex = .3) +
    geom_line(aes(y = I), cex = .3, alpha = .3, col = 'black') +
    geom_line(data = cbind(out, bci.out),
              aes(x = time, y = I), col = 'red') +
    geom_ribbon(data = cbind(out, bci.out),
                aes(x = time, ymin = lwr.ci, ymax = upr.ci),
                fill = 'red', alpha = .2) +
    
    geom_line(data = cbind(actual.out, bci.actual),
              aes(x = time, y = I), col = 'green') +
    geom_ribbon(data = cbind(actual.out, bci.actual),
                aes(x = time, ymin = lwr.ci, ymax = upr.ci),
                fill = 'green', alpha = .2) +

    labs(x = 'Time (t)', y = 'I(t)') +
  theme_classic()
  
}

p0 = do.call(grid.arrange, plts)



# plotting combined datasets for given neighbourhood and given pop size

# dataset using calibrated parameters
simdatfull2 = simdatfull |>
  mutate(
         pop = str_split(file, '/') |> lapply(FUN = \(x) x[2]) |> unlist() |> str_replace_all('p', '') |> as.numeric(),
         nb = str_split(file, '/') |> lapply(FUN = \(x) x[3]) |> unlist()
         ) |>
  select(-file)

simdatfull3 = simdatfull2 |>
  group_by(pop, nb, time) |>
  summarise(
    across(
      .cols = everything(),
      .fns = mean
    )
  )

# data using actual parameters
full.actual2 = full.actual |>
  mutate(
    pop = str_split(file, '/') |> lapply(FUN = \(x) x[2]) |> unlist() |> str_replace_all('p', '') |> as.numeric(),
    nb = str_split(file, '/') |> lapply(FUN = \(x) x[3]) |> unlist()
  ) |>
  select(-file)

full.actual3 = full.actual2 |>
  group_by(pop, nb, time) |>
  summarise(
    across(
      .cols = everything(),
      .fns = mean
    )
  )



fd_temp = abm_data |>
  mutate(
    pop = str_split(file, '/') |> lapply(FUN = \(x) x[2]) |> unlist() |> str_replace_all('p', '') |> as.numeric(),
    nb = str_split(file, '/') |> lapply(FUN = \(x) x[3]) |> unlist()
  ) |>
  select(-file) |>
  merge(simdatfull3 |> select(pop, nb, time, est, lwr.ci, upr.ci) |> mutate(type = 'calibrated'),
        by = c('time', 'pop', 'nb'), all.x = T) |>
  merge(full.actual3 |> select(pop, nb, time, est, lwr.ci, upr.ci) |> mutate(type = 'Actual'),
        by = c('time', 'pop', 'nb'), all.x = T) 


# combined plot for any one simulation run
p1=ggplot() + 
  geom_point(data = fd_temp, aes(x = time, y = I), cex = .7) +
  geom_line(data = simdatfull3, aes(x = time, y = I), col = 'red') +
  geom_ribbon(data = simdatfull3, aes(x = time, ymin = lwr.ci, ymax = upr.ci),
              fill = 'red', alpha = .2) +
  
  geom_line(data = full.actual3, aes(x = time, y = I), col = 'green') +
  geom_ribbon(data = full.actual3, aes(x = time, ymin = lwr.ci, ymax = upr.ci),
              fill = 'green', alpha = .2) +
  
  # facet_wrap(nb ~ pop, scales = 'free', ncol = 6) + 
  labs(x = 'Time', y = 'Number of infected individuals',
       subtitle = 'The green curves represent SIRS models fitted with actual parameters\nThe red curves represent SIRS models fitted with calibrated (MLE) parameters') +
  ggtitle('Calibration results for SIRS-ODE models') + 
  theme_classic()

# you can run this
# p0 = do.call(grid.arrange, plts)


small_plots <- wrap_plots(plts, ncol = 3)  # Adjust ncol to arrange plots as desired
layout <- (p1 | small_plots) + plot_layout(widths = c(1, 1))  # Adjust widths to control the relative size of p1
print(layout)

