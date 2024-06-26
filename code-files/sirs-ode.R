



# Libs --------------------------------------------------------------------

library(pacman)
p_load(dplyr, ggplot2, tidyr, deSolve, stringr, purrr)

# currently loading datasets from n = ?
current_n = 'p4000'


# Load data from ABM ------------------------------------------------------

data_files <- list.files(path = paste0('data/', current_n), pattern = "*.csv", full.names = TRUE)

#' Function to read and process a single CSV file
#'
#' @param file the file path and name and extension e.g. 'data/p500/sim.csv'
#'
#' @return a single file read
process_file <- function(file) {
  read.csv(file) |>
    select(time, S = Susceptible, I = Infected.infectious, R = Recovered) |>
    mutate(P = I / (S + I + R),
           file = file) |>
    select(S, I, R, time, P, file)
}

abm_data <- map_dfr(data_files, process_file)


# Functions --------------------------------------------------------------


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
  
  nlls <- -dbinom(obsDat$numPos, obsDat$numSamp, prob = simDat$P, log = T)
  return(sum(nlls))
}


# nllikelihood(par =  c(c = 1, p = .2, gamma = .4, xi = .5),
#              obsDat=myDat, abm_data = abm_data) ## loglikelihood of the true parameters (which we usually never know)



# Fitting the model by MLE ------------------------------------------------

for (i in 1:length(data_files)) {
  message('Running simulation: ', i, ' out of : ', length(data_files))
  
  # using only the chosen simulation
  filtered_data = abm_data |> 
    filter(file == data_files[i]) |>
    select(-file)
  
  # sampling from the ABM: Here we choose to use all the abm data
  myDat = sampleEpidemic(filtered_data)
  
  # optimizing 
  optim.vals <- optim(par = c(c = 3, p = .4, gamma = .05, xi = .1)
                      , nllikelihood
                      , obsDat = myDat
                      , abm_data = filtered_data
                      , control = list(trace = 1, maxit = 1000)
                      , method = "Nelder-Mead")
  
  pars = optim.vals$par |> data.frame() |> t()
  
  if (i == 1) fullpar = pars
  else fullpar = rbind(fullpar, pars)
  
}



# Plotting the parameter values 
fullpar2 = fullpar |>
  data.frame() |>
  setNames(c('Contact rate', 'Probability of \ninfection', 'Rate of \nrecovery', 'Rate of \nre-infection')) |>
  pivot_longer(everything())

truepars = data.frame(name = c('Contact rate', 'Probability of \ninfection', 'Rate of \nrecovery', 'Rate of \nre-infection'),
                      value = c(NA, .35, 1/13, 1/7))

ggplot(fullpar2) +
  geom_density(aes(x = value)) + 
  facet_wrap(~name, scales = "free") +
  geom_vline(data = truepars, aes(xintercept = value, group = name), col = 'blue', size = 1) + 
  labs(title = 'Parameter values for N = 500', x = 'Parameter', y = 'Density') +
  theme_bw(base_line_size = 0) +
  theme(panel.spacing = unit(1, "lines"))



# Plotting one occurence
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
             parms = fullpar[i, ])
  
  # Convert to data frame for easier handling
  out <- as.data.frame(out)
  
  # plottiong to confirm
  plts[[i]] = ggplot(filtered_data, aes(x = time)) + 
    geom_point(aes(y = I), cex = .3) + 
    geom_line(data = out, aes(x = time, y = I), col = 'red') + 
    theme_classic()
  
}

do.call(grid.arrange, plts)

fullpar = fullpar |>
  data.frame()
rownames(fullpar) = 1:nrow(fullpar)
write.csv(fullpar, paste0('output/sirs-ode/', current_n, '-pars.csv'), row.names = F)











