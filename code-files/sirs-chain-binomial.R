


# Libs --------------------------------------------------------------------

library(pacman)
p_load(dplyr, ggplot2, tidyr, deSolve, stringr)

# Load data from ABM ------------------------------------------------------

data_files <- list.files(path = paste0('data'), pattern = "*.csv", full.names = TRUE, recursive = T)
data_files = data_files[101:120]
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




# Fitting an SIR model by MLE ---------------------------------------------

#' chain binomial SIR function
#'
#' @param time gives the time we want our simulation to run from
#' @param state gives the state values for S, I, R
#' @param parameters include: prob_infection, prob_recovery, prob_reinfection, contact_rate
#'
#' @return list with dS, dI, dR
sirs_chainbinomial <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    N <- S + I + R
    
    prob_infect_full = 1 - (1 - c*(I/N)*prob_infect)**c
    
    # Differential equations
    infected_individuals = rbinom(size = last(S), n = 1, prob = prob_infect_full)
    recovered_individuals = rbinom(size = last(I), n = 1, prob = prob_recover)
    reinfected_individuals = rbinom(size = last(R), n = 1, prob = prob_reinfect)
    
    dS = reinfected_individuals - infected_individuals
    dI = infected_individuals - recovered_individuals
    dR = recovered_individuals - reinfected_individuals
    
    list(c(dS, dI, dR))
  })
  
}



#' The function for generating data given time points
#'
#' @param time gives the time we want our simulation to run from
#' @param state gives the state values for S, I, R
#' @param parameters include: prob_infection, prob_recovery, prob_reinfection, contact_rate
#'
#' @return dataset with time series of S, I, R
solve_chainbinomial <- function(time, initial_state, parameters) {
  states = as.list(initial_state)
  N <- last(states$S) + last(states$I) + last(states$R)
  
  
  
  prob_infect = 1 - (1 - parameters[4]*(last(states$I)/N)*parameters[1])**parameters[4] |> as.numeric()
  prob_recover = parameters[2] |> as.numeric()
  prob_reinfect = parameters[3] |> as.numeric()
  
  for (i in time) {
    
    # Differential equations
    infected_individuals = rbinom(size = last(states$S), n = 1, prob = prob_infect)
    recovered_individuals = rbinom(size = last(states$I), n = 1, prob = prob_recover)
    reinfected_individuals = rbinom(size = last(states$R), n = 1, prob = prob_reinfect)
    
    states$S = c(states$S, last(states$S) - infected_individuals + reinfected_individuals)
    states$I = c(states$I, last(states$I) + infected_individuals - recovered_individuals)
    states$R = c(states$R, last(states$R) + recovered_individuals - reinfected_individuals)
  }
  simDat <- list2DF(states) |>
    as.data.frame() |>
    mutate(P = I / (S + I + R))
  return(simDat)
}

# trial
# sirs_chainbinomial(time = 1:200, state = c(S=1000, I = 90, R = 30), parameters = parms)


#' Function for sampling from the ABM world dataset: returns either full or reduced sample
#'
#' @param simDat The data from the ABM 
#'
#' @return a dataset with time, positive cases, sample size, prevalence, lower and upper CI
sampleEpidemic <- function(simDat) {
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

# sample code
# myDat = sampleEpidemic(abm_data)



#' Function for computing the negative log likelihood
#'
#' @param parms the parameters of the SIRS model
#' @param obsDat the observation data i.e. a sample from the ABM data
#' @param abm_data # the actual ABM data
#'
#' @return numeric value: the sum of the negative log-likelihoods
nllikelihood <- function(parms, obsDat=myDat, abm_data = abm_data) {
  
  initial_state = (abm_data)[1, 1:3] |> unlist()
  time = seq(from = 1, to = max(abm_data$time), by = 1)
  
  simDat = solve_chainbinomial(time = time, 
                               initial_state = initial_state,
                               parameters = parms) |>
    mutate(P = I / (S + I + R))
  
  nlls <- -dbinom(obsDat$numPos, round(obsDat$numSamp, 0), prob = simDat$P, log = T)
  return(sum(nlls))
}

# sample code
# nllikelihood(parms = c(prob_infect = 0.3, prob_recover = 0.5, prob_reinfect = 0.4, c=2),
#              obsDat=myDat,
#              abm_data = abm_data) ## loglikelihood of the true parameters (which we usually never know)



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
  optim.vals <- optim(par = c(prob_infect=.3, prob_recover=.4, prob_reinfect=.3, c = 2)
                      , nllikelihood
                      , obsDat = myDat
                      , abm_data = filtered_data
                      , control = list(trace = 0, maxit = 1000)
                      , method = "Nelder-Mead")
  
  pars = optim.vals$par |> data.frame() |> t()
  
  if (i == 1) fullpar = pars
  else fullpar = rbind(fullpar, pars)
  
}


# creating the datasets showing population size and neighbourhood type
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
  mutate(`Probability of \ninfection` = ifelse(`Probability of \ninfection` < 0, 0, `Probability of \ninfection`))

fullpar2 |>
  pivot_longer(-c(`Population size`, `Neighbourhood`))

# boxplot
fullpar2 |> ggplot() + 
  geom_boxplot(aes(x = as.factor(`Population size`),
                   y = `Probability of \ninfection`,
                   fill = Neighbourhood)) + 
  geom_hline(aes(yintercept = .4), col = 'blue') + 
  labs(title = 'Probability of infection over population size ranges',
       x = 'Population size', y = 'Infection probability',
       fill = 'Neighbourhood type: ') +
  theme_bw(base_line_size = 0) + 
  theme(legend.position = 'bottom')

# density
fullpar2 |> ggplot() + 
  geom_density(aes(#x = as.factor(`Population size`),
    x = `Probability of \ninfection`,
    fill = as.factor(`Population size`)),
    alpha = .2) + 
  facet_wrap(~Neighbourhood) +
  geom_vline(aes(xintercept = .4), col = 'blue') + 
  labs(title = 'Probability of infection over population size ranges',
       x = 'Population size', y = 'Infection probability',
       fill = 'Neighbourhood type: ') +
  theme_bw(base_line_size = 0) + 
  theme(legend.position = 'bottom')









# Plotting the parameter values 
fullpar2 = fullpar |>
  data.frame() |>
  setNames(c('Contact rate', 'Probability of \ninfection', 'Rate of \nrecovery', 'Rate of \nre-infection')) |>
  pivot_longer(everything())

truepars = data.frame(name = c('Contact rate', 'Probability of \ninfection', 'Rate of \nrecovery', 'Rate of \nre-infection'),
                      value = c(NA, .3, 1/15, 1/5))

ggplot(fullpar2) +
  geom_density(aes(x = value)) + 
  facet_wrap(~name, scales = "free") +
  geom_vline(data = truepars, aes(xintercept = value, group = name), col = 'blue', size = 1) + 
  labs(title = 'Parameter values for N = 1,000', x = 'Parameter', y = 'Density') +
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
  out <- solve_chainbinomial(initial_state = (filtered_data)[1, 1:3] |> unlist(),
                             time = 1:nrow(filtered_data),
                             parameters = fullpar[i, ]) |>
    data.frame() |>
    mutate(time = row_number(),
           N = S + I + R)
  
  # getting binomial confidence intervals using I and N
  bci = with(out, DescTools::BinomCI(x = I, n = N, method = 'wilson') |>
               data.frame() |>
               mutate(
                 across(
                   .cols = everything(),
                   .fns = \(x) x * N
                 )
               )
  )
  
  
  plts[[i]] = ggplot(filtered_data, aes(x = time)) + 
    geom_point(aes(y = I), cex = .3) + 
    geom_line(data = cbind(out, bci),
              aes(x = time, y = I), col = 'red') +
    geom_ribbon(data = cbind(out, bci),
                aes(x = time, ymin = lwr.ci, ymax = upr.ci), 
                fill = 'red', alpha = .4) +
  theme_classic()
  
}

do.call(grid.arrange, plts)

# storing the results of the parameters
fullpar = fullpar |>
  data.frame()
rownames(fullpar) = 1:nrow(fullpar)
write.csv(fullpar, paste0('output/sirs-chain-binomial/', current_n, '-pars.csv'), row.names = F)

