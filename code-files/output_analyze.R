

# This is the file containing code for analyzing parameters ---------------
# we only use this file to analyze the contents of output folder

data_files <- list.files(path = 'output/sirs-ode', pattern = "*.csv", full.names = TRUE)

params <- map_dfr(data_files, \(file) {
  read.csv(file) |> 
    mutate(N = str_extract_all(file, '[0-9]') |> 
             stri_join_list())
})


# Visuaization ------------------------------------------------------------

# a visualization of a series of beta over values of N
params |>
  select(p, N) |>
  mutate(N = as.numeric(N)) |>
  ggplot() + 
  geom_boxplot(aes(x = factor(N), y = p)) +
  geom_hline(aes(yintercept = .3), col = 'blue') + 
  scale_y_continuous(labels = scales::percent_format()) + 
  labs(title = 'The probability of infection plot over N',
       x = 'N (Initial population size)', y = 'Probability of infection',
       subtitle = 'The blue line denotes the true probability of infection used in simulations') +
  theme_bw(base_line_size = 0) +
  theme(panel.spacing = unit(1, "lines"),
        legend.position = 'bottom') 
  
# a visualization of the density of betas for a given N
params |>
  select(p, N) |>
  mutate(N = as.numeric(N)) |> 
  ggplot() +
  geom_density(aes(x = p, fill = factor(N)), alpha = .2) + 
  scale_x_continuous(labels = scales::percent_format()) + 
  # facet_wrap(~N, scales = "free") +
  # geom_vline(data = truepars, aes(xintercept = value, group = name), col = 'blue', size = 1) + 
  labs(title = 'The density of P:', x = 'Probability of infection', y = 'Density',
       subtitle = 'The blue line denotes the true probability of infection used in simulations',
       fill = 'N: ') +
  theme_bw(base_line_size = 0) +
  theme(panel.spacing = unit(1, "lines"),
        legend.position = 'bottom')

