

# This is the file containing code for analyzing parameters ---------------
# we only use this file to analyze the contents of output folder

data_files <- list.files(path = 'output/sirs-ode', pattern = "*.csv", full.names = TRUE)

params <- map_dfr(data_files, \(file) {
  read.csv(file) |> 
    mutate(N = str_extract_all(file, '[0-9]') |> 
             stri_join_list())
})

params = params |> mutate(partype = 'ode')

actual_params <- data.frame(c = rep(NA, 6),
                            p = c(.3, .3, .35, .4, .4, .4),
                            gamma = c(1/15, 1/15, 1/13, 1/15, 1/15, 1/15),
                            xi = c(1/10, 1/5, 1/7, 1/7, 1/7, 1/7),
                            N = c(500, (1:5)*1000),
                            partype = 'abm')

dd = rbind(params, actual_params)


# Visuaization ------------------------------------------------------------

# a visualization of a series of beta over values of N
dd |>
  filter(partype == 'ode') |>
  select(p, N) |>
  mutate(N = as.numeric(N)) |>
  ggplot() + 
  geom_boxplot(aes(x = factor(N), y = p)) +
  
  # adding the ABM actual points
  geom_point(data = dd |> filter(partype == 'abm') |> mutate(N = as.numeric(N)), 
             aes(x = factor(N), y = p, group = factor(N)), col = 'blue', cex = 4) +
  
  scale_y_continuous(labels = scales::percent_format()) + 
  labs(title = 'The probability of infection plot over a range of N',
       x = 'N (Initial population size)', y = 'Probability of infection',
       subtitle = 'The blue points denotes the true probability of infection \nused in generating ABM simulations') +
  theme_bw(base_line_size = 0) +
  theme(panel.spacing = unit(1, "lines"),
        legend.position = 'bottom') 
  
# a visualization of the density of betas for a given N
params |>
  select(p, N) |>
  mutate(N = as.numeric(N)) |> 
  ggplot() +
  geom_density(aes(x = p, fill = (N), group = N), alpha = .2) + 
  scale_x_continuous(labels = scales::percent_format()) + 
  
  scale_fill_gradient(low = "lightblue", high = "blue4") + 
  
  # facet_wrap(~N, scales = "free") +
  # geom_vline(data = truepars, aes(xintercept = value, group = name), col = 'blue', size = 1) + 
  labs(title = 'The density of P:', x = 'Probability of infection', y = 'Density',
       subtitle = 'The blue line denotes the true probability of infection used in simulations',
       fill = 'N: ') +
  theme_bw(base_line_size = 0) +
  theme(panel.spacing = unit(1, "lines"),
        legend.position = 'bottom')


# faceted plot
dd |>
  filter(partype == 'ode') |>
  select(p, N) |>
  mutate(N = factor(N, levels = c(500, (1:5)*1e3))) |> 
  ggplot() +
  geom_density(aes(x = p), alpha = .2, fill = 'gray', col = 'black') + 
  scale_x_continuous(labels = scales::percent_format()) + 
  
  facet_wrap(~N, scales = "free") +
  geom_vline(data = dd %>% filter(partype == 'abm') |> mutate(N = factor(N, levels = c(500, (1:5)*1e3))),
             aes(xintercept = p), color = 'blue', size = 1) + 
  
  labs(title = 'The density of P:', x = 'Probability of infection', y = 'Density',
       subtitle = 'The blue line denotes the true probability of infection used in simulations',
       fill = 'N: ') +
  theme_bw(base_line_size = 0) +
  theme(panel.spacing = unit(1, "lines"),
        legend.position = 'bottom')




