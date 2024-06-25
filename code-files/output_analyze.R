

# This is the file containing code for analyzing parameters ---------------

data_files <- list.files(path = 'output/sirs-ode', pattern = "*.csv", full.names = TRUE)

#' Function to read and process a single CSV file
#'
#' @param file the file path and name and extension e.g. 'data/p500/sim.csv'
#'
#' @return a single file read
process_file <- function(file) read.csv(file) |> mutate(N = str_extract_all(file, '[0-9]') |> stri_join_list())

params <- map_dfr(data_files, process_file)


# Visuaization ------------------------------------------------------------

params |>
  select(p, N) |>
  mutate(N = as.numeric(N)) |>
  ggplot() + 
  geom_boxplot(aes(x = factor(N), y = p)) +
  geom_hline(aes(yintercept = .3), col = 'blue') + 
  labs(title = 'The probability of infection plot over N',
       x = 'N (Initial population size)', x = 'Probability of infection',
       subtitle = 'The blue line denotes the true probability of infection used in simulations') +
  theme_classic() 
  

params |>
  select(p, N) |>
  mutate(N = as.numeric(N)) |> 
  ggplot() +
  geom_density(aes(x = p, fill = factor(N)), alpha = .2) + 
  # facet_wrap(~N, scales = "free") +
  # geom_vline(data = truepars, aes(xintercept = value, group = name), col = 'blue', size = 1) + 
  labs(title = 'The density of P:', x = 'Probability of infection', y = 'Density',
       subtitle = 'The blue line denotes the true probability of infection used in simulations',
       fill = 'N: ') +
  theme_bw(base_line_size = 0) +
  theme(panel.spacing = unit(1, "lines"),
        legend.position = 'bottom')

