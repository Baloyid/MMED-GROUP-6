

# Load necessary libraries
library(sf)
library(spdep)
library(ggplot2)
library(dplyr)



# Uni neighbourhood -------------------------------------------------------


# Create a simple grid of spatial polygons (e.g., 3x3 grid)
n <- 3
grid <- st_make_grid(st_bbox(c(xmin = 0, xmax = n, ymin = 0, ymax = n)), n = n) %>%
  st_as_sf() %>%
  mutate(id = 1:(n*n))

# Compute Queen contiguity matrix
nb <- poly2nb(grid, queen = TRUE)

# Find the middle box id (assuming n is odd)
middle_id <- (n^2 + 1) / 2

# Get neighbors of the middle box
middle_neighbors <- nb[[middle_id]]

# Shade the middle box and its neighbors
grid <- grid %>%
  mutate(fill_color = case_when(
    id == middle_id ~ "darkblue", # Middle box
    id %in% middle_neighbors ~ "white", # Neighboring boxes
    TRUE ~ "white" # Other boxes
  ))

# Plot the grid and shade the specified regions
p1 = ggplot() +
  geom_sf(data = grid, aes(fill = fill_color), color = 'black') +
  scale_fill_identity() +
  theme_void() +
  labs(title = "Uni neighbourhood",
       subtitle = "")





# Rook neighbourhood ------------------------------------------------------


# Create a simple grid of spatial polygons (e.g., 3x3 grid)
n <- 3
grid <- st_make_grid(st_bbox(c(xmin = 0, xmax = n, ymin = 0, ymax = n)), n = n) %>%
  st_as_sf() %>%
  mutate(id = 1:(n*n))

# Compute Rook contiguity matrix
nb <- poly2nb(grid, queen = FALSE)

# Find the middle box id (assuming n is odd)
middle_id <- (n^2 + 1) / 2

# Get neighbors of the middle box
middle_neighbors <- nb[[middle_id]]

# Shade the middle box and its neighbors
grid <- grid %>%
  mutate(fill_color = case_when(
    id == middle_id ~ "darkblue", # Middle box
    id %in% middle_neighbors ~ "lightblue", # Neighboring boxes
    TRUE ~ "white" # Other boxes
  ))

# Plot the grid and shade the specified regions
p2 = ggplot() +
  geom_sf(data = grid, aes(fill = fill_color), color = 'black') +
  scale_fill_identity() +
  theme_void() +
  labs(title = "Four neighbourhood")





# Queen neighbourhood -----------------------------------------------------


# Load necessary libraries
library(sf)
library(spdep)
library(ggplot2)
library(dplyr)

# Create a simple grid of spatial polygons (e.g., 3x3 grid)
n <- 3
grid <- st_make_grid(st_bbox(c(xmin = 0, xmax = n, ymin = 0, ymax = n)), n = n) %>%
  st_as_sf() %>%
  mutate(id = 1:(n*n))

# Compute Queen contiguity matrix
nb <- poly2nb(grid, queen = TRUE)

# Find the middle box id (assuming n is odd)
middle_id <- (n^2 + 1) / 2

# Get neighbors of the middle box
middle_neighbors <- nb[[middle_id]]

# Shade the middle box and its neighbors
grid <- grid %>%
  mutate(fill_color = case_when(
    id == middle_id ~ "darkblue", # Middle box
    id %in% middle_neighbors ~ "lightblue", # Neighboring boxes
    TRUE ~ "white" # Other boxes
  ))

# Plot the grid and shade the specified regions
p3 = ggplot() +
  geom_sf(data = grid, aes(fill = fill_color), color = 'black') +
  scale_fill_identity() +
  theme_void() +
  labs(title = "Eight neighbourhood",
       subtitle = "")


grid.arrange(p1, p2, p3, nrow = 1)

