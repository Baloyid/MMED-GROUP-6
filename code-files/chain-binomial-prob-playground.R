


# This is a playground to test out how the various probability of  --------
# infection incorporating the infected play out


ii =  seq(from = 100, to = 900, by = 50)
c = 3
prob_infect = .3
N = 1000

# definition 1 for prob of infection
pp = 1 - (1 - c*(ii/N)*prob_infect)**c

# definition 2
pp = 1 - (1 - prob_infect)**(c*ii/N)

# definition 3
1 - exp(-c*prob_infect * ii/S0)

plot(ii, pp, type = 'o')

# visualizing for contacts differences
dd = data.frame(
  I = seq(from = 100, to = 1000, by = 50),
  N = 1000
) |>
  mutate(
    c1 = 1 - (1 - 1*(I/N)*.3)**1,
    c2 = 1 - (1 - 2*(I/N)*.3)**2,
    c3 = 1 - (1 - 3*(I/N)*.3)**3,
    c4 = 1 - (1 - 4*(I/N)*.3)**4,
    c5 = 1 - (1 - 5*(I/N)*.3)**5,
    
    across(.cols = starts_with('c'), .fns = \(x) ifelse(x > 1, 1, x))
    
    # c1 = 1 - (1 - .3)**(1*I/N),
    # c2 = 1 - (1 - .3)**(2*I/N),
    # c3 = 1 - (1 - .3)**(3*I/N),
    # c4 = 1 - (1 - .3)**(4*I/N),
    # c5 = 1 - (1 - .3)**(5*I/N)
    
    # c1 = exp(1 * .3 * I/(N - I)) |> pnorm(),
    # c2 = exp(2 * .3 * I/(N - I)) |> pnorm(),
    # c3 = exp(3 * .3 * I/(N - I)) |> pnorm(),
    # c4 = exp(4 * .3 * I/(N - I)) |> pnorm(),
    # c5 = exp(5 * .3 * I/(N - I)) |> pnorm()
  ) |>
  pivot_longer(-c(I, N))


ggplot(dd) + 
  geom_line(aes(x = I, y = value, col = name)) + 
  geom_point(aes(x = I, y = value, col = name)) + 
  theme_classic()
