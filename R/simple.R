# simple linear regression
#
# y ~ mx+n+epsilon
# y ~ N(mx+n, sigma)

# true data generating process
n <- 1
m <- 2
sigma <- 1

# done a series of experiments and repeated it 4 times
x <- rep.int(x = seq(from = 0, to = 20, by = 1),
             times = 4)
y <- sapply(x, function(i){
  rnorm(n = 1,
        mean = m * i + n,
        sd = sigma)
})

df <- dplyr::tibble(x = x,
                    y = y)

df <- dplyr::mutate(df, n_tilde = )

# prior


# likelihood

# sampler