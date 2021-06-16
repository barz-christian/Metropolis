Metropolis and Gibbs sampling
================
Christian Barz
25 5 2021

-   [1 Memo](#memo)
-   [2 True data generating process](#true-data-generating-process)
-   [3 model / theory](#model--theory)
    -   [3.1 model](#model)
    -   [3.2 derivation of the formulas for the
        likelihood](#derivation-of-the-formulas-for-the-likelihood)
-   [4 implementation](#implementation)
    -   [4.1 likelihood](#likelihood)
    -   [4.2 prior](#prior)
    -   [4.3 posterior](#posterior)
    -   [4.4 MH](#mh)
    -   [4.5 test](#test)
-   [5 diagnostics](#diagnostics)
    -   [5.1 traceplots warmup](#traceplots-warmup)
    -   [5.2 traceplots warmup discarded](#traceplots-warmup-discarded)

# 1 Memo

# 2 True data generating process

We generated data in a way one would do it in a lab, i.e. run an
experiment multiple times.

``` r
n <- 1
m <- 2
sigma <- 10.

# done a series of experiments and repeated it 4 times
x <- rep.int(x = seq(from = 0, to = 20, by = 1),
             times = 4)
y <- sapply(x, function(i){
  rnorm(n = 1,
        mean = m * i + n,
        sd = sigma)
})

# system matrix, design
df <- dplyr::tibble(x = x,
                    intercept = 1)
X <- as.matrix(df)
# model y = A*x, mit A hat hat eine konstante spalte 1 für den intercept
parameter <- c(m,n)
sd <- sigma
```

# 3 model / theory

## 3.1 model

We assume that our observation *y* depends linear on *x* and each
measurement has a normal distributed error, i.e. each observation can be
written as

*y*<sub>*i*</sub> = *m* ⋅ *x**i* + *n* + *ϵ*<sub>*i*</sub>,

where *ϵ*<sub>*i*</sub> ∼ *N*(0, *σ*<sup>2</sup>).

This can be written in the following equivalent form:

*y*<sub>*i*</sub> ∼ *N*(*m* ⋅ *x*<sub>*i*</sub> + *n*, *σ*<sup>2</sup>).

One benefit of the last notation is, that we directly see (the) 3
parameters of our model:

-   *m* the slope of the regression line,
-   *n* the intercept of the regression line,
-   *σ*<sup>2</sup> the variance of the error term (corresponds to the
    amount of error of each measurement).

We do a Bayesian approach to infer the distribution of the regression
parameters *n*, *m*. Hence we have to specify a prior on the parameters.
For this example we use:

**!!Include PRIORS here!!**

This way we directly see the hyperparameters of the model.

## 3.2 derivation of the formulas for the likelihood

We are going to derive step by step a term that is proportional to the
joint posterior distribution of *n*, *m* given *y*. But we do not only
want to derive a formula; **no** we also want to derive the sufficient
and necessary conditions for the formula, which are often not mentioned.

$$
\\begin{align}
p(m,n\|y) & =      & \\frac{p(y\|m,n) p(m,n)}{p(y)}\\\\ 
         &\\propto & p(y\|m,n)  p(m,n) \\\\
         &=       & p(y\|m,n)  p(m) p(n) 
\\end{align}
$$

Here the last line the last line requires that the parameters *m*, *n*
are statistically independent, i.e. the joint probability is equal to
the product of the marginal distributions:

*p*(*m*, *n*) = *p*(*m*)*p*(*n*)

In general the quantity *y* is a matrix, but for simplicity we restrict
to the case when *y* = (*y*<sub>1</sub>, …, *y*<sub>*n*</sub>) is a
vector. In this case the posterior is proportional to

$$
\\begin{align}
p(m,n\|y)
&    =    p(m,n\|y\_1,\\ldots,y\_n)\\\\
&\\propto  p(y\_1,\\ldots,y\_n\|m,n) p(m) p(n) \\\\
& =       \\prod\_{i=1}^{n} p(y\_i\|m,n) p(m) p(n)
\\end{align}
$$

In the last step we used the assumption that the measurements are
independently and identically distributed.

In order to avoid numerical underflow we apply the natural logarithm and
get:

### 3.2.1 log likelihood

$$
\\begin{align}
\\log p(y\|m,n) &= \\sum\_{i=1}^{n} \\log p(y\_i\|m,n)\\\\
              &= \\sum \\log p\_{normal}(y\_i-m \\cdot x\_i + n,\\sigma^2)
\\end{align}
$$

### 3.2.2 log prior

### 3.2.3 log posteriori

Because we have used conjugated priors we know the closed formof the
posterior and can compare the true posterior with the approximation we
get from our Metropolis implemenation.

# 4 implementation

## 4.1 likelihood

``` r
likelihood <- function(X,y, parameter, sd = 3){
  yhat <- X %*% parameter
 
  likelihood_per_data_point <- dnorm(y, 
                                     mean = yhat,
                                     sd = sd,
                                     log = TRUE)
  sum_of_likelihoods <- sum(likelihood_per_data_point)
  
  
  return(sum_of_likelihoods)
}
```

### 4.1.1 plot likelihood

## 4.2 prior

``` r
prior <- function(parameter){
  log_prob_parameter <- dnorm(parameter, sd = 3, log = TRUE)
  sum(log_prob_parameter)
}
```

## 4.3 posterior

``` r
posterior <- function(parameter, X,y){
  likelihood(X = X, y = y, parameter = parameter) + prior(parameter)
}
```

## 4.4 MH

**ingredients**

-   proposal distribution
-   random initial value

### 4.4.1 proposal function

``` r
proposal_function <- function(parameter){
  new_parameter <- parameter
  for (i in 1:NROW(parameter)) {
    new_parameter[i] <- rnorm(n = 1, 
                              mean = parameter[i], 
                              sd = 1)
  }
  
  new_parameter
}
```

### 4.4.2 acceptance

``` r
MH <- function(initial_value, iterations){
  chain <- array(dim = c(iterations+1, 2))
  
  chain[1,] <- initial_value
  for (i in 1:iterations) {
    
    candidate <- proposal_function(chain[i,])
    
    acceptance <- exp(posterior(X = X,
                                y= y,
                                parameter = candidate) - 
                        posterior(X = X, 
                                  y = y, 
                                  parameter = chain[i,])
                      )
   
    if(runif(1) < acceptance){
      chain[i+1,] <- candidate
    }else{
      chain[i+1,] <- chain[i,]
    }
  }
  return(chain)
}
```

## 4.5 test

recall our true data generating process *y* = *m* ⋅ *x* + *n* and the
values of the parameters

``` r
m
n
```

``` r
startvalue = c(10,10)
iterations <- 200
warmup_length <- round(iterations/10, digits = 0)

chain = MH(startvalue, iterations)

posterior_samples <- tibble(m = chain[,1],
             n = chain[,2]) %>%
  mutate(id = 1:nrow(.)) %>%
  mutate(warmup = ifelse(id <= warmup_length,TRUE,FALSE))

# plot parameters with warmup discarded
posterior_samples %>%
  filter(!warmup) %>%
  ggplot() + geom_histogram(aes(m))


posterior_samples %>%
  filter(!warmup) %>%
  ggplot() + geom_histogram(aes(n))
```

# 5 diagnostics

For simplicity we restrict to visual checks, precisely to traceplots.

## 5.1 traceplots warmup

Traceplot of the slope parameter `n`

``` r
posterior_samples %>%
  filter(warmup) %>%
  ggplot(aes(x=id, y = m)) + geom_point()
```

Traceplot of the intercept parameter `n`

``` r
posterior_samples %>%
  filter(warmup) %>%
  ggplot(aes(x=id, y = n)) + geom_point()
```

## 5.2 traceplots warmup discarded

``` r
posterior_samples %>%
  filter(!warmup) %>%
  ggplot(aes(x=id, y = m)) + geom_point()
```

``` r
posterior_samples %>%
  filter(!warmup) %>%
  ggplot(aes(x=id, y = m)) + geom_line()
```

``` r
posterior_samples %>%
  filter(!warmup) %>%
  ggplot(aes(x=id, y = n)) + geom_point()
```

``` r
posterior_samples %>%
  filter(!warmup) %>%
  ggplot(aes(x=id, y = n)) + geom_line()
```
