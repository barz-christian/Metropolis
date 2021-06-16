A Gibbs Sampler
================
Christian Barz
created: 2021-06-16, reviewed: 2021-06-16

-   [1 A bayesian regression](#a-bayesian-regression)
    -   [1.1 define conditional posterior
        probabilities](#define-conditional-posterior-probabilities)
    -   [1.2 gibbs sampler](#gibbs-sampler)
    -   [1.3 test](#test)

# 1 A bayesian regression

We consider the following regression

$$
\\begin{align}
y &= 
\\beta\_0 + \\beta\_1 x + \\epsilon\\\\
\\epsilon &\\sim 
N(0,1/\\tau)
\\end{align}
$$

To do a bayesian regression we need to set prior on the parameters. We
do this by using conditionally conjugated priors. A *conditionally
conjugated prior for a parameter *x** for a given likelihood gives a
posteriori when conditioned on all parameters

$$
\\begin{align}
\\beta\_i  & \\sim N(\\mu\_i,1/\\tau\_i)\\\\
\\tau & \\sim Gamma(a,b)
\\end{align}
$$

We assume that the hyperparameters
*μ*<sub>0</sub>, *μ*<sub>1</sub>, *τ*<sub>0</sub>, *τ*<sub>1</sub>, *a*, *b*
are known constants.

By bayes rule the posterior distribution is proportional to

*p*(*β*, *τ*\|*y*) ∝ *p*(*y*\|*β*, *τ*)*p*(*β*)*p*(*τ*)

To do Gibbs sampling, we need the conditional posteriori distributions
of the parameters.

$$
\\beta\_0\|x,y,\\beta\_1,\\tau \\sim N(\\frac{\\tau\_0\\mu\_0+\\tau\\sum\_i(y\_i-\\beta\_1x\_i)}{\\tau\_0+n\\tau},\\frac{1}{\\tau\_0+n\\tau})
$$
$$
\\beta\_1\|x,y,\\beta\_0,\\tau \\sim 
N(\\frac{
  \\tau\_1\\mu\_1 + 
  \\tau\\sum\_i(y\_i-\\beta\_0)x\_i 
}{
  \\tau\_1+\\tau\\sum\_ix\_i^2
  }
,
\\frac{1}{\\tau\_1+\\tau\\sum\_ix\_i^2}
)
$$

$$
\\tau\|\\beta\_0,\\beta\_1,a,b,x,y \\sim
Gamma(a + \\frac{N}{2}, b + \\sum\_i\\frac{(y\_i\\beta\_0-\\beta\_1x\_i)^2}{2})
$$

## 1.1 define conditional posterior probabilities

``` r
sample_beta0 <- function(y,x,beta1,tau,tau0,mu0){
  N <- NROW(y)
  precision <- 1 / (tau0 + N * tau)
  mittelwert <- (tau0 * mu0 + tau * sum(y - beta1 * x) ) / precision
  
  return(rnorm(n = 1, 
               mean = mittelwert,
               sd = 1/ sqrt(precision))
         )
}
```

``` r
sample_beta1 <- function(y,x,beta0,tau,tau1,mu1){
  N <- NROW(y)
  precision <- 1 / (tau1 + tau * sum(x * x))
  mittelwert <- (tau1 * mu1 + tau * sum((y-beta0) * x)) / precision
  
  return(rnorm(n = 1, 
               mean = mittelwert, 
               sd = 1/ sqrt(precision))
         )
}
```

``` r
sample_tau <- function(y,x,beta0, beta1, a,b){
  N <- NROW(y)
  a_posteriori <- a + N/2
  b_posteriori <- b + sum((y - beta0 - beta1 * x)^2)/2
  return(
    rgamma(n = 1, 
           shape = a_posteriori,
           scale = b_posteriori)
  )
}
```

## 1.2 gibbs sampler

``` r
gibbs <- function(y,x, 
                  initialvalue, 
                  #conditional_posteriors,
                  hyperparameter,
                  iterations){
  
  
  samples <- array(data = NA, 
                   dim = c(iterations,
                           length(initialvalue)))
  
  colnames(samples) <- names(initialvalue)
  
  samples[1,] <- initialvalue %>%
    purrr::reduce(cbind)
  
  for(i in 2:iterations){
    
    samples[i,"beta0"] <- sample_beta0(y = y, 
                                       x = x, 
                                       beta1 = samples[i-1,"beta1"],
                                       tau = samples[i-1,"tau"],
                                       tau0 = hyperparameter[["tau0"]],
                                       mu0 =  hyperparameter[["mu0"]])
    
    samples[i,"beta1"] <- sample_beta1(y = y, x = x,
                           beta0 = samples[i,"beta0"], 
                           tau = samples[i-1,"tau"], 
                           tau1 = hyperparameter[["tau1"]],
                           mu1 = hyperparameter[["mu1"]]
                            )
    
    samples[i,"tau"] <- sample_tau(y = y, x = x,
                                   beta0 = samples[i,"beta0"], 
                                   beta1 = samples[i,"beta1"],
                                   a = hyperparameter[["a"]], 
                                   b = hyperparameter[["b"]]
                                   )
  }
  
  return(samples)
}
```

## 1.3 test

### 1.3.1 synthetic data

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
```

### 1.3.2 initial values, hyperparameters

``` r
conditional_posteriors <- list(sample_beta0, 
                               sample_beta1,
                               sample_tau)
names(conditional_posteriors) <- c("beta0", "beta1", "tau")

initialvalue <- list(1,1,1)
names(initialvalue) <- c("beta0", "beta1", "tau")

hyperparameters <- list(0,1,0,1,2,2)
names(hyperparameters) <- c("mu0", "tau0",
                            "mu1", "tau1", 
                            "a", "b")
```

``` r
gibbs(y = y, x = x, 
      initialvalue = initialvalue,
      hyperparameter = hyperparameters, 
      iterations = 10)
```

    ## Warning in rnorm(n = 1, mean = mittelwert, sd = 1/sqrt(precision)): NAs
    ## produziert

    ## Warning in rgamma(n = 1, shape = a_posteriori, scale = b_posteriori): NAs
    ## produziert

    ## Warning in rnorm(n = 1, mean = mittelwert, sd = 1/sqrt(precision)): NAs
    ## produziert

    ## Warning in rnorm(n = 1, mean = mittelwert, sd = 1/sqrt(precision)): NAs
    ## produziert

    ## Warning in rgamma(n = 1, shape = a_posteriori, scale = b_posteriori): NAs
    ## produziert

    ## Warning in rnorm(n = 1, mean = mittelwert, sd = 1/sqrt(precision)): NAs
    ## produziert

    ## Warning in rnorm(n = 1, mean = mittelwert, sd = 1/sqrt(precision)): NAs
    ## produziert

    ## Warning in rgamma(n = 1, shape = a_posteriori, scale = b_posteriori): NAs
    ## produziert

    ## Warning in rnorm(n = 1, mean = mittelwert, sd = 1/sqrt(precision)): NAs
    ## produziert

    ## Warning in rnorm(n = 1, mean = mittelwert, sd = 1/sqrt(precision)): NAs
    ## produziert

    ## Warning in rgamma(n = 1, shape = a_posteriori, scale = b_posteriori): NAs
    ## produziert

    ## Warning in rnorm(n = 1, mean = mittelwert, sd = 1/sqrt(precision)): NAs
    ## produziert

    ## Warning in rnorm(n = 1, mean = mittelwert, sd = 1/sqrt(precision)): NAs
    ## produziert

    ## Warning in rgamma(n = 1, shape = a_posteriori, scale = b_posteriori): NAs
    ## produziert

    ## Warning in rnorm(n = 1, mean = mittelwert, sd = 1/sqrt(precision)): NAs
    ## produziert

    ## Warning in rnorm(n = 1, mean = mittelwert, sd = 1/sqrt(precision)): NAs
    ## produziert

    ## Warning in rgamma(n = 1, shape = a_posteriori, scale = b_posteriori): NAs
    ## produziert

    ## Warning in rnorm(n = 1, mean = mittelwert, sd = 1/sqrt(precision)): NAs
    ## produziert

    ## Warning in rnorm(n = 1, mean = mittelwert, sd = 1/sqrt(precision)): NAs
    ## produziert

    ## Warning in rgamma(n = 1, shape = a_posteriori, scale = b_posteriori): NAs
    ## produziert

    ##              beta0          beta1           tau
    ##  [1,] 1.000000e+00   1.000000e+00  1.000000e+00
    ##  [2,] 7.842933e+04  -7.561141e+11  1.680671e+29
    ##  [3,] 1.506993e+75 -4.104855e+140 4.674853e+286
    ##  [4,]          Inf            NaN           NaN
    ##  [5,]          NaN            NaN           NaN
    ##  [6,]          NaN            NaN           NaN
    ##  [7,]          NaN            NaN           NaN
    ##  [8,]          NaN            NaN           NaN
    ##  [9,]          NaN            NaN           NaN
    ## [10,]          NaN            NaN           NaN
