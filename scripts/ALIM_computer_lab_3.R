#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Advanced Life insurance mathematics: #
#            PC session 3              #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### 1. Multiple-state models: examples chapter 8 handbook Dickson ####

#### 1.1 Permanent disability model: example 8.4 ####
a1 = 4e-4
a2 = 5e-4
b1 = 3.4674e-6
b2 = 7.5858e-5
c1 = 0.138155
c2 = 0.087498


#### 1.1.1 Using the built-in R function 'integrate' (adaptive quadrature) ####
p00 = function(t, age, a1, b1, c1, a2, b2, c2) {
  p1 = (a1 + a2) * t
  p2 = (b1 / c1) * exp(age * c1) * (exp(c1 * t) - 1)
  p3 = (b2 / c2) * exp(age * c2) * (exp(c2 * t) - 1)
  prob = exp(-(p1 + p2 + p3))
  return(prob)
}

p00(10, 60, a1, b1, c1, a2, b2, c2)

p11 = function(t, age, a1, b1, c1, a2, b2, c2) {
  p1 = a2 * t
  p2 = (b2 / c2) * exp(age * c2) * (exp(c2 * t) - 1)
  prob = exp(-(p1 + p2))
  return(prob)
}

mu01 = function(age, a1, b1, c1, a2, b2, c2) {
  return(a1 + b1 * exp(c1 * age))
}

# Defining the integrand
p01integrand = function(t, age, end, a1, b1, c1, a2, b2, c2) {
  p1 = p00(t, age, a1, b1, c1, a2, b2, c2)
  p2 = mu01(age + t, a1, b1, c1, a2, b2, c2)
  p3 = p11(end - t, age + t, a1, b1, c1, a2, b2, c2)
  prob = p1 * p2 * p3
  return(prob)
}

# Specifying this integrand in the integrate function
integrate(
  p01integrand,
  lower = 0,
  upper = 10,
  age = 60,
  end = 10,
  a1,
  b1,
  c1,
  a2,
  b2,
  c2
)

# Putting everything together in 1 function
f <- function(x, a1, a2, b1, b2, c1, c2) {
  p = exp(-((a1 + a2) * x + b1 / c1 * exp(60 * c1) * (exp(c1 * x) - 1) + 
              b2 / c2 * exp(60 * c2) * (exp(c2 * x) - 1))) * 
    exp(-(a2 * (10 - x) + b2 / c2 * exp((60 + x) * c2) * (exp(c2 * (10 - x)) - 1)))
  p * (a1 + b1 * exp(c1 * (60 + x)))
}

integrate(
  f,
  0,
  10,
  a1 = a1,
  a2 = a2,
  b1 = b1,
  b2 = b2,
  c1 = c1,
  c2 = c2
)

#### 1.1.2 Trapezoidal rule ####
h     = 1 / 12
a     = 0
b     = 10
age   = 60
grid  = seq(from = a, to = b, by = h)
ngrid = length(grid)

f =
  p00(grid, age, a1, b1, c1, a2, b2, c2) * mu01(age + grid, a1, b1, c1, a2, b2, c2) *
  p11(10 - grid, age + grid, a1, b1, c1, a2, b2, c2)

(inttrap = h * (0.5 * f[1] + sum(f[2:(ngrid - 1)]) + 0.5 * f[ngrid]))

#### 1.1.3 Simpson rule ####
n = (b - a) / (2 * h)

(intSimp =
    (h / 3) * (f[1] + 4 * sum(f[2 * 1:n]) + 2 * sum(f[2 * 1:(n - 1) + 1]) +
                 f[ngrid]))


#### 1.2 Disability model (allowing for recovery): example 8.5 ####
mu01 <- function(a1, b1, c1, x) a1 + b1 * exp(c1 * x)
mu10 <- function(a1, b1, c1, x) 0.1 * mu01(a1, b1, c1, x)
mu02 <- function(a2, b2, c2, x) a2 + b2 * exp(c2 * x)
mu12 <- mu02

## Kolmogorov equations 
thP00 <- function(tP00, tP01, h, x) 
  tP00 - h * tP00 * (mu01(a1, b1, c1, x) + mu02(a2, b2, c2, x)) + h * tP01 * mu10(a1, b1, c1, x)
thP01 <- function(tP00, tP01, h, x)
  tP01 - h * tP01 * (mu12(a2, b2, c2, x) + mu10(a1, b1, c1, x)) + h * tP00 * mu01(a1, b1, c1, x)

formals(thP00)          # h defies step size, x defines age
thP00(1, 0, 1 / 12, 60) # Table 8.1, second row, column tP^{00}_{60}

## Reproduction table 8.1
Results = matrix(NA, 12 * 10 + 1, 2)
Results[1, ] = c(thP00(1, 0, 0, 60), thP01(1, 0, 0, 60))
for(i in 1:(12 * 10)) {
  Results[i + 1, ] = c(thP00(Results[i, 1], Results[i, 2], 1 / 12, 60 + (i - 1) * 1/ 12),
                       thP01(Results[i, 1], Results[i, 2], 1 / 12, 60 + (i - 1) * 1 / 12))
}
head(Results)
tail(Results)

t = seq(0, 10, by = 1 / 12)
Table8.1 = cbind.data.frame(t = t,
                            mu01 = mu01(a1, b1, c1, 60 + t), 
                            mu02 = mu02(a2, b2, c2, 60 + t),
                            mu10 = mu10(a1, b1, c1, 60 + t),
                            mu12 = mu12(a2, b2, c2, 60 + t),
                            thP00 = Results[, 1],
                            thP01 = Results[, 2])
head(round(Table8.1, 5))
tail(round(Table8.1, 5))


#### 1.3 Premiums: example 8.6 ####

## Handbook Dickson: All values below have been calculated using the repeated Simpson's rule,  ##
##                   with h = 1 / 12, using table 8.1                                          ##
## --> Hence, using the function integrate here, you will arrive at slightly different results ##

#### 1.3.1 example 8.6 (a) #### 
i     = 0.05
delta = log(1 + i)
h     = 1/12
a     = 0
b     = 10
n     = (b - a) / (2 * h)
grid  = seq(0, 10, by = h)
ngrid = length(grid)

list2env(Table8.1[, c("mu02", "mu12", "thP00", "thP01")], envir = .GlobalEnv)

# discounting
dis      = numeric(ngrid)
dis      = exp(-delta * grid)

## EPV of premium income ##
thP00dis = dis * thP00

# using trapezium rule
(intprem =
    h * (0.5 * thP00dis[1] + sum(thP00dis[2:(ngrid - 1)]) + 0.5 * thP00dis[ngrid]))

# using Simpson's rule
(intprem =
    (h / 3) * (thP00dis[1] + 4 * sum(thP00dis[2 * 1:n]) + 2 * sum(thP00dis[2 * 1:(n - 1) + 1]) + thP00dis[ngrid]))


## EPV of sickness benefit ##
thP01dis = dis * thP01

# using trapezium rule
(intsick = h * (0.5 * thP01dis[1] + sum(thP01dis[2:(ngrid - 1)]) + 0.5 *
                  thP01dis[ngrid]))

# using Simpson's rule
(intsick = (h / 3) * (thP01dis[1] + 4 * sum(thP01dis[2 * 1:n]) + 2 * sum(thP01dis[2 * 1:(n - 1) + 1]) + thP01dis[ngrid]))

## EPV of death benefit ##
fdis = dis * (thP00 * mu02 + thP01 * mu12)

# using trapezium rule
(intdeath = h * (0.5 * fdis[1] + sum(fdis[2:(ngrid - 1)]) + 0.5 * fdis[ngrid]))

# using Simpson's rule
(intdeath = (h / 3) * (fdis[1] + 4 * sum(fdis[2 * 1:n]) + 2 * sum(fdis[2 * 1:(n - 1) + 1]) + fdis[ngrid]))

## Premium
(P = (20000 * intsick + 50000 * intdeath) / intprem)

#### 1.3.2 example 8.6 (b) ####

## EPV of premium income ##
(prem = 1 / 12 * sum(thP00dis[-ngrid]))

## EPV of sickness benefit ##
(sick = 1 / 12 * sum(thP01dis[-1]))

## EPV of death benefit ##
intdeath

## Premium ##
(P = (20000 * sick + 50000 * intdeath) / prem)
P / 12

#### 1.4 Policy values and Thiele's differential equation: example 8.7 ####
x = 40
n = 20
d = 0.04
B = 100000
S = 500000
P = 5500 # Change to 6000 to get solution for example 8.7 (c) (i) (2)

tV0 = 0
tV1 = 0
h   = 1/12
t   = 10

a1 = 4e-4
a2 = 5e-4
b1 = 3.4674e-6
b2 = 7.5858e-5
c1 = 0.138155
c2 = 0.087498
mu01 <- function(x) a1 + b1 * exp(c1 * (x))
mu10 <- function(x) 0.1 * mu01(x)
mu02 <- function(x) a2 + b2 * exp(c2 * (x))
mu12 <- mu02

for(i in seq_len(120 * 2)) {
  tV0[i + 1] = tV0[i] * (1 - d * h) - P * h + h * mu01(x + n - ((i - 1) * h)) * (tV1[i] - tV0[i]) + 
    h * mu02(x + n - ((i - 1) * h)) * (S - tV0[i])
  tV1[i + 1] = tV1[i] * (1 - d * h) + B * h + h * mu10(x + n - ((i - 1) * h)) * (tV0[i] - tV1[i]) + 
    h * mu12(x + n - ((i - 1) * h)) * (S - tV1[i])
}

## 8.7 (c) (i) (1) ##
tV0[120 + 1]
tV1[120 + 1]
tV0[120 * 2 + 1]

## Finding the equivalence premium by defining a function and using the function uniroot ##
f <- function(P, decades = 2) {
  x = 40
  n = 20
  d = 0.04
  B = 100000
  S = 500000
  
  tV0 = 0
  tV1 = 0
  h   = 1/12
  t   = 10
  
  a1 = 4e-4
  a2 = 5e-4
  b1 = 3.4674e-6
  b2 = 7.5858e-5
  c1 = 0.138155
  c2 = 0.087498
  mu01 <- function(x) a1 + b1 * exp(c1 * (x))
  mu10 <- function(x) 0.1 * mu01(x)
  mu02 <- function(x) a2 + b2 * exp(c2 * (x))
  mu12 <- mu02
  
  for(i in seq_len(120 * decades)) {
    tV0[i + 1] = tV0[i] * (1 - d * h) - P * h + h * mu01(x + n - ((i - 1) * h)) * (tV1[i] - tV0[i]) + 
      h * mu02(x + n - ((i - 1) * h)) * (S - tV0[i])
    tV1[i + 1] = tV1[i] * (1 - d * h) + B * h + h * mu10(x + n - ((i - 1) * h)) * (tV0[i] - tV1[i]) + 
      h * mu12(x + n - ((i - 1) * h)) * (S - tV1[i])
  }
  return(tV0[120 * decades + 1])
}

(P = uniroot(f, c(5e3, 6e3))$root)

## Note that function goes from t = n to t = decades.      ##
## Hence, to arrive at {0}_V^{0}, we use f(P, decades = 2) ##
## as n = 20.                                              ##

f(5500, decades = 1) # Policy value with t = 10 and P = 5500
f(5500, decades = 2) # Policy value with t = 0 and P = 5500

f(6e3, decades = 1)  # Policy value with t = 10, with P = 6000
f(6e3, decades = 2)  # Policy value with t = 0 and P = 600

f(P, decades = 2)    # Policy value at begin of the contract, with P = equivalence premium

## Using techniques of example 8.6 ##
a00 = 12.8535
a01 = 0.31593
A02 = 0.08521

(S * A02 + B * a01) / a00


#### 1.5 Multiple decrement models: example 8.8 ####

## (a) ##
tP00 <- function(t) exp(-0.002 * t - 0.00025 * t^2)
mu01 <- function(t) 0.002 + 0.0005 * (t - 50)
f1   <- function(t, delta = 0.025) exp(-delta * t) * tP00(t) * mu01(50 + t)
f2   <- function(t, delta = 0.025) exp(-delta * t) * tP00(t)

(Abar = integrate(f1, 0, 10)$value)
(abar = integrate(f2, 0, 10)$value)

2e5 * Abar / abar

## (b) ##
tP00 <- function(t) exp(-0.052 * t - 0.00025 * t^2)
(Abar = integrate(f1, 0, 10)$value)
(abar = integrate(f2, 0, 10)$value)

2e5 * Abar / abar


