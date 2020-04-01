#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Advanced Life insurance mathematics: #
#            PC session 2              #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### 0. Settings ####
# install.packages("tidyverse")
# install.packages("demography")
# install.packages("forecast")
library(tidyverse)
library(demography)
library(forecast)
# install.packages("rstudioapi")
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

# colors - I copied most of these from # https://github.com/edrubin/EC524W20
KULbg <- "#116E8A"

## ------------------------------------------------------------------------------------
library(readr)
Belgium_male_1x1 <- readr::read_delim(
                    file = "../data/Belgium_male_life_table_1x1.txt", 
                    delim = " ", col_types = "nnnnnnnnnn", 
                    col_names = FALSE, skip = 1)

names(Belgium_male_1x1) <- c("Year", "Age", "mx", "qx", "ax", "lx", "dx", "Lx", "Tx", "ex")

## ------------------------------------------------------------------------------------
names(Belgium_male_1x1)

## ------------------------------------------------------------------------------------
min(Belgium_male_1x1$Year)
max(Belgium_male_1x1$Year)

## ------------------------------------------------------------------------------------
min(Belgium_male_1x1$Age)
max(Belgium_male_1x1$Age)

## ------------------------------------------------------------------------------------
library(dplyr)
start_year <- 1950 
end_year <- max(Belgium_male_1x1$Year)

Belgium_male <- filter(Belgium_male_1x1, Year >= start_year)

Belgium_male %>% select(Year, Age, mx, qx, lx, dx, ex) 

#### 1. Least-squares approach: iterative approach ####

## Alpha is estimated immediately, beta and kappa iteratively ##

#### 1.1. Estimate alpha ####
X <- model.matrix(
          ~ as.factor(Belgium_male$Age) 
          - 1)
dim(X)
y <- log(Belgium_male$mx)

alpha_est_expl <- solve(
          crossprod(X)) %*% 
                   t(X) %*% y # DIY
alpha_est <- lm(
              y ~ -1 + 
              as.factor(Belgium_male$Age))$coef 
              # use lm(.) function
              # extract $coef from the fitted 
              # lm(.) object

## ------------------------------------------------------------------------------------
library(ggplot2)
data <- tibble(age = 0:110, fit = alpha_est)
g <- ggplot(data) + geom_point(aes(age, alpha_est)) + 
        geom_line(aes(age, alpha_est), col = "black") +
        theme_bw() +
        ggtitle("Belgium - males, 1950 - 2018") + 
      labs(y = bquote(hat(beta)[x]^"(1)")) 

g

#### 1.2. Estimate kappa's: initialization ####

## Create new response variable: Z = log(mx) - alpha ##
Z <- log(Belgium_male$mx) - 
     alpha_est[Belgium_male$Age - 
                 min(Belgium_male$Age) + 1]

X <- model.matrix(~ as.factor(Belgium_male$Year) - 1)

kappa_est_expl <- solve(crossprod(X)) %*% t(X) %*% Z

kappa_est <- lm(
              Z ~ -1 + 
              as.factor(Belgium_male$Year))$coef 

## ------------------------------------------------------------------------------------
library(ggplot2)
data <- tibble(year = start_year:end_year, fit = kappa_est)
g <- ggplot(data) + geom_point(aes(year, kappa_est)) + 
        geom_line(aes(year, kappa_est), col = "black") +
        theme_bw() +
        ggtitle("Belgium - males, 1950 - 2018, starting values") + 
      labs(y = bquote(hat(kappa)[t]^"(2)")) 

g

#### 1.3. Estimate betas: initialization ####
## create variable with kappa's as values ##
var_kappa <- kappa_est[Belgium_male$Year - 
                         min(Belgium_male$Year) + 1]

## : in model formula -> interaction term ##
X <- model.matrix(~ 
              as.factor(Belgium_male$Age):var_kappa 
              - 1)

beta_est_expl  <- solve(crossprod(X)) %*% t(X) %*% Z 

beta_est <- lm(
          Z ~ -1 + 
          as.factor(Belgium_male$Age):var_kappa)$coef 

## ------------------------------------------------------------------------------------
library(ggplot2)
data <- tibble(age = 0:110, fit = beta_est)
g <- ggplot(data) + geom_point(aes(age, beta_est)) + 
        geom_line(aes(age, beta_est), col = "black") +
        theme_bw() +
        ggtitle("Belgium - males, 1950 - 2018, starting values") + 
      labs(y = bquote(hat(beta)[x]^"(2)")) 

g

#### 1.4. Incorporate this in algorithm ####

converged = F
iter      = 1

while(!converged){  
  beta_est_old  = beta_est
  kappa_est_old = kappa_est
  
  # (2): estimate kappa's
  var_beta = beta_est[Belgium_male$Age - min(Belgium_male$Age) + 1]
  X        = model.matrix(~ as.factor(Belgium_male$Year):var_beta - 1)
  kappa_est = solve(crossprod(X)) %*% t(X) %*% Z

  # (3): estimate beta's
  var_kappa = kappa_est[Belgium_male$Year - min(Belgium_male$Year) + 1]
  X         = model.matrix(~ as.factor(Belgium_male$Age):var_kappa - 1)
  beta_est   = solve(crossprod(X)) %*% t(X) %*% Z 

  # stopping criterion
  converged = 
  max(abs(beta_est - beta_est_old) / abs(beta_est_old), abs(kappa_est - kappa_est_old) / abs(kappa_est_old)) < 1e-8
  iter = iter + 1
  if(iter %% 1e2 == 0)
    cat("\n\nIteration number", iter, "\n\n")
}


## ------------------------------------------------------------------------------------
library(ggplot2)
data <- tibble(age = 0:110, fit = beta_est, init = beta_est_expl)
g <- ggplot(data) + geom_point(aes(age, fit)) + 
        geom_line(aes(age, fit), col = "black") +
        geom_line(aes(age, init), col = "red") +
        theme_bw() +
        ggtitle("Belgium - males, 1950 - 2018, starting + final values") + 
        labs(y = bquote(hat(beta)[x]^"(2)")) 

g

## ------------------------------------------------------------------------------------
## apply constraints ##
beta_est_LS  = beta_est / sum(beta_est)
kappa_est_LS = (kappa_est - mean(kappa_est)) * sum(beta_est)
alpha_est_LS = alpha_est + beta_est * mean(kappa_est)


## ------------------------------------------------------------------------------------
sum(beta_est_LS)
sum(kappa_est_LS)

## ------------------------------------------------------------------------------------
library(ggplot2)

data_period <- tibble(year = start_year:end_year, fit = kappa_est, init = kappa_est_expl)
data_age <- tibble(age = 0:110, fit_alpha = alpha_est, fit_beta = beta_est)

g_1 <- ggplot(data_age) + geom_point(aes(age, fit_alpha)) + 
       geom_line(aes(age, fit_alpha), col = "black") +
       theme_bw() +
       ggtitle("Belgium - males, 1950 - 2018, final estimates") + 
       labs(y = bquote(hat(beta)[x]^"(1)")) 

g_2 <- ggplot(data_age) + geom_point(aes(age, fit_beta)) + 
       geom_line(aes(age, fit_beta), col = "black") +
       theme_bw() + ggtitle("") +
       labs(y = bquote(hat(beta)[x]^"(2)")) 

g_3 <- ggplot(data_period) + geom_point(aes(year, fit)) + 
       geom_line(aes(year, fit), col = "black") +
       theme_bw() + ggtitle("") +
       labs(y = bquote(hat(kappa)[t]^"(2)")) 

library(gridExtra)
grid.arrange(g_1, g_2, g_3, ncol = 2)

#### 2. Lee-carter: GLM formulation ####

#### 2.0 Data ####

library(demography)
country <- c("BEL", "Belgium")
user    <- "vuulenbak42@hotmail.com"
pw      <- "testEAA"
df      <- hmd.mx(country[1], user , pw , country[2])
years   <- 1970:max(df$year)
ages    <- 0:89

#### Creates subset of demogdata object
df  <- demography::extract.years(
                      df, years = years)
df  <- demography::extract.ages(
                      df, ages = ages, 
                      combine.upper = FALSE)

min(df$year)
max(df$year)
max(df$age)

## ------------------------------------------------------------------------------------
etx <- t(df$pop$male)
dim(etx)

dtx <- etx * t(df$rate$male)
dim(dtx)

## ------------------------------------------------------------------------------------
df <- expand.grid(Year = years, Age = ages) # age 0 across all years, age 1 ...
rates <- dtx/etx # years in rows, ages in columns
df$rates <- log(as.vector(rates)) # age 0 across all years, age 1 ...
names(df) <- c("Year", "Age", "logRate")

p <- ggplot(df, aes(x = Age, y = logRate, group = Year)) + 
     geom_line(aes(colour = Year), size = 1, linetype = 1) +
     scale_colour_gradientn(colours = rainbow(10)) +
     scale_x_continuous(breaks = seq(ages[1], tail(ages, 1) + 1, 10)) +
     theme_bw() + ylab(expression("log" ~ m[x])) + xlab("Age (x)") 
p

#### 2. Poisson likelihood and NR scheme ####

#### 2.1. Poisson GLM: using fit701 function ####

source('fitModels.R')

LCfit701 <- fit701(ages, years, etx, dtx, matrix(1, length(years), length(ages)))

## ------------------------------------------------------------------------------------
names(LCfit701)

## ------------------------------------------------------------------------------------
library(ggplot2)

data_period <- tibble(year = years, fit = LCfit701$kappa2)
data_age <- tibble(age = ages, fit_alpha = LCfit701$beta1, fit_beta = LCfit701$beta2)

g_1 <- ggplot(data_age) + geom_point(aes(age, fit_alpha)) + 
       geom_line(aes(age, fit_alpha), col = "black") +
       theme_bw() +
       ggtitle("Belgium - males, 1970 - 2018, Lee Carter, Poisson") + 
       labs(y = bquote(hat(beta)[x]^"(1)")) 

g_2 <- ggplot(data_age) + geom_point(aes(age, fit_beta)) + 
       geom_line(aes(age, fit_beta), col = "black") +
       theme_bw() + ggtitle("") +
       labs(y = bquote(hat(beta)[x]^"(2)")) 

g_3 <- ggplot(data_period) + geom_point(aes(year, fit)) + 
       geom_line(aes(year, fit), col = "black") +
       theme_bw() + ggtitle("") + 
       labs(y = bquote(hat(kappa)[t]^"(2)")) 

library(gridExtra)
grid.arrange(g_1, g_2, g_3, ncol = 2)

#### 2.2. Using self-written function ####
## Function ##
LCNRopt <- function(dxt, ext, eps = 1e-4, maxiter = 1e4) {
  mxt  = dxt / ext
  m    = ncol(ext)
  LL <- function(dxt, ext, Beta, Kappa) {
    Mat = matrix(NA, nrow(dxt), ncol(dxt))
    for(i in seq_len(nrow(dxt)))
      for(j in seq_len(ncol(dxt)))
        Mat[i, j] = dxt[i, j] * (Beta[i, 1] * Beta[i, 2] * Kappa[j]) - 
          ext[i, j] * exp(Beta[i, 1] * Beta[i, 2] * Kappa[j])
      sum(apply(Mat, 1, sum))
  }
  Beta  = cbind(apply(mxt, 1, function(x) sum(log(x))) / ncol(mxt), rep(1 / nrow(dxt), nrow(mxt)))
  Kappa = (m : 1) - (m + 1) / 2
  Conv  = F
  iter  = 0
  LogL  = NULL
  LogL[iter + 1] = LL(dxt, ext, Beta, Kappa)
  
  while(!Conv) {
    if((iter %% 1000) == 0)
      cat("\n\nIteration number", iter, "\n\n")
    for(i in seq_len(nrow(dxt))) {
      B0i = Beta[i, 1]
      B2i = Beta[i, 2]
      dxti = dxt[i, ]
      exti = ext[i, ]
      B1i = B0i - (sum(dxti - exti * exp(B0i + B2i * Kappa))) /
        - (sum(exti * exp(B0i + B2i * Kappa)))
      
      while(abs(B0i - B1i) > 0.01) {
        B0i = B1i
        B1i = B0i - (sum(dxti - exti * exp(B0i + B2i * Kappa))) /
          - (sum(exti * exp(B0i + B2i * Kappa)))
      }
      Beta[i, 1] = B1i
    }
    
    for(i in seq_len(ncol(dxt))) {
      B1 = Beta[, 1]
      B2 = Beta[, 2]
      dxti = dxt[, i]
      exti = ext[, i]
      K0i = Kappa[i] 
      K1i = K0i - (sum( (dxti -  exti * exp(B1 + B2 * K0i)) * B2 ))  / 
        - (sum(exti * exp(B1 + B2 * K0i) * B2^2))
      
      while(abs(K1i - K0i) > 0.01) {
        K0i = K1i
        K1i = K0i - (sum( (dxti -  exti * exp(B1 + B2 * K0i)) * B2 ))  / 
          - (sum(exti * exp(B1 + B2 * K0i) * B2^2))
      }
      Kappa[i] = K1i
    }
    
    SumB2     = sum(Beta[, 2])
    AvgKappa  = mean(Kappa)
    Kappa     = SumB2 * (Kappa - AvgKappa)
    Beta[, 2] = Beta[, 2] / SumB2
    Beta[, 1] = Beta[, 1] + Beta[, 2] * SumB2 * AvgKappa
    
    for(i in seq_len(nrow(dxt))) {
      B0i = Beta[i, 1]
      B2i = Beta[i, 2]
      dxti = dxt[i, ]
      exti = ext[i, ]
      
      B2.1i = B2i - (sum( (dxti - exti * exp(B0i + B2i * Kappa)) * Kappa )) /
        - (sum( exti * (exp(B0i + B2i * Kappa)) * Kappa^2) )
      
      while(abs(B2i - B2.1i) > 0.01) {
        B2i = B2.1i
        B2.1i = B2i - (sum( (dxti - exti * exp(B0i + B2i * Kappa)) * Kappa )) /
          - (sum( exti * (exp(B0i + B2i * Kappa)) * Kappa^2) )
      }
      
      Beta[i, 2] = B2.1i
    }
    
    iter = iter + 1 
    LogL[iter + 1] = LL(dxt, ext, Beta, Kappa)
    if(abs(LogL[iter + 1] - LogL[iter]) < eps)
      break
    if(iter > maxiter)
      break
  }
  if(iter > maxiter)
    warning("Maximum number of iterations exceeded.")
  return(list(Beta = Beta, Kappa = Kappa, LogL = LogL))
}


formals(LCNRopt)

LeeCarterNR = LCNRopt(dxt   = t(dtx),
                      ext   = t(etx), eps = 1e-4, maxiter = 2e3)
Results = list(
  x = ages,
  y = years,
  beta1  = LeeCarterNR$Beta[, 1],
  beta2  = LeeCarterNR$Beta[, 2],
  kappa2 = LeeCarterNR$Kappa 
)

library(ggplot2)

data_period <- tibble(year = years, fit = Results$kappa2)
data_age <- tibble(age = ages, fit_alpha = Results$beta1, fit_beta = Results$beta2)

g_1 <- ggplot(data_age) + geom_point(aes(age, fit_alpha)) + 
       geom_line(aes(age, fit_alpha), col = "black") +
       theme_bw() +
       ggtitle("Belgium - males, 1970 - 2018, Lee Carter, Poisson") + 
       labs(y = bquote(hat(beta)[x]^"(1)")) 

g_2 <- ggplot(data_age) + geom_point(aes(age, fit_beta)) + 
       geom_line(aes(age, fit_beta), col = "black") +
       theme_bw() + ggtitle("") +
       labs(y = bquote(hat(beta)[x]^"(2)")) 

g_3 <- ggplot(data_period) + geom_point(aes(year, fit)) + 
       geom_line(aes(year, fit), col = "black") +
       theme_bw() + ggtitle("") + 
       labs(y = bquote(hat(kappa)[t]^"(2)")) 

library(gridExtra)
grid.arrange(g_1, g_2, g_3, ncol = 2)

#### 2.3. Residual plots ####
grid <- expand.grid(period = years, age = ages)
grid$res <- as.vector(LCfit701$epsilon)
names(grid) <- c("Year", "Age", "Residual")


## ------------------------------------------------------------------------------------
head(grid)

## ------------------------------------------------------------------------------------
grid <- expand.grid(period = years, age = ages)
grid$res <- as.vector(LCfit701$epsilon)
names(grid) <- c("Year", "Age", "Residual")

p <- ggplot(grid, aes(x = Year, y = Age)) + geom_tile(aes(fill = Residual)) +
  scale_fill_gradientn(colours =  topo.colors(7)) +
  theme_bw() + theme(legend.position = "bottom",
                     legend.title = element_text(size = 15))
p

#### 2.4. Inspecting goodness of fit: observed mx versus fitted mx ####

age <- 25
rates <- dtx/etx
df <- tibble(Year = years, obs = rates[, age - min(ages) + 1], fit = exp(LCfit701$beta1[age - min(ages) + 1] + LCfit701$beta2[age - min(ages) + 1] * LCfit701$kappa2))

g_25 <- ggplot(df) + geom_point(aes(Year, obs), col = "black") +
     geom_line(aes(Year, fit), col = "black", linetype = "dashed") +
     theme_bw() + geom_text(x = 2010, y = 0.00125, label = "Age 25", size = 10) +
     ggtitle("Belgium - males, 1970 - 2018, Lee Carter, Poisson") + 
     labs(y = bquote(hat(m)[25,][t])) + labs(x = bquote(Year (t)))

age <- 45
rates <- dtx/etx
df <- tibble(Year = years, obs = rates[, age - min(ages) + 1], fit = exp(LCfit701$beta1[age - min(ages) + 1] + LCfit701$beta2[age - min(ages) + 1] * LCfit701$kappa2))

g_45 <- ggplot(df) + geom_point(aes(Year, obs), col = "black") +
     geom_line(aes(Year, fit), col = "black", linetype = "dashed") +
     theme_bw() + geom_text(x = 2010, y = 0.004, label = "Age 45", size = 10) +
     ggtitle("Belgium - males, 1970 - 2018, Lee Carter, Poisson") + 
     labs(y = bquote(hat(m)[45,][t])) + labs(x = bquote(Year (t)))

age <- 65
rates <- dtx/etx
df <- tibble(Year = years, obs = rates[, age - min(ages) + 1], fit = exp(LCfit701$beta1[age - min(ages) + 1] + LCfit701$beta2[age - min(ages) + 1] * LCfit701$kappa2))

g_65 <- ggplot(df) + geom_point(aes(Year, obs), col = "black") +
     geom_line(aes(Year, fit), col = "black", linetype = "dashed") +
     theme_bw() + geom_text(x = 2010, y = 0.03, label = "Age 65", size = 10) +
     ggtitle("Belgium - males, 1970 - 2018, Lee Carter, Poisson") + 
     labs(y = bquote(hat(m)[65,][t])) + labs(x = bquote(Year (t)))

age <- 85
rates <- dtx/etx
df <- tibble(Year = years, obs = rates[, age - min(ages) + 1], fit = exp(LCfit701$beta1[age - min(ages) + 1] + LCfit701$beta2[age - min(ages) + 1] * LCfit701$kappa2))

g_85 <- ggplot(df) + geom_point(aes(Year, obs), col = "black") +
     geom_line(aes(Year, fit), col = "black", linetype = "dashed") +
     theme_bw() + geom_text(x = 2010, y = 0.18, label = "Age 85", size = 10) +
     ggtitle("Belgium - males, 1970 - 2018, Lee Carter, Poisson") + 
     labs(y = bquote(hat(m)[85,][t])) + labs(x = bquote(Year (t)))

grid.arrange(g_25, g_45, g_65, g_85, ncol = 2)

#### 2.5. Forecasting with time series ####

library(forecast)

time_series <- Arima(LCfit701$kappa2, 
                     order = c(0, 1, 0), 
                     include.drift = TRUE)

time_series


## ------------------------------------------------------------------------------------
plot(forecast(time_series, level = c(80, 85, 95)))

#### 2.6. Forecasting with the model: longevity charts ####

source('simModels.R')

sim_LC = sim2001(
  xx = LCfit701$x,
  yy = LCfit701$y,
  beta1v = LCfit701$beta1,
  beta2v = LCfit701$beta2,
  kappa2v = LCfit701$kappa2,
  nsim = 10000,
  tmax = 50,
  nyears = length(years)
)

names(sim_LC)


## ------------------------------------------------------------------------------------
sim_LC$y

dim(sim_LC$dda)

sim_LC$dda[ , 1:50, 1]

LCfit701$kappa2[length(years)]

dim(sim_LC$qaa)


## ------------------------------------------------------------------------------------
# time series kappa_t
plot(
  LCfit701$y,
  LCfit701$kappa2,
  pch = 20,
  xlim = c(1970, 2068),
  ylim = range(c(
    range(LCfit701$kappa2), range(sim_LC$dda[, , ])
  )),
  main = bquote(paste("Projection of ", kappa[t])),
  sub = "ARIMA(0,1,0), 10000 sim",
  xlab = "Year (t)",
  ylab = bquote(kappa[t])
)

# fan chart
fan(sim_LC$y, sim_LC$dda[, , ], color = "red")

## ------------------------------------------------------------------------------------
age = 65
minage = min(ages)

plot(
  LCfit701$y,
  exp(LCfit701$beta1[age - minage + 1] + LCfit701$beta2[age - minage + 1] * LCfit701$kappa2),
  lwd = 2,
  col = "red",
  type = "l",
  ylim = c(0, 0.04),
  main = bquote(paste(
    "Belgium: Lee-Carter, projection of ", mu[65](t)
  )),
  xlim = c(1970, 2068),
  sub = "ARIMA(0,1,0), 10000 sim",
  xlab = "Year (t)",
  ylab = bquote(paste(mu[65](t)))
)
fan(sim_LC$y,
    exp(LCfit701$beta1[age - minage + 1] + LCfit701$beta2[age - minage + 1] * sim_LC$dda[, , ]),
    col = "red")
points(LCfit701$y, rates[, age - minage + 1], col = "black", pch = 20)


## ------------------------------------------------------------------------------------
# plot of q_{65,75,85}
ages_sel  = c(65, 75, 85)
nages = length(ages_sel)
color = c("red", "green", "blue")

# y axis is logarithmic
plot(
  c(1970, 2068),
  c(0.005, 0.2),
  type = "n",
  log = "y",
  lwd = 2,
  col = "black",
  main = bquote(paste(
    "Belgium: projection of ", q[65](t), ", ", q[75](t), ", ", q[85](t)
  )),
  sub = "ARIMA(0,1,0), 10000 sim",
  xlab = "Year (t)",
  ylab = bquote(paste(q[65](t), ", ", q[75](t), ", ", q[85](t)))
)

for (j in 1:nages) {
  lines(LCfit701$y, 1 - exp(-exp(
    LCfit701$beta1[ages_sel[j] - minage + 1] + LCfit701$beta2[ages_sel[j] - minage + 1] *
      LCfit701$kappa2
  )), col = color[j], lwd = 2)
  matlines(sim_LC$y, sim_LC$qaa[ages_sel[j] - minage + 1, , 1:1000], col =
             "grey")
}
p = seq(0.05, 0.95, 0.05)

# Construct the corresponding empirical quantiles
for(j in 1:nages) {
  q_int = matrix(data = 0,
                 nrow = 50,
                 ncol = length(p))
  for (i in 1:50)
    q_int[i, ] = quantile(sim_LC$qaa[ages_sel[j] - minage + 1, i, ], p = p)
  # highlight 0.05 and 0.95 quantile
  lines(2018:2067,
        q_int[1:50, 1],
        col = color[j],
        lty = 1,
        lwd = 2)
  lines(2018:2067,
        q_int[1:50, length(p)],
        col = color[j],
        lty = 1,
        lwd = 2)
}

## ------------------------------------------------------------------------------------
# Alternative: plot of q_{65,75,85} using fan function in simModels.R
plot(
  c(1970, 2058),
  c(0.005, 0.2),
  log = "y",
  main = bquote(paste(
    "Belgium: projection of ", q[65](t), ", ", q[75](t), ", ", q[85](t)
  )),
  sub = "ARIMA(0,1,0), 10000 sim",
  xlab = "Year (t)",
  ylab = bquote(paste(q[65](t), ", ", q[75](t), ", ", q[85](t))),
  type = "n"
)


for (j in 1:nages) {
  points(LCfit701$y, 1 - exp(-exp(
    LCfit701$beta1[ages_sel[j] - minage + 1] + LCfit701$beta2[ages_sel[j] - minage + 1] *
      LCfit701$kappa2
  )), col = color[j], pch = 20)
  fan(sim_LC$y, sim_LC$qaa[ages_sel[j] - minage + 1, ,], color = color[j])
}


