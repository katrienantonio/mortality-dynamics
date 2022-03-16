#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Advanced Life insurance mathematics: #
#            PC session 1              #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### 0. Settings ####
# install.packages("tidyverse")
# install.packages("demography")
library(tidyverse)
library(demography)
# install.packages("rstudioapi")
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

# colors - I copied most of these from # https://github.com/edrubin/EC524W20
KULbg <- "#116E8A"

#### 0.1 Example to read in data from Human Mortality database ####
? hmd.mx
User = "vuulenbak42@hotmail.com"
pw   = "testEAA"
Df   = hmd.mx("BEL", User , pw , "Belgium")

## -------------------------------------------------------------------------------
names(Df)
names(Df$pop)
names(Df$rate)

## -------------------------------------------------------------------------------
str(Df$rate$female)

## -------------------------------------------------------------------------------
Df$type

## -------------------------------------------------------------------------------
Df$label

## -------------------------------------------------------------------------------
min(Df$year)
max(Df$year)

## -------------------------------------------------------------------------------
min(Df$age)
max(Df$age)

## -------------------------------------------------------------------------------
Belgium_female_2018 = read.table(file = "../data/Belgium_female_2018.txt", header = TRUE)

## -------------------------------------------------------------------------------
head(Belgium_female_2018)

## -------------------------------------------------------------------------------
Belgium_male_2018 = read.table(file = "../data/Belgium_male_2018.txt", header = TRUE)

## -------------------------------------------------------------------------------
KULbg <- "#116E8A"
library(ggplot2)
g_male <- ggplot(Belgium_male_2018, aes(Age, log(qx))) + 
          geom_point(col = KULbg) + 
          geom_line(stat = "identity", col = KULbg) + 
          theme_bw() +
          ggtitle("Belgium - males, 2018") + 
          labs(y = bquote(ln(q[x])))
g_male

## -------------------------------------------------------------------------------
KULbg <- "#116E8A"
library(ggplot2)
g_fem <- ggplot(Belgium_female_2018, aes(Age, log(qx))) + 
            geom_point(col = KULbg) + 
            geom_line(stat = "identity", col = KULbg) + 
            theme_bw() +
            ggtitle("Belgium - females, 2018") + 
            labs(y = bquote(ln(q[x])))
g_fem

#### 1. Parametric models for mortality ####

#### 1.1. Makeham ####

#### 1.1.1 Ballegeer's method ####
f <- function(par, x, px) {
  th5 = par[2]
  th3 = par[1]
  th6 = par[3]
  sum((log(th5) + (th3 - 1) * th3^x * log(th6) - log(px))^2)
}

# Initial values
par.init = c(1.1049542, 0.9999124, 0.9998060)

# resulting parameters
(Ballegeer = optim(par.init, f, x = Belgium_male_2018$Age[41:92],
                   lower = c(1, -Inf, -Inf), upper = c(Inf, 1, 1),
                   method = "L-BFGS-B",
                   px = 1 - Belgium_male_2018$qx[41:92])$par)
theta_3B = Ballegeer[1]
theta_5B = Ballegeer[2]
theta_6B = Ballegeer[3]

# tranform parameters Ballegeer to Makeham parameters original parameterization
th1B = -log(Ballegeer[2])
th2B = -log(Ballegeer[3]) * log(Ballegeer[1])
th3B = Ballegeer[1]

c(th1B, th2B, th3B)

## -------------------------------------------------------------------------------
g_male <- 
  ggplot(Belgium_male_2018, aes(Age, log(qx))) + 
  geom_point(col = KULbg) + 
  geom_line(stat = "identity", col = KULbg) + 
  theme_bw() +
  ggtitle("Belgium - males, 2018") + 
  labs(y = bquote(ln(q[x])))

df_Ball <- data.frame(
                Age = 40:91, 
                Makeham = log(1 - theta_5B * theta_6B ^ ((theta_3B - 1) * theta_3B ^ (40:91)))
            )

g_male <- g_male + geom_line(data = df_Ball, 
                    aes(Age, Makeham), 
                    col = "red", lwd = 1.2) 

g_male

####  1.1.2 De Vylder's method ####
BinLL <- function(par, x, dx, lx) {
  px = exp(-(par[1] + par[3] * par[2]^x)) 
  -sum((lx - dx) * log(px) + dx * log(1 - px))
}

par.init = c(8.758055e-05, 1.104954, 2.036360e-05)
fit1     = optim(par.init, BinLL, x = Belgium_male_2018$Age[41:92], lx = Belgium_male_2018$lx[41:92], 
                 dx = Belgium_male_2018$dx[41:92], 
                 control = list(reltol = 1e-10))

theta_1DV = fit1$par[1] ; theta_3DV = fit1$par[2] ; theta_7DV = fit1$par[3]

theta_1M  = theta_1DV
theta_2M  = theta_7DV * log(theta_3DV) / (theta_3DV - 1)
theta_3M  = theta_3DV

c(theta_1M, theta_2M, theta_3M)

## -------------------------------------------------------------------------------
g_male <- 
  ggplot(Belgium_male_2018, aes(Age, log(qx))) + 
  geom_point(col = KULbg) + 
  geom_line(stat = "identity", col = KULbg) + 
  theme_bw() +
  ggtitle("Belgium - males, 2018") + 
  labs(y = bquote(ln(q[x])))

df_DV <- data.frame(
            Age = 40:91, 
            Makeham = log(1 - exp(-(theta_1DV + theta_7DV * theta_3DV ^ (40:91))))
            )

g_male <- g_male + geom_line(data = df_DV, 
                    aes(Age, Makeham), 
                    col = "red", lwd = 1.2) 

g_male

#### 1.1.4 Simulate data from Makeham's law ####
#### 1.1.4.1  First illustration ####
age  = 50

nsim = 100000
v1   = rexp(nsim, 1)
v2   = rexp(nsim, 1)
t1   = (log(theta_2M * theta_3M ^ age +
              v1 * log(theta_3M)) -
          log(theta_2M)) / log(theta_3M) - age
t2   = v2 / theta_1M
t    = pmin(t1, t2)
age_death    = t + age

## -------------------------------------------------------------------------------
# group the simulated lifetimes in integer classes and plot
d_sim <- data.frame(
            x = age:110,
            sim = as.numeric(table(factor(floor(age_death), 
                                   levels = age:110)))
         )

# and plot
g_sim_Makeham <- ggplot(d_sim) + 
                 geom_step(aes(x, sim)) + 
                 theme_bw() + ggtitle("Belgium - males, 2018") + 
                 labs(y = bquote(d[x]))

g_sim_Makeham

# now we add the theoretical results, using Makeham's law
grid = age:110
tgrid = grid - age
p = exp(-theta_1M * tgrid - theta_2M * theta_3M ^ (age) * (theta_3M ^ tgrid - 1) / log(theta_3M))
l = nsim * p
d = c(-diff(l), 0) # add zero to match length of grid

d_theo <- data.frame(
             x = age:110,
             theo = d
)

g_sim_Makeham <- g_sim_Makeham + 
                 geom_step(data = d_theo, aes(x, theo), 
                           color = "red")

g_sim_Makeham


## -------------------------------------------------------------------------------
head(Belgium_male_2018)

#### 1.1.4.2 Second illustration ####
age  = 0

nsim = 100000
v1   = rexp(nsim, 1)
v2   = rexp(nsim, 1)
t1   = (log(theta_2M * theta_3M ^ age +
              v1 * log(theta_3M)) -
          log(theta_2M)) / log(theta_3M) - age
t2   = v2 / theta_1M
t    = pmin(t1, t2)
age_death    = t + age

# group the simulated lifetimes in integer classes and plot
d_sim <- data.frame(
            x = age:110,
            sim = as.numeric(table(factor(floor(age_death), 
                                   levels = age:110)))
         )

# and plot
g_sim_Makeham <- ggplot(d_sim) + 
                 geom_step(aes(x, sim)) + 
                 theme_bw() + ggtitle("Belgium - males, 2018") + 
                 labs(y = bquote(d[x]))

# now we add the theoretical results, using Makeham's law
grid = age:110
tgrid = grid - age
p = exp(-theta_1M * tgrid - theta_2M * theta_3M ^ (age) * (theta_3M ^ tgrid - 1) / log(theta_3M))
l = nsim * p
d = c(-diff(l), 0) # add zero to match length of grid

d_theo <- data.frame(
             x = age:110,
             theo = d
)

g_sim_Makeham <- g_sim_Makeham + 
                 geom_step(data = d_theo, aes(x, theo), 
                           color = "red")

g_sim_Makeham <- g_sim_Makeham + 
                 geom_step(data = Belgium_male_2018, 
                           aes(Age, dx), 
                           color = "blue")

g_sim_Makeham

#### 1.1.5 Heligman & Pollard's law: parametric model ####

#### 1.1.5.1 Only infant mortality ####
HPLLInfMort <- function(par, dx, lx, x){
  theta_1 = par[1]
  theta_2 = par[2]
  theta_3 = par[3]
  
  qxHP1   = theta_1 ^ ((x + theta_2) ^ (theta_3))
  qxHP    = qxHP1 / (1 + qxHP1)
  pxHP    = 1 - qxHP
  
  -sum((lx - dx) * log(pxHP) + dx * log(qxHP))
}


## -------------------------------------------------------------------------------
fitHP1 = optim(c(1e-6, 1e-6, 1e-6), HPLLInfMort, 
               dx = Belgium_male_2018$dx[1:5], 
               lx = Belgium_male_2018$lx[1:5], x = Belgium_male_2018$Age[1:5])

fitHP2 = optim(c(1e-2, 1e-2, 1e-2), HPLLInfMort, 
               dx = Belgium_male_2018$dx[1:5], 
               lx = Belgium_male_2018$lx[1:5], x = Belgium_male_2018$Age[1:5])

fitHP1$value; fitHP2$value
# Note that this is -log-likelihood, so the lower, the better

## -------------------------------------------------------------------------------
PH.lnqx <- function(par, grid) {
  th1 = par[1]
  th2 = par[2]
  th3 = par[3]
  qx1 = th1^((grid + th2)^th3)
  log(qx1 / (1 + qx1))
}

g_male <- 
  ggplot(Belgium_male_2018, aes(Age, log(qx))) + 
  geom_point(col = KULbg) + 
  geom_line(stat = "identity", col = KULbg) + 
  theme_bw() +
  ggtitle("Belgium - males, 2018") + 
  labs(y = bquote(ln(q[x])))

df_HP <- data.frame(
          x = 0:110,
          HP1 = PH.lnqx(fitHP1$par, 0:110),
          HP2 = PH.lnqx(fitHP2$par, 0:110)
)

g_male <- g_male + geom_line(data = df_HP, aes(x, HP1), col = "red")

g_male <- g_male + geom_line(data = df_HP, aes(x, HP2), col = "purple")

g_male



## -------------------------------------------------------------------------------
theta_1HPInf = fitHP2$par[1]
theta_2HPInf = fitHP2$par[2]
theta_3HPInf = fitHP2$par[3]


## ----echo=FALSE-----------------------------------------------------------------
PH.lnqx <- function(par, grid) {
  th1 = par[1]
  th2 = par[2]
  th3 = par[3]
  qx1 = th1^((grid + th2)^th3)
  log(qx1 / (1 + qx1))
}

g_male <- 
  ggplot(Belgium_male_2018, aes(Age, log(qx))) + 
  geom_point(col = KULbg) + 
  geom_line(stat = "identity", col = KULbg) + 
  theme_bw() +
  ggtitle("Belgium - males, 2018") + 
  labs(y = bquote(ln(q[x])))

df_HP <- data.frame(
          x = 0:110,
          HP1 = PH.lnqx(fitHP1$par, 0:110),
          HP2 = PH.lnqx(fitHP2$par, 0:110)
)

g_male <- g_male + geom_line(data = df_HP, aes(x, HP2), col = "purple")
g_male


#### 1.1.5.2 only accident hump ####
HPLLAccid <- function(par, dx, lx, x) {
  theta_4 = par[1]
  theta_5 = par[2]
  theta_6 = par[3]
  
  qxHP1 = theta_4 * exp(-theta_5 * 
                          (log(x) - log(theta_6)) ^ 2)
  qxHP  = qxHP1 / (1 + qxHP1)
  pxHP  = 1 - qxHP
  
  -sum((lx - dx) * log(pxHP) + dx * log(qxHP))
}

par.init = c(0.0001, 10, 18)
par.init = c(0.01, 5, 18)


## -------------------------------------------------------------------------------
fit1 = optim(par.init,
             HPLLAccid,
             dx = Belgium_male_2018$dx[16:25],
             lx = Belgium_male_2018$lx[16:25],
             x = Belgium_male_2018$Age[16:25])

(theta_4HPAcc = fit1$par[1])
(theta_5HPAcc = fit1$par[2])
(theta_6HPAcc = fit1$par[3])

## -------------------------------------------------------------------------------
df_HPAcc <- data.frame(
      x = 0:110,
      qxHP1 = theta_4HPAcc * exp(-theta_5HPAcc * 
                  (log(0:110) - log(theta_6HPAcc)) ^ 2)
)

df_HPAcc$qxHP <- df_HPAcc$qxHP1 / (1 + df_HPAcc$qxHP1)

g_male <- g_male + geom_line(data = df_HPAcc, 
                             aes(x, log(qxHP)), 
                             col = "green") +
          ylim(-12, 0)

g_male

#### 1.1.5.3 only adult mortality ####
HPLLAdult <- function(par,dx,lx,x){
  theta_7 = par[1]
  theta_8 = par[2]
  
  qxHP1 = theta_7 * theta_8 ^ x
  qxHP  = qxHP1 / (1 + qxHP1)
  pxHP  = 1 - qxHP
  
  -sum((lx - dx) * log(pxHP) + dx * log(qxHP))
}

par.init = c(0.01,1)


## -------------------------------------------------------------------------------
fit1 = optim(par.init,
             HPLLAdult,
             dx = Belgium_male_2018$dx[25:99],
             lx = Belgium_male_2018$lx[25:99],
             x = Belgium_male_2018$Age[25:99])

(theta_7HPAd = fit1$par[1])
(theta_8HPAd = fit1$par[2])


## ----echo = FALSE---------------------------------------------------------------
df_HPAd <- data.frame(
      x = 0:110,
      qxHP1 = theta_7HPAd * theta_8HPAd ^ (0:110)
)

df_HPAd$qxHP <- df_HPAd$qxHP1 / (1 + df_HPAd$qxHP1)

g_male <- g_male + geom_line(data = df_HPAd, 
                             aes(x, log(qxHP)), 
                             col = "blue") +
          ylim(-12, 0)

g_male

#### 1.1.5.4 full likelihood ####
HPLL <- function(par,dx,lx,x){
  theta_1 = par[1]
  theta_2 = par[2]
  theta_3 = par[3]
  theta_4 = par[4]
  theta_5 = par[5]
  theta_6 = par[6]
  theta_7 = par[7]
  theta_8 = par[8]
  
qxHP1 = theta_1 ^ ((x + theta_2) ^ (theta_3)) + 
theta_4 * exp(-theta_5 * (log(x) - log(theta_6)) ^ 2) + 
theta_7 * theta_8 ^ x

qxHP  = qxHP1 / (1 + qxHP1)
pxHP  = 1 - qxHP
  
  -sum((lx - dx) * log(pxHP) + dx * log(qxHP))
}

par.init = c(theta_1HPInf,theta_2HPInf,theta_3HPInf,
             theta_4HPAcc,theta_5HPAcc,theta_6HPAcc,
             theta_7HPAd,theta_8HPAd)



## -------------------------------------------------------------------------------
fit1 = optim(par.init,
             HPLL,
             dx = Belgium_male_2018$dx[1:99],
             lx = Belgium_male_2018$lx[1:99],
             x = Belgium_male_2018$Age[1:99])

theta_1HP = fit1$par[1]
theta_2HP = fit1$par[2]
theta_3HP = fit1$par[3]
theta_4HP = fit1$par[4]
theta_5HP = fit1$par[5]
theta_6HP = fit1$par[6]
theta_7HP = fit1$par[7]
theta_8HP = fit1$par[8]


## ----echo = FALSE---------------------------------------------------------------
df_HP <- data.frame(
  x = 0:110,
  qxHP1 = theta_1HP ^ ((0:110 + theta_2HP) ^ (theta_3HP)) + 
    theta_4HP * exp(-theta_5HP *(log(0:110) - log(theta_6HP)) ^ 2) + 
    theta_7HP * theta_8HP ^ (0:110)
)

df_HP$qxHP <- df_HP$qxHP1 / (1 + df_HP$qxHP1)

g_male <- g_male + geom_line(data = df_HP, 
                             aes(x, log(qxHP)), 
                             col = "red") +
          ylim(-12, 0)

g_male

#### 1.1.6 Heligman & Pollard's law: simulating random lifetimes ####
age = 0
x   = 0:110
p   = 1 / (1 + theta_1HP ^ ((x + theta_2HP) ^ 
                              (theta_3HP)) +
             theta_4HP *
             exp(-theta_5HP * (log(x) - 
                            log(theta_6HP)) ^ 2) + 
             theta_7HP * (theta_8HP ^ x))
n   = 110 - age

pmat = cumprod(p[(age + 1):length(p)])


## -------------------------------------------------------------------------------
nsim     = 100000
age_death = numeric(nsim)
u        = runif(nsim, 0, 1)
t1       = numeric(nsim)
t2       = numeric(nsim)

## -------------------------------------------------------------------------------
for(j in 1:nsim) {
  if (u[j] > pmat[1]) {
    t1[j] = 0
  }
  else{
    i = 1
    flag = FALSE
    while (i <= (n - 1) & !flag) {
      if (pmat[i] >= u[j] & u[j] > pmat[i + 1]) {
        t1[j] = i
        flag = TRUE
      }
      i = i + 1
    }
  }
  t = t1[j]
  if (t == 0) {
    t2[j] = log(u[j]) / (log(p[t + 1]))
  }
  if (t > 0) {
    t2[j] = (log(u[j]) - log(pmat[t])) / log(p[t + age + 1])
  }
  age_death[j] = age + t1[j] + t2[j]
}

## -------------------------------------------------------------------------------
# group the simulated lifetimes in integer classes and plot
d_sim <- data.frame(
            x = age:110,
            sim = as.numeric(table(factor(floor(age_death), 
                                   levels = age:110)))
         )

# and plot
g_sim_HP <- ggplot(d_sim) + 
                 geom_step(aes(x, sim)) + 
                 theme_bw() + ggtitle("Belgium - males, 2018") + 
                 labs(y = bquote(d[x]))

# now we add the theoretical results, using H&P law
l = nsim * c(1, pmat)
d = -diff(l)

d_theo <- data.frame(
             x = age:110,
             theo = d
)

g_sim_HP <- g_sim_HP + 
                 geom_step(data = d_theo, aes(x, theo), 
                           color = "red")

g_sim_HP <- g_sim_HP + 
                 geom_step(data = Belgium_male_2018, 
                           aes(Age, dx), 
                           color = "blue")

g_sim_HP



