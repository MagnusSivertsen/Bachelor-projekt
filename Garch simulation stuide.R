library(qrmtools)
library(sn)
library("rugarch")
library("fGarch")
library("MASS")
library("expectreg")
library("ismev")
library("lmom")
library("QRM")
library("skewt")
library(ggplot2)
library(dplyr)
alpha <- 0.025
n <- 250

# Funktioner, som udregner Var og ES.
Var_garch <- function(sigma,mu,qmodel){
  Var <- numeric(length(sigma))
  for (i in 1:length(sigma)){
    Var[i] <- mu[i] - sigma[i]*qmodel
  }
  return(Var)
}

ES_garch <- function(sigma,mu,ESmodel){
  ES <- numeric(length(sigma))
  for (i in 1:length(sigma)){
    ES[i] <- mu[i] - sigma[i]*ESmodel
  }
  return(ES)
}


Var_st <- function(df){-qstd(alpha,nu=df)}

ES_st <- function(df){(dstd((-1)*Var_st(df),nu=df)*(df+(-Var_st(df))^2)/(alpha*(df-1)))}

ES_t_standard <- function(alpha, df) {
  c <- sqrt(df / (df - 2))
  q_t <- qt(alpha, df)
  q_std <- q_t / c
  f_std_t <- function(z) {
    t_val <- z * c
    dt(t_val, df) * c  
  }
  
  integrand <- function(z) z * f_std_t(z)
  integral_val <- integrate(integrand, lower = -Inf, upper = q_std)$value
  ES <- (1 / alpha) * integral_val
  return(-ES)
}

# denne funktion udregner power ved de forskellige test.
PowerF <- function(matrix,Siglevel){
  quantile1 <- quantile(matrix[,1],Siglevel)
  powers <- numeric(ncol(matrix)-1)
  for (i in 2:ncol(matrix)){
    powers[i-1] <- ecdf(matrix[,i])(quantile1)
  }
  return(powers)
}


# denne funktion udregner hver test størrelse
Z_2 <- function(data,Var,ES){
  I_t <- ifelse(data < Var, 1, 0)
  sum(((data*I_t)/(n*alpha*-ES)))+1
}


Z_3 <- function(X,mu,sigma,nu){
  N <- length(X)
  Ta <- floor(n*alpha)
  
  Ut <- pdist("std", X, mu, sigma, shape=nu)
  
  EST <- numeric(N)
  for (i in 1:N){
    Y <- qdist("std", Ut, mu[i], sigma[i], shape=nu)
    EST[i] <- -1/Ta*sum(sort(Y)[1:Ta])
  }
  
  betamu <- integrate(function(p){pbeta(1-p,n-Ta,Ta)}, lower = 0, upper = 1)$value*mu
  betasigma <- integrate(function(p){pbeta(1-p,n-Ta,Ta)*qdist("std", p, shape=nu)}, 
                         lower = 0, upper = 1)$value*sigma
  
  EV <- -n/Ta*(betamu+betasigma)
  
  -mean(EST/EV)+1
}



Z_ES <- function(data,Var,ES){
  sum((ES-Var-(1/alpha)*pmin(data-Var,0))/ES)/n
}


Var_ex <- function(data,Var)
  -sum(data < -Var)

#Funktion som genererer en sti fra AR(1)-GARCH(1,1) model

simDGP <- function(N,mu,ar1,omega,al,be,Z_dist,seed){
  
  burn=1000  
  n=N+burn # desired length series plus the burn-in period
  
  # Series of innovations
  set.seed(seed=seed)
  x=Z_dist
  # plot(density(x))
  
  # AR(1)-GARCH(1,1) filter
  r=vector(mode="logical", length=n) # process, representing negated log-returns
  mu.t=vector(mode="logical", length=n)
  sig.t=vector(mode="logical", length=n)
  
  # starting values
  r[1]<-0; mu.t[1]<-0; sig.t[1]<-0
  
  for (i in 2:n)
  {
    mu.t[i]=mu+ar1*r[i-1]
    sig.t[i]=sqrt(omega+al*(r[i-1]-mu.t[i-1])^2+be*sig.t[i-1]^2)
    r[i]=mu.t[i]+sig.t[i]*x[i]
  }
  
  simdat = r[-(1:burn)]
  mut = mu.t[-(1:burn)]
  sigt = sig.t[-(1:burn)]
  # plot(simdat, type="l", ylab="")
  
  return(list(xt=simdat,mut=mut,sigt=sigt,burn = c(mut[burn],sigt[burn])))
}


#Udregning af fordelingerne af hver test under h_0 og h_1 ved forskellige df


start_time <- Sys.time()
M=1000
Var_array_garch1 <- matrix(0,M,4)
Z_2_array_garch1 <- matrix(0,M,4)
Z_3_array_garch1 <- matrix(0,M,4)
Z_ES_array_garch1<-matrix(0,M,4)
s<- ES_t_integral(alpha,100)
for (i in 1:M){
  data1 <- simDGP(250,-0.05,0.3,0.01,0.1,0.85,rstd(1250,nu =100),i)
  Var1 <- Var_garch(data1$sigt,data1$mut,Var_st(100))
  ES1 <- ES_garch(data1$sigt,data1$mut,s)
  Z_2_array_garch1[i,1] <- Z_2(data1$xt,Var1,ES1)
  Z_3_array_garch1[i,1] <- Z_3(data1$xt,data1$mut,data1$sigt,100)
  Z_ES_array_garch1[i,1] <- Z_ES(data1$xt,Var1,ES1)
  Var_array_garch1[i,1] <- Var_ex(data1$xt,-Var1)
}

df <- c(10,5,3)
for (j in 1:length(df)){
  for (i in 1:M){
    data1 <- simDGP(250,-0.05,0.3,0.01,0.1,0.85,rstd(1250,nu=df[j]),i)
    Var1 <- Var_garch(data1$sigt,data1$mut,Var_st(100))
    ES1 <- ES_garch(data1$sigt,data1$mut,s)
    Z_2_array_garch1[i,j+1] <- Z_2(data1$xt,Var1,ES1)
    Z_3_array_garch1[i,j+1] <- Z_3(data1$xt,data1$mut,data1$sigt,100)
    Z_ES_array_garch1[i,j+1] <- Z_ES(data1$xt,Var1,ES1)
    Var_array_garch1[i,j+1] <- Var_ex(data1$xt,-Var1)
  }
}
end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)

#Funktion som simulerer sti og alternative sigma_t ved modelmisspecifikation, af forholdet mellem alpha og beta

simDGP_2 <- function(N,mu,ar1,omega,al,al_prime,be,be_prime,Z_dist,seed){
  
  burn=1000  
  n=N+burn # desired length series plus the burn-in period
  
  # Series of innovations
  set.seed(seed=seed)
  x=Z_dist
  # plot(density(x))
  
  # AR(1)-GARCH(1,1) filter
  r=vector(mode="logical", length=n) # process, representing negated log-returns
  mu.t=vector(mode="logical", length=n)
  sig.t=vector(mode="logical", length=n)
  sig.t_prime = vector(mode="logical", length=n)
  # starting values
  r[1]<-0; mu.t[1]<-0; sig.t[1]<-0; sig.t_prime[1] <- 0;
  
  for (i in 2:n)
  {
    mu.t[i]=mu+ar1*r[i-1]
    sig.t[i]=sqrt(omega+al*(r[i-1]-mu.t[i-1])^2+be*sig.t[i-1]^2)
    r[i]=mu.t[i]+sig.t[i]*x[i]
    sig.t_prime[i] = sqrt(omega+al_prime*(r[i-1]-mu.t[i-1])^2+be_prime*sig.t_prime[i-1]^2)
  }
  
  simdat = r[-(1:burn)]
  mut = mu.t[-(1:burn)]
  sigt = sig.t[-(1:burn)]
  sigt_prime = sig.t_prime[-(1:burn)]
  # plot(simdat, type="l", ylab="")
  
  return(list(xt=simdat,mut=mut,sigt=sigt,sigt_prime = sigt_prime))
}

#Udregning af fordelingerne af hver test under h_0 og h_1 ved forskellig alpha og beta


start_time <- Sys.time()
M=1000
Z_2_array_garch2 <- matrix(0,M,4)
Z_2_array_garch2[,1] <- Z_2_array_garch1[,1]
Z_3_array_garch2 <- matrix(0,M,4)
Z_3_array_garch2[,1] <- Z_ES_array_garch1[,1]
Z_ES_array_garch2 <- matrix(0,M,4)
Z_ES_array_garch2[,1] <- Z_ES_array_garch1[,1]
alphas <- c(0.12,0.15,0.2)
for (j in 1:length(alphas)){
  for (i in 1:M){
    data1 <- simDGP_2(250,-0.05,0.3,0.01,alphas[j],0.1,0.95-alphas[j],0.85,rstd(1250,nu=100),i)
    Var1 <- Var_garch(data1$sigt_prime,data1$mut,Var_st(100))
    ES1 <- ES_garch(data1$sigt_prime,data1$mut,ES_t_integral(alpha,100))
    Z_2_array_garch2[i,j+1] <- Z_2(data1$xt,Var1,ES1)
    Z_3_array_garch2[i,j+1] <- Z_3(data1$xt,data1$mut,data1$sigt_prime,100)
    Z_ES_array_garch2[i,j+1] <- Z_ES(data1$xt,Var1,ES1)
  }
}
end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)

#Udregning af fordelingerne af hver test under h_0 og h_1 ved forskellig alphas og betas
start_time <- Sys.time()
M=1000
Z_2_array_garch4 <- matrix(0,M,4)
Z_2_array_garch4[,1] <- Z_2_array_garch1[,1]
Z_3_array_garch4 <- matrix(0,M,4)
Z_3_array_garch4[,1] <- Z_3_array_garch1[,1]
Z_ES_array_garch4 <- matrix(0,M,4)
Z_ES_array_garch4[,1] <- Z_3_array_garch1[,1]
alphas <- c(0.3,0.4,0.5)
for (j in 1:length(alphas)){
  for (i in 1:M){
    data1 <- simDGP_2(250,-0.05,0.3,0.01,alphas[j],0.1,0.95-alphas[j],0.85,rt(1250,100),i)
    Var1 <- Var_garch(data1$sigt_prime,data1$mut,Var_st(100))
    ES1 <- ES_garch(data1$sigt_prime,data1$mut,ES_st(100))
    Z_2_array_garch4[i,j+1] <- Z_2(data1$xt,Var1,ES1)
    Z_3_array_garch4[i,j+1] <- Z_3(data1$xt,data1$mut,data1$sigt_prime,100)
    Z_ES_array_garch4[i,j+1] <- Z_ES(data1$xt,Var1,ES1)
  }
}
end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)
#Udregning af fordelingerne af hver test under h_0 og h_1 ved forskellig alphas og betas


start_time <- Sys.time()
M=1000
Z_2_array_garch5 <- matrix(0,M,4)
Z_2_array_garch5[,1] <- Z_2_array_garch1[,1]
Z_3_array_garch5 <- matrix(0,M,4)
Z_3_array_garch5[,1] <- Z_3_array_garch1[,1]
Z_ES_array_garch5 <- matrix(0,M,4)
Z_ES_array_garch5[,1] <- Z_3_array_garch1[,1]
alphas <- c(0.05,0.03,0.01)
for (j in 1:length(alphas)){
  for (i in 1:M){
    data1 <- simDGP_2(250,-0.05,0.3,0.01,alphas[j],0.1,0.95-alphas[j],0.85,rt(1250,100),i)
    Var1 <- Var_garch(data1$sigt_prime,data1$mut,Var_st(100))
    ES1 <- ES_garch(data1$sigt_prime,data1$mut,ES_st(100))
    Z_2_array_garch5[i,j+1] <- Z_2(data1$xt,Var1,ES1)
    Z_3_array_garch5[i,j+1] <- Z_3(data1$xt,data1$mut,data1$sigt_prime,100)
    Z_ES_array_garch5[i,j+1] <- Z_ES(data1$xt,Var1,ES1)
  }
}
end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)


# Genererer sti ved misspec af omega

simDGP_3 <- function(N,mu,ar1,omega,omega_prime,al,be,Z_dist,seed){
  
  burn=1000  
  n=N+burn # desired length series plus the burn-in period
  
  # Series of innovations
  set.seed(seed=seed)
  x=Z_dist
  # plot(density(x))
  
  # AR(1)-GARCH(1,1) filter
  r=vector(mode="logical", length=n) # process, representing negated log-returns
  mu.t=vector(mode="logical", length=n)
  sig.t=vector(mode="logical", length=n)
  sig.t_prime2 = vector(mode="logical", length=n)
  # starting values
  r[1]<-0; mu.t[1]<-0; sig.t[1]<-0; sig.t_prime2[1] <- 0;
  
  for (i in 2:n)
  {
    mu.t[i]=mu+ar1*r[i-1]
    sig.t[i]=sqrt(omega+al*(r[i-1]-mu.t[i-1])^2+be*sig.t[i-1]^2)
    r[i]=mu.t[i]+sig.t[i]*x[i]
    sig.t_prime2[i] = sqrt(omega_prime+al*(r[i-1]-mu.t[i-1])^2+be*sig.t_prime2[i-1]^2)
  }
  
  simdat = r[-(1:burn)]
  mut = mu.t[-(1:burn)]
  sigt = sig.t[-(1:burn)]
  sigt_prime2 = sig.t_prime2[-(1:burn)]
  # plot(simdat, type="l", ylab="")
  
  return(list(xt=simdat,mut=mut,sigt=sigt,sigt_prime2 = sigt_prime2))
}

#Udregning af fordelingerne af hver test under h_0 og h_1 ved forskellig omega

M=1000
Z_2_array_garch3 <- matrix(0,M,4)
Z_2_array_garch3[,1] <- Z_2_array_garch1[,1]
Z_3_array_garch3 <- matrix(0,M,4)
Z_3_array_garch3[,1] <- Z_3_array_garch1[,1]
Z_ES_array_garch3 <- matrix(0,M,4)
Z_ES_array_garch3[,1] <- Z_3_array_garch1[,1]
omegas <- c(0.015,0.03,0.045)
for (j in 1:length(omegas)){
  for (i in 1:M){
    data1 <- simDGP_3(250,-0.05,0.3,omegas[j],0.01,0.1,0.85,rstd(1250,nu =100),i)
    Var1 <- Var_garch(data1$sigt_prime2,data1$mut,Var_st(100))
    ES1 <- ES_garch(data1$sigt_prime2,data1$mut,s)
    Z_2_array_garch3[i,j+1] <- Z_2(data1$xt,Var1,ES1)
    Z_3_array_garch3[i,j+1] <- Z_3(data1$xt,data1$mut,data1$sigt_prime2,100)
    Z_ES_array_garch3[i,j+1] <- Z_ES(data1$xt,Var1,ES1)
  }
}

# Udregning af power for alle test og misspec

print("df")
PowerF(Z_2_array_garch1,0.05);PowerF(Z_3_array_garch1,0.05);PowerF(Z_ES_array_garch1,0.05)
print("alpha")
PowerF(Z_2_array_garch2,0.05);PowerF(Z_3_array_garch2,0.05);PowerF(Z_ES_array_garch2,0.05)
print("omega")
PowerF(Z_2_array_garch3,0.05);PowerF(Z_3_array_garch3,0.05);PowerF(Z_ES_array_garch3,0.05)


end_time1 <- Sys.time()
runtime1 <- end_time1 - start_time1
print(runtime1)

# genererer sti med misspec af persistens.

simDGP_4 <- function(N,mu,ar1,omega,omega_prime,al,al_prime,be,be_prime,Z_dist,seed){
  
  burn=1000  
  n=N+burn # desired length series plus the burn-in period
  
  # Series of innovations
  set.seed(seed=seed)
  x=Z_dist
  # plot(density(x))
  
  # AR(1)-GARCH(1,1) filter
  r=vector(mode="logical", length=n) # process, representing negated log-returns
  mu.t=vector(mode="logical", length=n)
  sig.t=vector(mode="logical", length=n)
  sig.t_prime2 = vector(mode="logical", length=n)
  # starting values
  r[1]<-0; mu.t[1]<-0; sig.t[1]<-0; sig.t_prime2[1] <- 0;
  
  for (i in 2:n)
  {
    mu.t[i]=mu+ar1*r[i-1]
    sig.t[i]=sqrt(omega+al*(r[i-1]-mu.t[i-1])^2+be*sig.t[i-1]^2)
    r[i]=mu.t[i]+sig.t[i]*x[i]
    sig.t_prime2[i] = sqrt(omega_prime+al_prime*(r[i-1]-mu.t[i-1])^2+be_prime*sig.t_prime2[i-1]^2)
  }
  
  simdat = r[-(1:burn)]
  mut = mu.t[-(1:burn)]
  sigt = sig.t[-(1:burn)]
  sigt_prime2 = sig.t_prime2[-(1:burn)]
  # plot(simdat, type="l", ylab="")
  
  return(list(xt=simdat,mut=mut,sigt=sigt,sigt_prime2 = sigt_prime2))
}
#Udregning af fordelingerne af hver test under h_0 og h_1 ved forskellig persistens

M=1000
Z_2_array_garch10 <- matrix(0,M,5)
Z_3_array_garch10 <- matrix(0,M,5)
Z_ES_array_garch10 <- matrix(0,M,5)

d_par <- c(0.8/0.9,0.85/0.9,0.95/0.9,0.99/0.9) 
alphas1 <- d_par * 0.1
betas1 <- d_par *0.8
E_sig <- 0.01/(1-0.8-0.1)
omegas1 <- E_sig*(1-betas1-alphas1)
for (i in 1:M){
  data1 <- simDGP(250,-0.05,0.3,0.01,0.1,0.8,rstd(1250,nu =100),i)
  Var1 <- Var_garch(data1$sigt,data1$mut,Var_st(100))
  ES1 <- ES_garch(data1$sigt,data1$mut,s)
  Z_2_array_garch10[i,1] <- Z_2(data1$xt,Var1,ES1)
  Z_3_array_garch10[i,1] <- Z_3(data1$xt,data1$mut,data1$sigt,100)
  Z_ES_array_garch10[i,1] <- Z_ES(data1$xt,Var1,ES1)
}



for (j in 1:4){
  for (i in 1:M){
    data1 <- simDGP_4(250,-0.05,0.3,omegas1[j],0.01,alphas1[j],0.1,betas1[j],0.8,rstd(1250,nu =100),i)
    Var1 <- Var_garch(data1$sigt_prime2,data1$mut,Var_st(100))
    ES1 <- ES_garch(data1$sigt_prime2,data1$mut,s)
    Z_2_array_garch10[i,j+1] <- Z_2(data1$xt,Var1,ES1)
    Z_3_array_garch10[i,j+1] <- Z_3(data1$xt,data1$mut,data1$sigt_prime2,100)
    Z_ES_array_garch10[i,j+1] <- Z_ES(data1$xt,Var1,ES1)
  }
}

# power plot af misspec af alpha og beta

x_vals <- c(0.05,0.03,0.01,0.1,0.12,0.15,0.2,0.3,0.4,0.5)

data <- data.frame(
  alpha = rep(x_vals, 3),
  power = c( c((PowerF(Z_2_array_garch5,0.05)),0.05,PowerF(Z_2_array_garch2,0.05),PowerF(Z_2_array_garch4,0.05)), c((PowerF(Z_3_array_garch5,0.05)),0.05,PowerF(Z_3_array_garch2,0.05),PowerF(Z_3_array_garch4,0.05)), c((PowerF(Z_ES_array_garch5,0.05)),0.05,PowerF(Z_ES_array_garch2,0.05),PowerF(Z_ES_array_garch4,0.05))),
  test = rep(c("Z_2", "Z3","Z_ES"), each = 10)
)

ggplot(data, aes(x = alpha, y = power, color = test)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_continuous(name = "Alpha", limits = c(0.01, 0.5)) +
  scale_y_continuous(name = "Power", limits = c(0, 0.3)) +
  geom_vline(xintercept = 0.1, linetype = "solid", color = "black", size = 1) +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  labs(color = "Test", title = "")

# Power plot af misspec af persistens. 

x_vals <- c(d_par * 0.9,0.9)

data <- data.frame(
  alpha = rep(x_vals, 3),
  power = c( c((PowerF(Z_2_array_garch10,0.05)),0.05), c((PowerF(Z_3_array_garch10,0.05)),0.05), c((PowerF(Z_ES_array_garch10,0.05)),0.05)),
  test = rep(c("Z_2", "Z_3","Z_ES"), each = 5)
)


ggplot(data, aes(x = alpha, y = power, color = test)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_continuous(name = "Persistance", limits = c(0.8, 1)) +
  scale_y_continuous(name = "Power", limits = c(0, 0.2)) +
  geom_vline(xintercept = 0.9, linetype = "solid", color = "black", size = 1) +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  labs(color = "Test", title = "")

# Udregning af power igen

print("df")
PowerF(Z_2_array_garch1,0.05);PowerF(Z_3_array_garch1,0.05);PowerF(Z_ES_array_garch1,0.05)
print("alpha")
PowerF(Z_2_array_garch2,0.05);PowerF(Z_3_array_garch2,0.05);PowerF(Z_ES_array_garch2,0.05)
print("omega")
PowerF(Z_2_array_garch3,0.05);PowerF(Z_3_array_garch3,0.05);PowerF(Z_ES_array_garch3,0.05)



simulated_series7 <- simDGP(250,-0.05,0.3,0.01,0.1,0.85,rstd(1250,nu=20),2)
simulated_series8 <- simDGP(250,-0.05,0.3,0.01,0.1,0.85,rstd(1250,20),2)
Alpha = 0.025
Var2 <- ES_garch(simulated_series7$sigt,simulated_series7$mut,Var_st(20))

df <- data.frame(
  Time = seq_along(simulated_series7$xt),
  X_t = simulated_series7$xt,
  ES = Var2
)

ggplot(df, aes(x = Time)) +
  geom_line(aes(y = X_t, color = "X_t")) +
  geom_line(aes(y = ES, color = "ES")) +
  scale_color_manual(values = c("X_t" = "blue", "ES" = "red")) +
  theme_minimal() +
  labs(y = "Value", color = "Series", title = "X_t and Expected Shortfall Over Time")




end_time1 <- Sys.time()
runtime1 <- end_time1 - start_time1
print(runtime1)