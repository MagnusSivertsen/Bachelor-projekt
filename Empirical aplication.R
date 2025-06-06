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
library("quantmod")
alpha <- 0.025
n <- 250
set.seed(1964)

# Funktioner til at udregne Var og ES
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

ES_normal <- dnorm(qnorm(alpha))/alpha


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
# funktion, som simulerer en garch sti, som afhænger af de refittede modeller.

simDGP2 <- function(N,fit,seed){
  
  burn=1000  
  n=N+burn # desired length series plus the burn-in period
  
  fit1 <- fit[[1]]
  par <- coef(fit1)
  mu <- par["mu"];ar1 <- par["ar1"];omega <- par["omega"];alpha <- par["alpha1"];beta <- par["beta1"]
  set.seed(seed=seed)
  # Series of innovations
  x=rstd(burn,nu=par["shape"])
  
  r=numeric(n) # process, representing negated log-returns
  mu.t=numeric(n)
  sig.t=numeric(n)
  
  # starting values
  r[1]<-0; mu.t[1]<-0; sig.t[1]<-0
  
  for (i in 2:burn)
  {
    mu.t[i]=mu+ar1*r[i-1]
    sig.t[i]=sqrt(omega+alpha*(r[i-1]-mu.t[i-1])^2+beta*sig.t[i-1]^2)
    r[i]=mu.t[i]+sig.t[i]*x[i]
  }
  chunk_size <- N/(length(fit))
  for (k in 1:length(fit)){
    par <- coef(fit[[k]])
    mu <- par["mu"];ar1 <- par["ar1"];omega <- par["omega"];alpha <- par["alpha1"];beta <- par["beta1"]
    x1 <- rstd(chunk_size, 0, 1, nu = par["shape"])
    
    start_idx <- burn + 1 + (k - 1) * chunk_size
    end_idx <- burn + k * chunk_size
    
    for (i in start_idx:end_idx) {
      mu.t[i] <- mu + ar1 * r[i - 1]
      sig.t[i] <- sqrt(omega + alpha * (r[i - 1] - mu.t[i - 1])^2 + beta * sig.t[i - 1]^2)
      r[i] <- mu.t[i] + sig.t[i] * x1[i - start_idx + 1]
      
    }
  }
  
  simdat = r[-(1:burn)]
  mut = mu.t[-(1:burn)]
  sigt = sig.t[-(1:burn)]
  # plot(simdat, type="l", ylab="")
  
  return(list(xt=simdat,mut=mut,sigt=sigt,burn = c(mut[burn],sigt[burn])))
}



#Samme blot for normalfordelte innovationer


simDGP2norm <- function(N,fit,seed){
  
  burn=1000  
  n=N+burn # desired length series plus the burn-in period
  
  fit1 <- fit[[1]]
  par <- coef(fit1)
  mu <- par["mu"];ar1 <- par["ar1"];omega <- par["omega"];alpha <- par["alpha1"];beta <- par["beta1"]
  set.seed(seed=seed)
  # Series of innovations
  x=rnorm(burn,0,1)
  
  r=numeric(n) # process, representing negated log-returns
  mu.t=numeric(n)
  sig.t=numeric(n)
  
  # starting values
  r[1]<-0; mu.t[1]<-0; sig.t[1]<-0
  
  for (i in 2:burn)
  {
    mu.t[i]=mu+ar1*r[i-1]
    sig.t[i]=sqrt(omega+alpha*(r[i-1]-mu.t[i-1])^2+beta*sig.t[i-1]^2)
    r[i]=mu.t[i]+sig.t[i]*x[i]
  }
  chunk_size <- N/(length(fit))
  for (k in 1:length(fit)){
    par <- coef(fit[[k]])
    mu <- par["mu"];ar1 <- par["ar1"];omega <- par["omega"];alpha <- par["alpha1"];beta <- par["beta1"]
    x1 <- rnorm(chunk_size, 0, 1)
    
    start_idx <- burn + 1 + (k - 1) * chunk_size
    end_idx <- burn + k * chunk_size
    
    for (i in start_idx:end_idx) {
      mu.t[i] <- mu + ar1 * r[i - 1]
      sig.t[i] <- sqrt(omega + alpha * (r[i - 1] - mu.t[i - 1])^2 + beta * sig.t[i - 1]^2)
      r[i] <- mu.t[i] + sig.t[i] * x1[i - start_idx + 1]
      
    }
  }
  
  simdat = r[-(1:burn)]
  mut = mu.t[-(1:burn)]
  sigt = sig.t[-(1:burn)]
  # plot(simdat, type="l", ylab="")
  
  return(list(xt=simdat,mut=mut,sigt=sigt,burn = c(mut[burn],sigt[burn])))
}

#Backtest funktionerne

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

# Denne funktion inhenter aktiekurser for en given ticker

Importlog <- function(ticker,start,finish) {
  FinData <- getSymbols(ticker,src = "yahoo",from = start, to = finish, auto.assign = FALSE)
  LogR <- dailyReturn(Ad(FinData)[paste0(start, "::", finish)], type = "log")
  return(LogR)
}

# denne funktion tager en liste af it, som input, og tager Out of sample data, som input. 
# Det som funktionen gør er, at den tager mu_t og sigma_t fra fit, og rekursivt udregner mu_t+1 og sigma_t+1
# Dette gør den rekursivt indtil den når en dag, hvor der bliver refittet. Derefter skifter den parametre
# på AR-Garch modellen og kører videre. Hver dag, bruger den en observation/realisation fra OOS_data, 
# til at udregne næste dags mu og sigma.

garch_ahead <- function(fit,OOSdata,NumYear){
  a <- length(fit)
  n <- as.integer((250*NumYear)/a)
  
  OOSdata <- head(as.numeric(OOSdata),250*NumYear)
  mu_vector <- numeric(250*NumYear)
  sigma_vector <- numeric(250*NumYear)
  
  for (k in 1:(length(fit))){
    fit1 <- fit[[k]]
    xt_minus_1 <- tail(fit1@model$modeldata$data,1)
    sigma_fit <- tail(sigma(fit1),1)
    resid_fit <- tail(residuals(fit1),1)
    par <- coef(fit1)
    mu <- par["mu"];ar1 <- par["ar1"];omega <- par["omega"];alpha <- par["alpha1"];beta <- par["beta1"]
    
    for (i in 1:n){
      mu_t <- mu + ar1 * xt_minus_1
      sigma_t <- sqrt(omega + alpha * resid_fit^2+beta*sigma_fit^2)
      mu_vector[i+n*(k-1)] <- mu_t
      sigma_vector[i+n*(k-1)] <- sigma_t
      
      xt_minus_1 <- OOSdata[i+n*(k-1)]
      sigma_fit <- sigma_t
      resid_fit <- OOSdata[i+n*(k-1)] - mu_t
      
    }
  }  
  return(list(mu=mu_vector,sigma=sigma_vector))
}

#Denne funktion tager in-sample data, model specifikation, refit frekvens, længde af rolling window, og hvilket år, som skal være OOS.
# Det, som funktionen gør er, at den opsætter indeksne for de rullende vindue. Derefter fitter den en garchmodel
# som skal bruges til at forecaste Var og ES indtil næste refit. Så forskyder den vinduet og fitter igen.
# Output er en liste af de Garch modeller, som bliver fittet for  hvert rolling window. 
  
  MW.fit2 <- function(datas,spec,freq,lenWin,OOSYear,NumYear){
    xt <- as.numeric(datas)
    xtOOS <- xt[OOSYear]
    fit <- list()
    moving_window1 <- list()
    # definering af moving window
    
    years <- sort(unique(format(index(datas), "%Y"))) ; end_year <- years[length(years) - NumYear] ; end_date <- last(index(datas[format(index(datas), "%Y") == end_year])) ; 
    end_index <- which(index(datas) == end_date)
    for (i in 1:((250*NumYear)/freq)){
      start_date <- index(datas[end_index-lenWin + (i-1)*(freq)])
      end_date <- index(datas[end_index+ (i-1)*(freq)])
      moving_window <- datas[paste0(start_date, "/", end_date)]
      moving_window1[[i]] <- moving_window
      fit[[i]] <- ugarchfit(spec,moving_window,solver= "hybrid")
    }
    return(list(fit,moving_window1))
  }
  
  tickers <- c("AAPL", "MSFT", "AMZN", "GOOGL", "GOOG", "BRK-B", "NVDA", "META", "TSLA", "UNH",
               "JNJ", "V", "XOM", "JPM", "PG", "MA", "HD", "CVX", "LLY", "MRK",
               "ABBV", "PEP", "KO", "AVGO", "COST", "WMT", "BAC", "DIS", "ADBE", "PFE",
               "NFLX", "TMO", "INTC", "CSCO", "VZ", "ABT", "CRM", "ACN", "MCD", "DHR",
               "NKE", "TXN", "LIN", "ORCL", "QCOM", "PM", "NEE", "MDT", "UPS", "UNP")
  
# Denne funktion importere data, fitter garchmodeller med MW.fit2, og derefter forecaster VaR og ES.
# Den udregner så den realiseret Z(X) og derefter simulerer under h_0 for hver test, og outputter den 
# realiserede Z(x) og fordelingerne af Z_2 og Z_ES under h_0,
# den udregner også p-værdien for begge test.
  
Emp_test4 <- function(tick,M,freq,NumYear){
  data2 <- Importlog(tick,"2013-01-01","2025-01-01")
  years <- sort(unique(format(index(data2), "%Y"))) ; end_year <- years[length(years) - NumYear]
  end_date <- last(index(data2[format(index(data2), "%Y") == end_year])) ; 
    end_index <- which(index(data2) == end_date)
  spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T),distribution.model = "std")
  test1 <- MW.fit2(data2,spec,freq,2000,"2024",NumYear)
  OOS <- garch_ahead(test1[[1]],data2[(end_index+1):(end_index+(250*NumYear))],NumYear)
    
  Var_1 <- numeric(250*NumYear)
  ES_1 <- numeric(250*NumYear)
    
  for (k in 1:((250*NumYear)/freq)){
    Var_1[(1+(k-1)*freq):(freq*k)] <- Var_garch(OOS[[2]][(1+(k-1)*freq):(freq*k)], OOS[[1]][(1+(k-1)*freq):(freq*k)],Var_st(coef(test1[[1]][[k]])["shape"]))
    ES_1[(1+(k-1)*freq):(freq*k)] <- ES_garch(OOS[[2]][(1+(k-1)*freq):(freq*k)], OOS[[1]][(1+(k-1)*freq):(freq*k)],ES_t_standard(alpha,coef(test1[[1]][[k]])["shape"]))
    }
  
    
      
  
  Z_2emp <- Z_2(data2[(end_index+1):(end_index+(250*NumYear))],Var_1,ES_1)
  Z_ESemp <- Z_ES(data2[(end_index+1):(end_index+(250*NumYear))],Var_1,ES_1)
      
    
    
    
  Z_ES_test_array <- matrix(0,M,1)
  Z_2_test_array <- matrix(0,M,1)
  for (i in 1:M){
    data1 <- simDGP2(250*NumYear,test1[[1]],i)
    Var_2 <- numeric(250*NumYear)
    ES_2 <- numeric(250*NumYear)
      for (k in 1:((250*NumYear)/freq)){
        Var_2[(1+(k-1)*freq):(freq*k)] <- Var_garch(data1$sigt[(1+(k-1)*freq):(freq*k)], data1$mut[(1+(k-1)*freq):(freq*k)],Var_st(coef(test1[[1]][[k]])["shape"]))
        ES_2[(1+(k-1)*freq):(freq*k)] <- ES_garch(data1$sigt[(1+(k-1)*freq):(freq*k)], data1$mut[(1+(k-1)*freq):(freq*k)],ES_t_standard(alpha,coef(test1[[1]][[k]])["shape"]))
    }
      
  
    Z_ES_test_array[i,1] <- Z_ES(data1$xt,Var_2,ES_2)
    Z_2_test_array[i,1] <- Z_2(data1$xt,Var_2,ES_2)
  }
    
  p_Z_2 <- sum(Z_2_test_array < Z_2emp)/M
  mean_Z_2 <- mean(Z_2_test_array)
    
  p_Z_ES <- sum(Z_ES_test_array < Z_ESemp)/M
  mean_Z_ES <- mean(Z_ES_test_array)
    
  df <- data.frame(
    p_value = c(p_Z_2,p_Z_ES),
    mean = c(mean_Z_2,mean_Z_ES))
  rownames(df) <- c("Z_2","Z_ES")
  return(list(df,Z_2_test_array,Z_ES_test_array,Z_2emp,Z_ESemp))

}
# Dette er samme funktion som Emptest4 blot med normalfordelte innovationer i garch modellen.

ES_normal <- dnorm(qnorm(alpha))/alpha
Emp_test5 <- function(tick,M,freq,NumYear){
  data2 <- Importlog(tick,"2013-01-01","2025-01-01")
  years <- sort(unique(format(index(data2), "%Y"))) ; end_year <- years[length(years) - NumYear]
  end_date <- last(index(data2[format(index(data2), "%Y") == end_year])) ; 
  end_index <- which(index(data2) == end_date)
  spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T),distribution.model = "norm")
  test1 <- MW.fit2(data2,spec,freq,2500,"2024",NumYear)
  OOS <- garch_ahead(test1[[1]],data2[(end_index+1):(end_index+(250*NumYear))],NumYear)
  
  Var_1 <- numeric(250*NumYear)
  ES_1 <- numeric(250*NumYear)
  
  for (k in 1:((250*NumYear)/freq)){
    Var_1[(1+(k-1)*freq):(freq*k)] <- Var_garch(OOS[[2]][(1+(k-1)*freq):(freq*k)], OOS[[1]][(1+(k-1)*freq):(freq*k)],-qnorm(alpha,0,1))
    ES_1[(1+(k-1)*freq):(freq*k)] <- ES_garch(OOS[[2]][(1+(k-1)*freq):(freq*k)], OOS[[1]][(1+(k-1)*freq):(freq*k)],ES_normal)
  }

  
  
  
  Z_2emp <- Z_2(data2[(end_index+1):(end_index+(250*NumYear))],Var_1,ES_1)
  Z_ESemp <- Z_ES(data2[(end_index+1):(end_index+(250*NumYear))],Var_1,ES_1)
  
  
  
  
  Z_ES_test_array <- matrix(0,M,1)
  Z_2_test_array <- matrix(0,M,1)
  for (i in 1:M){
    data1 <- simDGP2norm(250*NumYear,test1[[1]],i)
    Var_2 <- numeric(250*NumYear)
    ES_2 <- numeric(250*NumYear)
    for (k in 1:((250*NumYear)/freq)){
      Var_2[(1+(k-1)*freq):(freq*k)] <- Var_garch(data1$sigt[(1+(k-1)*freq):(freq*k)], data1$mut[(1+(k-1)*freq):(freq*k)],-qnorm(alpha,0,1))
      ES_2[(1+(k-1)*freq):(freq*k)] <- ES_garch(data1$sigt[(1+(k-1)*freq):(freq*k)], data1$mut[(1+(k-1)*freq):(freq*k)],ES_normal)
    }
    
    
    Z_ES_test_array[i,1] <- Z_ES(data1$xt,Var_2,ES_2)
    Z_2_test_array[i,1] <- Z_2(data1$xt,Var_2,ES_2)
  }
  
  p_Z_2 <- sum(Z_2_test_array < Z_2emp)/M
  mean_Z_2 <- mean(Z_2_test_array)
  
  p_Z_ES <- sum(Z_ES_test_array < Z_ESemp)/M
  mean_Z_ES <- mean(Z_ES_test_array)
  
  df <- data.frame(
    p_value = c(p_Z_2,p_Z_ES),
    mean = c(mean_Z_2,mean_Z_ES))
  rownames(df) <- c("Z_2","Z_ES")
  return(list(df,Z_2_test_array,Z_ES_test_array,Z_2emp,Z_ESemp))

}

#Resten af scriptet er kørsler af Emp_test4 og 5 med forskellige parametre, for alle aktier i tickers.


Z_ES_emp_array_250_250 <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_ES_emp_array_250_250) <- tickers

Z_2_emp_array_250_250 <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_ES_emp_array_250_250) <- tickers


start_time <- Sys.time()

for (i in 1:length(tickers)){
  cool <- Emp_test4(tickers[i],1000,250,1)[[1]]
  Z_ES_emp_array_250_250[i,1] <- cool[2,1]
  Z_ES_emp_array_250_250[i,2] <- cool[2,2]
  Z_2_emp_array_250_250[i,1] <- cool[1,1]
  Z_2_emp_array_250_250[i,2] <- cool[1,2]
  cat("Iteration", i, "\n")
}

end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)

Z_ES_emp_array_250_250[,1,drop=FALSE]
Z_2_emp_array_250_250[,1,drop=FALSE]

sum(Z_ES_emp_array_250_250[,1]<0.025)/length(tickers)
sum(Z_2_emp_array_250_250[,1]<0.025)/length(tickers)

Z_ES_emp_array_125_250 <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_ES_emp_array_125_250) <- tickers

Z_2_emp_array_125_250 <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_2_emp_array_125_250) <- tickers



start_time <- Sys.time()

for (i in 1:length(tickers)){
  cool <- Emp_test4(tickers[i],1000,125,1)[[1]]
  Z_ES_emp_array_125_250[i,1] <- cool[2,1]
  Z_ES_emp_array_125_250[i,2] <- cool[2,2]
  Z_2_emp_array_125_250[i,1] <- cool[1,1]
  Z_2_emp_array_125_250[i,2] <- cool[1,2]
  cat("Iteration", i, "\n")
}

end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)

Z_ES_emp_array_125_250[,1,drop=FALSE]
Z_2_emp_array_125_250[,1,drop=FALSE]
sum(Z_ES_emp_array_125_250[,1]<0.025)/length(tickers)
sum(Z_2_emp_array_125_250[,1]<0.025)/length(tickers)

Z_ES_emp_array_50_250 <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_ES_emp_array_50_250) <- tickers

Z_2_emp_array_50_250 <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_2_emp_array_50_250) <- tickers




start_time <- Sys.time()

for (i in 1:length(tickers)){
  cool <- Emp_test4(tickers[i],1000,50,1)[[1]]
  Z_ES_emp_array_50_250[i,1] <- cool[2,1]
  Z_ES_emp_array_50_250[i,2] <- cool[2,2]
  Z_2_emp_array_50_250[i,1] <- cool[1,1]
  Z_2_emp_array_50_250[i,2] <- cool[1,2]
  cat("Iteration", i, "\n")
}

end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)

Z_ES_emp_array_50_250[,1,drop=FALSE]
Z_2_emp_array_50_250[,1,drop=FALSE]
sum(Z_ES_emp_array_50_250[,1]<0.025)/length(tickers)
sum(Z_2_emp_array_50_250[,1]<0.025)/length(tickers)

Z_ES_emp_array_5_250 <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_ES_emp_array_5_250) <- tickers

Z_2_emp_array_5_250 <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_2_emp_array_5_250) <- tickers




start_time <- Sys.time()

for (i in 1:length(tickers)){
  cool <- Emp_test4(tickers[i],1000,5,1)[[1]]
  Z_ES_emp_array_5_250[i,1] <- cool[2,1]
  Z_ES_emp_array_5_250[i,2] <- cool[2,2]
  Z_2_emp_array_5_250[i,1] <- cool[1,1]
  Z_2_emp_array_5_250[i,2] <- cool[1,2]
  cat("Iteration", i, "\n")
}

end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)

Z_ES_emp_array_5_250[,1,drop=FALSE]
Z_2_emp_array_5_250[,1,drop=FALSE]
sum(Z_ES_emp_array_5_250[,1]<0.025)/length(tickers)
sum(Z_2_emp_array_5_250[,1]<0.025)/length(tickers)


Z_ES_emp_array_250_500 <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_ES_emp_array_250_500) <- tickers

Z_2_emp_array_250_500 <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_ES_emp_array_250_500) <- tickers


start_time <- Sys.time()

for (i in 1:length(tickers)){
  cool <- Emp_test4(tickers[i],1000,250,2)[[1]]
  Z_ES_emp_array_250_500[i,1] <- cool[2,1]
  Z_ES_emp_array_250_500[i,2] <- cool[2,2]
  Z_2_emp_array_250_500[i,1] <- cool[1,1]
  Z_2_emp_array_250_500[i,2] <- cool[1,2]
  cat("Iteration", i, "\n")
}

end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)

Z_ES_emp_array_250_500[,1,drop=FALSE]
Z_2_emp_array_250_500[,1,drop=FALSE]
sum(Z_ES_emp_array_250_500[,1]<0.025)/length(tickers)
sum(Z_2_emp_array_250_500[,1]<0.025)/length(tickers)

Z_ES_emp_array_125_500 <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_ES_emp_array_125_500) <- tickers

Z_2_emp_array_125_500 <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_2_emp_array_125_500) <- tickers



start_time <- Sys.time()

for (i in 1:length(tickers)){
  cool <- Emp_test4(tickers[i],1000,125,2)[[1]]
  Z_ES_emp_array_125_500[i,1] <- cool[2,1]
  Z_ES_emp_array_125_500[i,2] <- cool[2,2]
  Z_2_emp_array_125_500[i,1] <- cool[1,1]
  Z_2_emp_array_125_500[i,2] <- cool[1,2]
  cat("Iteration", i, "\n")
}

end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)

Z_ES_emp_array_125_500[,1,drop=FALSE]
Z_2_emp_array_125_500[,1,drop=FALSE]
sum(Z_ES_emp_array_125_500[,1]<0.025)/length(tickers)
sum(Z_2_emp_array_125_500[,1]<0.025)/length(tickers)

Z_ES_emp_array_50_500 <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_ES_emp_array_50_500) <- tickers

Z_2_emp_array_50_500 <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_2_emp_array_50_500) <- tickers




start_time <- Sys.time()

for (i in 1:length(tickers)){
  cool <- Emp_test4(tickers[i],1000,50,2)[[1]]
  Z_ES_emp_array_50_500[i,1] <- cool[2,1]
  Z_ES_emp_array_50_500[i,2] <- cool[2,2]
  Z_2_emp_array_50_500[i,1] <- cool[1,1]
  Z_2_emp_array_50_500[i,2] <- cool[1,2]
  cat("Iteration", i, "\n")
}

end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)

Z_ES_emp_array_50_500[,1,drop=FALSE]
Z_2_emp_array_50_500[,1,drop=FALSE]
sum(Z_ES_emp_array_50_500[,1]<0.025)/length(tickers)
sum(Z_2_emp_array_50_500[,1]<0.025)/length(tickers)

Z_ES_emp_array_5_500 <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_ES_emp_array_5_500) <- tickers

Z_2_emp_array_5_500 <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_2_emp_array_5_500) <- tickers




start_time <- Sys.time()

for (i in 1:length(tickers)){
  cool <- Emp_test4(tickers[i],1000,5,2)[[1]]
  Z_ES_emp_array_5_500[i,1] <- cool[2,1]
  Z_ES_emp_array_5_500[i,2] <- cool[2,2]
  Z_2_emp_array_5_500[i,1] <- cool[1,1]
  Z_2_emp_array_5_500[i,2] <- cool[1,2]
  cat("Iteration", i, "\n")
}

end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)

Z_ES_emp_array_5_500[,1,drop=FALSE]
Z_2_emp_array_5_500[,1,drop=FALSE]
sum(Z_ES_emp_array_5_500[,1]<0.025)/length(tickers)
sum(Z_2_emp_array_5_500[,1]<0.025)/length(tickers)

Z_ES_emp_array_250_500n <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_ES_emp_array_250_500n) <- tickers

Z_2_emp_array_250_500n <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_2_emp_array_250_500n) <- tickers




start_time <- Sys.time()

for (i in 1:length(tickers)){
  cool <- Emp_test5(tickers[i],1000,250,2)[[1]]
  Z_ES_emp_array_250_500n[i,1] <- cool[2,1]
  Z_ES_emp_array_250_500n[i,2] <- cool[2,2]
  Z_2_emp_array_250_500n[i,1] <- cool[1,1]
  Z_2_emp_array_250_500n[i,2] <- cool[1,2]
  cat("Iteration", i, "\n")
}

end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)

Z_ES_emp_array_250_500n[,1,drop=FALSE]
Z_2_emp_array_250_500n[,1,drop=FALSE]
sum(Z_ES_emp_array_250_500n[,1]<0.025)/length(tickers)
sum(Z_2_emp_array_250_500n[,1]<0.025)/length(tickers)

Z_ES_emp_array_125_500n <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_ES_emp_array_125_500n) <- tickers

Z_2_emp_array_125_500n <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_2_emp_array_125_500n) <- tickers




start_time <- Sys.time()

for (i in 1:length(tickers)){
  cool <- Emp_test5(tickers[i],1000,125,2)[[1]]
  Z_ES_emp_array_125_500n[i,1] <- cool[2,1]
  Z_ES_emp_array_125_500n[i,2] <- cool[2,2]
  Z_2_emp_array_125_500n[i,1] <- cool[1,1]
  Z_2_emp_array_125_500n[i,2] <- cool[1,2]
  cat("Iteration", i, "\n")
}

end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)

Z_ES_emp_array_125_500n[,1,drop=FALSE]
Z_2_emp_array_125_500n[,1,drop=FALSE]
sum(Z_ES_emp_array_125_500n[,1]<0.025)/length(tickers)
sum(Z_2_emp_array_125_500n[,1]<0.025)/length(tickers)

Z_ES_emp_array_50_500n <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_ES_emp_array_50_500n) <- tickers

Z_2_emp_array_50_500n <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_2_emp_array_50_500n) <- tickers




start_time <- Sys.time()

for (i in 1:length(tickers)){
  cool <- Emp_test5(tickers[i],1000,50,2)[[1]]
  Z_ES_emp_array_50_500n[i,1] <- cool[2,1]
  Z_ES_emp_array_50_500n[i,2] <- cool[2,2]
  Z_2_emp_array_50_500n[i,1] <- cool[1,1]
  Z_2_emp_array_50_500n[i,2] <- cool[1,2]
  cat("Iteration", i, "\n")
}

end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)

Z_ES_emp_array_50_500n[,1,drop=FALSE]
Z_2_emp_array_50_500n[,1,drop=FALSE]
sum(Z_ES_emp_array_50_500n[,1]<0.025)/length(tickers)
sum(Z_2_emp_array_50_500n[,1]<0.025)/length(tickers)

Z_ES_emp_array_5_500n <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_ES_emp_array_5_500n) <- tickers

Z_2_emp_array_5_500n <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_2_emp_array_5_500n) <- tickers



start_time <- Sys.time()

for (i in 1:length(tickers)){
  cool <- Emp_test5(tickers[i],1000,5,2)[[1]]
  Z_ES_emp_array_5_500n[i,1] <- cool[2,1]
  Z_ES_emp_array_5_500n[i,2] <- cool[2,2]
  Z_2_emp_array_5_500n[i,1] <- cool[1,1]
  Z_2_emp_array_5_500n[i,2] <- cool[1,2]
  cat("Iteration", i, "\n")
}

end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)



Z_ES_emp_array_5_500n[,1,drop=FALSE]
Z_2_emp_array_5_500n[,1,drop=FALSE]
sum(Z_ES_emp_array_5_500n[,1]<0.025)/length(tickers)
sum(Z_2_emp_array_5_500n[,1]<0.025)/length(tickers)



xxx




Z_ES_emp_array_250_250n <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_ES_emp_array_250_250n) <- tickers

Z_2_emp_array_250_250n <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_2_emp_array_250_250n) <- tickers






start_time <- Sys.time()

for (i in 1:length(tickers)){
  cool <- Emp_test5(tickers[i],1000,250,1)[[1]]
  Z_ES_emp_array_250_250n[i,1] <- cool[2,1]
  Z_ES_emp_array_250_250n[i,2] <- cool[2,2]
  Z_2_emp_array_250_250n[i,1] <- cool[1,1]
  Z_2_emp_array_250_250n[i,2] <- cool[1,2]
  cat("Iteration", i, "\n")
}

end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)



Z_ES_emp_array_250_250n[,1,drop=FALSE]
Z_2_emp_array_250_250n[,1,drop=FALSE]
sum(Z_ES_emp_array_250_250n[,1]<0.025)/length(tickers)
sum(Z_2_emp_array_250_250n[,1]<0.025)/length(tickers)



Z_ES_emp_array_125_250n <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_ES_emp_array_125_250n) <- tickers

Z_2_emp_array_125_250n <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_2_emp_array_125_250n) <- tickers






start_time <- Sys.time()

for (i in 1:length(tickers)){
  cool <- Emp_test5(tickers[i],1000,125,1)[[1]]
  Z_ES_emp_array_125_250n[i,1] <- cool[2,1]
  Z_ES_emp_array_125_250n[i,2] <- cool[2,2]
  Z_2_emp_array_125_250n[i,1] <- cool[1,1]
  Z_2_emp_array_125_250n[i,2] <- cool[1,2]
  cat("Iteration", i, "\n")
}

end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)



Z_ES_emp_array_125_250n[,1,drop=FALSE]
Z_2_emp_array_125_250n[,1,drop=FALSE]
sum(Z_ES_emp_array_125_250n[,1]<0.025)/length(tickers)
sum(Z_2_emp_array_125_250n[,1]<0.025)/length(tickers)



Z_ES_emp_array_50_250n <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_ES_emp_array_50_250n) <- tickers

Z_2_emp_array_50_250n <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_2_emp_array_50_250n) <- tickers






start_time <- Sys.time()

for (i in 1:length(tickers)){
  cool <- Emp_test5(tickers[i],1000,50,1)[[1]]
  Z_ES_emp_array_50_250n[i,1] <- cool[2,1]
  Z_ES_emp_array_50_250n[i,2] <- cool[2,2]
  Z_2_emp_array_50_250n[i,1] <- cool[1,1]
  Z_2_emp_array_50_250n[i,2] <- cool[1,2]
  cat("Iteration", i, "\n")
}

end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)



Z_ES_emp_array_50_250n[,1,drop=FALSE]
Z_2_emp_array_50_250n[,1,drop=FALSE]
sum(Z_ES_emp_array_50_250n[,1]<0.025)/length(tickers)
sum(Z_2_emp_array_50_250n[,1]<0.025)/length(tickers)




Z_ES_emp_array_5_250n <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_ES_emp_array_5_250n) <- tickers

Z_2_emp_array_5_250n <- data.frame(
  p_value = numeric(length(tickers)),
  mean = numeric(length(tickers))
)
rownames(Z_2_emp_array_5_250n) <- tickers






start_time <- Sys.time()

for (i in 1:length(tickers)){
  cool <- Emp_test5(tickers[i],1000,5,1)[[1]]
  Z_ES_emp_array_5_250n[i,1] <- cool[2,1]
  Z_ES_emp_array_5_250n[i,2] <- cool[2,2]
  Z_2_emp_array_5_250n[i,1] <- cool[1,1]
  Z_2_emp_array_5_250n[i,2] <- cool[1,2]
  cat("Iteration", i, "\n")
}

end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)



Z_ES_emp_array_5_250n[,1,drop=FALSE]
Z_2_emp_array_5_250n[,1,drop=FALSE]
sum(Z_ES_emp_array_5_250n[,1]<0.025)/length(tickers)
sum(Z_2_emp_array_5_250n[,1]<0.025)/length(tickers)




