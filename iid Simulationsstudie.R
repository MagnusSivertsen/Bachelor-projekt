
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

T <- 250
Alpha <- 0.025
Beta <- 0.01

# Funktioner til at udregne Var og ES for t fordeling og normaliseret t fordeling

Var_st <- function(df){-qt(Alpha,df)}

ES_st <- function(df){(dt((-1)*Var_st(df),df)*(df+(-Var_st(df))^2)/(Alpha*(df-1)))}

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
Var_std <- function(df){-qstd(Alpha,0,1,df)}

# Funktioner til at udregne de forskellige tests.

Z_1 <- function(data,Var,ES){
  I_t <- as.numeric(data < -Var)
  (sum((data*I_t)/ES)/sum(I_t))+1
}
Z_2 <- function(data,Var,ES){
  I_t <- as.numeric(data < -Var)
  sum((data*I_t/(T*Alpha*ES)))+1
}
Z_3 <- function(data,df1){
  TOP <- (-1/floor(T*Alpha))*sum(sort(data)[1:floor(T*Alpha)])
  BOTTOM <- -T/floor(T*Alpha)*integrate(function(p) pbeta(1-p,T-floor(T*Alpha),floor(T*Alpha))*qt(p,df=df1),lower=0,upper=1)$value
  (-TOP/BOTTOM)+1
}
Z_3_2 <- function(data,df1,BOTTOM){
  TOP <- (-1/floor(T*Alpha))*sum(sort(data)[1:floor(T*Alpha)])
  (-TOP/BOTTOM)+1
}
BOTTOMFNC <- function(df1){
  -T/floor(T*Alpha)*integrate(function(p) pbeta(1-p,T-floor(T*Alpha),floor(T*Alpha))*qt(p,df=df1),lower=0,upper=1)$value
}
BOTTOMFNC2 <- function(df1){
  -T/floor(T*Alpha)*integrate(function(p) pbeta(1-p,T-floor(T*Alpha),floor(T*Alpha))*qstd(p,0,1,nu=df1),lower=0,upper=1)$value
}
Var_ex <- function(data,Var)
  -sum(data < -Var)
Z_ES <- function(data,Var,ES){
  sum((ES-(Var-1/Alpha*pmin(data+Var,0)))/ES)/T
}

# funktion til at udregne power  

PowerF <- function(matrix,Siglevel){
  quantile1 <- quantile(matrix[,1],Siglevel)
  powers <- numeric(ncol(matrix)-1)
  for (i in 2:ncol(matrix)){
    powers[i-1] <- ecdf(matrix[,i])(quantile1)
  }
  return(powers)
}



# Sektion xx udregning normal 

set.seed(1964)
M <- 100000 
df_1 <-c(100,10,5,3)
Z2_array_1 = matrix(0,M,4) 
Z3_array_1 = matrix(0,M,4)
Var_array_1 = matrix(0,M,4)
Z_ES_array_1 = matrix(0,M,4)
bot1 <- BOTTOMFNC(100)
for (i in 1:4){
  for (k in 1:M){
    data1 <- rt(T,df=df_1[i])
    Z2_array_1[k,i] <- Z_2(data1,Var_st(df_1[1]),ES_st(df_1[1]))
    Z3_array_1[k,i] <- Z_3_2(data1,100,bot1)
    Var_array_1[k,i] <- Var_ex(data1,-qt(Beta,df_1[1]))
    Z_ES_array_1[k,i] <- Z_ES(data1,Var_st(df_1[1]),ES_st(df_1[1]))
  }
}

set.seed(1964)
df_2 <-c(10,5,3)
Z2_array_2 <- matrix(0,M,3) 
Z3_array_2 <- matrix(0,M,3)
Var_array_2 = matrix(0,M,3)
Z_ES_array_2 = matrix(0,M,3)
bot2 <- BOTTOMFNC(10)
for (i in 1:3){
  for (k in 1:M){
    data2 <- rt(T,df=df_2[i])
    Z2_array_2[k,i] <- Z_2(data2,Var_st(df_2[1]),ES_st(df_2[1]))
    Z3_array_2[k,i] <- Z_3_2(data2,10,bot2)
    Var_array_2[k,i] <- Var_ex(data2,-qt(Beta,df_2[1]))
    Z_ES_array_2[k,i] <- Z_ES(data2,Var_st(df_2[1]),ES_st(df_2[1]))
  }
}


# Sektion xx uregning normaliseret t

set.seed(1964)
M <- 100000 
df_1 <-c(100,10,5,3)
Z2_array_1.2 = matrix(0,M,4) 
Z3_array_1.2 = matrix(0,M,4)
Var_array_1.2 = matrix(0,M,4)
Z_ES_array_1.2 = matrix(0,M,4)
bot1.2 <- BOTTOMFNC2(100)
for (i in 1:4){
  for (k in 1:M){
    data1 <- rstd(T,nu=df_1[i])
    Z2_array_1.2[k,i] <- Z_2(data1,Var_std(df_1[1]),ES_t_standard(Alpha,df_1[1]))
    Z3_array_1.2[k,i] <- Z_3_2(data1,100,bot1.2)
    Var_array_1.2[k,i] <- Var_ex(data1,-qstd(Beta,0,1,df_1[1]))
    Z_ES_array_1.2[k,i] <- Z_ES(data1,Var_std(df_1[1]),ES_t_standard(Alpha,df_1[1]))
  }
}

set.seed(1964)
df_2 <-c(10,5,3)
Z2_array_2.2 <- matrix(0,M,3) 
Z3_array_2.2 <- matrix(0,M,3)
Var_array_2.2 = matrix(0,M,3)
Z_ES_array_2.2 = matrix(0,M,3)
bot2.2 <- BOTTOMFNC2(10)
for (i in 1:3){
  for (k in 1:M){
    data2 <- rstd(T,nu=df_2[i])
    Z2_array_2.2[k,i] <- Z_2(data2,Var_std(df_2[1]),ES_t_standard(Alpha,df_2[1]))
    Z3_array_2.2[k,i] <- Z_3_2(data2,10,bot2.2)
    Var_array_2.2[k,i] <- Var_ex(data2,-qstd(Beta,0,1,df_2[1]))
    Z_ES_array_2.2[k,i] <- Z_ES(data2,Var_std(df_2[1]),ES_t_standard(Alpha,df_2[1]))
  }
}

# Sektion xx udregning fix var 


M <- 100000
Z2_array_5 <- matrix(0,M,4) 
Z3_array_5 <- matrix(0,M,4)
Z2_array_6 <- matrix(0,M,3) 
Z3_array_6 <- matrix(0,M,3)
Z1_array_1 <- matrix(0,M,4)
Z1_array_2 <- matrix(0,M,3)
Var_array_5 <- matrix(0,M,4)
Var_array_6 <- matrix(0,M,3)
Z_ES_array_5 <- matrix(0,M,4)
Z_ES_array_6 <- matrix(0,M,3)

for (i in 1:4){
  for (k in 1:M){
    data1 <- rt(T,df=df_1[i]) + qt(0.025,df_1[1])-qt(0.025,df_1[i])
    Z2_array_5[k,i] <- Z_2(data1,Var_st(df_1[1]),ES_st(df_1[1]))
    Z3_array_5[k,i] <- Z_3_2(data1,100,bot1)
    Z1_array_1[k,i] <- Z_1(data1,Var_st(df_1[1]),ES_st(df_1[1]))
    Var_array_5[k,i] <- Var_ex(data1,-qt(Beta,100))
    Z_ES_array_5[k,i]<- Z_ES(data1,Var_st(df_1[1]),ES_st(df_1[1]))
  }
}


for (i in 1:3){
  for (k in 1:M){
    data2 <- rt(T,df=df_2[i]) + qt(0.025,df_2[1])-qt(0.025,df_2[i])
    Z2_array_6[k,i] <- Z_2(data2,Var_st(df_2[1]),ES_st(df_2[1]))
    Z3_array_6[k,i] <- Z_3_2(data2,10,bot2)
    Z1_array_2[k,i] <- Z_1(data2,Var_st(df_2[1]),ES_st(df_2[1]))
    Var_array_6[k,i] <- Var_ex(data2,-qt(Beta,10))
    Z_ES_array_6[k,i]<- Z_ES(data2,Var_st(df_2[1]),ES_st(df_2[1]))
  }
}

# Sektion xx udregning fix var med normaliseret t.
set.seed(1964)
M <- 100000 
df_1 <-c(100,10,5,3)
Z2_array_5.2 = matrix(0,M,4) 
Z3_array_5.2 = matrix(0,M,4)
Var_array_5.2 = matrix(0,M,4)
Z_ES_array_5.2 = matrix(0,M,4)
Z1_array_1.2 <- matrix(0,M,4)
Z1_array_2.2 <- matrix(0,M,3)
bot1.2 <- BOTTOMFNC2(100)
for (i in 1:4){
  for (k in 1:M){
    data1 <- rstd(T,nu=df_1[i]) + qstd(0.025,0,1,df_1[1])-qstd(0.025,0,1,df_1[i])
    
    Z2_array_5.2[k,i] <- Z_2(data1,Var_std(df_1[1]),ES_t_standard(Alpha,df_1[1]))
    Z3_array_5.2[k,i] <- Z_3_2(data1,100,bot1.2)
    Var_array_5.2[k,i] <- Var_ex(data1,-qstd(Beta,0,1,df_1[1]))
    Z_ES_array_5.2[k,i] <- Z_ES(data1,Var_std(df_1[1]),ES_t_standard(Alpha,df_1[1]))
    Z1_array_1.2[k,i] <- Z_1(data1,Var_std(df_1[1]),ES_t_standard(Alpha,df_1[1]))
  }
}

set.seed(1964)
df_2 <-c(10,5,3)
Z2_array_6.2 <- matrix(0,M,3) 
Z3_array_6.2 <- matrix(0,M,3)
Var_array_6.2 = matrix(0,M,3)
Z_ES_array_6.2 = matrix(0,M,3)
bot2.2 <- BOTTOMFNC2(10)
for (i in 1:3){
  for (k in 1:M){
    data2 <- rstd(T,nu=df_2[i]) + qstd(0.025,0,1,df_2[1])-qstd(0.025,0,1,df_2[i])
    Z2_array_6.2[k,i] <- Z_2(data2,Var_std(df_2[1]),ES_t_standard(Alpha,df_2[1]))
    Z3_array_6.2[k,i] <- Z_3_2(data2,10,bot2.2)
    Var_array_6.2[k,i] <- Var_ex(data2,-qstd(Beta,0,1,df_2[1]))
    Z_ES_array_6.2[k,i] <- Z_ES(data2,Var_std(df_2[1]),ES_t_standard(Alpha,df_2[1]))
    Z1_array_2.2[k,i] <- Z_1(data2,Var_std(df_2[1]),ES_t_standard(Alpha,df_2[1]))
  }
}

# Udregning af power for sektion xx
print("Z_2")
PowerF(Z2_array_1,0.041)
PowerF(Z2_array_1,0.104)
PowerF(Z2_array_2,0.04)
PowerF(Z2_array_2,0.106)
print("Z_3")
PowerF(Z3_array_1,0.041)
PowerF(Z3_array_1,0.104)
PowerF(Z3_array_2,0.04)
PowerF(Z3_array_2,0.106)
print("Z_ES")
PowerF(Z_ES_array_1,0.041)
PowerF(Z_ES_array_1,0.104)
PowerF(Z_ES_array_2,0.04)
PowerF(Z_ES_array_2,0.106)
print("VaR")
PowerF(Var_array_1,0.041)
PowerF(Var_array_1,0.104)
PowerF(Var_array_2,0.04)
PowerF(Var_array_2,0.106)

# udregning af power for sektion xx standard t

print("Z_2")
PowerF(Z2_array_1.2,0.044)
PowerF(Z2_array_1.2,0.112)
PowerF(Z2_array_2.2,0.044)
PowerF(Z2_array_2.2,0.11)
print("Z_3")
PowerF(Z3_array_1.2,0.044)
PowerF(Z3_array_1.2,0.112)
PowerF(Z3_array_2.2,0.044)
PowerF(Z3_array_2.2,0.11)
print("Z_ES")
PowerF(Z_ES_array_1.2,0.044)
PowerF(Z_ES_array_1.2,0.112)
PowerF(Z_ES_array_2.2,0.044)
PowerF(Z_ES_array_2.2,0.11)
print("VaR")
PowerF(Var_array_1.2,0.044)
PowerF(Var_array_1.2,0.112)
PowerF(Var_array_2.2,0.044)
PowerF(Var_array_2.2,0.11)

# udregning af power for sektion xx fix var 

print("Z_1")
PowerF(na.omit(Z1_array_1),0.041)
PowerF(na.omit(Z1_array_1),0.111)
PowerF(na.omit(Z1_array_2),0.042)
PowerF(na.omit(Z1_array_2),0.114)
print("Z_2")
PowerF(Z2_array_5,0.041)
PowerF(Z2_array_5,0.111)
PowerF(Z2_array_6,0.042)
PowerF(Z2_array_6,0.114)
print("Z_3")
PowerF(Z3_array_5,0.041)
PowerF(Z3_array_5,0.111)
PowerF(Z3_array_6,0.042)
PowerF(Z3_array_6,0.114)
print("Z_ES")
PowerF(Z_ES_array_5,0.041)
PowerF(Z_ES_array_5,0.111)
PowerF(Z_ES_array_6,0.042)
PowerF(Z_ES_array_6,0.114)
print("VaR")
PowerF(Var_array_5,0.041)
PowerF(Var_array_5,0.111)
PowerF(Var_array_6,0.042)
PowerF(Var_array_6,0.114)

# udregning af power for sektion xx fixed var standard t

print("Z_1")
PowerF(Z1_array_1.2,0.041)
PowerF(Z1_array_1.2,0.111)
PowerF(Z1_array_2.2,0.042)
PowerF(Z1_array_2.2,0.114)
print("Z_2")
PowerF(Z2_array_5.2,0.041)
PowerF(Z2_array_5.2,0.111)
PowerF(Z2_array_6.2,0.042)
PowerF(Z2_array_6.2,0.114)
print("Z_3")
PowerF(Z3_array_5.2,0.041)
PowerF(Z3_array_5.2,0.111)
PowerF(Z3_array_6.2,0.042)
PowerF(Z3_array_6.2,0.114)
print("Z_ES")
PowerF(Z_ES_array_5.2,0.041)
PowerF(Z_ES_array_5.2,0.111)
PowerF(Z_ES_array_6.2,0.042)
PowerF(Z_ES_array_6.2,0.114)
print("VaR")
PowerF(Var_array_5.2,0.041)
PowerF(Var_array_5.2,0.111)
PowerF(Var_array_6.2,0.042)
PowerF(Var_array_6.2,0.114)

# udregning af mean for Z_2 og Z_ES i sektion xx, og kritiske værdier

EStt <- ES_st(5)
sigma_values <- seq(0.5, 1.5, length.out = b)
mu_values <- -EStt + sigma_values * EStt
DiffVar <- -(mu_values + sigma_values*qt(0.025,df=5))

Z_ES_array_fix_ES1.2 <- matrix(0,M,b) 
Z_2_array_fix_ES1.2 <- matrix(0,M,b)

for (i in 1:b){
  for (k in 1:M){
    data1 <- rt(T,5)
    Z_2_array_fix_ES1.2[k,i] <- Z_2(data1,DiffVar[i],EStt)
    Z_ES_array_fix_ES1.2[k,i] <- Z_ES(data1,DiffVar[i],EStt)
  }
}

mean_Z_ES1.2 <- matrix(0,b,2)
for (i in 1:b){
  mean_Z_ES1.2[i,1] <- (DiffVar[i]+qt(0.025,5))/DiffVar[i]
  mean_Z_ES1.2[i,2] <- mean(Z_ES_array_fix_ES1.2[,i]) 
}  

mean_Z_21.2 <- matrix(0,b,2)
for (i in 1:b){
  mean_Z_21.2[i,1] <- (DiffVar[i]+qt(0.025,5))/DiffVar[i]
  mean_Z_21.2[i,2] <- mean(Z_2_array_fix_ES1.2[,i]) 
}  


crit_val_Z_ES1.2 <- quantile(Z_ES_array_fix_ES1.2[,((b-1)/2)+1],0.05)
crit_val_Z_21.2 <- quantile(Z_2_array_fix_ES1.2[,((b-1)/2)+1],0.05)


crit_val_Z_ES1.2
crit_val_Z_21.2
# plot xx i sektion xx
df_Z2 <- data.frame(x = mean_Z_21.2[,1], y = mean_Z_21.2[,2], test = "Z2")
df_ES <- data.frame(x = mean_Z_ES1.2[,1], y = mean_Z_ES1.2[,2], test = "ZES")


df_all <- rbind(df_Z2, df_ES)


plot2 <- ggplot(df_all, aes(x = x, y = y, color = test)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = -0.7434122, linetype = "dotted", color = "#F8766D", size = 1.2) +
  geom_hline(yintercept = -0.3215843, linetype = "dotted", color = "#00BFC4", size = 1.2) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 1) +
  coord_cartesian(xlim = c(-0.3, 0.2), ylim = c(-1, 0.3)) +
  labs(title = "Sensitiviteten af Z_ES og Z_2 når ES er estimeret korrekt og VaR er forkert",
       x = "(v-VaR)/v",y="") +
  theme_minimal()
plot2

# udregning af Z_ES og Z_2 ved e = 0.8ES
set.seed(1964)
M <- 10000
b <- 41
ESnormal <- ESnorm(1-0.025,0,1)
sigma_values1 <- seq(-1, 1, length.out = b)
mu_values1 <- -ESnormal*0.8 + sigma_values1 * ESnormal

DiffVar1 <- -(mu_values1 + sigma_values1*qnorm(0.025,0,1))

Z_ES_array_fix_ES2 <- matrix(0,M,b)
Z_2_array_fix_ES2 <- matrix(0,M,b)

for (i in 1:b){
  for (k in 1:M){
    data1 <- rnorm(T,0,1)
    Z_2_array_fix_ES2[k,i] <- Z_2(data1,DiffVar1[i],ESnormal*0.8)
    Z_ES_array_fix_ES2[k,i] <- Z_ES(data1,DiffVar1[i],ESnormal*0.8)
  }
}
Null_array <- matrix(0,M,2)
for (i in 1:M){
  data1 <- rnorm(T,0,1)
  Null_array[i,2] <- Z_ES(data1,-qnorm(0.025,0,1),ESnormal)
  Null_array[i,1] <- Z_2(data1,-qnorm(0.025,0,1),ESnormal)
}

mean_Z_ES2 <- matrix(0,b,2)
for (i in 1:b){
  mean_Z_ES2[i,1] <- (DiffVar1[i]+qnorm(0.025,0,1))/DiffVar1[i]
  mean_Z_ES2[i,2] <- mean(Z_ES_array_fix_ES2[,i]) 
}  

mean_Z_22 <- matrix(0,b,2)
for (i in 1:b){
  mean_Z_22[i,1] <- (DiffVar1[i]+qnorm(0.025,0,1))/DiffVar1[i]
  mean_Z_22[i,2] <- mean(Z_2_array_fix_ES2[,i]) 
} 


crit_val_Z_ES2 <- quantile(Null_array[,2],0.05)
crit_val_Z_22 <- quantile(Null_array[,1],0.05)



crit_val_Z_22
crit_val_Z_ES2

# Plot af dette resultat

df_Z22 <- data.frame(x = mean_Z_22[,1], y = mean_Z_22[,2], test = "Z2")
df_ES2 <- data.frame(x = mean_Z_ES2[,1], y = mean_Z_ES2[,2], test = "ZES")
df_all2 <- rbind(df_Z22, df_ES2)


plot1 <- ggplot(df_all2, aes(x = x, y = y, color = test)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = -0.690929, linetype = "dotted", color = "#F8766D", size = 1.2) +
  geom_hline(yintercept = -0.1623956 , linetype = "dotted", color = "#00BFC4", size = 1.2) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 1) +  
  coord_cartesian(xlim = c(-0.3, 0.05), ylim = c(-1.6, 0.05)) +
  labs(title = "Sensitiviteten af Z_ES og Z_2 når e = 0.8ES og VaR er forkert" ,
       x = "(v-VaR)/v",y="") +
  theme_minimal()

# Udregning af kritiske værdier for forskellig location og df for Z_ES

df_12 <- c(100,10,5,3)
critmatrix0 <- matrix(0,M,4)
critmatrix1 <- matrix(0,M,4)
critmatrix_1 <- matrix(0,M,4)

for (i in 1:4){
  for (k in 1:M){
    data1 <- rt(250,df_12[i])
    data2 <- rt(250,df_12[i]) + 1 
    data3 <- rt(250,df_12[i]) - 1
    critmatrix0[k,i] <- Z_ES(data1,Var_st(df_12[i]),ES_st(df_12[i]))
    critmatrix1[k,i] <- Z_ES(data2,Var_st(df_12[i])+1,ES_st(df_12[i])+1)
    critmatrix_1[k,i] <- Z_ES(data3,Var_st(df_12[i])-1,ES_st(df_12[i])-1)
    
  }
}




c(quantile(critmatrix0[,1],0.05),quantile(critmatrix0[,2],0.05),quantile(critmatrix0[,3],0.05),quantile(critmatrix0[,4],0.05))
c(quantile(critmatrix1[,1],0.05),quantile(critmatrix1[,2],0.05),quantile(critmatrix1[,3],0.05),quantile(critmatrix1[,4],0.05))
c(quantile(critmatrix_1[,1],0.05),quantile(critmatrix_1[,2],0.05),quantile(critmatrix_1[,3],0.05),quantile(critmatrix_1[,4],0.05))

c(quantile(critmatrix0[,1],0.0001),quantile(critmatrix0[,2],0.0001),quantile(critmatrix0[,3],0.0001),quantile(critmatrix0[,4],0.0001))
c(quantile(critmatrix1[,1],0.0001),quantile(critmatrix1[,2],0.0001),quantile(critmatrix1[,3],0.0001),quantile(critmatrix1[,4],0.0001))
c(quantile(critmatrix_1[,1],0.0001),quantile(critmatrix_1[,2],0.0001),quantile(critmatrix_1[,3],0.0001),quantile(critmatrix_1[,4],0.0001))

# Udregning af kritiske værdier for forskellig location og df for Z_2


df_12 <- c(100,10,5,3)
critmatrix00 <- matrix(0,M,4)
critmatrix11 <- matrix(0,M,4)
critmatrix_11 <- matrix(0,M,4)

for (i in 1:4){
  for (k in 1:M){
    data1 <- rt(250,df_12[i])
    data2 <- rt(250,df_12[i]) + 1 
    data3 <- rt(250,df_12[i]) - 1
    critmatrix00[k,i] <- Z_2(data1,Var_st(df_12[i]),ES_st(df_12[i]))
    critmatrix11[k,i] <- Z_2(data2,Var_st(df_12[i])-1,ES_st(df_12[i])-1)
    critmatrix_11[k,i] <- Z_2(data3,Var_st(df_12[i])+1,ES_st(df_12[i])+1)
    
  }
}

c(quantile(critmatrix00[,1],0.05),quantile(critmatrix00[,2],0.05),quantile(critmatrix00[,3],0.05),quantile(critmatrix00[,4],0.05))
c(quantile(critmatrix11[,1],0.05),quantile(critmatrix11[,2],0.05),quantile(critmatrix11[,3],0.05),quantile(critmatrix11[,4],0.05))
c(quantile(critmatrix_11[,1],0.05),quantile(critmatrix_11[,2],0.05),quantile(critmatrix_11,0.05),quantile(critmatrix_11[,4],0.05))

c(quantile(critmatrix00[,1],0.0001),quantile(critmatrix00[,2],0.0001),quantile(critmatrix00[,3],0.0001),quantile(critmatrix00[,4],0.0001))
c(quantile(critmatrix11[,1],0.0001),quantile(critmatrix11[,2],0.0001),quantile(critmatrix11[,3],0.0001),quantile(critmatrix11[,4],0.0001))
c(quantile(critmatrix_11[,1],0.0001),quantile(critmatrix_11[,2],0.0001),quantile(critmatrix_11,0.0001),quantile(critmatrix_11[,4],0.0001))

