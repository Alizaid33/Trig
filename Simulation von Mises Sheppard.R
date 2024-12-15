library("Directional")
library("circular")
library("openxlsx")
library(ggplot2)
library(tidyr)

GroupedVaR<-function(n,k,Mu, kappa,p){
  
  x <- rvonmises(n, Mu, kappa )
  
  
  # Grouping
  group <- seq(min(x) - 0.1, max(x) + 0.1, length = k+1)
  
  
  midpoints <- (group[-(k + 1)] + group[-1]) / 2
  
  # The length of actual group
  hI<-(group[-1]-group[-(k + 1)])
  h<-hI[1]  # Group width based on actual data range
  
  y <- cut(x, breaks = group, length = k+1)
  group <- matrix( c( group[1], rep(group[2:(k)], each = 2), group[k+1]), ncol = 2, byrow = TRUE)
  fi <- as.vector( table(y) )
  
  
  # GroupRp<-function(p,x,fi,midpoints,h){
  
  # The variance for raw data
  VRaw<-1-trigonometric.moment(x, p , center = FALSE, control.circular = list())$rho
  
  # p trignometric moments for grouped data 
  
  Si<-fi*sin(p*midpoints)
  Ci<-fi*cos(p*midpoints)
  
  # uncorrected mean resultant length
  GrRp<-sqrt((sum(Si)/n)^2+(sum(Ci)/n)^2)
  
  # uncorrected variance for grouped data
  VarGroup<-1-GrRp
  
  # Corrected variance circular Sheppard
  CGVar<-1-((p*h)/(2*sin((p*h)/2)))*GrRp
  
  # Corrected variance non -circular Sheppard
  SGrVar<-VarGroup-((h^2)/12)
  
  Variances<-cbind(VRaw, VarGroup, CGVar,SGrVar)
  
  
  colnames(Variances)<-c("Variance of Raw","Uncorrected Variance","Circular Correction","Non_Cir Correction")
  return(list(Variances= Variances))
  
}

GroupedVarSimu<-function(n,k,Mu, kappa,p){
  
  Simu<-1000
  
  Estimations<-matrix(0,nrow=Simu, ncol=4)
  
  for(i in 1:Simu){
    
    Estimations[i,]<-GroupedVaR(n=n,k=k,Mu=pi, kappa=kappa,p)$Variances
  }
  TrueMean1=mean(Estimations[,1])
  
  
  # Calculate the bias for the Variance of grouped data
  BiasUncoVar <- abs(mean(Estimations[,2] - TrueMean1))
  BiasCircCoVar <- abs(mean(Estimations[,3] - TrueMean1))
  BiasNonCircCoVar <- abs(mean(Estimations[,4] - TrueMean1))
  
  # Calculate the Mean Squared Error (MSE) for the Variance of grouped data
  
  MSEUncoVar <- mean(((Estimations[,2]) - TrueMean1)^2)
  MSECircCoVar <- mean(((Estimations[,3]) - TrueMean1)^2)
  MSENonCircCoVar <- mean(((Estimations[,4]) - TrueMean1)^2)
  
  
  
  Indicators<-cbind("BiasUncoVar"=BiasUncoVar,"BiasCircCoVar"=BiasCircCoVar,"BiasNonCircCoVar"=BiasNonCircCoVar,"MSEUncoVar"=MSEUncoVar,"MSECircCoVar"=MSECircCoVar,"MSENonCircCoVar"=MSENonCircCoVar)
  
  return(list(Indicators = Indicators))
}



set.seed(2024)

# Define parameter combinations for n, k, and kappa
n_values <- c(30, 50, 70, 100, 150, 200)
k_values <- c(5, 7, 10, 15)
kappa_values <- c(0.2, 0.5, 1, 2, 5, 10)



# Generate all combinations of n, k, kappa
param_combinations <- expand.grid(n = n_values, k = k_values, kappa = kappa_values)


# Initialize an empty list to store the results
results_list <- list()

simulation_results<-matrix(0,nrow=nrow(param_combinations), ncol=6)

# Run the simulation for each combination of n, k, and kappa
for (i in 1:nrow(param_combinations)) {
  n <- param_combinations$n[i]
  k <- param_combinations$k[i]
  kappa <- param_combinations$kappa[i]
  
  
  # Run the simulations for this combination of parameters
  
  simulation_results[i,] <- GroupedVarSimu(n = n, k = k, Mu = pi, kappa = kappa,p)$Indicators
}

results_df <- cbind(param_combinations, simulation_results)
colnames(results_df) <- c("n","k","kappa","BiasUncoVar","BiasCircCoVar","BiasNonCircCoVar","MSEUncoVar","MSECircCoVar","MSENonCircCoVar")

data<-as.data.frame(results_df)
# Filter the data for n = 100 and k = 7
filtered_data <- data[data$n == 100 & data$k == 7, ]

# Reshape the data into long format (melt the multiple GBiasNorm columns)
long_data <- filtered_data %>%
  gather(key = "Metric", value = "Value", BiasUncoVar,BiasCircCoVar,BiasNonCircCoVar)
#gather(key = "Metric", value = "Value", BiasUncoVar,BiasCircCoVar,BiasNonCircCoVar,MSEUncoVar,MSECircCoVar,MSENonCircCoVar)
# Plot the data using ggplot2
aa<-ggplot(long_data, aes(x = kappa, y = Value, color = Metric)) +
  # geom_point(size = 0) +               # Scatter plot with points
  geom_line(aes(group = Metric), size = 1) +  # Add lines connecting points
  #labs(title = "Bias of estimated variances for grouped data, n= 100 and K=7",
  labs(x = "Kappa",
       y = "Bias") +
  theme_minimal() +
  theme(legend.title = element_blank())

ggsave("myplotVariance.pdf",aa)

# Filter the data for kappa = 2 and k = 7
filtered_data <- data[data$kappa == 2 & data$k == 7, ]

# Reshape the data into long format (melt the multiple GBiasNorm columns)
long_data <- filtered_data %>%
  gather(key = "Metric", value = "Value", BiasUncoVar,BiasCircCoVar,BiasNonCircCoVar)
#gather(key = "Metric", value = "Value", BiasUncoVar,BiasCircCoVar,BiasNonCircCoVar,MSEUncoVar,MSECircCoVar,MSENonCircCoVar)
# Plot the data using ggplot2
aa<-ggplot(long_data, aes(x = n, y = Value, color = Metric)) +
  # geom_point(size = 0) +               # Scatter plot with points
  geom_line(aes(group = Metric), size = 1) +  # Add lines connecting points
  #labs(title = "Bias of estimated variances for grouped data, kappa= 5 and K=7",
  labs(x = "n",
       y = "Bias") +
  theme_minimal() +
  theme(legend.title = element_blank())

ggsave("VarianceSampleSize.pdf",aa)



# Write the results to an Excel file
write.xlsx(results_df, "SimuVMVaeiances.xlsx")





