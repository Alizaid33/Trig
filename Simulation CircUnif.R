
library("Directional")
library("circular")
library("openxlsx")

GroupedR<-function(n,k,p){
  
  x<-rcircularuniform(n, control.circular=list())
  
  
  # The resultant length for raw data
  a<-trigonometric.moment(x, p , center = FALSE, control.circular = list())$rho
  
  # Grouping
  
  group <- seq(min(x) - 0.1, max(x) + 0.1, length = k+1)
  
  
  midpoints <- (group[-(k + 1)] + group[-1]) / 2
  
  # The length of actual group
  h<-(group[-1]-group[-(k + 1)])
  hB<-h[1]
  
  y <- cut(x, breaks = group, length = k+1)
  group <- matrix( c( group[1], rep(group[2:(k)], each = 2), group[k+1]), ncol = 2, byrow = TRUE)
  
  fi <- as.vector( table(y) )
  Si<-fi*sin(p*midpoints)
  Ci<-fi*cos(p*midpoints)
  
  # Frequency distribution of grouped data
  FrqDist<-cbind(group,midpoints, fi,Si,Ci)
  
  # uncorrected mean resultant length
  GrR<-sqrt((sum(Si)/n)^2+(sum(Ci)/n)^2)  
  
  
  # Subinterval length
  
  hPi=(2*pi)/(k)     # Group width based on data range =2 pi
  
  h=hB               # Group width based on actual data range
  
  
  # Corrected resultant length based on hPi
  CofRhPi<-((p*hPi)/(2*sin((p*hPi)/2)))*GrR
  
  # Corrected resultant length based on h
  CofR2<-((p*h)/(2*sin((p*h)/2)))*GrR
  
  Results<-cbind("Raw Data"=a,"Uncorrected"=GrR, "CofR2"=CofR2, "CofRhPi"=CofRhPi)
  
  return(list(Results = Results))
  
}


GroupedRSimu<-function(n,k, p){
  
  Simu<-1000
  
  Estimations<-matrix(0,nrow=Simu, ncol=4)
  
  for(i in 1:Simu){
    
    Estimations[i,]<-GroupedR(n=n,k=k,p)$Results
  }
  TrueMean=mean(Estimations[,1])
  # Calculate the bias for each column
  bias_uncorrected <- abs(mean(Estimations[,2] - Estimations[,1]))
  bias_cofR2 <- abs(mean(Estimations[,3] - Estimations[,1]))
  bias_cofRhPi <- abs(mean(Estimations[,4] - Estimations[,1]))
  
  # Calculate the Mean Squared Error (MSE) for each column
  mse_uncorrected <- mean(((Estimations[,2]) - Estimations[,1])^2)
  mse_cofR2 <- mean(((Estimations[,3]) - Estimations[,1])^2)
  mse_cofRhPi<- mean(((Estimations[,4]) - Estimations[,1])^2)
  
  Indicators<-cbind("TrueMean"=TrueMean,"bias_uncorrected"=bias_uncorrected,"bias_cofR2"=bias_cofR2,"bias_cofRhPi"=bias_cofRhPi,"mse_uncorrected"=mse_uncorrected,"mse_cofR2"=mse_cofR2,"mse_cofRhPi"=mse_cofRhPi)
  
  return(list(Indicators = Indicators))
}
set.seed(2024)

# Define parameter combinations for n, k, and kappa
n_values <- c(30, 50, 70, 100, 150, 200)
k_values <- c(5, 7, 10, 15)

p=2   # moment order

# Generate all combinations of n, k, kappa
param_combinations <- expand.grid(n = n_values, k = k_values)


# Initialize an empty list to store the results
results_list <- list()

simulation_results<-matrix(0,nrow=nrow(param_combinations), ncol=7)

# Run the simulation for each combination of n, k, and kappa
for (i in 1:nrow(param_combinations)) {
  n <- param_combinations$n[i]
  k <- param_combinations$k[i]
  
  
  # Run the simulations for this combination of parameters
  
  simulation_results[i,] <- GroupedRSimu(n = n, k = k,p)$Indicators
}

results_df <- cbind(param_combinations, simulation_results)
colnames(results_df) <- c("n","k","TrueMean","bias_uncorrected","bias_cofR2","bias_cofRhPi","mse_uncorrected","mse_cofR2","mse_cofRhPi")


# Write the results to an Excel file
write.xlsx(results_df, "SimulationUniformP=2.xlsx")


