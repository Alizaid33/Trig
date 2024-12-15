
library("Directional")
library("circular")
library("openxlsx")
library(ggplot2)
library(tidyr)
library(dplyr)

GroupedR<-function(n,k,Mu, kappa){
  
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
  print(fi)
  
  GroupRp<-function(p,x,fi,midpoints){
    
    # The resultant length for raw data
    Rp<-trigonometric.moment(x, p , center = FALSE, control.circular = list())$rho
    
    # p trignometric moments for grouped data 
    
    Si<-fi*sin(p*midpoints)
    Ci<-fi*cos(p*midpoints)
    
    # uncorrected mean resultant length
    GrRp<-sqrt((sum(Si)/n)^2+(sum(Ci)/n)^2)
    
    # Corrected resultant length based on h
    CGrRp<-((p*h)/(2*sin((p*h)/2)))*GrRp
    
    Estimates<-cbind(Rp, GrRp,CGrRp)
    
    list(Estimates=Estimates)
  }
  
  
  
  GroupR1<- GroupRp(1,x,fi,midpoints)$Estimates
  GroupR2<- GroupRp(2,x,fi,midpoints)$Estimates
  GroupR3<- GroupRp(3,x,fi,midpoints)$Estimates
  GroupR4<- GroupRp(4,x,fi,midpoints)$Estimates
  
  
  Norm1<- GroupR1
  Norm2<-(GroupR2+1)/2
  Norm3<-(GroupR3+3*GroupR1)/4
  Norm4<-(GroupR4+4*GroupR2+3)/8
  
  NORMS<-cbind(Norm1,Norm2,Norm4,Norm4)
  # G_Norm1= Grouped Norm 1,  CG_Norm1= Corrected Grouped Norm 1
  
  colnames(NORMS)<-c("Row_Norm1","G_Norm1","CG_Norm1","Row_Norm2","G_Norm2","CG_Norm2","Row_Norm3","G_Norm3","CG_Norm3","Row_Norm4","G_Norm4","CG_Norm4")
  return(list( NORMS=NORMS))
  
}


GroupedRSimu<-function(n,k,Mu, kappa){
  
  Simu<-1000
  
  Estimations<-matrix(0,nrow=Simu, ncol=12)
  
  for(i in 1:Simu){
    
    Estimations[i,]<-GroupedR(n=n,k=k,Mu=pi, kappa=kappa)$NORMS
  }
  TrueMean1=mean(Estimations[,1])
  TrueMean2=mean(Estimations[,4])
  TrueMean3=mean(Estimations[,7])
  TrueMean4=mean(Estimations[,10])
  
  # Calculate the bias for the Norm of grouped data
  GBiasNorm1 <- abs(mean(Estimations[,2] - TrueMean1))
  CGBiasNorm1 <- abs(mean(Estimations[,3] - TrueMean1))
  
  GBiasNorm2 <- abs(mean(Estimations[,5] - TrueMean2))
  CGBiasNorm2 <- abs(mean(Estimations[,6] - TrueMean2))
  
  GBiasNorm3 <- abs(mean(Estimations[,8] - TrueMean3))
  CGBiasNorm3 <- abs(mean(Estimations[,9] - TrueMean3))
  
  GBiasNorm4 <- abs(mean(Estimations[,11] - TrueMean4))
  CGBiasNorm4 <- abs(mean(Estimations[,12] - TrueMean4))
  
  
  # Calculate the Mean Squared Error (MSE) for the Norm of grouped data
  
  GMSENorm1 <- mean(((Estimations[,2]) - TrueMean1)^2)
  CGMSENorm1 <- mean(((Estimations[,3]) - TrueMean1)^2)
  
  GMSENorm2 <- mean(((Estimations[,5]) - TrueMean2)^2)
  CGMSENorm2 <- mean(((Estimations[,6]) - TrueMean2)^2)
  
  GMSENorm3 <- mean(((Estimations[,8]) - TrueMean3)^2)
  CGMSENorm3 <- mean(((Estimations[,9]) - TrueMean3)^2)
  
  GMSENorm4 <- mean(((Estimations[,11]) - TrueMean4)^2)
  CGMSENorm4 <- mean(((Estimations[,12]) - TrueMean4)^2)
  
  Indicators<-cbind("GBiasNorm1"=GBiasNorm1,"GBiasNorm2"=GBiasNorm2,"GBiasNorm3"=GBiasNorm3,"GBiasNorm4"=GBiasNorm4,"CGBiasNorm1"=CGBiasNorm1,"CGBiasNorm2"=CGBiasNorm2,"CGBiasNorm3"=CGBiasNorm3,"CGBiasNorm4"=CGBiasNorm4,"GMSENorm1"=GMSENorm1,"GMSENorm2"=GMSENorm2,"GMSENorm3"=GMSENorm3,"GMSENorm4"=GMSENorm4,"CGMSENorm1"=CGMSENorm1,"CGMSENorm2"=CGMSENorm2,"CGMSENorm3"=CGMSENorm3,"CGMSENorm4"=CGMSENorm4)
  
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

simulation_results<-matrix(0,nrow=nrow(param_combinations), ncol=16)

# Run the simulation for each combination of n, k, and kappa
for (i in 1:nrow(param_combinations)) {
  n <- param_combinations$n[i]
  k <- param_combinations$k[i]
  kappa <- param_combinations$kappa[i]
  
  
  # Run the simulations for this combination of parameters
  
  simulation_results[i,] <- GroupedRSimu(n = n, k = k, Mu = pi, kappa = kappa)$Indicators
}

results_df <- cbind(param_combinations, simulation_results)
colnames(results_df) <- c("n","k","kappa","GBiasNorm1","GBiasNorm2","GBiasNorm3","GBiasNorm4","CGBiasNorm1","CGBiasNorm2","CGBiasNorm3","CGBiasNorm4","GMSENorm1","GMSENorm2","GMSENorm3","GMSENorm4","CGMSENorm1","CGMSENorm2","CGMSENorm3","CGMSENorm4")

data<-as.data.frame(results_df)
# Filter the data for n = 100 and k = 5
filtered_data <- data[data$n == 100 & data$k == 7, ]

# Reshape the data into long format (melt the multiple GBiasNorm columns)
long_data <- filtered_data %>%
  gather(key = "Metric", value = "Value", 
         GBiasNorm1, GBiasNorm2, GBiasNorm3, GBiasNorm4,
         CGBiasNorm1, CGBiasNorm2, CGBiasNorm3, CGBiasNorm4)

# Create a new column to group the metrics
long_data$MetricGroup <- case_when(
  long_data$Metric %in% c("GBiasNorm1", "CGBiasNorm1") ~ "Norm 1",
  long_data$Metric %in% c("GBiasNorm2", "CGBiasNorm2") ~ "Norm 2",
  long_data$Metric %in% c("GBiasNorm3", "CGBiasNorm3") ~ "Norm 3",
  long_data$Metric %in% c("GBiasNorm4", "CGBiasNorm4") ~ "Norm 4",
  TRUE ~ "Other"
)

# Create the plot using ggplot2 with faceting by 'MetricGroup'
BBB <- ggplot(long_data, aes(x = kappa, y = Value, color = Metric)) +
  geom_line(aes(group = Metric), size = 1) +  # Add lines connecting points
  facet_wrap(~ MetricGroup, scales = "free", ncol = 2) +  # Facet into 2 columns
 # labs(title = "Scatter plot of GBiasNorm and CGBiasNorm values",
  labs(x = "Kappa",
       y = "Bias Norm Values") +
  theme_minimal() +
  theme(legend.title = element_blank())

# Display the plot
print(BBB)

ggsave("Normplot.pdf",BBB)

# Write the results to an Excel file
write.xlsx(results_df, "SimuVMNORMS.xlsx")


