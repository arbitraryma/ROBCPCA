library(MASS)
library(clue)
library(fungible)
library(reservoirnet)
japanese_vowels <- generate_data(
  dataset = "japanese_vowels",
  one_hot_encode = TRUE,
  repeat_targets = FALSE,
  reload = FALSE,
  n_timesteps,
  tau = 17,
  a = 0.2,
  b = 0.1,
  n = 10,
  x0 = 1.2,
  h = 1
)

adjacjv <- rep(0,100)
adjacrobjv <- rep(0,100)
adjacroblagjv <- rep(0,100)
adjacroblag2jv <- rep(0,100)
adjacroblag12jv <- rep(0,100)
accuracyjv <- rep(0,100)
accuracyrobjv <- rep(0,100)
accuracyroblagjv <- rep(0,100)
accuracyroblag2jv <- rep(0,100)
accuracyroblag12jv <- rep(0,100)

true_labels <- rep(0,length(japanese_vowels[["japanese_vowels"]][["Y_test"]]))
for (i in 1:length(japanese_vowels[["japanese_vowels"]][["Y_test"]])) {
  true_labels[i] <- which(japanese_vowels[["japanese_vowels"]][["Y_test"]][[i]] == 1)
}



for (i in 1:100) {
  

b <- cpca(japanese_vowels[["japanese_vowels"]][["X_test"]],9)
adjacjv[i] <- adjustedRandIndex(b,true_labels)
accuracyjv[i] <- true_accuracy(b,true_labels)
 


a <- robcpca(japanese_vowels[["japanese_vowels"]][["X_test"]],9)
adjacrobjv[i] <- adjustedRandIndex(a,true_labels)
accuracyrobjv[i] <- true_accuracy(a,true_labels)
 

d <- robcpca_lag1(japanese_vowels[["japanese_vowels"]][["X_test"]],9)
adjacroblagjv[i] <- adjustedRandIndex(d,true_labels)
accuracyroblagjv[i] <- true_accuracy(d,true_labels)

e <- robcpca_lag_2(japanese_vowels[["japanese_vowels"]][["X_test"]],9)
adjacroblag2jv[i] <- adjustedRandIndex(e,true_labels)
accuracyroblag2jv[i] <- true_accuracy(e,true_labels)

f <- robcpca_lag_12(japanese_vowels[["japanese_vowels"]][["X_test"]],9)
adjacroblag12jv[i] <- adjustedRandIndex(f,true_labels)
accuracyroblag12jv[i] <- true_accuracy(f,true_labels)
}

 
mean(accuracyjv)
# 0.3321892

mean(accuracyrobjv)
# 0.4412703
 
mean(accuracyroblagjv)
#  0.3662162
mean(accuracyroblag2jv)
#  0.3831892
mean(accuracyroblag12jv)
# 0.4854595




round(quantile(accuracyjv),2)
# 0%       25%       50%       75%      100% 
# 0.24 0.30 0.34 0.36 0.43 

round(quantile(accuracyrobjv),2)
# 0%       25%       50%       75%      100% 
# 0.18 0.40 0.45 0.48 0.57 


round(quantile(accuracyroblagjv),2)
# 0.26 0.34 0.37 0.39 0.45 


round(quantile(accuracyroblag2jv),2)
# 0.27 0.35 0.39 0.41 0.47 

round(quantile(accuracyroblag12jv),2)
# 0.36 0.45 0.48 0.52 0.63 



sd(accuracyjv) # 0.04
sd(accuracyrobjv) # 0.06
sd(accuracyroblagjv) #  0.04
sd(accuracyroblag2jv) # 0.04
sd(accuracyroblag12jv) # 0.05

t.test(accuracyjv,accuracyroblagjv,paired = TRUE) #  < 2.2e-16


# plot 




data <- data.frame(
  values = c(data$accuracyroblag12jv, accuracyroblag2jv, accuracyroblagjv, accuracyrobjv, accuracyjv),
  group = factor(rep(c("ROBCPCA012", "ROBCPCA02", "ROBCPCA01","ROBCPCA0", "CPCA"), each = length(accuracyjv)))
)
# Set the desired order of the groups
data$group <- factor(data$group, levels = c("CPCA", "ROBCPCA0","ROBCPCA01", "ROBCPCA02", "ROBCPCA012"))

# Create the boxplot with custom colors and ordered groups
ggplot(data, aes(x = group, y = values, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("CPCA" = "salmon", "ROBCPCA0" = "darkseagreen4","ROBCPCA01" = "goldenrod", "ROBCPCA02" = "cyan", "ROBCPCA012" = "orchid")) +
  labs(  x = "Methods", y = "True Accuracy") +
  theme_minimal()




data <- read_xlsx('/Users/maz0b/Desktop/ROBCPCA/code/Japanese_Vowl.xlsx')
# Extract the columns into individual vectors
accuracyjv <- data$tanato
accuracyrobjv <- data$tanatorob
accuracyroblagjv <- data$tanatoroblag
accuracyroblag2jv <- data$tanatoroblag2
accuracyroblag12jv <- data$tanatoroblag12

# Combine data into a data frame in long format
data_long <- data.frame(
  values = c(accuracyjv, accuracyrobjv, accuracyroblagjv, accuracyroblag2jv, accuracyroblag12jv),
  group = factor(rep(c("CPCA", "ROBCPCA0", "ROBCPCA01", "ROBCPCA02", "ROBCPCA012"), 
                     each = length(accuracyjv)))
)

# Set the desired order of the groups
data_long$group <- factor(data_long$group, levels = c("CPCA", "ROBCPCA0", "ROBCPCA01", "ROBCPCA02", "ROBCPCA012"))

ggplot(data_long, aes(x = group, y = values, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("CPCA" = "salmon", 
                               "ROBCPCA0" = "darkseagreen4",
                               "ROBCPCA01" = "goldenrod", 
                               "ROBCPCA02" = "cyan", 
                               "ROBCPCA012" = "orchid")) +
  labs(x = "Methods", y = "True Accuracy") +
  theme_minimal() +
  theme(
    text = element_text(size = 14),          # Set overall font size
    axis.text = element_text(size = 14),     # Axis tick labels
    axis.title = element_text(size = 14),    # Axis titles
    legend.text = element_text(size = 14),   # Legend text
    legend.title = element_text(size = 14),  # Legend title
    plot.title = element_text(size = 14) ,    # Plot title
    legend.key.size = unit(1, 'cm')
  )


