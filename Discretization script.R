# Discretization
# Discretization refers to the process of converting continuous data or variables into discrete categories or intervals. 
# This is often done to simplify the analysis, make data more manageable, or to prepare it for specific statistical or machine learning

# Select a dataset in /Datasets
Dataset <- t(read.delim("Datasets/expr_GSE139078.txt", header=T, sep=" ")) # Example

i=1 # i = first col, the cols are the genes.
numcol <- ncol(Dataset) # number of columns in the entire dataset

Dataset_BN = matrix(ncol = numcol,
                      nrow = nrow(Dataset)) # Matrix, that will store the discretized dataset.

for (i in i:numcol) {
  Dataset_BN[,i] <- ifelse(Dataset[,i] < quantile(Dataset[,i], probs = 0.25),
                                       "A",
                                       ifelse(Dataset[,i] < quantile(Dataset[,i], probs = 0.5),
                                              "B",
                                              ifelse(Dataset[,i] < quantile(Dataset[,i], probs = 0.75),
                                                     "C", "D")))
} 
# The first quantile <25%, will receive the letter A, between 25% and 50% the letter B, 50% to 75% the letter C and >75% the letter D.

