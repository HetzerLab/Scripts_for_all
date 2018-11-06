# Let's do the differential expression in R. 
# The example below is for using edgeR, you can also do this in RStudio server. 
# For that type http://tardis.salk.edu:8787 into a web browser and login with your credentials.
# Navigate to the R_analysis folder and set that as your working directory.
# If not using RStudio, type R in the console below to start R in the command line.
# R
setwd("~/RNAseq_example/R_analysis")
library(edgeR)
library(readr)
countTableraw <- read_delim("countTableraw.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(countTableraw)
newcolnames <- colnames(countTableraw)
newcolnames[1] <- "TranscriptID"
newcolnames[9:10] <- c("A1_37C", "A2_43C")
colnames(countTableraw) <- newcolnames
colnames(countTableraw)

x <- as.matrix(countTableraw[,9:10]) # Change to indicate the cols for your samples
group <- factor(c(1,2),labels = c("Temp37C","Temp42C")) # Change to indicate your sample groups
genes <- countTableraw[,1:8]
y <- DGEList(counts = x, group = group, genes = genes)

y <- calcNormFactors(y)
design <- model.matrix(~group)

# Uncomment and run one of the two below:
# y <- estimateDisp(y,design) # If you do not have replicates this will not work. Use the commented function below
# y$common.dispersion <- 0.05 # If you have replicates it's best to use the function above.

de.com <- exactTest(y)
results <- topTags(de.com,n = length(x[,1]))

# write the output to a text file
results <- results$table
write.table(results,file="outputFile_edgeR.txt",sep="\t",row.names = FALSE)