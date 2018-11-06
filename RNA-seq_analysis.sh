# RNA-seq start up analysis

# Log in to the server
# Replace salk ID with yours.

ssh salkID@tardis.salk.edu

# Use screen so you can shut down the computer if you'd like and still leave things running
# Use screen -rd to get back to the same place.

# Downloading your RNA-seq files to the server:
# Replace the URLs with the ones for your files (likely in an email from Nasun?)
# For this example I'm using 2 files that Robbie and Hannah provided
# File 1 22RC @37C: http://igc2.snl.salk.edu/illumina/runs/160412-SE75/Hetzer_Robbie/AR1/A1_S1_R1_001.fastq.gz
# File 2 22RC @43C: http://igc2.snl.salk.edu/illumina/runs/160412-SE75/Hetzer_Robbie/AR2/A2_S2_R1_001.fastq.gz

# I'm also creating a folder to download them into:

mkdir RNAseq_example
cd RNAseq_example

curl -O http://igc2.snl.salk.edu/illumina/runs/160412-SE75/Hetzer_Robbie/AR1/A1_S1_R1_001.fastq.gz
curl -O http://igc2.snl.salk.edu/illumina/runs/160412-SE75/Hetzer_Robbie/AR2/A2_S2_R1_001.fastq.gz

# Nasun does the fastqc for our files automaticaly, so we can skip that part.
# If that's not the case, uncomment the section below and run the following:

# fastqc A1_S1_R1_001.fastq.gz A2_S2_R1_001.fastq.gz

# Let's organize the files a bit:

mkdir fastq_files
mkdir fastqc_results

mv *.fastq.gz fastq_files
mv *fastqc.* fastqc_results

# Assuming the fastqc results looked ok, we move ahead with the alignment.

mkdir alignments
for file in fastq_files/*.fastq.gz
do
STAR --readFilesCommand zcat --genomeDir /data/genomes/STARgenomes/hg19.star --runThreadN 24 --readFilesIn ${file} --outFileNamePrefix alignments/${file#fastq_files} --genomeLoad LoadAndKeep
done

cd alignments

# Let's organize the files a bit better:

mkdir Logs
for file in *Log.*
do
mv $file Logs/
done

mkdir SJs
for file in *.gzSJ.*
do
mv $file SJs/
done

mkdir sam_files
for file in *.sam*
do
mv $file sam_files/
done

#Next we sort the resulting sam files and convert them to bam.

cd sam_files
for file in *.sam*
do
samToSortedBam $file
done

# We can use the sorted bam files to create bedgraphs.
# These bedgraph files can be uploaded to the UCSC genome broweser for vizualization.

for file in *Aligned.out.bam
do
bamToBedGraph hg19 $file
done

# Time to organize our files a bit more again.
cd ..
mkdir bam_files
mkdir bedgraphs

cd sam_files

for file in *.bam
do
mv $file ~/RNAseq_example/alignments/bam_files/
done

for file in *.bedGraph.gz
do
mv $file ~/RNAseq_example/alignments/bedgraphs/
done

# To move ahead with the RNA-seq analysis using HOMER we need to create tag directories.

for file in *Aligned.out.sam
do
makeTagDirectory ${file%fastq.gzAligned.out.sam}TagDir $file -format sam -flip
done

# Now we agregate the reads into a tabular format.
# Replate the name of my TagDirs below with yours.
# If you have more than 2 groups of samples you plan to compare order here is important. Your controls should be the first samples.

analyzeRepeats.pl rna hg19 -count exons -condenseGenes -raw -d A1_S1_R1_001.TagDir/ A2_S2_R1_001.TagDir/ > countTableraw.txt

# Finally we do differential expression analysis.
# This will call DEseq2 in R and it need the raw count values to work best (not FPKM)
# If you have more than 2 types of samples the command below will compare the first sample kind against all the others.
# Replace the name of my samples with those appropriate for yours. The sample names must be in the same order as what you provided in the function above.
# One quick note, this is what you usually will run assuming you have replicates of your RNA-seq data. IF YOU DO NOT HAVE REPLICATES this is likely not going to give you any genes at all (very conservative estimate of variance).
# You can than use a less conservative function call using edgeR, see the second function below and use that.
# I have used both here just as an example for comparison. But usually you will choose the one that best matches your data.

#getDiffExpression.pl countTableraw.txt Temp37C Temp43C > diffOutputDEseq2.txt 

#getDiffExpression.pl countTableraw.txt Temp37C Temp43C -edgeR -simpleNorm -dispersion 0.05 > diffOutputedgeR.txt 

# A little more reorganizing:
cd ../..
mkdir R_analysis
mv /home/jcapitanio/RNAseq_example/alignments/sam_files/countTableraw.txt /home/jcapitanio/RNAseq_example/R_analysis/countTableraw.txt
cd R_analysis

# Let's do the differential expression in R. 
# The example below is for using edgeR, you can also do this in RStudio server. 
# For that type http://tardis.salk.edu:8787 into a web browser and login with your credentials.
# Navigate to the R_analysis folder and set that as your working directory.
# If not using RStudio, type R in the console below to start R in the command line.
R
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
quit()
# You should now have a table with your genes and FDR values indicating if they were differentially expressed.
# You can download that with filezilla and move ahead with your analysis.




