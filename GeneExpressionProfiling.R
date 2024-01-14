#1
library(BioStudies)
library(SummarizedExperiment)

dir.create("data") 

getBio("E-MEXP-3396", path = "data") 


sdrf <- read.table("data/E-MEXP-3396.sdrf.txt", header = TRUE, as.is = TRUE,
                   row.names = 1, sep = "\t")

library(oligo)

cel_files <- sdrf$Array.Data.File
rawData <- read.celfiles(filenames=paste0("data/E-MEXP-3396.sdrf.txt",cel_files))

#Normalize data
eset <- rma(rawData)

#Creating a factor for cell type
group <- factor(sub(pattern = " .+",
                    replacement = "",
                    x = sdrf$Factor.Value.CELL_TYPE.))

#Creating a factor for the pairs
Pairs <- factor(sub(pattern = ".+_.+_",
                    replacement = "", x = sdrf$Characteristics.Individual.))


boxplot(exprs(eset),
        col = rep(c("red", "blue"), each =5),
        cex.axis = 0.6,
        las = 3,
        main = "Treatment") #Inspect the RNA-data for normalization


par(mfrow = c(1,2))
hist(exprs(eset))
qqnorm(exprs(eset))
qqline(exprs(eset), lwd=1, col="red")


library(limma)
#Create desing matrix
studyDesign = data.frame(Group = group, Pairs = Pairs)
studyDesign

designMatrix <- model.matrix(~ 0 + Group + Pairs, data = studyDesign)
designMatrix

colnames(designMatrix)[1:2] <- c("Artrial", "Ventricular")
designMatrix

#Fit linear model
fit <- lmFit(exprs(eset), designMatrix)

#Create contrast matrix
contrastMatrix <- makeContrasts(VentricularVsArtrial = Ventricular - Artrial,
                                levels = designMatrix)
contrastMatrix

#Perform testing
fit <- eBayes(fit, robust = TRUE)

#Filter results
results <- topTable(fit, coef = "VentricularVsArtrial",,number=nrow(eset),
                    p.value = 0.05, lfc = 1)
results[1:6,1:6]
