#### Timepoint 24h
Sys.unsetenv("R_LIBS_USER")
dir.create("RLibrary")
.libPaths()
.libPaths(paste(getwd(), "RLibrary", sep="/"))
setRepositories()
library(tidyverse)# Provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # Package for getting Kallisto results into R
library(ensembldb) # Helps deal with ensembl

# Looking for the name of the available Theobroma Criollo cacao datasets in BioMart
library(biomaRt)
cac.anno <- useMart(biomart = "plants_mart", dataset = "tccriollo_eg_gene", host="https://plants.ensembl.org")
cac.attributes <- listAttributes(cac.anno)
view(cac.attributes)
Tx.cac <- getBM(attributes=c('ensembl_transcript_id',
                               'ensembl_gene_id', 'description'),
                  mart = cac.anno)

targets3 <- read_tsv("ca_study.txt") 
path3 <- file.path("T24h",targets3$sample, "abundance.tsv") # Set file paths to mapped data
all(file.exists(path3)) 
Tx.cac <- as_tibble(Tx.cac)
Tx.cac <- dplyr::rename(Tx.cac, target_id = ensembl_transcript_id, 
                          gene_name = ensembl_gene_id)
Tx.cac <- dplyr::select(Tx.cac, "target_id", "gene_name")
Txcac_gene <- tximport(path3, 
                     type = "kallisto", 
                     tx2gene = Tx.cac, 
                     txOut = FALSE, # Determines whether the data is represented at transcript or gene level
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)

#### Step 2 ####
library(edgeR)
library(matrixStats)
library(cowplot)

myTPM3 <- Txcac_gene$abundance # For transcripts: Txi_transcript$abundance
myCounts3 <- Txcac_gene$counts
colSums(myTPM3)
colSums(myCounts3)

sampleLabels3 <- targets3$sample


myDGEList3 <- DGEList(Txcac_gene$counts)
log4.cpm <- cpm(myDGEList3, log=TRUE) 

log4.cpm.df <- as_tibble(log4.cpm, rownames = "geneID")
colnames(log4.cpm.df) <- c("geneID", sampleLabels3)
log4.cpm.df.pivot <- pivot_longer(log4.cpm.df, # Dataframe to be pivoted
                                  cols = path2:cont4, # Column names to be stored as a SINGLE variable
                                  names_to = "samples", # Name of that new variable (column)
                                  values_to = "expression") # Name of new variable (column) storing all the values (data)

p1_cac <- ggplot(log4.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log3 expression", x = "sample",
       title="Log3 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

cpm3 <- cpm(myDGEList3)
table(rowSums(myDGEList3$counts==0)==10)
keepers2 <- rowSums(cpm3>1)>=4 # User defined
myDGEList3.filtered <- myDGEList3[keepers2,]

log4.cpm.filtered <- cpm(myDGEList3.filtered, log=TRUE)
log4.cpm.filtered.df <- as_tibble(log4.cpm.filtered, rownames = "geneID")
colnames(log4.cpm.filtered.df) <- c("geneID", sampleLabels3)
log4.cpm.filtered.df.pivot <- pivot_longer(log4.cpm.filtered.df, # Dataframe to be pivoted
                                           cols = path2:cont4, # Column names to be stored as a SINGLE variable
                                           names_to = "samples", # Name of that new variable (column)
                                           values_to = "expression") # Name of new variable (column) storing all the values (data)

p2_cac <- ggplot(log4.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log3 expression", x = "sample",
       title="Log3 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

myDGEList3.filtered.norm <- calcNormFactors(myDGEList3.filtered, method = "TMM")
log4.cpm.filtered.norm <- cpm(myDGEList3.filtered.norm, log=TRUE)
log4.cpm.filtered.norm.df <- as_tibble(log4.cpm.filtered.norm, rownames = "geneID")
colnames(log4.cpm.filtered.norm.df) <- c("geneID", sampleLabels3)
log4.cpm.filtered.norm.df.pivot <- pivot_longer(log4.cpm.filtered.norm.df, # Dataframe to be pivoted
                                                cols = path2:cont4, # Column names to be stored as a SINGLE variable
                                                names_to = "samples", # Name of that new variable (column)
                                                values_to = "expression") # Name of new variable (column) storing all the values (data)

p3_cac <- ggplot(log4.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log3 expression", x = "sample",
       title="Log3 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

plot_grid(p1_cac, p2_cac, p3_cac, labels = c('A', 'B', 'C'), label_size = 12)

### Step3 ###
library(DT) # Interactive and searchable tables of our GSEA results
library(gt)
library(plotly)

group3 <- targets3$group
group3 <- factor(group3)

pca.res3 <- prcomp(t(log4.cpm.filtered.norm), scale.=F, retx=T)
pc.var3 <- pca.res3$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per3 <- round(pc.var3/sum(pc.var3)*100, 1) 
pca.res.df3 <- as_tibble(pca.res3$x)
pca.plot3 <- ggplot(pca.res.df3) +
  aes(x=PC1, y=PC2, label=sampleLabels3, color = group3) +
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per3[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per3[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggplotly(pca.plot3)
mydata.df3 <- mutate(log4.cpm.filtered.norm.df,
                     healthy.AVG = (cont1 + cont2 + cont3 + cont4)/4, 
                     disease.AVG = (path2 + path3 + path4 + path5)/4,
                     # Now make columns comparing each of the averages above that you're interested in
                     LogFC = (disease.AVG - healthy.AVG)) %>% 
  mutate_if(is.numeric, round, 2)

datatable(mydata.df3[,c(1,10:12)], 
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         #dom = "Blfrtip", 
                         #buttons = c("copy", "csv", "excel"),
                         lengthMenu = c("10", "25", "50", "100")))

### Step 4 ###
library(limma)

group3 <- factor(targets3$group)  # With this function the vector group2 of our dataframe(tibble) targets2 is encoded as a factor
design3 <- model.matrix(~0 + group3)  # The new factor group2 is used to expand it to a matrix according to the appearance of disease and healty in the 8 samples. It is then annotated as 1.
colnames(design3) <- levels(group3) # The column names are healthy and disease and specified as levels

v.DEGList3.filtered.norm <- voom(myDGEList3.filtered.norm, design3, plot = T) # The counts are saved as filtered.norm to the design2, as well as the weights and transcripts of the sample 
fit3 <- lmFit(v.DEGList3.filtered.norm, design3) # The arrays are fitted by weighted or generalized least squares -> so disease and healthy
contrast.matrix3 <- makeContrasts(infection = disease - healthy,      # Specify the contrasts between the disease and healthy -> disease is 1 and healthy is -1
                                  levels=design3)                      # Create a contrast matrix

fits3 <- contrasts.fit(fit3, contrast.matrix3)                        # Compute estimated coefficients and standard errors for given set of contrasts. The input is the lmfit and contrast matrix. Finally, it re-orientated the fitted model object from the coefficients to 'disease' and 'healthy'.
ebFit3 <- eBayes(fits3)                                               # Compute moderated t-statistics, F-statistic and log-odds of our fitted data
myTopHits3 <- topTable(ebFit3, adjust ="BH", coef=1, number=40000, sort.by="logFC") # Extract the top-ranked genes from our linear model fit object and sort them by the (absolute) coefficient representing the log-fold-change
myTopHits.df3 <- myTopHits3 %>%    # Creating a dataframe(tibble) with the columns:'logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B'. It's interesting that the B-statistics is below 0; but a B-statistic of zero corresponds to a 50-50 chance that the gene is differentially expressed.
  as_tibble(rownames = "geneID")
gt(myTopHits.df3) # View t-statistics; B-statistics and the p-values

vplot3 <- ggplot(myTopHits.df3) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=1) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +      
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  annotate("rect", xmin = 1, xmax = 5, ymin = -log10(0.01), ymax = 5.5, alpha=.3, fill="#BE684D") +
  annotate("rect", xmin = -1, xmax = -5, ymin = -log10(0.01), ymax = 5.5, alpha=.3, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "Theobroma Criollo cacao",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

ggplotly(vplot3)                # The blue square visualizes the down regulated differentially expressed genes and the red square the up regulated ones. None of the transcripts are significant.


results3 <- decideTests(ebFit3, method="global", adjust.method="BH", p.value = 0.01, lfc = 1)     # Checking which genes are significantly differentially expressed for each contrast from the fitted object. No significant differentially expressed genes. The p.value represents the alpha significance; which is 0.01. 
colnames(v.DEGList3.filtered.norm$E) <- sampleLabels3       # Specifing sampleLables according to pathogen1 to control4
diffGenes3 <- v.DEGList3.filtered.norm$E[results3[,1] !=0,] # Taking from v.DEGList3.filtered.norm the table E and only the results from the first row which are unequal 0
diffGenes.df3 <- as_tibble(diffGenes3, rownames = "geneID")  # Putting diffGenes3 into a tibble
datatable(diffGenes.df3,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in Theobroma Criollo cacao',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:9), digits=2)  # Displaying the datatable with the expression values (E displays the numeric matrix of normalized expression values on the log2 scale)of the pathogen1 till control4, which are above 0. But here not a single one is revealed since they are not above significance = 0.01 and lfc = 1.

### Step 5 ###

myTopHits3 <- topTable(ebFit3, adjust ="BH", coef=1, number=50, sort.by="logFC")
print((rownames(myTopHits3)))

# Conducting the GO enrichment analysis via gprofiler webpage, no results are shown for the TopHits3.
# So the results of the same strain at different timepoints doesn't vary significantly here either.
