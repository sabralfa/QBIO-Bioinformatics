Sys.unsetenv("R_LIBS_USER")
dir.create("RLibrary")
.libPaths()
.libPaths(paste(getwd(), "RLibrary", sep="/"))
setRepositories()
library(tidyverse) #provides a collection of packages which help in data science and provide some methods for analysis 
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #gives access to ensembl and biomart
library(edgeR) #allows to conduct statistical analysis specifically for differential gene expression 
library(matrixStats) #function working on rows and columns of a matrix to conduct basic statistics 
library(cowplot) #improvement of ggplot to produce publication-quality graphs 
library(DT) # DataTables:data objects created in R (matrices or data frames) displayed as tables on HTML pages
library(plotly) # make interactive graphics from graphs of ggplot2 
library(gt) #makes presentable and easy to read tables 
library(limma) # mainly used for  fitting a broad class of statistical models called “linear models” 
library(biomaRt)
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots
library(qusage) # Quantitative Set Analysis for Gene Expression
library(heatmaply)

#### Step 1 ####

# Looking for the name of the available Theobroma Criollo cacao datasets in BioMart
cacao.anno <- useMart(biomart = "plants_mart", dataset = "tccriollo_eg_gene", host="https://plants.ensembl.org") # importing the annotation file for Cacao 
cacao.attributes <- listAttributes(cacao.anno)  # saving all obtained attributes obtained from the annotation file 
view(cacao.attributes)
Tx.cacao <- getBM(attributes=c('ensembl_transcript_id',       # querying the annotation file only for the specified attributes and saving the respectively belonging data 
                             'ensembl_gene_id', 'description'),
                mart = cacao.anno)

targets2 <- read_tsv("cacao_study.txt") # importing study design file as “targets2”
path2 <- file.path(targets2$sample, "abundance.tsv") # Setting the path to the respective abundance files created by Kallisto 
all(file.exists(path2)) 
Tx.cacao <- as_tibble(Tx.cacao) # converting into tibble format 
Tx.cacao <- dplyr::rename(Tx.cacao, target_id = ensembl_transcript_id, # renaming attributes for later accession 
                        gene_name = ensembl_gene_id)
Tx.cacao <- dplyr::select(Tx.cacao, "target_id", "gene_name")   # for further analysis and accession of the data it is important that the attributes are stored in this exact order in the data frame 

# here the Kallisto data is imported and the data counts at gene level are obtained and stored in this dataframe. Further analysis can now be conducted on this basis  
Txc_gene <- tximport(path2, 
                     type = "kallisto", 
                     tx2gene = Tx.cacao, 
                     txOut = FALSE, #determines whether your data represented at transcript or gene level
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)

#### Step 2 ####

myTPM2 <- Txc_gene$abundance # for transcripts: Txi_transcript$abundance
myCounts2 <- Txc_gene$counts
colSums(myTPM2)
colSums(myCounts2)

sampleLabels2 <- targets2$sample  # this list keeps the samples names for labelling use afterwards 
myDGEList2 <- DGEList(Txc_gene$counts)   # Storing all counts for the respective transcripts in a DGEList which is a compatible format for gene expression in R
log3.cpm <- cpm(myDGEList2, log=TRUE)  # Transforming the counts in counts per million with an additional log call to transform the data in a better visible format. This makes analysis and specifically visualisation easier.

log3.cpm.df <- as_tibble(log3.cpm, rownames = "geneID") # A dataframe is made from the DGEList where rows are respective genes and columns refer to samples we took. Counts per million are stored as data.
colnames(log3.cpm.df) <- c("geneID", sampleLabels2)
log3.cpm.df.pivot <- pivot_longer(log3.cpm.df, # The pivot call rearranges the dataframe into a longer one (having more rows, less columns). The resulting dataframe has stored each transcript 8 times to refer to the respective samples and the according expression of it. 
                                  cols = pathogen1:control4, # column names to be stored as a SINGLE variable
                                  names_to = "samples", 
                                  values_to = "expression") 

# Plotting the just created dataframe and storing it as p1 for later usage. The plot shows to be not descriptive enough thus further transformations need to be performed
p1_c <- ggplot(log3.cpm.df.pivot) +
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
# The data needs to be filtered. Thus in this step all transcripts which have no belonging counts and thus are not expressed in any of the samples are filtered out of the dataframe.  
cpm2 <- cpm(myDGEList2)
table(rowSums(myDGEList2$counts==0)==10)
keepers2 <- rowSums(cpm2>1)>=4 #user defined
myDGEList2.filtered <- myDGEList2[keepers2,]

# Conversion of the filtered DGEList into a dataframe with gene Ids in the rows and sample names as columns as done before 
log3.cpm.filtered <- cpm(myDGEList2.filtered, log=TRUE)
log3.cpm.filtered.df <- as_tibble(log3.cpm.filtered, rownames = "geneID")
colnames(log3.cpm.filtered.df) <- c("geneID", sampleLabels2)
# Using the pivot call on the filtered dataframe for transformation in the same manner as for the unfiltered data
log3.cpm.filtered.df.pivot <- pivot_longer(log3.cpm.filtered.df, 
                                           cols = pathogen1:control4, 
                                           names_to = "samples", 
                                           values_to = "expression") 
# Plotting the filtered data and storing it under p2 for further usage. Still envisioning the dataframe reveals that normalisation is missing. 
p2_c <- ggplot(log3.cpm.filtered.df.pivot) +
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

myDGEList2.filtered.norm <- calcNormFactors(myDGEList2.filtered, method = "TMM") # This function is normalising the data into effective library sizes
# Transforming the normalised data into a dataframe with the rows referring to the gene ids and columns to the samples given
log3.cpm.filtered.norm <- cpm(myDGEList2.filtered.norm, log=TRUE)
log3.cpm.filtered.norm.df <- as_tibble(log3.cpm.filtered.norm, rownames = "geneID")
colnames(log3.cpm.filtered.norm.df) <- c("geneID", sampleLabels2)
# Using the pivot call on the normalised data in the same manner as in the steps beforehand.
log3.cpm.filtered.norm.df.pivot <- pivot_longer(log3.cpm.filtered.norm.df, 
                                                cols = pathogen1:control4, 
                                                names_to = "samples", 
                                                values_to = "expression") 

# Plotting the now filtered and normalised data as a violin plot. No further changes need to be induced. Assigning the plot to p3. 
p3_c <- ggplot(log3.cpm.filtered.norm.df.pivot) +
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

# Plotting all three generated plots into one picture to have direct comparison on how the respective methods performed transformed the data and possibly improved visualisation of it. 
plot_grid(p1_c, p2_c, p3_c, labels = c('A', 'B', 'C'), label_size = 12)



#### Step 3 ####

group2 <- targets2$group
group2 <- factor(group2)

distance_cacao <- dist(t(log3.cpm.filtered.norm), method = "maximum") # Transforms the dataframe into a data matrix in order to implement hierarchical clustering
# Cluster dendrogram
clusters_cacao <- hclust(distance_cacao, method = "average") 
plot(clusters_cacao, labels=sampleLabels2)


pca.res2 <- prcomp(t(log3.cpm.filtered.norm), scale.=F, retx=T)
summary(pca.res2) # prints variance summary for all principal components
# PCA attempts to put emphasis on variation and extract strong patterns in a provided dataset.
# PC indicates principal components, and PC1 indicates the axis that is associated with the most variation of data, while PC2 
# represents the axis, which the data alongside portrays the second most variation, and etc.


pc.var2 <- pca.res2$sdev^2  # sdev^2 refers to the eigenvalues of the PCA 
pc.per2 <- round(pc.var2/sum(pc.var2)*100, 1) # use eigenvalues to calculate the percentage of variance linked to each PC, in order to obtain information about the relative value instead of the absolute value pc.per2

# in a plot where PC1 and PC2 is plotted against each other,it is visualised how much each sample (8 in total) contributes to PC1 and PC2
pca.res.df2 <- as_tibble(pca.res2$x)
pca.plot2 <- ggplot(pca.res.df2) +
  aes(x=PC1, y=PC2, label=sampleLabels2, color = group2) +
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per2[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per2[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggplotly(pca.plot2)

# Create a PCA 'small multiples' chart ----
# Visualise all the PCA loadings; the impact of each sample on each principal component
pca.res.df2 <- pca.res2$x[,1:8] %>% 
  as_tibble() %>%
  add_column(sample = sampleLabels2,
             group = group2)
pca.pivot2 <- pivot_longer(pca.res.df2, # Selected dataframe
                                cols = PC1:PC8, # Column names stored as a single variable
                                names_to = "PC", # Name of new variable mentioned
                                values_to = "loadings") # Name of new variable that stores the data
ggplot(pca.pivot2) +
  aes(x=sample, y=loadings, fill=group) + 
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()
# PC1 spans over the whole data set, and the coverage of the data decreases as the PC value becomes higher
# The dispersion for PC8 is too small and rounded up to zero, that it is not portrayed in the graph

# Use dplyr 'mutate' function to add new columns
mydata.df2 <- mutate(log3.cpm.filtered.norm.df,
                    healthy.AVG = (control1 + control2 + control3 + control4)/4, 
                    disease.AVG = (pathogen1 + pathogen2 + pathogen3 + pathogen4)/4,
                    # Now make columns comparing each of the averages above that you're interested in
                    LogFC = (disease.AVG - healthy.AVG)) %>% 
  mutate_if(is.numeric, round, 2)

# Using arrange to sort rows based on the values of the indicated column
# Using select to view the columns of interest and exclude the others
data.sort2 <- mydata.df2 %>%
  dplyr::arrange(desc(LogFC)) %>% 
  dplyr::select(geneID, LogFC)

# Use DT package to make an interactive table  ----
datatable(mydata.df2[,c(1,10:12)], # for  geneID, healthy.AVG, disease.AVG, LogFC
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         #dom = "Blfrtip", 
                         #buttons = c("copy", "csv", "excel"),
                         lengthMenu = c("10", "25", "50", "100")))

# Use plotly to make an interactive scatter plot  -----
# First make a ggplot 
plot2 <- ggplot(mydata.df2) + 
  aes(x=healthy.AVG, y=disease.AVG,
      text = paste("Symbol:", geneID)) + # Include the information of the geneID for every point, visualised by mouseover tooltip
  geom_point(shape=16, size=1) +
  ggtitle("disease vs. healthy") +
  theme_bw()

# Then convert this ggplot into an interactive plot
ggplotly(plot2)

### Step 4 ###

group2 <- factor(targets2$group)  # This function encodes the vector group2 as factor from the dataframe targets2
design2 <- model.matrix(~0 + group2)   # Using the new factor group2 to expand it to a matrix according to the appearance of disease and healthy in the 8 samples. It is then annotated as 1.
colnames(design2) <- levels(group2) # the column names are healthy and disease and are specified as levels 


v.DEGList2.filtered.norm <- voom(myDGEList2.filtered.norm, design2, plot = T) # The counts added are saved as filtered.norm to the design2, as well as the weights and transcripts of the sample 
fit2 <- lmFit(v.DEGList2.filtered.norm, design2) # Fit the arrays by weighted or generalized least squares -> so disease and healthy
contrast.matrix2 <- makeContrasts(infection = disease - healthy,      # The contrast between the disease and healthy is specified-> disease is 1 and healthy is -1 
                                 levels=design2)                      # Create a contrast matrix


fits2 <- contrasts.fit(fit2, contrast.matrix2)                        # In this call estimated coefficients and standard errors for given set of contrasts are estimated. The input is the lmfit and contrast matrix. Finally, it re-orientated the fitted model object from the coefficients to 'disease' and 'healthy'.
ebFit2 <- eBayes(fits2)                                               # Compute moderated t-statistics, F-statistic and log-odds of the fitted data
myTopHits2 <- topTable(ebFit2, adjust ="BH", coef=1, number=40000, sort.by="logFC") # Extraction of the top-ranked genes from the linear model fit object and sorting them by the (absolute) coefficient representing the log-fold-change
myTopHits.df2 <- myTopHits2 %>%    # creating a dataframe(tibble) with the columns:'logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B'. It's interesting that the B-statistics is below 0; but a B-statistic of zero corresponds to a 50-50 chance that the gene is differentially expressed.
  as_tibble(rownames = "geneID")
gt(myTopHits.df2) # The t-statistics; B-statistics and the p-values are portrayed here 


vplot2 <- ggplot(myTopHits.df2) +
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

ggplotly(vplot2)                # The blue square is presenting the own regulated differentially expressed genes and the red square the up regulated ones. None of the transcripts are significant.


results2 <- decideTests(ebFit2, method="global", adjust.method="BH", p.value = 0.01, lfc = 1)     # Checking which genes are significantly differentially expressed for each contrast from our fitted object. We don't have any significant differentially expressed genes. The p.value here represents the alpha significance; which is chosen here 0.05. 
colnames(v.DEGList2.filtered.norm$E) <- sampleLabels2       # Specifing sampleLables according to pathogen1 to control4
diffGenes2 <- v.DEGList2.filtered.norm$E[results2[,1] !=0,] # Taking from v.DEGList2.filtered.norm the table E and only the results from the first row which are unequal 0
diffGenes.df2 <- as_tibble(diffGenes2, rownames = "geneID")  # Now we transform diffGenes2 into tibble
datatable(diffGenes.df2,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in Theobroma Criollo cacao',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:9), digits=2)  # Displaying the datatable with the expression values (E displays the numeric matrix of normalized expression values on the log2 scale)of the pathogen1 till control4, which are above 0. But here not a single one is revealed since they are not above significance = 0.01 and lfc = 1.


### Step 5 ###

# Carry out GO enrichment using gProfiler2 ----
# Use topTable result to pick the top genes for carrying out a Gene Ontology (GO) enrichment analysis
myTopHits2 <- topTable(ebFit2, adjust ="BH", coef=1, number=50, sort.by="logFC")
print((rownames(myTopHits2)))

### The following steps were conducted on the webpage, since in Rstudio the analysis didn't work. 
### The results of the GO enrichment analysis on the webpage were actually non-existing.
### Thus, a heatmap and plotting the results is not necessary and has not been done.

### The following couldn't be conducted since no results were obtained from GO enrichment analysis.

