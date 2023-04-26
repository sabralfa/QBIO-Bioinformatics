Sys.unsetenv("R_LIBS_USER")
dir.create("RLibrary")
.libPaths()
.libPaths(paste(getwd(), "RLibrary", sep="/"))
setRepositories()
library(tidyverse)# provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #helps deal with ensembl

# we looked for our Theobroma Criollo cacao in the view available datasets in the BioMart
library(biomaRt)
cacao.anno <- useMart(biomart = "plants_mart", dataset = "tccriollo_eg_gene", host="https://plants.ensembl.org")
cacao.attributes <- listAttributes(cacao.anno)
view(cacao.attributes)
Tx.cacao <- getBM(attributes=c('ensembl_transcript_id',
                             'ensembl_gene_id', 'description'),
                mart = cacao.anno)

targets2 <- read_tsv("cacao_study.txt")# read in your study design
targets2
path2 <- file.path(targets2$sample, "abundance.tsv") # set file paths to your mapped data
path2
all(file.exists(path2)) 
Tx.cacao <- as_tibble(Tx.cacao)
Tx.cacao <- dplyr::rename(Tx.cacao, target_id = ensembl_transcript_id, 
                        gene_name = ensembl_gene_id)
Tx.cacao <- dplyr::select(Tx.cacao, "target_id", "gene_name")
Txc_gene <- tximport(path2, 
                     type = "kallisto", 
                     tx2gene = Tx.cacao, 
                     txOut = FALSE, #determines whether your data represented at transcript or gene level
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)

#### Step 2 ####
library(edgeR)
library(matrixStats)
library(cowplot)

myTPM2 <- Txc_gene$abundance # for transcripts: Txi_transcript$abundance
myCounts2 <- Txc_gene$counts
colSums(myTPM2)
colSums(myCounts2)

sampleLabels2 <- targets2$sample
#myTPM.stats <- transform(myTPM2, 
#                         SD=rowSds(myTPM2), 
#                         AVG=rowMeans(myTPM2),
#                         MED=rowMedians(myTPM2))

myDGEList2 <- DGEList(Txc_gene$counts)
log3.cpm <- cpm(myDGEList2, log=TRUE) 

log3.cpm.df <- as_tibble(log3.cpm, rownames = "geneID")
colnames(log3.cpm.df) <- c("geneID", sampleLabels2)
log3.cpm.df.pivot <- pivot_longer(log3.cpm.df, # dataframe to be pivoted
                                  cols = pathogen1:control4, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

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

cpm2 <- cpm(myDGEList2)
table(rowSums(myDGEList2$counts==0)==10)
keepers2 <- rowSums(cpm2>1)>=4 #user defined
myDGEList2.filtered <- myDGEList2[keepers2,]

log3.cpm.filtered <- cpm(myDGEList2.filtered, log=TRUE)
log3.cpm.filtered.df <- as_tibble(log3.cpm.filtered, rownames = "geneID")
colnames(log3.cpm.filtered.df) <- c("geneID", sampleLabels2)
log3.cpm.filtered.df.pivot <- pivot_longer(log3.cpm.filtered.df, # dataframe to be pivoted
                                           cols = pathogen1:control4, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)

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

myDGEList2.filtered.norm <- calcNormFactors(myDGEList2.filtered, method = "TMM")
log3.cpm.filtered.norm <- cpm(myDGEList2.filtered.norm, log=TRUE)
log3.cpm.filtered.norm.df <- as_tibble(log3.cpm.filtered.norm, rownames = "geneID")
colnames(log3.cpm.filtered.norm.df) <- c("geneID", sampleLabels2)
log3.cpm.filtered.norm.df.pivot <- pivot_longer(log3.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = pathogen1:control4, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

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

plot_grid(p1_c, p2_c, p3_c, labels = c('A', 'B', 'C'), label_size = 12)



#### Step 3 ####

library(DT) #interactive and searchable tables of our GSEA results
library(gt)
library(plotly)

group2 <- targets2$group
group2 <- factor(group2)

pca.res2 <- prcomp(t(log3.cpm.filtered.norm), scale.=F, retx=T)
pc.var2 <- pca.res2$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per2 <- round(pc.var2/sum(pc.var2)*100, 1) 
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

mydata.df2 <- mutate(log3.cpm.filtered.norm.df,
                    healthy.AVG = (control1 + control2 + control3 + control4)/4, 
                    disease.AVG = (pathogen1 + pathogen2 + pathogen3 + pathogen4)/4,
                    #now make columns comparing each of the averages above that you're interested in
                    LogFC = (disease.AVG - healthy.AVG)) %>% #note that this is the first time you've seen the 'pipe' operator
  mutate_if(is.numeric, round, 2)

datatable(mydata.df2[,c(1,10:12)], 
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



group2 <- factor(targets2$group)  # with this function we want to encode the vector group2 of our dataframe(tibble) targets2 as a factor
design2 <- model.matrix(~0 + group2)  # we are using the new factor group2 to expand it to a matrix according to the appearance of disease and healty in the 8 samples. It is then annotated as 1.
colnames(design2) <- levels(group2) # the column names are healthy and disease and we specifiy them know as levels

v.DEGList2.filtered.norm <- voom(myDGEList2.filtered.norm, design2, plot = T) # we now add the counts saved as filtered.norm to the design2, as well as the weights and transcripts of the sample 
fit2 <- lmFit(v.DEGList2.filtered.norm, design2) # now we fit the arrays by weighted or generalized least squares -> so disease and healthy
contrast.matrix2 <- makeContrasts(infection = disease - healthy,      # here we specify the contrasts between the disease and healthy -> disease is 1 and healthy is -1
                                 levels=design2)                      # and we also create a contrast matrix

fits2 <- contrasts.fit(fit2, contrast.matrix2)                        # here we want to compute estimated coefficients and standard errors for given set of contrasts. The input is the lmfit and contrast matrix. Finally, it re-orientated the fitted model object from the coefficients to 'disease' and 'healthy'.
ebFit2 <- eBayes(fits2)                                               # we know compute moderated t-statistics, F-statistic and log-odds of our fitted data
myTopHits2 <- topTable(ebFit2, adjust ="BH", coef=1, number=40000, sort.by="logFC") # with this we extract the top-ranked genes from our linear model fit object and sort them by the (absolute) coefficient representing the log-fold-change
myTopHits.df2 <- myTopHits2 %>%    # creating a dataframe(tibble) with the columns:'logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B'. It's interesting that the B-statistics is below 0; but a B-statistic of zero corresponds to a 50-50 chance that the gene is differentially expressed.
  as_tibble(rownames = "geneID")
gt(myTopHits.df2) # here we can see our t-statistics; B-statistics and the p-values

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

ggplotly(vplot2)                # So the blue square is giving us the own regulated differentially expressed genes and the red square the up regulated ones. For out p adjusted values we have pretty normal differentially expressed genes. Apparently none of our transcripts are significant.


results2 <- decideTests(ebFit2, method="global", adjust.method="BH", p.value = 0.2, lfc = 1)     # checking which genes are significantly differentially expressed for each contrast from our fitted object. We don't have any significant differentially expressed genes. The p.value here represents the alpha significance; which normally is 0.01 in medical statistics. 
colnames(v.DEGList2.filtered.norm$E) <- sampleLabels2       #specifing sampleLables according to pathogen1 to control4
diffGenes2 <- v.DEGList2.filtered.norm$E[results2[,1] !=0,] #taking from v.DEGList2.filtered.norm the table E and only the results from the first row and only the ones which unequal 0
diffGenes.df2 <- as_tibble(diffGenes2, rownames = "geneID")  # now we put diffGenes2 into a tibble
datatable(diffGenes.df2,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in Theobroma Criollo cacao',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:9), digits=2)  # just displaying the datatable with the expression values (E displays the numeric matrix of normalized expression values on the log2 scale)of the pathogen1 till control4, which are above 0. But here not a single one is revealed since they are not above significance = 0.01 and lfc = 1.

### Step 5 ###


#install.packages('gplots')
library(gplots) #for heatmaps
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots



# Carry out GO enrichment using gProfiler2 ----
# use topTable result to pick the top genes for carrying out a Gene Ontology (GO) enrichment analysis
myTopHits2 <- topTable(ebFit2, adjust ="BH", coef=1, number=50, sort.by="logFC")
myTopHits2
ebFit2
# use the 'gost' function from the gprofiler2 package to run GO enrichment analysis
gost.res2 <- gost(rownames(myTopHits2), organism = "tccriollo", correction_method = "fdr")
gost.res2
# produce an interactive manhattan plot of enriched GO terms
gostplot(gost.res2, interactive = T, capped = F) #set interactive=FALSE to get plot for publications
# produce a publication quality static manhattan plot with specific GO terms highlighted
# rerun the above gostplot function with 'interactive=F' and save to an object 'mygostplot'
publish_gostplot(
  mygostplot, #your static gostplot from above
  highlight_terms = c("GO:0034987"),
  filename = NULL,
  width = NA,
  height = NA)

#you can also generate a table of your gost results
# publish_gosttable(
#   gost.res,
#   highlight_terms = NULL,
#   use_colors = TRUE,
#   show_columns = c("source", "term_name", "term_size", "intersection_size"),
#   filename = NULL,
#   ggplot=TRUE)
# now repeat the above steps using only genes from a single module from the step 6 script, by using `rownames(myModule)`
# what is value in breaking up DEGs into modules for functional enrichment analysis?


# Competitive GSEA using CAMERA----
# for competitive tests the null hypothesis is that genes in the set are, at most, as often differentially expressed as genes outside the set
# first let's create a few signatures to test in our enrichment analysis
mySig <- rownames(myTopHits) #if your own data is from mice, wrap this in 'toupper()' to make gene symbols all caps
mySig2 <- sample((rownames(v.DEGList.filtered.norm$E)), size = 50, replace = FALSE)
collection <- list(real = mySig, fake = mySig2)
# now test for enrichment using CAMERA
camera.res <- camera(v.DEGList.filtered.norm$E, collection, design, contrast.matrix[,1]) 
camera.df <- as_tibble(camera.res, rownames = "setName")
camera.df

# Self-contained GSEA using ROAST----
# remember that for self-contained the null hypothesis is that no genes in the set are differentially expressed
mroast(v.DEGList.filtered.norm$E, collection, design, contrast=1) #mroast adjusts for multiple testing

# now repeat with an actual gene set collection
# camera requires collections to be presented as a list, rather than a tibble, so we must read in our signatures using the 'getGmt' function
broadSet.C2.ALL <- getGmt("c2.cp.kegg.v2023.1.Hs.symbols.gmt", geneIdType=SymbolIdentifier())
# Option for plants: goto http://structuralbiology.cau.edu.cn/PlantGSEA/download.php
# Download the gene set of interest for the species of interest (e.g. Ara_KEGG.txt)
# broadSet.C2.ALL <- getGmt("Ara_KEGG.txt",geneIdType=SymbolIdentifier())
# broadSet.C2.ALL <- geneIds(broadSet.C2.ALL)

#extract as a list
broadSet.C2.ALL <- geneIds(broadSet.C2.ALL)
camera.res <- camera(v.DEGList.filtered.norm$E, broadSet.C2.ALL, design, contrast.matrix[,1]) 
camera.df <- as_tibble(camera.res, rownames = "setName")
camera.df

# filter based on FDR and display as interactive table
camera.df <- filter(camera.df, FDR<=0.01)

datatable(camera.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2,4,5), digits=2)

#as before, add a variable that maps up/down regulated pathways with phenotype
camera.df <- camera.df %>%
  mutate(phenotype = case_when(
    Direction == "Up" ~ "disease",
    Direction == "Down" ~ "healthy"))

#easy to filter this list based on names of signatures using 'str_detect'
#here is an example of filtering to return anything that has 'CD8' or 'CYTOTOX' in the name of the signature
camera.df.sub <- camera.df %>%
  dplyr::filter(str_detect(setName, "CD8|CYTOTOX"))

# graph camera results as bubble chart 
ggplot(camera.df[1:25,], aes(x=phenotype, y=setName)) + 
  geom_point(aes(size=NGenes, color = Direction, alpha=-log10(FDR))) +
  theme_bw()

# Single sample GSEA using the GSVA package----
# the GSVA package offers a different way of approaching functional enrichment analysis.  
# A few comments about the approach:
# In contrast to most GSE methods, GSVA performs a change in coordinate systems,
# transforming the data from a gene by sample matrix to a gene set (signature) by sample matrix. 
# this allows for the evaluation of pathway enrichment for each sample.
# the method is both non-parametric and unsupervised
# bypasses the conventional approach of explicitly modeling phenotypes within enrichment scoring algorithms. 
# focus is therefore placed on the RELATIVE enrichment of pathways across the sample space rather than the absolute enrichment with respect to a phenotype. 
# however, with data with a moderate to small sample size (< 30), other GSE methods that explicitly include the phenotype in their model are more likely to provide greater statistical power to detect functional enrichment.

# be aware that if you choose a large MsigDB file here, this step may take a while
GSVA.res.C2CP <- gsva(v.DEGList.filtered.norm$E, #your data
                      broadSet.C2.ALL, #signatures
                      min.sz=3, max.sz=500, #criteria for filtering gene sets
                      mx.diff=FALSE,
                      method="gsva") #options for method are "gsva", ssgsea', "zscore" or "plage"

# Apply linear model to GSVA result
# now using Limma to find significantly enriched gene sets in the same way you did to find diffGenes
# this means you'll be using topTable, decideTests, etc
# note that you need to reference your design and contrast matrix here
fit.C2CP <- lmFit(GSVA.res.C2CP, design)
ebFit.C2CP <- eBayes(fit.C2CP)

# use topTable and decideTests functions to identify the differentially enriched gene sets
topPaths.C2CP <- topTable(ebFit.C2CP, adjust ="BH", coef=1, number=50, sort.by="logFC")
res.C2CP <- decideTests(ebFit.C2CP, method="global", adjust.method="BH", p.value=0.05, lfc=0)
# the summary of the decideTests result shows how many sets were enriched in induced and repressed genes in all sample types
summary(res.C2CP)

# pull out the gene sets that are differentially enriched between groups
diffSets.C2CP <- GSVA.res.C2CP[res.C2CP[,1] !=0,]
head(diffSets.C2CP)
dim(diffSets.C2CP)



# Create a heatmap of differentially expressed genes ----

heatmaply(diffSets.C2CP, 
          #dendrogram = "row",
          xlab = "Samples", ylab = "KEGG pathways", 
          main = "Responsive KEGG pathways in cutaneous leishmaniasis",
          scale = "column",
          margins = c(60,100,40,20),
          grid_color = "white",
          grid_width = 0.0000001,
          titleX = T,
          titleY = T,
          hide_colorbar = TRUE,
          branches_lwd = 0.1,
          label_names = c("Gene", "Sample:", "Value"),
          fontsize_row = 5, fontsize_col = 5,
          labCol = colnames(diffSets.C2CP),
          labRow = rownames(diffSets.C2CP),
          heatmap_layers = theme(axis.line=element_blank())
)


### essentials ###
gost2.res_up <- gost(rownames(myModule_up), organism = "Theobroma cacao (Theobroma cacao Belizian Criollo B97-61/B2)", correction_method = "fdr")
gostplot(gost2.res_up, interactive = T, capped = F) #set interactive=FALSE to get plot for publications
gost.res_down <- gost(rownames(myModule_down), organism = "Theobroma cacao (Theobroma cacao Belizian Criollo B97-61/B2)", correction_method = "fdr")
gostplot(gost.res_down, interactive = T, capped = T) #set interactive=FALSE to get plot for publications

hs_gsea_c2 <- msigdbr(species = "Theobroma cacao (Theobroma cacao Belizian Criollo B97-61/B2)", # change depending on species your data came from
                      category = "H") %>% # choose your msigdb collection of interest
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature 

# Now that you have your msigdb collections ready, prepare your data
# grab the dataframe you made in step3 script
# Pull out just the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
mydata.df.sub2 <- dplyr::select(mydata.df2, geneID, LogFC)
mydata.gsea2 <- mydata.df.sub2$LogFC
names(mydata.gsea2) <- as.character(mydata.df.sub2$geneID)
mydata.gsea2 <- sort(mydata.gsea2, decreasing = TRUE)

# run GSEA using the 'GSEA' function from clusterProfiler
myGSEA.res2 <- GSEA(mydata.gsea2, TERM2GENE=hs_gsea_c2, verbose=FALSE)
myGSEA.df2 <- as_tibble(myGSEA.res2@result)

# view results as an interactive table
datatable(myGSEA.df2, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(3:10), digits=2)
# create enrichment plots using the enrichplot package
gseaplot2(myGSEA.res2, 
          geneSetID = 47, #can choose multiple signatures to overlay in this plot
          pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
          title = myGSEA.res2$Description[47]) #can also turn off this title

# add a variable to this result that matches enrichment direction with phenotype
myGSEA.df2 <- myGSEA.df2 %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "disease",
    NES < 0 ~ "healthy"))

# create 'bubble plot' to summarize y signatures across x phenotypes
ggplot(myGSEA.df2[1:20,], aes(x=phenotype, y=ID)) + 
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()


