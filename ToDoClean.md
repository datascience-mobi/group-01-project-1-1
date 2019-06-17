# Data Cleanup

* Remove data which are not brain cancer from annotation
```
brain.anno = allDepMapData$annotation[allDepMapData$annotation$Subtype.Disease == "Glioblastoma", ]
```

* make cell line identifier be the rownames in annotation
```
rownames(brain.anno) = brain.anno$DepMap_ID
brain.anno = brain.anno[, -1]
```

* remove irrelevant cell lines from the other tables as well
```
cell.lines = dput(rownames(brain.anno))
exp.clean = allDepMapData$expression[, -which(colnames(allDepMapData$expression) %!in% cell.lines)]
copy.clean = allDepMapData$copynumber[, -which(colnames(allDepMapData$copynumber) %!in% cell.lines)]
relevant.mutations = subset(allDepMapData$mutation, names(allDepMapData$mutation) %in% cell.lines)
ceres.clean = allDepMapData$kd.ceres[, -which(colnames(allDepMapData$kd.ceres) %!in% cell.lines)]
prob.clean = allDepMapData$kd.prob[, -which(colnames(allDepMapData$kd.prob) %!in% cell.lines)]
    # library "operators" is needed for "%!in%"
```
* remove NA and values out of range (data contain a lot of 0 and negative values(esp. copynumber where x >=0))
```
sum(is.na(ceres.clean)) //0
sum(is.na(prob.clean)) //0
sum(is.na(copy.clean)) //1780
sum(is.na(exp.clean)) //0
sum(is.na(relevant.mutations)) //0, some matrices in this list have some though
copy.clean = copy.clean[-which(apply(copy.clean, 1, function(x) {sum(is.na(x))}) > 0), ]
sum(is.na(copy.clean)) //0
```
  +  expression matrix is supposed to be > 0
  ```
  sum(exp.clean < 0) //0
  ```
  +  kd.prob matrix is supposed to be 0<=x<=1
  ``` 
  sum(prob.clean < 0) //0
  sum(prob.clean > 1) //0
  ```
  +  check for values = 0
  
  ```
  sum(ceres.clean == 0) //28
  sum(prob.clean == 0) //80
  sum(copy.clean == 0) //0
  sum(exp.clean == 0) // 644705
  ```
* scan for outliers via boxplots and threshold or remove if necessary
```
boxplot(ceres.clean) //two values identified as outliers
ceres.transf[ceres.clean > 3] = 3 //thresholding
ceres.transf = ceres.clean[-which(ceres.clean > 3), ] //removing
```

* compare gene data availability between data sets
```
dim(copy.clean) == dim(exp.clean)
gene.data.co = c(rownames(copy.clean))
gene.data.ex = c(rownames(exp.clean))
exp.clean = exp.clean[-which(rownames(exp.clean) %!in% gene.data.co),]
copy.clean = copy.clean[-which(rownames(copy.clean) %!in% gene.data.ex),]
genes.clean = rownames(exp.clean)
ceres.clean = ceres.clean[-which(rownames(ceres.clean) %!in% genes.clean),]
prob.clean = prob.clean[-which(rownames(prob.clean) %!in% genes.clean),]
genes.clean = rownames(ceres.clean)
exp.clean = exp.clean[-which(rownames(exp.clean) %!in% genes.clean),]
copy.clean = copy.clean[-which(rownames(copy.clean) %!in% genes.clean),]
rm(gene.data.co)
rm(gene.data.ex)
```
  dimensions: 16970x28

* assign data formats, make all data nominal (annotation)
```
brain.anno$CCLE_Name = factor(brain.anno$CCLE_Name)
brain.anno$Aliases = factor(brain.anno$Aliases)
brain.anno$Primary.Disease = factor(brain.anno$Primary.Disease)
brain.anno$Subtype.Disease = factor(brain.anno$Subtype.Disease)
brain.anno$Subtype.Gender = factor(brain.anno$Gender)
brain.anno$Subtype.Source = factor(brain.anno$Source)
```

* order rows (gene names) and columns (cell lines) alphabetically
```
ceres.clean = ceres.clean[order(rownames(ceres.clean)), order(colnames(ceres.clean))]
prob.clean = prob.clean[order(rownames(prob.clean)), order(colnames(prob.clean))]
exp.clean = exp.clean[order(rownames(exp.clean)), order(colnames(exp.clean))]
copy.clean = copy.clean[order(rownames(copy.clean)), order(colnames(copy.clean))]
```
* most commonly mutated genes among relevant cells, Gen-Namen aus allen glioblastoma Zelllinien in eine Liste geben (old code)
```									
#common.genes = as.data.frame(table(c(relevant.mutations$`ACH-000036`$Hugo_Symbol,relevant.mutations$`ACH-000040`$Hugo_Symbol, relevant.mutations$`ACH-000075`$Hugo_Symbol, relevant.mutations$`ACH-000098`$Hugo_Symbol, relevant.mutations$`ACH-000208`$Hugo_Symbol, relevant.mutations$`ACH-000215`$Hugo_Symbol, relevant.mutations$`ACH-000137`$Hugo_Symbol, relevant.mutations$`ACH-000152`$Hugo_Symbol, relevant.mutations$`ACH-000591`$Hugo_Symbol, relevant.mutations$`ACH-000673`$Hugo_Symbol, relevant.mutations$`ACH-000231`$Hugo_Symbol, relevant.mutations$`ACH-000244`$Hugo_Symbol, relevant.mutations$`ACH-000760`$Hugo_Symbol, relevant.mutations$`ACH-000368`$Hugo_Symbol, relevant.mutations$`ACH-000376`$Hugo_Symbol, relevant.mutations$`ACH-000445`$Hugo_Symbol, relevant.mutations$`ACH-000464`$Hugo_Symbol, relevant.mutations$`ACH-000469`$Hugo_Symbol, relevant.mutations$`ACH-000479`$Hugo_Symbol, relevant.mutations$`ACH-000570`$Hugo_Symbol, relevant.mutations$`ACH-000571`$Hugo_Symbol, relevant.mutations$`ACH-000623`$Hugo_Symbol, relevant.mutations$`ACH-000631`$Hugo_Symbol, relevant.mutations$`ACH-000738`$Hugo_Symbol, relevant.mutations$`ACH-000756`$Hugo_Symbol, relevant.mutations$`ACH-000819`$Hugo_Symbol, relevant.mutations$`ACH-000128`$Hugo_Symbol, relevant.mutations$`ACH-000887`$Hugo_Symbol)))#
```

* mutation matrices of all glioblastoma cell lines combined into one																		
```
relevant.mutations.combi = do.call(rbind, lapply(which(names(allDepMapData$mutation) %in% cell.lines), function(a) allDepMapData$mutation[[a]]))
common.genes = as.matrix(table(c(relevant.mutations.combi$Hugo_Symbol)))

summary(common.genes)
```

common.genes  
 MT-ND5 :   32  //mitochondrially encoded NADH:ubiquinone oxidoreductase core subunit 5
 TTN    :   26  //titin, large abundant protein of striated muscle
 TP53   :   24  //
 MUC16  :   23  //mucin
 MT-CYB :   18  //mitochondrially encoded cytochrome b
 PTEN   :   17  //phosphatase and tensin homolog
 (Other):10096
 
```
rownames(common.genes) = common.genes$Var1
common.genes$Var1 = NULL
barplot(common.genes, beside = T, names.arg = rownames(common.genes), las = 2)
common.genes.c = subset(common.genes, common.genes$Freq >11)
common.genes.c = common.genes.c[c(6, 11, 10, 7, 4, 9, 3, 2, 1, 8),]
barplot_commongenes <- barplot(common.genes.c, beside = T, names.arg = rownames(common.genes.c), ylab = "Frequency", main = "Most common gene mutations", las = 2)
```

* combine the mutation matrices of all GBM cell lines

```
mutations.all = rbind(relevant.mutations$`ACH-000036`,relevant.mutations$`ACH-000040`, relevant.mutations$`ACH-000075`, relevant.mutations$`ACH-000098`, relevant.mutations$`ACH-000208`, relevant.mutations$`ACH-000215`, relevant.mutations$`ACH-000137`, relevant.mutations$`ACH-000152`, relevant.mutations$`ACH-000591`, relevant.mutations$`ACH-000673`, relevant.mutations$`ACH-000231`, relevant.mutations$`ACH-000244`, relevant.mutations$`ACH-000760`, relevant.mutations$`ACH-000368`, relevant.mutations$`ACH-000376`, relevant.mutations$`ACH-000445`, relevant.mutations$`ACH-000464`, relevant.mutations$`ACH-000469`, relevant.mutations$`ACH-000479`, relevant.mutations$`ACH-000570`, relevant.mutations$`ACH-000571`, relevant.mutations$`ACH-000623`, relevant.mutations$`ACH-000631`, relevant.mutations$`ACH-000738`, relevant.mutations$`ACH-000756`, relevant.mutations$`ACH-000819`, relevant.mutations$`ACH-000128`, relevant.mutations$`ACH-000887`$Hugo_Symbol)
```

* Extract the relevant columns "Hugo_Symbol", "IsDeleterious" and "DepMap_ID" (package dplyr)

```
mutations.all= select(mutations.all, c(2,21, 36))
```

* Determine the GBM cell lines, which contain our four most prominent driver mutations

```
list.cells = subset(mutations.all, mutations.all$Hugo_Symbol %in% rownames(common.genes.c))
unique(list.cells$DepMap_ID)
```

* Extract the cell lines which contain one specific driving mutation or contain one specific mutation **not**, respectively

For MT-ND5:
```
list.mtnd5 = unique(subset(mutations.all, mutations.all$Hugo_Symbol == "MT-ND5"))
cells.mtnd5 = c(list.mtnd5$DepMap_ID)
```

For non-MT-ND5:
```
list.non_mtnd5 = unique(subset(mutations.all, mutations.all$Hugo_Symbol != "MT-ND5"))
cells.non_mtnd5 = c(list.non_mtnd5$DepMap_ID)
```

For MUC16:
```
list.muc16 = unique(subset(mutations.all, mutations.all$Hugo_Symbol == "MUC16"))
cells.muc16 = c(list.muc16$DepMap_ID)
```

For non-MUC16:
```
list.non_muc16 = unique(subset(mutations.all, mutations.all$Hugo_Symbol != "MUC16"))
cells.non_muc16 = c(list.non_muc16$DepMap_ID)
```

For TP53:
```
list.tp53 = unique(subset(mutations.all, mutations.all$Hugo_Symbol == "TP53"))
cells.tp53 = c(list.tp53$DepMap_ID)
```

For non_TP53:
```
list.non_tp53 = unique(subset(mutations.all, mutations.all$Hugo_Symbol != "TP53"))
cells.non_tp53 = c(list.non_tp53$DepMap_ID)
```

For TTN:
```
list.ttn = unique(subset(mutations.all, mutations.all$Hugo_Symbol == "TTN"))
cells.ttn = c(list.ttn$DepMap_ID)
```

For non-TTN:
```
list.non_ttn = unique(subset(mutations.all, mutations.all$Hugo_Symbol != "TTN"))
cells.non_ttn = c(list.non_ttn$DepMap_ID)
```


* Convert expression, CN, CERES and probability matrix for each driving mutation, that it just contains the necessary cell lines

For MT-ND5:
```
exp.mtnd5 = exp.clean[,which(colnames(exp.clean) %in% cells.mtnd5)]
copy.mtnd5 = copy.clean[,which(colnames(copy.clean) %in% cells.mtnd5)]
ceres.mtnd5 = ceres.clean[,which(colnames(ceres.clean) %in% cells.mtnd5)]
prob.mtnd5 = prob.clean[,which(colnames(prob.clean) %in% cells.mtnd5)]
```

For non-MT-ND5:
```
exp.non_mtnd5 = exp.clean[,which(colnames(exp.clean) %in% cells.non_mtnd5)]
copy.non_mtnd5 = copy.clean[,which(colnames(copy.clean) %in% cells.non_mtnd5)]
ceres.non_mtnd5 = ceres.clean[,which(colnames(ceres.clean) %in% cells.non_mtnd5)]
prob.non_mtnd5 = prob.clean[,which(colnames(prob.clean) %in% cells.non_mtnd5)]
```

For MUC16:
```
exp.muc16 = exp.clean[,which(colnames(exp.clean) %in% cells.muc16)]
copy.muc16 = copy.clean[,which(colnames(copy.clean) %in% cells.muc16)]
ceres.muc16 = ceres.clean[,which(colnames(ceres.clean) %in% cells.muc16)]
prob.muc16 = prob.clean[,which(colnames(prob.clean) %in% cells.muc16)]
```

For non-MUC16:
```
exp.non_muc16 = exp.clean[,which(colnames(exp.clean) %in% cells.non_muc16)]
copy.non_muc16 = copy.clean[,which(colnames(copy.clean) %in% cells.non_muc16)]
ceres.non_muc16 = ceres.clean[,which(colnames(ceres.clean) %in% cells.non_muc16)]
prob.non_muc16 = prob.clean[,which(colnames(prob.clean) %in% cells.non_muc16)]
```

For TP53:
```
exp.tp53 = exp.clean[,which(colnames(exp.clean) %in% cells.tp53)]
copy.tp53 = copy.clean[,which(colnames(copy.clean) %in% cells.tp53)]
ceres.tp53 = ceres.clean[,which(colnames(ceres.clean) %in% cells.tp53)]
prob.tp53 = prob.clean[,which(colnames(prob.clean) %in% cells.tp53)]
```

For non-TP53:
```
exp.non_tp53 = exp.clean[,which(colnames(exp.clean) %in% cells.non_tp53)]
copy.non_tp53 = copy.clean[,which(colnames(copy.clean) %in% cells.non_tp53)]
ceres.non_tp53 = ceres.clean[,which(colnames(ceres.clean) %in% cells.non_tp53)]
prob.non_tp53 = prob.clean[,which(colnames(prob.clean) %in% cells.non_tp53)]
```

For TTN:
```
exp.ttn = exp.clean[,which(colnames(exp.clean) %in% cells.ttn)]
copy.ttn = copy.clean[,which(colnames(copy.clean) %in% cells.ttn)]
ceres.ttn = ceres.clean[,which(colnames(ceres.clean) %in% cells.ttn)]
prob.ttn = prob.clean[,which(colnames(prob.clean) %in% cells.ttn)]
```

For non-TTN:
```
exp.non_ttn = exp.clean[,which(colnames(exp.clean) %in% cells.non_ttn)]
copy.non_ttn = copy.clean[,which(colnames(copy.clean) %in% cells.non_ttn)]
ceres.non_ttn = ceres.clean[,which(colnames(ceres.clean) %in% cells.non_ttn)]
prob.non_ttn = prob.clean[,which(colnames(prob.clean) %in% cells.non_ttn)]
```

* Determine mean expression/CN/CERES/probability of all genes over the cell lines, which contain one/or not one specific driving mutation

For MT-ND5:
```
mtnd5.exp.mean = as.matrix(c(rowMeans(exp.mtnd5)))
mtnd5.copy.mean = as.matrix(c(rowMeans(copy.mtnd5)))
mtnd5.ceres.mean = as.matrix(c(rowMeans(ceres.mtnd5)))
mtnd5.prob.mean = as.matrix(c(rowMeans(prob.mtnd5)))
```

For non_MT-ND5:
```
non_mtnd5.exp.mean = as.matrix(c(rowMeans(exp.non_mtnd5)))
non_mtnd5.copy.mean = as.matrix(c(rowMeans(copy.non_mtnd5)))
non_mtnd5.ceres.mean = as.matrix(c(rowMeans(ceres.non_mtnd5)))
non_mtnd5.prob.mean = as.matrix(c(rowMeans(prob.non_mtnd5)))
```

For MUC16:
```
muc16.exp.mean = as.matrix(c(rowMeans(exp.muc16)))
muc16.copy.mean = as.matrix(c(rowMeans(copy.muc16)))
muc16.ceres.mean = as.matrix(c(rowMeans(ceres.muc16)))
muc16.prob.mean = as.matrix(c(rowMeans(prob.muc16)))
```
For non_MUC16:
```
non_muc16.exp.mean = as.matrix(c(rowMeans(exp.non_muc16)))
non_muc16.copy.mean = as.matrix(c(rowMeans(copy.non_muc16)))
non_muc16.ceres.mean = as.matrix(c(rowMeans(ceres.non_muc16)))
non_muc16.prob.mean = as.matrix(c(rowMeans(prob.non_muc16)))
```

For TP53:
```
tp53.exp.mean = as.matrix(c(rowMeans(exp.tp53)))
tp53.copy.mean = as.matrix(c(rowMeans(copy.tp53)))
tp53.ceres.mean = as.matrix(c(rowMeans(ceres.tp53)))
tp53.prob.mean = as.matrix(c(rowMeans(prob.tp53)))
```

For non_TP53:
```
non_tp53.exp.mean = as.matrix(c(rowMeans(exp.non_tp53)))
non_tp53.copy.mean = as.matrix(c(rowMeans(copy.non_tp53)))
non_tp53.ceres.mean = as.matrix(c(rowMeans(ceres.non_tp53)))
non_tp53.prob.mean = as.matrix(c(rowMeans(prob.non_tp53)))
```

For TTN:
```
ttn.exp.mean = as.matrix(c(rowMeans(exp.ttn)))
ttn.copy.mean = as.matrix(c(rowMeans(copy.ttn)))
ttn.ceres.mean = as.matrix(c(rowMeans(ceres.ttn)))
ttn.prob.mean = as.matrix(c(rowMeans(prob.ttn)))
```

For non_TTN:
```
non_ttn.exp.mean = as.matrix(c(rowMeans(exp.non_ttn)))
non_ttn.copy.mean = as.matrix(c(rowMeans(copy.non_ttn)))
non_ttn.ceres.mean = as.matrix(c(rowMeans(ceres.non_ttn)))
non_ttn.prob.mean = as.matrix(c(rowMeans(prob.non_ttn)))
```


# Data Visualization

* Boxplots for whole expression and CN matrix

```{r}
boxplot_expression <- boxplot(exp.clean, ylab ="Expression level", main = "Distribution of expression", par(las =2))
boxplot_CN <- boxplot(copy.clean, ylab = "Copy number", main = "Distribution of copy number", par(las=2))
```
* Boxplot for MT-ND5 mean expression and copy number
 
```
boxplot_mtnd5_exp <- boxplot(mtnd5.exp.mean, ylab = "Expression level", main = "Mean expression of all genes containing MT-ND5 as DM")
boxplot_mtnd5_CN <- boxplot(mtnd5.copy.mean, ylab = "Copy number", main = "Mean copy number of all genes containing MT-ND5 as DM")
```

* Heatmaps of whole mutation/CN/expression/CERES and probability matrices

Expression matrix (few overexpressed, but a lot deleted genes):
```
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(exp.clean, show_rownames = F)
setHook("grid.newpage", NULL, "replace")
grid.text("celllines", y=-0.07, gp=gpar(fontsize=16))
grid.text("genes", x=-0.07, rot=90, gp=gpar(fontsize=16))
```

Copy number matrix:
```
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(copy.clean, show_rownames = F)
setHook("grid.newpage", NULL, "replace")
grid.text("celllines", y=-0.07, gp=gpar(fontsize=16))
grid.text("genes", x=-0.07, rot=90, gp=gpar(fontsize=16))
```

CERES matrix (few genes, which seems to be essential):
```
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(ceres.clean, show_rownames = F)
setHook("grid.newpage", NULL, "replace")
grid.text("celllines", y=-0.07, gp=gpar(fontsize=16))
grid.text("genes", x=-0.07, rot=90, gp=gpar(fontsize=16))
```

Probability matrix (few genes which seem to be essential):
```
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(prob.clean, show_rownames = F)
setHook("grid.newpage", NULL, "replace")
grid.text("celllines", y=-0.07, gp=gpar(fontsize=16))
grid.text("genes", x=-0.07, rot=90, gp=gpar(fontsize=16))
```

# Data reduction

## Additional cleanup beforehand 

* Remove rows containing zero in matrices

```
sum(exp.clean == 0)
exp.clean.w0 <- exp.clean[!(apply(exp.clean, 1, function(y) any(y == 0))),]
sum(exp.clean.w0 == 0)

sum(copy.clean == 0)
copy.clean.w0 <- copy.clean[!(apply(copy.clean, 1, function(y) any(y == 0))),]
sum(copy.clean.w0 == 0)

sum(ceres.clean == 0)
ceres.clean.w0 <- ceres.clean[!(apply(ceres.clean, 1, function(y) any(y == 0))),]
sum(ceres.clean.w0 == 0)

sum(prob.clean == 0)
prob.clean.w0 <- prob.clean[!(apply(prob.clean, 1, function(y) any(y == 0))),]
sum(prob.clean.w0 == 0)
```

* compare gene data availability between data sets
```
dim(copy.clean.w0) == dim(exp.clean.w0)
gene.data.ex = c(rownames(exp.clean.w0))
copy.clean.w0 = copy.clean.w0[-which(rownames(copy.clean.w0) %!in% gene.data.ex),]
genes.clean.w0 = rownames(exp.clean.w0)
ceres.clean.w0 = ceres.clean.w0[-which(rownames(ceres.clean.w0) %!in% genes.clean.w0),]
prob.clean.w0 = prob.clean.w0[-which(rownames(prob.clean.w0) %!in% genes.clean.w0),]
genes.clean.w0 = rownames(prob.clean.w0)
exp.clean.w0 = exp.clean.w0[-which(rownames(exp.clean.w0) %!in% genes.clean.w0),]
copy.clean.w0 = copy.clean.w0[-which(rownames(copy.clean.w0) %!in% genes.clean.w0),]
ceres.clean.w0 = ceres.clean.w0[-which(rownames(ceres.clean.w0) %!in% genes.clean.w0),]
rm(gene.data.ex)
```
  dimensions: 11545x28

## Principal component analysis

* Expression matrix

```
pca_exp <- prcomp(t(exp.clean.w0), scale = TRUE) //matrix must be transposed, because prcomp() expects samples to be rows and genes to be columns; if we dont do it, we get a graph how genes are related to each other
plot(pca_exp$x[,1], pca_exp$x[,2]) //PC1 plotted against PC2
pca_exp_var <- pca_exp$sdev^2 //calculate variation each PC accounts for
pca_exp_var_per <- round(pca_exp_var/sum(pca_exp_var)*100, 1) //convert it into percentage
barplot(pca_exp_var_per, main = "Proportion of variance", xlab = "Principal component", ylab = "Percent variation")

pca_data_exp <- data.frame(Sample = rownames(pca_exp$x), X = pca_exp$x[,1], Y = pca_exp$x[,2]) //create data frame for ggplot
pca_data_exp
pca_exp_plot <- ggplot(data = pca_data_exp, aes(x=X, y=Y, label = Sample)) //tell ggplot which data to use
+ geom_point()   // plot dots, with geom.text() sample names are plotted
+ xlab(paste("PC1 -", pca_exp_var_per[1], "%", sep = ""))  // add labels with percentage
+ ylab(paste("PC2 -", pca_exp_var_per[2], "%", sep="")) 
+ theme_bw()  // make background white, without its grey (idc what to use, just an idea)                                                                            
+ ggtitle("PCA expression") // add title
pca_exp_plot

loading_scores_exp <- pca_exp$rotation[,1]      // determine which genes have largest effect on where the samples are plotted 
gene_scores_exp <- abs(loading_scores_exp)          // make values absolute
gene_score_exp_ranked <- sort(gene_scores_exp, decreasing = TRUE)
top_10_genes_exp <- names(gene_score_exp_ranked[1:10])
top_10_genes_exp
pca_exp$rotation[top_10_genes_exp,1]
```

* Copy number matrix

```
pca_copy <- prcomp(t(copy.clean.w0), scale = TRUE) 
plot(pca_copy$x[,1], pca_copy$x[,2]) 
pca_copy_var <- pca_copy$sdev^2 
pca_copy_var_per <- round(pca_copy_var/sum(pca_copy_var)*100, 1) 
barplot(pca_copy_var_per, main = "Proportion of variance", xlab = "Principal component", ylab = "Percent variation")

pca_data_copy <- data.frame(Sample = rownames(pca_copy$x), X = pca_copy$x[,1], Y = pca_copy$x[,2]) 
pca_data_copy
pca_copy_plot <- ggplot(data = pca_data_copy, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC1 -", pca_copy_var_per[1], "%", sep = ""))  + ylab(paste("PC2 -", pca_copy_var_per[2], "%", sep="")) + theme_bw() + ggtitle("PCA copy number") 
pca_copy_plot

loading_scores_copy <- pca_copy$rotation[,1]      
gene_scores_copy <- abs(loading_scores_copy)          
gene_score_copy_ranked <- sort(gene_scores_copy, decreasing = TRUE)
top_10_genes_copy <- names(gene_score_copy_ranked[1:10])
top_10_genes_copy
pca_copy$rotation[top_10_genes_copy,1]
```

* CERES matrix

```
pca_ceres <- prcomp(t(ceres.clean.w0), scale = TRUE) 
plot(pca_ceres$x[,1], pca_ceres$x[,2]) 
pca_ceres_var <- pca_ceres$sdev^2 
pca_ceres_var_per <- round(pca_ceres_var/sum(pca_ceres_var)*100, 1) 
barplot(pca_ceres_var_per, main = "Proportion of variance", xlab = "Principal component", ylab = "Percent variation")

pca_data_ceres <- data.frame(Sample = rownames(pca_ceres$x), X = pca_ceres$x[,1], Y = pca_ceres$x[,2]) 
pca_data_ceres
pca_ceres_plot <- ggplot(data = pca_data_ceres, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC1 -", pca_ceres_var_per[1], "%", sep = ""))  + ylab(paste("PC2 -", pca_ceres_var_per[2], "%", sep="")) + theme_bw() + ggtitle("PCA CERES score")
pca_ceres_plot

loading_scores_ceres <- pca_ceres$rotation[,1]      
gene_scores_ceres <- abs(loading_scores_ceres)          
gene_score_ceres_ranked <- sort(gene_scores_ceres, decreasing = TRUE)
top_10_genes_ceres <- names(gene_score_ceres_ranked[1:10])
top_10_genes_ceres
pca_ceres$rotation[top_10_genes_ceres,1]
```

* Probability matrix

```
pca_prob <- prcomp(t(prob.clean.w0), scale = TRUE) 
plot(pca_prob$x[,1], pca_prob$x[,2]) 
pca_prob_var <- pca_prob$sdev^2 
pca_prob_var_per <- round(pca_prob_var/sum(pca_prob_var)*100, 1) 
barplot(pca_prob_var_per, main = "Proportion of variance", xlab = "Principal component", ylab = "Percent variation")

pca_data_prob <- data.frame(Sample = rownames(pca_prob$x), X = pca_prob$x[,1], Y = pca_prob$x[,2]) 
pca_data_prob
pca_prob_plot <- ggplot(data = pca_data_prob, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC1 -", pca_prob_var_per[1], "%", sep = ""))  + ylab(paste("PC2 -", pca_prob_var_per[2], "%", sep="")) + theme_bw() + ggtitle("PCA probability score")
pca_prob_plot

loading_scores_prob <- pca_prob$rotation[,1]      
gene_scores_prob <- abs(loading_scores_prob)          
gene_score_prob_ranked <- sort(gene_scores_prob, decreasing = TRUE)
top_10_genes_prob <- names(gene_score_prob_ranked[1:10])
top_10_genes_prob
pca_prob$rotation[top_10_genes_prob,1]
```

## Further steps:

* Interprete the proportion of variances, maybe plot other PC´s against each other
* Interprete plots plus additional plots
* loading scores for more than ten genes (Attention: in these cases the loading scores just refer to the PC1!!!!)
* loading scores for PC2, PC3, .... for each matrix
* Determine which possible 2nd site targets to look at (how many?) for further analyzing
* Make plots more beautiful