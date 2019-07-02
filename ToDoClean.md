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

```
library(pheatmap)
library(grid)
```

Expression matrix (few overexpressed, but a lot deleted genes):

```
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(exp.clean, show_rownames = F)
setHook("grid.newpage", NULL, "replace")
grid.text("celllines", y=-0.07, gp=gpar(fontsize=16))
grid.text("genes", x=-0.07, rot=90, gp=gpar(fontsize=16))
```

```
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(exp.clean.w0, show_rownames = F)
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

```
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(copy.clean.w0, show_rownames = F)
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

```
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(ceres.clean.w0, show_rownames = F)
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
```
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(prob.clean.w0, show_rownames = F)
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

### Expression matrix

```
pca_exp <- prcomp(t(exp.clean.w0), scale = TRUE) //matrix must be transposed, because prcomp() expects samples to be rows and genes to be columns; if we dont do it, we get a graph how genes are related to each other
plot(pca_exp$x[,1], pca_exp$x[,2]) //PC1 plotted against PC2
pca_exp_var <- pca_exp$sdev^2 //calculate variation each PC accounts for
pca_exp_var_per <- round(pca_exp_var/sum(pca_exp_var)*100, 1) //convert it into percentage
barplot(pca_exp_var_per, main = "Proportion of variance", xlab = "Principal component", ylab = "Percent variation")
```

* Plotting PC's against each other

*   PC1 against PC2
```{r}
pca_data_exp_12 <- data.frame(Sample = rownames(pca_exp$x), X = pca_exp$x[,1], Y = pca_exp$x[,2]) 
pca_data_exp_12
pca_exp_plot_12 <- ggplot(data = pca_data_exp_12, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC1 -", pca_exp_var_per[1], "%", sep = ""))  + ylab(paste("PC2 -", pca_exp_var_per[2], "%", sep="")) + theme_bw() + ggtitle("PCA expression score")
pca_exp_plot_12
```

*   PC1 against PC3

```{r}
pca_data_exp_13 <- data.frame(Sample = rownames(pca_exp$x), X = pca_exp$x[,1], Y = pca_exp$x[,3]) 
pca_data_exp_13
pca_exp_plot_13 <- ggplot(data = pca_data_exp_13, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC1 -", pca_exp_var_per[1], "%", sep = ""))  + ylab(paste("PC3 -", pca_exp_var_per[3], "%", sep="")) + theme_bw() + ggtitle("PCA expression score")
pca_exp_plot_13
```

* PC1 against PC4:

```{r}
pca_data_exp_14 <- data.frame(Sample = rownames(pca_exp$x), X = pca_exp$x[,1], Y = pca_exp$x[,4]) 
pca_data_exp_14
pca_exp_plot_14 <- ggplot(data = pca_data_exp_14, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC1 -", pca_exp_var_per[1], "%", sep = ""))  + ylab(paste("PC4 -", pca_exp_var_per[4], "%", sep="")) + theme_bw() + ggtitle("PCA expression score")
pca_exp_plot_14
```

* PC1 against PC5:

```{r}
pca_data_exp_15 <- data.frame(Sample = rownames(pca_exp$x), X = pca_exp$x[,1], Y = pca_exp$x[,5]) 
pca_data_exp_15
pca_exp_plot_15 <- ggplot(data = pca_data_exp_15, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC1 -", pca_exp_var_per[1], "%", sep = ""))  + ylab(paste("PC5 -", pca_exp_var_per[5], "%", sep="")) + theme_bw() + ggtitle("PCA expression score")
pca_exp_plot_15
```

* PC2 against PC3:

```{r}
pca_data_exp_23 <- data.frame(Sample = rownames(pca_exp$x), X = pca_exp$x[,2], Y = pca_exp$x[,3]) 
pca_data_exp_23
pca_exp_plot_23 <- ggplot(data = pca_data_exp_23, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC2 -", pca_exp_var_per[2], "%", sep = ""))  + ylab(paste("PC3 -", pca_exp_var_per[3], "%", sep="")) + theme_bw() + ggtitle("PCA expression score")
pca_exp_plot_23
```

* PC2 against PC4:

```{r}
pca_data_exp_24 <- data.frame(Sample = rownames(pca_exp$x), X = pca_exp$x[,2], Y = pca_exp$x[,4]) 
pca_data_exp_24
pca_exp_plot_24 <- ggplot(data = pca_data_exp_24, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC2 -", pca_exp_var_per[2], "%", sep = ""))  + ylab(paste("PC4 -", pca_exp_var_per[4], "%", sep="")) + theme_bw() + ggtitle("PCA expression score")
pca_exp_plot_24
```

* PC2 against PC5:

```{r}
pca_data_exp_25 <- data.frame(Sample = rownames(pca_exp$x), X = pca_exp$x[,2], Y = pca_exp$x[,5]) 
pca_data_exp_25
pca_exp_plot_25 <- ggplot(data = pca_data_exp_25, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC2 -", pca_exp_var_per[2], "%", sep = ""))  + ylab(paste("PC5 -", pca_exp_var_per[5], "%", sep="")) + theme_bw() + ggtitle("PCA expression score")
pca_exp_plot_25
```

* PC3 against PC4:

```{r}
pca_data_exp_34 <- data.frame(Sample = rownames(pca_exp$x), X = pca_exp$x[,3], Y = pca_exp$x[,4]) 
pca_data_exp_34
pca_exp_plot_34 <- ggplot(data = pca_data_exp_34, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC3 -", pca_exp_var_per[3], "%", sep = ""))  + ylab(paste("PC4 -", pca_exp_var_per[4], "%", sep="")) + theme_bw() + ggtitle("PCA expression score")
pca_exp_plot_34
```

* PC3 against PC5:

```{r}
pca_data_exp_35 <- data.frame(Sample = rownames(pca_exp$x), X = pca_exp$x[,3], Y = pca_exp$x[,5]) 
pca_data_exp_35
pca_exp_plot_35 <- ggplot(data = pca_data_exp_35, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC3 -", pca_exp_var_per[3], "%", sep = ""))  + ylab(paste("PC5 -", pca_exp_var_per[5], "%", sep="")) + theme_bw() + ggtitle("PCA expression score")
pca_exp_plot_35
```

* PC4 against PC5:

```{r}
pca_data_exp_45 <- data.frame(Sample = rownames(pca_exp$x), X = pca_exp$x[,4], Y = pca_exp$x[,5]) 
pca_data_exp_45
pca_exp_plot_45 <- ggplot(data = pca_data_exp_45, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC4 -", pca_exp_var_per[4], "%", sep = ""))  + ylab(paste("PC5 -", pca_exp_var_per[5], "%", sep="")) + theme_bw() + ggtitle("PCA expression score")
pca_exp_plot_45
```

* Determine genes with highest variances


* Loading scores PC1

```{r}
loading_scores_exp_1 <- pca_exp$rotation[,1]      
gene_scores_exp_1 <- abs(loading_scores_exp_1)          
gene_score_exp_ranked_1 <- sort(gene_scores_exp_1, decreasing = TRUE)
top_10_genes_exp_1 <- names(gene_score_exp_ranked_1[1:10])
top_10_genes_exp_1
pca_exp$rotation[top_10_genes_exp_1,1]
```

* Loading scores PC2

```{r}
loading_scores_exp_2 <- pca_exp$rotation[,2]      
gene_scores_exp_2 <- abs(loading_scores_exp_2)          
gene_score_exp_ranked_2 <- sort(gene_scores_exp_2, decreasing = TRUE)
top_10_genes_exp_2 <- names(gene_score_exp_ranked_2[1:10])
top_10_genes_exp_2
pca_exp$rotation[top_10_genes_exp_2,2]
```

```{r}
pca_data_exp <- data.frame(Sample = rownames(pca_exp$x), X = pca_exp$x[,1], Y = pca_exp$x[,2]) 
pca_data_exp
pca_exp_plot <- ggplot(data = pca_data_exp, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC1 -", pca_exp_var_per[1], "%", sep = ""))  + ylab(paste("PC2 -", pca_exp_var_per[2], "%", sep="")) + theme_bw() + ggtitle("PCA expression score")
pca_exp_plot

loading_scores_exp <- pca_exp$rotation[,1]      
gene_scores_exp <- abs(loading_scores_exp)          
gene_score_exp_ranked <- sort(gene_scores_exp, decreasing = TRUE)
top_10_genes_exp <- names(gene_score_exp_ranked[1:10])
top_10_genes_exp
pca_exp$rotation[top_10_genes_exp,1]
```

* Loading scores PC3

```{r}
loading_scores_exp_3 <- pca_exp$rotation[,3]      
gene_scores_exp_3 <- abs(loading_scores_exp_3)          
gene_score_exp_ranked_3 <- sort(gene_scores_exp_3, decreasing = TRUE)
top_10_genes_exp_3 <- names(gene_score_exp_ranked_3[1:10])
top_10_genes_exp_3
pca_exp$rotation[top_10_genes_exp_3,3]
```

* Loading scores PC4

```{r}
loading_scores_exp_4 <- pca_exp$rotation[,4]      
gene_scores_exp_4 <- abs(loading_scores_exp_4)          
gene_score_exp_ranked_4 <- sort(gene_scores_exp_4, decreasing = TRUE)
top_10_genes_exp_4 <- names(gene_score_exp_ranked_4[1:10])
top_10_genes_exp_4
pca_exp$rotation[top_10_genes_exp,4]
```

Loading scores PC5

```{r}
loading_scores_exp_5 <- pca_exp$rotation[,5]      
gene_scores_exp_5 <- abs(loading_scores_exp_5)          
gene_score_exp_ranked_5 <- sort(gene_scores_exp_5, decreasing = TRUE)
top_10_genes_exp_5 <- names(gene_score_exp_ranked_5[1:10])
top_10_genes_exp_5
pca_exp$rotation[top_10_genes_exp,5]
```

### Copy number matrix

Barplot with proportional variances

```{r}
pca_copy <- prcomp(t(copy.clean.w0), scale = TRUE)
plot(pca_copy$x[,1], pca_copy$x[,2])
pca_copy_var <- pca_copy$sdev^2
pca_copy_var_per <- round(pca_copy_var/sum(pca_copy_var)*100, 1)
barplot(pca_copy_var_per, main = "Proportion of variance", xlab = "Principal component", ylab = "Percent variation")
```

Plotting PC's against each other

Plot PC1 against PC2:
```{r}
pca_data_copy_12 <- data.frame(sample = rownames(pca_copy$x), X = pca_copy$x[,1], Y = pca_copy$x[,2])
pca_data_copy_12
pca_copy_plot_12 <- ggplot(data = pca_data_copy_12, aes(x=X, y=Y, label = sample)) + geom_point()   + xlab(paste("PC1 -", pca_copy_var_per[1], "%", sep = ""))  + ylab(paste("PC2 -", pca_copy_var_per[2], "%", sep="")) + theme_bw() + ggtitle("PCA copy score")
pca_copy_plot_12
```

Plot PC1 against PC3:

```{r}
pca_data_copy_13 <- data.frame(sample = rownames(pca_copy$x), X = pca_copy$x[,1], Y = pca_copy$x[,3])
pca_data_copy_13
pca_copy_plot_13 <- ggplot(data = pca_data_copy_13, aes(x=X, y=Y, label = sample)) + geom_point()   + xlab(paste("PC1 -", pca_copy_var_per[1], "%", sep = ""))  + ylab(paste("PC3 -", pca_copy_var_per[3], "%", sep="")) + theme_bw() + ggtitle("PCA copy score")
pca_copy_plot_13
```
Plot PC1 against PC4:

```{r}
pca_data_copy_14 <- data.frame(sample = rownames(pca_copy$x), X = pca_copy$x[,1], Y = pca_copy$x[,4])
pca_data_copy_14
pca_copy_plot_14 <- ggplot(data = pca_data_copy_14, aes(x=X, y=Y, label = sample)) + geom_point()   + xlab(paste("PC1 -", pca_copy_var_per[1], "%", sep = ""))  + ylab(paste("PC4 -", pca_copy_var_per[4], "%", sep="")) + theme_bw() + ggtitle("PCA copy score")
pca_copy_plot_14
```
Plot PC1 against PC5:

```{r}
pca_data_copy_15 <- data.frame(sample = rownames(pca_copy$x), X = pca_copy$x[,1], Y = pca_copy$x[,5])
pca_data_copy_15
pca_copy_plot_15 <- ggplot(data = pca_data_copy_15, aes(x=X, y=Y, label = sample)) + geom_point()   + xlab(paste("PC1 -", pca_copy_var_per[1], "%", sep = ""))  + ylab(paste("PC5 -", pca_copy_var_per[5], "%", sep="")) + theme_bw() + ggtitle("PCA copy score")
pca_copy_plot_15
```

Plot PC2 against PC3:

```{r}
pca_data_copy_23 <- data.frame(sample = rownames(pca_copy$x), X = pca_copy$x[,2], Y = pca_copy$x[,3])
pca_data_copy_23
pca_copy_plot_23 <- ggplot(data = pca_data_copy_23, aes(x=X, y=Y, label = sample)) + geom_point()   + xlab(paste("PC2 -", pca_copy_var_per[2], "%", sep = ""))  + ylab(paste("PC3 -", pca_copy_var_per[3], "%", sep="")) + theme_bw() + ggtitle("PCA copy score")
pca_copy_plot_23
```

Plot PC2 against PC4:

```{r}
pca_data_copy_24 <- data.frame(sample = rownames(pca_copy$x), X = pca_copy$x[,2], Y = pca_copy$x[,4])
pca_data_copy_24
pca_copy_plot_24 <- ggplot(data = pca_data_copy_24, aes(x=X, y=Y, label = sample)) + geom_point()   + xlab(paste("PC2 -", pca_copy_var_per[2], "%", sep = ""))  + ylab(paste("PC4 -", pca_copy_var_per[4], "%", sep="")) + theme_bw() + ggtitle("PCA copy score")
pca_copy_plot_24
```

Plot PC2 against PC5:

```{r}
pca_data_copy_25 <- data.frame(sample = rownames(pca_copy$x), X = pca_copy$x[,2], Y = pca_copy$x[,5])
pca_data_copy_25
pca_copy_plot_25 <- ggplot(data = pca_data_copy_25, aes(x=X, y=Y, label = sample)) + geom_point()   + xlab(paste("PC2 -", pca_copy_var_per[2], "%", sep = ""))  + ylab(paste("PC5 -", pca_copy_var_per[5], "%", sep="")) + theme_bw() + ggtitle("PCA copy score")
pca_copy_plot_25
```

Plot PC3 against PC4:

```{r}
pca_data_copy_34 <- data.frame(sample = rownames(pca_copy$x), X = pca_copy$x[,3], Y = pca_copy$x[,4])
pca_data_copy_34
pca_copy_plot_34 <- ggplot(data = pca_data_copy_34, aes(x=X, y=Y, label = sample)) + geom_point()   + xlab(paste("PC3 -", pca_copy_var_per[3], "%", sep = ""))  + ylab(paste("PC4 -", pca_copy_var_per[4], "%", sep="")) + theme_bw() + ggtitle("PCA copy score")
pca_copy_plot_34
```

Plot PC3 against PC5:

```{r}
pca_data_copy_35 <- data.frame(sample = rownames(pca_copy$x), X = pca_copy$x[,3], Y = pca_copy$x[,5])
pca_data_copy_35
pca_copy_plot_35 <- ggplot(data = pca_data_copy_35, aes(x=X, y=Y, label = sample)) + geom_point()   + xlab(paste("PC3 -", pca_copy_var_per[3], "%", sep = ""))  + ylab(paste("PC5 -", pca_copy_var_per[5], "%", sep="")) + theme_bw() + ggtitle("PCA copy score")
pca_copy_plot_35
```

Plot PC4 against PC5: 

```{r}
pca_data_copy_45 <- data.frame(sample = rownames(pca_copy$x), X = pca_copy$x[,4], Y = pca_copy$x[,5])
pca_data_copy_45
pca_copy_plot_45 <- ggplot(data = pca_data_copy_45, aes(x=X, y=Y, label = sample)) + geom_point()   + xlab(paste("PC4 -", pca_copy_var_per[4], "%", sep = ""))  + ylab(paste("PC5 -", pca_copy_var_per[5], "%", sep="")) + theme_bw() + ggtitle("PCA copy score")
pca_copy_plot_45
```
**Determine genes with highest variances

Loading scores PC1

```{r}
loading_scores_copy_1 <- pca_copy$rotation[,1]
gene_scores_copy_1 <- abs(loading_scores_copy_1)
gene_score_copy_ranked_1 <- sort(gene_scores_copy_1, decreasing = TRUE)
top_10_genes_copy_1 <- names(gene_score_copy_ranked_1[1:10])
top_10_genes_copy_1
```

Loading scores PC2

```{r}
loading_scores_copy_2 <- pca_copy$rotation[,2]
gene_scores_copy_2 <- abs(loading_scores_copy_2)
gene_score_copy_ranked_2 <- sort(gene_scores_copy_2, decreasing = TRUE)
top_10_genes_copy_2 <- names(gene_score_copy_ranked_2[1:10])
top_10_genes_copy_2
```

Loading scores PC3

```{r}
loading_scores_copy_3 <- pca_copy$rotation[,3]
gene_scores_copy_3 <- abs(loading_scores_copy_3)
gene_score_copy_ranked_3 <- sort(gene_scores_copy_3, decreasing = TRUE)
top_10_genes_copy_3 <- names(gene_score_copy_ranked_3[1:10])
top_10_genes_copy_3
```

Loading scores PC4

```{r}
loading_scores_copy_4 <- pca_copy$rotation[,4]
gene_scores_copy_4 <- abs(loading_scores_copy_4)
gene_score_copy_ranked_4 <- sort(gene_scores_copy_4, decreasing = TRUE)
top_10_genes_copy_4 <- names(gene_score_copy_ranked_4[1:10])
top_10_genes_copy_4
```

Loading scores PC5

```{r}
loading_scores_copy_5 <- pca_copy$rotation[,5]
gene_scores_copy_5 <- abs(loading_scores_copy_5)
gene_score_copy_ranked_5 <- sort(gene_scores_copy_5, decreasing = TRUE)
top_10_genes_copy_5 <- names(gene_score_copy_ranked_5[1:10])
top_10_genes_copy_5
```


### CERES matrix

```
pca_ceres <- prcomp(t(ceres.clean.w0), scale = TRUE) 
plot(pca_ceres$x[,1], pca_ceres$x[,2]) 
pca_ceres_var <- pca_ceres$sdev^2 
pca_ceres_var_per <- round(pca_ceres_var/sum(pca_ceres_var)*100, 1) 
barplot(pca_ceres_var_per, main = "Proportion of variance", xlab = "Principal component", ylab = "Percent variation")
```
* Plotting PCs against each other

Plot PC1 against PC2:
```
pca_data_ceres_12 <- data.frame(Sample = rownames(pca_ceres$x), X = pca_ceres$x[,1], Y = pca_ceres$x[,2]) 
pca_data_ceres_12
pca_ceres_plot_12 <- ggplot(data = pca_data_ceres_12, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC1 -", pca_ceres_var_per[1], "%", sep = ""))  + ylab(paste("PC2 -", pca_ceres_var_per[2], "%", sep="")) + theme_bw() + ggtitle("PCA ceres score")
pca_ceres_plot_12
```

Plot PC1 against PC3:
```
pca_data_ceres_13 <- data.frame(Sample = rownames(pca_ceres$x), X = pca_ceres$x[,1], Y = pca_ceres$x[,3]) 
pca_data_ceres_13
pca_ceres_plot_13 <- ggplot(data = pca_data_ceres_13, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC1 -", pca_ceres_var_per[1], "%", sep = ""))  + ylab(paste("PC3 -", pca_ceres_var_per[3], "%", sep="")) + theme_bw() + ggtitle("PCA ceres score")
pca_ceres_plot_13
```

Plot PC1 against PC4:
```
pca_data_ceres_14 <- data.frame(Sample = rownames(pca_ceres$x), X = pca_ceres$x[,1], Y = pca_ceres$x[,4]) 
pca_data_ceres_14
pca_ceres_plot_14 <- ggplot(data = pca_data_ceres_14, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC1 -", pca_ceres_var_per[1], "%", sep = ""))  + ylab(paste("PC4 -", pca_ceres_var_per[4], "%", sep="")) + theme_bw() + ggtitle("PCA ceres score")
pca_ceres_plot_14
```

Plot PC1 against PC5:
```
pca_data_ceres_15 <- data.frame(Sample = rownames(pca_ceres$x), X = pca_ceres$x[,1], Y = pca_ceres$x[,5]) 
pca_data_ceres_15
pca_ceres_plot_15 <- ggplot(data = pca_data_ceres_15, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC1 -", pca_ceres_var_per[1], "%", sep = ""))  + ylab(paste("PC5 -", pca_ceres_var_per[5], "%", sep="")) + theme_bw() + ggtitle("PCA ceres score")
pca_ceres_plot_15
```

Plot PC2 against PC3:
```
pca_data_ceres_23 <- data.frame(Sample = rownames(pca_ceres$x), X = pca_ceres$x[,2], Y = pca_ceres$x[,3]) 
pca_data_ceres_23
pca_ceres_plot_23 <- ggplot(data = pca_data_ceres_23, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC2 -", pca_ceres_var_per[2], "%", sep = ""))  + ylab(paste("PC3 -", pca_ceres_var_per[3], "%", sep="")) + theme_bw() + ggtitle("PCA ceres score")
pca_ceres_plot_23
```

Plot PC2 against PC4:
```
pca_data_ceres_24 <- data.frame(Sample = rownames(pca_ceres$x), X = pca_ceres$x[,2], Y = pca_ceres$x[,4]) 
pca_data_ceres_24
pca_ceres_plot_24 <- ggplot(data = pca_data_ceres_24, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC2 -", pca_ceres_var_per[2], "%", sep = ""))  + ylab(paste("PC4 -", pca_ceres_var_per[4], "%", sep="")) + theme_bw() + ggtitle("PCA ceres score")
pca_ceres_plot_24
```

Plot PC2 against PC5:
```
pca_data_ceres_25 <- data.frame(Sample = rownames(pca_ceres$x), X = pca_ceres$x[,2], Y = pca_ceres$x[,5]) 
pca_data_ceres_25
pca_ceres_plot_25 <- ggplot(data = pca_data_ceres_25, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC2 -", pca_ceres_var_per[2], "%", sep = ""))  + ylab(paste("PC5 -", pca_ceres_var_per[5], "%", sep="")) + theme_bw() + ggtitle("PCA ceres score")
pca_ceres_plot_25
```

Plot PC3 against PC4:
```
pca_data_ceres_34 <- data.frame(Sample = rownames(pca_ceres$x), X = pca_ceres$x[,3], Y = pca_ceres$x[,4]) 
pca_data_ceres_34
pca_ceres_plot_34 <- ggplot(data = pca_data_ceres_34, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC3 -", pca_ceres_var_per[3], "%", sep = ""))  + ylab(paste("PC4 -", pca_ceres_var_per[4], "%", sep="")) + theme_bw() + ggtitle("PCA ceres score")
pca_ceres_plot_34
```

Plot PC3 against PC5:
```
pca_data_ceres_35 <- data.frame(Sample = rownames(pca_ceres$x), X = pca_ceres$x[,3], Y = pca_ceres$x[,5]) 
pca_data_ceres_35
pca_ceres_plot_35 <- ggplot(data = pca_data_ceres_35, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC3 -", pca_ceres_var_per[3], "%", sep = ""))  + ylab(paste("PC5 -", pca_ceres_var_per[5], "%", sep="")) + theme_bw() + ggtitle("PCA ceres score")
pca_ceres_plot_35
```

Plot PC4 against PC5:
```
pca_data_ceres_45 <- data.frame(Sample = rownames(pca_ceres$x), X = pca_ceres$x[,4], Y = pca_ceres$x[,5]) 
pca_data_ceres_45
pca_ceres_plot_45 <- ggplot(data = pca_data_ceres_45, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC4 -", pca_ceres_var_per[4], "%", sep = ""))  + ylab(paste("PC5 -", pca_ceres_var_per[5], "%", sep="")) + theme_bw() + ggtitle("PCA ceres score")
pca_ceres_plot_45
```
* Determine genes with highest variances

Loading scores PC1
```
loading_scores_ceres_1 <- pca_ceres$rotation[,1]      
gene_scores_ceres_1 <- abs(loading_scores_ceres_1)          
gene_score_ceres_ranked_1 <- sort(gene_scores_ceres_1, decreasing = TRUE)
top_10_genes_ceres_1 <- names(gene_score_ceres_ranked_1[1:10])
top_10_genes_ceres_1
pca_ceres$rotation[top_10_genes_ceres_1,1]
```

Loading scores PC2
```
loading_scores_ceres_2 <- pca_ceres$rotation[,2]      
gene_scores_ceres_2 <- abs(loading_scores_ceres_2)          
gene_score_ceres_ranked_2 <- sort(gene_scores_ceres_2, decreasing = TRUE)
top_10_genes_ceres_2 <- names(gene_score_ceres_ranked_2[1:10])
top_10_genes_ceres_2
pca_ceres$rotation[top_10_genes_ceres_2,2]
```

Loading scores PC3
```
loading_scores_ceres_3 <- pca_ceres$rotation[,3]      
gene_scores_ceres_3 <- abs(loading_scores_ceres_3)          
gene_score_ceres_ranked_3 <- sort(gene_scores_ceres_3, decreasing = TRUE)
top_10_genes_ceres_3 <- names(gene_score_ceres_ranked_3[1:10])
top_10_genes_ceres_3
pca_ceres$rotation[top_10_genes_ceres_3,3]
```

Loading scores PC4
```
loading_scores_ceres_4 <- pca_ceres$rotation[,4]      
gene_scores_ceres_4 <- abs(loading_scores_ceres_4)          
gene_score_ceres_ranked_4 <- sort(gene_scores_ceres_4, decreasing = TRUE)
top_10_genes_ceres_4 <- names(gene_score_ceres_ranked_4[1:10])
top_10_genes_ceres_4
pca_ceres$rotation[top_10_genes_ceres_4,4]
```
Loading scores PC5
```
loading_scores_ceres_5 <- pca_ceres$rotation[,5]      
gene_scores_ceres_5 <- abs(loading_scores_ceres_5)          
gene_score_ceres_ranked_5 <- sort(gene_scores_ceres_5, decreasing = TRUE)
top_10_genes_ceres_5 <- names(gene_score_ceres_ranked_5[1:10])
top_10_genes_ceres_5
pca_ceres$rotation[top_10_genes_ceres_5,5]
```


### Probability matrix

* Barplot with proportional variances

```
pca_prob <- prcomp(t(prob.clean.w0), scale = TRUE) 
plot(pca_prob$x[,1], pca_prob$x[,2]) 
pca_prob_var <- pca_prob$sdev^2 
pca_prob_var_per <- round(pca_prob_var/sum(pca_prob_var)*100, 1) 
barplot(pca_prob_var_per, main = "Proportion of variance", xlab = "Principal component", ylab = "Percent variation")
```

* Plotting PCs against each other

Plot PC1 against PC2:
```
pca_data_prob_12 <- data.frame(Sample = rownames(pca_prob$x), X = pca_prob$x[,1], Y = pca_prob$x[,2]) 
pca_data_prob_12
pca_prob_plot_12 <- ggplot(data = pca_data_prob_12, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC1 -", pca_prob_var_per[1], "%", sep = ""))  + ylab(paste("PC2 -", pca_prob_var_per[2], "%", sep="")) + theme_bw() + ggtitle("PCA probability score")
pca_prob_plot_12
```

Plot PC1 against PC3:
```
pca_data_prob_13 <- data.frame(Sample = rownames(pca_prob$x), X = pca_prob$x[,1], Y = pca_prob$x[,3]) 
pca_data_prob_13
pca_prob_plot_13 <- ggplot(data = pca_data_prob_13, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC1 -", pca_prob_var_per[1], "%", sep = ""))  + ylab(paste("PC3 -", pca_prob_var_per[3], "%", sep="")) + theme_bw() + ggtitle("PCA probability score")
pca_prob_plot_13
```

Plot PC1 against PC4:
```
pca_data_prob_14 <- data.frame(Sample = rownames(pca_prob$x), X = pca_prob$x[,1], Y = pca_prob$x[,4]) 
pca_data_prob_14
pca_prob_plot_14 <- ggplot(data = pca_data_prob_14, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC1 -", pca_prob_var_per[1], "%", sep = ""))  + ylab(paste("PC4 -", pca_prob_var_per[4], "%", sep="")) + theme_bw() + ggtitle("PCA probability score")
pca_prob_plot_14
```

Plot PC1 against PC5:
```
pca_data_prob_15 <- data.frame(Sample = rownames(pca_prob$x), X = pca_prob$x[,1], Y = pca_prob$x[,5]) 
pca_data_prob_15
pca_prob_plot_15 <- ggplot(data = pca_data_prob_15, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC1 -", pca_prob_var_per[1], "%", sep = ""))  + ylab(paste("PC5 -", pca_prob_var_per[5], "%", sep="")) + theme_bw() + ggtitle("PCA probability score")
pca_prob_plot_15
```
Plot PC2 against PC3:
```
pca_data_prob_23 <- data.frame(Sample = rownames(pca_prob$x), X = pca_prob$x[,2], Y = pca_prob$x[,3]) 
pca_data_prob_23
pca_prob_plot_23 <- ggplot(data = pca_data_prob_23, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC2 -", pca_prob_var_per[2], "%", sep = ""))  + ylab(paste("PC3 -", pca_prob_var_per[3], "%", sep="")) + theme_bw() + ggtitle("PCA probability score")
pca_prob_plot_23
```

Plot PC2 against PC4:
```
pca_data_prob_24 <- data.frame(Sample = rownames(pca_prob$x), X = pca_prob$x[,2], Y = pca_prob$x[,4]) 
pca_data_prob_24
pca_prob_plot_24 <- ggplot(data = pca_data_prob_24, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC2 -", pca_prob_var_per[2], "%", sep = ""))  + ylab(paste("PC4 -", pca_prob_var_per[4], "%", sep="")) + theme_bw() + ggtitle("PCA probability score")
pca_prob_plot_24
```

Plot PC2 against PC5:
```
pca_data_prob_25 <- data.frame(Sample = rownames(pca_prob$x), X = pca_prob$x[,2], Y = pca_prob$x[,5]) 
pca_data_prob_25
pca_prob_plot_25 <- ggplot(data = pca_data_prob_25, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC2 -", pca_prob_var_per[2], "%", sep = ""))  + ylab(paste("PC5 -", pca_prob_var_per[5], "%", sep="")) + theme_bw() + ggtitle("PCA probability score")
pca_prob_plot_25
```

Plot PC3 against PC4:
```
pca_data_prob_34 <- data.frame(Sample = rownames(pca_prob$x), X = pca_prob$x[,3], Y = pca_prob$x[,4]) 
pca_data_prob_34
pca_prob_plot_34 <- ggplot(data = pca_data_prob_34, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC3 -", pca_prob_var_per[3], "%", sep = ""))  + ylab(paste("PC4 -", pca_prob_var_per[4], "%", sep="")) + theme_bw() + ggtitle("PCA probability score")
pca_prob_plot_34
```

Plot PC3 against PC5:
```
pca_data_prob_35 <- data.frame(Sample = rownames(pca_prob$x), X = pca_prob$x[,3], Y = pca_prob$x[,5]) 
pca_data_prob_35
pca_prob_plot_35 <- ggplot(data = pca_data_prob_35, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC3 -", pca_prob_var_per[3], "%", sep = ""))  + ylab(paste("PC5 -", pca_prob_var_per[5], "%", sep="")) + theme_bw() + ggtitle("PCA probability score")
pca_prob_plot_35
```

Plot PC4 against PC5:
```
pca_data_prob_45 <- data.frame(Sample = rownames(pca_prob$x), X = pca_prob$x[,4], Y = pca_prob$x[,5]) 
pca_data_prob_45
pca_prob_plot_45 <- ggplot(data = pca_data_prob_45, aes(x=X, y=Y, label = Sample)) + geom_point()   + xlab(paste("PC4 -", pca_prob_var_per[4], "%", sep = ""))  + ylab(paste("PC5 -", pca_prob_var_per[5], "%", sep="")) + theme_bw() + ggtitle("PCA probability score")
pca_prob_plot_45
```

* Determine genes with highest variances

Loading scores PC1
```
loading_scores_prob_1 <- pca_prob$rotation[,1]      
gene_scores_prob_1 <- abs(loading_scores_prob_1)          
gene_score_prob_ranked_1 <- sort(gene_scores_prob_1, decreasing = TRUE)
top_10_genes_prob_1 <- names(gene_score_prob_ranked_1[1:10])
top_10_genes_prob_1
pca_prob$rotation[top_10_genes_prob_1,1]
```

Loading scores PC2
```
loading_scores_prob_2 <- pca_prob$rotation[,2]      
gene_scores_prob_2 <- abs(loading_scores_prob_2)          
gene_score_prob_ranked_2 <- sort(gene_scores_prob_2, decreasing = TRUE)
top_10_genes_prob_2 <- names(gene_score_prob_ranked_2[1:10])
top_10_genes_prob_2
pca_prob$rotation[top_10_genes_prob_2,2]
```

Loading scores PC3
```
loading_scores_prob_3 <- pca_prob$rotation[,3]      
gene_scores_prob_3 <- abs(loading_scores_prob_3)          
gene_score_prob_ranked_3 <- sort(gene_scores_prob_3, decreasing = TRUE)
top_10_genes_prob_3 <- names(gene_score_prob_ranked_3[1:10])
top_10_genes_prob_3
pca_prob$rotation[top_10_genes_prob_3,3]
```
Loading scores PC4
```
loading_scores_prob_4 <- pca_prob$rotation[,4]      
gene_scores_prob_4 <- abs(loading_scores_prob_4)          
gene_score_prob_ranked_4 <- sort(gene_scores_prob_4, decreasing = TRUE)
top_10_genes_prob_4 <- names(gene_score_prob_ranked_4[1:10])
top_10_genes_prob_4
pca_prob$rotation[top_10_genes_prob_4,4]
```

Loading scores PC5
```
loading_scores_prob_5 <- pca_prob$rotation[,5]      
gene_scores_prob_5 <- abs(loading_scores_prob_5)          
gene_score_prob_ranked_5 <- sort(gene_scores_prob_5, decreasing = TRUE)
top_10_genes_prob_5 <- names(gene_score_prob_ranked_5[1:10])
top_10_genes_prob_5
pca_prob$rotation[top_10_genes_prob_5,5]
```

## Further steps:

<<<<<<< HEAD
=======
* Interprete the proportion of variances, maybe plot other PC??s against each other
* Interprete plots plus additional plots
* loading scores for more than ten genes (Attention: in these cases the loading scores just refer to the PC1!!!!)
* loading scores for PC2, PC3, .... for each matrix
>>>>>>> 83a8f5219e2db8dca17addcd0a62522817a6a870
* Determine which possible 2nd site targets to look at (how many?) for further analyzing
* Make plots more beautiful



## "PCA" plus test (David's idea)

### Data cleanup

* Determine our four most prominent driving mutations

Already done before: TP53, TTN, MT-ND5, MUC16

* Mutation matrices with just the celllines containing one specific DM to determine all other genes, which are mutated too
  Since all other mutated genes are a few thousand, i extracted the ones, which had "is Deleterious = TRUE"

For TP53:
```
mutations.tp53 = rbind(relevant.mutations$`ACH-000036`,relevant.mutations$`ACH-000040`, relevant.mutations$`ACH-000098`, relevant.mutations$`ACH-000208`, relevant.mutations$`ACH-000215`, relevant.mutations$`ACH-000137`, relevant.mutations$`ACH-000152`, relevant.mutations$`ACH-000591`, relevant.mutations$`ACH-000673`, relevant.mutations$`ACH-000231`, relevant.mutations$`ACH-000368`, relevant.mutations$`ACH-000376`, relevant.mutations$`ACH-000445`, relevant.mutations$`ACH-000464`, relevant.mutations$`ACH-000469`, relevant.mutations$`ACH-000570`, relevant.mutations$`ACH-000571`, relevant.mutations$`ACH-000623`, relevant.mutations$`ACH-000673`, relevant.mutations$`ACH-000738`, relevant.mutations$`ACH-000756`, relevant.mutations$`ACH-000819`, relevant.mutations$`ACH-000128`)
mutations.tp53= select(mutations.tp53, c(2,21, 36))//dplyr packgage
mutations.tp53 = subset(mutations.tp53, mutations.tp53$isDeleterious == "TRUE")
```

For TTN:
```
mutations.ttn = rbind(relevant.mutations$`ACH-000040`, relevant.mutations$`ACH-000098`, relevant.mutations$`ACH-000231`, relevant.mutations$`ACH-000445`, relevant.mutations$`ACH-000464`, relevant.mutations$`ACH-000570`, relevant.mutations$`ACH-000623`, relevant.mutations$`ACH-000738`, relevant.mutations$`ACH-000756`, relevant.mutations$`ACH-000760`, relevant.mutations$`ACH-000819`, relevant.mutations$`ACH-000244`, relevant.mutations$`ACH-000479`, relevant.mutations$`ACH-000631`)
mutations.ttn= select(mutations.ttn, c(2,21, 36))
mutations.ttn = subset(mutations.ttn, mutations.ttn$isDeleterious == "TRUE")
```

For MT-ND5:
```
mutations.mtnd5 = rbind(relevant.mutations$`ACH-000040`, relevant.mutations$`ACH-000098`, relevant.mutations$`ACH-000075`, relevant.mutations$`ACH-000231`, relevant.mutations$`ACH-000244`, relevant.mutations$`ACH-000368`, relevant.mutations$`ACH-000445`, relevant.mutations$`ACH-000479`, relevant.mutations$`ACH-000570`, relevant.mutations$`ACH-000571`, relevant.mutations$`ACH-000591`, relevant.mutations$`ACH-000621`, relevant.mutations$`ACH-000738`, relevant.mutations$`ACH-000756`, relevant.mutations$`ACH-000819`)
mutations.mtnd5= select(mutations.mtnd5, c(2,21, 36))
mutations.mtnd5 = subset(mutations.mtnd5, mutations.mtnd5$isDeleterious == "TRUE")
```

For MUC16:
```
mutations.muc16 = rbind(relevant.mutations$`ACH-000036`, relevant.mutations$`ACH-000137`, relevant.mutations$`ACH-000152`, relevant.mutations$`ACH-000208`, relevant.mutations$`ACH-000231`, relevant.mutations$`ACH-000368`, relevant.mutations$`ACH-000376`, relevant.mutations$`ACH-000464`, relevant.mutations$`ACH-000571`, relevant.mutations$`ACH-000738`, relevant.mutations$`ACH-000756`, relevant.mutations$`ACH-000819`)
mutations.muc16= select(mutations.muc16, c(2,21, 36))
mutations.muc16 = subset(mutations.muc16, mutations.muc16$isDeleterious == "TRUE")
```

* Determine the remaining genes and extract them from the CERES matrix to continue working with just them

For TP53:
```
list.tp53.genes = unique(subset(mutations.tp53))
genes.tp53 <- c(list.tp53.genes$Hugo_Symbol)
ceres.tp53.genes <- ceres.clean[which(rownames(ceres.clean) %in% genes.tp53),]
```

For TTN:
```
list.ttn.genes = unique(subset(mutations.ttn))
genes.ttn <- c(list.ttn.genes$Hugo_Symbol)
ceres.ttn.genes <- ceres.clean[which(rownames(ceres.clean) %in% genes.ttn),]
```
For MT-ND5:
```
list.mtnd5.genes = unique(subset(mutations.mtnd5))
genes.mtnd5 <- c(list.mtnd5.genes$Hugo_Symbol)
ceres.mtnd5.genes <- ceres.clean[which(rownames(ceres.clean) %in% genes.mtnd5),]
```

For MUC16:
```
list.muc16.genes = unique(subset(mutations.muc16))
genes.muc16 <- c(list.muc16.genes$Hugo_Symbol)
ceres.muc16.genes <- ceres.clean[which(rownames(ceres.clean) %in% genes.muc16),]
```

### Determination of potential second site targets (SSTs) through correlation tests: Is there a correlation between DM and other mutations in the cell?

#### Check data distribution to determine which correlation test to apply
Parametric correlation tests (e.g. Pearson correlation) need normally distributed data
If data is not normally distributed apply non-parametric tests (e. g. Spearman correlatiion or Wilcoxon Rank sum test)

Methods to determine the type of distribution: qq-plots (visual determination), Shapiro-Wilk-test

* exemplary qq-plot of CERES scores of all genes of cell line ACH-000036

```{r}
qqnorm(ceres.tp53.genes$`ACH-000036`, main="QQ-Plot of CERES scores of all genes of cell line ACH-000036")
qqline(ceres.tp53.genes$`ACH-000036`, datax = FALSE, distribution = qnorm,
        probs = c(0.25, 0.75), qtype = 7)
```

Observation: qq-plot shows curve that deviates from a linear curve shape
Conclusion: CERES scores are not normally distributed
      
* Shapiro-Wilk-test

```{r}
shapiro.test(ceres.tp53.genes$`ACH-000036`) #exemplary for one cell line
```
* Result:

Shapiro-Wilk normality test

data:  ceres.tp53.genes$`ACH-000036`
W = 0.85504, p-value < 2.2e-16

* Interpretation: p-vylue < 0.05 --> CERES scores of all genes of cell line ACH-000036 are not normally distributed

use apply-function to determine distribution of all cell lines

* create a matrix that contains all mutated genes of cell lines

```{r}
ceres.allDM.genes=rbind(ceres.tp53.genes, ceres.ttn.genes, ceres.mtnd5.genes, ceres.muc16.genes)
```


```{r}
lapply(ceres.allDM.genes,shapiro.test)
```

* All cell lines show p-values < 0.05 -> no normal distribution

#### Perform non-parametric tests to determine the correlation of the mutated genes to the driving mutations in order to find potential SSTs

* Pearson correlation

Correlate the CERES scores of every gene across all cell lines belonging to one specific DM to the DM

* Transpose the matrices ceres.DM.genes to ceres.DM.genes_t to make cell lines rows and genes columns

```{r}
t(ceres.tp53.genes)->ceres.tp53.genes_t
t(ceres.ttn.genes)->ceres.ttn.genes_t
t(ceres.muc16.genes)->ceres.muc16.genes_t
t(ceres.mtnd5.genes)->ceres.mtnd5.genes_t

```

* Select columns (CERES scores of one gene across all cell lines) and calculate correlation to DM (one specific column)
  Perform for all genes of a matrix using lapply
  
* exemplary correlation test for one specific gene

```{r}
cor(ceres.tp53.genes_t$ABCA13, ceres.tp53.genes_t$TP53, method = "spearman")
```
#### Problem: cannot define ceres.tp53.genes_t$TP53 as vector; error: $ operator is invalid for atomic vectors -> how to solve? works fine for the untransposed matrix but problems for the transposed one

* solution: generation of a correlation matrix that includes the correlation coefficients between all genes in the matrix
```{r}
cor.ceres.tp53.genes_t<-cor(ceres.tp53.genes_t, method="spearman")
```

####Problem: resulting correlation matrix is very large and most of the data is not relevant for the correlation between the mutations and the DM because this function calculates the correlation coefficients between all genes but be need only one column that includes all correlation coefficients between all mutated genes and the DM

* solution: extract one specific column

???{r}
cor.ceres.tp.53.only<-cor.ceres.tp53.genes_t[1:734,642]
???
* cor.ceres.tp.53.only includes correlation coefficients of all genes to TP53
  for all other matrices analogue
  
* TTN

```{r}
cor.ceres.ttn.genes_t<-cor(ceres.ttn.genes_t, method="spearman")
cor.ceres.ttn.only<-cor.ceres.ttn.genes_t[1:585,520]
```

* MUC16

```{r}
cor.ceres.muc16.genes_t<-cor(ceres.muc16.genes_t, method="spearman")
cor.ceres.muc16.only<-cor.ceres.muc16.genes_t[1:478,251]
```
* MTND5

```{r}
cor.ceres.mtnd5.genes_t<-cor(ceres.mtnd5.genes_t, method="spearman")
cor.ceres.mtnd5.only<-cor.ceres.mtnd5.genes_t[1:496,???] \\`MTND5 does not appear in the matrix???
```

#### Correlation coefficients for all genes and DMs were calculated - next step: determination of p-values to get correlation matrices with significance levels (statistical statements)

* Hmisc package needs to be installed first

```{r}
install.packages("Hmisc")
library(Hmisc)
```

* function rcorr(x, type = c("pearson","spearman")) can be used to compute the significance levels for Spearman and Pearson correlation


* calculate significance levels for TP53

```{r}
sig.cor.ceres.tp53.genes_t<-rcorr(as.matrix(ceres.tp53.genes_t), type="spearman")
View(sig.cor.ceres.tp53.genes_t$r)
```

* exact meaning of r, n, P?

* extract relevant column:

```{r}
sig.cor.ceres.tp53.only<-sig.cor.ceres.tp53.genes_t$r[1:734, 642]

```

* repeat for remaining DMs

* TTN

```{r}
sig.cor.ceres.ttn.genes_t<-rcorr(as.matrix(ceres.ttn.genes_t), type="spearman")
View(sig.cor.ceres.ttn.genes_t$r)
sig.cor.ceres.ttn.only<-sig.cor.ceres.ttn.genes_t$r[1:585,520]
```

* MUC16

```{r}
sig.cor.ceres.muc16.genes_t<-rcorr(as.matrix(ceres.muc16.genes_t), type="spearman")
View(sig.cor.ceres.muc16.genes_t$r)
sig.cor.ceres.muc16.only<-sig.cor.ceres.muc16.genes_t$r[1:478,251]
```

* MTND5

```{r}
sig.cor.ceres.mtnd5.genes_t<-rcorr(as.matrix(ceres.mtnd5.genes_t), type="spearman")
View(sig.cor.ceres.mtnd5.genes_t$r)
sig.cor.ceres.mtnd5.only<-sig.cor.ceres.mtnd5.genes_t$r[1:496,???] \\see above
```
## Multiple linear regression

* Calculate mean of expression, copy number, ceres score and probability throughput all cell lines and combine all values in one matrix for multiple linear regression. At first using only matrices without 0 values.

```
mlr.mat = as.data.frame(cbind(rowMeans(exp.clean.w0), rowMeans(copy.clean.w0), rowMeans(ceres.clean.w0), rowMeans(prob.clean.w0)))
colnames(mlr.mat) = c("expression", "copynumber", "ceres", "probability")
```
* Perform multiple linear regression and look at results 
( R^2 values of copynumber and expression don't seem right.)
```
summary(lm(expression ~ ., data = mlr.mat)) //r2=0.1563
summary(lm(copynumber ~ ., data = mlr.mat)) //r2=0.003437
summary(lm(ceres ~ ., data = mlr.mat)) //r2=0.9352
summary(lm(probability ~ ., data = mlr.mat)) //r2=0.9358
```
* R^2 values for copynumber and expression of first multiple linear regression don't seem right. May have to do with matrices not containing 0 values?

* Do step one again, this time with matrices containing 0 values.
```
mlr.mat.0 = as.data.frame(cbind(rowMeans(exp.clean), rowMeans(copy.clean), rowMeans(ceres.clean), rowMeans(prob.clean)))
colnames(mlr.mat.0) = c("expression", "copynumber", "ceres", "probability")
```
* Perform multiple linear regression again. 
```
summary(lm(expression ~ ., data = mlr.mat.0)) //r2=0.1934
summary(lm(copynumber ~ ., data = mlr.mat.0)) //r2=0.00411
summary(lm(ceres ~ ., data = mlr.mat.0)) //r2=0.9243
summary(lm(probability ~ ., data = mlr.mat.0)) //r2=0.9249
```
* R^2 values for copynumber and expression got better, but not significantly. Values for ceres and probability got worse. 

# Testing our multiple linear regression model, (package "caTools" required)

#For the expression model:

*Splitting mlr.mat into the training set and test set
```
set.seed(123) #initialize the random numbers
split.exp = sample.split(mlr.mat$expression, SplitRatio = 0.8) #split the mlr.mat into 4/5 Training and 1/5 Testing expression values
split.exp = as.data.frame(split.exp) #converting split.exp to a data frame for further analysis
training_set_exp = subset(mlr.mat, split.exp == TRUE) #use labels to get training data
test_set_exp = subset(mlr.mat, split.exp == FALSE) #dim(test_set_exp) will give you 2309 --> 11545/5*1 = 2309 --> train/test split worked
``` 
Feature Scaling can be performed
```
training_set_exp = scale(training_set_exp)
test_set_exp = scale(test_set_exp)
```

*Fitting multiple linear regression to the Training set
```
exp.regressor = lm(expression ~ ., data = training_set_exp)
```

*Predicting the test set results

```
exp_pred = predict(exp.regressor, newdata = test_set_exp) #predict expression based on your testing data (data taht the model did NEVER see and highly useful to evaluate the performance of the model)
test_set_exp$Prediction = exp_pred #add predictions to mlr.mat 
test_set_exp #now a comparison of the Predictions (last column) with the real values for the expression (1st column) is possible 
```
*Do the same for other models (copynumber, ceres, probability)

#For copy number model
```
set.seed(123)
split.copy = sample.split(mlr.mat$copynumber, SplitRatio = 0.8)
split.copy = as.data.frame(split.copy)
training.set.copy = subset(mlr.mat, split.copy == TRUE)
test.set.copy = subset(mlr.mat, split.copy == FALSE)
dim(test.set.copy) 
```
dimension was 2021, which was not expected and does not equal one fifth of mlr.mat
--> continue with other models, come back to this later.

#For ceres model

*Splitting data into training and test sets
```
set.seed(123)
split.ceres = sample.split(mlr.mat$ceres, SplitRatio = 0.8)
training.set.ceres = subset(mlr.mat, split.ceres == TRUE)
test.set.ceres = subset(mlr.mat, split.ceres == FALSE)
dim(test.set.ceres) #2309 --> training/test split successful
```

*Scaling
```
training.set.ceres = scale(training.set.ceres)
test.set.ceres = scale(test.set.ceres)
```

*Fitting multiple linear regression to the Training set
```
ceres.regressor = lm(ceres ~ ., data = as.data.frame(training.set.ceres))
```

*Predicting test set results
```
ceres_pred = predict(ceres.regressor, newdata = test.set.ceres)
test.set.ceres$Prediction = ceres_pred
View(test.set.ceres)
```

#For probability model

*Splitting data into training/test sets
```
set.seed(123)
split.prob = sample.split(mlr.mat$probability, SplitRatio = 0.8)
training.set.prob = subset(mlr.mat, split.prob ==TRUE)
test.set.prob = subset(mlr.mat, split.prob == FALSE)
```

*Fitting multiple linear regression to the training set
```
prob.regressor = lm(probability ~., data = as.data.frame(training.set.prob))
```

*Predicting test set results
```
prob_pred = predict(prob.regressor, newdata = test.set.prob)
test.set.prob$Prediction = prob_pred
```

# Visualizing real values against test values

*For expression model
```
plot(density(test_set_exp$expression), col = "red")
lines(density(test_set_exp$Prediction), col = "blue")
```
*For ceres model
```
plot(density(test.set.ceres$ceres), col = "red")
lines(density(test.set.ceres$Prediction))
```

*For probability model
```
plot(density(test.set.prob$probability), col = "red")
lines(density(test.set.prob$Prediction))
```

#### Follow up: Interpretation of the p-values, maybe Wilcoxon Rank Sum test

### Determination of potential second site targets (SSTs) through correlation tests: Is there a correlation between DM and other mutations in the cell?

#### Check data distribution to determine which correlation test to apply
* Parametric correlation tests (e.g. Pearson correlation) need normally distributed data
* If data is not normally distributed apply non-parametric tests (e. g. Spearman correlation or Wilcoxon Rank sum test)

* Methods to determine the type of distribution: qq-plots (visual determination), Shapiro-Wilk-test

* Exemplary qq-plot of CERES scores of all genes of cell line ACH-000036

```{r}
qqnorm(ceres.tp53.genes$`ACH-000036`, main="QQ-Plot of CERES scores of all genes of cell line ACH-000036")
qqline(ceres.tp53.genes$`ACH-000036`, datax = FALSE, distribution = qnorm,
        probs = c(0.25, 0.75), qtype = 7)
```

* Observation: qq-plot shows curve that deviates from a linear curve shape
* Conclusion: CERES scores are not normally distributed
      
* Shapiro-Wilk-test

```{r}
shapiro.test(ceres.tp53.genes$`ACH-000036`) #exemplary for one cell line
```
* Result:

Shapiro-Wilk normality test

data:  ceres.tp53.genes$`ACH-000036`
W = 0.85504, p-value < 2.2e-16

* Interpretation: p-value < 0.05 --> CERES scores of all genes of cell line ACH-000036 are not normally distributed

use apply-function to determine distribution of all cell lines

* create a matrix that contains all mutated genes of cell lines

```{r}
ceres.allDM.genes=rbind(ceres.tp53.genes, ceres.ttn.genes, ceres.mtnd5.genes, ceres.muc16.genes)
```


```{r}
lapply(ceres.allDM.genes,shapiro.test)
```

* Observation: All cell lines show p-values < 0.05
* Conclusion: No normal distribution

#### Perform non-parametric tests to determine the correlation of the mutated genes to the driving mutations in order to find potential SSTs

* Spearman correlation

* Correlate the CERES scores of every gene across all cell lines belonging to one specific DM to the DM

* Transpose the matrices ceres.DM.genes to ceres.DM.genes_t to make cell lines rows and genes columns

```{r}
t(ceres.tp53.genes)->ceres.tp53.genes_t
t(ceres.ttn.genes)->ceres.ttn.genes_t
t(ceres.muc16.genes)->ceres.muc16.genes_t
t(ceres.mtnd5.genes)->ceres.mtnd5.genes_t

```

* Select columns (CERES scores of one gene across all cell lines) and calculate correlation to DM (one specific column)
  Perform for all genes of a matrix using lapply
  
* exemplary correlation test for one specific gene

```{r}
cor(ceres.tp53.genes_t$ABCA13, ceres.tp53.genes_t$TP53, method = "spearman")
```
#### Problem: Cannot define ceres.tp53.genes_t$TP53 as vector; error: $ operator is invalid for atomic vectors -> how to solve? works fine for the untransposed matrix but problems for the transposed one


* Solution: Generation of a correlation matrix that includes the correlation coefficients between all genes in the matrix
```{r}
cor.ceres.tp53.genes_t<-cor(ceres.tp53.genes_t, method="spearman")
```

####Problem: resulting correlation matrix is very large and most of the data is not relevant for the correlation between the mutations and the DM because this function calculates the correlation coefficients between all genes but be need only one column that includes all correlation coefficients between all mutated genes and the DM

* Solution: Extract one specific column

???{r}
cor.ceres.tp.53.only<-cor.ceres.tp53.genes_t[1:734,642]
???
* cor.ceres.tp.53.only includes correlation coefficients of all genes to TP53
  for all other matrices analogue
  
* TTN

```{r}
cor.ceres.ttn.genes_t<-cor(ceres.ttn.genes_t, method="spearman")
cor.ceres.ttn.only<-cor.ceres.ttn.genes_t[1:585,520]
```

* MUC16

```{r}
cor.ceres.muc16.genes_t<-cor(ceres.muc16.genes_t, method="spearman")
cor.ceres.muc16.only<-cor.ceres.muc16.genes_t[1:478,251]
```
* MTND5

```{r}
cor.ceres.mtnd5.genes_t<-cor(ceres.mtnd5.genes_t, method="spearman")
cor.ceres.mtnd5.only<-cor.ceres.mtnd5.genes_t[1:496,???] \\`MTND5 does not appear in the matrix???
```

#### Correlation coefficients for all genes and DMs were calculated - next step: determination of p-values to get correlation matrices with significance levels (statistical statements)

* Hmisc package needs to be installed first

```{r}
install.packages("Hmisc")
library(Hmisc)
```

<<<<<<< HEAD
* Function rcorr(x, type = c("pearson","spearman")) can be used to compute the significance levels for Spearman and Pearson correlation
=======
* function rcorr(x, type = c("pearson","spearman")) can be used to compute the significance levels for Spearman and Pearson correlation

>>>>>>> acde25c6f0558dd5272804d55d119e905257a012

* Calculate significance levels for TP53

```{r}
sig.cor.ceres.tp53.genes_t<-rcorr(as.matrix(ceres.tp53.genes_t), type="spearman")
View(sig.cor.ceres.tp53.genes_t$r)
```

* Exact meaning of r, n, P?

* Extract relevant column:

```{r}
sig.cor.ceres.tp53.only<-sig.cor.ceres.tp53.genes_t$r[1:734, 642]

```

* Repeat for remaining DMs

* TTN

```{r}
sig.cor.ceres.ttn.genes_t<-rcorr(as.matrix(ceres.ttn.genes_t), type="spearman")
View(sig.cor.ceres.ttn.genes_t$r)
sig.cor.ceres.ttn.only<-sig.cor.ceres.ttn.genes_t$r[1:585,520]
```

* MUC16

```{r}
sig.cor.ceres.muc16.genes_t<-rcorr(as.matrix(ceres.muc16.genes_t), type="spearman")
View(sig.cor.ceres.muc16.genes_t$r)
sig.cor.ceres.muc16.only<-sig.cor.ceres.muc16.genes_t$r[1:478,251]
```

* MTND5

```{r}
sig.cor.ceres.mtnd5.genes_t<-rcorr(as.matrix(ceres.mtnd5.genes_t), type="spearman")
View(sig.cor.ceres.mtnd5.genes_t$r)
sig.cor.ceres.mtnd5.only<-sig.cor.ceres.mtnd5.genes_t$r[1:496,???] \\see above
```
#### Follow up: Interpretation of the p-values, maybe Wilcoxon Rank Sum test if needed --> asked David, Spearman correlation is sufficient

* P-value: is the probability of obtaining a larger (one-sided upper tail), a smaller (one-sided lower tail) or more extreme value (two-sided or two-tailed) value of the test statistics if H0 is valid --> low p-value: H0 is invalid with a high probability
* We need a significance level alpha to evaluate our p-values in order to identify potential SSTs --> p < ??: H0 is invalid and can be rejected, the observed effect is significant, H1 is statistically proven
* ?? = 0.05 is often used as a standard value --> we will apply the same value to evaluate our correlation results and to evaluate H0 and H1

#### Formulation of our H0 and H1 hypothesis:

* H0: The correlation of a mutated gene and the DM is not significant --> no potential SST
* H1: The correlation of a mutated gene and the DM is significant --> potential SST

* Define a threshold for the p-values to identify all the genes that are correlated to the DMs with a high probability

* TP53

```{r}
sig.cor.ceres.tp53.only[abs(sig.cor.ceres.tp53.only) > 0.05]<- NA \\ only significant p-values remain --> all genes that are NA are no potential SSTs; all genes where H0 cannot be rejected are defined as NA
sst.tp53<-na.omit(sig.cor.ceres.tp53.only) \\ returns a matrix that contains all potential SSTs that we can compare to literature afterwards
```

* TTN

```{r}
sig.cor.ceres.ttn.only[abs(sig.cor.ceres.ttn.only) > 0.05]<- NA
sst.ttn<-na.omit(sig.cor.ceres.ttn.only)
```

* MUC16

```{r}
sig.cor.ceres.muc16.only[abs(sig.cor.ceres.muc16.only) > 0.05] <- NA
SST.muc16<-na.omit(sig.cor.ceres.muc16.only)
```

* MTND-5

```{r}
sig.cor.ceres.mtnd5.only[abs(sig.cor.ceres.mtnd5.only) > 0.05] <- NA
SST.mtnd5<-na.omit(sig.cor.ceres.mtnd5.only)
```

#### Follow up: Comparison to literature

* TP53: literature: top plausible SSL pairs with significant impact on overall survival in GBM: 
    SLC1A5 --> our data:
    PLK1 --> our data:
    MDM2: all genes except MDM2 are overexpressed if TP53 is mutated, MDM2 is overexpressed by TP53 wt ??? mutually exclusive --> our data:
    DBF4: overexpressed in context of GBM
    TCP1: overexpressed in context of GBM when TP53 is mutated --> our data: 
    IDH1 mutations are significantly correlated with the drastic overexpression of several known GBM survival genes --> our data: 
    AKT3 overexpression was found in a significant fraction of breast and prostate cancers and has been reported as a possible oncogene and a     potential glioma survival gene; drastic overexpression of AKT3 is exclusively and significantly associated with IDH1 mutation
    
#### Further research and interpretation needed \\ will do it during the following week
