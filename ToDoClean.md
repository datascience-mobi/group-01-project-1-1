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
* most commonly mutated genes among relevant cells

```
common.genes = as.data.frame(table(c(relevant.mutations$`ACH-000036`$Hugo_Symbol,relevant.mutations$`ACH-000040`$Hugo_Symbol, relevant.mutations$`ACH-000075`$Hugo_Symbol, relevant.mutations$`ACH-000098`$Hugo_Symbol, relevant.mutations$`ACH-000208`$Hugo_Symbol, relevant.mutations$`ACH-000215`$Hugo_Symbol, relevant.mutations$`ACH-000137`$Hugo_Symbol, relevant.mutations$`ACH-000152`$Hugo_Symbol, relevant.mutations$`ACH-000591`$Hugo_Symbol, relevant.mutations$`ACH-000673`$Hugo_Symbol, relevant.mutations$`ACH-000231`$Hugo_Symbol, relevant.mutations$`ACH-000244`$Hugo_Symbol, relevant.mutations$`ACH-000760`$Hugo_Symbol, relevant.mutations$`ACH-000368`$Hugo_Symbol, relevant.mutations$`ACH-000376`$Hugo_Symbol, relevant.mutations$`ACH-000445`$Hugo_Symbol, relevant.mutations$`ACH-000464`$Hugo_Symbol, relevant.mutations$`ACH-000469`$Hugo_Symbol, relevant.mutations$`ACH-000479`$Hugo_Symbol, relevant.mutations$`ACH-000570`$Hugo_Symbol, relevant.mutations$`ACH-000571`$Hugo_Symbol, relevant.mutations$`ACH-000623`$Hugo_Symbol, relevant.mutations$`ACH-000631`$Hugo_Symbol, relevant.mutations$`ACH-000738`$Hugo_Symbol, relevant.mutations$`ACH-000756`$Hugo_Symbol, relevant.mutations$`ACH-000819`$Hugo_Symbol, relevant.mutations$`ACH-000128`$Hugo_Symbol, relevant.mutations$`ACH-000887`$Hugo_Symbol)))
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
common.genes.c = subset(common.genes, common.genes$Freq >11)