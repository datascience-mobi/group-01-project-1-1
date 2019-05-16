# To Do

## Data overview

## Data cleanup

1.  remove irrelevant columns:

<!-- end list -->

  - data which are not brain
        cancer
    
        allDepMapData$annotation = allDepMapData$annotation[allDepMapData$annotation$Primary.Disease == "Brain Cancer", ]

  - (remove source, gender, aliases?)

<!-- end list -->

2.  remove NA values (data contain a lot of 0 and negative values(esp. copynumber where x >=0))
  copynumber
  allDepMapData$copynumber = allDepMapData$copynumber[-which(apply(allDepMapData$copynumber, 1, function(x){sum(is.na(x))}) > 0),]

<!-- end list -->

  - (remove rows with missing values?)

<!-- end list -->

3.  assign data formats, make all data nominal
    (annotation)

<!-- end list -->

    allDepMapData$annotation$DepMap_ID = factor(allDepMapData$annotation$DepMap_ID)
    allDepMapData$annotation$CCLE_Name = factor(allDepMapData$annotation$CCLE_Name)
    allDepMapData$annotation$Aliases = factor(allDepMapData$annotation$Aliases)
    allDepMapData$annotation$Primary.Disease = factor(allDepMapData$annotation$Primary.Disease)
    allDepMapData$annotation$Subtype.Disease = factor(allDepMapData$annotation$Subtype.Disease)
    allDepMapData$annotation$Subtype.Gender = factor(allDepMapData$annotation$Gender)
    allDepMapData$annotation$Subtype.Source = factor(allDepMapData$annotation$Source)

4.  make cell line identifier be the rownames in
    annotation?

<!-- end list -->

    rownames(allDepMapData$annotation) = allDepMapData$annotation$DepMap_ID
    allDepMapData$annotation = allDepMapData$annotation[, -1]

5.  remove irrelevant cell lines from the other tables as well

<!-- end list -->

    cell.lines = dput(rownames(allDepMapData$annotation))
    allDepMapData$expression = allDepMapData$expression[, -which(colnames(allDepMapData$expression) %!in% cell.lines)]
    allDepMapData$copynumber = allDepMapData$copynumber[, -which(colnames(allDepMapData$copynumber) %!in% cell.lines)]
    relevant.mutations = subset(allDepMapData$mutation, names(allDepMapData$mutation) %in% cell.lines)
    # library "operators" is needed for "%!in%"


summary.anno.sub = summary(allDepMapData$annotation$Subtype.Disease)

compare gene data availability between data sets
dim(allDepMapData$copynumber) == dim(allDepMapData$expression)

gene.data.co = c(rownames(allDepMapData$copynumber))
gene.data.ex = c(rownames(allDepMapData$expression))
allDepMapData$expression = allDepMapData$expression[-which(rownames(allDepMapData$expression) %!in% gene.data.co),]
allDepMapData$copynumber = allDepMapData$copynumber[-which(rownames(allDepMapData$copynumber) %!in% gene.data.ex),]

order rows (gene names) and columns (cell lines) alphabetically

allDepMapData$expression = allDepMapData$expression[order(rownames(allDepMapData$expression)), order(colnames((allDepMapData$expression)))]
allDepMapData$copynumber = allDepMapData$copynumber[order(rownames(allDepMapData$copynumber)), order(colnames((allDepMapData$copynumber)))]

combine expression and copynumber into one data set?