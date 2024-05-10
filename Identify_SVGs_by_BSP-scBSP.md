# Identifying Spatially Variable Genes using Big-Small Patch (BSP) algorithm and sparse matrix implementation as scBSP

## BSP Algorithm

Big-small patch (BSP) is a granularity-guided, data-driven, and parameter-free model for identifying spatial variable genes in 2D and 3D high-throughput spatial transcriptomics data.

![BSP](flowchart.png)

Original implementation of BSP in python is available at https://github.com/juexinwang/BSP


## scBSP - A Fast Tool for Single-Cell Spatially Variable Genes Identifications on Large-Scale Spatially Resolved Transcriptomics Data


This package utilizes a granularity-based dimension-agnostic tool, single-cell big-small patch (scBSP), implementing **sparse matrix** operation and KD-tree/balltree method for distance calculation, for the identification of spatially variable genes on
large-scale data. A corresponding Python library is available at [https://pypi.org/project/scbsp](https://pypi.org/project/scbsp/).


# System Requirement
Tested on MacOS Sonoma version 14.4.1 with R version 4.3.1, XXX 

# Installation
This package can be installed on R CRAN
```
install.packages("scBSP")
```

# Example 1: Breast Cancer Data
Load the scBSP package and Breast cancer data set, which can be downloaded [here](https://github.com/juexinwang/Tutorial_DahShu2024/blob/master/data/Layer2_BC_Count.rds).

    library('scBSP')
    load("./Layer2_BC_Count.rds")
     
View the expression count matrix rawcount, each row denotes a gene and each column represents a cell/spot.

```
rawcount[1:5,1:5]

    17.907x4.967   18.965x5.003   18.954x5.995    17.846x5.993 20.016x6.019
GAPDH   1   7   5   1   2
USP4    1   0   0   0   0
MAPKAPK2    1   1   0   0   1
CPEB1   0   0   0   0   0
LANCL2  0   0   0   0   0
```

## Prepare the input: coordinates and gene expression

Extract the coordinates from the rawdata
```
info <- cbind.data.frame(x=as.numeric(sapply(strsplit(colnames(rawcount),split="x"),"[",1)),
                         y=as.numeric(sapply(strsplit(colnames(rawcount),split="x"),"[",2)))
rownames(info) <- colnames(rawcount)
Coords=info[,1:2]
```

View the coordinates of x and y on 2-D space

```
head(Coords)
                  x     y
17.907x4.967 17.907 4.967
18.965x5.003 18.965 5.003
18.954x5.995 18.954 5.995
17.846x5.993 17.846 5.993
20.016x6.019 20.016 6.019
20.889x6.956 20.889 6.956
```

View the dimension of coordinates

```
dim(Coords)
[1] 251   2
```

## Preprocessing
Excluding low expressed genes
```
Filtered_ExpMat <- SpFilter(rawcount)
```

(Optional) For this data, we need to explicitly convert the data format to Matrix
```
Filtered_ExpMat <- Matrix::Matrix(Filtered_ExpMat, sparse = TRUE)
```

Check the dimension of the input data 
```
dim(Filtered_ExpMat)
[1] 11769   251
```

## Computing p-values
The inputs are the coordinates and the expresson matrix, scBSP calculates the p-values
```
P_values <- scBSP(Coords, Filtered_ExpMat)
```

## Output the final results
```
head(P_values)
  GeneNames     P_values
1     GAPDH 2.704573e-09
2      USP4 2.452562e-01
3  MAPKAPK2 4.177737e-03
4     CPEB1 7.656481e-01
5    LANCL2 7.272278e-01
6      MCL1 9.308317e-06
```

Final results are sorted and filtered if P_values<0.05 
```
results <- P_values[P_values$P_values<0.05,]
results <- results[order(results$P_values),]
head(results)
     GeneNames     P_values
326     SPINT2 0.000000e+00
261      POSTN 6.439294e-15
2013    COL1A1 6.550316e-15
1347       B2M 7.660539e-15
2370       FN1 1.521006e-14
1072    EEF1A1 1.356693e-13

dim(results)
[1] 1765    2
```




# Example 2: HDST data of mouse hippocampus

HDST provides a subcellular resolution spatial trianscriptomics data, which contains XX spots and XX genes with density XXX. The data can be downloaded at [here](https://github.com/juexinwang/Tutorial_DahShu2024/blob/master/data/CN24_D1_unmodgtf_filtered_red_ut_HDST_final_clean.rds)




# Example 2: SlideSeq V2 data on Mouse Olfactory Bulb

install.packages("remotes")
remotes::install_github("satijalab/seurat-data")

library(Seurat)
library(SeuratData)
library(SeuratObject)
library(peakRAM)

AvailableData()

InstallData("ssHippo")
slide.seq <- LoadData("ssHippo")
data_extracted <- scBSP::LoadSpatial(slide.seq)
ExpMatrix_Filtered <- scBSP::SpFilter(data_extracted$ExpMatrix, Threshold = 1)
P_values <- scBSP::scBSP(data_extracted$Coords, ExpMatrix_Filtered)

scBSP_mem <- peakRAM({
  P_values <- scBSP::scBSP(data_extracted$Coords, ExpMatrix_Filtered)
})

# Cite
Wang, J., Li, J., Kramer, S.T. et al. Dimension-agnostic and granularity-based spatially variable gene identification using BSP. Nat Commun 14, 7367 (2023). https://doi.org/10.1038/s41467-023-43256-5

# References:
1. https://github.com/juexinwang/BSP/
2. https://github.com/CastleLi/scBSP/
3. Stahl, P. L. et al. Visualization and analysis of gene expression in tissue sections by spatial transcriptomics. Science 353, 78-82, 2016
4. https://github.com/mssanjavickovic/3dst

