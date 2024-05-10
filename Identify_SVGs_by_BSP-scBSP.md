# Identifying Spatially Variable Genes using Big-Small Patch (BSP) algorithm and sparse matrix implementation as scBSP

## BSP Algorithm

Big-small patch (BSP) is a granularity-guided, data-driven, and parameter-free model for identifying spatial variable genes in 2D and 3D high-throughput spatial transcriptomics data.

![BSP](flowchart.png)

Original implementation of BSP in python is available at https://github.com/juexinwang/BSP


## scBSP - A Fast Tool for Single-Cell Spatially Variable Genes Identifications on Large-Scale Spatially Resolved Transcriptomics Data


This package utilizes a granularity-based dimension-agnostic tool, single-cell big-small patch (scBSP), implementing **sparse matrix** operation and KD-tree/balltree method for distance calculation, for the identification of spatially variable genes on
large-scale data. A corresponding Python library is available at [https://pypi.org/project/scbsp](https://pypi.org/project/scbsp/). The R source code is available at https://github.com/CastleLi/scBSP
. The Python source code is available at https://github.com/YQ-Wang/scBSP


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
HDST provides a subcellular resolution spatial trianscriptomics data, which contains 181,367 spots and 13,243 genes with density XXX. The data can be downloaded at [here](https://github.com/juexinwang/Tutorial_DahShu2024/blob/master/data/CN24_D1_unmodgtf_filtered_red_ut_HDST_final_clean.rds)

    library('scBSP')
    load("./CN24_D1_unmodgtf_filtered_red_ut_HDST_final_clean.rds")

View the expression count matrix as a sparse matrix, each row denotes a gene and each column represents a cell/spot.

```
sp_count[1:5,1:5]

5 x 5 sparse Matrix of class "dgCMatrix"
        1000x100 1000x103 1000x113 1000x114 1000x116
Rcn2           .        .        .        .        .
Mycbp2         .        .        .        .        .
mt-Rnr2        .        .        .        .        .
Mprip          .        .        .        .        .
Mroh1          .        .        .        .        .
```

## Prepare the input: coordinates and gene expression

Extract the coordinates
```
info <- cbind.data.frame(x=as.numeric(sapply(strsplit(colnames(sp_count),split="x"),"[",1)),
                         y=as.numeric(sapply(strsplit(colnames(sp_count),split="x"),"[",2)))
rownames(info)  <- colnames(sp_count)
Coords        <- as.matrix(info)
```

View the coordinates of x and y on 2-D space

```
head(Coords)
            x   y
1000x100 1000 100
1000x103 1000 103
1000x113 1000 113
1000x114 1000 114
1000x116 1000 116
1000x142 1000 142
```

View the dimension of coordinates

```
dim(Coords)
[1] 181367      2
```

## Preprocessing
Removing mitochondrial genes
```
mt_idx      <- grep("mt-",rownames(sp_count))
if(length(mt_idx)!=0){
    sp_count    <- sp_count[-mt_idx,]
}
```

Excluding low expressed genes
```
Filtered_ExpMat <- SpFilter(sp_count)
```

Check the dimension of the input data 
```
dim(Filtered_ExpMat)
[1]  13209 181367
```

## Computing p-values
Calculate Spatially Variable Genes with this sparse matrix using scBSP.
```
P_values <- scBSP(Coords, Filtered_ExpMat)
```

## Recording computation time
Install peakRAM package to record computational resources
```
install.packages("peakRAM")
```
or

```
install.packages("remotes")
remotes::install_github("tpq/peakRAM")
```

Record computational resources, it takes seconds on a laptop.
```
library(peakRAM)
peakRAM(P_values <- scBSP(Coords, Filtered_ExpMat))

Normalizing the expression matrix
Normalizing the coords matrix
Calculating p-values
                            Function_Call Elapsed_Time_sec Total_RAM_Used_MiB Peak_RAM_Used_MiB
1 P_values<-scBSP(Coords,Filtered_ExpMat)            2.877                0.1         505559907
```

## Output the final results

Final results are sorted and filtered if P_values<0.05 
```
results <- P_values[P_values$P_values<0.05,]
results <- results[order(results$P_values),]
head(results)
    GeneNames     P_values
24    Gm42418 0.000000e+00
70      Cmss1 0.000000e+00
74     Camk1d 3.330669e-16
49       Gphn 7.771561e-16
64 CT010467.1 1.975087e-13
43       Cdk8 1.445843e-12

dim(results)
[1] 1829    2
```

# Example 3: SlideSeq V2 data on Mouse Olfactory Bulb

Use SlideSeq V2 data with higher density.

## Install data using SeuratData
```
install.packages("remotes")
remotes::install_github("satijalab/seurat-data")

library(Seurat)
library(SeuratData)
library(SeuratObject)
library(peakRAM)
```

Check the available data
```
AvailableData()
```

Install Slide-seq v2 dataset of mouse hippocampus
```
InstallData("ssHippo")
```

## load data and preprocessing
```
slide.seq <- LoadData("ssHippo")
data_extracted <- scBSP::LoadSpatial(slide.seq)
ExpMatrix_Filtered <- scBSP::SpFilter(data_extracted$ExpMatrix, Threshold = 1)
```

Check the dimensions
```
dim(data_extracted$Coords)
[1] 53173     2

dim(ExpMatrix_Filtered)
[1] 23243 53173
```

Record the computation time
```
peakRAM({ P_values <- scBSP::scBSP(data_extracted$Coords,ExpMatrix_Filtered)})

Normalizing the expression matrix
Normalizing the coords matrix
Calculating p-values
                                                       Function_Call Elapsed_Time_sec Total_RAM_Used_MiB Peak_RAM_Used_MiB
1 {P_values<-scBSP::scBSP(data_extracted$Coords,ExpMatrix_Filtered)}           15.572                  0         456867078
```
We can see it takes more time than HDST data for they have different density.

# Cite
Wang, J., Li, J., Kramer, S.T. et al. Dimension-agnostic and granularity-based spatially variable gene identification using BSP. Nat Commun 14, 7367 (2023). https://doi.org/10.1038/s41467-023-43256-5

# References:
1. https://github.com/juexinwang/BSP/
2. https://github.com/CastleLi/scBSP/
3. Stahl, P. L. et al. Visualization and analysis of gene expression in tissue sections by spatial transcriptomics. Science 353, 78-82, 2016
4. https://github.com/mssanjavickovic/3dst

