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

## extract the coordinates from the rawdata
```
info <- cbind.data.frame(x=as.numeric(sapply(strsplit(colnames(rawcount),split="x"),"[",1)),
                         y=as.numeric(sapply(strsplit(colnames(rawcount),split="x"),"[",2)))
rownames(info) <- colnames(rawcount)
Coords=info[,1:2]
```

# Excluding low expressed genes
```
Filtered_ExpMat <- SpFilter(rawcount)
Filtered_ExpMat <- Matrix::Matrix(Filtered_ExpMat, sparse = TRUE)
```
Check the dimension of the input data, 
```
dim(Coords)
dim(Filtered_ExpMat)
```

# Computing p-values
P_values <- scBSP(Coords, Filtered_ExpMat)
Test the spatially expressed pattern genes.

# Output the final results, i.e., combined p-values, adjusted p-values, etc.

head(P_values)

  GeneNames     P_values
1     GAPDH 2.704573e-09
2      USP4 2.452562e-01
3  MAPKAPK2 4.177737e-03
4     CPEB1 7.656481e-01
5    LANCL2 7.272278e-01
6      MCL1 9.308317e-06

# Example 2:

# Cite
Wang, J., Li, J., Kramer, S.T. et al. Dimension-agnostic and granularity-based spatially variable gene identification using BSP. Nat Commun 14, 7367 (2023). https://doi.org/10.1038/s41467-023-43256-5

# References:
1. https://github.com/juexinwang/BSP/
2. https://github.com/CastleLi/scBSP/
3. Stahl, P. L. et al. Visualization and analysis of gene expression in tissue sections by spatial transcriptomics. Science 353, 78-82, 2016
4. https://github.com/mssanjavickovic/3dst


### System Requirements
* Python 3.7+
* scikit-learn
* numpy
* pandas

### Installation
Tested on Windows 10, Ubuntu 16.04, CentOS 7, MacOS Monterey version 12.4, and MacOS M1 Pro Ventura 13.2.1.

### Quick Start

Place your spatial transcriptomic data as a folder under data/ folder. MOB (2D ST mouse olfactory from Stahl et al.) and 3Dsim are provided as the tutorial usage.

## Example 1: 2D spatial transcriptomics of ST mouse olfactor
```
python BSP.py --datasetName MOB --spaLocFilename Rep11_MOB_spa.csv --expFilename Rep11_MOB_count.csv
```

This step will load location and expression files individually under data/MOB/ folder, and generate MOB_P_values.csv in the project folder, where each row corresponds to each gene, each gene name with the inferred pvalue.

If use beta distribution:
```
python BSP.py --datasetName MOB --spaLocFilename Rep11_MOB_spa.csv --expFilename Rep11_MOB_count.csv --fitDist beta --adjustP
```

User can also output top-quantile genes regardless the p-values using argument ```--empirical```, and manually define quantiles by ```--quantiles```. 
```
python BSP.py --datasetName MOB --spaLocFilename Rep11_MOB_spa.csv --expFilename Rep11_MOB_count.csv --empirical --quantiles 0.05
```

## Example 2: 3D spatial transcriptomics from simulation
```
python BSP.py --inputDir data/3Dsim/  --for3DTag --useDirTag 
```
This step will load all location and expression combined files under data/3Dsim/ folder, and generate Pattern_1_P_values.csv in the project folder, where each row corresponds to each gene, each gene name with the inferred pvalue.

Both 2D and 3D examples should be finished in several seconds. On a MacOS M1 Pro Ventura 13.2.1, example 1 takes ~1 seconds, example 2 takes less than 1 seconds.

### Support data formats
1. Use Coordinates file and Expression file with single study (as example 1)
* Coordinates file: Row as spots, Column as x,y (for 2D), x,y,z (for 3D)
* Expression file: Row as spots, Column as genes

2. Input with single .csv file (as example 2)
* Rows as spots, Columns as 3D Coordinates ("x","y","z") or 2D Coordinates ("x","y")+ Genes

## Usage

```
# Creating coords and expression matrix
Coords <- expand.grid(1:100,1:100, 1:3)
RandFunc <- function(n) floor(10 * stats::rbeta(n, 1, 5))
Raw_Exp <- Matrix::rsparsematrix(nrow = 10^4, ncol = 3*10^4, density = 0.0001, rand.x = RandFunc)

# Excluding low expressed genes
Filtered_ExpMat <- SpFilter(Raw_Exp)
rownames(Filtered_ExpMat) <- paste0("Gene_", 1:nrow(Filtered_ExpMat))

# Computing p-values
P_values <- scBSP(Coords, Filtered_ExpMat)

```
