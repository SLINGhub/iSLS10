## iSLS10 Workshop Lipidomics Data Analysis

# Sphingolipids and T2DM Risk

## Manage, analyse and interpret lipidomics data using R

Summary

This repository contains datasets and R codes used for 9th International Singapore Lipid Symposium [iSLS10](https://sling.sg/news-events/isls/) workshop held Wednesday, March 9, 2022. We will be using the data published in *Chew et al., 2019, JCI Insight 5(13)* [10.1172/jci.insight.126925](https://doi.org/10.1172/jci.insight.126925) as an example dataset for this workshop.

In the first part, we will import and clean a lipidomics dataset and de-identified metadata into R, look into quality control aspects, convert lipid names to current nomenclature, and re-structure the data for subsequent statistical analysis covered in the second part. We will also explore how to efficiently conduct descriptive statistics annd make plots for each analyte in the dataset using R.

In the second part, we will inspect the overall data trends in both sample meta data and lipidomics data via visualization. We will study the correlation structure between sphingolipids and other competing risk factors collected at baseline, and test associations of individual lipids with the risk of DM incidence, using logistic regression (binary outcome) and Cox regression analysis (time-to-event analysis). We will also explore options for multi-variable modelling, such as multi-variable Cox regression with group structured LASSO penalty and popular machine learning methods (e.g. random forest or boosting).

## Download and run the R scripts

Download the R Project containing the scripts and data used in this workshop from this repository.

![](images/Screenshot%202022-03-05%20141709.png){width="312"}

Following R packages and external software should be installed on your R system before running the scripts (see also comments in the R scripts of Part 1 and 2):

-   Part 1:

    -   CRAN packages: `here`, `tidyverse`, `ggrepel`

    -   Optional: `rgoslin` (used to convert lipid nomenclatures). Only available via github (<https://github.com/lifs-tools/rgoslin>). Note that on Windows, you will need an installation of [rtools](https://cran.r-project.org/bin/windows/Rtools/), and on macOS of [XCode](https://apps.apple.com/sg/app/xcode/id497799835?mt=12) in order to install `rgoslin.`

-   Part 2:

    -   CRAN packages: `scales`, `gplots`, `survival, glmnet, grpreg, randomForest, huge`

    -   Bioconductor packages: [`mixOmics`](https://bioconductor.org/packages/release/bioc/html/mixOmics.html), [`qvalue`](https://bioconductor.org/packages/release/bioc/html/qvalue.html) , [`impute`](https://www.bioconductor.org/packages/release/bioc/html/impute.html)

    -   Cytoscape (<https://cytoscape.org/download.html>)

## Authors

-   Bo Burla - [Singapore Lipidomics Incubator \@ National University of Singapore](https://sling.sg)

-   Hyungwon Choi - [Computational & Statistical Systems Biology Laboratory \@ National University of Singapore](https://www.cssblab.org)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

-   Wee Siong Chew and Shanshan Ji
