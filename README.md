RNAAgeCalc
===========

`RNAAgeCalc` is a multi-tissue transcriptionals age calculator. 

The `RNAAgeCalc` paper can be found [here](https://doi.org/10.1371/journal.pone.0237006).   
The Bioconductor package can be found [here](https://bioconductor.org/packages/release/bioc/html/RNAAgeCalc.html).    
The Bioconductor package vignette can be found [here](https://bioconductor.org/packages/release/bioc/vignettes/RNAAgeCalc/inst/doc/RNAAge-vignette.html).    
An interactive Shiny App can be found [here](https://xuren2120.shinyapps.io/RNAAgeCalcshiny/).    
The Shiny App demo can be found [here](http://www.ams.sunysb.edu/~pfkuan/RNAAgeCalc/instructions.html).    
The Python version of RNAAgeCalc can be found [here](https://pypi.org/project/racpy).

Package installation
------------
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager")    
BiocManager::install("RNAAgeCalc")    
```

Citation
------------
Ren, X., & Kuan, P. F. (2020). RNAAgeCalc: A multi-tissue transcriptional age calculator. PloS one, 15(8), e0237006. 

@article{ren2020rnaagecalc,    
title={RNAAgeCalc: A multi-tissue transcriptional age calculator},    
author={Ren, Xu and Kuan, Pei-Fen},    
journal={PLOS ONE},    
volume={15},    
number={8},    
pages={1-21},    
year={2020},    
publisher={Public Library of Science}    
}
