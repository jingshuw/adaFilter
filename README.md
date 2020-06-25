# adaFilter

This is an R package for finding mutliple replicating signals under the Partial Conjunction framework. Our method adaFilter contains two procedures: adaFilter Bonferroni and adaFilter BH that can effciently identify signals that replicate in at least r studies


## Installation
```{r}
library(devtools)

install_github("jingshuw/adaFilter")
```

## Example
adaFilter starts with a p-value matrix of size M \times n where M is the number of 
