# adaFilter

This is an R package for finding mutliple replicating signals under the Partial Conjunction framework. Our method adaFilter contains two procedures: adaFilter Bonferroni and adaFilter BH that can effciently identify signals that replicate in at least r studies


## Installation
```{r}
library(devtools)

install_github("jingshuw/adaFilter")
```

## Example
AdaFilter starts with a p-value matrix of size M * n where M is the number of hypotheses in one study and n is the number of studies. In some high-throughput genetic experiments, it is common that some hypotheses have missing p-values in some studies. AdaFilter allows missing values in the p-value matrix (which are simply NA values).

Testing for replicability, we reject a hypothesis only when the individual hypotheses are nonnull in at least r studies. For instance, for genetic data, we reject a gene only when there is signal in at least r studies. Here, r is a user-specified replicability level, and should be at least 2.

As an example, we simulate p-values for M = 1000 hypotheses in n = 4 studies, and set 10% hypotheses to be true partial conjunction non-nulls with r = 4. In other words, there are 100 hypotheses that are non-null in both studies, and are the true signals that we want to find.
```{r}
library(adaFilter)
set.seed(1)
data <- GenPMat(M = 1000, n = 4, r = 4)
head(data$pvalue.mat)
```
Then, we can run adaFilter to find signals that replicate in both two studies and control FDR at level 0.05. AdaFilter returns a data frame of 5 columns, storing the decision of whether rejecting a hypothesis or not, adaFilter adjusted p-values, the selection and filtering p-values and the adaFilter adjustment numbers (for more details, see reference). 
```{r}
result <- adaFilter(data$pvalue.mat, r = 4)
head(result)
## show indices of true positives
print(which(result$decision == 1 & data$truth.pc == 1))
```
AdaFilter finds 77 true positives.

In contrast, we can compare adaFilter with three other standard multiple testing procedures, that apply BH on three types of partial conjunction p-values. In this example of n = r= 2, the three types of partial conjunction p-values become the same, which are the maximum of the two p-values for each hypothesis. We can run one of them for comparison with adaFilter.
```{r}
result <- ClassicalMTP(data$pvalue.mat, r = 2, alpha = 0.05)
print(which(result$decision))
```
 
