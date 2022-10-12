# inactiveXX: an R package for inferring X-inactivation (Xi) status at the single cell level from single transcriptomes

Infer the X-inactivation status from your single cell transcriptomic data.

## Installation

```
devtools::install_github('constantAmateur/inactiveXX')
```

## Basic usage

The basic steps are:

1. Call heterozygous SNPs
2. Get counts at heterozygous SNPs in single cells, on the X chromosome
3. Infer Xi status from this data
4. Downstream analysis using Xi states.

A rough example of this in code would be:

```
#Find het SNPs by looking at expression at 1k genomes SNPs across all single cells
hSNPs = hetSNPsFromRNA(bams10X,refGenome10X,outputs=sprintf('%s_scRNA_1kSNPs_XCnts.tsv',names(bams10X)))
#Filter just to usable het SNPs on X
XCnts = filterCountsX(hSNPs)
#Run inferance
fit = inferInactiveX(XCnts,nParallel=8)
#Visualise the fit
plotSolutitons(fit)
```

## More detailed instructions

### Step 1 - call heterozygous SNPS

Unless you have some other source of genotype information (e.g. whole genome sequencing), the first step will count reads at sites of 1,000 genomes SNPs on the X chromosome across all cells.  Any SNP that has evidence of both genotypes is then accepted as heterozygous and taken forward.  As counting reads from BAM files is a slow process, this should be done in parallel where possible (by setting `nParellel=8` or similar), and it is a good idea to save the output for future use (by setting the `outputs` parameter).  As well as the BAM file containing the single cell transcriptomic counts, it is necessary to provide the reference genome used to create the BAM file.

It is possible to adjust the default cut-offs that determine when a SNP should be called heterozygous, by setting `errRate` and `alpha`, but don't do this unless you know what you're doing. 

Note, if you obtained heterozygous SNPs from some other source than the scRNA BAM file, you also need to first get the counts at those SNPs in each cell.  Do this by calling `getAllelicExpression` using the heterozygous SNPs you have generated from another source. 

Note also that you should always process **all** BAM files relating to the same individual at the same time, by providing them as a vector to `bams10X`.  The more cells from an individual are available, the more accurate Xi inference (and heterozygous SNP inference) will be.

Typically, the code for this step look something like this:

```
hSNPs = hetSNPsFromRNA(bams10X,refGenome10X,outputs=sprintf('%s_scRNA_1kSNPs_XCnts.tsv',names(bams10X)))
```

### Step 2 - Get counts from heterozygous SNPs on X

Having obtained heterozygous SNPs, these now need to be filtered to the usable ones.  This consists of removing anything with too few observations to be useful, with a high likelihood of escaping Xi, or with direct evidence that both alleles are expressed in a single cell (i.e., we see evidence of Xi escape in the data itself).  We also typically exclude anything that is not mapped to a gene.

Typical code:

```
XCnts = filterCountsX(hSNPs)
```

### Step 3 - Infer Xi status

We now have the data we need to fit the EM model and estimate the Xi status.  This is done by running `inferInactiveX`, but it is worth understanding a few things about how the model works and how to troubleshoot it before blindly proceeding.

`inferInactiveX` will run an EM algorithm to find the X chromosome genotype and Xi status of cells using `nStarts` (1,000 by default) random starts.  As this can be time consuming, it is useful to fit the `nStart` EM models in parallel.  However, due to Rs terrible memory handling, the memory requirement for `nParellel` parallel threads is roughly `nParallel` times bigger.  As such, you may need to use a much lower `nParallel` than the number of cores available, particularly for fits with many cells.

Each EM iteration is itself a multi-part process involving iterating to converge (or until `maxIter` is reached) on smoothed versions of the likelihood space, a process known as deterministic annealing.  This process should help the EM method to converge upon a global maximum rather than a local one.  The details of this multi-part process can be controlled with the `tauInit`, `betaStart` and `betaFac` parameters, but you shouldn't touch them.

Once all `nStart` models have been fit, an extra check is made to ensure that the distribution of results is consistent with the model having found a global maximum.  That is, when a global maximum has been found, we expect that the best fitting `tauDiffFrac` (default 10%) models will all have the same Xi skew (tau).  If these are found to vary by more than `tauDiffThresh` (5% by default), then it is likely that the global maximum has not been found, probably because the Xi skew is more extreme than it is possible to estimate with the available data (or you screwed up and are trying to estimate Xi on a male).  In this case, an error usually occurs, but this can be converted to a warning if `tauDiffWarnOnly` is set.

Finally, the best fit is summarised in various ways in the output.  The primary summary is the probability of each cell being in each Xi state.  `inferInactiveX` will give probabilities of each state, which are converted into either 'Maternal', 'Paternal', or unassigned based on the thresholds set by `logitCut` and `pCut`.  Note that 'Maternal' and 'Paternal' are arbitrary.

Typical code: 

```
fit = inferInactiveX(XCnts,nParallel=8)
```

It is highly recommended to visualise the distribution of the `nStarts` solutions, to visually check that there is a cluster of best fitting models.  This can be done by running `plotSolutions`.  Note that because it is not possible to determine which X copy is maternal, we expect the solutions to be symmetrically distributed around an Xi skew of 50.


### Step 4 - Downstream analysis.

Although there is lots of extra information in the fit object, the main information is the mapping between cell ID and Xi state.  You may of course do whatever you like with this information, but `inactiveXX` contains a few additional functions to help with common uses.  These are briefly described here.

#### Deviations in Xi skew within one individual

The Xi state of a group of cells in an individual is determined by the Xi skew of the founder cells, that is the set of cells present when Xi was first initiated, and the inheritance pattern that relates those founder cells to the observed cells.  On average, a random collection of cells should mirror the founder Xi skew.  Therefore, significant deviations from the founder Xi skew are a sign of a population bottleneck of some kind for those cells, either due to local restriction, clonal expansion, or other factors.

The function `estCelltypeSkew` can be used to create a [miloR](https://github.com/MarioniLab/miloR) object, which tests for regions with significant deviation from the founder Xi skew.  This is done by assuming that the founder Xi skew is approximately equal to the average Xi skew across all cells, although any founder Xi skew (tau) can be specified using the `tau` parameter.  Neighbourhoods around cells are then tested for significant deviation from the founder skew, using the milo framework.

In order to do this, the table of counts for the same cells that Xi state has been estimated needs to be provided.  The most convenient way to do this is by providing a [Seurat](https://satijalab.org/seurat/) object, as all the meta-data can then be passed through to the resulting milo object output by `estCelltypeSkew`.  When providing a Seurat object, the columns of meta-data to keep need to be specified using the `resultsPassthrough` parameter.

One oddity of the resulting milo object is that the `logFC` measure will, by default, actually represent a linear deviation from the founder Xi skew.  That is, if the founder Xi skew is 0.6, a logFC value of 0.2 represents a value of 0.6+0.2 = 0.8.  This behaviour can be disabled by setting `logFC=TRUE`.

Where this function will tend to go wrong, is if the Seurat object (or table of counts) provided to `toc` has cellIDs in a different format to the inactiveXX fit object provided to `fit`.  

A typical code usage would be something like:

```
milo = estCelltypeSkew(srat,fit,resultsPassthrough='annot')
plotDAbeeswarm(milo$res,group.by='annot')
```

#### Population bottleneck size estimation

If you have Xi skew values from the same cell types from multiple individuals, this can be used to estimate the effective population bottleneck size that these cells have passed through.  The assumption here is that for each individual, the same number of founder cells contribute to the observed cells.  If this is true, then the variance in the distribution of Xi skews across individuals will depend on the number of founder cells contributing.

To perform this estimation, use the `estPopSize` function.  The main input to this function is a list of `inactiveXX` fit objects produced by `inferInactiveX`, provided to the `fits` parameter.  

If you wish to calculate the population size of particular groupings of cells, such as cell types, you need to specify the cellIDs belonging to each cell group using the `cellsToUse` function.  This is done by constructing a list, named by the cell groupings, with cellIDs indicating which cells across all samples belong to that cell type.  For example, if you have two cell types A and B, you would specify `cellsToUse` like this:

```
cellsToUse = list(typeA = c("4602STDY6976422_AAACGGGTCCGCAGTG", "4602STDY6976422_AAAGCAACATTTCACT", 
"4602STDY6976422_AAAGCAATCAAGGCTT", "4602STDY6976422_AAAGCAATCAGCAACT", 
"4602STDY6976422_AAAGTAGAGAGCTGCA", "4602STDY6976422_AAATGCCTCTAGAGTC"
), typeB = c("4602STDY6976422_AAATGCCGTCCAGTAT", "4602STDY6976422_AACACGTAGCCCTAAT", 
"4602STDY6976422_AACACGTCATGGTAGG", "4602STDY6976422_AACCATGAGGGCTCTC", 
"4602STDY6976422_AACCGCGCAACCGCCA", "4602STDY6976422_AACCGCGCATTAGCCA"
))
```

This function will output a data.frame with estimates of the population size for each cell group/type (or all cells if cellsToUse=NULL), with confidence intervals around those estimates that are specified using the `distQuants` parameter.

#### Correlation between cell types

The other analysis you can do when you have Xi information from multiple individuals, is to calculate how correlated cell types Xi skew are with one another using the `estCorrelation` function.  The input parameters are similar to `estPopSize`, with fits and cell groupings specified in the same way.

By default, the correlation in raw Xi skews is calculated.  However, as basically all cells are derived from the same set of founder cells, all cells tend to be correlated with one another.  It is often more useful to look at which cell types have correlated deviation from the founder Xi skew.  That is, which cell types consistently move away from the founder Xi skew in the same direction.  This can be done by setting `relativeTo` to `median` or `tau`.  `median` will cause the founder Xi to be calculated as the median value across all cell types (the preferred option), `tau` will simply use the global average in the fit object.  It is also possible to manually specify a fixed value to compare against (e.g. 0.5), by leaving `relativeTo='fixed'` and setting `fixedVal`.

The other consideration is that the most accurate estimate of the uncertainty due to the number of individuals available, is complex and slow to calculate.  By default, a faster but less accurate approximation is used, but if the uncertainty needs to be known to high accuracy, consider setting `exactConf=TRUE` to force the non-approximate uncertainty estimate.

The function returns a data.frame summarising the correlation between cell types and the uncertainty in these estimates.  Usually the most useful way to understand these results is to visualise them using the `plotCorrelations` function, which produces a heatmap of the correlations and their uncertainties.



