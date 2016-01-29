Version 0.2.2:  January 29, 2016
-------------------------------------------------------------------------------

General:

+ Updated license from CC BY-NC-SA 3.0 to CC BY-NC-SA 4.0.
+ Internal changes to conform to CRAN policies.

Amino Acid Analysis:

+ Fixed bug where arguments for the `aliphatic()` function were not being
  passed through the ellipsis argument of `aminoAcidProperties()`.
+ Improved amino acid analysis vignette.
+ Added check for correctness of amino acids sequences to `aminoAcidProperties()`.
+ Renamed `AA_TRANS` to `ABBREV_AA`.

Diversity:

+ Added evenness and bootstrap standard deviation to `rarefyDiversity()` 
  output.

Lineage:

+ Added `ExampleTrees` data with example output from `buildPhylipLineage()`.


Version 0.2.1:  December 18, 2015
-------------------------------------------------------------------------------

General:

+ Removed plyr dependency.
+ Added dplyr, lazyeval and stringi dependencies.
+ Added strict requirement for igraph version >= 1.0.0.
+ Renamed `getDNADistMatrix()` and `getAADistMatrix()` to `getDNAMatrix` and 
  `getAAMatrix()`, respectively.
+ Added `getSeqMatrix()` which calculates a pairwise distance matrix for a set 
  of sequences.
+ Modified default plot sizing to be more appropriate for export to PDF 
  figures with 7-8 inch width.
+ Added `multiggplot()` function for performing multiple panel plots.

Amino Acid Analysis:

+ Migrated amino acid property analysis from Change-O CTL to alakazam. 
  Includes the new functions `gravy()`, `bulk()`, `aliphatic()`, `polar()`, 
  `charge()`, `countPatterns()` and `aminoAcidProperties()`.

Annotation:

+ Added support for unusual TCR gene names, such as 'TRGVA*01'.
+ Added removal of 'D' label (gene duplication) from gene names when parsed 
  with `getSegment()`, `getAllele()`, `getGene()` and `getFamily()`.  May be 
  disabled by providing the argument `strip_d=FALSE`.
+ Added `countGenes()` to tabulate V(D)J allele, gene and family usage.

Diversity:

+ Added several functions related to analysis of clone size distributions, 
  including `countClones()`, `estimateAbundance()` and `plotAbundance()`.
+ Renamed `resampleDiversity()` to `rarefyDiversity()` and changed many of
  the internals. Bootstrapping is now performed on an inferred complete
  relative abundance distribution.
+ Added support for inclusion of copy number in clone size determination
  within `rarefyDiversity()` and `testDiversity()`.
+ Diversity scores and confiderence intervals within `rarefyDiversity()`
  and `testDiversity()` are now calculated using the mean and standard 
  deviation of the bootstrap realizations, rather than the median and
  upper/lower quantiles.
+ Added ability to add counts to the legend in `plotDiversityCurve()`.


Version 0.2.0:  June 15, 2015
-------------------------------------------------------------------------------

Initial public release.

General:

+ Added citations for the `citation("alakazam")` command.


Version 0.2.0.beta-2015-05-30:  May 30, 2015
-------------------------------------------------------------------------------

Lineage:

+ Added more error checking to `buildPhylipLineage()`.


Version 0.2.0.beta-2015-05-26:  May 26, 2015
-------------------------------------------------------------------------------

Lineage:

+ Fixed issue where `buildPhylipLineage()` would hang on R 3.2 due to R change 
  request PR#15508.


Version 0.2.0.beta-2015-05-05:  May 05, 2015
-------------------------------------------------------------------------------

Prerelease for review.