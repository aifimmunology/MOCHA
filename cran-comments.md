## R CMD check results

There were no ERRORs* or WARNINGs. 

  Bioconductor does not yet build and check packages for R version 4.3; see
    https://bioconductor.org/install
  
  *This error occurred with winbuilder but may be ignored according to https://github.com/r-hub/rhub/issues/471

There were 2 NOTEs:

Possibly misspelled words in DESCRIPTION:
  ATAC (29:44, 37:27)

Package has a FOSS license but eventually depends on the following
package which restricts use:
  CNEr

We use the GPL v3 License. CNEr on Bioconductor uses the GPL v2 license and
has dependencies with GPL v3 in Bioconductor.

❯ checking package dependencies ... NOTE
  Imports includes 30 non-default packages.
  Importing from so many packages makes the package vulnerable to any of
  them becoming unavailable.  Move as many as possible to Suggests and
  use conditionally.

MOCHA is a comprehensive analysis package for the fundamental steps of downstream analysis with scATAC data - having many dependencies is not unexpected. Most dependencies are very popular or are found in Bioconductor, and are well maintained per Bioconductor guidelines. 

* This is a new release.
