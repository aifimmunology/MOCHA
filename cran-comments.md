## R CMD check results

There were no ERRORs or WARNINGs. 

There was 1 NOTE:

❯ checking package dependencies ... NOTE
  Imports includes 31 non-default packages.
  Importing from so many packages makes the package vulnerable to any of
  them becoming unavailable.  Move as many as possible to Suggests and
  use conditionally.

MOCHA is a comprehensive analysis package for the fundamental steps of downstream analysis with scATAC data - having many dependencies is not unexpected. Most dependencies are very popular or are found in Bioconductor, and are well maintained per Bioconductor guidelines. 

* This is a new release.
