## R CMD check results

There were no ERRORs or WARNINGs. 

There were 5 NOTEs:

❯ checking CRAN incoming feasibility ... [13s] NOTE
  
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    scATAC (29:32, 37:15)
    translational (37:30)
  
  Package has a FOSS license but eventually depends on the following
  package which restricts use:
    CNEr
  
  Suggests or Enhances not in mainstream repositories:
    ArchR
  
  Found the following (possibly) invalid URLs:
    URL: https://github.com/aifimmunology/MOCHA/actions/workflows/R-CMD-check-manual.yml
      From: README.md
      Status: 404
      Message: Not Found

'scATAC' and 'translational' are field-specific terms.

'CNEr' is a Bioconductor package and a dependency of 'TFBSTools' which has
a GPL-2 open-source license.
      
This GitHub repository is currently private, but will be made public when 
the manuscript for MOCHA is released.

ArchR is an optional package available in GitHub. 
MOCHA is designed to allow the use of ArchR as a starting point, though it is 
not required.

❯ checking package dependencies ... NOTE
  Imports includes 31 non-default packages.
  Importing from so many packages makes the package vulnerable to any of
  them becoming unavailable.  Move as many as possible to Suggests and
  use conditionally.

MOCHA is a comprehensive analysis package for the fundamental steps of downstream analysis with scATAC data - having many dependencies is not unexpected. Most dependencies are very popular or are found in Bioconductor, and are well maintained per Bioconductor guidelines. 

* This is a new release.
