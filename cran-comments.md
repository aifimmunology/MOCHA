This is a resubmission. To reduce runtime of tests, skip_on_cran() statements 
were added to all tests that are dependent on external test data.

## R CMD check results

There were no ERRORs or WARNINGs. 

There were 2 NOTEs:
  1)
  Suggests or Enhances not in mainstream repositories:
    ArchR
  Availability using Additional_repositories specification:
    ArchR   yes   https://imran-aifi.github.io/drat
  
  2)
  Packages suggested but not available for checking:
    'ArchR', 'TxDb.Hsapiens.UCSC.hg38.refGene'
  
  TxDb is available in Bioconductor.
  