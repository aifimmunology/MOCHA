## Resubmission
This is a resubmission. In this version I have:

* Explained acronyms in DESCRIPTION (ATAC)

* Added references describing zero-inflated statistical methods to the DESCRIPTION

* Added missing return value tags to documentation for GRangesToString.Rd and StringsToGRanges.Rd

* Decreased the size of the package tarball to 1.9MB

* Removed foo:::f access to unexported objects in examples

* In examples requiring database packages: Replaced, where possible, dontrun{} with donttest{} and run those examples conditionally

* Replaced print() with message() or warning() where applicable, and make all messages suppressable with a 'verbose' argument

* Replaced default path writing in tests with tempdir()
  
## R CMD check results

There were no ERRORs or WARNINGs. 

There were 2 NOTEs:
  1)
  Possibly misspelled words in DESCRIPTION:
    ATAC (29:50, 37:45)
    Ghazanfar (39:3)
    Pimentel (40:3)
    Spearman's (40:46)
    Su (39:27)
    Transposase (30:3)
    al (39:37)
    et (39:34)

  Suggests or Enhances not in mainstream repositories:
    ArchR
  Availability using Additional_repositories specification:
    ArchR   yes   https://imran-aifi.github.io/drat
  
  2)
  Packages suggested but not available for checking:
    'ArchR', 'TxDb.Hsapiens.UCSC.hg38.refGene'
  
  TxDb is available in Bioconductor.
  
* This is a new release.