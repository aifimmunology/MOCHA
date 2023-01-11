MOChA: Model-based single-cell Open Chromatin Analysis
------------------------------------------------------------------------

MOCHA is an R package containing a novel single-cell peak-calling algorithm that leverages single-cell information to determine whether a particular genomic region is open by calculating two measures of intensities, and using these to call peaks via a hierarchical model.

------------------------------------------------------------------------

### Table of Contents

-   [Installation](#installation)
-   [Overview](#overview)
-   [Contact](#contact)
-   [License](#license)


-----------------------------------------------------------------------



## <a name="installation"></a> Installation
Install from binaries (stable release on CRAN):
  
    install.packages("MOCHA")
    
Install from source:

    devtools::install_github("aifimmunology/MOCHA")

Install a specific development branch from source:

    devtools::install_github("aifimmunology/MOCHA", ref = "your_branch_name")

## <a name="overview"></a> Usage Overview

Please view the example usage found in the vignette found under
`vignettes/COVID-walkthrough.html`.

The example usage demonstrates this workflow: 

![](man/figures/workflow.svg)

## <a name="contact"></a> Contact

To contact the developers on issues and feature requests, please contact us via the discussions tab for feature requests, or open issues for any bugs.

## <a name="license"></a> License

MOCHA follows the Allen Institute Software License - full information about the license can be found on the LICENSE file.
